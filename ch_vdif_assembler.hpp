//
// ch_vdif_assembler: a module for analysis of CHIME high-speed data.
//
// Defines an "assembler" class which handles nuisance issues generic
// to high-speed analysis.
//
//   1. Packets arrive out of order, data from different frequencies can
//      arrive at slightly different times, packets can get dropped.
//
//      The assembler keeps data in a ring buffer and presents the outside
//      world with a regular array, with dropped packets represented by mask/
//
//   2. 32-bit FPGA counts wrap around every 3 hours.  The assembler transparently
//      converts to 64-bit timestamps which are wraparound free.
//
//   3. Noise source "edge" every 21 seconds.  The assembler has flags to trim
//      data near the edge, and to split data on noise source edges.
//
//   4. Zeros / railed values.  The way these are handled in the low-level CHIME
//      software may change soon, but for now there are some flags in the
//      assembler which may help.
//
//   5. RFI mask.  The assembler has flags to apply this automatically.
//
// To use this libarary, you:
//
//   - Construct an object of type vdif_assembler, which implements some
//     generic logic for handling the above nuisance issues.
// 
//   - Define a subclass of the virtual base class vdif_processor which
//     contains task-specific data processing logic, and register the
//     processor with the assembler (vdif_assembler::register_processor())
//
//   - Create an object of type vdif_stream, which represents the high-speed
//     data stream to be processed.  Right now the only subclasses of vdif_stream
//     are vdif_file (representing a single file) or vdif_acquisition
//     (representing multiple files).
//
//     To get a list of filenames in an acqusition on moose, you may find
//     the utility show-moose-acquisitions.py useful.
//    

#ifndef _CH_VDIF_ASSEMBLER_HPP
#define _CH_VDIF_ASSEMBLER_HPP

#include <stdint.h>
#include <vector>
#include <complex>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace ch_vdif_assembler {
#if 0
}; // pacify emacs c-mode
#endif

// forward declarations...
struct vdif_header;
struct vdif_stream;
struct vdif_processor;


//
// @mask should be an array of length vdif_assembler::nfreq (=1024)
// @version=0 gets the most recent mask version
//
// This function fills the mask array with 0's (representing bad channels)
// or 1's (good channels).
// 
extern void get_rfi_mask(int *mask, int version);


//
// The 32-bit FPGA counts used as timestamps "wrap around" every 3 hours.
//
// This helper class converts a stream of unsigned 32-bit FPGA counts to
// a stream of signed 64-bit timestamps without wraparound.  It correctly
// handles out-of-order timestamps as long as the "jitter" isn't more than
// 2^31 timestamps (which is an unrealistic case).
//
// WARNING: timestamps are always represented by 'int64_t', and it is a bug to
// convert them to 'int'.  Be careful since it's easy to do this by accident!
//
class timestamp_unwrapper {
private:
    bool     initialized;
    int64_t  last_timestamp;

public:
    timestamp_unwrapper() : initialized(false), last_timestamp(0) { }

    inline int64_t unwrap(uint32_t x)
    {
	if (!initialized) {
	    initialized = true;
	    last_timestamp = (int64_t)((uint64_t)x);   // sets upper 32 bits to zero
	    return last_timestamp;
	}

	int32_t delta = (int32_t)x - (int32_t)last_timestamp;
	last_timestamp += (int64_t)delta;
	return last_timestamp;
    }
};


// -------------------------------------------------------------------------------------------------
//
// The main vdif_assembler class.  
//
// To process real-time data, construct an vdif_assembler, register some processors,
// construct an vdif_stream, and call vdif_assembler::run()!
//
// The current implementation of vdif_assembler could be improved (too much copying) 
// but I don't think I'll need to change the core interface!
//

class vdif_assembler {
public:
    //
    // CHIME frequencies are represented as pairs (ifreq_major, ifreq_minor) 
    // where (0 <= ifreq_major < nfreq major) and (0 <= ifreq_minor < nfreq_minor).
    // The channel index is ifreq = ifreq_major * nfreq_major + ifreq_minor.
    // Channel index 0 is 800 MHz, and channel index (nfreq-1) is 400 MHz.
    //
    // Each row of the vdif_file contains row_nt time samples, for all 8 minor 
    // frequencies at a fixed major frequency and polarization.  Each time sample 
    // is a (4+4) bit complex number, represented as an offset-encoded 8-bit unsigned
    // integer.  The row header contains {ifreq_major, polarization, initial_time}
    // for the frame (see definition of struct vdif_header below)
    //
    static const int row_nt = 625;
    static const int header_nbytes = 32;
    static const int nfreq_major = 128;
    static const int nfreq_minor = 8;
    static const int nfreq = nfreq_major * nfreq_minor;

protected:
    int      buf_nt;       // defines window size in units of fpga counts
    int64_t  buf_t0;       // current buffer position
    bool     initialized;
    bool     empty;
    int      rfi_mask_version;
    double   trim_frac;
    int64_t  trim_cutoff;  // derived from trim_frac

    timestamp_unwrapper tsu;

    std::vector<int> rfi_mask;

    // Sliding buffer
    std::vector<uint8_t>  vdif_flag;      // shape (nfreq, 2, 2*buf_nt)
    std::vector<uint8_t>  vdif_data;      // shape (nfreq, 2, 2*buf_nt)

    std::vector< boost::shared_ptr<vdif_processor> >  registered_processors;

public:
    // I made some of these boolean flags public
    bool     mask_zeros;
    bool     mask_rails;
    bool     mask_rfi;
    bool     span_frames;
    bool     allow_drops;

    //
    // Arguments to the vdif_assembler constructor are as follows.
    //
    //   @mask_zeros: Currently a (0+0j) electric field value is sometimes a "real"
    //      zero, but sometimes represents a dropped packet!  There is a long-term
    //      fix for this planned, but in the meantime the caller has the option to
    //      either treat (0+0j) as zero, or mask it.  It's my understanding that
    //      this type of packet drop is rare, so mask_zeros=false is the default.
    //
    //   @mask_rails: If set, then electric field measurements where either the
    //      real or imaginary part is "railed" will be masked.
    //
    //   @mask_rfi: If set, then one of the hardcoded RFI masks will be applied
    //      (see ch_rfi_assembler::get_rfi_mask())
    //
    //   @rfi_mask_version: Determines version of the RFI mask; only meaningful if
    //      mask_rfi=true.  Setting rfi_mask_version=0 uses the most current mask.
    //
    //   @trim_frac: Currently there is a big level change in the data every 21.47
    //      seconds when the noise source switches on/off.  We mask a small fraction
    //      of the data near the level changes to remove transient effects.
    //
    //   @span_frames: By default (span_frames=false), calls to vdif_processor::process_data()
    //      are split up so that they never span a noise source switch event.  This simplifies
    //      detrending for example.  Set this flag to true if you don't want this behavior.
    //
    //      To determine which noise source frame a given timestamp belongs to, do:
    //          int64_t frame = (fpga_count % (1 << 23));
    //
    //   @buf_nt: Choosing a larger value here may have the benefit of reducing
    //      the number of dropped packets, at the expense of memory footprint,
    //      cache efficiency, and latency.  The default value (64K) means that
    //      the memory footprint of the assembler will be ~2GB, the latency will
    //      be ~0.1 sec, and packet drops almost never happen.
    //    
    //   @allow_drops: Set to false if you want to fail with an error if the assembler
    //      needs to drop packets.  (I don't recommend ever setting this.)
    //
    vdif_assembler(bool mask_zeros=false, bool mask_rails=true, bool mask_rfi=false,
		  int rfi_mask_version=0, double trim_frac=1.0e-5, bool span_frames=false,
		  int buf_nt=65536, bool allow_drops=true);

    void register_processor(const boost::shared_ptr<vdif_processor> &p);
    void run(vdif_stream &vs);

    //
    // The processing logic works as follows.  Whenever a chunk of data is
    // ready, the vdif_assembler calls
    //
    //   vdif_processor::process_data()
    //
    // To get the data, the vdif_processor calls one or more of the following 
    // routines back in the vdif_assembler.  Note that all arrays are 1D-packed.
    //
    //   vdif_assembler::get_mask()
    //      fills shape-(nfreq,2,nt) array
    //      the middle index is the polarization, either 0 or 1
    //      each element is guaranteed to be 0 (=bad) or 1 (=good)
    //
    //   vdif_assembler::get_raw_data()
    //      fills shape-(nfreq,2,nt) array, middle index is polarization
    //      each element is an 8-bit offset-encoded electric field value
    //      elements where mask=0 are guaranteed to equal 0x88 (which is 0+0j offset-encoded)
    //
    //   vdif_assembler::get_complex_data()
    //      fills shape-(nfreq,2,nt) array, middle index is polarization
    //      each element is a complex electric field value, without any encoding
    //      elements where mask=0 are guaranteed to equal 0+0j
    //
    //   vdif_assembler::get_visibilities()
    //      fills shape-(nfreq,4,nt) array
    //      packing order for the middle index is
    //         (X^*X)  (Y^*Y)  (Re X^*Y)  (Im X^*Y)
    //      masked elements are guaranteed to equal 0
    //
    // When vdif_processor::process_data() is called, the vdif_assembler gives the vdif_processor 
    // a range [t0,t0+nt) of timestamps which are ready.
    //
    // IMPORTANT: Usually these ranges will be contiguous between calls, e.g.
    //   [t0,t0+nt)   [t0+nt,t0+2*nt)   [t0+2*nt,t0+3*nt)   ...
    // but the vdif_processor should not assume that this!  If there is a temporary 
    // interruption in data stream, then a timestamp gap will appear.
    //
    // Each of the vdif_assembler::get_*() callbacks takes a timestamp range [t0,t0+nt).
    // This will usually be the full range of timestamps passed to vdif_processor::process_data(),
    // but the vdif processor is free to specify a subrange or read in multiple chunks.
    //
    void get_mask(uint8_t *out, int64_t t0, int nt);
    void get_raw_data(uint8_t *out, int64_t t0, int nt);
    void get_complex_data(std::complex<float> *out, int64_t t0, int nt);
    void get_visibilities(float *out, int64_t t0, int nt);

    // Decode (4+4)-bit offset encoding
    static inline void offset_decode(int &re, int &im, uint8_t byte)
    {
	re = (int)((byte & 0xf0) >> 4) - 8;
	im = (int)(byte & 0x0f) - 8;
    }

    
    
    // helpers for run(), public for cython's benefit
    int _run(vdif_stream &vs);
    void _advance(int nt);
    void _allocate();
    void _deallocate();

    // more catering to cython
    void get_complex_data(float __complex__ *out, int64_t t0, int nt);
    int64_t _get_buf_t0() { return buf_t0; }
};


// -------------------------------------------------------------------------------------------------
//
// "Processor" class
//
// This is just a wrapper for user-defined behavior!
//
// For an example of a processor class written in C++, see waterfall_plotter.cpp


struct vdif_processor {
    //
    // Called by the vdif_assembler whenever a chunk of assembled data is ready.
    //
    // This routine will probably "call back" vdif_assembler::get_*() to get
    // data it needs.  See the long comment in class vdif_assembler for more
    // details.
    //
    virtual void process_data(vdif_assembler &a, int64_t t0, int nt) = 0;

    // Called when the input stream reaches EOF
    virtual void finalize() = 0;

    virtual ~vdif_processor() { }
};


// -------------------------------------------------------------------------------------------------
//
// Virtual base class which represents an input stream to the vdif_assembler
//
// Different flavors of input stream are represented by different subclasses:
//
//   vdif_file: single file on disk
//   vdif_acquisition: multiple files on disk, read in parallel from different threads for speed
//   vdif_rt_stream: real-time network capture (currently in an "alpha" state)


struct vdif_stream {
    //
    // Reads frame header at current position in stream.
    //
    virtual vdif_header read_header() const = 0;

    //
    // Reads offset-encoded data into a shape-(8,625) array with 
    // _column-major_ ordering
    // 
    virtual void read_data(uint8_t *out) const = 0;

    //
    // Reads offset-encoded data into a shape-(8,625) array with 
    // row-major ordering and specified memory stride >= 625 for 
    // the row (frequency) axis.
    //
    virtual void read_data_strided(uint8_t *out, int stride) const = 0;

    // Advance to next frame in stream
    virtual void advance() = 0;

    // Is the stream at EOF?
    virtual bool eof() const = 0;

    virtual ~vdif_stream() { }
};


//
// Each row of the vdif file contains an (8,625) array containing
// electric field measurements at 8 minor frequencies and 625 time
// samples, at a fixed polarization and major frequency.  The "header"
// contains information which identifies the row.
//
struct vdif_header {
    int        freq_major;   // 0 <= freq_major < 128
    int        pol;          // 0 <= pol < 2
    uint32_t   fpga_count;   // 32-bit timestamp of first sample in row
};


//
// vdif_file is a subclass of the virtual base class vdif_stream.
// It defines a stream consisting of a single file.
//
struct vdif_file : public vdif_stream {
    std::string filename;
    int nrows;         // number of frames in file, derived from file size (always 50000 in practice?)
    int current_row;   // current position in stream

    boost::shared_array<uint8_t> data;

    vdif_file(const std::string &filename);

    // Devirtualize vdif_stream
    virtual void advance();
    virtual bool eof() const;
    virtual void read_data(uint8_t *out) const;
    virtual void read_data_strided(uint8_t *out, int stride) const;
    virtual vdif_header read_header() const;
    virtual ~vdif_file() { }

    // returns empty pointer if exception gets thrown anywhere
    static boost::shared_ptr<vdif_file> try_to_read(const std::string &filename);
};


//
// Helper class for vdif_acquisition: context for one I/O thread.  (We use
// multiple threads to read the acquisition files so that we can read from
// multiple physical devices in parallel.)
//
struct vdif_acq_context {    
    // read-only after construction, so not protected by mutex
    std::vector<std::string> filename_list;
    int nfiles;

    // data exchange between IO and assembler threads
    int consumer_fileid;
    int producer_fileid;
    boost::shared_ptr<vdif_file> curr_fptr;
    pthread_mutex_t mutex;
    pthread_cond_t cond_file_produced;
    pthread_cond_t cond_file_consumed;
    
    pthread_t thread;

    vdif_acq_context(const std::vector<std::string> &filename_list);

    // Called by producer (IO) thread
    void add_file(const boost::shared_ptr<vdif_file> &f);
    void wait_for_consumer();

    // Called by consumer (assembler) thread
    boost::shared_ptr<vdif_file> get_file();
};


//
// vdif_acquisition is a subclass of the virtual base class vdif_stream.
//
// It defines a stream consisting of a sequence of files.
//
// The script 'show-moose-acquisitions.py' may be useful for generating
// filename lists for the vdif_acquistion constructor!
//
class vdif_acquisition : public vdif_stream {
protected:
    std::vector< std::string > filename_list;
    std::vector< boost::shared_ptr<vdif_acq_context> > context_list;
    int nthreads;
    int nfiles;

    boost::shared_ptr<vdif_file> curr_fptr;
    int curr_fileid;

    void _construct(const std::vector<std::string> &filename_list, int nthreads);
    void _get_fptr_from_iothread(int it);

public:
    // Construct from list of filenames
    vdif_acquisition(const std::vector<std::string> &filename_list, int nthreads=10);

    // Construct from file which contains a list of filenames (one per line)
    vdif_acquisition(const std::string &list_filename, int nthreads=10);

    // Devirtualize vdif_stream
    virtual void advance();
    virtual bool eof() const;
    virtual void read_data(uint8_t *out) const;
    virtual void read_data_strided(uint8_t *out, int stride) const;
    virtual vdif_header read_header() const;
    virtual ~vdif_acquisition() { }
};


//
// Helper class for vdif_rt_stream.
//
struct vdif_rt_context
{
    vdif_rt_context(int ring_buffer_size = 65536);
    
    int nring;
    uint8_t *packet0;
    int64_t packet_stride;
    
    // packet index currently being processed by producer/consumer
    int64_t ix_producer;
    int64_t ix_consumer;
    bool eof;

    pthread_mutex_t mutex;
    pthread_cond_t cond_packet_produced;
    pthread_cond_t cond_packet_consumed;

    std::vector<uint8_t> allocated_slab;
    pthread_t io_thread;

    uint8_t *advance_producer();
    uint8_t *advance_consumer();
    void set_eof();
};


//
// vdif_rt_stream: reperesents a real-time network stream
//
struct vdif_rt_stream : public vdif_stream
{
    boost::shared_ptr<vdif_rt_context> context;
    uint8_t *current_packet;

    vdif_rt_stream(int ring_buffer_size = 65536);

    // Devirtualize vdif_stream
    virtual void advance();
    virtual bool eof() const;
    virtual void read_data(uint8_t *out) const;
    virtual void read_data_strided(uint8_t *out, int stride) const;
    virtual vdif_header read_header() const;
    virtual ~vdif_rt_stream() { }
};


// -------------------------------------------------------------------------------------------------
//
// Some processing classes that have been useful so far


// waterfall_plotter.cpp
extern boost::shared_ptr<vdif_processor> make_waterfall_plotter(const std::string &outdir);

// rfi_histogrammer.cpp
extern boost::shared_ptr<vdif_processor> make_rfi_histogrammer(const std::string &output_hdf5_filename, bool ref_flag);


}  // namespace ch_vdif_assembler

#endif  // _CH_VDIF_ASSEMBLER_HPP

/*
 * Local variables:
 *  c-basic-offset: 4
 * End:
 */
