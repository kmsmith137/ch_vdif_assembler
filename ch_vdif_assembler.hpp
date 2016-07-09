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
// To use this library, you:
//
//   - Construct an object of type vdif_assembler, which implements some
//     generic logic for handling the above nuisance issues.
// 
//   - Define a subclass of the virtual base class vdif_processor which
//     contains task-specific data processing logic, and register the
//     processor with the assembler (assembler.register_processor()).
//
//   - Create an object of type vdif_stream, which represents the high-speed
//     data stream to be processed.  Currently the following subclasses are
//     implemented:
//
//         vdif_network_stream:    real-time network capture
//         vdif_file_stream:       previous network capture saved on disk
//         vdif_simulated_stream:  for testing (currently doesn't simulate much)
//
//     To get a list of filenames in an acqusition on moose, you may find
//     the utility show-moose-acquisitions.py useful.
//
//   - Call assembler.run()
//
//
// Multiple processors can be registered, and will be run in separate threads.
// For example we could try to run a waterfall_plotter, an FRB search, and a GPU
// pulsar backend at the same time.  The data can also be saved to disk in real
// time in vdif format, or buffered and saved to disk if a processor calls
// vdif_assembler::trigger().
//
// Currently, each assembler instance may only process one stream, and you'll
// get an error if you try to call run() or register_processor() after the assembler
// has finished.  I may generalize this later.


#ifndef _CH_VDIF_ASSEMBLER_HPP
#define _CH_VDIF_ASSEMBLER_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <unistd.h>
#include <string>
#include <vector>
#include <memory>
#include <complex>
#include <iostream>
#include <stdexcept>

#include "ch_vdif_assembler_kernels.hpp"


// Branch predictor hint
#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif

//
// xassert(): like assert, but throws an exception in order to work smoothly with python.
// We also print the exception to stdout, so that we see it regardless of which thread threw it.
//
#ifndef xassert
#define xassert(cond) xassert2(cond, __LINE__)
#define xassert2(cond,line) \
    do { \
        if (_unlikely(!(cond))) { \
	    const char *msg = "Assertion '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")\n"; \
	    std::cout << msg << std::flush; \
	    throw std::runtime_error(msg); \
	} \
    } while (0)
#endif


namespace ch_vdif_assembler {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// Constants


namespace constants {
    static const int chime_nfreq = 1024;
    static const int header_nbytes = 32;
    static const int timestamps_per_packet = 625;
    static const int timestamps_per_frame = 1 << 23;    // cadence of noise source
    static const int packet_nbytes = 5032;              // = header_nbytes + 8 * timestamps_per_packet
    static const int num_disks = 10;
    static const int default_abuf_size = 4;
    static const int default_assembler_nt = 65536;
    static const int cache_line_size = 64;
};


// -------------------------------------------------------------------------------------------------
//
// Helper routines
//
// One thing that would be generally useful is a thread-safe stream class for logging.
// Right now the output from different threads can be interleaved, which can be confusing.
// In the meantime, it helps is to write each line as a single string, i.e. instead of
//   cout << "The value of i is: " << i << endl;
//
// do either
//   cout << (string("The value of i is: ") << to_string(i) << "\n") << flush;
//
// or
//   stringstream ss;
//   ss << "The value of i is: " << i << "\n";
//   cout << ss.str() << flush;


struct noncopyable
{
    noncopyable() { }
    noncopyable(const noncopyable &) = delete;
    noncopyable& operator=(const noncopyable &) = delete;
};


struct vdif_stream;
struct vdif_processor;
struct assembled_chunk;
struct assembled_chunk_pool;
struct assembler_killer;
struct assembler_nerve_center;

// Creates directory, but doesn't throw an exception if it already exists
extern void xmkdir(const std::string &dirname);

extern bool is_empty_dir(const std::string &dirname);

template<typename T>
inline std::string to_string(const T &x)
{
    std::stringstream ss;
    ss << x;
    return ss.str();
}


// -------------------------------------------------------------------------------------------------
//
// The assembler class
//
// FIXME some optional argments supported in v1 are not supported yet in v2 (e.g. trim_frac).
// I'll reinstate them soon.


struct vdif_assembler
{
    std::shared_ptr<assembler_nerve_center> nc;
    std::shared_ptr<assembler_killer> killer;

    //
    // Note: when the constructor is called, "assembler" and "disk_writer" threads
    // are automatically spawned.  The 'write_to_disk' argument determines whether
    // data is written to disk at the same time it's processed.  
    // 
    // If write_to_disk is set to false, but a large value of rbuf_size is chosen,
    // the assembler will ring-buffer the data, and write to disk if a processor
    // calls trigger().
    //
    vdif_assembler(bool write_to_disk=false, 
		   int rbuf_size=constants::num_disks, 
		   int abuf_size=constants::default_abuf_size, 
		   int assembler_nt=constants::default_assembler_nt);

    ~vdif_assembler();

    // Each call to register_processor() spawns one processing thread.
    void register_processor(const std::shared_ptr<vdif_processor> &p);
    
    //
    // Spawns one or more I/O threads, and waits for assembler to finish.
    // One quirk of the current code is: run() will still process the data
    // even if no processing threads are registered and write_to_disk=false!
    //
    void run(const std::shared_ptr<vdif_stream> &s);

    // The asynchronous version of run()
    void start_async(const std::shared_ptr<vdif_stream> &s);
    void wait_until_end();
};


// -------------------------------------------------------------------------------------------------
//
// vdif_stream: a virtual base class which represents a stream of unassembled data.
// Currently four subclasses are defined:
//
//    - network_stream: a real-time capture arriving over the network
//    - file_stream: a previous capture which has been saved to disk
//    - sim_stream: a simulated 6.4 Gpbs capture (right now not much is simulated, but this will change)
//    - timing_stream: a simulated capture which times the assembler and registered processors


struct vdif_stream {
    bool is_realtime;
    
    vdif_stream(bool is_realtime_) : is_realtime(is_realtime_) { }
    virtual ~vdif_stream() { }
    
    virtual void spawn_threads(const std::shared_ptr<assembler_nerve_center> &nc) = 0;
};


//
// The streams are implemented as factory functions returning pointers, in order
// to avoid polluting the .hpp files with details of their implementations.
//
extern std::shared_ptr<vdif_stream> make_file_stream(const std::string &filelist_filename);

extern std::shared_ptr<vdif_stream> make_file_stream(const std::vector<std::string> &filename_list);

extern std::shared_ptr<vdif_stream> make_network_stream();

extern std::shared_ptr<vdif_stream> make_simulated_stream(double gbps=6.4, double nsec=60.0);

extern std::shared_ptr<vdif_stream> make_timing_stream(int npackets_per_chunk, int nchunks);


// -------------------------------------------------------------------------------------------------
//
// vdif_processor: represents a processing task which runs on assembled chunks.  For a reference
// example, see waterfall_plotter.cpp.
//
// Each processor is run in its own thread, which is spawned when vdif_assembler::register_processor()
// is called.  Each processing thread ends when end-of-stream is reached, an unrecoverable error
// occurs, or if the processor throws an exception.


struct vdif_processor : noncopyable {
    std::string name;
    bool is_critical;

    pthread_mutex_t mutex;
    bool runflag;

    //
    // If a processor is 'critical', then an exception thrown in the processor will kill the entire
    // assembler.  Otherwise, when the processor dies, the assembler (and any other registered processors)
    // keeps running until end-of-stream.
    //
    vdif_processor(const std::string &name, bool is_critical=false);
    virtual ~vdif_processor() { }

    bool is_running();
    void set_running();

    //
    // These are the member functions which must be instantiated in order to define a vdif_processor.
    //
    // The assembler calls process_chunk() multiple times, presenting the processor with the appearance
    // of a uniform sequence of assembled data.  When the stream ends, the assembler calls finalize().
    // See below for details on how to use the assembled_chunk!
    //
    virtual void process_chunk(const std::shared_ptr<assembled_chunk> &a) = 0;
    virtual void finalize() = 0;
};


//
// Some processors which are currently implemented, in waterfall_plotter.cpp and rfi_histogrammer.cpp.
//
// We export factory functions returning pointers, in order to avoid polluting the .hpp files with details
// of the processor implementations.
//
// waterfall_plotter.cpp is a good reference for implementing a vdif_processor.
//
std::shared_ptr<vdif_processor> make_waterfall_plotter(const std::string &outdir, bool is_critical=false);

std::shared_ptr<vdif_processor> make_intensity_beam(const std::string &acqdir);

std::shared_ptr<vdif_processor> make_rfi_histogrammer(const std::string &output_hdf5_filename, bool is_critical=false, bool ref_flag=false);


// -------------------------------------------------------------------------------------------------
//
// assembled_chunk: from the perspective of the vdif_processor, this is the class which contains
// the high-speed data.
//
// Each assembled chunk corresponds to range of timestamps [t0,t0+nt), where t0 is a 64-bit
// wraparound-free timestamp. 
//
// WARNING 1: timestamps are always represented by 'int64_t', and it is a bug to convert to 'int'/
// Be careful since it's easy to do this by accident!
//
// WARNING 2: Usually these ranges will be contiguous between calls, e.g.
//   [t0,t0+nt)   [t0+nt,t0+2*nt)   [t0+2*nt,t0+3*nt)   ...
// but the vdif_processor should not assume that this!  If there is a temporary 
// interruption in data stream, then a timestamp gap will appear.
//
// The data is represented as a uint8_t array of shape (constants::chime_nfreq, 2, nt)
// where the middle index is polarization.  Each uint8_t is an offset-encoded (4+4)-bit
// complex number, see offset_decode() below for the conversion to a pair (re,im) where
// re,im are integers between -7 and 7 inclusive.
//
// WARNING 3: Some entries in the data array will be 0x00, which is a special value
// indicating "missing data".  Calling offset_decode() on these entries will return
// (-8-8j) which is probably not what you want.  Handling missing data is an important
// aspect of the vdif_processor since it happens all the time.  If a GPU correlator node 
// is down, which is a frequent occurrence, then some frequencies will be "all missing".
// There are also routine packet loss events on second-timescales which result in some
// high-speed samples being flagged as missing data.
//
// When developing a new vdif_processor, I suggest starting by testing for missing data
// with (if (byte != 0x00) ...) and using offset_decode() to convert to complex.  However
// it's unlikely that a processor which uses this approach will be fast enough to run
// on a real-time network capture!  Unfortunately it seems to be necessary to use
// assembly language kernels.  
//
// You can see an example below in sum16_auto_correlations(),  where I started with a 
// "reference" kernel which is easy to understand, and then wrote an equivalent assembly 
// language kernel which is fast enough for production use.  I left the reference kernel
// in the code as a way of documenting what the assembly language kernel actually does.
// My plan is to keep doing this on an ad hoc basis as new vdif_processors are developed,
// until we end up with a complete set of kernels, which hide the complexity of assembly
// language in member functions of assembled_chunk.
//
// WARNING 4: For real-time network captures, assembly language kernels will probably
// be necessary, feel free to email me if you need one!
//
// A final note: for processors which only need to use downsampled intensity data
// (not high-speed baseband data) the simplest approach is to use the helper class
// downsampled_intensity, see below!


struct assembled_chunk : noncopyable {
    std::shared_ptr<assembled_chunk_pool> pool;

    //
    // Offset-encoded raw data; array of shape (constants::chime_nfreq, 2, nt) 
    // This data is read simultaneously by multiple threads, so don't modify it!
    //
    // "Invalid" entries, e.g. data which never arrived due to dead correlator
    // nodes or dropped packets, are represented by (uint8_t)0.
    //
    const uint8_t *const buf;

    // Timestamp range
    const int64_t t0;
    const int nt;
    
    // Used internally by the assembler, modifying this field will probably crash or deadlock!
    int pcount;

    assembled_chunk(const std::shared_ptr<assembled_chunk_pool> &pool, int64_t t0);
    ~assembled_chunk();

    // These methods are a little slow and intended for unit testing
    bool is_zero() const;
    bool is_equal(const assembled_chunk &a) const;

    //
    // Should probably never be used in "production" code (an assembly language kernel will
    // probably be necessary) but useful during development.
    //
    static inline void offset_decode(int &re, int &im, uint8_t byte)
    {
	re = (int)((byte & 0xf0) >> 4) - 8;
	im = (int)(byte & 0x0f) - 8;
    }

    //
    // This kernel is used in the waterfall plotter.  It returns the sum of all non-missing
    // visibilities |E^2| over 16 timestamps, and a count of the number which were non-missing.
    //
    static inline void sum16_auto_correlations_reference(int &sum, int &count, const uint8_t *buf)
    {
	int re, im;
	sum = count = 0;

	for (int i = 0; i < 16; i++) {
	    if (buf[i] != 0) {
		offset_decode(re, im, buf[i]);
		sum += (re*re + im*im);
		count += 1;
	    }
	}
    }

    // This is the equivalent assembly language version, which is much faster
    static inline void sum16_auto_correlations(int &sum, int &count, const uint8_t *buf)
    {
	// Defined in ch_vdif_assembler_kernels.hpp
	_sum16_auto_correlations(sum, count, buf);
    }

    //
    // This kernel fills two shape-(nfreq,2,nt) arrays:
    //   efield: complex electric field values, with missing data represented by (0+0j)
    //   mask: either 0 or 1 to indicate missing vs non-missing data
    //
    // FIXME: this kernel needs an assembly language version but I haven't written
    // it yet!  (It's called from cython)
    // 
    inline void fill_efield_array_reference(std::complex<float> *efield, int *mask)
    {
	int arr_size = constants::chime_nfreq * 2 * this->nt;
	int re, im;

	for (int i = 0; i < arr_size; i++) {
	    if (buf[i] != 0) {
		offset_decode(re, im, buf[i]);
		efield[i] = std::complex<float> (float(re), float(im));
		mask[i] = 1;
	    }
	    else {
		efield[i] = std::complex<float> (0.0, 0.0);
		mask[i] = 0;
	    }
	}
    }

    //
    // This kernel fills two shape-(nfreq,2,nt) arrays:
    //   vis: real auto correlations, with missing data represented by zeros
    //   mask: either 0 or 1 to indicate missing vs non-missing data
    //
    // FIXME: this kernel needs an assembly language version but I haven't written
    // it yet!  (It's used in rfi_histogrammer.cpp)
    // 
    inline void fill_auto_correlations_reference(float *vis, int *mask)
    {
	int arr_size = constants::chime_nfreq * 2 * this->nt;
	int re, im;

	for (int i = 0; i < arr_size; i++) {
	    if (buf[i] != 0) {
		offset_decode(re, im, buf[i]);
		vis[i] = (float)(re*re + im*im);
		mask[i] = 1;
	    }
	    else {
		vis[i] = 0.0;
		mask[i] = 0;
	    }
	}
    }	

    // also intended for testing
    static std::shared_ptr<assembled_chunk> make_random(const std::shared_ptr<assembled_chunk_pool> &pool, int64_t min_allowed_t0);
};


// -------------------------------------------------------------------------------------------------
//
// A helper class for vdif_processors which process downsampled intensity data rather than high-speed
// baseband.  Rather than handle the assembled_chunk directly, the vdif_processor can pass it to
// downsampled_intensity::process_chunk(), and fill a buffer containing all downsampled intensities
// derived from the assembled_chunk.


struct downsampled_intensity {
    const int nt_downsample;  // downsampling factor between baseband and intensity
    const double dt_sample;   // downsampled sample length, int seconds

    bool initialized;
    int64_t initial_t0;       // downsampled
    int64_t curr_chunk_t0;    // downsampled
    int64_t curr_chunk_nt;    // downsampled

    ssize_t nt_alloc;         // must be at least curr_chunk_nt
    std::vector<float> intensity_buf;
    std::vector<float> weights_buf;

    downsampled_intensity(int nt_downsample);

    //
    // When process_chunk() returns, the following members of 'struct downsampled_intensity' will be initialized
    //    curr_chunk_t0
    //    curr_chunk_nt
    //    intensity_buf     shape (constants::chime_nfreq, 2, curr_chunk_nt)   note that the stride is curr_chunk_nt, not nt_alloc
    //    weights_buf       shape (constants::chime_nfreq, 2, curr_chunk_nt)
    //
    void process_chunk(const std::shared_ptr<assembled_chunk> &a);
    
    // A slow reference version of process_chunk(), only used for testing
    void process_chunk_reference(const std::shared_ptr<assembled_chunk> &a);
};


}  // namespace ch_vdif_assembler

#endif // _CH_VDIF_ASSEMBLER_HPP
