#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/epoll.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <errno.h>
#include <fcntl.h>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <boost/make_shared.hpp>

#include "ch_vdif_assembler.hpp"

using namespace std;
using namespace boost;

#ifndef xassert
#define xassert(cond) xassert2(cond, __LINE__)
#define xassert2(cond,line) do { if (!(cond)) throw std::runtime_error("Assevdifon '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")"); } while (0)
#endif

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


// forward declarations
static void *acq_iothread_main(void *opaque_arg);
static void *rt_iothread_main(void *opaque_arg);
static void rt_iothread_main2(vdif_rt_context *context);


// std::vector doesn't provide a member function which guarantees deallocation!
template<typename T> 
static inline void deallocate(std::vector<T> &v)
{
    std::vector<T> w;
    v.swap(w);
}


// -------------------------------------------------------------------------------------------------
//
// VDIF assembler.
//
// We note some important invariants of the vdif_assembler data structure which must be preserved!
//
// Masked elements must be set to (0+0j), or explicitly with the offset-encoding:
//   If vdif_flag[i]=0, then vdif_data[i]=0x88
//
// If empty=true, then vdif_flag[:]=0 and vdif_data[:]=0x88


vdif_assembler::vdif_assembler(bool mask_zeros_, bool mask_rails_, bool mask_rfi_, int rfi_mask_version_,
			       double trim_frac_, bool span_frames_, int buf_nt_, bool allow_drops_)
{
    this->buf_nt = buf_nt_;
    this->buf_t0 = 0;
    this->initialized = false;
    this->empty = true;
    this->mask_zeros = mask_zeros_;
    this->mask_rails = mask_rails_;
    this->mask_rfi = mask_rfi_;
    this->span_frames = span_frames_;
    this->allow_drops = allow_drops_;
    this->rfi_mask_version = rfi_mask_version_;
    this->trim_frac = trim_frac_;

    xassert(buf_nt >= 8192);  // a reasonble minimum value (?)
    xassert(trim_frac >= 0.0 && trim_frac <= 0.1);

    // Samples will be cut if ((fpga_count + trim_cutoff) % 2^23) < (2*trim_cutoff)
    this->trim_cutoff = (int)(0.5 * trim_frac * (1<<23));
    xassert(trim_cutoff >= 0 && trim_cutoff < (1<<22));

    rfi_mask.resize(vdif_assembler::nfreq, 0);
    get_rfi_mask(&rfi_mask[0], rfi_mask_version);

    cerr << "vdif_assembler initialized\n";

    //
    // Show constructor arguments which differ from their default values
    //
    if (mask_zeros)
	cerr << "   mask_zeros=true\n";
    if (!mask_rails)
	cerr << "   mask_rails=false\n";
    if (!mask_rfi)
	cerr << "   mask_rfi=false\n";
    if (rfi_mask_version != 0)
	cerr << "   rfi_mask_version=" << rfi_mask_version << "\n";
    if (trim_frac != 1.0e-4)
	cerr << "   trim_frac=" << trim_frac << ", trim_cutoff=" << trim_cutoff << "\n";
    if (span_frames)
	cerr << "   span_frames=true\n";
    if (buf_nt != 65535)
	cerr << "   buf_nt=" << buf_nt << "\n";
    if (!allow_drops)
	cerr << "   allow_drops=false\n";
}


void vdif_assembler::register_processor(const shared_ptr<vdif_processor> &p)
{
    registered_processors.push_back(p);
}


void vdif_assembler::run(vdif_stream &vs)
{
    cerr << "vdif_assembler: start!\n";
    this->_allocate();

    for (;;) {
	int nt = this->_run(vs);

	if (nt == 0)
	    break;

	for (unsigned int i = 0; i < registered_processors.size(); i++)
	    registered_processors[i]->process_data(*this, buf_t0, nt);

	this->_advance(nt);
    }

    for (unsigned int i = 0; i < registered_processors.size(); i++)
	registered_processors[i]->finalize();

    this->_deallocate();
    cerr << "vdif_assembler: end\n";
}


//
// _run(): this routine contains the main logic of the vdif_assembler.
//
// Reads data from the vdif_stream until the buffer gets full, then returns 
// the number of time samples that the processors should read from the
// buffer.  
//
int vdif_assembler::_run(vdif_stream &vs)
{
    if (!initialized) {
	// 
	// First call to _run(), at beginning of input stream.
	// We put the packet at the buffer midpoint, so that we have 
	// headroom for reordered packets.
	//
	vdif_header h = vs.read_header();
	this->buf_t0 = tsu.unwrap(h.fpga_count) - buf_nt;
	this->initialized = true;
	this->empty = true;
    }

    //
    // Outer loop over packets read from stream
    //
    for (;;) {

	if (vs.eof()) {
	    //
	    // We've reached the end of the input stream, but the processors
	    // still need to process whatever data remains in the buffer.
	    //
	    // FIXME: currently we return (2*buf_nt), telling the processors
	    // to read the entire buffer at once.  This can result in a one-time
	    // doubling of the memory footprint, since the usual return value
	    // is (buf_nt).  
	    //
	    // Note that the "big read" of size (2*buf_nt) has the side effect
	    // of setting the 'empty' flag (see _advance() below), which means
	    // that the next call to _run() will return 0, triggering the exit
	    // path in run().  This logic for reaching the exit path is a little
	    // convoluted!
	    //
	    // To do: clean up the above (plan is to define a buf_t1 field and
	    // track the "live" range of the buffer more intelligently)
	    //

	    if (empty)
		return 0;

	    int ret = 2*buf_nt;

	    if (!span_frames) {
		int fmax = (1<<23) - (buf_t0 % (1<<23));
		ret = min(ret, fmax);
	    }

	    return ret;
	}

	vdif_header h = vs.read_header();

	// Beginning of data
	int64_t data_t0 = tsu.unwrap(h.fpga_count);

	if (data_t0 < buf_t0) {
	    //
	    // If we get here, the packet extends past the beginning of the buffer, so
	    // we need to drop it.  This can happen if packets are badly reordered.
	    // 
	    // FIXME: For now we drop packets silently; it would be good to add some 
	    // informational messages such as periodically printing the drop rate. 
	    // We also drop the entire packet if any part of it precedes the buffer 
	    // start time; this could be improved but it's not a high priority.
	    //
	    if (!allow_drops)
		throw runtime_error("vdif_assembler: badly reordered packets, either increase buffer size or set allow_drops=true");

	    vs.advance();
	    continue;
	}

	int64_t buf_t1 = buf_t0 + 2*buf_nt;                  // end of buffer
	int64_t data_t1 = data_t0 + vdif_assembler::row_nt;   // end of data

	if (data_t1 > buf_t1) {
	    // If we get here, the packet extends past the end of the buffer.

	    if (!empty) {
		//
		// Typical case, processors just need to extract some data to make room.
		// We need to tell the processors how much data to extract.
		//

		int ret = data_t1 - buf_t1;    // minimum which must be extracted to make room for packet
		ret = max(ret, buf_nt);        // ... but always extract at least buf_nt for efficiency
		ret = min(ret, 2*buf_nt);      // ... but don't try to extract past end of buffer of course!

		if (!span_frames) {
		    int fmax = (1<<23) - (buf_t0 % (1<<23));
		    ret = min(ret, fmax);
		}

		return ret;
	    }

	    // 
	    // Since the buffer is empty, we can just shift it to accommodate the packet.
	    // As above, we put the packet at the buffer midpoint, so that we have room
	    // for reordered packets which may come later.
	    //
	    buf_t0 = data_t0 - buf_nt;
	    buf_t1 = buf_t0 + 2*buf_nt;
	    
	    //
	    // Catch an unlikely corner case where a "mega-packet" is comparable or 
	    // larger than the buffer size.  This should never happen and presumably
	    // indicates data corruption somewhere, so we throw an exception so that
	    // a human can sort things out.  (If this corner case occurs and we don't
	    // trap it here, then we'll segfault later!)
	    //
	    if (data_t1 > buf_t1)
		throw runtime_error("Mega-packet detected, this should never happen!  See ch_vdif_assembler.cpp");
	}

	//
	// If we get here, we're ready to read the packet into the buffer.
	//   n = initial offset in data (or flag) buffer, in bytes
	//   s = minor frequency "stride" in data (or flag) buffer, in bytes
	//
	// i.e. the read is to memory locations
	//   vdif_data[n + i*s + j]     where (0 <= i < nfreq_minor) and (0 <= j < row_nt)
	//
	int n = (2*h.freq_major + h.pol) * (2*buf_nt) + (data_t0 - buf_t0);   // initial offset in data/flag buffer
	int s = (2*vdif_assembler::nfreq_major) * (2*buf_nt);                  // frequency "stride" in data/flag buffer

	//
	// Read data and update mask.
	//
	vs.read_data_strided(&vdif_data[n], s);

	for (int i = 0; i < vdif_assembler::nfreq_minor; i++)
	    memset(&vdif_flag[i*s+n], 0x01, vdif_assembler::row_nt);

	//
	// Remaining logic in this routine handles the zero, rail, and frame trim masks.
	// 
	// Reminder: be careful to preserve the data structure invariant 
	//    flag[i]==0 => data[i]==0x88
	//
	// FIXME: the logic below could be made more cache-efficient
	//

	if (mask_zeros) {
	    for (int i = 0; i < vdif_assembler::nfreq_minor; i++) {
		for (int j = 0; j < vdif_assembler::row_nt; j++) {
		    if (vdif_data[n+i*s+j] == 0x88)   // (0+0j) offset-encodes to 0x88
			vdif_flag[n+i*s+j] = 0;
		}
	    }
	}
	
	if (mask_rails) {
	    for (int i = 0; i < vdif_assembler::nfreq_minor; i++) {
		for (int j = 0; j < vdif_assembler::row_nt; j++) {
		    uint8_t d = vdif_data[n+i*s+j];

		    //
		    // d is "railed" if its real or complex part is +/-7 or 8.
		    // This is equivalent to the following test on the offset-encoded value.
		    //
		    if (((d & 0x0e) == 0) || ((d & 0x0f) == 0x0f) || ((d & 0xe0) == 0) || ((d & 0xf0) == 0xf0)) {
			vdif_data[n+i*s+j] = 0x88;
			vdif_flag[n+i*s+j] = 0;
		    }
		}
	    }
	}	    
    
	// Frame-trimming mask is ((u+j) % (1<<23)) < 2*trim_cutoff
	int u = (data_t0 + trim_cutoff) % (1<<23);

	// A fast test which usually skips the noise source cutiff
	if ((u < 2*trim_cutoff) || ((u + vdif_assembler::row_nt) >= (1<<23))) {
	    for (int j = 0; j < vdif_assembler::row_nt; j++) {
		int v = (u+j) % (1<<23);
		if (v >= 2*trim_cutoff) 
		    continue;
		
		for (int i = 0; i < vdif_assembler::nfreq_minor; i++) {
		    vdif_data[n+i*s+j] = 0x88;
		    vdif_flag[n+i*s+j] = 0;
		}
	    }
	}

	this->empty = false;
	vs.advance();
    }
}


//
// Helper for vdif_assembler::advance()
//
// @p points to a buffer of shape (m,n).  We advance the buffer by k steps 
// along the second axis, and intiailize new values to 'c':
//      p[:,0:(n-k)] = p[:,k:n]
//      p[:,(n-k):n] = c
//
static inline void _advance_helper(uint8_t *p, int k, int m, int n, int c)
{
    k = min(k, n);

    for (int i = 0; i < m; i++) {
	memmove(&p[i*n], &p[i*n+k], n-k);
	memset(&p[i*n+n-k], c, k);
    }
}


void vdif_assembler::_advance(int nt)
{
    xassert(nt >= 0);

    if (empty || nt==0)
	return;

    _advance_helper(&vdif_flag[0], nt, 2*vdif_assembler::nfreq, 2*buf_nt, 0);
    _advance_helper(&vdif_data[0], nt, 2*vdif_assembler::nfreq, 2*buf_nt, 0x88);
    buf_t0 += nt;

    if (nt >= 2*buf_nt)
	empty = true;
}


void vdif_assembler::_allocate()
{
    this->initialized = false;
    this->empty = true;
    this->buf_t0 = 0;
    
    vdif_flag.resize(nfreq * 2 * (2*buf_nt), 0);
    vdif_data.resize(nfreq * 2 * (2*buf_nt), 0x88);
}


void vdif_assembler::_deallocate()
{
    this->initialized = false;
    this->empty = true;
    this->buf_t0 = 0;
    
    deallocate(vdif_flag);
    deallocate(vdif_data);
}


// fills array of shape (nfreq, 2, nt)
void vdif_assembler::get_mask(uint8_t *out, int64_t t0, int nt)
{
    xassert(t0 >= buf_t0);
    xassert(t0+nt <= buf_t0 + 2*buf_nt);
    
    for (int i = 0; i < 2*vdif_assembler::nfreq; i++)
	memcpy(&out[i*nt], &vdif_flag[i*2*buf_nt + (t0-buf_t0)], nt);
}

// fills array of shape (nfreq, 2, nt)
void vdif_assembler::get_raw_data(uint8_t *out, int64_t t0, int nt)
{
    xassert(t0 >= buf_t0);
    xassert(t0+nt <= buf_t0 + 2*buf_nt);
    
    for (int i = 0; i < 2*vdif_assembler::nfreq; i++)
	memcpy(&out[i*nt], &vdif_data[i*2*buf_nt + (t0-buf_t0)], nt);
}

// fills array of shape (nfreq, 2, nt)
void vdif_assembler::get_complex_data(complex<float> *out, int64_t t0, int nt)
{
    xassert(t0 >= buf_t0);
    xassert(t0+nt <= buf_t0 + 2*buf_nt);

    const uint8_t *d = &vdif_data[0];

    for (int i = 0; i < 2*vdif_assembler::nfreq; i++) {
	for (int j = 0; j < nt; j++) {
	    int s = i*(2*buf_nt) + (t0-buf_t0) + j;

	    int re, im;
	    vdif_assembler::offset_decode(re, im, d[s]);
	    out[i*nt+j] = complex<float>(re, im);
	}
    }
}

// fills array of shape (nfreq, 4, nt)
void vdif_assembler::get_visibilities(float *out, int64_t t0, int nt)
{
    xassert(t0 >= buf_t0);
    xassert(t0+nt <= buf_t0 + 2*buf_nt);

    const uint8_t *d = &vdif_data[0];

    for (int i = 0; i < vdif_assembler::nfreq; i++) {
	for (int j = 0; j < nt; j++) {
	    int xre, xim, yre, yim;
	    vdif_assembler::offset_decode(xre, xim, d[(2*i)*(2*buf_nt) + (t0-buf_t0) + j]);
	    vdif_assembler::offset_decode(yre, yim, d[(2*i+1)*(2*buf_nt) + (t0-buf_t0) + j]);

	    int xauto = xre*xre + xim*xim;
	    int yauto = yre*yre + yim*yim;
	    int cross_re = xre*yre + xim*yim;   // Re(x^* y)
	    int cross_im = xre*yim - xim*yre;   // Im(x^* y)

	    out[(4*i)*nt + j] = (float)xauto;
	    out[(4*i+1)*nt + j] = (float)yauto;
	    out[(4*i+2)*nt + j] = (float)cross_re;
	    out[(4*i+3)*nt + j] = (float)cross_im;
	}
    }
}


void vdif_assembler::get_complex_data(float __complex__ *out, int64_t t0, int nt)
{
    this->get_complex_data(reinterpret_cast<std::complex<float> *> (out), t0, nt);
}


// -------------------------------------------------------------------------------------------------


vdif_file::vdif_file(const string &filename_)
    : filename(filename_), nrows(0), current_row(0)

{
    static const int io_block_size = 8192;

    do {
	// Use stringstream to write informational message non-incrementally,
	// to decrease chance of interleaving in multithreaded case
	stringstream ss;
	ss << "reading " << filename << "\n";
	cerr << ss.str().c_str();
    } while (0);

    struct stat s;
    int err = stat(filename.c_str(), &s);

    if (err < 0) {
	stringstream ss;
	ss << filename << ": stat() failed: " << strerror(errno);
	throw runtime_error(ss.str().c_str());
    }
    
    int nbytes_per_row = vdif_assembler::header_nbytes + vdif_assembler::nfreq_minor * vdif_assembler::row_nt;

    if (s.st_size % nbytes_per_row) {
	cerr << filename << ": warning: file size (=" << s.st_size << ")"
	     << " is not divisible by bytes_per_row (=" << nbytes_per_row << "),"
	     << " truncating file\n";
	// we now fall through here, instead of throwing an exception...
    }

    this->nrows = s.st_size / nbytes_per_row;
    xassert(nrows > 0);  // zero-length files create inconvenient corner cases in e.g. vdif_acquisition::advance()

    // Current implementation reads the whole file up-front!
    this->data = shared_array<uint8_t> (new uint8_t[s.st_size]);

    int pos = 0;
    int sz = s.st_size;
    int fd = open(filename.c_str(), O_RDONLY);

    if (fd < 0) {
	stringstream ss;
	ss << filename << ": open() failed: " << strerror(errno);
	throw runtime_error(ss.str().c_str());
    }

    while (pos < sz) {
	int count = min(sz - pos, io_block_size);
	int ret = read(fd, &data[pos], count);

	if (ret <= 0) {
	    stringstream ss;
	    ss << filename << ": error in read(), or unexpected end-of-file: " << strerror(errno);
	    throw runtime_error(ss.str().c_str());
	}
	
	pos += ret;
    }

    close(fd);
}


// virtual member function (devirtualizing base class vdif_stream)
void vdif_file::advance()
{
    if (current_row < nrows)
	current_row++;
}

// virtual member function (devirtualizing base class vdif_stream)
void vdif_file::read_data(uint8_t *out) const
{
    xassert(current_row < nrows);
    xassert(data != NULL);

    int data_nbytes = vdif_assembler::nfreq_minor * vdif_assembler::row_nt;
    int data_location = (current_row * data_nbytes) + ((current_row+1) * vdif_assembler::header_nbytes);
    memcpy(out, &data[data_location], data_nbytes);
}


// virtual member function (devirtualizing base class vdif_stream)
void vdif_file::read_data_strided(uint8_t *out, int stride) const
{
    xassert(current_row < nrows);
    xassert(data != NULL);

    int data_nbytes = vdif_assembler::nfreq_minor * vdif_assembler::row_nt;
    int data_location = (current_row * data_nbytes) + ((current_row+1) * vdif_assembler::header_nbytes);
    const uint8_t *q = &data[data_location];

    for (int it = 0; it < vdif_assembler::row_nt; it++)
	for (int f = 0; f < vdif_assembler::nfreq_minor; f++)
	    out[f*stride + it] = q[it*vdif_assembler::nfreq_minor + f];
}


bool vdif_file::eof() const
{
    return (current_row >= nrows);
}


// virtual member function (devirtualizing base class vdif_stream)
vdif_header vdif_file::read_header() const
{
    xassert(current_row < nrows);
    xassert(data != NULL);

    int data_nbytes = vdif_assembler::nfreq_minor * vdif_assembler::row_nt;
    int header_location = current_row * (data_nbytes + vdif_assembler::header_nbytes);

    const uint32_t *p = reinterpret_cast<const uint32_t *> (&data[header_location]);
    uint32_t pol = (p[3] & 0x0000fffff);
    uint32_t link_id = (p[3] & 0x000f0000) >> 16;
    uint32_t slot_id = (p[3] & 0x01f00000) >> 20;
    uint32_t fpga_count = p[5];

    // pathfinder
    xassert(link_id >= 0 && link_id < 8);
    xassert(slot_id >= 0 && slot_id < 16);

    vdif_header ret;
    ret.fpga_count = fpga_count;
    ret.freq_major = slot_id + 16*link_id;
    ret.pol = pol;

    return ret;
}


// static member function
shared_ptr<vdif_file> vdif_file::try_to_read(const string &filename_)
{
    try {
	return boost::make_shared<vdif_file> (filename_);
    }
    catch (...) {
	return shared_ptr<vdif_file>();
    }
}


// -------------------------------------------------------------------------------------------------


vdif_acq_context::vdif_acq_context(const vector<string> &filename_list_)
    : filename_list(filename_list_), nfiles(filename_list_.size()), 
      consumer_fileid(-1), producer_fileid(0)
{
    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond_file_produced, NULL);
    pthread_cond_init(&cond_file_consumed, NULL);
}


void vdif_acq_context::add_file(const shared_ptr<vdif_file> &f)
{
    pthread_mutex_lock(&mutex);

    while (curr_fptr)
	pthread_cond_wait(&cond_file_consumed, &mutex);

    curr_fptr = f;
    producer_fileid++;
    pthread_cond_signal(&cond_file_produced);

    pthread_mutex_unlock(&mutex);	    
}


void vdif_acq_context::wait_for_consumer()
{
    pthread_mutex_lock(&mutex);

    while (curr_fptr)
	pthread_cond_wait(&cond_file_consumed, &mutex);

    pthread_mutex_unlock(&mutex);
}


shared_ptr<vdif_file> vdif_acq_context::get_file()
{
    pthread_mutex_lock(&mutex);

    consumer_fileid++;
    while (producer_fileid <= consumer_fileid)
	pthread_cond_wait(&cond_file_produced, &mutex);

    shared_ptr<vdif_file> ret = curr_fptr;
    curr_fptr = shared_ptr<vdif_file> ();
    pthread_cond_signal(&cond_file_consumed);

    pthread_mutex_unlock(&mutex);
    return ret;
}


// -------------------------------------------------------------------------------------------------


vdif_acquisition::vdif_acquisition(const vector<string> &filename_list_, int nthreads_)
{
    this->_construct(filename_list_, nthreads_);
}


vdif_acquisition::vdif_acquisition(const std::string &list_filename, int nthreads_)
{
    //
    // Very simple parsing which expects one filename per line, with no extra whitespace
    //
    // FIXME could be improved, but low priority
    //
    ifstream f(list_filename.c_str());
    vector<string> filename_list_;
    string line;
    
    if (f.fail()) {
	stringstream s;
	s << "vdif_acquistion: couldn't open file '" << list_filename << "'";
	throw runtime_error(s.str().c_str());
    }

    while (getline(f, line))
	filename_list_.push_back(line);
    
    this->_construct(filename_list_, nthreads_);
}


void vdif_acquisition::_construct(const vector<string> &filename_list_, int nthreads_)
{
    xassert(context_list.size() == 0);   // detects double call to _construct()

    this->filename_list = filename_list_;
    this->nthreads = nthreads_;
    this->nfiles = filename_list_.size();
    this->curr_fptr = shared_ptr<vdif_file> ();
    this->curr_fileid = -1;             // so that advance() will read first file

    xassert(nthreads > 0);
    xassert(nfiles > 0);

    if (nthreads > nfiles) {
	cerr << "vdif_acquisition: number of threads (=" << nthreads << ")"
	     << " will be reduced to number of files (=" << nfiles << ")\n";
	nthreads = nfiles;
    }

    if (nthreads > 1)
	cerr << "vdif_acqusition: spawning " << nthreads << " IO threads\n";

    context_list.resize(nthreads);

    for (int it = 0; it < nthreads; it++) {
	// round robin assignment (as assumed in vdif_acquisition::advance())
	vector<string> filename_sublist;
	for (int f = it; f < nfiles; f += nthreads)
	    filename_sublist.push_back(filename_list[f]);

	context_list[it] = boost::make_shared<vdif_acq_context> (filename_sublist);
	
	// bare pointer to shared pointer!
	shared_ptr<vdif_acq_context> *p = new shared_ptr<vdif_acq_context> (context_list[it]);

	int err = pthread_create(&context_list[it]->thread, NULL, acq_iothread_main, reinterpret_cast<void *> (p));

	if (err != 0) {
	    cerr << "vdif_acquisition: pthread_create() failed!\n";
	    delete p;
	    throw runtime_error("vdif_acquisition: pthread_create() failed!\n");
	}
    }

    cerr << "vdif_acquisition: entering first blocking call (waiting on iothread 0)\n";
    this->advance();
    cerr << "vdif_acquisition: first blocking call returned\n";
}


void vdif_acquisition::advance()
{
    if (curr_fptr) {
	curr_fptr->advance();
	if (!curr_fptr->eof())
	    return;
    }

    if (curr_fileid >= nfiles)
	return;

    for (;;) {
	// try to advance curr_fileid -> curr_fileid+1
	curr_fptr = shared_ptr<vdif_file>();
	curr_fileid++;

	if (curr_fileid >= nfiles)
	    return;

	int thread_id = curr_fileid % nthreads;   // round-robin
	curr_fptr = context_list[thread_id]->get_file();

	// get_file() can return empty pointer, if there is an error in the IO thread
	if (curr_fptr) {
	    cerr << "processing " << filename_list[curr_fileid] << "\n";
	    return;
	}

	// an error occured, skip this file and go to the next
	cerr << filename_list[curr_fileid] << ": error reading file, this one will be skipped\n";
    }
}


bool vdif_acquisition::eof() const
{
    return (curr_fileid >= nfiles);
}


void vdif_acquisition::read_data(uint8_t *out) const
{
    xassert(curr_fptr);
    curr_fptr->read_data(out);
}

void vdif_acquisition::read_data_strided(uint8_t *out, int stride) const
{
    xassert(curr_fptr);
    curr_fptr->read_data_strided(out, stride);
}

vdif_header vdif_acquisition::read_header() const
{
    xassert(curr_fptr);
    return curr_fptr->read_header();
}


static void *acq_iothread_main(void *opaque_arg)
{
    // The opaque argument is a bare pointer to a shared pointer
    shared_ptr<vdif_acq_context> *pp = reinterpret_cast<shared_ptr<vdif_acq_context> *> (opaque_arg);
    shared_ptr<vdif_acq_context> p = *pp;
    delete pp;

    for (int i = 0; i < p->nfiles; i++) {
	p->add_file(vdif_file::try_to_read(p->filename_list[i]));
	p->wait_for_consumer();   // reduces memory footprint
    }

    return NULL;
}


// -------------------------------------------------------------------------------------------------


vdif_rt_context::vdif_rt_context(int ring_buffer_size)
{
    // Using less than this would be insane
    xassert(ring_buffer_size >= 1024);

    this->nring = ring_buffer_size;
    this->packet_stride = vdif_assembler::header_nbytes + vdif_assembler::nfreq_minor * vdif_assembler::row_nt;
    this->allocated_slab.resize(nring * packet_stride, 0);
    this->packet0 = &allocated_slab[0];

    this->ix_producer = -1;   // first call to advance_producer() will return packet 0
    this->ix_consumer = -1;   // first call to advance_consumer() will block and return packet 0
    this->eof = false;

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond_packet_produced, NULL);
    pthread_cond_init(&cond_packet_consumed, NULL);
}


uint8_t *vdif_rt_context::advance_producer()
{
    pthread_mutex_lock(&mutex);

    ix_producer++;
    pthread_cond_signal(&cond_packet_produced);

    while (ix_consumer <= ix_producer - nring)
	pthread_cond_wait(&cond_packet_consumed, &mutex);

    pthread_mutex_unlock(&mutex);
    
    // OK if this is outside lock
    return packet0 + packet_stride * (ix_producer % nring);
}


uint8_t *vdif_rt_context::advance_consumer()
{
    pthread_mutex_lock(&mutex);

    ix_consumer++;
    pthread_cond_signal(&cond_packet_consumed);

    for (;;) {
	if (ix_producer > ix_consumer) {
	    pthread_mutex_unlock(&mutex);

	    // OK if this is outside lock
	    return packet0 + packet_stride * (ix_consumer % nring);
	}

	if (eof) {
	    pthread_mutex_unlock(&mutex);
	    return NULL;
	}

	pthread_cond_wait(&cond_packet_produced, &mutex);
    }
}


void vdif_rt_context::set_eof()
{
    pthread_mutex_lock(&mutex);
    
    eof = true;
    pthread_cond_signal(&cond_packet_produced);

    pthread_mutex_unlock(&mutex);
}



// -------------------------------------------------------------------------------------------------


vdif_rt_stream::vdif_rt_stream(int ring_buffer_size)
{
    this->context = boost::make_shared<vdif_rt_context> (ring_buffer_size);
    
    // Bare pointer to a shared pointer (IO thread will get reference and call delete)
    shared_ptr<vdif_rt_context> *p = new shared_ptr<vdif_rt_context> (context);

    cerr << "vdif_assembler spawning and waiting for IO thread\n";

    int err = pthread_create(&context->io_thread, NULL, rt_iothread_main, reinterpret_cast<void *> (p));

    if (err != 0) {
	cerr << "vdif_assembler: pthread_create() failed!\n";
	delete p;
	throw runtime_error("vdif_assembler: pthread_create() failed!\n");
    }
    
    this->current_packet = context->advance_consumer();

    if (current_packet)
	cerr << "vdif_assembler is receiving packets from IO thread!\n";
    else
	cerr << "vdif_assembler did not receive any packets from IO thread\n";
}


void vdif_rt_stream::advance()
{
    this->current_packet = context->advance_consumer();
}


bool vdif_rt_stream::eof() const
{
    return (current_packet == NULL);
}


//
// FIXME clean up some cut-and-paste with vdif_file
//
vdif_header vdif_rt_stream::read_header() const
{
    xassert(current_packet != NULL);

    const uint32_t *p = reinterpret_cast<const uint32_t *> (current_packet);
    uint32_t pol = (p[3] & 0x0000fffff);
    uint32_t link_id = (p[3] & 0x000f0000) >> 16;
    uint32_t slot_id = (p[3] & 0x01f00000) >> 20;
    uint32_t fpga_count = p[5];

    vdif_header ret;
    ret.fpga_count = fpga_count;
    ret.freq_major = slot_id + 16*link_id;
    ret.pol = pol;

    return ret;
}


void vdif_rt_stream::read_data(uint8_t *out) const
{
    xassert(current_packet != NULL);

    const uint8_t *in = current_packet + vdif_assembler::header_nbytes;
    memcpy(out, in, vdif_assembler::nfreq_minor * vdif_assembler::row_nt);
}


void vdif_rt_stream::read_data_strided(uint8_t *out, int stride) const
{    
    xassert(current_packet != NULL);

    const uint8_t *in = current_packet + vdif_assembler::header_nbytes;

    for (int it = 0; it < vdif_assembler::row_nt; it++)
	for (int f = 0; f < vdif_assembler::nfreq_minor; f++)
	    out[f*stride + it] = in[it*vdif_assembler::nfreq_minor + f];
}


// -------------------------------------------------------------------------------------------------


static void *rt_iothread_main(void *opaque_arg)
{
    cerr << "vdif_assembler IO thread starting\n";

    // The opaque argument is a bare pointer to a shared pointer
    shared_ptr<vdif_rt_context> *pp = reinterpret_cast<shared_ptr<vdif_rt_context> *> (opaque_arg);
    shared_ptr<vdif_rt_context> p = *pp;
    delete pp;

    try {
	rt_iothread_main2(p.get());
    } catch (...) {
	cerr << "vdif_assembler: network IO thread crashed\n";
    }

    // called even if an exception is thrown
    p->set_eof();

    cerr << "vdif_assembler IO thread terminating\n";
    return NULL;
}


static void rt_iothread_main2(vdif_rt_context *context)
{
    static const int max_events = 100;
    struct epoll_event events[max_events];
    struct epoll_event ev;

    int epoll_fd = epoll_create(10);

    if (epoll_fd < 0) {
	cerr << "vdif_assembler: epoll_create() failed in network IO thread\n";
	throw runtime_error("epoll_create failed");
    }

    int sock_fd = socket(AF_INET, SOCK_DGRAM | SOCK_NONBLOCK, IPPROTO_UDP);

    if (sock_fd < 0) {
	cerr << "vdif_assembler: socket() failed in network IO thread\n";
	throw runtime_error("socket failed");
    }

    struct sockaddr_in server_address;
    memset(&server_address, 0, sizeof(server_address));

    server_address.sin_family = AF_INET;
    inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
    server_address.sin_port = htons(10251);

    if (bind(sock_fd, (struct sockaddr *) &server_address, sizeof(server_address)) < 0) {
	cerr << "vdif_assembler: bind() failed in network IO thread\n";
	throw runtime_error("bind failed");
    }

    int n = 256 * 1024 * 1024;
    if (setsockopt(sock_fd, SOL_SOCKET, SO_RCVBUF, (void *) &n, sizeof(n)) < 0) {
	cerr << "vdif_assembler: setsockopt() failed in network IO thread\n";
	throw runtime_error("setsockopt failed");
    }

    ev.events = EPOLLIN;
    ev.data.fd = sock_fd;
    if (epoll_ctl(epoll_fd, EPOLL_CTL_ADD, sock_fd, &ev) < 0) {
	cerr << "vdif_assembler: epoll_ctl() failed in network IO thread\n";
	throw runtime_error("epoll_ctl() failed");
    }

    uint8_t *buf = context->advance_producer();
    int packet_size = vdif_assembler::header_nbytes + vdif_assembler::nfreq_minor * vdif_assembler::row_nt;

    for (;;) {
        int num_events = epoll_wait(epoll_fd, events, max_events, -1);

        if (num_events < 0) {
	    cerr << "vdif_assembler: epoll_wait() failed in network IO thread\n";
	    throw runtime_error("epoll_wait failed");
	}
	 
        for (int i = 0; i < num_events; i++) {
            if (events[i].data.fd != sock_fd)
                continue;

            ssize_t bytes_read = read(sock_fd, buf, packet_size);

	    // FIXME silently drop?
            if (bytes_read != packet_size)
		continue;

	    buf = context->advance_producer();
        }
    }
    
    // FIXME how to exit gracefully and call set_eof()?
}



// -------------------------------------------------------------------------------------------------
//
// RFI mask


//
// Adds a range of bad frequency channels to a channel mask
//   @mask = an array of length vdif_assembler::nfreq, containing 0's or 1's
//   @freq_min, @freq_max should be in MHz
//
// Note that the mask uses CHIME indexing: channel 0 is the highest frequency 
// (800 MHz) and channel 1023 is the lowest (400 MHz)
//
static void extend_rfi_mask(int *mask, double freq_min, double freq_max)
{
    assert(freq_min <= freq_max);

    // Convert freq_min, freq_max to channel index
    // Note freq_min and freq_max are swapped here since CHIME indices run 800 -> 400
    int chan0 = (int)((800. - freq_max) / 400. * vdif_assembler::nfreq);       // round down 
    int chan1 = (int)((800. - freq_min) / 400. * vdif_assembler::nfreq) + 1;   // round up

    // Clamp to CHIME range
    chan0 = max(chan0, 0);
    chan1 = min(chan1, vdif_assembler::nfreq);

    for (int i = chan0; i < chan1; i++)
	mask[i] = 0;
}


//
// @version = 0 gets the most recent version.
//
// (Currently, only one version of the mask is supported, so the 'version'
// argument is just a hook for future extensions!  However, I plan to propose
// a few RFI mask soon, since the current mask is out-of-date and the RFI
// environment at the site has changed.
//
void get_rfi_mask(int *mask, int version)
{
    for (int i = 0; i < vdif_assembler::nfreq; i++)
	mask[i] = 1;

    if ((version == 0) || (version == 1)) {
	// ch_util November 2015 vintage
	extend_rfi_mask(mask, 439.26, 441.20);
	extend_rfi_mask(mask, 449.81, 450.97);
	extend_rfi_mask(mask, 455.28, 456.05);
	extend_rfi_mask(mask, 469.34, 470.11);
	extend_rfi_mask(mask, 483.40, 484.96);
	extend_rfi_mask(mask, 488.09, 494.33);
	extend_rfi_mask(mask, 499.81, 507.22);
	extend_rfi_mask(mask, 529.89, 536.52);
	extend_rfi_mask(mask, 542.00, 554.49);
	extend_rfi_mask(mask, 565.83, 572.85);
	extend_rfi_mask(mask, 577.93, 584.57);
	extend_rfi_mask(mask, 599.81, 600.58);
	extend_rfi_mask(mask, 695.12, 696.28);
	extend_rfi_mask(mask, 740.04, 745.11);
	extend_rfi_mask(mask, 746.29, 756.05);
	extend_rfi_mask(mask, 799.81, 800.19);
    }
    else
	throw runtime_error("get_rfi_mask: version not recognized");
}


}  // namespace ch_vdif_assembler

/*
 * Local variables:
 *  c-basic-offset: 4
 * End:
 */
