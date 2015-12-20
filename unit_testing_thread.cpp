#include <cstring>
#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// unit_test_buffer


unit_test_buffer::unit_test_buffer(int capacity_)
    : capacity(capacity_)
{
    xassert(capacity > 0);

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_produced, NULL);
    this->ix0 = 0;
    this->ix1 = 0;
    this->buf.resize(capacity);
    this->producer_exit_flag = false;
}


unit_test_buffer::~unit_test_buffer()
{
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&cond_produced);
}


int unit_test_buffer::get_size() const
{
    pthread_mutex_lock(&lock);
    int ret = ix1-ix0;
    pthread_mutex_unlock(&lock);
    return ret;
}


shared_ptr<assembled_chunk> unit_test_buffer::get_chunk()
{
    pthread_mutex_lock(&lock);
    
    for (;;) {
	if (ix0 < ix1) {
	    shared_ptr<assembled_chunk> ret = buf[ix0 % capacity];
	    buf[ix0 % capacity] = shared_ptr<assembled_chunk> ();  // drop reference to allocated memory
	    ix0++;
	    
	    pthread_mutex_unlock(&lock);
	    xassert(ret);
	    return ret;
	}
	
	if (producer_exit_flag) {
	    pthread_mutex_unlock(&lock);
	    return shared_ptr<assembled_chunk> ();
	}

	pthread_cond_wait(&cond_produced, &lock);
    }
}


void unit_test_buffer::put_chunk(const shared_ptr<assembled_chunk> &chunk)
{
    xassert(chunk);
    pthread_mutex_lock(&lock);

    if (producer_exit_flag) {
	pthread_mutex_unlock(&lock);
	throw runtime_error("unit_test_buffer::put_chunk() called after producer_exit()");
    }
    
    if (ix1 >= ix0 + capacity) {
	pthread_mutex_unlock(&lock);
	throw runtime_error("unit_test_buffer capacity exceeded");
    }

    buf[ix1 % capacity] = chunk;
    ix1++;

    pthread_cond_broadcast(&cond_produced);
    pthread_mutex_unlock(&lock);
}


void unit_test_buffer::producer_exit()
{
    pthread_mutex_lock(&lock);
    producer_exit_flag = true;
    pthread_mutex_unlock(&lock);
}


// -------------------------------------------------------------------------------------------------
//
// unit_test_stream


struct unit_test_stream_thread : public thread_base 
{
    shared_ptr<assembler_nerve_center> nc;
    shared_ptr<unit_test_buffer> ubuf;
    int nchunks;
    int asm_nt;

    unit_test_stream_thread(const shared_ptr<assembler_nerve_center> &nc_, const shared_ptr<unit_test_buffer> &ubuf_, int nchunks_, int assembler_nt_)
	: thread_base("unit_test stream thread"),
	  nc(nc_), ubuf(ubuf_), nchunks(nchunks_), asm_nt(assembler_nt_)
    {
	xassert(nc);
	xassert(asm_nt >= 8192);
	xassert(ubuf);
	xassert(nchunks > 0);
    }

    virtual ~unit_test_stream_thread() { }

    // stream_start() gets called prior to spawning thread, but thread is responsible for calling stream_end()
    virtual void thread_body()
    {
	static const int packets_per_chunk = 50000;
	
	// kill assembler if we throw an exception somewhere
	assembler_killer killer(nc, "unit_test simulation thread threw exception");
	
	// starting time (in fpga counts) of the stream
	const uint64_t tstart = asm_nt;

	// max size of random offset (in fpga counts)
	const uint64_t tjitter = (int)(0.9 * asm_nt) - constants::timestamps_per_packet;

	// fpga counts per packet
	const double dt_dpacket = 625./256.;

	// reference assembler state
	uint64_t asm_t0 = tstart;
	shared_ptr<assembled_chunk> asm_chunk0 = make_shared<assembled_chunk> (asm_t0, asm_nt);
	shared_ptr<assembled_chunk> asm_chunk1 = make_shared<assembled_chunk> (asm_t0 + asm_nt, asm_nt);

	for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	    shared_ptr<vdif_chunk> chunk = make_shared<vdif_chunk> (packets_per_chunk, ichunk);
	    chunk->size = chunk->capacity;
	    chunk->set_zero();

	    uint8_t *packet0 = chunk->buf;
	    
	    for (int ipacket = 0; ipacket < packets_per_chunk; ipacket++) {
		uint32_t *header = reinterpret_cast<uint32_t *> (packet0 + ipacket*constants::packet_nbytes);
		uint8_t *data = packet0 + ipacket*constants::packet_nbytes + constants::header_nbytes;
		
		// randomly generate frequency and polarization
		int freq_major = randint(0, 128);
		int pol = randint(0, 2);

		// randomly generate starting time
		uint64_t tpacket = tstart;
		tpacket += (uint64_t)((ichunk*packets_per_chunk+ipacket) * dt_dpacket);
		tpacket += randint(0, tjitter);

		// write vdif header
		int link_id = freq_major / 16;
		int slot_id = freq_major % 16;
		header[3] = (link_id << 16) | (slot_id << 20) | pol;
		header[5] = (uint32_t) tpacket;

		// randomly generate data, and put it in both the vdif_chunk and the assembled_chunks
		for (int i = 0; i < constants::timestamps_per_packet; i++) {
		    // advance assembler
		    if (tpacket + i >= asm_t0 + 2*asm_nt) {
			ubuf->put_chunk(asm_chunk0);
			asm_chunk0 = asm_chunk1;
			asm_chunk1 = make_shared<assembled_chunk> (asm_t0 + 2*asm_nt, asm_nt);
			asm_t0 += asm_nt;
		    }

		    xassert(tpacket + i >= asm_t0);
		    xassert(tpacket + i < asm_t0 + 2*asm_nt);

		    // locate timestamp (tpacket + i) in assembler buffer
		    int s = (2*freq_major + pol) * asm_nt + (tpacket + i - asm_t0);
		    const uint8_t *pp = ((tpacket + i) >= (asm_t0 + asm_nt)) ? (asm_chunk1->buf + s - asm_nt) : (asm_chunk0->buf + s);
		    uint8_t *p = const_cast<uint8_t *> (pp);

		    for (int j = 0; j < 8; j++) {
			uint8_t d = randint(0, 256);
			data[8*i + j] = d;
			p[256*asm_nt*j] = d;
		    }
		}
	    }
	
	    stringstream ss;
	    ss << "unit test: chunk=" << ichunk << "/" << nchunks << "\n";
	    cout << ss.str() << flush;

	    nc->check_alive();
	    nc->stream_put_chunk(chunk, timer);
	}

	ubuf->put_chunk(asm_chunk0);
	ubuf->put_chunk(asm_chunk1);
	ubuf->producer_exit();

	nc->stream_end();
	killer.let_live();
    }
};


struct unit_test_stream : public vdif_stream
{
    shared_ptr<unit_test_buffer> ubuf;
    int nchunks;
    int assembler_nt;


    unit_test_stream(const shared_ptr<unit_test_buffer> &ubuf_, int nchunks_, int assembler_nt_)
	: vdif_stream(false),    // is_realtime=false
	  ubuf(ubuf_), nchunks(nchunks_), assembler_nt(assembler_nt_)
    {
	xassert(ubuf);
	xassert(nchunks > 0);
	xassert(assembler_nt > 0);
    }

    virtual ~unit_test_stream() { }

    virtual void spawn_threads(const shared_ptr<assembler_nerve_center> &nc)
    {
	xassert(nc);
	nc->check_alive();
	spawn_thread<unit_test_stream_thread> (nc, ubuf, nchunks, assembler_nt);
    }
};


shared_ptr<vdif_stream> make_unit_test_stream(const shared_ptr<unit_test_buffer> &ubuf, int nchunks, int assembler_nt)
{
    return make_shared<unit_test_stream> (ubuf, nchunks, assembler_nt);
}


// -------------------------------------------------------------------------------------------------
//
// unit_test_processor


struct unit_test_processor : public vdif_processor
{
    shared_ptr<unit_test_buffer> ubuf;


    unit_test_processor(const shared_ptr<unit_test_buffer> &ubuf_)
	: vdif_processor("unit_test processor", true),    // is_critical=true
	  ubuf(ubuf_)
    {
	xassert(ubuf);
    }

    virtual ~unit_test_processor() { }

    
    virtual void process_chunk(const shared_ptr<assembled_chunk> &a)
    {
	if (a->is_zero())
	    return;

	shared_ptr<assembled_chunk> b;

	for (;;) {
	    b = ubuf->get_chunk();
	    xassert(b);
	    xassert(b->nt == a->nt);
	    
	    if (b->t0 == a->t0) {
		xassert(b->is_equal(*a));
		return;
	    }

	    xassert(b->t0 < a->t0);
	    xassert(b->is_zero());
	}
    }


    virtual void finalize()
    {
	for (;;) {
	    shared_ptr<assembled_chunk> b = ubuf->get_chunk();
	    if (!b)
		return;

	    xassert(b->is_zero());
	}

	cout << "unit test: pass\n" << flush;
    }
};


shared_ptr<vdif_processor> make_unit_test_processor(const shared_ptr<unit_test_buffer> &ubuf)
{
    return make_shared<unit_test_processor> (ubuf);
}


}   // namespace ch_vdif_assembler
