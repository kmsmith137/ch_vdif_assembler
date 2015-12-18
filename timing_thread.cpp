#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


struct timing_thread : public thread_base 
{
    shared_ptr<assembler_nerve_center> nc;
    int packets_per_chunk;
    int nchunks;

    timing_thread(const shared_ptr<assembler_nerve_center> &nc_, int packets_per_chunk_, int nchunks_)
	: thread_base("timing thread"),
	  nc(nc_), packets_per_chunk(packets_per_chunk_), nchunks(nchunks_)
    {
	xassert(packets_per_chunk >= 8192);
	xassert(nchunks >= 4);
    }

    virtual ~timing_thread() { }

    // stream_start() gets called prior to spawning thread, but thread is responsible for calling stream_end()
    virtual void thread_body()
    {
	static const int packet_size = constants::packet_nbytes;

	double assumed_gbps = 6.4;
	double packets_per_sec = assumed_gbps / (8.0 * packet_size / pow(2.,30.));
	double fpga_counts_per_packet = 1.0 / packets_per_sec / 2.5e-6;
	
	// kill assembler if we throw an exception somewhere
	assembler_killer killer(nc, "timing thread threw exception");
	
	// FIXME cut-and-paste with sim_thread here; could define helper function or common base class
	for (int i = 0; i < nchunks; i++) {
	    struct timeval tv0 = get_time();

	    shared_ptr<vdif_chunk> chunk = make_shared<vdif_chunk> (packets_per_chunk, i);
	    chunk->size = chunk->capacity;
	    chunk->set_zero();
	    
	    uint8_t *packet0 = chunk->buf;

	    for (int j = 0; j < packets_per_chunk; j++) {
		uint8_t *packet = packet0 + j*packet_size;
		uint32_t *header = reinterpret_cast<uint32_t *> (packet);
		
		int freq_major = randint(0, 128);
		int pol = randint(0, 2);
		
		int link_id = freq_major / 16;
		int slot_id = freq_major % 16;
		header[3] = (link_id << 16) | (slot_id << 20) | pol;
		
		uint64_t t0 = (uint64_t) randint(100000, 150000);
		t0 += (uint64_t)((i*packets_per_chunk+j) * fpga_counts_per_packet);
		header[5] = (uint32_t) t0;
	    }

	    double dt = time_diff(tv0, get_time());
	    double instantaneous_gbps = assumed_gbps * packets_per_chunk / packets_per_sec / dt;
	
	    stringstream ss;
	    ss << this->name << ": chunk=" << i << "/" << nchunks << ", instantaneous_gpbs=" << instantaneous_gbps << "\n";
	    cout << ss.str() << flush;

	    nc->stream_put_chunk(chunk, timer);
	}

	nc->stream_end();
	killer.let_live();
    }
};


struct timing_stream : public vdif_stream
{
    int npackets_per_chunk;
    int nchunks;

    timing_stream(int npackets_per_chunk_, int nchunks_)
	: vdif_stream(false),   // is_realtime = false
	  npackets_per_chunk(npackets_per_chunk_), nchunks(nchunks_)
    {
	xassert(npackets_per_chunk >= 8192);
	xassert(nchunks >= 4);
    }

    virtual ~timing_stream() { }

    virtual void spawn_threads(const shared_ptr<assembler_nerve_center> &nc)
    {
	xassert(nc);
	nc->check_alive();
	spawn_thread<timing_thread> (nc, npackets_per_chunk, nchunks);
    }
};
    

shared_ptr<vdif_stream> make_timing_stream(int npackets_per_chunk, int nchunks)
{
    return make_shared<timing_stream> (npackets_per_chunk, nchunks);
}


}   // namespace ch_vdif_assembler
