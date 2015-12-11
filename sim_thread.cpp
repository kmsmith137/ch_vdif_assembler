#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


struct sim_thread : public thread_base 
{
    shared_ptr<assembler_nerve_center> nc;
    double gbps;
    double nsec;

    sim_thread(const shared_ptr<assembler_nerve_center> &nc_, double gbps_, double nsec_)
	: thread_base("sim thread"),
	  nc(nc_), gbps(gbps_), nsec(nsec_)
    {
	xassert(nc);
	xassert(gbps > 0.0);
	xassert(nsec > 0.0);
    }

    virtual ~sim_thread() { }

    // stream_start() gets called prior to spawning thread, but thread is responsible for calling stream_end()
    virtual void thread_body()
    {
	static const int packets_per_chunk = 50000;
	static const int packet_size = constants::packet_nbytes;
	
	// kill assembler if we throw an exception somewhere
	assembler_killer killer(nc, "sim thread threw exception");

	double packets_per_sec = gbps / (8.0 * packet_size / pow(2.,30.));   // note "giga" is 2^30, not 10^9
	double fpga_counts_per_packet = 1.0 / packets_per_sec / 2.5e-6;
	double seconds_per_chunk = packets_per_chunk / packets_per_sec;
	int nchunks = (int)(nsec / seconds_per_chunk) + 1;
	
	cout << (this->name + string(": packets_per_sec=") + to_string(packets_per_sec)) << endl;
	cout << (this->name + string(": fpga_counts_per_packet=") + to_string(fpga_counts_per_packet)) << endl;
	cout << (this->name + string(": seconds_per_chunk=") + to_string(seconds_per_chunk)) << endl;
	cout << (this->name + string(": nchunks=") + to_string(nchunks)) << endl;
	
	struct timeval tv0 = get_time();

	for (int i = 0; i < nchunks; i++) {
	    double t = time_diff(tv0, get_time());
	    double target = i * seconds_per_chunk;
	    
	    if (t < target) {
		int usec = (int)(1.0e6 * (target-t));
		usec = max(usec,0);  // paranoia
		usleep(usec);
		t = time_diff(tv0, get_time());
	    }

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
	
	    stringstream ss;
	    ss << this->name << ": chunk=" << i << "/" << nchunks << ", t=" << t << ", target=" << target << "\n";
	    cout << ss.str() << flush;

	    nc->stream_put_chunk(chunk, timer);
	}

	nc->stream_end();
	killer.let_live();
    }
};


struct sim_stream : public vdif_stream
{
    double gpbs;
    double nsec;

    sim_stream(double gpbs_, double nsec_)
	: vdif_stream(true),   // is_realtime = true
	  gpbs(gpbs_), nsec(nsec_)
    {
	xassert(gpbs > 0);
	xassert(nsec > 0);
    }

    virtual ~sim_stream() { }

    virtual void spawn_threads(const shared_ptr<assembler_nerve_center> &nc)
    {
	xassert(nc);
	nc->check_alive();
	spawn_thread<sim_thread> (nc, gpbs, nsec);
    }
};


shared_ptr<vdif_stream> make_simulated_stream(double gbps, double nsec)
{
    return make_shared<sim_stream> (gbps, nsec);
}


}   // namespace ch_vdif_assembler
