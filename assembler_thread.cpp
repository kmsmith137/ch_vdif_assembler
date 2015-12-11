#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


//
// A helper class for the assembler thread
//
struct _assembler_fragment {
    int s;             // array stride in output
    int n;             // number of timestamps in fragment
    uint8_t *dst;
    const uint8_t *src;

    _assembler_fragment(int assembler_nt) 
	: s(constants::chime_nfreq/4 * assembler_nt),
	  n(0), dst(NULL), src(NULL) 
    { }

    inline void prefetch()
    {
	if (n > 0) {
	    // Prefetch to L1 cache (_MM_HINT_T0)
	    _mm_prefetch(dst, _MM_HINT_T0);
	    _mm_prefetch(dst+s, _MM_HINT_T0);
	    _mm_prefetch(dst+2*s, _MM_HINT_T0);
	    _mm_prefetch(dst+3*s, _MM_HINT_T0);
	    _mm_prefetch(dst+4*s, _MM_HINT_T0);
	    _mm_prefetch(dst+5*s, _MM_HINT_T0);
	    _mm_prefetch(dst+6*s, _MM_HINT_T0);
	    _mm_prefetch(dst+7*s, _MM_HINT_T0);
	}
    }

    inline void assemble()
    {
	// ch_vdif_assembler_kernels.hpp
	_assemble128(dst, s, src, n);
    }
};


// -------------------------------------------------------------------------------------------------


struct assembler_thread : public thread_base {
    shared_ptr<assembler_nerve_center> nc;
    const int assembler_nt;

    shared_ptr<assembled_chunk> chunk0;
    shared_ptr<assembled_chunk> chunk1;

    //
    // The prefetching logic below requires that chunks wait here for one iteration
    // of the packet procesing loop, before being passed to the processing threads.
    //
    shared_ptr<assembled_chunk> inactive_chunk0;
    shared_ptr<assembled_chunk> inactive_chunk1;

    int64_t chunk_t1;             // always equal to chunk1->t0
    uint8_t *buf0;                // always equal to chunk0->buf
    uint8_t *buf1;                // always equal to chunk1->buf
    bool have_inactive_chunks;    // always equal to (inactive_chunk0 || inactive_chunk1)
    bool dropflag;                // just used to throttle a warning message

    timestamp_unwrapper ts_unwrapper;


    assembler_thread(const shared_ptr<assembler_nerve_center> &nc_)
	: thread_base("assembler thread"),
	  nc(nc_), assembler_nt(nc_->get_assembler_nt()),
	  chunk_t1(-(1<<30)),   // see below for explanation of this sentinel value
	  buf0(NULL), buf1(NULL), have_inactive_chunks(false), dropflag(false)
    {
	xassert(assembler_nt % 64 == 0);
	xassert(assembler_nt >= 8192);     // using less than 8K is probably a bad idea!
    }

    virtual ~assembler_thread() { }

    inline void analyze_packet(_assembler_fragment &frag0, _assembler_fragment &frag1, const uint8_t *packet)
    {
	static const int packet_nt = constants::timestamps_per_packet;

	const uint32_t *header = reinterpret_cast<const uint32_t *> (packet);    
	uint32_t pol = (header[3] & 0x0000ffff);
	uint32_t link_id = (header[3] & 0x000f0000) >> 16;
	uint32_t slot_id = (header[3] & 0x01f00000) >> 20;
	uint32_t fpga_count = header[5];
	
	//
	// Fast test for bad packet using bitwise operators.
	//
	// Pathfinder-specific: test is ((link_id < 8) && (slot_id < 16) && (pol < 2)).
	//
	// We should never get bad packets, but processing them "blind" would cause memory
	// corruption, so I included this test even though it's in a speed-critical part
	// of the code.
	//
	uint32_t bad_packet = (header[3] & 0x0108fffe);
	if (_unlikely(bad_packet))
	    throw runtime_error("vdif_assembler: bad packet");
	
	int64_t packet_t0 = ts_unwrapper.unwrap(fpga_count);
	int64_t packet_t1 = packet_t0 + packet_nt;    
	
	if (_unlikely(packet_t0 < chunk_t1 - assembler_nt)) {
	    // drop packet
	    if (!dropflag) {
		cout << "warning: assembler dropping packets due to large out-of-order arrival times\n" << flush;
		dropflag = true;
	    }

	    frag0.n = 0;
	    frag1.n = 0;
	    return;
	}

	if (_unlikely(packet_t1 > chunk_t1 + assembler_nt))
	    this->advance_buffers(packet_t1);

	int freq_major = slot_id + 16*link_id;

	int dt = packet_t0 - chunk_t1;                       // packet arrival time relative to buf1
	int i0 = (2*freq_major + pol) * assembler_nt + dt;   // array offset of beginning of packet in buf1
	int m = min(max(-dt,0), packet_nt);                  // overlap between packet and buf0
	
	frag0.n = m;
	frag0.dst = buf0 + (i0 + assembler_nt);
	frag0.src = packet + constants::header_nbytes;

	frag1.n = packet_nt - m;
	frag1.dst = buf1 + (i0 + m);
	frag1.src = packet + constants::header_nbytes + 8*m;
    }


    void thread_body()
    {
	// Kill assembler if we throw an exception somewhere
	assembler_killer killer(nc, "assembler thread threw exception");
	
	_assembler_fragment frag00(assembler_nt);
	_assembler_fragment frag01(assembler_nt);
	_assembler_fragment frag10(assembler_nt);
	_assembler_fragment frag11(assembler_nt);
	
	//
	// Initialize chunk_t1 to a huge negative sentinel value.  This
	// ensures that the first call to analyze_packet() allocates buffers.
	//
	this->chunk_t1 = -(1 << 30);

	for (;;) {
	    shared_ptr<vdif_chunk> packets = nc->assembler_get_chunk(timer);

	    if (!packets)
		break;  // end of stream
	
	    int npackets = packets->size;
	    
	    if (npackets == 0)
		continue;
	    
	    // Analyze first packet
	    this->dropflag = false;
	    uint8_t *packet0 = packets->buf;
	    analyze_packet(frag00, frag01, packet0);
	    
	    if (_unlikely(have_inactive_chunks))
		this->process_inactive_chunks();
	    
	    // Note npackets-1 here
	    for (int ipacket = 0; ipacket < npackets-1; ipacket++) {
		// Invariant: at top of loop, packet i has been analyzed in {frag00,frag01}
		
		// Analyze and prefetch packet (i+1)
		analyze_packet(frag10, frag11, packet0 + (ipacket+1) * constants::packet_nbytes);
		frag10.prefetch();
		frag11.prefetch();
		
		// Assemble packet[i]
		frag00.assemble();
		frag01.assemble();
		
		if (_unlikely(have_inactive_chunks))
		    this->process_inactive_chunks();
		
		// Preserve loop invariant
		frag00 = frag10;
		frag01 = frag11;
	    }
	    
	    // Assemble last packet (already analyzed)
	    frag00.assemble();
	    frag01.assemble();
	    
	    if (_unlikely(have_inactive_chunks))
		this->process_inactive_chunks();
	}

	// FIXME revisit if we introduce stream pausing/unpausing
	nc->set_non_realtime();

	this->process_all_chunks();

	nc->assembler_end();
	killer.let_live();
    }


    void advance_buffers(int64_t packet_t1)
    {
	xassert(!inactive_chunk0);
	xassert(!inactive_chunk1);
	xassert(!have_inactive_chunks);
	xassert(packet_t1 > chunk_t1 + assembler_nt);
	
	if (packet_t1 <= chunk_t1 + 2*assembler_nt) {
	    chunk_t1 += assembler_nt;
	    
	    inactive_chunk0 = chunk0;
	    chunk0 = chunk1;
	    chunk1 = make_shared<assembled_chunk> (chunk_t1, assembler_nt);
	}
	else {
	    xassert(packet_t1 >= 0);  // FIXME not guaranteed but extremely likely
	    chunk_t1 = packet_t1 - (packet_t1 % assembler_nt);
	    
	    inactive_chunk0 = chunk0;
	    inactive_chunk1 = chunk1;
	    chunk0 = make_shared<assembled_chunk> (chunk_t1-assembler_nt, assembler_nt);
	    chunk1 = make_shared<assembled_chunk> (chunk_t1, assembler_nt);
	}

	xassert(packet_t1 >= chunk_t1);
	xassert(packet_t1 <= chunk_t1 + assembler_nt);
	
	// const_cast: the assembler is allowed to modify this buffer but no one else
	this->buf0 = const_cast<uint8_t *> (chunk0->buf);
	this->buf1 = const_cast<uint8_t *> (chunk1->buf);
	this->have_inactive_chunks = (!!inactive_chunk0) || (!!inactive_chunk1);
    }


    void process_inactive_chunks()
    {
	if (inactive_chunk0) {
	    nc->assembler_put_chunk(inactive_chunk0, timer);
	    inactive_chunk0 = shared_ptr<assembled_chunk> ();
	}

	if (inactive_chunk1) {
	    nc->assembler_put_chunk(inactive_chunk1, timer);
	    inactive_chunk1 = shared_ptr<assembled_chunk> ();
	}

	have_inactive_chunks = false;
    }


    void process_all_chunks()
    {
	this->process_inactive_chunks();

	if (chunk0) {
	    nc->assembler_put_chunk(chunk0, timer);
	    chunk0 = shared_ptr<assembled_chunk> ();
	}

	if (chunk1) {
	    nc->assembler_put_chunk(chunk1, timer);
	    chunk1 = shared_ptr<assembled_chunk> ();
	}

	chunk_t1 = -(1 << 30);
	buf0 = NULL;
	buf1 = NULL;
    }
};


void spawn_assembler_thread(const std::shared_ptr<assembler_nerve_center> &nc)
{
    xassert(nc);
    nc->check_alive();
    spawn_thread<assembler_thread> (nc);
}


}  // namespace ch_vdif_assembler
