//
// The assembler_nerve_center data structure contains concurrent buffering logic
// involving various threads: stream IO threads, assembler thread, disk writer
// threads, processing threads.
//

#include "ch_vdif_assembler_internals.hpp"

//
// This macro only makes sense in member functions of assembler_nerve_center,
// and must be called with the lock held!
//
#define assert_locked(cond) assert_locked2(cond, __LINE__)

#define assert_locked2(cond,line) \
    do { \
        if (_unlikely(!(cond))) { \
	    this->_kill("Assertion '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")"); \
	    this->_unlock(); \
	    cout << (string(this->killmsg) + "\n") << flush; \
	    throw runtime_error(this->killmsg); \
	} \
    } while (0)


//
// The lock(), unlock(), wait(), and broadcast() macros were my attempt at
// general-purpose concurrency debugging output, but they aren't very useful
// in their current form.  I'll replace them with something better later!
//
// FIXME would it be a useful cleanup to have a lock_releaser object on the stack?
//
// FIXME something else to think about: logic to release buffers "early" if we
// know they'll never be used?
//
#if THREAD_DEBUG < 2

#define lock() this->_lock()
#define unlock() this->_unlock()
#define wait(cond,timer) this->_wait(cond,timer)
#define broadcast(cond) pthread_cond_broadcast(&(cond))

#else

#define lock() lock_debug(__LINE__)
#define lock_debug(line) \
    do { \
	cout << ("lock(): " __FILE__ ":" __STRING(line) "\n") << flush;	\
        this->_lock(); \
    } while(0)

#define unlock() unlock_debug(__LINE__)
#define unlock_debug(line) \
    do { \
	cout << ("unlock(): " __FILE__ ":" __STRING(line) "\n") << flush; \
        this->_unlock(); \
    } while(0)

#define wait(cond,timer) wait_debug(cond,timer,__LINE__)
#define wait_debug(cond,timer,line) \
    do { \
	cout << ("wait(): " __FILE__ ":" __STRING(line) "\n") << flush; \
        this->_wait(cond,timer); \
    } while(0)

#define broadcast(cond) broadcast_debug(cond,__LINE__)
#define broadcast_debug(cond,line) \
    do { \
	cout << ("broadcast(): " __FILE__ ":" __STRING(line) "\n") << flush; \
	pthread_cond_broadcast(&cond); \
    } while(0)

#endif  // THREAD_DEBUG


using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


assembler_nerve_center::assembler_nerve_center(bool write_to_disk, int rbuf_size_, int abuf_size_, int assembler_nt_)
    : assembler_nt(assembler_nt_), is_writing_to_disk(write_to_disk), rbuf_size(rbuf_size_), abuf_size(abuf_size_)
{
    // Reasonable ranges?
    xassert(assembler_nt >= 8192);
    xassert(assembler_nt <= 262144);
    xassert(assembler_nt % 64 == 0);   // cache, simd alignment
    xassert(rbuf_size >= constants::num_disks);
    xassert(rbuf_size <= 80);  // 20 GB, 25 sec, would like to increase but moose would need more memory
    xassert(abuf_size >= 2);
    xassert(abuf_size <= 4);

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond_abuf_produced, NULL);
    pthread_cond_init(&cond_abuf_consumed, NULL);
    pthread_cond_init(&cond_done, NULL);

    for (int ithread = 0; ithread < constants::num_disks; ithread++) {
	pthread_cond_init(&cond_rbuf_produced[ithread], NULL);
	pthread_cond_init(&cond_rbuf_consumed[ithread], NULL);
    }

    startflag = false;
    stream_done = false;
    assembler_done = false;
    processors_done = false;
    disk_writers_done = false;

    killflag = false;
    killmsg = "[should never see this]";
    
    // placeholder values, will be initialized in stream_start()
    is_realtime = false;
    is_writing_to_disk = false;

    ndrops_assembler = 0;
    ndrops_disk_writer = 0;

    rbuf.resize(rbuf_size);
    rbuf_ix = 0;
    rbuf_iasm = 0;
    rbuf_itrigger = 0;

    for (int ithread = 0; ithread < constants::num_disks; ithread++)
	rbuf_idisk[ithread] = ithread;
    
    abuf.resize(abuf_size);
    abuf_ix = 0;

    num_processors = 0;
}


assembler_nerve_center::~assembler_nerve_center()
{
    cout << "assembler object destroyed\n" << flush;

    if (ndrops_assembler > 0) {
	double drop_rate = (double)ndrops_assembler / (double)(rbuf_ix);

	cout << ("  !!! assembler: ndrops=" + to_string(ndrops_assembler) 
		 + ", drop_rate=" + to_string(drop_rate) + " !!!\n") << flush;
    }

    if (ndrops_disk_writer > 0) {
	double drop_rate = (double)ndrops_disk_writer / (double)(rbuf_ix);

	cout << ("  !!! assembler: ndrops=" + to_string(ndrops_disk_writer) 
		 + ", drop_rate=" + to_string(drop_rate) + " !!!\n") << flush;
    }
}

// called with mutex locked
void assembler_nerve_center::_kill(const char *killmsg_)
{
    killflag = true;
    killmsg = killmsg_ ? killmsg_ : "internal error";

    for (int ithread = 0; ithread < constants::num_disks; ithread++) {
	broadcast(cond_rbuf_produced[ithread]);
	broadcast(cond_rbuf_consumed[ithread]);
    }

    broadcast(cond_abuf_produced);
    broadcast(cond_abuf_consumed);
    broadcast(cond_done);
}


void assembler_nerve_center::_test_for_processors_done()
{
    assert_locked(!processors_done);

    if (assembler_done && !num_processors) {
	processors_done = true;

	// Wake up any waiting disk writer threads so they can exit
	for (int ithread = 0; ithread < constants::num_disks; ithread++)
	    broadcast(cond_rbuf_produced[ithread]);
    }
}


void assembler_nerve_center::check_alive()
{
    lock();  // throws exception if killflag is set
    unlock();
}    


void assembler_nerve_center::kill_assembler(const char *killmsg_)
{
    // not assembler_nerve_center::lock()
    pthread_mutex_lock(&mutex);

    if (killflag) {
	// already dead, nothing to do
	pthread_mutex_unlock(&mutex);
	return;
    }

    _kill(killmsg_);
    pthread_mutex_unlock(&mutex);

    cout << (string(killmsg) + ", killing assembler\n") << flush;
    cout << "it may take ~30 secs for all threads to die (this should be fixed soon)\n" << flush;
}


void assembler_nerve_center::wait_until_end()
{
    thread_timer timer_unused;

    lock();

    while (!disk_writers_done)
	wait(cond_done, timer_unused);

    unlock();
}

//
// Note: we currently don't define the counterpart assembler_buf::set_realtime(),
// but if this is ever implemented, it will need to wake up threads that are blocked
// waiting for data!
//
void assembler_nerve_center::set_non_realtime()
{
    lock();
    this->is_realtime = false;
    unlock();
}


void assembler_nerve_center::trigger()
{
    lock();

    int i0 = max(rbuf_itrigger, rbuf_ix - rbuf_size);

    for (int ix = i0; ix < rbuf_ix; ix++) {
	shared_ptr<vdif_chunk> chunk = rbuf[ix % rbuf_size];

	if (chunk->is_on_disk || chunk->want_on_disk)
	    continue;  // no change

	// trigger!
	chunk->want_on_disk = true;
	
	// back up disk_writer thread if necessary
	int ithread = ix % constants::num_disks;
	rbuf_idisk[ithread] = min(rbuf_idisk[ithread], ix);

	broadcast(cond_rbuf_produced[ithread]);
    }

    rbuf_itrigger = rbuf_ix;
    unlock();
}


void assembler_nerve_center::stream_start(bool is_realtime_)
{
    lock();

    if (startflag) {
	unlock();
	throw runtime_error("running the same assembler twice is not currently supported");
    }

    startflag = true;
    is_realtime = is_realtime_;
    bool nothing_to_do = (num_processors == 0) && (!is_writing_to_disk || !is_realtime);

    unlock();

    if (nothing_to_do)
	cout << "assembler: data will be processed to end-of-stream, even though no processors are registered!\n" << flush;
}


void assembler_nerve_center::stream_put_chunk(const std::shared_ptr<vdif_chunk> &chunk, thread_timer &timer)
{
    bool ddrop, adrop;
    int ix = chunk->seq_id;
    int ithread = ix % constants::num_disks;
    int ithread_prev = (ix + constants::num_disks - 1) % constants::num_disks;

    // important, since many routines assume no empty pointers in the buffer
    xassert(chunk);
    
    lock();
    assert_locked(startflag);

    while (rbuf_ix < ix) {
	assert_locked(!stream_done);
	wait(cond_rbuf_produced[ithread_prev], timer);
    }

    for (;;) {
	assert_locked(!stream_done);

	// catches case where seq_id is behind the current buffer location
	assert_locked(rbuf_ix == ix);

	shared_ptr<vdif_chunk> old = rbuf[ix % rbuf_size];
	ddrop = old && !old->is_on_disk && old->want_on_disk;
	adrop = (rbuf_iasm <= ix - rbuf_size);

	if (is_realtime)
	    break;

	if (!adrop && !ddrop)
	    break;

	wait(cond_rbuf_consumed[ithread], timer);
    }

    // Note: want_on_disk set here!
    chunk->want_on_disk = this->is_writing_to_disk;
    
    rbuf[ix % rbuf_size] = chunk;
    rbuf_ix++;
    
    broadcast(cond_rbuf_produced[ithread]);
    unlock();

    if (ddrop) {
	cout << "  !!! disk_writer threads running too slow to keep up with stream, some data was lost !!!\n";
	ndrops_disk_writer++;
    }

    if (adrop) {
	cout << "  !!! assembler running too slow to keep up with stream, some data was lost !!!\n";
	ndrops_assembler++;
    }
}


void assembler_nerve_center::stream_end()
{
    lock();

    assert_locked(startflag);
    assert_locked(!stream_done);

    stream_done = true;

    for (int ithread = 0; ithread < constants::num_disks; ithread++)
	broadcast(cond_rbuf_produced[ithread]);

    unlock();
}


shared_ptr<vdif_chunk> assembler_nerve_center::disk_writer_get_chunk(int ithread, thread_timer &timer)
{
    static const int nd = constants::num_disks;
    xassert(ithread >= 0 && ithread < nd);

    lock();

    for (;;) {
	int ix = rbuf_idisk[ithread];

	// Bookkeeping of drops happens in stream_put_chunk()
	if (ix < rbuf_ix - rbuf_size) {
	    ix = max(rbuf_ix - rbuf_size, 0);
	    ix += (ithread - (ix % nd) + nd) % nd;   // round up to make congruent to ithread (mod nd)
	}

	assert_locked(ix % nd == ithread);

	while (ix < rbuf_ix) {
	    shared_ptr<vdif_chunk> chunk = rbuf[ix % rbuf_size];
	    ix += nd;

	    if (!chunk->is_on_disk && chunk->want_on_disk) {
		rbuf_idisk[ithread] = ix;
		chunk->is_on_disk = true;
		broadcast(cond_rbuf_consumed[ithread]);
		
		unlock();
		return chunk;
	    }
	}

	rbuf_idisk[ithread] = ix;

	if (!processors_done) { 
	    //
	    // It's only safe to return if processor_done is set.  Testing stream_done
	    // isn't sufficient, since a processor might still call trigger().
	    //
	    // When we return from pthread_cond_wait(), the value of rbuf_idisk[ithread]
	    // may have changed due to a trigger.  This loop handles this correctly!
	    //
	    wait(cond_rbuf_produced[ithread], timer);
	    continue;
	}

	//
	// If we get here, this disk_writer is done.
	// First determine whether all the disk writers are done.
	//
	disk_writers_done = true;
	for (int ithread = 0; ithread < constants::num_disks; ithread++) {
	    if (rbuf_idisk[ithread] < rbuf_ix) {
		disk_writers_done = false;
		break;
	    }
	}

	if (disk_writers_done)
	    broadcast(cond_done);

	unlock();
	return shared_ptr<vdif_chunk> ();
    }
}


shared_ptr<vdif_chunk> assembler_nerve_center::assembler_get_chunk(thread_timer &timer)
{
    lock();

    for (;;) {
	assert_locked(!assembler_done);

	// Bookkeeping of drops happens in stream_put_chunk()
	rbuf_iasm = max(rbuf_iasm, rbuf_ix - rbuf_size);

	if (rbuf_iasm < rbuf_ix) {
	    int ithread = rbuf_iasm % constants::num_disks;
	    broadcast(cond_rbuf_consumed[ithread]);

	    shared_ptr<vdif_chunk> ret = rbuf[rbuf_iasm % rbuf_size];
	    rbuf_iasm++;

	    unlock();
	    return ret;
	}

	if (stream_done) {
	    unlock();
	    return shared_ptr<vdif_chunk> ();
	}

	int ithread = rbuf_iasm % constants::num_disks;
	wait(cond_rbuf_produced[ithread], timer);
    }
}


void assembler_nerve_center::assembler_put_chunk(const shared_ptr<assembled_chunk> &chunk, thread_timer &timer)
{
    // important, since many routines assume no empty pointers in the buffer
    xassert(chunk);

    lock();

    for (;;) {
	assert_locked(!assembler_done);
	
	if (is_realtime)
	    break;

	shared_ptr<assembled_chunk> old = abuf[abuf_ix % abuf_size];

	if (!old || (old->pcount >= num_processors))
	    break;

	wait(cond_abuf_consumed, timer);
    }

    // Bookkeeping of drops is done in producer_get_chunk()
    abuf[abuf_ix % abuf_size] = chunk;
    abuf_ix++;

    broadcast(cond_abuf_produced);
    unlock();
}


void assembler_nerve_center::assembler_end()
{
    lock();

    assert_locked(stream_done);
    assert_locked(!assembler_done);
    assembler_done = true;

    _test_for_processors_done();
    broadcast(cond_abuf_produced);
    unlock();
}


void assembler_nerve_center::processor_start()
{
    lock();

    if (processors_done) {
	unlock();
	throw runtime_error("registering processors on a finished assembler is not currently supported");
    }

    num_processors++;
    unlock();
}


shared_ptr<assembled_chunk> assembler_nerve_center::processor_get_chunk(int &ichunk, int &ndrops, thread_timer &timer)
{
    ndrops = 0;

    lock();
    assert_locked(ichunk <= abuf_ix);

    for (;;) {
	int imin = max(abuf_ix - abuf_size, 0);

	if (ichunk < 0)
	    ichunk = imin;
	
	if (ichunk < imin) {
	    ndrops += (imin - ichunk);
	    ichunk = imin;
	}

	if (ichunk < abuf_ix)
	    break;

	if (assembler_done) {
	    unlock();
	    return shared_ptr<assembled_chunk> ();
	}

	wait(cond_abuf_produced, timer);
    }

    shared_ptr<assembled_chunk> ret = abuf[ichunk % abuf_size];
    ichunk++;

    ret->pcount++;
    broadcast(cond_abuf_consumed);
    unlock();

    return ret;
}


void assembler_nerve_center::processor_end(int ichunk)
{
    lock();

    assert_locked(!processors_done);
    assert_locked(ichunk <= abuf_ix);

    int imin = max(abuf_ix - abuf_size, 0);
    
    for (int i = imin; i < ichunk; i++) {
	shared_ptr<assembled_chunk> chunk = abuf[i % abuf_size];
	assert_locked(chunk->pcount > 0);
	chunk->pcount--;
    }

    assert_locked(num_processors > 0);
    num_processors--;

    _test_for_processors_done();

    broadcast(cond_abuf_consumed);
    unlock();
}


}   // namespace ch_vdif_assembler
