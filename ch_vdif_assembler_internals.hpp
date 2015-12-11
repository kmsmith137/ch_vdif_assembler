#ifndef _CH_VDIF_ASSEMBLER_INTERNALS_HPP
#define _CH_VDIF_ASSEMBLER_INTERNALS_HPP

#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <iomanip>

#include "ch_vdif_assembler.hpp"

// Set to 0 to run in normal mode
// Set to 1 to see thread spawn/exit events
// Set to 2 to see lock/unlock/wait/broadcast events (verbose!)
#define THREAD_DEBUG 0

namespace ch_vdif_assembler {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// Miscellaneous helper routines


inline double time_diff(const struct timeval &tv1, const struct timeval &tv2)
{
    return (tv2.tv_sec - tv1.tv_sec) + 1.0e-6 * (tv2.tv_usec - tv1.tv_usec);
}

inline struct timeval get_time()
{
    struct timeval ret;
    if (gettimeofday(&ret, NULL) < 0)
	throw std::runtime_error("gettimeofday() failed");
    return ret;
}

// FIXME switch to thread-safe rng (currently doesn't matter since only one thread ever generates random numbers)
inline double uniform_rand()
{
    return (rand() + 0.5) / (RAND_MAX + 1.0);
}

inline int randint(int lo, int hi)
{
    int ret = lo + (int)((hi-lo)*uniform_rand());
    ret = std::max(ret, lo);  // should be redundant
    ret = std::min(ret, hi);  // should be redundant
    return ret;
}

// sometimes useful for multithreaded debugging
inline void deadlock()
{
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_lock(&mutex);

    pthread_cond_t cond;
    pthread_cond_init(&cond, NULL);
    pthread_cond_wait(&cond, &mutex);
}


extern std::string make_dataset_name();
extern std::string make_data_dir(const std::string &dataset_name, int disk_id);


// -------------------------------------------------------------------------------------------------
//
// Thread utilities


extern void spawn_assembler_thread(const std::shared_ptr<assembler_nerve_center> &nc);
extern void spawn_disk_writer_thread(const std::shared_ptr<assembler_nerve_center> &nc, const std::string &data_dir, int ithread);
extern void spawn_processing_thread(const std::shared_ptr<assembler_nerve_center> &nc, const std::shared_ptr<vdif_processor> &p);


struct thread_timer {
    double total_waiting_time;
    double total_running_time;

    struct timeval tv0_run;
    struct timeval tv0_wait;

    thread_timer() : total_waiting_time(0), total_running_time(0) { }

    void start_running()  { tv0_run = get_time(); }
    void stop_running()   { total_running_time += time_diff(tv0_run, get_time()); }
    void start_waiting()  { tv0_wait = get_time(); }
    void stop_waiting()   { total_waiting_time += time_diff(tv0_wait, get_time()); }

    double busyfrac()     { return 1.0 - total_waiting_time / total_running_time; }
};


struct thread_base {
    std::string name;
    pthread_t pthread;

    // Not protected by a mutex; only accessed by thread itself
    thread_timer timer;

    thread_base(const std::string &name);
    virtual ~thread_base() { }

    virtual void thread_body() = 0;

    //
    // Assumes thread_base has been allocated with 'new'
    // If the thread fails to spawn, calls 'delete' and throws an exception.
    // Thus the caller can assume that the pointer is always deleted.
    //
    static void _spawn(thread_base *p);

    static void *_pthread_main(void *arg);
};


//
// A one-liner for spawning threads:
//   spawn_thread<thread_subclass> (arg1, arg2, ...);
//
template<typename T, typename... Args>
void spawn_thread(Args... args)
{
    T *t = new T(args...);
    thread_base::_spawn(t);
}


// -------------------------------------------------------------------------------------------------
//
// Some core classes


struct vdif_chunk : noncopyable {
    uint8_t *buf;
    int capacity;  // in packets
    int size;      // in packets

    int seq_id;
    bool is_on_disk;
    bool want_on_disk;

    // Warning: this constructor allocates the buffer but doesn't initialize it
    // In both constructors, it's important that the seq_ids go 0,1,2,etc.  (See stream_put_chunk below)
    // This constructor sets is_on_disk=false
    vdif_chunk(int npackets, int seq_id);

    // Read from file (this constructor sets is_on_disk=true)
    vdif_chunk(const std::string &filename, int seq_id);

    void set_zero();
    void write(const std::string &filename);

    ~vdif_chunk();
};


struct assembler_killer {
    std::shared_ptr<assembler_nerve_center> nc;
    const char *killmsg;

    assembler_killer();
    assembler_killer(const std::shared_ptr<assembler_nerve_center> &nc_, const char *killmsg_);
    ~assembler_killer();

    void set_victim(const std::shared_ptr<assembler_nerve_center> &nc_, const char *killmsg_);
    void let_live();
};


//
// The 32-bit FPGA counts used as timestamps "wrap around" every 3 hours.
//
// This helper class converts a stream of unsigned 32-bit FPGA counts to
// a stream of signed 64-bit timestamps without wraparound.  It correctly
// handles out-of-order timestamps as long as the "jitter" isn't more than
// 2^31 timestamps (which is an unrealistic case).
//
class timestamp_unwrapper {
private:
    int64_t last_timestamp;

public:
    timestamp_unwrapper()
    {
	// this initial value will mean that the first timestamp is in the range [0,2^32-1]
	last_timestamp = int64_t(1) << 31;
    }

    inline int64_t unwrap(uint32_t x)
    {
	int32_t delta = (int32_t)x - (int32_t)last_timestamp;
	last_timestamp += (int64_t)delta;
	return last_timestamp;
    }
};


// -------------------------------------------------------------------------------------------------
//
// The monster: assembler_nerve_center


class assembler_nerve_center : noncopyable
{
protected:
    //
    // Might fine-grain later
    // Should broadcast (not signal) all these conditionals
    //
    pthread_mutex_t mutex;
    pthread_cond_t cond_rbuf_produced[constants::num_disks];
    pthread_cond_t cond_rbuf_consumed[constants::num_disks];
    pthread_cond_t cond_abuf_produced;
    pthread_cond_t cond_abuf_consumed;
    pthread_cond_t cond_done;

    //
    // Assembler state
    //
    // Right now, the assembler can only be run on one stream!
    // However, the data structures are written with asychronous
    // operation in mind, so this will be easy (I hope!) to generalize 
    // later once we understand the typical use cases.
    //
    bool startflag;
    bool stream_done;        // stream_end() called
    bool assembler_done;     // assembler_end() called
    bool processors_done;    // assembler_end() called and last processor detached
    bool disk_writers_done;
    
    // Set if disaster strikes (e.g. exception in assembler thread)
    bool killflag;
    const char *killmsg;

    int assembler_nt;
    bool is_writing_to_disk;
    bool is_realtime;

    int ndrops_assembler;
    int ndrops_disk_writer;
    
    // ring buffer containing unasembled chunks
    std::vector<std::shared_ptr<vdif_chunk> > rbuf;
    int rbuf_size;
    int rbuf_ix;
    int rbuf_iasm;                           // index of first unassembled chunk
    int rbuf_itrigger;                       // index of first chunk which is a candidate for trigger()
    int rbuf_idisk[constants::num_disks];    // index of first chunk which is a candidate for writing to disk

    // ring buffer containing assembled chunks
    std::vector<std::shared_ptr<assembled_chunk> > abuf;
    int abuf_size;
    int abuf_ix;

    int num_processors;


    inline void _lock()
    {
	pthread_mutex_lock(&mutex);

	if (killflag) {
	    pthread_mutex_unlock(&mutex);
	    throw std::runtime_error(killmsg);
	}
    }

    // Member functions beginning with underscore (except _lock()) are called with the lock held
    inline void _unlock()
    {
	pthread_mutex_unlock(&mutex);
    }

    inline void _wait(pthread_cond_t &cond, thread_timer &timer) 
    {
	timer.start_waiting();
	pthread_cond_wait(&cond, &mutex);
	timer.stop_waiting();

	if (killflag) {
	    pthread_mutex_unlock(&mutex);
	    throw std::runtime_error(killmsg);
	}
    }

    void _kill(const char *killmsg);
    void _test_for_processors_done();

public:
    assembler_nerve_center(bool write_to_disk, int rbuf_size, int abuf_size, int assembler_nt);
    ~assembler_nerve_center();

    inline int get_assembler_nt() const { return assembler_nt; }

    // High-level control
    void check_alive();
    void kill_assembler(const char *killmsg);
    void wait_until_end();
    void set_non_realtime();
    void trigger();

    //
    // Called by stream I/O threads.
    //
    // In stream_put_chunk(), chunk->seq_id is used to determine placement in the buffer.
    // If this is ahead of the current buffer location, the call will block, expecting another thread
    // to place the intermediates.  If behind the current buffer location, the assembler will be killed.
    //
    void stream_start(bool is_realtime);
    void stream_put_chunk(const std::shared_ptr<vdif_chunk> &chunk, thread_timer &timer);
    void stream_end();

    // Called by disk writer threads
    std::shared_ptr<vdif_chunk> disk_writer_get_chunk(int ithread, thread_timer &timer);

    // Called by assembler threads
    std::shared_ptr<vdif_chunk> assembler_get_chunk(thread_timer &timer);
    void assembler_put_chunk(const std::shared_ptr<assembled_chunk> &chunk, thread_timer &timer);
    void assembler_end();

    // Called by processor threads
    void processor_start();
    std::shared_ptr<assembled_chunk> processor_get_chunk(int &ichunk, int &ndrops, thread_timer &timer);
    void processor_end(int ichunk);
};


// -------------------------------------------------------------------------------------------------
//
// A lowish-level interface to the assembler's output
//
// C++ processors will probably want to use the higher-level interface in class 'vdif_processor'.
//
// (Python processors use the lower-level interface, but this is hidden in 
//  ch_vdif_assembler.vdif_processor.run())
//


struct processor_handle : noncopyable {
    std::string name;
    std::shared_ptr<assembler_nerve_center> nc;

    int ichunk;  // next chunk to retrieve
    int ndrops;
    int nprocessed;

    // constructor registers processor
    processor_handle(const std::string &name, const std::shared_ptr<assembler_nerve_center> &nc);

    // destructor unregisters processor
    ~processor_handle();

    // Called by processor threads
    std::shared_ptr<assembled_chunk> get_next_chunk(thread_timer &timer);
};


}  // namespace ch_vdif_assembler


#endif // _CH_VDIF_ASSEMBLER_INTERNALS_HPP
