#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


// p->name, with no possibility of segfault
inline string get_name(const shared_ptr<vdif_processor> &p)
{
    xassert(p);
    return p->name;
}

static void throw_rerun_exception()
{
    static const char *msg = "currently, a vdif_processor can't be run twice, you'll need to re-construct and start from scratch\n";
    cout << msg << flush;
    throw runtime_error(msg);
}


// -------------------------------------------------------------------------------------------------


struct processing_thread : public thread_base {
    shared_ptr<assembler_nerve_center> nc;
    shared_ptr<vdif_processor> processor;


    processing_thread(const shared_ptr<assembler_nerve_center> &nc_, const shared_ptr<vdif_processor> &processor_)
	: thread_base(get_name(processor_) + " thread"),
	  nc(nc_), processor(processor_)
    {
	xassert(nc);
	xassert(processor);
	
	if (processor->is_running())
	    throw_rerun_exception();
    }

    virtual ~processing_thread() { }
	
    void thread_body()
    {
	processor->set_running();

	processor_handle ph(processor->name, nc);

	assembler_killer killer;
	if (processor->is_critical)
	    killer.set_victim(nc, "exception thrown in critical processor");

	for (;;) {
	    shared_ptr<assembled_chunk> chunk = ph.get_next_chunk(timer);
	    
	    if (!chunk) {
		processor->finalize();
		break;
	    }

	    processor->process_chunk(chunk);
	}

	killer.let_live();
    }
};


// -------------------------------------------------------------------------------------------------


vdif_processor::vdif_processor(const string &name_, bool is_critical_)
    : name(name_), is_critical(is_critical_), runflag(false)
{
    pthread_mutex_init(&mutex, NULL);
}


bool vdif_processor::is_running()
{
    pthread_mutex_lock(&mutex);
    bool ret = runflag;
    pthread_mutex_unlock(&mutex);
    
    return ret;
}


void vdif_processor::set_running()
{
    pthread_mutex_lock(&mutex);
    
    if (runflag) {
	pthread_mutex_unlock(&mutex);
	throw_rerun_exception();
    }

    runflag = true;
    pthread_mutex_unlock(&mutex);
}


// -------------------------------------------------------------------------------------------------


void spawn_processing_thread(const shared_ptr<assembler_nerve_center> &nc, const shared_ptr<vdif_processor> &p)
{
    xassert(nc);
    nc->check_alive();
    spawn_thread<processing_thread> (nc, p);
}


}   // namespace ch_vdif_assembler
