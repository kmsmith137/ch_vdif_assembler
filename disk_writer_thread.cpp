#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


struct disk_writer_thread : public thread_base {
    shared_ptr<assembler_nerve_center> nc;
    string outdir;
    int thread_id;
    bool mkdir_called;

    disk_writer_thread(const shared_ptr<assembler_nerve_center> &nc_, const string &outdir_, int thread_id_)
	: thread_base("disk writer thread " + to_string(thread_id_)),
	  nc(nc_), outdir(outdir_), thread_id(thread_id_), mkdir_called(false)
    { 
	xassert(thread_id >= 0);
	xassert(thread_id < constants::num_disks);
    }

    virtual ~disk_writer_thread() { }


    // Devirtualize thread_base
    virtual void thread_body()
    {
	// Kill the assembler if we throw an exception before the end
	assembler_killer k(nc, "disk_writer thread threw exception");
    
	for (;;) {
	    shared_ptr<vdif_chunk> data = nc->disk_writer_get_chunk(thread_id, timer);
	    if (!data)
		break;

	    if (!mkdir_called) {
		xmkdir(outdir);
		cout << (name + string(": created directory ") + outdir + "\n") << flush;
		mkdir_called = true;
	    }
	    
	    // Make filename
	    stringstream ss;
	    ss << "/" << setfill('0') << setw(7) << data->seq_id << ".dat";
	    
	    string basename = ss.str();
	    string filename = outdir + basename;
	    
	    data->write(filename);
	    cout << (name + string(": wrote ") + filename + "\n") << flush;
	}
	
	k.let_live();
    }	
};


void spawn_disk_writer_thread(const shared_ptr<assembler_nerve_center> &nc, const string &data_dir, int ithread)
{
    xassert(nc);
    nc->check_alive();
    spawn_thread<disk_writer_thread> (nc, data_dir, ithread);
}


}   // namespace ch_vdif_assembler
