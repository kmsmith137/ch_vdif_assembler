#include <fstream>
#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


struct disk_reader_thread : public thread_base {
    shared_ptr<assembler_nerve_center> nc;
    vector<string> filename_list;
    int thread_id;
    int nfiles_loc;
    int nfiles_tot;

    disk_reader_thread(const shared_ptr<assembler_nerve_center> &nc_, const vector<string> &filename_list_, int thread_id_, int nfiles_tot_)
	: thread_base("disk reader thread " + to_string(thread_id_)),
	  nc(nc_), filename_list(filename_list_), thread_id(thread_id_), nfiles_loc(filename_list_.size()),  nfiles_tot(nfiles_tot_)
    { 
	xassert(thread_id >= 0);
	xassert(thread_id < constants::num_disks);
	xassert(nfiles_tot > 0);
	xassert(nfiles_loc == (nfiles_tot + constants::num_disks - thread_id - 1) / constants::num_disks);
    }

    virtual ~disk_reader_thread() { }

    // stream_start() gets called prior to spawning thread, but thread is responsible for calling stream_end()
    virtual void thread_body()
    {
	// kill assembler if we throw an exception somewhere
	assembler_killer killer(nc, "disk reader thread threw exception");

	for (int ifile = 0; ifile < nfiles_loc; ifile++) {
	    // round-robin file assignment assumed here
	    int seq_id = thread_id + ifile * constants::num_disks;

	    string filename = filename_list[ifile];
	    shared_ptr<vdif_chunk> chunk = make_shared<vdif_chunk> (filename, seq_id);
	    cout << (name + string(": read ") + filename) << endl;

	    nc->stream_put_chunk(chunk, timer);
	}

	int ichunk_last = thread_id + (nfiles_loc-1) * constants::num_disks;
	if (ichunk_last == nfiles_tot-1)
	    nc->stream_end();

	killer.let_live();
    }	
};


struct file_stream : public vdif_stream
{
    vector<string> filename_list;
    int nfiles;

    file_stream(const vector<string> &filename_list_)
	: vdif_stream(false),    // is_realtime = false
	  filename_list(filename_list_),
	  nfiles(filename_list_.size())
    {
	xassert(nfiles > 0);
    }

    virtual ~file_stream() { }


    virtual void spawn_threads(const shared_ptr<assembler_nerve_center> &nc)
    {
	xassert(nc);
	nc->check_alive();

	// spawn threads with round robin file assignment, as assumed above
	for (int it = 0; it < constants::num_disks; it++) {
	    vector<string> filename_sublist;
	    for (int f = it; f < nfiles; f += constants::num_disks)
		filename_sublist.push_back(filename_list[f]);

	    spawn_thread<disk_reader_thread> (nc, filename_sublist, it, filename_list.size());
	}
    }
};


static void parse_file_list(vector<string> &filename_list, const string &list_filename)
{
    filename_list.clear();

    //
    // Very simple parsing which expects one filename per line, with no extra whitespace
    //
    // FIXME could be improved, but low priority
    //
    ifstream f(list_filename.c_str());
    string line;
    
    if (f.fail()) {
	string err = list_filename + ": couldn't open file\n";
	cout << err << flush;
	throw runtime_error(err.c_str());
    }

    while (getline(f, line))
	filename_list.push_back(line);

    if (filename_list.size() == 0) {
	string err = list_filename + ": empty file\n";
	cout << err << flush;
	throw runtime_error(err.c_str());
    }
}


shared_ptr<vdif_stream> make_file_stream(const string &filelist_filename)
{
    vector<string> filename_list;
    parse_file_list(filename_list, filelist_filename);

    return make_file_stream(filename_list);
}


shared_ptr<vdif_stream> make_file_stream(const vector<string> &filename_list)
{
    return make_shared<file_stream> (filename_list);
}


}   // namespace ch_vdif_assembler
