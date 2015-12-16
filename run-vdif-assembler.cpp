#include <cstring>
#include "ch_vdif_assembler_internals.hpp"

using namespace std;
using namespace ch_vdif_assembler;


static void usage()
{
    cerr << "Usage: run-vdif-assembler [FLAGS]\n"
	 << "where FLAGS include the following\n"
	 << "    -w waterfall_plot_outdir           to make waterfall plots\n"
	 << "    -r rfi_histogram_outfile.hdf5      to make rfi histograms\n"
	 << "    -R rfi_histogram_outfile.hdf5      to make rfi histograms with a reference implementation\n"
	 << "    -f file_list.txt                   to run on a disk capture\n"
	 << "    -n                                 to run on a real-time network capture\n"
	 << "    -s                                 to run on a simulated network capture (6.4 Gbps, 60 sec)\n"
	 << "    -S num_seconds                     to run on a simulated network capture (6.4 Gpbs, specified duration)\n"
	 << "    -t                                 to run in \"timing mode\": reported running time will be determined by the slowest thread\n"
	 << "    -d                                 to save stream on disk (will no-op if already running on a disk capture)\n"
	 << "    -m                                 to run a concurrent \"mischief thread\" which memcpy's between two 0.5 GB buffers\n"
	 << "\n"
	 << "You may find the script show-moose-acqusitions.py useful for making file lists\n"
	 << "\n"
	 << "Suggested usage:\n"
	 << "    show-moose-acquisitions.py 41537 > filelist_41537.txt\n"
	 << "    run-vdif-assembler -w waterfall_41537 -f filelist_41537.txt\n"
	 << "    index-vdif-waterfalls.py waterfall_41537\n"
	 << "\n"
	 << "The index-vdif-waterfalls.py script makes an HTML summary page with clickable thumbnails\n";

    exit(2);
}


// -------------------------------------------------------------------------------------------------
//
// mischief thread


// these globals are initialized in main(), before spawning the thread
static pthread_t mischief_thread;
static pthread_mutex_t mischief_mutex;
static bool mischief_killflag;


static void *mischief_pthread_main(void *arg)
{
    static const int nbytes = (1 << 29);  // 0.5 GB
    
    char *buf1 = new char[nbytes];
    char *buf2 = new char[nbytes];

    memset(buf1, 0, nbytes);
    memset(buf2, 0, nbytes);
    
    cout << "mischief thread running\n" << flush;
    
    // caller initializes mutex, conditional
    for (;;) {
	pthread_mutex_lock(&mischief_mutex);

	if (mischief_killflag) {
	    pthread_mutex_unlock(&mischief_mutex);
	    cout << "mischief thread done\n" << flush;
	    return NULL;
	}

	pthread_mutex_unlock(&mischief_mutex);
	memcpy(buf1, buf2, nbytes);
	memcpy(buf2, buf1, nbytes);
    }
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    // only used in timing thread (-t)
    static const int timing_npackets_per_chunk = 50000;
    static const int timing_nchunks = 128;

    bool is_timing = false;
    bool write_to_disk = false;
    bool mischief_flag = false;
    shared_ptr<vdif_stream> stream;
    shared_ptr<vdif_processor> waterfall_plotter;
    shared_ptr<vdif_processor> rfi_histogrammer;
    vector<shared_ptr<vdif_processor> > processors;

    //
    // Low-tech command line parsing.
    //
    // Note that when constructing processors, we always set is_critical=true.
    // The assumption is that if we're running this standalone driver program, we
    // just want to crash if a processing thread throws an exception.  In the 
    // realtime system, we may want to do someting more sophisticated.
    //
    int pos = 1;
    while (pos < argc) {
	if ((strlen(argv[pos]) != 2) || (argv[pos][0] != '-'))
	    usage();

	char cs = argv[pos][1];

	// Parse switches with no positional argument

	if (cs == 'n') {
	    if (stream)
		usage();
	    stream = make_network_stream();
	    pos++;
	    continue;
	}

	if (cs == 'd') {
	    write_to_disk = true;
	    pos++;
	    continue;
	}

	if (cs == 's') {
	    if (stream)
		usage();
	    stream = make_simulated_stream();
	    pos++;
	    continue;
	}

	if (cs == 't') {
	    if (stream)
		usage();
	    stream = make_timing_stream(timing_npackets_per_chunk, timing_nchunks);
	    is_timing = true;
	    pos++;
	    continue;
	}

	if (cs == 'm') {
	    if (mischief_flag)
		usage();
	    mischief_flag = true;
	    pos++;
	    continue;
	}

	// Parse switches with one positional argument

	if ((pos+1 >= argc) || (argv[pos+1][0] == '-'))
	    usage();

	const char *positional_arg = argv[pos+1];

	if (cs == 'f') {
	    if (stream)
		usage();
	    stream = make_file_stream(positional_arg);
	    pos += 2;
	    continue;
	}

	if (cs == 'S') {
	    if (stream)
		usage();
	    double nsec = atof(positional_arg);
	    cout << "simulated stream duration: " << nsec << " seconds\n";
	    xassert(nsec > 0.0);
	    stream = make_simulated_stream(6.4, nsec);
	    pos += 2;
	    continue;
	}

	if (cs == 'w') {
	    if (waterfall_plotter)
		usage();
	    waterfall_plotter = make_waterfall_plotter(positional_arg, true);
	    processors.push_back(waterfall_plotter);
	    pos += 2;
	    continue;
	}

	if (cs == 'r') {
	    if (rfi_histogrammer)
		usage();
	    rfi_histogrammer = make_rfi_histogrammer(positional_arg, true, false);
	    processors.push_back(rfi_histogrammer);
	    pos += 2;
	    continue;
	}

	if (cs == 'R') {
	    if (rfi_histogrammer)
		usage();
	    rfi_histogrammer = make_rfi_histogrammer(positional_arg, true, true);
	    processors.push_back(rfi_histogrammer);
	    pos += 2;
	    continue;
	}

	// unrecognized switch
	usage();
    }

    if (!stream)
	usage();
    if (!is_timing && !write_to_disk && !processors.size())
	usage();    // nothing to do
    
    vdif_assembler assembler(write_to_disk);

    for (unsigned int i = 0; i < processors.size(); i++)
	assembler.register_processor(processors[i]);

    if (mischief_flag) {
	pthread_mutex_init(&mischief_mutex, NULL);
	mischief_killflag = false;
	
	if (pthread_create(&mischief_thread, NULL, mischief_pthread_main, NULL) < 0) {
	    cout << (string("couldn't spawn mischief thread: ") + strerror(errno) + "\n") << flush;
	    exit(1);
	}
    }

    struct timeval tv0 = get_time();

    assembler.run(stream);

    if (is_timing) {
	// note: "giga" is 2^30, not 10^9
	double dt = time_diff(tv0, get_time());
	double gbytes = (double)timing_nchunks * (double)timing_npackets_per_chunk * (double)constants::packet_nbytes / pow(2.,30.);
	double gbps = 8.0 * gbytes / dt;

	stringstream ss;
	ss << "vdif-assembler: processed " << gbytes << " GB at " << gbps << " Gbps\n";
	cout << ss.str() << flush;
    }

    if (mischief_flag) {
	pthread_mutex_lock(&mischief_mutex);
	mischief_killflag = true;
	pthread_mutex_unlock(&mischief_mutex);

	if (pthread_join(mischief_thread, NULL) != 0) {
	    cout << (string("couldn't join mischief thread: ") + strerror(errno) + "\n") << flush;
	    exit(1);
	}
    }

    return 0;
}
