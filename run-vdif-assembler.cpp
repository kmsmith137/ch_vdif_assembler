#include <cstring>
#include "ch_vdif_assembler.hpp"

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
	 << "    -n                                 to run on a real-time network capture (in 'alpha')\n"
	 << "\n"
	 << "You may find the script show-moose-acqusitions.py useful for making file lists\n"
	 << "\n"
	 << "Suggested usage:\n"
	 << "    show-moose-acquisitions.py 41537 > filelist_41537.txt\n"
	 << "    run-vdif-assembler -w waterfall_41537 -f filelist_41537.txt\n"
	 << "    ./index-thumbnails.py waterfall_41537\n"
	 << "\n"
	 << "The last script makes an HTML summary page with clickable thumbnails\n";

    exit(2);
}


int main(int argc, char **argv)
{
    string waterfall_plot_outdir;
    string rfi_histogram_outfile;
    string file_list;
    bool rfi_rflag = false;
    bool nflag = false;

    //
    // Low-tech command line parsing
    //
    int pos = 1;
    while (pos < argc) {
	if ((strlen(argv[pos]) != 2) || (argv[pos][0] != '-'))
	    usage();

	if (argv[pos][1] == 'n') {
	    if (nflag || !file_list.empty())
		usage();
	    nflag = true;
	    pos++;
	    continue;
	}

	if ((pos+1 >= argc) || (argv[pos+1][0] == '-'))
	    usage();

	if (argv[pos][1] == 'w') {
	    if (!waterfall_plot_outdir.empty())
		usage();
	    waterfall_plot_outdir = argv[pos+1];
	}
	else if (argv[pos][1] == 'r') {
	    if (!rfi_histogram_outfile.empty())
		usage();
	    rfi_histogram_outfile = argv[pos+1];
	}
	else if (argv[pos][1] == 'R') {
	    if (!rfi_histogram_outfile.empty())
		usage();
	    rfi_histogram_outfile = argv[pos+1];
	    rfi_rflag = true;
	}
	else if (argv[pos][1] == 'f') {
	    if (nflag || !file_list.empty())
		usage();
	    file_list = argv[pos+1];
	}
	else
	    usage();

	pos += 2;
    }

    if (!nflag && file_list.empty())
	usage();
    if (waterfall_plot_outdir.empty() && rfi_histogram_outfile.empty())
	usage();

    //
    // Construct assembler object
    //
    vdif_assembler assembler;

    //
    // Construct processor objects
    //
    if (!waterfall_plot_outdir.empty())
	assembler.register_processor(make_waterfall_plotter(waterfall_plot_outdir));
    if (!rfi_histogram_outfile.empty())
	assembler.register_processor(make_rfi_histogrammer(rfi_histogram_outfile, rfi_rflag));

    //
    // Construct stream object and run assembler
    //
    if (!file_list.empty()) {
	vdif_acquisition stream(file_list);
	assembler.run(stream);
    }
    else {
	// real-time network capture
	vdif_rt_stream stream;
	assembler.run(stream);
    }

    return 0;
}
