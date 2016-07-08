#include <cstring>
#include "ch_vdif_assembler_internals.hpp"

#ifdef HAVE_CH_FRB_IO
#include <ch_frb_io.hpp>
#endif

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


#ifndef HAVE_CH_FRB_IO

shared_ptr<vdif_processor> make_intensity_beam(const string &acqdir)
{
    throw runtime_error("make_intensity_beam() was called, but this version of ch_vdif_assembler was compiled without libch_frb_io");
}

#else  // HAVE_CH_FRB_IO


struct intensity_beam : public vdif_processor {
    // hardcoded for now, may make a parameter later
    static constexpr int nt_downsample = 512;
    static constexpr int nt_maxfile = 16384;

    string acqdir;
    downsampled_intensity ds_int;
    std::unique_ptr<ch_frb_io::intensity_hdf5_ofile> curr_ofile;

    intensity_beam(const string &acqdir);
    virtual ~intensity_beam() { }

    virtual void process_chunk(const shared_ptr<assembled_chunk> &a) override;
    virtual void finalize() override { }
};


intensity_beam::intensity_beam(const string &acqdir_) :
    vdif_processor("intensity_beam", true),
    acqdir(acqdir_),
    ds_int(nt_downsample)
{ }


// virtual override
void intensity_beam::process_chunk(const shared_ptr<assembled_chunk> &a)
{
    ds_int.process_chunk(a);

    if (!curr_ofile) {
	ssize_t fileid = (ds_int.curr_chunk_t0 - ds_int.initial_t0) / 1000;
	if ((fileid < 0) || (fileid > 99999999))
	    throw runtime_error("intensity_beam: internal error: bad fileid");

	// Following a convention used throughout CHIMEFRB, we use a filename
	// of the form NNNNNNNN.h5, where N=[0,9].

	stringstream ss;
	ss << acqdir << "/" << setfill('0') << setw(8) << fileid << ".h5";
	
	string filename = ss.str();
	int nfreq = constants::chime_nfreq;
	vector<string> pol = { "XX", "YY" };
	double freq0_MHz = 800.0;
	double freq1_MHz = 400.0;
	double dt_sample = ds_int.dt_sample;
	ssize_t ipos0 = ds_int.curr_chunk_t0;
	double time0 = ds_int.curr_chunk_t0 * dt_sample;
	int bitshuffle = 3;    // mandatory compression

	this->curr_ofile = make_unique<ch_frb_io::intensity_hdf5_ofile> (filename, nfreq, pol, freq0_MHz, freq1_MHz, dt_sample, ipos0, time0, bitshuffle);
    }

    curr_ofile->append_chunk(ds_int.curr_chunk_nt, &ds_int.intensity_buf[0], &ds_int.weights_buf[0], ds_int.curr_chunk_t0);

    if (curr_ofile->curr_nt >= nt_maxfile)
	this->curr_ofile = unique_ptr<ch_frb_io::intensity_hdf5_ofile> ();  // emptying this pointer flushes and writes file
}


shared_ptr<vdif_processor> make_intensity_beam(const string &acqdir)
{
    xmkdir(acqdir);
    
    if (!is_empty_dir(acqdir))
	throw runtime_error("make_intensity_beam(): fatal: acquisition directory '" + acqdir + "' is nonempty");

    return make_shared<intensity_beam> (acqdir);
}

#endif  // HAVE_CH_FRB_IO


}   // namespace ch_vdif_assembler
