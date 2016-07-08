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

    bool initialized;
    int64_t first_timestamp;
    int64_t last_fileid;
    std::unique_ptr<ch_frb_io::intensity_hdf5_ofile> curr_ofile;
    
    ssize_t nt_alloc;
    std::vector<float> intensity_buf;
    std::vector<float> weights_buf;

    intensity_beam(const string &acqdir);
    virtual ~intensity_beam() { }

    virtual void process_chunk(const shared_ptr<assembled_chunk> &a) override;
    virtual void finalize() override { }
};


intensity_beam::intensity_beam(const string &acqdir_) :
    vdif_processor("intensity_beam", true),
    acqdir(acqdir_), initialized(false), first_timestamp(0), nt_alloc(0)
{ }


// virtual override
void intensity_beam::process_chunk(const shared_ptr<assembled_chunk> &a)
{
    ssize_t nt_chunk_hires = a->nt;
    ssize_t nt_chunk_lores = nt_chunk_hires / nt_downsample;
    const int nfreq = constants::chime_nfreq;

    if ((nt_chunk_hires <= 0) || (nt_chunk_hires != nt_downsample * nt_chunk_lores))
	throw runtime_error("intensity_beam: internal error: bad assembled_chunk size");

    if (!initialized) {
	this->first_timestamp = a->t0;
	this->last_fileid = -1;
	this->initialized = true;
    }

    if (!curr_ofile) {
	int64_t fileid = (a->t0 - first_timestamp) / 400000;
	if ((fileid < 0) || (fileid > 99999999))
	    throw runtime_error("intensity_beam: internal error: bad timestamp");
	if (fileid <= last_fileid)
	    throw runtime_error("intensity_beam: internal error: bad timestamp");

	// Following a convention used throughout CHIMEFRB, we use a filename
	// of the form NNNNNNNN.h5, where N=[0,9].

	stringstream ss;
	ss << acqdir << "/" << setfill('0') << setw(8) << fileid << ".h5";
	
	string filename = ss.str();
	vector<string> pol = { "XX", "YY" };
	double freq0_MHz = 800.0;
	double freq1_MHz = 400.0;
	double dt_sample = 2.56e-6 * nt_downsample;
	ssize_t ipos0 = a->t0;
	double time0 = a->t0 * dt_sample;
	int bitshuffle = 3;    // mandatory compression

	this->curr_ofile = make_unique<ch_frb_io::intensity_hdf5_ofile> (filename, nfreq, pol, freq0_MHz, freq1_MHz, dt_sample, ipos0, time0, bitshuffle);
	this->last_fileid = fileid;
    }

    if (nt_alloc < nt_chunk_lores) {
	nt_alloc = nt_chunk_lores;
	intensity_buf.resize(nfreq * nt_alloc);
	weights_buf.resize(nfreq * nt_alloc);
    }
    
    float *intensityp = &intensity_buf[0];
    memset(intensityp, 0, nfreq * nt_chunk_lores * sizeof(float));

    float *weightsp = &weights_buf[0];
    memset(weightsp, 0, nfreq * nt_chunk_lores * sizeof(float));
    
    // assumed in loop below (sum16_auto_correlations)
    static_assert(nt_downsample % 16 == 0, "intensity_beam::nt_downsample not divisible by 16");

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int ipol = 0; ipol < 2; ipol++) {
	    // pointer to time series at given (ifreq, ipol)
	    const uint8_t *p = a->buf + (2*ifreq+1)*nt_chunk_hires;

	    // loop over lores time samples
	    for (int it = 0; it < nt_chunk_lores; it++) {
		int acc_sum = 0;
		int acc_count = 0;

		// loop over hires time samples which fall in the given lores sample
		for (int it2 = it*nt_downsample; it2 < (it+1)*nt_downsample; it2 += 16) {
		    int sum, count;
		    assembled_chunk::sum16_auto_correlations(sum, count, p+it2);
		    acc_sum += sum;
		    acc_count += count;
		}

		intensityp[ifreq*nt_chunk_lores + it] += (float)acc_sum;
		weightsp[ifreq*nt_chunk_lores + it] += (float)acc_count;
	    }
	}
    }

    for (int i = 0; i < nfreq*nt_chunk_lores; i++) {
	intensityp[i] /= max(weightsp[i], (float)1.0);  // convert weight*intensity -> intensity (avoiding divide-by-zero)
	weightsp[i] /= (2*nt_downsample);               // normalize to max weight 1
    }

    curr_ofile->append_chunk(nt_chunk_lores, intensityp, weightsp, a->t0);

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
