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


// Helper class.
struct downsampled_intensity {
    const int nt_downsample;  // downsampling factor between baseband and intensity
    const double dt_sample;   // downsampled sample length, int seconds

    bool initialized;
    int64_t initial_t0;
    int64_t curr_chunk_t0;   // downsampled
    int64_t curr_chunk_nt;    // downsampled

    ssize_t nt_alloc;
    std::vector<float> intensity_buf;
    std::vector<float> weights_buf;

    downsampled_intensity(int nt_downsample);

    //
    // When process_chunk() returns, the following members of 'struct downsampled_intensity' will be initialized
    //    curr_chunk_t0     (in downsampled samples)
    //    curr_chunk_nt
    //    intensity_buf     shape (constants::chime_nfreq, 2, curr_chunk_nt)
    //    weights_buf       shape (constants::chime_nfreq, 2, curr_chunk_nt)
    //
    void process_chunk(const std::shared_ptr<assembled_chunk> &a);
};


downsampled_intensity::downsampled_intensity(int nt_downsample_) :
    nt_downsample(nt_downsample_), 
    dt_sample(2.56e-6 * nt_downsample_),
    initialized(false),
    nt_alloc(0)
{
    // Note: it's assumed in downsampled_intensity::process_chunk() that nt_downsample is a multiple of 16
    if (nt_downsample < 16 || nt_downsample > 65536)
	throw runtime_error("downsampled_intensity: currently nt_downsample must be between 16 and 65536");
    if (!is_power_of_two(nt_downsample))
	throw runtime_error("downsampled_intensity: currently nt_downsample must be a power of two");
}


void downsampled_intensity::process_chunk(const shared_ptr<assembled_chunk> &a)
{
    const int nfreq = constants::chime_nfreq;

    ssize_t t0_chunk_hires = a->t0;
    ssize_t nt_chunk_hires = a->nt;
    ssize_t t0_chunk_lores = t0_chunk_hires / nt_downsample;
    ssize_t nt_chunk_lores = nt_chunk_hires / nt_downsample;

    if (t0_chunk_hires < 0)
	throw runtime_error("downsampled_intensity: internal_error: assembled_chunk::t0 was negative?!");
    if (nt_chunk_hires <= 0)
	throw runtime_error("downsampled_intensity: internal_error: assembled_chunk::nt was non-positive?!");

    if (nt_chunk_hires != nt_downsample * nt_chunk_lores)
	throw runtime_error("downsampled_intensity: fatal: assembled_chunk::nt is not divisible by nt_downsample");
    if (t0_chunk_hires != nt_downsample * t0_chunk_lores)
	throw runtime_error("downsampled_intensity: fatal: assembled_chunk::t0 is not divisible by nt_downsample");

    if (!initialized) {
	this->initial_t0 = t0_chunk_lores;
	this->curr_chunk_t0 = t0_chunk_lores;
	this->curr_chunk_nt = 0;
	this->initialized = true;
    }

    if (t0_chunk_lores < curr_chunk_t0 + curr_chunk_nt)
	throw runtime_error("downsampled_intensity: internal error: assembled_chunks are overlapping or have non-monotonic timestamps");

    this->curr_chunk_nt = nt_chunk_lores;
    this->curr_chunk_t0 = t0_chunk_lores;

    if (nt_alloc < nt_chunk_lores) {
	nt_alloc = nt_chunk_lores;
	intensity_buf.resize(nfreq * 2 * nt_alloc);
	weights_buf.resize(nfreq * 2 * nt_alloc);
    }
    
    float *intensityp = &intensity_buf[0];
    float *weightsp = &weights_buf[0];

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int ipol = 0; ipol < 2; ipol++) {
	    // pointer to time series at given (ifreq, ipol)
	    const uint8_t *p = a->buf + (2*ifreq+ipol)*nt_chunk_hires;

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

		int i = (2*ifreq+ipol)*nt_chunk_lores + it;
		intensityp[i] = (float)acc_sum / (float)max(acc_count,1);   // max() avoids divide-by-zero
		weightsp[i] = (float)acc_count / (float)(2*nt_downsample);  // normalize to max weight 1
	    }
	}
    }    
}


// -------------------------------------------------------------------------------------------------




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
