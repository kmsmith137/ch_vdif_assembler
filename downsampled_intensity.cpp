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
		weightsp[i] = (float)acc_count / (float)nt_downsample;      // normalize to max weight 1
	    }
	}
    }    
}


void downsampled_intensity::process_chunk_reference(const std::shared_ptr<assembled_chunk> &a)
{
    // The first part of this routine is cut-and-paste from process_chunk()
    // FIXME define common helper function

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

    // The logic below deviates from process_chunk()

    vector<float> vis_hires(nfreq * 2 * nt_chunk_hires, 0.0);
    vector<int> mask_hires(nfreq * 2 * nt_chunk_hires, 0);
    a->fill_auto_correlations_reference(&vis_hires[0], &mask_hires[0]);

    memset(&intensity_buf[0], 0, nfreq * 2 * nt_chunk_lores * sizeof(float));
    memset(&weights_buf[0], 0, nfreq * 2 * nt_chunk_lores * sizeof(float));    

    for (int ii = 0; ii < 2*nfreq; ii++) {    // (ifreq,ipol) pairs
	for (int it_hires = 0; it_hires < nt_chunk_hires; it_hires++) {
	    int it_lores = it_hires / nt_downsample;   // integer division
	    int idst = ii*nt_chunk_lores + it_lores;
	    int isrc = ii*nt_chunk_hires + it_hires;

	    if (mask_hires[isrc]) {
		intensity_buf[idst] += vis_hires[isrc];
		weights_buf[idst] += 1.0;
	    }
	}
    }

    for (int i = 0; i < 2*nfreq*nt_chunk_lores; i++) {   // (ifreq,ipol,it_lores) triples
	intensity_buf[i] /= max(weights_buf[i], (float)1.);
	weights_buf[i] /= (float)nt_downsample;
    }
}


}   // namespace ch_vdif_assembler
