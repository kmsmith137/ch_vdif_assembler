//
// rfi_histogrammer: this is work in progress toward real-time RFI analysis.
//
// The code is not cleaned up or commented systematically, so it's not a good
// reference for writing vdif_processors.
//

#include <cstring>
#include <stdexcept>

#include <hdf5.h>
#include "ch_vdif_assembler.hpp"

using namespace std;

// backwards compatibility hacks for HDF5 1.6
#if H5_VERS_MINOR == 6
#  define H5Aiterate1 H5Aiterate
#  define H5Acreate1 H5Acreate
#  define H5Dopen1 H5Dopen
#  define H5Dcreate1 H5Dcreate
#  define H5Gcreate1 H5Gcreate
#  define H5Eset_auto1 H5Eset_auto
#  define H5Ewalk1 H5Ewalk
#  define H5E_error1_t H5E_error_t
#  define H5Gopen1 H5Gopen
#endif

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


struct histogram {
    int ifreq;
    int pol;
    int p;
    int n;
    double dx;
    vector<int64_t> counts;

    histogram(int ifreq, int pol, int p, int n, double dx);

    inline void add(double val)
    {
	int i = (int)(min(val,n*dx) / dx);
	i = max(i, 0);
	i = min(i, n-1);
	counts[i] += 1.0;
    }

    inline void add(double num1, double den1, double num2, double den2, double thresh, int64_t bin, int p)
    {
	// skip if noise source boundary is straddled
	if ((bin % (1 << (23-p))) == 0)
	    return;

	// skip if not enough samples
	if ((den1 < thresh) || (den2 < thresh))
	    return;

	double val = fabs(num1/den1 - num2/den2);
	this->add(val);
    }
};


histogram::histogram(int ifreq_, int pol_, int p_, int n_, double dx_)
    : ifreq(ifreq_), pol(pol_), p(p_), n(n_), dx(dx_), counts(n_,0)
{
    if (n <= 0)
	throw runtime_error("bad value of n passed to histogram constructor");
    if (dx <= 0)
	throw runtime_error("bad value of dx passed to histogram constructor");
}


// -------------------------------------------------------------------------------------------------


struct histogram_set {
    static const int pmin = 5;
    static const int pmax = 18;

    //
    // buf is a shape-(pmax+1, 2, 2) array
    //   first index = p
    //   second index = time index in a 2-sample buffer
    //   third index = num/den
    //
    vector<double> buf;

    // length (pmax+1)
    vector<histogram> histograms;

    bool empty;
    int64_t tcurr;

    histogram_set(int ifreq, int pol);

    // not written for speed, intended as a reference
    void add_sample(int64_t t, double vis);

    // written to be fast, unit tested by comparing to add_sample()
    void add_samples(int64_t t0, int nt, const float *vis, const int *counts, bool ref_flag=false);

    // helper for add_samples()
    void _add_binned_samples(int p, int64_t b0, int64_t b1, const double *vis_binned, const double *counts_binned);

    void finalize();
};


histogram_set::histogram_set(int ifreq, int pol)
    : buf(4*(pmax+1),0), empty(true), tcurr(0)
{
    // dummy histograms
    for (int p = 0; p < pmin; p++)
	histograms.push_back(histogram(ifreq,pol,p,1,1));

    for (int p = pmin; p <= pmax; p++) {
	// reasonable defaults?
	int n = 512;
	double dx = 1. / pow(2,p/2.);
	histograms.push_back(histogram(ifreq,pol,p,n,dx));
    }
}


void histogram_set::add_sample(int64_t t, double vis)
{	
    if (empty) {
	tcurr = t;
	empty = false;
    }

    for (int p = pmin; p <= pmax; p++) {
	int64_t pp = (1 << p);
	double thresh = 0.5*(pp+1);

	int64_t b = t/pp;
	int64_t bcurr = tcurr/pp;

	if (b <= bcurr) {
	    buf[4*p+2] += vis;
	    buf[4*p+3] += 1.0;
	    continue;
	}
	
	histograms[p].add(buf[4*p], buf[4*p+1], buf[4*p+2], buf[4*p+3], thresh, bcurr, p);
	
	if (b == bcurr+1) {
	    buf[4*p] = buf[4*p+2];
	    buf[4*p+1] = buf[4*p+3];
	}
	else {
	    buf[4*p] = 0.0;
	    buf[4*p+1] = 0.0;
	}
	
	buf[4*p+2] = vis;
	buf[4*p+3] = 1.0;
    }

    tcurr = t;
}


// assumes b0 < b1
void histogram_set::_add_binned_samples(int p, int64_t b0, int64_t b1, const double *vis_binned, const double *counts_binned)
{
    int64_t pp = 1 << p;
    int64_t bcurr = tcurr/pp;
    double thresh = 0.5*(pp+1);

    double *bufp = &buf[4*p];
    histogram &histp = histograms[p];

    if (b0 <= bcurr) {
	// update current bin
	bufp[2] += vis_binned[0];
	bufp[3] += counts_binned[0];
	vis_binned++;
	counts_binned++;
	b0 = bcurr + 1;
    }
    
    if (b1 <= bcurr+1)
	return;

    // if we get here, then bcurr < b0 < b1, i.e. we have data beyond the current bin
    histp.add(bufp[0], bufp[1], bufp[2], bufp[3], thresh, bcurr, p);

    if (b1 == bcurr+2) {
	// if we get here, then (b0,b1)=(bcurr+1,bcurr+2), i.e. we're just extending by one partial bin
	bufp[0] = bufp[2];
	bufp[1] = bufp[3];
	bufp[2] = vis_binned[0];
	bufp[3] = counts_binned[0];
	return;
    }

    if (b0 == bcurr+1) {
	// if we get here, then we have bin (bcurr+1) and it is complete
	histp.add(bufp[2], bufp[3], vis_binned[0], counts_binned[0], thresh, bcurr+1, p);
    }

    // OK if this loop is empty
    for (int64_t b = b0; b <= b1-3; b++)
	histp.add(vis_binned[b-b0], counts_binned[b-b0], vis_binned[b-b0+1], counts_binned[b-b0+1], thresh, b+1, p);

    bufp[0] = (b1 >= b0+2) ? vis_binned[b1-b0-2] : 0.0;
    bufp[1] = (b1 >= b0+2) ? counts_binned[b1-b0-2] : 0.0;
    bufp[2] = vis_binned[b1-b0-1];
    bufp[3] = counts_binned[b1-b0-1];
}


void histogram_set::add_samples(int64_t t0, int nt, const float *vis_arr, const int *flag_arr, bool ref_flag)
{
    if (ref_flag) {
	for (int i = 0; i < nt; i++) {
	    if (flag_arr[i])
		this->add_sample(t0+i, vis_arr[i]);
	}
	return;
    }

    if (empty) {
	tcurr = t0;
	empty = false;
    }

    //
    // First step: bin data to width 2^pmin
    //
    int64_t pp = (1 << pmin);
    int64_t b0 = t0/pp;              // start bin in data
    int64_t b1 = (t0+nt-1)/pp + 1;   // end bin in data
    int64_t nb = b1 - b0;            // number of bins spanned by data

    // allocate temp arrays
    vector<double> scratch(4*nb, 0);
    double *vis_binned = &scratch[0];
    double *counts_binned = &scratch[nb];
    double *vis_tmp = &scratch[2*nb];
    double *counts_tmp = &scratch[3*nb];

    for (int64_t t = t0; t < t0+nt; t++) {
	int64_t b = t/pp;
	vis_binned[b-b0] += vis_arr[t-t0];
	counts_binned[b-b0] += flag_arr[t-t0];
    }

    this->_add_binned_samples(pmin, b0, b1, vis_binned, counts_binned);

    for (int p = pmin+1; p <= pmax; p++) {
	//
	// Rebin (p-1) -> p
	//
	int64_t new_b0 = b0/2;
	int64_t new_b1 = (b1-1)/2 + 1;

	memset(vis_tmp, 0, (new_b1-new_b0) * sizeof(double));
	memset(counts_tmp, 0, (new_b1-new_b0) * sizeof(double));

	for (int64_t b = b0; b < b1; b++) {
	    vis_tmp[(b/2)-new_b0] += vis_binned[b-b0];
	    counts_tmp[(b/2)-new_b0] += counts_binned[b-b0];
	}

	std::swap(vis_binned, vis_tmp);
	std::swap(counts_binned, counts_tmp);
	b0 = new_b0;
	b1 = new_b1;

	//
	// Bin samples!
	//
	this->_add_binned_samples(p, b0, b1, vis_binned, counts_binned);
    }

    this->tcurr = t0+nt-1;
}


void histogram_set::finalize()
{
    if (empty)
	return;

    for (int p = pmin; p <= pmax; p++) {
	int64_t pp = (1 << p);
	double thresh = 0.5*(pp+1);
	histograms[p].add(buf[4*p], buf[4*p+1], buf[4*p+2], buf[4*p+3], thresh, tcurr/pp, p);
    }
}


// -------------------------------------------------------------------------------------------------


struct rfi_histogrammer : public vdif_processor {    
    vector<histogram_set> histograms;   // shape (nfreq, pol)
    string output_hdf5_filename;
    bool ref_flag;

    // If @ref_flag is specified, then the reference implementation will be used (slow but simple code)
    rfi_histogrammer(const string &output_hdf5_filename, bool is_critical, bool ref_flag);
    virtual ~rfi_histogrammer() { }

    void write_hdf5_file(const string &filename) const;

    // Devirtualize vdif_assembler_callback
    virtual void initialize() { }
    virtual void process_chunk(const shared_ptr<assembled_chunk> &a);
    virtual void finalize();
};


rfi_histogrammer::rfi_histogrammer(const string &output_hdf5_filename_, bool is_critical_, bool ref_flag_)
    : vdif_processor("rfi_histogrammer", is_critical_),
      output_hdf5_filename(output_hdf5_filename_), 
      ref_flag(ref_flag_)
{ 
    for (int ifreq = 0; ifreq < 1024; ifreq++)
	for (int pol = 0; pol < 2; pol++)
	    histograms.push_back(histogram_set(ifreq, pol));

    if (ref_flag)
	cerr << "Note: using rfi_histogrammer reference implementation\n";
}


void rfi_histogrammer::write_hdf5_file(const string &filename) const
{
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (file_id < 0) {
        cerr << "Fatal: couldn't create HDF5 file " << filename << endl;
        exit(1);
    }

    hid_t group_id = H5Gopen1(file_id, ".");
    xassert(group_id >= 0);

    for (int p = histogram_set::pmin; p <= histogram_set::pmax; p++) {
	// use single precision floating-point to allow large dynamic range while saving disk space
	int nbins = this->histograms[0].histograms[p].n;
	vector<float> all_counts(histograms.size() * nbins, 0);

	for (unsigned int i = 0; i < histograms.size(); i++) {
	    const histogram &h = this->histograms[i].histograms[p];
	    xassert(h.n == nbins);

	    for (int j = 0; j < nbins; j++)
		all_counts[i*nbins + j] = (float)h.counts[j];
	}

	vector<hsize_t> shape(3);
	shape[0] = 1024;   // frequency
	shape[1] = 2;      // polarization
	shape[2] = nbins;  
	
	hid_t dataspace_id = H5Screate(H5S_SIMPLE);
	xassert(dataspace_id >= 0);

	int ret = H5Sset_extent_simple(dataspace_id, 3, &shape[0], &shape[0]);
	xassert(ret >= 0);

	stringstream s;
	s << "BIN" << (1 << p);
	string dataset_name = s.str();
	
	hid_t dataset_id = H5Dcreate1(group_id, dataset_name.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
	xassert(dataset_id >= 0);

	ret = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &all_counts[0]);
	xassert(ret >= 0);
	
	cerr << filename << ": wrote shape (1024,2," << nbins << ") array\n";
    }

    H5Fclose(file_id);
    cerr << filename << ": closed file\n";
}


void rfi_histogrammer::process_chunk(const shared_ptr<assembled_chunk> &chunk)
{
    int nt = chunk->nt;
    int64_t t0 = chunk->t0;

    vector<float> vis(2048 * nt, 0);
    vector<int> mask(2048 * nt, 0);
    chunk->fill_auto_correlations_reference(&vis[0], &mask[0]);

    for (int ifreq = 0; ifreq < 1024; ifreq++) {
	histograms[2*ifreq].add_samples(t0, nt, &vis[2*ifreq*nt], &mask[2*ifreq*nt], this->ref_flag);
	histograms[2*ifreq+1].add_samples(t0, nt, &vis[(2*ifreq+1)*nt], &mask[(2*ifreq+1)*nt], this->ref_flag);
    }
}


// empty for now
void rfi_histogrammer::finalize()
{
    for (unsigned int i = 0; i < histograms.size(); i++)
	histograms[i].finalize();

    this->write_hdf5_file(output_hdf5_filename);
}


shared_ptr<vdif_processor> make_rfi_histogrammer(const string &output_hdf5_filename, bool is_critical, bool ref_flag)
{
    return make_shared<rfi_histogrammer> (output_hdf5_filename, is_critical, ref_flag);
}


}  // namespace ch_vdif_assembler

/*
 * Local variables:
 *  c-basic-offset: 4
 * End:
 */
