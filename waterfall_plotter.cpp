#include <errno.h>
#include <sys/stat.h>
#include <cstring>
#include <string>
#include <stdexcept>
#include <boost/make_shared.hpp>
#include <boost/noncopyable.hpp>

#include <png.h>
#include "ch_vdif_assembler.hpp"

using namespace std;
using namespace boost;

#ifndef xassert
#define xassert(cond) xassert2(cond, __LINE__)
#define xassert2(cond,line) do { if (!(cond)) throw std::runtime_error("Assevdifon '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")"); } while (0)
#endif

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// png_writer: a thin wrapper class around libpng for writing images (nothing in this class
// has anything to do with vdif data)


struct png_writer : boost::noncopyable {
    png_structp png_ptr;
    png_infop info_ptr;
    FILE *fp;

    png_writer();
    ~png_writer();

    //
    // @rgb should be an array of shape (ny,nx,3) where the last index is "rgb".
    // The first index (y) runs from top to bottom (not bottom to top) in the image.
    // 
    void write(const string &filename, png_byte *rgb, int nx, int ny);
    void deallocate();
};


png_writer::png_writer()
    : png_ptr(NULL), info_ptr(NULL), fp(NULL)
{ }

png_writer::~png_writer()
{
    this->deallocate();
}


void png_writer::deallocate()
{
    if (fp) {
	fclose(fp);
	fp = NULL;
    }
    
    if (info_ptr) {
	png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	info_ptr = NULL;
    }
    
    if (png_ptr) {
	png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	png_ptr = NULL;
    }
}

void png_writer::write(const string &filename, png_byte *rgb, int nx, int ny)
{
    // just in case any state remains from a previous call
    this->deallocate();

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    xassert(png_ptr);

    info_ptr = png_create_info_struct(png_ptr);
    xassert(info_ptr);

    fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        cerr << "png_writer: couldn't open file " << filename << endl;
	this->deallocate();
        return;
    }
    
    if (setjmp(png_jmpbuf(png_ptr))) {
	cerr << "png_writer: libpng internal failure\n";
	this->deallocate();
	return;
    }

    png_init_io(png_ptr, fp);

    png_set_IHDR(png_ptr, info_ptr, nx, ny,
		 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    for (int i = 0; i < ny; i++)
	png_write_row(png_ptr, &rgb[3*i*nx]);

    png_write_end(png_ptr, NULL);

    this->deallocate();
    cerr << "wrote " << filename << "\n";
}


// -------------------------------------------------------------------------------------------------
//
// This class really represents a pair of images, one for each polarization


struct waterfall_image {
    int nfreq;
    int nt_bins;

    // These buffers have shape (nfreq, 2, nt_bins) where the second index is polarization
    std::vector<double> acc_num;
    std::vector<double> acc_den;

    png_writer pw;
    
    // Construct empty image
    waterfall_image(int nfreq, int nt_bins);

    // Construct image by downgrading existing image (e.g. to make html thumbnail)
    waterfall_image(const waterfall_image &w, int new_nfreq, int new_nt_bins);
    
    //
    // This is the main computational routine.
    //
    // It takes a chunk of data from the vdif_assembler, and absorbs it into the
    // waterfall plot.  
    //
    // Arugments t0,nt have the same meanings as in vdif_processor::process_data().
    //
    // The 'vis' and 'mask' arrays are packed as in vdif_assembler::get_visibilities()
    // and vdif_assembler::get_mask(), i.e. vis is an array of shape (nfreq, 4, nt)
    // where the middle index is {xauto,yauto,recross,imcross}, and mask is an array
    // of shape (nfreq, 2, nt) where the middle index is polarization.
    //
    void absorb_data(int64_t t0, int nt, const float *vis, const uint8_t *mask);

    void clear();
    void subtract_medians();

    // The mask here has shape (nfreq, nt_bins) and applies to a single polarization
    double compute_variance_outside_mask(int pol, const int *mask);

    // Generates a mask internally using iterative clipping
    double compute_clipped_variance(int pol);
    
    void write_image(const string &filename, int pol, double dv);
    void write_image(const string &filename, int pol);
};


// Construct empty image
waterfall_image::waterfall_image(int nfreq_, int nt_bins_)
    : nfreq(nfreq_), nt_bins(nt_bins_),
      acc_num(nfreq_ * 2 * nt_bins_, 0),
      acc_den(nfreq_ * 2 * nt_bins_, 0)
{ 
    xassert(nfreq > 0);
    xassert(nt_bins > 0);
}


// Construct image by downgrading existing image (e.g. to make html thumbnail)
waterfall_image::waterfall_image(const waterfall_image &w, int new_nfreq, int new_nt_bins)
    : nfreq(new_nfreq), nt_bins(new_nt_bins),
      acc_num(new_nfreq * 2 * new_nt_bins, 0),
      acc_den(new_nfreq * 2 * new_nt_bins, 0)
{
    xassert(nfreq > 0);
    xassert(nt_bins > 0);

    xassert(w.nfreq % nfreq == 0);
    xassert(w.nt_bins % nt_bins == 0);

    int ff = w.nfreq / nfreq;
    int bb = w.nt_bins / nt_bins;

    for (int f = 0; f < w.nfreq; f++) {
	for (int pol = 0; pol < 2; pol++) {
	    int isrc = (2*f + pol) * w.nt_bins;
	    int idst = (2*(f/ff) + pol) * nt_bins;

	    for (int b = 0; b < w.nt_bins; b++) {
		acc_num[idst + (b/bb)] += w.acc_num[isrc + b];
		acc_den[idst + (b/bb)] += w.acc_den[isrc + b];
	    }
	}
    }
}


//
// Since this routine is the rate-limiting step of the waterfall_assembler,
// it's written with a lot of optimization
//
void waterfall_image::absorb_data(int64_t t0_, int nt, const float *vis, const uint8_t *mask)
{
    // This routine assumes we're running at native CHIME frequency resolution
    xassert(nfreq == vdif_assembler::nfreq);

    // This routine assumes that the number of time bins evenly divides the CHIME frame count
    xassert((1 << 23) % nt_bins == 0);

    // Number of timestamps per bin
    int tt = (1<<23) / nt_bins;

    // Shift t0 to beginning of frame
    int t0 = (t0_ % (1 << 23));

    // Bin range
    int b0 = t0/tt;
    int b1 = ((t0+nt-1) / tt) + 1;

    // This should be guaranteed but never hurts to check
    xassert(b1 <= nt_bins);

    for (int f = 0; f < nfreq; f++) {
	for (int pol = 0; pol < 2; pol++) {
	    const float *vis_row = &vis[(4*f+pol)*nt];
	    const uint8_t *mask_row = &mask[(2*f+pol)*nt];

	    for (int b = b0; b < b1; b++) {
		// range of time indices falling in bin b
		int j0 = max(b*tt-t0, 0);
		int j1 = min((b+1)*tt-t0, nt);
		
		float num = 0.0;
		int den = 0;
		
		for (int j = j0; j < j1; j++) {
		    num += vis_row[j];
		    den += mask_row[j];
		}

		acc_num[f*2*nt_bins + pol*nt_bins + b] += num;
		acc_den[f*2*nt_bins + pol*nt_bins + b] += den;
	    }
	}
    }
}


void waterfall_image::subtract_medians()
{
    std::vector<double> buf(nt_bins);
    int nvals;

    // loop over (freq,pol pairs)
    for (int i = 0; i < 2*nfreq; i++) {
	nvals = 0;
	for (int j = 0; j < nt_bins; j++) {
	    if (acc_den[i*nt_bins+j] > 99.5)
		buf[nvals++] = acc_num[i*nt_bins+j] / acc_den[i*nt_bins+j];
	}
	
	if (!nvals)
	    continue;

	std::nth_element(buf.begin(), buf.begin() + (nvals/2), buf.end());
	double median = buf[nvals/2]; 

	for (int j = 0; j < nt_bins; j++)
	    acc_num[i*nt_bins+j] -= median * acc_den[i*nt_bins+j];
    }
}


void waterfall_image::clear()
{
    memset(&acc_num[0], 0, acc_num.size() * sizeof(acc_num[0]));
    memset(&acc_den[0], 0, acc_den.size() * sizeof(acc_den[0]));
}


//
// @mask assumed 0 or 1
//
double waterfall_image::compute_variance_outside_mask(int pol, const int *mask)
{
    double num2 = 0.0;
    double den2 = 0.0;

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int b = 0; b < nt_bins; b++) {
	    if (mask[ifreq*nt_bins + b]) {
		double num = acc_num[ifreq*(2*nt_bins) + (pol*nt_bins) + b];
		double den = acc_den[ifreq*(2*nt_bins) + (pol*nt_bins) + b];
		num2 += (num*num)/(den*den);
		den2 += 1.0;
	    }
	}
    }

    return (den2 > 0.0) ? (num2/den2) : 0.0;
}


double waterfall_image::compute_clipped_variance(int pol)
{
    vector<int> mask(nfreq * nt_bins, 0);

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int b = 0; b < nt_bins; b++) { 
	    double den = acc_den[ifreq*(2*nt_bins) + (pol*nt_bins) + b];
	    mask[ifreq*nt_bins + b] = (den > 99.5);
	}
    }

    double vvar = compute_variance_outside_mask(pol, &mask[0]);

    if (vvar <= 0.0)
	return 0.0;

    for (int n = 0; n < 3; n++) {
	// extend mask by clipping at 2 sigma
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int b = 0; b < nt_bins; b++) {		
		double num = acc_num[ifreq*(2*nt_bins) + (pol*nt_bins) + b];
		double den = acc_den[ifreq*(2*nt_bins) + (pol*nt_bins) + b];
		
		int ok = (num*num < 4*vvar*den*den);
		mask[ifreq*nt_bins + b] &= ok;
	    }
	}

	// recompute variance with extended mask
	double vvar2 = compute_variance_outside_mask(pol, &mask[0]);
	if (vvar2 <= 0.0)
	    return vvar;

	vvar = vvar2;
    }

    return vvar;
}


void waterfall_image::write_image(const string &filename, int pol, double dv)
{
    vector<png_byte> rgb(3 * nfreq * nt_bins, 0);

    if (dv <= 0.0) {
	pw.write(filename, &rgb[0], nt_bins, nfreq);
	return;
    }

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int b = 0; b < nt_bins; b++) {
	    double num = acc_num[ifreq*(2*nt_bins) + (pol*nt_bins) + b];
	    double den = acc_den[ifreq*(2*nt_bins) + (pol*nt_bins) + b];

	    if (den > 99.5) {
		int val = 128 * (num/den/dv) + 128;
		val = max(val, 0);
		val = min(val, 255);

		rgb[ifreq*(3*nt_bins) + 3*b] = (unsigned int)val;
		rgb[ifreq*(3*nt_bins) + 3*b + 2] = (unsigned int)(255-val);
	    }
	}
    }

    pw.write(filename, &rgb[0], nt_bins, nfreq);
}

void waterfall_image::write_image(const string &filename, int pol)
{
    double vvar = this->compute_clipped_variance(pol);
    double dv = 2.5 * sqrt(max(vvar,0.0));
    this->write_image(filename, pol, dv);
}


// -------------------------------------------------------------------------------------------------


struct waterfall_plotter : public vdif_processor
{
    std::string outdir;
    std::string imgdir;

    int64_t curr_frame;
    bool curr_frame_initialized;

    waterfall_image curr_img;

    waterfall_plotter(const std::string &outdir);

    // Devirtualize vdif_processor
    virtual void process_data(vdif_assembler &assembler, int64_t t0, int nt);
    virtual void finalize();
    virtual ~waterfall_plotter() { }

    void write_images();
    static void make_dir(const std::string &dirname);
};


waterfall_plotter::waterfall_plotter(const string &outdir_)
    : outdir(outdir_),       
      imgdir(outdir_ + string("/img")),
      curr_frame(0),
      curr_frame_initialized(false),
      curr_img(vdif_assembler::nfreq, 1024)
{
    make_dir(outdir);
    make_dir(imgdir);
}


//
// This routine defines the behavior of the waterfall_plotter class.
//
// It calls assembler.get_mask() and assembler.get_visibilities() to get
// the current chunk of data, and absorbs it into the waterfall plot.
//
void waterfall_plotter::process_data(vdif_assembler &assembler, int64_t t0, int nt)
{
    if (assembler.span_frames) {
	cerr << "waterfall plotter only works if assembler.span_frames is set to false\n";
	throw runtime_error("waterfall plotter only works if assembler.span_frames is set to false");
    }

    int64_t frame = t0 / (1<<23);

    if (!curr_frame_initialized) {
	curr_frame = frame;
	curr_frame_initialized = true;
    }
    else if (curr_frame != frame) {
	// Frame complete, write and move on to the next
	this->write_images();
	curr_frame = frame;
	curr_img.clear();
    }

    vector<uint8_t> mask(vdif_assembler::nfreq * 2 * nt, 0);
    assembler.get_mask(&mask[0], t0, nt);

    vector<float> vis(vdif_assembler::nfreq * 4 * nt, 0);
    assembler.get_visibilities(&vis[0], t0, nt);

    curr_img.absorb_data(t0, nt, &vis[0], &mask[0]);
}


void waterfall_plotter::finalize()
{
    this->write_images();
}


void waterfall_plotter::write_images()
{
    curr_img.subtract_medians();

    waterfall_image thumbnail(curr_img, 256, 256);
    thumbnail.subtract_medians();
    
    for (int pol = 0; pol < 2; pol++) {
	stringstream ss1;
	ss1 << imgdir << "/full_pol" << pol << "_frame" << curr_frame << ".png";
	curr_img.write_image(ss1.str(), pol);

	stringstream ss2;
	ss2 << imgdir << "/thumbnail_pol" << pol << "_frame" << curr_frame << ".png";
	thumbnail.write_image(ss2.str(), pol);
    }
}


void waterfall_plotter::make_dir(const string &dirname)
{
    int err = mkdir(dirname.c_str(), 0777);

    if (!err)
	return;

    if (errno != EEXIST) {
	stringstream ss;
	ss << "couldn't create directory " << dirname << ": " << strerror(errno);
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }
    
    struct stat s;
    err = stat(dirname.c_str(), &s);
    
    if (err < 0) {
	stringstream ss;
	ss << "couldn't stat file " << dirname << ": " << strerror(errno);
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }

    if (!S_ISDIR(s.st_mode)) {
	stringstream ss;
	ss << "couldn't create directory " << dirname << ": file already exists and is not a directory";
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }
}


shared_ptr<vdif_processor> make_waterfall_plotter(const string &outdir)
{
    return boost::make_shared<waterfall_plotter> (outdir);
}


}  // namespace ch_vdif_assembler

/*
 * Local variables:
 *  c-basic-offset: 4
 * End:
 */
