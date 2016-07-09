// FIXME this code could be streamlined significantly by using the downsampled_intensity helper class

#include <cmath>
#include <cstring>
#include <algorithm>
#include <png.h>

#include "ch_vdif_assembler.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// png_writer: a thin wrapper class around libpng for writing images (nothing in this class
// has anything to do with vdif data)


struct png_writer : noncopyable {
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
    // close file first
    if (fp) {
	fclose(fp);
	fp = NULL;
    }
    
    //
    // The libpng documentation doesn't really explain how to clean up its state.
    //
    // I dug into the source code and decided that one call to png_destroy_write_struct()
    // destroys both the png_ptr and the info_ptr, and does the right thing if only the
    // png_ptr is non-NULL.
    //

    if (png_ptr || info_ptr) {
	png_destroy_write_struct(&png_ptr, &info_ptr);
	png_ptr = NULL;
	info_ptr = NULL;
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
    cout << "wrote " << filename << endl;
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
    
    // This is the main computational routine.
    void absorb_chunk(const assembled_chunk &chunk);

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
void waterfall_image::absorb_chunk(const assembled_chunk &a)
{
    // Number of timestamps per waterfall bin
    int tt = constants::timestamps_per_frame / this->nt_bins;

    // These assumptions simplify the code, but could be removed if necessary
    xassert(this->nfreq == constants::chime_nfreq);
    xassert(constants::timestamps_per_frame % this->nt_bins == 0);
    xassert(a.t0 % tt == 0);
    xassert(a.nt % tt == 0);
    xassert(tt % 16 == 0);

    int b0 = (a.t0 / tt) % nt_bins;    // first waterfall bin in assembled chunk
    int nb = (a.nt / tt);              // number of waterfall bins in assembled chunk

    for (int f = 0; f < nfreq; f++) {
	for (int pol = 0; pol < 2; pol++) {
	    for (int b = 0; b < nb; b++) {
		// Points to the data which will be accumulated into waterfall bin (f,p,b0+b)
		const uint8_t *buf = &a.buf[(2*f+pol) * a.nt + b*tt];
		
		int num = 0;
		int den = 0;
		int sum, count;

		for (int j = 0; j < tt; j += 16) {
		    assembled_chunk::sum16_auto_correlations(sum, count, buf+j);
		    num += sum;
		    den += count;
		}		    

		acc_num[f*2*nt_bins + pol*nt_bins + b+b0] += (float)num;
		acc_den[f*2*nt_bins + pol*nt_bins + b+b0] += (float)den;
	    }
	}
    }
}


void waterfall_image::subtract_medians()
{
    vector<double> buf(nt_bins);
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

    waterfall_plotter(const std::string &outdir, bool is_critical);
    virtual ~waterfall_plotter() { }

    // Devirtualize vdif_processor
    virtual void process_chunk(const shared_ptr<assembled_chunk> &chunk);
    virtual void initialize();
    virtual void finalize();

    void write_images();
};


waterfall_plotter::waterfall_plotter(const string &outdir_, bool is_critical)
    : vdif_processor("waterfall plotter", is_critical), 
      outdir(outdir_),       
      imgdir(outdir_ + string("/img")),
      curr_frame(0),
      curr_frame_initialized(false),
      curr_img(constants::chime_nfreq, 1024)
{
    xmkdir(outdir);
    xmkdir(imgdir);
}


//
// This routine defines the behavior of the waterfall_plotter class.
//
// It calls assembler.get_mask() and assembler.get_visibilities() to get
// the current chunk of data, and absorbs it into the waterfall plot.
//
void waterfall_plotter::process_chunk(const shared_ptr<assembled_chunk> &chunk)
{
    int64_t frame = chunk->t0 / constants::timestamps_per_frame;

    if (!curr_frame_initialized) {
	curr_frame = frame;
	curr_frame_initialized = true;
	curr_img.clear();
    }
    else if (curr_frame != frame) {
	// Frame complete, write and move on to the next
	this->write_images();
	curr_frame = frame;
	curr_img.clear();
    }

    curr_img.absorb_chunk(*chunk);
}


void waterfall_plotter::initialize()
{
    curr_frame_initialized = false;
    curr_img.clear();
}

void waterfall_plotter::finalize()
{
    this->write_images();
    curr_img.clear();
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

shared_ptr<vdif_processor> make_waterfall_plotter(const string &outdir, bool is_critical)
{
    return make_shared<waterfall_plotter> (outdir, is_critical);
}


}   // namespace ch_vdif_assembler
