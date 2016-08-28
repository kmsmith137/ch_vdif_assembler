#include "ch_vdif_assembler.hpp"
#include "ch_vdif_assembler_dspsr.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


// static global varibles
const int dspsr_handle::nfreq = constants::chime_nfreq;
const int dspsr_handle::nt_chunk = constants::default_assembler_nt;
const double dspsr_handle::sampling_rate_Hz = 1.0 / constants::dt_fpga;


struct dspsr_handle_implementation : public dspsr_handle
{
    shared_ptr<vdif_stream> stream;

    dspsr_handle_implementation(const string &filelist_filename)
    {
	this->stream = make_file_stream(filelist_filename);
	cerr << "file_stream constructed! " << filelist_filename << endl;
    }

    virtual ~dspsr_handle_implementation() { }
};


// static member function
dspsr_handle *dspsr_handle::make(const string &filelist_filename)
{
    return new dspsr_handle_implementation(filelist_filename);
}


}   // namespace ch_vdif_assembler
