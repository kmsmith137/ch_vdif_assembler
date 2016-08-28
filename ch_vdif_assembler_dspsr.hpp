//
// We sometimes use the vdif_assembler to supply data to the 'dspsr'
// pulsar data processing software.  One technical nuisance issue is
// that dspsr is compiled using a C++03 compiler, whereas ch_vdif_assembler
// assumes C++11.
//
// This header file is C++03 compliant and defines a virtual base class
// ch_vdif_assembler::dspsr_handle which can be used from the dspsr code,
// and a factory function dspsr_handle::make().  The factory function (and
// corresponding non-virtual subclass) are implemented in dspsr_handle.cpp
// where C++11 features can be used.
//

#ifndef _CH_VDIF_ASSEMBLER_DSPSR_HPP
#define _CH_VDIF_ASSEMBLER_DSPSR_HPP

// Note: no C++11 includes (including ch_vdif_assembler.hpp!)
#include <string>

namespace ch_vdif_assembler {
#if 0
}; // pacify emacs c-mode
#endif


struct dspsr_handle {
    static const int nfreq;
    static const int nt_chunk;
    static const double sampling_rate_Hz;
    
    // 
    // When the dspsr_handle is constructed, curr_chunk_ix is set to (-1) and curr_data is set to NULL.
    // In the first call to advance(), the stream is started, curr_chunk_ix is set to 0, and curr_data
    // points to the (4+4)-bit data.  Subsequent calls to advance() will modify curr_chunk_ix and curr_data
    // appropriately until end-of-stream is reached, at which point curr_data is set to NULL.
    // 
    int64_t curr_chunk_ix;
    const uint8_t *curr_data;

    virtual void advance() = 0;
    virtual ~dspsr_handle() { }
    
    static dspsr_handle *make(const std::string &filelist_filename);
};


}  // namespace ch_vdif_assembler

#endif // _CH_VDIF_ASSEMBLER_DSPSR_HPP
