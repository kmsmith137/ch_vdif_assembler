#include "ch_vdif_assembler_dspsr.hpp"
#include "ch_vdif_assembler_internals.hpp"

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
    shared_ptr<vdif_assembler> assembler;
    shared_ptr<processor_handle> processor;
    shared_ptr<vdif_stream> stream;
    shared_ptr<assembled_chunk> chunk_reference;
    thread_timer tm;  // not actually used anywhere, but the processor_handle interface requires it

    dspsr_handle_implementation(const string &filelist_filename)
    {
	this->curr_chunk_ix = -1;
	this->curr_data = nullptr;
	this->assembler = make_shared<vdif_assembler> ();
	this->processor = make_shared<processor_handle> ("dspsr", assembler->nc);
	this->stream = make_file_stream(filelist_filename);
	cerr << "ch_vdif_assembler::dspsr_handle:  dspsr_handle constructed\n";
    }

    virtual ~dspsr_handle_implementation() { }

    virtual void advance()
    {
	if (curr_chunk_ix < 0) {
	    cerr << "ch_vdif_assembler::dspsr_handle:  starting assembler\n";
	    assembler->start_async(stream);
	}
	else if (!curr_data)
	    throw runtime_error("ch_vdif_assembler::dspsr_handle::advance() called after end-of-stream");
	    
	this->curr_chunk_ix++;
	this->curr_data = nullptr;
	cerr << "ch_vdif_assembler::dspsr_handle: calling processor->get_next_chunk()\n";
	this->chunk_reference = processor->get_next_chunk(tm);

	// Dropped chunks should never occur, since we haven't implemented a real-time search yet.
	// FIXME: what happens if there is a big gap between timestamps?
	if (processor->ndrops > 0)
	    throw runtime_error("ch_vdif_assembler: dropped chunk?!");

	if (chunk_reference) {
	    cerr << ("ch_vdif_assembler::dspsr_handle: got chunk " + to_string(curr_chunk_ix) + "\n");
	    this->curr_data = chunk_reference->buf;
	    return;
	}

	if (curr_chunk_ix == 0)
	    throw runtime_error("ch_vdif_assembler: empty stream?!");

	cerr << "ch_vdif_assembler::dspsr_handle: end of stream\n";
	return;
    }
};


// static member function
dspsr_handle *dspsr_handle::make(const string &filelist_filename)
{
    return new dspsr_handle_implementation(filelist_filename);
}


}   // namespace ch_vdif_assembler
