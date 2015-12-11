//
// This file defines some extra helper classes and routines which are
// useful for the python interface, but aren't used in the C++ interface.
//

#ifndef _CH_VDIF_ASSEMBLER_CYTHON_HPP
#define _CH_VDIF_ASSEMBLER_CYTHON_HPP

#include "ch_vdif_assembler_internals.hpp"


namespace ch_vdif_assembler {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// cython_stream: just a wrapper around shared_ptr<vdif_stream>, to avoid confusing cython
//

struct cython_stream {
    std::shared_ptr<vdif_stream> p;
    cython_stream(const std::shared_ptr<vdif_stream> &p_) : p(p_) { xassert(p); }
};

inline cython_stream *cython_file_stream(const std::vector<std::string> &filename_list)
{
    std::shared_ptr<vdif_stream> p = make_file_stream(filename_list);
    return new cython_stream(p);
}

inline cython_stream *cython_simulated_stream(double gbps, double nsec)
{
    std::shared_ptr<vdif_stream> p = make_simulated_stream(gbps, nsec);
    return new cython_stream(p);
}

inline cython_stream *cython_network_stream()
{
    std::shared_ptr<vdif_stream> p = make_network_stream();
    return new cython_stream(p);
}


// -------------------------------------------------------------------------------------------------
//
// cpp_processor: wrapper around shared_ptr<vdif_processor>, to avoid confusing cython
//


struct cpp_processor {
    std::shared_ptr<vdif_processor> p;
    cpp_processor(const std::shared_ptr<vdif_processor> &p_) : p(p_) { xassert(p); }
};

inline cpp_processor *cpp_waterfall_plotter(const std::string &outdir, bool is_critical)
{
    std::shared_ptr<vdif_processor> p = make_waterfall_plotter(outdir, is_critical);
    return new cpp_processor(p);
}


// -------------------------------------------------------------------------------------------------
//
// cython_assembled_chunk: wrapper around shared_ptr<assembled_chunk>, to avoid confusing cython
//

struct cython_assembled_chunk {
    std::shared_ptr<assembled_chunk> p;
    int64_t t0;
    int nt;

    cython_assembled_chunk(const std::shared_ptr<assembled_chunk> &p_) 
	: p(p_) 
    {
	xassert(p);
	this->t0 = p->t0;
	this->nt = p->nt;
    }

    // FIXME how to get a "float complex" pointer from cython?
    inline void fill_efield(void *efield_hack, int32_t *mask)
    {
	if (!efield_hack || !mask)
	    throw std::runtime_error("NULL pointer passed to fill_efield()");

	std::complex<float> *efield = reinterpret_cast<std::complex<float> *> (efield_hack);
	p->fill_efield_array_reference(efield, mask);
    }
};


// -------------------------------------------------------------------------------------------------
//
// cython_assembler


struct cython_assembler {
    vdif_assembler a;
    std::shared_ptr<processor_handle> python_processor;

    cython_assembler(bool write_to_disk, int rbuf_size, int abuf_size, int assembler_nt)
	: a(write_to_disk, rbuf_size, abuf_size, assembler_nt)
    { }

    void register_cpp_processor(cpp_processor *processor)
    {
	xassert(processor);
	xassert(processor->p);
	a.register_processor(processor->p);
    }

    void register_python_processor()
    {
	if (python_processor)
	    throw std::runtime_error("double call to cython_assembler::register_python_processor");
	python_processor = std::make_shared<processor_handle> ("python processor", a.nc);
    }

    // can return NULL
    cython_assembled_chunk *get_next_python_chunk()
    {
	if (!python_processor)
	    throw std::runtime_error("cython_assembler::get_next_python_chunk() called, but no python processor was registered");

	// We don't bother collecting timing statistics for the python processor
	thread_timer timer_unused;

	std::shared_ptr<assembled_chunk> chunk = this->python_processor->get_next_chunk(timer_unused);

	if (!chunk)
	    return NULL;

	return new cython_assembled_chunk(chunk);
    }

    // this seemed like a good idea, before calling wait_until_end()
    void unregister_python_processor()
    {
	python_processor = std::shared_ptr<processor_handle> ();
    }

    void start_async(cython_stream *s)
    {
	xassert(s);
	xassert(s->p);
	a.start_async(s->p);
    }

    void wait_until_end()
    {
	if (python_processor)
	    throw std::runtime_error("cython_assembler::wait_until_end() called with python processor registered");

	a.wait_until_end();
    }
};


}  // namespace ch_vdif_assembler

#endif // _CH_VDIF_ASSEMBLER_CYTHON_HPP
