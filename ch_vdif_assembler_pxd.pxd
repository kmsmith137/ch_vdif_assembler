from libc.stdint cimport int32_t, int64_t

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from "ch_vdif_assembler_cython.hpp" namespace "ch_vdif_assembler::constants":
    cdef int chime_nfreq
    cdef int timestamps_per_frame
    cdef int num_disks


cdef extern from "ch_vdif_assembler_cython.hpp" namespace "ch_vdif_assembler":
    cdef cppclass cython_stream:
        pass

    cdef cppclass cpp_processor:
        pass    # opaque to cython


    cdef cppclass cython_assembled_chunk:
        int64_t t0
        int nt

        void fill_efield(float complex *efield, int32_t *mask) except +


    cdef cppclass cython_assembler:
        cython_assembler(bool write_to_disk, int rbuf_size, int abuf_size, int assembler_nt) except +

        void register_cpp_processor(cpp_processor *p) except +
        void register_python_processor() except +
        void unregister_python_processor() except +
        void start_async(cython_stream *s) except +
        void wait_until_end() except +

        cython_assembled_chunk *get_next_python_chunk() except +


    # Factory functions
    cython_stream *cython_file_stream(const vector[string] &filename_list) except +
    cython_stream *cython_simulated_stream(double gbps, double nsec) except +
    cython_stream *cython_network_stream() except +
    cpp_processor *cpp_waterfall_plotter(const string &outdir, bool is_critical) except +

