from libc.stdint cimport uint8_t, int64_t

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.complex cimport complex


cdef extern from "ch_vdif_assembler.hpp" namespace "ch_vdif_assembler":
    cdef cppclass vdif_stream:
        pass

    cdef cppclass vdif_acquisition(vdif_stream):
        vdif_acquisition(const vector[string] &filename_list, int nthreads) except +


    cdef cppclass vdif_assembler:
        vdif_assembler(bool mask_zeros, bool mask_rails, bool mask_rfi,
                       int rfi_mask_version, double trim_frac, bool span_frames,
                       int buf_nt, bool allow_drops) except +

        int _run(vdif_stream &vs) except +
        void _advance(int nt) except +
        void _allocate() except +
        void _deallocate() except +

        void get_mask(uint8_t *out, int64_t t0, int nt) except +
        void get_complex_data(float complex *out, int64_t t0, int nt) except +
        void get_visibilities(float *out, int64_t t0, int nt) except +

        int64_t _get_buf_t0() except +


    cdef void get_rfi_mask(int *mask, int version) except +
