import numpy as np
cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport int64_t

cimport _ch_vdif_assembler_c

# FIXME fix some hardcoded 1024's below


cdef class vdif_acquisition:
    cdef _ch_vdif_assembler_c.vdif_acquisition *_p

    def __cinit__(self):
        self._p = NULL

    def __init__(self, filename_list, nthreads=10):
        self._p = new _ch_vdif_assembler_c.vdif_acquisition(filename_list, nthreads)

    def __dealloc__(self):
        if self._p != NULL:
            del self._p
            self._p = NULL


cdef class _vdif_assembler:
    cdef _ch_vdif_assembler_c.vdif_assembler *_p

    def __cinit__(self):
        self._p = NULL

    def __init__(self, mask_zeros=False, mask_rails=True, mask_rfi=True, rfi_mask_version=0, trim_frac=1.0e-5, span_frames=False, buf_nt=65536, allow_drops=True):
        self._p = new _ch_vdif_assembler_c.vdif_assembler(mask_zeros, mask_rails, mask_rfi, rfi_mask_version, trim_frac, span_frames, buf_nt, allow_drops)
        self.callbacks = [ ]

    def __dealloc__(self):
        if self._p != NULL:
            del self._p
            self._p = NULL

    def _run(self, vdif_acquisition vs):
        assert self._p != NULL
        assert vs._p != NULL
        return self._p._run(vs._p[0])

    def _advance(self, int nt):
        assert self._p != NULL
        self._p._advance(nt)

    def _allocate(self):
        assert self._p != NULL
        self._p._allocate()

    def _deallocate(self):
        assert self._p != NULL
        self._p._deallocate()

    def _get_buf_t0(self):
        assert self._p != NULL
        return self._p._get_buf_t0()


    def get_mask(self, int64_t t0, int nt):
        assert self._p != NULL

        cdef np.ndarray[np.uint8_t,ndim=3,mode='c'] a = np.zeros((1024,2,nt),dtype=np.uint8)
        self._p.get_mask(&a[0,0,0], t0, nt)
        return a


    def get_complex_data(self, int64_t t0, int nt):
        assert self._p != NULL

        cdef np.ndarray[np.complex64_t,ndim=3,mode='c'] a = np.zeros((1024,2,nt),dtype=np.complex64)
        self._p.get_complex_data(&a[0,0,0], t0, nt)
        return a


    def get_visibilities(self, int64_t t0, int nt):
        assert self._p != NULL

        cdef np.ndarray[np.float32_t,ndim=3,mode='c'] a = np.zeros((1024,4,nt),dtype=np.float32)
        self._p.get_visibilities(&a[0,0,0], t0, nt)
        return a


#def get_rfi_mask(version=0):
#    cdef np.ndarray[np.int,ndim=1,mode='c'] a = np.zeros(1024, dtype=np.int)
#    _ch_vdif_assembler_c.get_rfi_mask(&a[0], version)
