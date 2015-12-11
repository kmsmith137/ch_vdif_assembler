from libc.stdint cimport int32_t, int64_t
from libcpp.vector cimport vector
from libcpp cimport bool

import numpy as np
cimport numpy as np

cimport ch_vdif_assembler_pxd


##############################################  Constants  #########################################


chime_nfreq = ch_vdif_assembler_pxd.chime_nfreq
timestamps_per_frame = ch_vdif_assembler_pxd.timestamps_per_frame
num_disks = ch_vdif_assembler_pxd.num_disks


##############################################  Streams  ###########################################


cdef class stream:
    cdef ch_vdif_assembler_pxd.cython_stream *p

    def __cinit__(self):
        self.p = NULL

    def __dealloc__(self):
        if self.p != NULL:
            del self.p
            self.p = NULL


def make_file_stream(filename_list):
    ret = stream()
    ret.p = ch_vdif_assembler_pxd.cython_file_stream(filename_list)
    return ret

def make_network_stream():
    ret = stream()
    ret.p = ch_vdif_assembler_pxd.cython_network_stream()
    return ret

def make_simulated_stream(gbps, nsec):
    ret = stream()
    ret.p = ch_vdif_assembler_pxd.cython_simulated_stream(gbps, nsec)
    return ret


############################################  cpp_processor  #######################################


cdef class cpp_processor:
    cdef ch_vdif_assembler_pxd.cpp_processor *p

    def __cinit__(self):
        self.p = NULL

    def __dealloc__(self):
        if self.p != NULL:
            del self.p
            self.p = NULL


def cpp_waterfall_plotter(outdir, is_critical):
    ret = cpp_processor()
    ret.p = ch_vdif_assembler_pxd.cpp_waterfall_plotter(outdir, is_critical)
    return ret


###########################################  assembled_chunk  ######################################


cdef class assembled_chunk:
    cdef ch_vdif_assembler_pxd.cython_assembled_chunk *p

    def __cinit__(self):
        self.p = NULL

    def __dealloc__(self):
        if self.p != NULL:
            del self.p
            self.p = NULL

    def get_data(self):
        if self.p == NULL:
            return None

        t0 = self.p[0].t0
        nt = self.p[0].nt

        cdef np.ndarray[np.complex64_t,ndim=3,mode='c'] efield = np.zeros((chime_nfreq,2,nt),dtype=np.complex64)
        cdef np.ndarray[np.int32_t,ndim=3,mode='c'] mask = np.zeros((chime_nfreq,2,nt),dtype=np.int32)

        self.p[0].fill_efield(&efield[0,0,0], &mask[0,0,0])
        return (t0, nt, efield, mask)


##############################################  Assembler  #########################################


cdef class assembler:
    cdef ch_vdif_assembler_pxd.cython_assembler *p

    def __cinit__(self, bool write_to_disk, int rbuf_size, int abuf_size, int assembler_nt):
        self.p = new ch_vdif_assembler_pxd.cython_assembler(write_to_disk, rbuf_size, abuf_size, assembler_nt)

    def __dealloc__(self):
        if self.p != NULL:
            del self.p
            self.p = NULL

    def register_cpp_processor(self, cpp_processor processor):
        self.p[0].register_cpp_processor(processor.p)

    def register_python_processor(self):
        self.p[0].register_python_processor()

    def get_next_python_chunk(self):
        cdef ch_vdif_assembler_pxd.cython_assembled_chunk *a = self.p[0].get_next_python_chunk()

        if a == NULL:
            return None
 
        ret = assembled_chunk()
        ret.p = a
        return ret

    def unregister_python_processor(self):
        self.p[0].unregister_python_processor()

    def start_async(self, stream s):
        self.p[0].start_async(s.p)

    def wait_until_end(self):
        self.p[0].wait_until_end()
