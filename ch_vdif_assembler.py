"""
ch_vdif_assembler: a module for analysis of CHIME high-speed data.

This is the python version, which is too slow for real-time processing, 
but useful for prototyping.  For "production" you probably want to use
the C++ version.

See toy-python-assembler.py for an example.  The basic steps are:

     - Construct an object of type vdif_assembler.
 
     - Define a "processor" class and register it with the assembler.
       This should be a subclass of class vdif_processor, see below.
       See class python_waterfall_plotter below for an example.

     - Register the processor with vdif_assembler.register_processor().
       You can also register one more more C++ processors.  Right now the
       only C++ processor which is exported via cython is cpp_waterfall_plotter
       (see below, or example usage in toy-python-assembler.py).  You can
       run any number of C++ processors, and either 0 or 1 python processors,
       in parallel (each processor will run in a separate threads).

     - Create an object of type vdif_stream, which represents the high-speed
       data stream to be processed.  Since python will probably be too slow
       for real-time use, the stream will probably be one of:
          make_file_stream() - pass a list of vdif filenames
          moose_acquisition() - pass the name of an acqusition on moose, e.g. '41537'

     - Call vdif_assembler.run()

For more documentation, see the C++ source file ch_vdif_assembler.hpp, which 
is better commented!
"""


import os
import re
import sys
import errno
import numpy as np

import ch_vdif_assembler_cython


class constants:
    chime_nfreq = ch_vdif_assembler_cython.chime_nfreq                      # 1024
    timestamps_per_frame = ch_vdif_assembler_cython.timestamps_per_frame    # 2^23 (cadence of noise source)
    num_disks = ch_vdif_assembler_cython.num_disks                          # 10 (moose)


class assembler(object):
    def __init__(self, write_to_disk=False, rbuf_size=constants.num_disks, abuf_size=4, assembler_nt=65536):
        self._assembler = ch_vdif_assembler_cython.assembler(write_to_disk, rbuf_size, abuf_size, assembler_nt)
        self.python_processors = []


    def register_processor(self, p):
        if isinstance(p, ch_vdif_assembler_cython.cpp_processor):
            self._assembler.register_cpp_processor(p)   # register C++ processor (this actually spawns a processing thread)
        elif not isinstance(p, processor):
            raise RuntimeError('Argument to assembler.register_processor() must be either an object of class ch_vdif_assembler.processor, or a C++ processor (e.g. returned by make_waterfall_plotter)')
        self.python_processors.append(p)


    def run(self, stream):
        if not self.python_processors:
            self._assembler.start_async(stream)
            self._assembler.wait_until_end()
            return

        self._assembler.register_python_processor()

        processors_byte_data = [p.byte_data for p in self.python_processors]
        need_byte_data = True in processors_byte_data
        need_complex_data = False in processors_byte_data

        try:
            self._assembler.start_async(stream)

            while True:
                chunk = self._assembler.get_next_python_chunk()
                if chunk is None:
                    break
                if need_byte_data:
                    byte_data = chunk.get_byte_data()
                if need_complex_data:
                    complex_data = chunk.get_data()
                for p in self.python_processors:
                    if p.byte_data:
                        t0, nt, efield = byte_data
                        mask = None
                    else:
                        t0, nt, efield, mask = complex_data
                    p.process_chunk(t0, nt, efield, mask)

            for p in self.python_processors:
                p.finalize()

        finally:
            self._assembler.unregister_python_processor()

        self._assembler.wait_until_end()


########################################  Streams  ########################################


def make_file_stream(filename_list):
    """Returns a stream object from a list of vdif_filenames."""
    return ch_vdif_assembler_cython.make_file_stream(filename_list)


def moose_acquisition(acq_name):
    """
    Helper function which returns a stream object corresponding to one acqusition on moose.

    The acq_name can be either a complete acquisition name, e.g. '20150910T041537Z_chime_beamformed',
    or a uniquely determining substring, e.g. '41537'.  For a list of all acquisitions, you may find
    the script show-moose-acquisitions.py useful.
    """
    filename_list = moose_filename_list(acq_name)
    return ch_vdif_assembler_cython.make_file_stream(filename_list)

def make_network_stream():
    """
    Returns a network stream, but the python processor will probably see a lot of missing
    chunks since it will be too slow to keep up!
    """
    return ch_vdif_assembler_cython.make_network_stream()


def make_simulated_stream(gbps=6.4, nsec=60.0):
    """
    Returns a simualted 6.4 Gpbs stream, but the python processor will probably see a lot of missing
    chunks since it will be too slow to keep up!
    """
    return ch_vdif_assembler_cython.make_simulated_stream(gbps, nsec)


########################################  Processors  ########################################


def cpp_waterfall_plotter(outdir, is_critical=False):
    """This is the C++ version of the waterfall plotter, exported via cython."""

    # Create directory if it doesn't already exist
    try:
        os.makedirs(outdir)
    except OSError, exc:
        if (exc.errno != errno.EEXIST) or not os.path.isdir(outdir):
            raise

    return ch_vdif_assembler_cython.cpp_waterfall_plotter(outdir, is_critical)


class processor(object):
    """
    To define a python processor, you subclass this base class.

    When the assembler runs, it will call process_chunk() with a sequence of chunks, represented
    by a (t0,nt,efield,mask) quadruple.

    Each chunk corresponds to range of timestamps [t0,t0+nt), where t0 is a 64-bit wraparound-free
    timestamp.

    WARNING 1: Usually these ranges will be contiguous between calls, e.g.
        [t0,t0+nt)   [t0+nt,t0+2*nt)   [t0+2*nt,t0+3*nt)   ...
    but the vdif_processor should not assume that this!  If there is a temporary 
    interruption in data stream, then a timestamp gap will appear.

    The 'efield' arg is a shape (nfreq,2,nt) complex array with electric field values, where
    the middle index is polarziation.  Missing data is represented by (0+0j).  The 'mask' arg
    is a shape (nfreq,2,nt) integer array which is 0 for missing data, and 1 for non-missing.

    For subclasses with attribute `byte_data = True`, 'efield' is a shape
    (nfreq,2,nt) byte array, as in the C++ versions. 'mask' is None.

    WARNING 2: Handling missing data is an important aspect of the vdif_processor since it 
    happens all the time.  If a GPU correlator node is down, which is a frequent occurrence, 
    then some frequencies will be "all missing".  There are also routine packet loss events 
    on second-timescales which result in some high-speed samples being flagged as missing data.
    """

    byte_data = False

    def process_chunk(self, t0, nt, efield, mask):
        if mask is None:
            print 'process_chunk called! t0=%s nt=%s efield (%s,%s)' % (t0, nt, efield.dtype, efield.shape)
        else:
            print 'process_chunk called! t0=%s nt=%s efield (%s,%s) mask (%s,%s)' % (t0, nt, efield.dtype, efield.shape, mask.dtype, mask.shape)

    def finalize(self):
        pass


####################################################################################################
#
# Some classes and helper functions for inventory-ing acquisitions on moose.
#
# See also the script show-moose-acquisitions.py


class moose_inventory(object):
    def __init__(self):
        self.topdirs = [ ('/drives/E/%d' % i) for i in xrange(10) ]
        self.subdirs = set()

        for t in self.topdirs:
            for s in os.listdir(t):
                if not s.endswith('_chime_beamformed'):
                    continue
                if len(os.listdir(os.path.join(t,s))) == 0:
                    continue   # skip empty "acquisitions"
                self.subdirs.add(s)

        self.subdirs = sorted(self.subdirs)


    def __iter__(self):
        for s in self.subdirs:
            yield s


    def match(self, name):
        candidate_subdirs = [ s for s in self.subdirs if s.find(name) >= 0 ]

        if len(candidate_subdirs) == 0:
            raise RuntimeError("%s: no matching acquisition found" % name)

        if len(candidate_subdirs) > 1:
            raise RuntimeError("%s: multiple matching acquisitions found: %s" % (name, sorted(candidate_subdirs)))

        return candidate_subdirs.pop()


    def get_filename_list(self, subdir):
        subdir = self.match(subdir)   # disambiguate
        filename_hash = { }           # hash integer index -> filename
    
        for t in self.topdirs:
            d = os.path.join(t, subdir)
            if not os.path.exists(d):
                continue

            for f in os.listdir(d):
                df = os.path.join(d,f)

                m = re.match(r'(\d+)\.dat', f)
                if not m:
                    print >>sys.stderr, "Warning: stray file '%s' found in acquisition directory" % df
                    continue

                i = int(m.group(1))
                if filename_hash.has_key(i):
                    print >>sys.stderr, "Warning: redundant acquisition files exist? (%s and %s)" % (filename_hash[i], df)
                    continue

                filename_hash[i] = df

        if len(filename_hash) == 0:
            return [ ]

        filename_list = [ ]
        nmissing = 0

        for i in xrange(max(filename_hash.keys()) + 1):
            if not filename_hash.has_key(i):
                nmissing += 1
                continue

            filename_list.append(filename_hash.pop(i))

        if nmissing > 0:
            print >>sys.stderr, 'Warning: in acquisition %s, %d files were missing' % (subdir, nmissing)

        return filename_list


def moose_filename_list(name):
    mi = moose_inventory()
    subdir = mi.match(name)
    return mi.get_filename_list(subdir)


####################################################################################################
#
# python waterfall plotter
#
# I originally prototyped this in python, then re-implemented in C++ because the python version was
# too slow.  I left the python version intact since it's a good reference for implementing processors.
#
# This means that there are currently two waterfall plotters, the C++ version which is faster, or
# and the python version which is simpler code.


try:
    import PIL.Image
except:
    print >>sys.stderr, 'ch_vdif_assembler: import PIL.Image failed; everything will still work except python_waterfall_plotter'


class python_waterfall_plotter(processor):
    def __init__(self, outdir):
        self.outdir = outdir
        self.imgdir = os.path.join(outdir, 'img')
        self.nt_bins = 1024
        self.curr_frame = None

        # Create output dirs, but don't throw an exception if they already exist
        try:
            os.makedirs(self.imgdir)
        except OSError, exc:
            if (exc.errno != errno.EEXIST) or not os.path.isdir(self.imgdir):
                raise

    
    def process_chunk(self, t0, nt, efield, mask):
        assert efield.shape == (constants.chime_nfreq, 2, nt)
        assert efield.dtype == np.complex64

        # Square electric field to get visibilities
        vis = efield.real**2 + efield.imag**2

        frame = t0 // constants.timestamps_per_frame

        if self.curr_frame is None:
            self.acc_num = np.zeros((constants.chime_nfreq, 2, self.nt_bins))
            self.acc_den = np.zeros((constants.chime_nfreq, 2, self.nt_bins))
            self.curr_frame = frame

        elif self.curr_frame != frame:
            self.write_images()
            self.acc_num = np.zeros((constants.chime_nfreq, 2, self.nt_bins))
            self.acc_den = np.zeros((constants.chime_nfreq, 2, self.nt_bins))
            self.curr_frame = frame

        # Number of timestamps per bin
        tt = constants.timestamps_per_frame // self.nt_bins

        # These assumptions simplify the code, but could be removed if necessary
        assert constants.timestamps_per_frame % self.nt_bins == 0
        assert t0 % tt == 0
        assert nt % tt == 0

        b0 = (t0 // tt) % self.nt_bins   # first waterfall bin in assembled chunk
        nb = (nt // tt)                  # number of waterfall bins in assembled chunk

        for b in xrange(nb):
            (i, j) = (b*tt, (b+1)*tt)
            self.acc_num[:,:,b+b0] += np.sum(vis[:,:,i:j], axis=2)
            self.acc_den[:,:,b+b0] += np.sum(mask[:,:,i:j], axis=2)


    def finalize(self):
        if self.curr_frame is not None:
            self.write_images()


    def write_images(self):
        self.acc_num = self.subtract_medians(self.acc_num, self.acc_den)

        thumbnail_num = self.downgrade(self.acc_num, 256, 256)
        thumbnail_den = self.downgrade(self.acc_den, 256, 256)

        for pol in [0,1]:
            filename1 = '%s/full_pol%d_frame%d.png' % (self.imgdir, pol, self.curr_frame)
            filename2 = '%s/thumbnail_pol%d_frame%d.png' % (self.imgdir, pol, self.curr_frame)

            self.write_image(filename1, self.acc_num[:,pol,:], self.acc_den[:,pol,:])
            self.write_image(filename2, thumbnail_num[:,pol,:], thumbnail_den[:,pol,:])

        self.curr_frame = None


    def subtract_medians(self, num, den):
        """This helper routine applies median subtraction along the last array axis."""

        assert num.shape == den.shape

        nt = num.shape[-1]
        num_rs = np.reshape(num, (-1,nt))
        den_rs = np.reshape(den, (-1,nt))

        mask_rs = (den_rs > 99.5)
        vals_rs = num_rs / np.maximum(den_rs, 1.0)    # avoid division by zero
        
        nrs = num_rs.shape[0]
        med = np.zeros(nrs)

        for i in xrange(nrs):
            t = np.extract(mask_rs[i,:], vals_rs[i,:])
            if len(t) >= 1:
                med[i] = np.median(t)

        ret_rs = num_rs - (med[:,np.newaxis] * den_rs)
        return np.reshape(ret_rs, num.shape)


    def downgrade(self, img, new_nx, new_ny):
        """
        This is called on a "polarized" image of shape (old_nx, 2, old_ny).
        It returns a downgraded image of shape (new_nx, 2, new_ny).
        """

        assert img.shape == (img.shape[0], 2, img.shape[-1])
        (old_nx, old_ny) = (img.shape[0], img.shape[-1])

        assert old_nx % new_nx == 0
        assert old_ny % new_ny == 0

        ret = np.reshape(img, (new_nx, old_nx//new_nx, 2, new_ny, old_ny//new_ny))
        ret = np.sum(ret, axis=4)
        ret = np.sum(ret, axis=1)
        return ret


    def write_image(self, filename, num, den):
        assert (num.shape == den.shape) and (num.ndim == 2)

        rgb = np.zeros((num.shape[0], num.shape[1] ,3), dtype=np.uint8)
        vvar = self.compute_clipped_variance(num, den)

        if vvar > 0.0:
            # plot saturates at 2.5 sigma
            dv = 2.5 * max(vvar,0.0)**0.5

            vis = num / np.maximum(den, 1.0)
            mask = (den > 99.5)

            vis = 128*vis/dv + 128.
            vis = np.array(vis, dtype=np.int)
            vis = np.maximum(vis, 0)
            vis = np.minimum(vis, 255)
            
            rgb[:,:,0] = np.where(mask, vis, 0)
            rgb[:,:,2] = np.where(mask, 255-vis, 0)
        
        #
        # PIL has weird plotting conventions: the x-axis of the plot is array index 1,
        # and the y-axis of the plot is array index 0 with a sign flip.  This turns out
        # to be exactly what we want, since we have time in array index 1, and we have
        # frequency channel (with the CHIME sign flip) in array index 0.  So we can just
        # pass the RGB array without any transposes or sign flips.
        #
        img = PIL.Image.fromarray(rgb)
        img.save(filename)
        print >>sys.stderr, 'wrote %s' % filename


    def compute_clipped_variance(self, num, den):
        assert (num.shape == den.shape) and (num.ndim == 2)
        
        mask = (den > 99.5)
        vvar = self.compute_variance_outside_mask(num, den, mask)

        if vvar <= 0.0:
            return 0.0
    
        for n in xrange(3):
            # extend mask by clipping at 2 sigma
            mask = np.logical_and(mask, num*num < 4*vvar*den*den)

            # recompute variance with extended mask
            vvar2 = self.compute_variance_outside_mask(num, den, mask)

            if vvar2 <= 0.0:
                return vvar

            vvar = vvar2
        
        return vvar


    def compute_variance_outside_mask(self, num, den, mask):
        msum = np.sum(mask)
        if msum <= 0.5:
            return 0.0

        # assume denominator is > 0 outside mask (FIXME assert this?)
        x = mask * (num*num) / np.maximum(den*den, 1.0)
        return np.sum(x) / msum
