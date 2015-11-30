"""
ch_vdif_assembler: a module for analysis of CHIME high-speed data.

This is the python version, which is too slow for real-time processing, 
but useful for prototyping.  For "production" you probably want to use
the C++ version.

To use this library, you:

     - Construct an object of type vdif_assembler, which implements some
       generic logic for handling the above nuisance issues.
 
     - Define a "processor" class and register it with the assembler.
       This could be as minimal as defining a class with the following 
       methods defined:

          self.process_data(asm, t0, nt)
          self.finalize()    # optional

     - Create an object of type vdif_stream, which represents the high-speed
       data stream to be processed.  Right now the only option is vdif_acqu of vdif_stream
       are vdif_file (representing a single file) or vdif_acquisitio
       (representing multiple files).

       To get a list of filenames in an acqusition on moose, you may find
       the utility show-moose-acquisitions.py useful.

For more information, see the waterfall_plotter example below, or see
the C++ source files, which have more documentation.
"""

import os
import re
import sys
import time
import errno
import numpy as np

# cython extension classes
import ch_vdif_assembler_c

# note: constructor syntax is vdif_acquistion(filename_list, nthreads)
# where nthreads=10 is a good choice on moose
from ch_vdif_assembler_c import vdif_acquisition


class vdif_assembler(ch_vdif_assembler_c._vdif_assembler):
    """
    Note: the following methods are defined in the cython base class:
    
       self.get_mask(t0, nt)          -> returns shape (nfreq,2,nt) uint8-array
       self.get_complex_data(t0, nt)  -> returns shape (nfreq,2,nt) complex array
       self.get_visibilities(t0, nt)  -> returns shape (nfreq,4,nt) real array (middle index is xauto, yauto, recross, imcross)
    """

    def __init__(self, mask_zeros=False, mask_rails=True, mask_rfi=True, rfi_mask_version=0, trim_frac=1.0e-5, span_frames=False, buf_nt=65536, allow_drops=True):
        ch_vdif_assembler_c._vdif_assembler.__init__(self, mask_zeros, mask_rails, mask_rfi, rfi_mask_version, trim_frac, span_frames, buf_nt, allow_drops)
        self.processors = [ ]


    def register_processor(self, x):
        self.processors.append(x)


    def run(self, vdif_stream):
        #
        # Mimic vdif_assembler::run() from C++ code
        #

        sys.stderr.write('vdif_assembler: start!\n')
        self._allocate()
        
        while True:
            nt = self._run(vdif_stream)

            if nt == 0:
                break

            for x in self.processors:
                x.process_data(self, self._get_buf_t0(), nt)

            self._advance(nt)

        for x in self.processors:
            if hasattr(x, 'finalize'):
                x.finalize()

        self._deallocate()
        sys.stderr.write('vdif_assembler: end\n')


####################################################################################################


class moose_inventory:
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


def moose_acquisition(name, nthreads=10):
    mi = moose_inventory()
    subdir = mi.match(name)
    filename_list = mi.get_filename_list(subdir)

    return ch_vdif_assembler_c.vdif_acquisition(filename_list, nthreads)


####################################################################################################
#
# Everything below this needs cleanup; ignore for now


class vis_frame_analyzer:
    """
    Analyzes 21-sec (2^23 FPGA sample) "frames"

    This is a vdif_assembler processor, i.e. it defines the method
       self.process_data(asm, t0, nt)

    This routine does frame accumulation, and defers the processing to the
    following method, which must be defined by the subclass:
        self.process_frame(vis, counts, s)

    where @vis, @counts are real arrays of shape (nfreq_bins, 2, nt)
    and @s is the frame index (starting from an arbitrary value; subclass
    shouldn't assume that s=0 in the first frame).

    The subclass may also define
        self.finalize()

    Note that we only process auto correlations for now.
    """


    def __init__(self, nfreq_bins, nt_bins, trim_frac=0.001):
        assert (nfreq_bins > 0) and (1024 % nfreq_bins == 0)
        assert (nt_bins > 0)
        assert (trim_frac > 0.0) and (trim_frac < 0.5)

        self.nfreq_bins = nfreq_bins
        self.nt_bins = nt_bins
        self.trim_frac = trim_frac
        
        self.vis_accum = np.zeros((nfreq_bins,2,nt_bins))
        self.count_accum = np.zeros((nfreq_bins,2,nt_bins), dtype=np.int)
        self.scurr = 0
        self.empty = True

        # range of unpadded t-values in each frame
        self.up0 = int(0.5 * trim_frac * 2**23)
        self.up1 = 2**23 - self.up0
        assert self.up1 - self.up0 >= nt_bins


    def _rebin_freqs(self, arr):
        """
        A little helper routine to rebin in frequency.
            Input @arr is an array of shape (1024, ...)
            Returns an array of shape (self.nfreq_bins, ...).
        """

        assert arr.shape[0] == 1024

        tmp_shape = (self.nfreq_bins, -1) + arr.shape[1:]
        tmp = np.reshape(arr, tmp_shape)
        return np.sum(tmp, axis=1)


    def process_data(self, asm, t0, nt):
        data = asm.get_complex_data(t0, nt)     # shape (nfreq, 2, nt) complex
        data = data.real**2 + data.imag**2      # complex -> real (electric field -> visibility)
        data = self._rebin_freqs(data)

        # note conversion to 'int' here, to avoid 8-bit overflow
        counts = np.array(asm.get_mask(t0,nt), dtype=np.int)   # shape (nfreq, 2, nt)
        counts = self._rebin_freqs(counts)
        assert data.shape == counts.shape == (self.nfreq_bins, 2, nt)

        #
        # Handle case where data spans a frame boundary
        #
        # This loop is correct for a data block which spans N frames,
        # although N=1 and N=2 are usually the only possibilities.
        #
        while nt > 0:
            t1 = ((t0 // 2**23) + 1) * 2**23    # start of next frame
            nt_frame = min(t1-t0, nt)
            self.process_data2(data[:,:,:nt_frame], counts[:,:,:nt_frame], t0, nt_frame)

            data = data[:,:,nt_frame:]
            counts = counts[:,:,nt_frame:]
            t0 += nt_frame
            nt -= nt_frame


    def process_data2(self, vis, counts, t0, nt):
        """
        Helper for process_data(): process one chunk of data, which has been binned
        in frequency but not time, and guaranteed not to span a frame boundary.
        """

        assert vis.shape == counts.shape == (self.nfreq_bins, 2, nt)

        s0 = t0 // 2**23
        u0 = t0 % 2**23
        s1 = (t0 + nt - 1) // 2**23
        assert s0 == s1

        up0 = self.up0
        up1 = self.up1
        nt_bins = self.nt_bins

        if self.empty:
            print >>sys.stderr, 'Frame %s: start' % s0
            self.scurr = s0
            self.empty = False
        elif s0 < self.scurr:
            raise RuntimeError("internal error in vis_frame_analyzer: timestamps running backwards")
        elif s0 > self.scurr:
            print >>sys.stderr, 'Frame %s: end' % self.scurr
            self.vis_accum /= np.maximum(self.count_accum,1)
            self.process_frame(self.vis_accum, self.count_accum, self.scurr)
            print >>sys.stderr, 'Frame %s: start' % s0
            self.vis_accum[:,:,:] = 0
            self.count_accum[:,:,:] = 0
            self.scurr = s0

        #
        # Remaining code in this routine does the accumulation.
        # First we compute the range of t-bins corresponding to this block.
        # The map from t -> t_bin is given by b = ((t-up0) * nt_bins) // (up1-up0).
        #
        b0 = ((u0 - up0) * nt_bins) // (up1-up0)
        b1 = ((u0+nt-1 - up0) * nt_bins) // (up1-up0) + 1
        (b0, b1) = (max(b0,0), min(b1,nt_bins))    # needed if trim_frac > 0

        for b in xrange(b0, b1):
            # The map from t_bin -> t is given by t = b * (up1-up0) // nt_bins + up0
            i = (b * (up1-up0)) // nt_bins + up0-u0
            j = ((b+1) * (up1-up0)) // nt_bins + up0-u0
            (i, j) = (max(i,0), min(j,nt))

            self.vis_accum[:,:,b] += np.sum(vis[:,:,i:j], axis=2)
            self.count_accum[:,:,b] += np.sum(counts[:,:,i:j], axis=2)


class vis_waterfall_plotter(vis_frame_analyzer):
    """
    A subclass of vis_frame_analyzer which makes waterfall plots of high-res 
    data.  This is generally useful for understanding what's happening in a 
    capture.  The waterfall plots are linked into a top-level file 
    ${outdir}/index.html with clickable thumbnails.
    """

    def __init__(self, outdir, trim_frac=0.001):
        vis_frame_analyzer.__init__(self, 1024, 1024, trim_frac=trim_frac)

        self.outdir = outdir
        self.imgdir = os.path.join(outdir, 'img')
        # self.rfi_mask = get_rfi_mask()
        self.slist = [ ]

        for d in [ self.outdir, self.imgdir ]:
            try:
                os.makedirs(d)
            except OSError, exc:
                if (exc.errno != errno.EEXIST) or not os.path.isdir(d):
                    raise

        filename = os.path.join(self.outdir, 'index.html')
        self.f_html = open(filename, 'w')
        print >>self.f_html, '<html><head><title>run 2132</title></head>'
        print >>self.f_html, '<body>'
        print >>self.f_html, '<table cellspacing="10">'


    def process_frame(self, vis, counts, s):
        assert vis.shape == counts.shape == (1024, 2, 1024)
        self.slist.append(s)

        for pol in [0,1]:
            for is_thumbnail in [False,True]:
                self._write_image(vis, counts, s, pol, is_thumbnail)

        print >>self.f_html, '<tr> <td> %s' % s

        for pol in [0,1]:
            full_filename = self._img_basename(s, pol, is_thumbnail=False, rfi_flag=False)
            thumbnail_filename = self._img_basename(s, pol, is_thumbnail=True, rfi_flag=False)
            print >>self.f_html, '     <td> <a href="img/%s"><img src="img/%s"></a>' % (full_filename, thumbnail_filename)


    def _make_image(self, vis, counts, s, pol, is_thumbnail):
        """
        Helper for write_image(): returns RGB array with shape (n,n,3) and dtype uint8.

        The rfi mask is not applied.  To apply it to the image, do e.g.
           rgb[:,:,:] *= self.rfi_mask[:,np.newaxis,np.newaxis]
        """

        vis = vis[:,pol,:]
        counts = counts[:,pol,:]
        mask = np.array(counts > 100, dtype=np.int)

        # Per-frequency median subtraction outside mask
        for ifreq in xrange(1024):
            v = np.extract(mask[ifreq,:], vis[ifreq,:])
            if len(v) >= 2:
                vis[ifreq,:] -= np.median(v)


        if is_thumbnail:
            # (1024,1024) -> (256,256)
            vis = np.reshape(counts*vis, (256, 4, 256, 4))
            vis = np.sum(vis, axis=3)
            vis = np.sum(vis, axis=1)

            counts = np.reshape(counts, (256, 4, 256, 4))
            counts = np.sum(counts, axis=3)
            counts = np.sum(counts, axis=1)

            vis /= np.maximum(counts,1)
            mask = np.array(counts > 100, dtype=np.int)


        #
        # The next bit of code computes vvar, a clipped scalar
        # variance for the entire 2d visibility array.
        #

        vvar = np.sum(mask * vis**2) / max(np.sum(mask), 1)

        if vvar <= 0:
            print >>sys.stderr, 'Frame %d: warning: not enough data to make plot' % s
            return np.zeros((vis.shape[0], vis.shape[1] ,3), dtype=np.uint8)

        # clip mask
        cmask = np.copy(mask)

        for n in xrange(5):
            # try to compute variance with smaller clip mask...
            cmask = np.logical_and(cmask, np.abs(vis) < 2*vvar**0.5)
            vvar2 = np.sum(cmask * vis**2) / max(np.sum(cmask), 1)
        
            if vvar2 <= 0:
                # ... if we fail, that's ok
                print >>sys.stderr, 'Frame %d: warning: clipping logic went awry' % s                
                break
 
            # ... if we succeed, update variance and iterate
            vvar = vvar2

        #
        # Now we fill the RGB image, using vvar to set the plot scale
        #

        dv = 2.5 * vvar**0.5   # plot saturates at 2.5 sigma
        vis = 128*vis/dv + 128.
        vis = np.array(vis, dtype=np.int)
        vis = np.maximum(vis, 0)
        vis = np.minimum(vis, 255)

        rgb = np.zeros((vis.shape[0], vis.shape[1] ,3), dtype=np.uint8)
        rgb[:,:,0] = np.where(mask, vis, 0)
        rgb[:,:,2] = np.where(mask, 255-vis, 0)

        #
        # PIL has weird plotting conventions: the x-axis of the plot is array index 1,
        # and the y-axis of the plot is array index 0 with a sign flip.  This turns out
        # to be exactly what we want, since we have time in array index 1, and we have
        # frequency channel (with the CHIME sign flip) in array index 0.  So we can just
        # return the RGB array without any transposes or sign flips.
        #
        return rgb


    def _img_basename(self, s, pol, is_thumbnail, rfi_flag):
        prefix = 'thumbnail' if is_thumbnail else 'full'
        suffix = 'm' if rfi_flag else ''
        basename = '%s_pol%d_frame%d%s.png' % (prefix, pol, s, suffix)
        return basename


    def _write_image(self, vis, counts, s, pol, is_thumbnail):
        import PIL.Image

        rgb = self._make_image(vis, counts, s, pol, is_thumbnail)

        filename = os.path.join(self.imgdir, self._img_basename(s, pol, is_thumbnail, rfi_flag=False))
        img = PIL.Image.fromarray(rgb)
        img.save(filename)
        print >>sys.stderr, 'wrote %s' % filename
        
        if is_thumbnail:
            return

        #rgb[:,:,:] *= self.rfi_mask[:,np.newaxis,np.newaxis]
        #filename = os.path.join(self.imgdir, self._img_basename(s, pol, is_thumbnail, rfi_flag=True))
        #img = PIL.Image.fromarray(rgb)
        #img.save(filename)
        #print >>sys.stderr, 'wrote %s' % filename


    def finalize(self):
        print >>self.f_html, '</table></body></html>'

