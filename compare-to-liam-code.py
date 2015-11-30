#!/usr/bin/env python
#
# unit-tests the ch_vdif_assembler, by comparing it to Liam's BeamReader code


import os
import sys
import numpy as np
import numpy.random as rand
import ch_vdif_assembler

# Liam code
import ReadBeamform


if len(sys.argv) != 2:
    print >>sys.stderr, 'Usage: compare-to-liam-code.py <acq_name>'
    print >>sys.stderr, '    acq_name can be a complete acquisition name, e.g. 20150910T041537Z_chime_beamformed'
    print >>sys.stderr, '    or a uniquely determining substring, e.g. 41537'
    sys.exit(2)

acq_name = sys.argv[1]

# rand.seed(123)

inv = ch_vdif_assembler.moose_inventory()
filename_list = inv.get_filename_list(acq_name)


def br_read_file(filename):
    """
    Read vdif file using Liam code and return (header,data), where

        pol = header[:,0]       # pol (0 or 1)
        link_id = header[:,1]
        slot_id = header[:,2]         # node number
        frame_ix = header[:,3]
        times_j2000 = header[:,4]
        fpga_count = header[:,5]

    The frequency index is (slot_id + 16 * link_id + 128 * frame_number).

    The data is an array of shape (N, 625, 8, 2) where the first index is 
    vdif row, second is fpga count, third is frame_number, last index is 
    real/complex.
    """

    sys.stderr.write('reading %s...' % filename)
    RB = ReadBeamform.ReadBeamform()
    (header, data) = RB.read_file_dat(filename)
    sys.stderr.write('done\n')

    data = np.reshape(data, (data.shape[0], 625, 8, 2))
    return (header, data)


print >>sys.stderr, 'Reading first and last file, in order to determine FPGA count range'

# note this only works if no wraparound
fpga_min = br_read_file(filename_list[0])[0][0,5]
fpga_max = br_read_file(filename_list[-1])[0][-1,5]

# widen range slightly, to test corner cases
fpga_min -= 10000
fpga_max += 10000

if fpga_min > fpga_max:
    print >>sys.stderr, "This acquisition contains FPGA count wraparound, which isn't supported in this unit test! (FIXME)"
    sys.exit(1)

# make a random set of "comparison sites", i.e. (freq, pol, time) triples where the two codes are to be compared
nt = 10000
sites = np.zeros((1024, 2, nt), dtype=np.int)

for i in xrange(1024):
    for j in xrange(2):
        # generate random set of sorted distinct integers between fpga_min and fpga_max
        v = rand.randint(fpga_min, fpga_max-nt, size=nt)
        v = np.sort(v)
        v += np.arange(nt)
        sites[i,j,:] = v

# these arrays will store complex electric field values at the sites, as computed by the two codes
kms = np.zeros((1024, 2, nt), dtype=np.complex)
liam = np.zeros((1024, 2, nt), dtype=np.complex)

# initialize to some nominal values which will mean "no data"
kms[:,:,:] = 1.0e10 + 1.0e10j
liam[:,:,:] = 1.0e10 + 1.0e10j

print >>sys.stderr, "\nReading files with Liam's BeamReader code"

for filename in filename_list:
    (file_header, file_data) = br_read_file(filename)

    for (row_header, row_data) in zip(file_header, file_data):
        assert row_header.shape == (6,)
        assert row_data.shape == (625, 8, 2)

        pol = row_header[0]
        link_id = row_header[1]
        slot_id = row_header[2]
        fpga0 = row_header[5]

        for frame_ix in xrange(8):
            freq = slot_id + 16 * link_id + 128 * frame_ix
            
            # determine site range for this row
            i = np.searchsorted(sites[freq,pol,:], fpga0)
            j = np.searchsorted(sites[freq,pol,:], fpga0+625)

            for site_ix in xrange(i,j):
                fpga_site = sites[freq,pol,site_ix]
                t_ix = fpga_site - fpga0
                assert 0 <= t_ix < 625

                liam[freq,pol,site_ix] = row_data[t_ix,frame_ix,0] + row_data[t_ix,frame_ix,1]*1j


class kms_reader:
    def __init__(self, sites, output_arr):
        self.sites = sites
        self.output = output_arr

        assert sites.shape == (1024, 2, self.sites.shape[-1])
        assert output_arr.shape == (1024, 2, self.sites.shape[-1])


    def process_data(self, asm, t0, nt):
        mask = asm.get_mask(t0, nt)
        assert mask.shape == (1024, 2, nt)

        data = asm.get_complex_data(t0, nt)
        assert data.shape == (1024, 2, nt)

        for freq in xrange(1024):
            for pol in xrange(2):
                # determine site range for this row
                i = np.searchsorted(self.sites[freq,pol,:], t0)
                j = np.searchsorted(self.sites[freq,pol,:], t0+nt)

                for site_ix in xrange(i,j):
                    fpga_site = self.sites[freq,pol,site_ix]
                    t_ix = fpga_site - t0
                    assert 0 <= t_ix < nt

                    if mask[freq,pol,t_ix]:
                        self.output[freq,pol,site_ix] = data[freq,pol,t_ix]
                    else:
                        assert data[freq,pol,t_ix] == 0+0j


print >>sys.stderr, "\nReading files with ch_vdif_assembler"

# an assembler which keeps everything
asm = ch_vdif_assembler.vdif_assembler(mask_zeros=False, mask_rails=False, mask_rfi=False, trim_frac=0, buf_nt=2**17, allow_drops=False)

asm.register_processor(kms_reader(sites,kms))

acq = ch_vdif_assembler.vdif_acquisition(filename_list)
asm.run(acq)

(fi,pi,ti) = np.where(np.abs(kms-liam) > 0.1)
n = len(fi)

if n == 0:
    print 'Success!'
    sys.exit(0)

# Just print first failure
(fi,pi,ti) = (fi[0],pi[0],ti[0])
print 'At fi=%s, pi=%s, ti=%s:  kms=%s  liam=%s' % (fi, pi, ti, kms[fi,pi,ti], liam[fi,pi,ti])
