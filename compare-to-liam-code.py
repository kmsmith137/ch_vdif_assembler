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

#rand.seed(123)


if len(sys.argv) != 2:
    print >>sys.stderr, 'Usage: compare-to-liam-code.py <filename_list.txt>'
    print >>sys.stderr, 'You may find the script show-moose-acqusitions.py useful for making file lists'
    sys.exit(2)

filename_list = [ f.strip() for f in open(sys.argv[1]) ]

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

                re = row_data[t_ix,frame_ix,0]
                im = row_data[t_ix,frame_ix,1]

                # emulate masking logic in ch_vdif_assembler
                if (re < -7.5) and (im < -7.5):
                    continue
                
                liam[freq,pol,site_ix] = re + im*1j


class kms_reader(ch_vdif_assembler.processor):
    def __init__(self, sites, output_arr):
        self.sites = sites
        self.output = output_arr

        assert sites.shape == (1024, 2, self.sites.shape[-1])
        assert output_arr.shape == (1024, 2, self.sites.shape[-1])


    def process_chunk(self, t0, nt, data, mask):
        assert data.shape == (1024, 2, nt)
        assert mask.shape == (1024, 2, nt)

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

# 2^17 was large enough to avoid drops in the example I used
assembler = ch_vdif_assembler.assembler(assembler_nt = 2**17)
assembler.register_processor(kms_reader(sites,kms))

s = ch_vdif_assembler.make_file_stream(filename_list)
assembler.run(s)

(fi,pi,ti) = np.where(np.abs(kms-liam) > 0.1)
n = len(fi)

if n == 0:
    print 'Success!'
    sys.exit(0)

# Just print first failure
(fi,pi,ti) = (fi[0],pi[0],ti[0])
fpga = sites[fi,pi,ti]
print 'At fi=%s, pi=%s, ti=%s, fpga=%s:  kms=%s  liam=%s' % (fi, pi, ti, fpga, kms[fi,pi,ti], liam[fi,pi,ti])
