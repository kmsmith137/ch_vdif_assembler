#!/usr/bin/env python
#
# Usage:
#
# With no command-line arguments,
#   show-moose-acquisitions.py
#
# will show a list of all acquisitions on moose.  With command-line arguments,
#   show-moose-acquisitions.py [acq_name1] [acq_name2] ...
#
# will print a list of filenames (one filename per line) contained in the
# acquisitions.  Each acq_name can be a complete acquisition name, e.g.
# 20150910T041537Z_chime_beamformed, or a uniquely determining substring,
# e.g. 41537.  (I often redirect the output to a file, in order to supply
# a filename list to another program.)
#


import sys
import ch_vdif_assembler

mi = ch_vdif_assembler.moose_inventory()

if len(sys.argv) > 1:
    filename_list = [ ]

    for acq_name in sys.argv[1:]:
        subdir = mi.match(acq_name)
        filename_list.extend(mi.get_filename_list(subdir))

    # output in separate loop so warnings will be shown separately
    for filename in filename_list:
        print filename

else:
    # list of (subdir, nfiles) pairs
    acquisition_list = [ ]

    for subdir in mi:
        n = len(mi.get_filename_list(subdir))
        acquisition_list.append((subdir,n))

    # output in separate loop so warnings will be shown separately
    for (subdir,n) in acquisition_list:
        print '%s: %d files' % (subdir, n)
