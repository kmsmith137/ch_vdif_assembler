#!/usr/bin/env python

import sys


if len(sys.argv) != 2:
    print >>sys.stderr, 'Usage: toy-python-assembler.py <acq_name>'
    print >>sys.stderr, '    acq_name can be a complete acquisition name, e.g. 20150910T041537Z_chime_beamformed'
    print >>sys.stderr, '    or a uniquely determining substring, e.g. 41537'
    print >>sys.stderr
    print >>sys.stderr, 'This script will run both the C++ and python versions of the waterfall_plotter.'
    print >>sys.stderr, 'Output files will appear in the waterfall_cpp and waterfall_python directories.'
    print >>sys.stderr
    print >>sys.stderr, '(If you compare the running time with run-vdif-assembler -w, you\'ll see how'
    print >>sys.stderr, ' much faster the C++ code is than python.  But if you compare the waterfall plotter'
    print >>sys.stderr, ' code in waterfall_plotter.cpp and ch_vdif_assembler.py, you\'ll see how much'
    print >>sys.stderr, ' more convenient python is than C++!)'
    print >>sys.stderr
    print >>sys.stderr, 'Suggested usage:'
    print >>sys.stderr, '   toy-python-assembler.py 41537'
    print >>sys.stderr, '   index-vdif-waterfalls.py waterfall_cpp_41537'
    print >>sys.stderr, '   index-vdif-waterfalls.py waterfall_python_41537'
    print >>sys.stderr, 'The index-vdif-waterfalls script makes an HTML summary page with clickable thumbnails'

    sys.exit(2)


acq_name = sys.argv[1]

import ch_vdif_assembler

# Construct the stream object
# ch_moose_acquisition() is a helper function which finds files on moose and returns a stream object to read them
stream = ch_vdif_assembler.moose_acquisition(acq_name)

# Construct an assembler object (this spawns some worker threads)
assembler = ch_vdif_assembler.assembler()

# Construct a C++ processor exported with cython (the only one currently exported is the waterfall_plotter)
outdir = 'waterfall_cpp_' + acq_name
p = ch_vdif_assembler.cpp_waterfall_plotter(outdir, is_critical=True)
assembler.register_processor(p)

# Construct a python processor (the reference waterfall plotter in ch_vdif_assembler.py)
outdir = 'waterfall_python_' + acq_name
p = ch_vdif_assembler.python_waterfall_plotter(outdir)
assembler.register_processor(p)

# In just a few minutes, you will have your beautiful waterfall plots, you lucky devil
assembler.run(stream)
