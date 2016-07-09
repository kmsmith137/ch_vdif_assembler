ch_vdif_assembler: a multithreaded, real-time acqusition system for CHIME high-speed data.

This version (v2) is fast enough to keep up with a real-time 
6.4 Gbps network capture, but using assembly language kernels 
seems to be necessary!  (see ch_vdif_assembler_kernels.hpp).

Unforuntately it's likely that as new processing tasks are
added, it may be necessary to add new assembly language kernels
if we want to run in real time.  My plan is to keep doing this
on an ad hoc basis until we hopefully converge to a set of
kernels which is general enough to cover all practical cases.

This version also supports running multiple processing tasks,
which will be run in different threads.  It can also stream
data to disk, or buffer data and save it to disk if one of the
processing tasks sets a trigger.


### INSTALLATION

- You'll need the following prerequisites.  Boost is no longer
  needed but you'll need a newish C++ compiler with C++11 support.
    - python
    - numpy
    - cython
    - libpng
    - libhdf5

- Create a file Makefile.local as described in comments in 
  Makefile.local.example.  You can probably start with one
  of the examples in `site/` (or just symlink one of these
  examples to `./Makefile.local`)

- You should now be able to compile with `make all install`.

- If you want to run some unit tests:
```
./test-timestamp-unwrapper
./test-kernels
./test-downsampled-intensity
./run-vdif-assembler -u
```

- To try it out, run the executable `./run-vdif-assembler`
  and follow the instructions!

- A performance puzzle: the assembly-language-kernel-enabled
  assembler runs 5 times faster on moose (3.5 GHz Haswell-E,
  gcc 4.4.7) than on my desktop (2.4 GHz Haswell, gcc 4.8.3)!
  This seems worth understanding and fixing, but I haven't gotten
  to the bottom of it yet.  In the meantime, if you use moose
  you should be fine, but you may find poor performance on other
  machines.  For a performance test, do ./run-vdif-assembler -t.
  This gives ~13 Gbps on moose and 2.5 Gbps on my desktop, where
  6.4 Gbps is the minimum needed to keep up with a real-time
  network capture.


### WRITING EXTENSIONS

To integrate C++ code with ch_vdif_assembler, you'll want to 
define a subclass of the virtual base class vdif_processor.  For 
a lot more detail, see comments in ch_vdif_assembler.hpp or
an example C++ class in waterfall_plotter.cpp.

There's also a python interface which will probably never be fast 
enough for realtime processing, but may be useful for prototyping 
code using captured data on disk.
