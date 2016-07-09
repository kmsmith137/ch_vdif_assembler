- v3:
  - downsampled_intensity helper class, for vdif_processors which process downsampled intensity data (rather than high-speed baseband)
  - intensity_beam processor thread, for writing downsampled intensity data to disk in CHIMEFRB hdf5 format (`run-vdif-assembler -i`)
  - fix some minor bitrotting (e.g. `/drives/E` -> `/drives/G` on moose
  - This update was focused on the shortest path to fixing the "two-value" bug which was affecting the CHIMEFRB search.
    As a result there is a big loose end: cythonizing the downsampled_intensity class and using it consistently in more places
    (e.g. both Python and C++ waterfall plotters should use it)

- v2: effective initial release

- v1: useless crap
