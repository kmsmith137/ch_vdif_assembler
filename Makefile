# Makefile.local must define the following variables
#   BINDIR   install dir for executables
#   LIBDIR   install dir for C++ libraries
#   INCDIR   install dir for C++ headers
#   PYDIR    install dir for python/cython modules
#   CPP      g++ compiler command line (must support c++0x, and -lpthread or -pthread if necessary)
#
# See Makefile.local.example for an example.

include Makefile.local

ifndef CPP
$(error Fatal: Makefile.local must define CPP variable)
endif

ifndef LIBDIR
$(error Fatal: Makefile.local must define LIBDIR variable)
endif

ifndef INCDIR
$(error Fatal: Makefile.local must define INCDIR variable)
endif

ifndef BINDIR
$(error Fatal: Makefile.local must define BINDIR variable)
endif

ifndef PYDIR
$(error Fatal: Makefile.local must define PYDIR variable)
endif


####################################################################################################


INCFILES=ch_vdif_assembler.hpp ch_vdif_assembler_internals.hpp ch_vdif_assembler_kernels.hpp
BINFILES=run-vdif-assembler
LIBFILES=libch_vdif_assembler.so
LIBCYTHON=ch_vdif_assembler_cython.so
PYMODULES=ch_vdif_assembler.py
SCRIPTS=show-moose-acquisitions.py index_vdif_waterfalls.py
TESTBINFILES=test-timestamp-unwrapper test-kernels time-kernels peek-at-kernels

OFILES=assembler_nerve_center.o assembler_thread.o disk_reader_thread.o disk_writer_thread.o misc.o network_thread.o processing_thread.o rfi_histogrammer.o sim_thread.o timing_thread.o unit_testing_thread.o waterfall_plotter.o


all: $(BINFILES) $(LIBFILES) $(LIBCYTHON) $(TESTBINFILES)

cython: $(LIBCYTHON)

%.o: %.cpp $(INCFILES)
	$(CPP) -c -o $@ $<

%_cython.cpp: %_cython.pyx ch_vdif_assembler_pxd.pxd ch_vdif_assembler_cython.hpp $(INCFILES)
	cython --cplus $<

libch_vdif_assembler.so: $(OFILES)
	$(CPP) -o $@ -shared $^

ch_vdif_assembler_cython.so: ch_vdif_assembler_cython.cpp libch_vdif_assembler.so
	$(CPP) -shared -o $@ $< -lch_vdif_assembler -lhdf5 -lpng

run-vdif-assembler: run-vdif-assembler.o libch_vdif_assembler.so
	$(CPP) -o $@ $< -lch_vdif_assembler -lhdf5 -lpng

test-timestamp-unwrapper: test-timestamp-unwrapper.cpp ch_vdif_assembler_internals.hpp
	$(CPP) -o $@ $<

test-kernels: test-kernels.cpp ch_vdif_assembler_kernels.hpp
	$(CPP) -o $@ $<

time-kernels: time-kernels.cpp ch_vdif_assembler_kernels.hpp
	$(CPP) -o $@ $<

peek-at-kernels: peek-at-kernels.cpp ch_vdif_assembler_kernels.hpp
	$(CPP) -o $@ $<

test: $(TESTBINFILES)
	for f in $(TESTBINFILES); do ./$$f; done

install: $(INCFILES) $(BINFILES) $(LIBFILES) $(LIBCYTHON)
	mkdir -p $(INCDIR) $(LIBDIR) $(BINDIR) $(PYDIR)
	cp -f $(INCFILES) $(INCDIR)/
	cp -f $(LIBFILES) $(LIBDIR)/
	cp -f $(BINFILES) $(SCRIPTS) $(BINDIR)/
	cp -f $(LIBCYTHON) $(PYMODULES) $(PYDIR)/

#cp -f run-vdif-assembler show-moose-acquisitions.py $(BINDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	for f in $(LIBFILES); do rm -f $(LIBDIR)/$$f; done
	for f in $(BINFILES) $(SCRIPTS); do rm -f $(BINDIR)/$$f; done
	for f in $(LIBCYTHON) $(PYMODULES); do rm -f $(PYDIR)/$$f; done

clean:
	rm -f *~ *.o *_cython.cpp *.so *.pyc $(BINFILES) $(TESTBINFILES)
