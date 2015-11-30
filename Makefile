# Makefile.local must define the following variables
#   BINDIR   install dir for executables
#   LIBDIR   install dir for C++ libraries
#   INCDIR   install dir for C++ headers
#   PYDIR    install dir for python/cython modules
#   CPP      g++ compiler command line (including -pthread or -lpthread if necessary)
#
# See Makefile.local.example for an example.

include Makefile.local

all: libch_vdif_assembler.so ch_vdif_assembler_c.so run-vdif-assembler

%.o: %.cpp ch_vdif_assembler.hpp
	$(CPP) -c -o $@ $<

libch_vdif_assembler.so: ch_vdif_assembler.o waterfall_plotter.o rfi_histogrammer.o
	$(CPP) -o $@ -shared $^

ch_vdif_assembler_c.cpp: ch_vdif_assembler_c.pyx _ch_vdif_assembler_c.pxd ch_vdif_assembler.hpp
	cython --cplus $<

ch_vdif_assembler_c.so: ch_vdif_assembler_c.cpp libch_vdif_assembler.so
	$(CPP) -shared -o $@ $< -L. -lch_vdif_assembler

run-vdif-assembler: run-vdif-assembler.o libch_vdif_assembler.so
	$(CPP) -o $@ $< -L. -lch_vdif_assembler -lhdf5 -lpng

install: libch_vdif_assembler.so ch_vdif_assembler_c.so
	cp -f ch_vdif_assembler.hpp $(INCDIR)/ch_vdif_assembler.hpp
	cp -f libch_vdif_assembler.so $(LIBDIR)/libch_vdif_assembler.so
	cp -f ch_vdif_assembler_c.so $(PYDIR)/ch_vdif_assembler_c.so
	cp -f ch_vdif_assembler.py $(PYDIR)/ch_vdif_assembler.py
	cp -f run-vdif-assembler show-moose-acquisitions.py $(BINDIR)/

clean:
	rm -f *~ *.o *.so *.pyc ch_vdif_assembler_c.cpp rfi-histogram
