## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) H2-freq.out \
			 *.prc fldiff.report $(EXTRA_DIST)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/H2-freq/input.dat .
	cp -f $(top_srcdir)/tests/H2-freq/posinp .
	$(run_parallel) $(top_builddir)/src/frequencies > H2-freq.out
	python $(top_srcdir)/tests/fldiff.py --bigdft H2-freq.out $(top_srcdir)/tests/Refs/H2-freq.out.ref | tee fldiff.report
