## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H-spin.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .


check: $(PPFILES) posinp input.dat
	cp -f $(top_srcdir)/tests/H-spin/input.dat .
	cp -f $(top_srcdir)/tests/H-spin/posinp .
	$(run_parallel) $(top_builddir)/src/cluster > H-spin.out
	python $(top_srcdir)/tests/fldiff.py H-spin.out $(top_srcdir)/tests/Refs/H-spin.out.ref | tee fldiff.report
