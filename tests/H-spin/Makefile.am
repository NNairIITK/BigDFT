## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H-spin.out posout* *.prc \
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .


check: $(PPFILES) posinp.xyz input.dft
	cp -f $(top_srcdir)/tests/H-spin/input.dft .
	cp -f $(top_srcdir)/tests/H-spin/posinp.xyz .
	$(run_parallel) $(top_builddir)/src/cluster > H-spin.out
	python $(top_srcdir)/tests/fldiff.py --bigdft H-spin.out $(top_srcdir)/tests/Refs/H-spin.out.ref | tee fldiff.report
