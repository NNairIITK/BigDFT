## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz 

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* SiBulk.out posout* *.prc \
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES)
	cp -f $(top_srcdir)/tests/SiBulk/input.dft .
	cp -f $(top_srcdir)/tests/SiBulk/posinp.xyz .
	$(run_parallel) $(top_builddir)/src/cluster > SiBulk.out
	python $(top_srcdir)/tests/fldiff.py --bigdft SiBulk.out $(top_srcdir)/tests/Refs/SiBulk.out.ref | tee fldiff.report
