## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp 

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* SiBulk.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES)
	cp -f $(top_srcdir)/tests/SiBulk/input.dat .
	cp -f $(top_srcdir)/tests/SiBulk/posinp .
	$(run_parallel) $(top_builddir)/src/cluster > SiBulk.out
	python $(top_srcdir)/tests/fldiff.py SiBulk.out $(top_srcdir)/tests/Refs/SiBulk.out.ref | tee fldiff.report
