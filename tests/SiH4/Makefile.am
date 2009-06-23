## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	input.geopt \
	posinp.xyz 

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) SiH4.out geopt.mon posout* *.prc\
			fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH-K/$@ .

check: $(PPFILES)
	cp -f $(top_srcdir)/tests/SiH4/input.dft .
	cp -f $(top_srcdir)/tests/SiH4/input.geopt .
	cp -f $(top_srcdir)/tests/SiH4/posinp.xyz .
	$(run_parallel) $(top_builddir)/src/cluster > SiH4.out
	python $(top_srcdir)/tests/fldiff.py --bigdft SiH4.out $(top_srcdir)/tests/Refs/SiH4.out.ref | tee fldiff.report
