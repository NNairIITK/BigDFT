## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	GPU.config

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) GPU.out posout* *.prc \
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/GPU/input.dft .
	cp -f $(top_srcdir)/tests/GPU/posinp.xyz .
	cp -f $(top_srcdir)/tests/GPU/GPU.config .
	$(top_builddir)/src/cluster > GPU.out
	python $(top_srcdir)/tests/fldiff.py --bigdft GPU.out $(top_srcdir)/tests/Refs/GPU.out.ref | tee fldiff.report
