## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	input.geopt \
	posinp.xyz \
	input.geopt

PPFILES = \
	psppar.Na \
	psppar.Cl

CLEANFILES = $(PPFILES) NaCl.out posout* geopt.mon *.prc\
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
psppar.Na psppar.Cl:
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES)
	cp -f $(top_srcdir)/tests/NaCl/input.dft .
	cp -f $(top_srcdir)/tests/NaCl/input.geopt .
	cp -f $(top_srcdir)/tests/NaCl/posinp.xyz .
	$(run_parallel) $(top_builddir)/src/cluster > NaCl.out
	python $(top_srcdir)/tests/fldiff.py --bigdft NaCl.out $(top_srcdir)/tests/Refs/NaCl.out.ref | tee fldiff.report
