## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz 


PPFILES = \
	psppar.H \
	psppar.O

CLEANFILES = $(PPFILES) inguess.new inguess.dat *.pot *.pot_avg* *xyz H2O-slab.out posout* \
	*prc fldiff.report $(EXTRA_DIST) $(PPFILES)	

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH-K/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/H2O-slab/input.dft .
	cp -f $(top_srcdir)/tests/H2O-slab/posinp.xyz .
	$(run_parallel) $(top_builddir)/src/cluster > H2O-slab.out
	python $(top_srcdir)/tests/fldiff.py --bigdft H2O-slab.out $(top_srcdir)/tests/Refs/H2O-slab.out.ref | tee fldiff.report
