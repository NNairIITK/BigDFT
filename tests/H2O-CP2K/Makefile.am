## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp \
	gaubasis.dat \
	gaucoeff.dat 

PPFILES = \
	psppar.H \
	psppar.O

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H2O-CP2K.out posout* *.prc \
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/H2O-CP2K/input.dat .
	cp -f $(top_srcdir)/tests/H2O-CP2K/posinp .
	cp -f $(top_srcdir)/tests/H2O-CP2K/gaucoeff.dat .
	cp -f $(top_srcdir)/tests/H2O-CP2K/gaubasis.dat .
	$(run_parallel) $(top_builddir)/src/cluster > H2O-CP2K.out
	python $(top_srcdir)/tests/fldiff.py H2O-CP2K.out $(top_srcdir)/tests/Refs/H2O-CP2K.out.ref | tee fldiff.report
