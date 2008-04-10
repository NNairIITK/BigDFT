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
			 fldiff.report

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES) posinp input.dat gaucoeff.dat gaubasis.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/H2O-CP2K/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/H2O-CP2K/posinp ; \
	fi
	if test ! -f gaucoeff.dat ; then \
	  ln -s $(top_srcdir)/tests/H2O-CP2K/gaucoeff.dat ; \
	fi
	if test ! -f gaubasis.dat ; then \
	  ln -s $(top_srcdir)/tests/H2O-CP2K/gaubasis.dat ; \
	fi

	$(run_parallel) $(top_builddir)/src/cluster > H2O-CP2K.out
	python $(top_srcdir)/tests/fldiff.py H2O-CP2K.out $(top_srcdir)/tests/Refs/H2O-CP2K.out.ref | tee fldiff.report
