## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	gaubasis.dat \
	gaucoeff.dat \
	H2O-CP2K.out.ref

PPFILES = \
	psppar.H \
	psppar.O

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H2O-CP2K.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

H2O-CP2K.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/gaucoeff.dat ; \
	  ln -fs $(srcdir)/gaubasis.dat ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f gaucoeff.dat ; \
	  rm -f gaubasis.dat ; \
        fi

check: H2O-CP2K.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft H2O-CP2K.out \
		$(top_srcdir)/tests/DFT/H2O-CP2K/H2O-CP2K.out.ref | tee fldiff.report
