## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	input.abscalc \
	posinp.xyz \
	gproje_Fe_1s_1

PPFILES = \
	psppar.Fe

CLEANFILES = $(PPFILES) Xabs*.out posout* *.pot *.prc coeffs_chebychev cheb_spectra_199 cheb_spectra_399 \
			 fldiff.report  del.del cheb_spectra_* alpha_beta*

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

XabsCheb.out: $(top_builddir)/src/abscalc $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/input.abscalc ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	$(top_builddir)/src/abscalc > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f input.abscalc ; \
	  rm -f gproje_Fe_1s_1 ; \
	fi

XabsLanc.out: $(top_builddir)/src/abscalc $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft input.dft  ; \
	  ln -fs $(srcdir)/input.abscalcLanc input.abscalc ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	$(top_builddir)/src/abscalc > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f input.abscalc ; \
	  rm -f gproje_Fe_1s_1 ; \
	fi


Xabsb2B.out: $(top_builddir)/src/abscalc $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/inputGenPot.dft input.dft  ; \
	  ln -fs $(srcdir)/input.abscalc  input.abscalc; \
	  ln -fs $(srcdir)/posinpGenPot.xyz   posinp.xyz ; \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	$(top_builddir)/src/abscalc > del.del
	if ! test x"$(srcdir)" = x"." ; then \
	  mv -f local_potentialb2B.cube  b2B_xanes.cube ; \
	  cp posinp.xyz  b2B_xanes.xyz ; \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f input.abscalc ; \
	  rm -f gproje_Fe_1s_1 ; \
	fi
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft input.dft  ; \
	  ln -fs $(srcdir)/inputb2B.abscalc  input.abscalc; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	$(top_builddir)/src/abscalc > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f input.abscalc ; \
	  rm -f gproje_Fe_1s_1 ; \
	  rm -f  b2B_xanes.xyz ; \
	  rm -f  b2B_xanes.cube ; \
	fi


check: XabsCheb.out  XabsLanc.out Xabsb2B.out
	cat XabsCheb.out  XabsLanc.out  Xabsb2B.out >  XabsChebLanc.out 
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft XabsChebLanc.out $(top_srcdir)/tests/Refs/XabsCheb.out.ref | tee fldiff.report
