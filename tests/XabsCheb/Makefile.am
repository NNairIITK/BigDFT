## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	inputRef.dft \
	inputRef.abscalc \
	posinpRef.xyz \
	gproje_Fe_1s_1 \
	XabsCheb.out.ref \
	inputGenPot.dft	\
	posinpGenPot.xyz 

PPFILES = \
	psppar.Fe

CLEANFILES = $(PPFILES) Xabs*.out posout* *.pot *.prc *.cube *.report \
	coeffs_chebychev cheb_spectra_199 cheb_spectra_399 \
	fldiff.report  cheb_spectra_* alpha_beta* \
	alphabeta* local_potential* potx.dat

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

XabsCheb.out: $(top_builddir)/src/abscalc $(PPFILES)
	@if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	@ln -fs $(srcdir)/inputRef.abscalc input.abscalc
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	$(top_builddir)/src/abscalc > $@
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc
	@if ! test x"$(srcdir)" = x"." ; then \
	  rm -f gproje_Fe_1s_1 ; \
	fi

XabsLanc.out: $(top_builddir)/src/abscalc $(PPFILES)
	@if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	@ln -fs $(srcdir)/input.abscalcLanc input.abscalc
	$(top_builddir)/src/abscalc > $@
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc
	@if ! test x"$(srcdir)" = x"." ; then \
	  rm -f gproje_Fe_1s_1 ; \
	fi


Xabsb2B.out: $(top_builddir)/src/abscalc $(PPFILES)
	@if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	@ln -fs $(srcdir)/inputGenPot.dft input.dft
	@ln -fs $(srcdir)/posinpGenPot.xyz posinp.xyz
	@ln -fs $(srcdir)/inputRef.abscalc input.abscalc
	$(top_builddir)/src/abscalc > Xabsdel.out
	mv -f local_potentialb2B.cube  b2B_xanes.cube
	cp posinp.xyz  b2B_xanes.xyz
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	@ln -fs $(srcdir)/inputb2B.abscalc input.abscalc
	$(top_builddir)/src/abscalc > $@
	rm -f  b2B_xanes.xyz
	rm -f  b2B_xanes.cube
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc
	@if ! test x"$(srcdir)" = x"." ; then \
	  rm -f gproje_Fe_1s_1 ; \
	fi
#dependencies for make check in parallel
Xabsb2B.out: XabsLanc.out
XabsLanc.out: XabsCheb.out

check: XabsCheb.out  XabsLanc.out Xabsb2B.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft XabsCheb.out \
		$(top_srcdir)/tests/XabsCheb/XabsCheb.out.ref | tee XabsCheb.out.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft XabsLanc.out \
		$(top_srcdir)/tests/XabsCheb/XabsLanc.out.ref | tee XabsLanc.out.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft Xabsdel.out \
		$(top_srcdir)/tests/XabsCheb/Xabsdel.out.ref | tee Xabsdel.out.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft Xabsb2B.out \
		$(top_srcdir)/tests/XabsCheb/Xabsb2B.out.ref | tee Xabsb2B.out.report
