## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	input.abscalc \
	posinp.xyz \
	gproje_Fe_1s_1

PPFILES = \
	psppar.Fe

CLEANFILES = $(PPFILES) XabsCheb.out posout* *.pot *.prc coeffs_chebychev cheb_spectra_199 cheb_spectra_399 \
			 fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

XabsCheb.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$srcdir" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/input.abscalc ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/gproje_Fe_1s_1 ; \
	fi
	$(top_builddir)/src/cluster > $@
	if ! test x"$srcdir" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f input.abscalc ; \
	  rm -f gproje_Fe_1s_1 ; \
	fi

check: XabsCheb.out
	python $(top_srcdir)/tests/fldiff.py --bigdft XabsCheb.out $(top_srcdir)/tests/Refs/XabsCheb.out.ref | tee fldiff.report
