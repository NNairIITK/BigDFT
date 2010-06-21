## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	input.geopt \
	posinp.ascii \
	NH3.out.ref

PPFILES = \
	psppar.H \
	psppar.N

CLEANFILES = $(PPFILES) NH3.out geopt.mon posout* *.prc\
			fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH-K/$@

NH3.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.ascii ; \
	  ln -fs $(srcdir)/input.geopt ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.ascii ; \
	  rm -f input.geopt ; \
        fi

check: NH3.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft NH3.out \
		$(top_srcdir)/tests/DFT/NH3/NH3.out.ref | tee fldiff.report
