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
			 fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

NaCl.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/input.geopt ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f input.geopt ; \
        fi

check: NaCl.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=2.e-8 NaCl.out $(top_srcdir)/tests/Refs/NaCl.out.ref | tee fldiff.report
