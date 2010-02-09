## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz 

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) Davidson-SiH4.out posout* *.pot *.prc\
			 fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

Davidson-SiH4.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
        fi
	$(run_parallel) $(top_builddir)/src/cluster > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
        fi

check: Davidson-SiH4.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.e-9 --bigdft Davidson-SiH4.out $(top_srcdir)/tests/Refs/Davidson-SiH4.out.ref | tee fldiff.report
