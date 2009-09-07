## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	input.geopt \
	posinp.ascii

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) SiH4.out geopt.mon posout* *.prc\
			fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH-K/$@

SiH4.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$srcdir" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.ascii ; \
	  ln -fs $(srcdir)/input.geopt ; \
        fi
	$(run_parallel) $(top_builddir)/src/cluster > $@
	if ! test x"$srcdir" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.ascii ; \
	  rm -f input.geopt ; \
        fi

check: SiH4.out
	python $(top_srcdir)/tests/fldiff.py --bigdft SiH4.out $(top_srcdir)/tests/Refs/SiH4.out.ref | tee fldiff.report
