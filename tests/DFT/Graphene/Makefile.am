## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	input.kpt \
	posinp.ascii \
	Graphene.out.ref

PPFILES = \
	psppar.C

CLEANFILES = $(PPFILES) Graphene.out *.prc\
			fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

Graphene.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.ascii ; \
	  ln -fs $(srcdir)/input.kpt ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.ascii ; \
	  rm -f input.kpt ; \
        fi

check: Graphene.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft Graphene.out \
		   $(top_srcdir)/tests/DFT/Graphene/Graphene.out.ref | tee fldiff.report
