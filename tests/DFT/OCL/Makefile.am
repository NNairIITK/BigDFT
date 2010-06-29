## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	OCL.out.ref

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) OCL.out posout* *.prc \
			 fldiff.report

#dependencies 
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

GPU.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
        fi
	$(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
        fi

check: GPU.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft OCL.out \
		$(top_srcdir)/tests/DFT/OCL/OCL.out.ref | tee fldiff.report
