## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	GPU.config

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) GPU.out posout* *.prc \
			 fldiff.report

#dependencies 
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

GPU.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/GPU.config ; \
        fi
	$(top_builddir)/src/cluster > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f GPU.config ; \
        fi

check: GPU.out
	python $(top_srcdir)/tests/fldiff.py --bigdft GPU.out $(top_srcdir)/tests/Refs/GPU.out.ref | tee fldiff.report
