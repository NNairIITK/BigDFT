## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) H2-freq.out \
			 *.prc fldiff.report $(EXTRA_DIST)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .

H2-freq.out: $(top_builddir)/src/frequencies $(PPFILES) 
	if ! test x"$srcdir" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
        fi
	$(run_parallel) $(top_builddir)/src/frequencies > $@
	if ! test x"$srcdir" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
        fi

check : H2-freq.out
	python $(top_srcdir)/tests/fldiff.py --bigdft H2-freq.out $(top_srcdir)/tests/Refs/H2-freq.out.ref | tee fldiff.report
