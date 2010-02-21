## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.freq \
	input.dft \
	posinp.xyz

EXTRA_CLEAN = hessian.dat frequencies.xyz frequencies.res

PPFILES = psppar.H

CLEANFILES = H2-freq.out *.prc fldiff.report $(EXTRA_CLEAN) $(PPFILES)

#dependencies
$(PPFILES):
	@cp $(top_srcdir)/utils/PSPfiles/$@ .

H2-freq.out: $(top_builddir)/src/frequencies $(PPFILES) 
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    ln -fs $(srcdir)/$$i ; \
	  done ; \
        fi
	$(run_parallel) $(top_builddir)/src/frequencies > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi
	@for i in $(EXTRA_CLEAN) ; do \
	  rm -f $$i ; \
	done

check : H2-freq.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft H2-freq.out $(top_srcdir)/tests/Refs/H2-freq.out.ref | tee fldiff.report
