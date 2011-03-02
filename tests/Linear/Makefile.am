## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = input.dft  Linear.out.ref  orbitalsValues  posinp.xyz  psppar.C  psppar.H

CLEANFILES = Linear.out fort* *.prc fldiff.report $(EXTRA_CLEAN) 

#dependencies
Linear.out: $(top_builddir)/src/bigdft
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    ln -fs $(srcdir)/$$i ; \
	  done ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi
	@for i in $(EXTRA_CLEAN) ; do \
	  rm -f $$i ; \
	done

check : Linear.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=1.e-8 Linear.out \
		$(top_srcdir)/tests/Linear/Linear.out.ref | tee fldiff.report
