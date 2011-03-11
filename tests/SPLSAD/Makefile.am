## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = input.dft  input.splsad  ll_input.dft  posinp2.xyz  posinp.xyz  psppar.H  psppar.N SPLSAD.out.ref

EXTRA_CLEAN = 

PPFILES =

CLEANFILES = SPLSAD.out *.prc fort* vogt* nogt* anchorpoints???.xyz path???.xyz fldiff.report $(EXTRA_CLEAN) $(PPFILES)


#Check also the restart system (frequencies.res file)
#by doing twice run
SPLSAD.out: $(top_builddir)/src/splsad  
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    ln -fs $(srcdir)/$$i ; \
	  done ; \
        fi
	$(run_parallel) $(top_builddir)/src/splsad > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi
	@for i in $(EXTRA_CLEAN) ; do \
	  rm -f $$i ; \
	done

check : SPLSAD.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=1.e-9 SPLSAD.out \
		$(top_srcdir)/tests/SPLSAD/SPLSAD.out.ref | tee fldiff.report
