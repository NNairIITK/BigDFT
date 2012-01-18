## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = input.dft input.inter C6H6.nnkp BigDFT2Wannier.out.ref posinp.xyz psppar.C psppar.H virtuals.etsf wavefunction.etsf

CLEANFILES = BigDFT2Wannier.out *amn *mmn *_default fort* *.prc fldiff.report $(EXTRA_CLEAN) 

#dependencies
BigDFT2Wannier.out: $(top_builddir)/src/bigdft
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    ln -fs $(srcdir)/$$i ; \
	  done ; \
        fi
	$(run_parallel) $(top_builddir)/src/BigDFT2Wannier > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi
	@for i in $(EXTRA_CLEAN) ; do \
	  rm -f $$i ; \
	done

check : BigDFT2Wannier.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=1.e-8 BigDFT2Wannier.out \
		$(top_srcdir)/tests/BigDFT2Wannier/BigDFT2Wannier.out.ref | tee fldiff.report
