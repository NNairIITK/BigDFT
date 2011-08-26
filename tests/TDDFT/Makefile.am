## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = input.dft input.tddft posinp.xyz TDDFT.out.ref

EXTRA_CLEAN = 

PPFILES = psppar.N

CLEANFILES = TDDFT.out *.prc fort* fldiff.report *mon $(EXTRA_CLEAN) $(PPFILES)

#dependencies
$(PPFILES):
	@cp $(top_srcdir)/utils/PSPfiles/$@ .


TDDFT.out: $(top_builddir)/src/bigdft $(PPFILES) 
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

check : TDDFT.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=2.e-9 TDDFT.out \
		$(top_srcdir)/tests/TDDFT/TDDFT.out.ref | tee fldiff.report
