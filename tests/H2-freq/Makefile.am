## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.freq \
	input.dft \
	posinp.xyz \
	H2-freq.out.ref \
	H2-freq.out.2.ref

EXTRA_CLEAN = hessian.dat frequencies.xyz frequencies.res

PPFILES = psppar.H

CLEANFILES = H2-freq.out H2-freq.out *.prc fldiff.report $(EXTRA_CLEAN) $(PPFILES)

#dependencies
$(PPFILES):
	@cp $(top_srcdir)/utils/PSPfiles/$@ .

#Check also the restart system (frequencies.res file)
#by doing twice run
H2-freq.out: $(top_builddir)/src/frequencies $(PPFILES) 
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    ln -fs $(srcdir)/$$i ; \
	  done ; \
        fi
	$(run_parallel) $(top_builddir)/src/frequencies > $@
	$(run_parallel) $(top_builddir)/src/frequencies > $@.2
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi
	@for i in $(EXTRA_CLEAN) ; do \
	  rm -f $$i ; \
	done

check : H2-freq.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=2.e-10 H2-freq.out \
		$(top_srcdir)/tests/H2-freq/H2-freq.out.ref | tee fldiff.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=2.e-10 H2-freq.out.2 \
		$(top_srcdir)/tests/H2-freq/H2-freq.out.2.ref | tee fldiff.2.report
