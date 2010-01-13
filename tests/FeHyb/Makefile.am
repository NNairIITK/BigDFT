## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	occup.dat

PPFILES = \
	psppar.Fe

CLEANFILES = $(PPFILES) FeHyb.out posout* *.pot *.prc fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

$(EXTRA_DIST):
	ln -fs $(srcdir)/$@


FeHyb.out: $(top_builddir)/src/cluster $(PPFILES)
	@if ! test x"$(srcdir)" = x"." ; then \
	  for file in $(EXTRA_DIST); do \
	    ln -fs $(srcdir)/$$file ; \
	  done ; \
	fi
	$(run_parallel) $(top_builddir)/src/cluster > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for file in $(EXTRA_DIST); do \
	    rm -f $$file ; \
	  done ; \
	fi

check: FeHyb.out
	python $(top_srcdir)/tests/fldiff.py --bigdft FeHyb.out $(top_srcdir)/tests/Refs/FeHyb.out.ref | tee fldiff.report
