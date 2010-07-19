## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	input.geopt \
	velocities.ref \
	H2O-slab.out.ref


PPFILES = \
	psppar.H \
	psppar.O

CLEANFILES = $(PPFILES) inguess.new inguess.dat *.pot *.pot_avg* H2O-slab.out posmd* \
	*prc fldiff.report velocities.xyz

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH-K/$@

H2O-slab.out: $(top_builddir)/src/bigdft $(PPFILES)
	@if ! test x"$(srcdir)" = x"." ; then \
	  for file in $(EXTRA_DIST); do \
	    ln -fs $(srcdir)/$$file ; \
	  done ; \
        fi
	cp -f velocities.ref velocities.xyz
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for file in $(EXTRA_DIST); do \
	    rm -f $$file ; \
	  done ; \
        fi

check: H2O-slab.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft H2O-slab.out \
		$(top_srcdir)/tests/DFT/H2O-slab/H2O-slab.out.ref | tee fldiff.report
