## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.kpt \
	write.dft \
	read.dft \
	posinp.xyz \
	write.out.ref \
	read.out.ref

EXTRA_CLEAN = wavefunction.etsf

PPFILES = psppar.O

CLEANFILES = write.out read.out *.prc fldiff.*.report geopt.mon input.*_default $(EXTRA_CLEAN) $(PPFILES)

#dependencies
$(PPFILES):
	@cp $(top_srcdir)/utils/PSPfiles/HGH-K/$@ .

in:
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    cp -p $(srcdir)/$$i . ; \
	  done ; \
        fi

clean.in:
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi

#Check read and write by doing twice run
write.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft input.dft
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@rm -f input.dft
	$(MAKE) clean.in

read.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft 
	$(MAKE) in
	@cp read.dft input.dft
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@rm -f input.dft
	$(MAKE) clean.in

check : write.out write.out.ref read.out read.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft write.out \
		$(srcdir)/write.out.ref | tee fldiff.write.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft read.out \
		$(srcdir)/read.out.ref | tee fldiff.read.report
