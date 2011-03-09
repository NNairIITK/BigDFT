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

CLEANFILES = write.out read.out *.prc fldiff.*.report $(EXTRA_CLEAN) $(PPFILES)

#dependencies
$(PPFILES):
	@cp $(top_srcdir)/utils/PSPfiles/HGH-K/$@ .

#Check read and write by doing twice run
write.out: $(top_builddir)/src/bigdft $(PPFILES) 
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    cp -p $(srcdir)/$$i . ; \
	  done ; \
        fi
	@cp write.dft input.dft
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi

read.out: $(top_builddir)/src/bigdft $(PPFILES) 
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    cp -p $(srcdir)/$$i . ; \
	  done ; \
        fi
	@cp read.dft input.dft
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi

check : write.out read.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft write.out \
		$(srcdir)/write.out.ref | tee fldiff.write.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft read.out \
		$(srcdir)/read.out.ref | tee fldiff.read.report
