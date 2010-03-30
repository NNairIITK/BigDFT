## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	SiBulk.out.ref

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* SiBulk.out \
	posout* *.prc fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

SiBulk.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
        fi
	$(run_parallel) $(top_builddir)/src/cluster > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
        fi

check: SiBulk.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=5e-10 SiBulk.out \
		$(top_srcdir)/tests/DFT/SiBulk/SiBulk.out.ref | tee fldiff.report
