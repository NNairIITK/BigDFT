## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	H2-spin.out.ref

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H2-spin.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES): 
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

H2-spin.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
        fi

check: H2-spin.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft H2-spin.out \
		$(top_srcdir)/tests/DFT/H2-spin/H2-spin.out.ref | tee fldiff.report
