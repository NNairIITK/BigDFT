## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H-spin.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES): 
	ln -fs $(top_srcdir)/utils/PSPfiles/HGH/$@

H-spin.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
        fi
	$(run_parallel) $(top_builddir)/src/cluster > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
        fi

check: H-spin.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft H-spin.out $(top_srcdir)/tests/Refs/H-spin.out.ref | tee fldiff.report
