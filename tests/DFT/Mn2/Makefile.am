## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	moments   \
	input.dft \
	posinp.xyz \
	Mn2.memguess.ref \
	Mn2.out.ref

PPFILES = \
	psppar.Mn

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* Mn2.out posout* \
			 *.prc Mn2.memguess fldiff.report memguess.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

Mn2.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/moments ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f moments ; \
        fi

Mn2.memguess: $(top_builddir)/src/memguess $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/moments ; \
        fi
	$(top_builddir)/src/memguess 1 > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f moments ; \
        fi

check: Mn2.memguess Mn2.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft Mn2.memguess \
		$(top_srcdir)/tests/DFT/Mn2/Mn2.memguess.ref | tee memguess.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft Mn2.out \
		$(top_srcdir)/tests/DFT/Mn2/Mn2.out.ref | tee fldiff.report
