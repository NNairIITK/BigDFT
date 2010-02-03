## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	moments   \
	input.dft \
	posinp.xyz

PPFILES = \
	psppar.Mn

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* Mn2.out posout* \
			 *.prc Mn2.memguess fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

Mn2.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/moments ; \
        fi
	$(run_parallel) $(top_builddir)/src/cluster > $@
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
	python $(top_srcdir)/tests/fldiff.py --bigdft Mn2.memguess $(top_srcdir)/tests/Refs/Mn2.memguess.ref | cat
	python $(top_srcdir)/tests/fldiff.py --bigdft Mn2.out $(top_srcdir)/tests/Refs/Mn2.out.ref | tee fldiff.report
