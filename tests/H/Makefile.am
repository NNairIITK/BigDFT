## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H.out posout* \
			 *.prc H.memguess memguess.report fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

H.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
    fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
    fi

H.memguess: $(top_builddir)/src/memguess $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
    fi
	$(top_builddir)/src/memguess 1 > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
    fi

check: H.memguess H.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft H.memguess $(top_srcdir)/tests/Refs/H.memguess.ref | tee memguess.report
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft H.out $(top_srcdir)/tests/Refs/H.out.ref | tee fldiff.report
