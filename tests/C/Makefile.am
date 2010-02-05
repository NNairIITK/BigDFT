## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	occup.dat \
	posinp.xyz

PPFILES = \
	psppar.C

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* C.out posout* *.out *.prc\
			 fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

inguess.dat:
	ln -fs $(top_srcdir)/tests/inguess.dat

C.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
	  ln -fs $(srcdir)/occup.dat ; \
        fi
	$(run_parallel) $(top_builddir)/src/cluster > C.out
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
	  rm -f occup.dat ; \
        fi

check: inguess.dat C.out
	python $(top_srcdir)/tests/fldiff.py --bigdft C.out $(top_srcdir)/tests/Refs/C.out.ref | tee fldiff.report
