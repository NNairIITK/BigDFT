## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	posinp.xyz \
	Ca2.out.ref

PPFILES = \
	psppar.Ca

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* Ca2.out posout* *.prc\
			 fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

inguess.dat:
	ln -fs $(top_srcdir)/tests/DFT/inguess.dat

Ca2.out: $(top_builddir)/src/bigdft $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/posinp.xyz ; \
        fi
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f input.dft ; \
	  rm -f posinp.xyz ; \
        fi

#Use 4.0e-10 to test fldiff.py and also difference betwen 1 (gfortran) and many cores (g95) (TD)
check: inguess.dat Ca2.out
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft --discrepancy=4.e-10 Ca2.out \
		$(top_srcdir)/tests/DFT/Ca2/Ca2.out.ref | tee fldiff.report
