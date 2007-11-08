## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	occup.dat \
	posinp

PPFILES = \
	psppar.C_lda 

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* C.out posout* *.out *.prc\
			 fldiff.report

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln -s $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat occup.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/C/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/C/posinp ; \
	fi
	if test ! -f occup.dat ; then \
	  ln -s $(top_srcdir)/tests/C/occup.dat ; \
	fi
	$(run_parallel) $(top_builddir)/src/cluster > C.out
	python $(top_srcdir)/tests/fldiff.py C.out $(top_srcdir)/tests/Refs/C.out.ref | tee fldiff.report
