## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.Na_lda \
	psppar.Cl_lda

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* NaCl.out posout* *.prc\
			 fldiff.report

#dependencies
$(PPFILES):
	ln $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat
	if test ! -f input.dat ; then \
	  ln $(top_srcdir)/tests/NaCl/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln $(top_srcdir)/tests/NaCl/posinp ; \
	fi
	$(run_parallel) $(top_builddir)/src/cluster > NaCl.out
	python $(top_srcdir)/tests/fldiff.py NaCl.out $(top_srcdir)/tests/Refs/NaCl.out.ref | tee fldiff.report
