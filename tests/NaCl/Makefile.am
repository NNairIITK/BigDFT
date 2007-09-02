## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.Na_lda \
	psppar.Cl_lda

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* *.out posout* *.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln -s $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/NaCl/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/NaCl/posinp ; \
	fi
	$(top_builddir)/src/cluster > NaCl.out
	perl $(top_srcdir)/tests/fldiff.pl NaCl.out $(top_srcdir)/tests/Refs/NaCl.out.ref | cat
