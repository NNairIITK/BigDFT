## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.H_lda 

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* *.out posout* *.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

#inguess.dat:
#	ln -s $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/H-spin/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/H-spin/posinp ; \
	fi
	$(top_builddir)/src/cluster > H-spin.out
	perl $(top_srcdir)/tests/fldiff.pl H-spin.out $(top_srcdir)/tests/Refs/H-spin.out.ref | cat
