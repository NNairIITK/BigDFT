## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.Ca_lda

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* *.out posout* *.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln -s $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/Ca2/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/Ca2/posinp ; \
	fi
	$(top_builddir)/src/cluster > Ca2.out
	python $(top_srcdir)/tests/fldiff.py Ca2.out $(top_srcdir)/tests/Refs/Ca2.out.ref | cat
