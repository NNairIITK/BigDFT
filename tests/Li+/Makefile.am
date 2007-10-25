## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.Li_lda 

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* *.out posout* malloc *.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

#inguess.dat:
#	ln -s $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/Li+/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/Li+/posinp ; \
	fi
	$(top_builddir)/src/cluster > Li+.out
	python $(top_srcdir)/tests/fldiff.py Li+.out $(top_srcdir)/tests/Refs/Li+.out.ref | cat
