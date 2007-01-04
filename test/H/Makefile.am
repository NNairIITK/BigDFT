## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.H_lda 

CLEANFILES = $(PPFILES) inguess.dat wavefunction* *.out posout* malloc time.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln -s $(top_srcdir)/test/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/test/H/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/test/H/posinp ; \
	fi
	$(top_builddir)/src/cluster > H.out
	diff H.out $(top_srcdir)/test/Refs/H.out.ref | cat
