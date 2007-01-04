## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.C_lda 

CLEANFILES = $(PPFILES) inguess.dat wavefunction* *.out posout* malloc *.out time.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln -s $(top_srcdir)/test/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/test/C/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/test/C/posinp ; \
	fi
	$(top_builddir)/src/cluster > C.out
	diff C.out $(top_srcdir)/test/Refs/C.out.ref | cat
