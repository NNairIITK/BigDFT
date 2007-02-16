## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp 

PPFILES = \
	psppar.H_lda \
	psppar.Si_lda

CLEANFILES = $(PPFILES) inguess.dat wavefunction* *.out posout* malloc time.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln -s $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/SiH4/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/SiH4/posinp ; \
	fi
	$(top_builddir)/src/cluster > SiH4.out
	perl ../fldiff.pl SiH4.out $(top_srcdir)/tests/Refs/SiH4.out.ref | cat
