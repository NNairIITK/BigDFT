## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	occup.dat \
	posinp

PPFILES = \
	psppar.C_lda 

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* *.out posout* *.out *.prc

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	ln -s $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/C/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/C/posinp ; \
	fi
	$(top_builddir)/src/cluster > C.out
	perl $(top_srcdir)/tests/fldiff.pl C.out $(top_srcdir)/tests/Refs/C.out.ref | cat
