## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp \
	SiH4.out.ref

PPFILES = \
	psppar.H_lda \
	psppar.Si_lda

CLEANFILES = $(PPFILES) inguess.dat wavefunction* *.out posout* malloc time.prc

#dependencies
$(PPFILES):
	ln -s ../../utils/PSPfiles/$@ .

inguess.dat:
	ln -s ../inguess.dat .

check: $(PPFILES) inguess.dat
	../../src/cluster > SiH4.out
	diff SiH4.out ../Refs/SiH4.out.ref | cat
