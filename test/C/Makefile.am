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
	ln -s ../../utils/PSPfiles/$@ .

inguess.dat:
	ln -s ../inguess.dat .

check: $(PPFILES) inguess.dat
	../../src/cluster > C.out
	diff C.out ../Refs/C.out.ref | cat
