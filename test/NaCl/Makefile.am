## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.Na_lda \
	psppar.Cl_lda

CLEANFILES = $(PPFILES) inguess.dat wavefunction* *.out posout* malloc time.prc

#dependencies
$(PPFILES):
	ln -s ../../utils/PSPfiles/$@ .

inguess.dat:
	ln -s ../inguess.dat .

check: $(PPFILES) inguess.dat
	../../src/cluster > NaCl.out
	diff NaCl.out ../Refs/NaCl.out.ref | cat
