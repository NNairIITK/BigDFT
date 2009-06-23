## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dft \
	occup.dat \
	posinp.xyz

PPFILES = \
	psppar.C

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* C.out posout* *.out *.prc\
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	cp $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat
	cp -f $(top_srcdir)/tests/C/input.dft .
	cp -f $(top_srcdir)/tests/C/posinp.xyz .
	cp -f $(top_srcdir)/tests/C/occup.dat .
	$(run_parallel) $(top_builddir)/src/cluster > C.out
	python $(top_srcdir)/tests/fldiff.py --bigdft C.out $(top_srcdir)/tests/Refs/C.out.ref | tee fldiff.report
