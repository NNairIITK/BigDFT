## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.Na \
	psppar.Cl

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* NaCl.out posout* *.prc\
			 fldiff.report

#dependencies
psppar.Na:
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .

psppar.Cl:
	cp $(top_srcdir)/utils/PSPfiles/$@ .

inguess.dat:
	cp $(top_srcdir)/tests/inguess.dat .

check: $(PPFILES) inguess.dat
	cp -f $(top_srcdir)/tests/NaCl/input.dat .
	cp -f $(top_srcdir)/tests/NaCl/posinp .
	$(run_parallel) $(top_builddir)/src/cluster > NaCl.out
	python $(top_srcdir)/tests/fldiff.py NaCl.out $(top_srcdir)/tests/Refs/NaCl.out.ref | tee fldiff.report
