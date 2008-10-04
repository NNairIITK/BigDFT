## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H.out posout* \
			 *.prc H.memguess fldiff.report $(EXTRA_DIST)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/H/input.dat .
	cp -f $(top_srcdir)/tests/H/posinp .
	$(top_builddir)/src/memguess 1 > H.memguess
	python $(top_srcdir)/tests/fldiff.py --bigdft H.memguess $(top_srcdir)/tests/Refs/H.memguess.ref | cat
	$(run_parallel) $(top_builddir)/src/cluster > H.out
	python $(top_srcdir)/tests/fldiff.py --bigdft H.out $(top_srcdir)/tests/Refs/H.out.ref | tee fldiff.report
