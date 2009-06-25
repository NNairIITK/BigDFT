## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	moments   \
	input.dft \
	posinp.xyz

PPFILES = \
	psppar.Mn

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* Mn2.out posout* \
			 *.prc Mn2.memguess fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/Mn2/input.dft .
	cp -f $(top_srcdir)/tests/Mn2/posinp.xyz .
	cp -f $(top_srcdir)/tests/Mn2/moments .
	$(top_builddir)/src/memguess 1 > Mn2.memguess
	python $(top_srcdir)/tests/fldiff.py --bigdft Mn2.memguess $(top_srcdir)/tests/Refs/Mn2.memguess.ref | cat
	$(run_parallel) $(top_builddir)/src/cluster > Mn2.out
	python $(top_srcdir)/tests/fldiff.py --bigdft Mn2.out $(top_srcdir)/tests/Refs/Mn2.out.ref | tee fldiff.report
