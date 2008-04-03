## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	moments   \
	input.dat \
	posinp

PPFILES = \
	psppar.Mn

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* Mn.out posout* \
			 *.prc Mn.memguess fldiff.report 

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/Mn2/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/Mn2/posinp ; \
	fi
	$(top_builddir)/src/memguess 1 > Mn.memguess
	python $(top_srcdir)/tests/fldiff.py Mn.memguess $(top_srcdir)/tests/Refs/Mn.memguess.ref | cat
	$(run_parallel) $(top_builddir)/src/cluster > Mn.out
	python $(top_srcdir)/tests/fldiff.py Mn.out $(top_srcdir)/tests/Refs/Mn.out.ref | tee fldiff.report
