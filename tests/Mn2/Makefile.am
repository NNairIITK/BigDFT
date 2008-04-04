## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	moments   \
	input.dat \
	posinp

PPFILES = \
	psppar.Mn

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* Mn2.out posout* \
			 *.prc Mn2.memguess fldiff.report 

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
	if test ! -f moments ; then \
	  ln -s $(top_srcdir)/tests/Mn2/moments ; \
	fi
	$(top_builddir)/src/memguess 1 > Mn2.memguess
	python $(top_srcdir)/tests/fldiff.py Mn2.memguess $(top_srcdir)/tests/Refs/Mn2.memguess.ref | cat
	$(run_parallel) $(top_builddir)/src/cluster > Mn2.out
	python $(top_srcdir)/tests/fldiff.py Mn2.out $(top_srcdir)/tests/Refs/Mn2.out.ref | tee fldiff.report
