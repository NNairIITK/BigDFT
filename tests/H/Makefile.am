## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.H

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* H.out posout* \
			 *.prc H.memguess fldiff.report 

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/H/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/H/posinp ; \
	fi
	$(top_builddir)/src/memguess 1 > H.memguess
	python $(top_srcdir)/tests/fldiff.py H.memguess $(top_srcdir)/tests/Refs/H.memguess.ref | cat
	$(run_parallel) $(top_builddir)/src/cluster > H.out
	python $(top_srcdir)/tests/fldiff.py H.out $(top_srcdir)/tests/Refs/H.out.ref | tee fldiff.report
