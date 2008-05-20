## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp 

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* SiH4.out posout* *.prc\
			 fldiff.report

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/HGH-K/$@ .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/SiH4/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/SiH4/posinp ; \
	fi
	$(run_parallel) $(top_builddir)/src/cluster > SiH4.out
	python $(top_srcdir)/tests/fldiff.py SiH4.out $(top_srcdir)/tests/Refs/SiH4.out.ref | tee fldiff.report
