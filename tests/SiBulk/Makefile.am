## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp 

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* SiBulk.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/SiBulk/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/SiBulk/posinp ; \
	fi

	$(run_parallel) $(top_builddir)/src/cluster > SiBulk.out
	python $(top_srcdir)/tests/fldiff.py SiBulk.out $(top_srcdir)/tests/Refs/SiBulk.out.ref | tee fldiff.report
