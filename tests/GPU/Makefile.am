## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp 

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) GPU.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES):
	ln -s $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln -s $(top_srcdir)/tests/GPU/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln -s $(top_srcdir)/tests/GPU/posinp ; \
	fi

	$(top_builddir)/src/cluster > GPU.out
	python $(top_srcdir)/tests/fldiff.py GPU.out $(top_srcdir)/tests/Refs/GPU.out.ref | tee fldiff.report
