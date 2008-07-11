## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp 

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) Davidson-SiH4.out posout* *.pot *.prc\
			 fldiff.report

#dependencies
$(PPFILES):
	ln $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln $(top_srcdir)/tests/Davidson-SiH4/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln $(top_srcdir)/tests/Davidson-SiH4/posinp ; \
	fi
	$(run_parallel) $(top_builddir)/src/cluster > Davidson-SiH4.out
	python $(top_srcdir)/tests/fldiff.py Davidson-SiH4.out $(top_srcdir)/tests/Refs/Davidson-SiH4.out.ref | tee fldiff.report
