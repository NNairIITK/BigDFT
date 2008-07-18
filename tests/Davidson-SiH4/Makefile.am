## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp 

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) Davidson-SiH4.out posout* *.pot *.prc\
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES)
	cp -f $(top_srcdir)/tests/Davidson-SiH4/input.dat .
	cp -f $(top_srcdir)/tests/Davidson-SiH4/posinp .
	$(run_parallel) $(top_builddir)/src/cluster > Davidson-SiH4.out
	python $(top_srcdir)/tests/fldiff.py Davidson-SiH4.out $(top_srcdir)/tests/Refs/Davidson-SiH4.out.ref | tee fldiff.report
