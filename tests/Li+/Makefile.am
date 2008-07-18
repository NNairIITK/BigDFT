## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp

PPFILES = \
	psppar.Li_lda 

CLEANFILES = $(PPFILES) inguess.new inguess.dat wavefunction* Li+.out posout* malloc *.prc\
			 fldiff.report $(EXTRA_DIST) $(PPFILES)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .


check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/Li+/input.dat .
	cp -f $(top_srcdir)/tests/Li+/posinp .
	$(run_parallel) $(top_builddir)/src/cluster > Li+.out
	python $(top_srcdir)/tests/fldiff.py Li+.out $(top_srcdir)/tests/Refs/Li+.out.ref | tee fldiff.report
