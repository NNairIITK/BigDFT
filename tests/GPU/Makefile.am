## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.dat \
	posinp \
	cpu_gpu_aff.config

PPFILES = \
	psppar.Si

CLEANFILES = $(PPFILES) GPU.out posout* *.prc \
			 fldiff.report

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/GPU/input.dat .
	cp -f $(top_srcdir)/tests/GPU/posinp .
	cp -f $(top_srcdir)/tests/GPU/cpu_gpu_aff.config .
	$(top_builddir)/src/cluster > GPU.out
	python $(top_srcdir)/tests/fldiff.py GPU.out $(top_srcdir)/tests/Refs/GPU.out.ref | tee fldiff.report
