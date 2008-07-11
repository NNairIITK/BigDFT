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
	ln $(top_srcdir)/utils/PSPfiles/HGH/$@ .

check: $(PPFILES) posinp input.dat
	if test ! -f input.dat ; then \
	  ln $(top_srcdir)/tests/GPU/input.dat ; \
	fi
	if test ! -f posinp ; then \
	  ln $(top_srcdir)/tests/GPU/posinp ; \
	fi
	if test ! -f cpu_gpu_aff.config ; then \
	  ln $(top_srcdir)/tests/GPU/cpu_gpu_aff.config ; \
	fi


	$(top_builddir)/src/cluster > GPU.out
	python $(top_srcdir)/tests/fldiff.py GPU.out $(top_srcdir)/tests/Refs/GPU.out.ref | tee fldiff.report
