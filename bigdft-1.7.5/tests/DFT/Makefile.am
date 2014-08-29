## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  cubic \
	  postSCF 

if USE_OCL
else
  SUBDIRS += linear
endif

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif


EXTRA_DIST = \
	README \
	inguess.dat

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
