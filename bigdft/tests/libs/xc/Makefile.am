## Process this file with automake to produce Makefile.in

SHORT_CHECK=
LONG_CHECK=xc.out.yaml xc.ref.yaml xc.report.yaml

include $(srcdir)/../check.mk

if HAVE_PYYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif


SUBDIRS = 

EXTRA_DIST =

CLEANFILES = xc.out xc.tmp.* xc.report

check_PROGRAMS = test_xc

AM_FCFLAGS = @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

test_xc_SOURCES = test_xc.f90 
test_xc_LDADD = -lbigdft-1

check: xc.report.yaml

xc.out.yaml: test_xc
	-$(run_parallel) ./test_xc > $@

xc.report.yaml: xc.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(pythondir)/fldiff_yaml.py -r $(srcdir)/xc.ref.yaml -d $^ -o $@ \
			-t $(top_srcdir)/tests/tols-BigDFT.yaml --label=xc
