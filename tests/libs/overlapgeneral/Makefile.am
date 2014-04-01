## Process this file with automake to produce Makefile.in

EXTRA_DIST = test-overlapgeneral.out.ref.yaml

check_PROGRAMS = driver

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

AM_FCFLAGS = -I$(top_builddir)/includes -I. @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@  @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@ -lrt
endif

driver_SOURCES = driver.f90
driver_LDADD = -lbigdft-1 @LINALG_LIBS@


CLEANFILES = test-overlapgeneral.out.yaml 

test-overlapgeneral.out.yaml: driver
	cp $(srcdir)/inputdata.fake $(PWD)/inputdata.fake
	$(run_parallel) ./driver

check: test-overlapgeneral.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/test-overlapgeneral.ref.yaml -d test-overlapgeneral.out.yaml -t $(abs_srcdir)/tols.yaml \
	        --output fldiff.report.yaml

report: 
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
