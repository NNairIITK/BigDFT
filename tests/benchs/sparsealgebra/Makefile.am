## Process this file with automake to produce Makefile.in

EXTRA_DIST = 

check_PROGRAMS = smatmul

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

smatmul_SOURCES = smatmul.f90
smatmul_LDADD = -lbigdft-1 @LINALG_LIBS@

CLEANFILES = *.dat *.yaml report

test-smatmul.out.yaml: smatmul
	cp -f $(srcdir)/mat_6_36.dat $(PWD)
	$(run_parallel) ./smatmul --filename=mat_6_36.dat > test-smatmul.out.yaml


check: $(srcdir)/test-smatmul.ref.yaml test-smatmul.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/test-smatmul.ref.yaml -d test-smatmul.out.yaml -t $(abs_srcdir)/tols.yaml \
	        --output fldiff.report.yaml


report: 
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
