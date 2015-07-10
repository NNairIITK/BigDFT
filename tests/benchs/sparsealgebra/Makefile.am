## Process this file with automake to produce Makefile.in

EXTRA_DIST = test-smatmul.ref.yaml mat_114_8986.dat  mat_1202_77462.dat  mat_600_33620.dat tols.yaml

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
	cp -f $(srcdir)/mat_114_8986.dat $(PWD)/matrix.dat
	$(run_parallel) ./smatmul --filename=matrix.dat --nit=10 > test-smatmul.out.yaml
	cp -f $(srcdir)/mat_600_33620.dat $(PWD)/matrix.dat
	$(run_parallel) ./smatmul --filename=matrix.dat --nit=100 --verbosity=1 >> test-smatmul.out.yaml
	cp -f $(srcdir)/mat_1202_77462.dat  $(PWD)/matrix.dat
	$(run_parallel) ./smatmul --filename=matrix.dat --nit=20 --verbosity=1 >> test-smatmul.out.yaml


check: test-smatmul.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/test-smatmul.ref.yaml -d test-smatmul.out.yaml -t $(abs_srcdir)/tols.yaml \
	        --output fldiff.report.yaml


report: 
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
