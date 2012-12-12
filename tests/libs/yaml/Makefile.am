## Process this file with automake to produce Makefile.in

OTHERS = yaml_test.f90

EXTRA_DIST = yaml.ref.yaml $(OTHERS)

check_PROGRAMS = yaml_test

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

AM_FCFLAGS = -I$(top_builddir)/src/modules @MPI_INCLUDE@

yaml_test_SOURCES = yaml_test.f90
yaml_test_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@


CLEANFILES = yaml_test.out.yaml

yaml_test.out.yaml: yaml_test
	$(run_parallel) ./yaml_test > yaml_test.out.yaml

check: yaml_test.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/yaml_test.ref.yaml -d yaml_test.out.yaml --output fldiff.report.yaml

