## Process this file with automake to produce Makefile.in

OTHERS = yaml_test.f90

EXTRA_DIST = yaml_test.ref.yaml $(OTHERS)

check_PROGRAMS = yaml_test

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

AM_FCFLAGS = -I$(top_builddir)/src/modules @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

yaml_test_SOURCES = yaml_test.f90
yaml_test_LDADD = -lbigdft-1


CLEANFILES = yaml_test.out.yaml

yaml_test.out.yaml: yaml_test
	./yaml_test > yaml_test.out.yaml

check: yaml_test.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/yaml_test.ref.yaml -d yaml_test.out.yaml --output fldiff.report.yaml

