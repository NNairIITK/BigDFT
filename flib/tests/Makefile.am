## Process this file with automake to produce Makefile.in

OTHERS = yaml_test.f90 yaml_test2.f90

EXTRA_DIST = yaml_test.ref.yaml $(OTHERS)

check_PROGRAMS = yaml_test yaml_test2

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

AM_FCFLAGS = -I../src
AM_LDFLAGS = -L$(top_builddir)/flib/src
#LIBS += -lflib -lrt @LIB_YAML_LIBS@

if USE_MPI
mpi_source =
else
mpi_source = $(top_srcdir)/flib/src/MPIfake.f90
endif

if BUILD_LIB_PSOLVER
# Static case always
flib_libraries = $(top_builddir)/flib/src/libflib-1.a
LIBS += @LIB_YAML_LIBS@
else
if BUILD_DYNAMIC_LIBS
flib_libraries = -lflib-1
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/flib/src
else
flib_libraries = $(top_builddir)/flib/src/libflib-1.a
LIBS += @LIB_YAML_LIBS@
endif
endif

yaml_test_SOURCES = yaml_test.f90 errs.f90 yamlout.f90 dicts.f90 dynmem.f90 yaml_invoice.f90 utls.f90 $(mpi_source)
yaml_test_LDADD = $(flib_libraries)
yaml_test2_SOURCES = yaml_test2.f90 errs.f90 yamlout.f90 dicts.f90 dynmem.f90 yaml_invoice.f90 utls.f90 $(mpi_source)
yaml_test2_LDADD = $(flib_libraries)

CLEANFILES = yaml_test.out.yaml yaml_test.report.yaml

yaml_test.o: $(top_builddir)/flib/src/libflib-1.a
yaml_test2.o: $(top_builddir)/flib/src/libflib-1.a

yaml_test.out.yaml: yaml_test yaml_test2
	./yaml_test > yaml_test.out.yaml
	./yaml_test2 >> yaml_test.out.yaml

yaml_test.report.yaml: $(top_srcdir)/tests/fldiff_yaml.py $(top_srcdir)/tests/tols-BigDFT.yaml $(abs_srcdir)/yaml_test.ref.yaml yaml_test.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py --label yaml_test -t $(top_srcdir)/tests/tols-BigDFT.yaml 	\
		-r $(abs_srcdir)/yaml_test.ref.yaml -d yaml_test.out.yaml --output $@

check: yaml_test.report.yaml

