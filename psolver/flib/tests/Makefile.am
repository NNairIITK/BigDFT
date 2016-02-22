## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	yaml_hl.cfg \
	yaml_hl.py \
	tols.yaml \
	yaml_hl_remarks.cfg \
	fldiff_yaml.py.in \
	report.py.in \
	yaml_test.ref.yaml $(OTHERS)

OTHERS = yaml_test.f90

check_PROGRAMS = yaml_test

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

CLEANFILES = yaml_test.out.yaml yaml_test.report.yaml

yaml_test.o: $(top_builddir)/flib/src/libflib-1.a

yaml_test.out.yaml: yaml_test
	./yaml_test > yaml_test.out.yaml

yaml_test.report.yaml: $(top_builddir)/tests/fldiff_yaml.py $(top_srcdir)/tests/tols-BigDFT.yaml $(abs_srcdir)/yaml_test.ref.yaml yaml_test.out.yaml
	python $(top_builddir)/tests/fldiff_yaml.py --label yaml_test -t $(top_srcdir)/tests/tols-BigDFT.yaml 	\
		-r $(abs_srcdir)/yaml_test.ref.yaml -d yaml_test.out.yaml --output $@

check: yaml_test.report.yaml
