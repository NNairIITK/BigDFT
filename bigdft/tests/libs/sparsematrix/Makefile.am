## Process this file with automake to produce Makefile.in

SHORT_CHECK = test-ccs.report.yaml
LONG_CHECK = $(SHORT_CHECK)

include $(srcdir)/../check.mk

EXTRA_DIST = test-ccs.ref.yaml \
	matrix1.dat matrix2.dat

check_PROGRAMS = driver driver_ccs

if HAVE_PYYAML
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

driver_ccs_SOURCES = utilities.f90 driver_ccs.f90
driver_ccs_LDADD = -lbigdft-1 @LINALG_LIBS@

CLEANFILES = test-ccs.out.yaml test-ccs.report.yaml test-ccs.report

test-ccs.report.yaml: $(srcdir)/test-ccs.ref.yaml test-ccs.out.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-ccs.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-ccs -o $@ 
##	python $(top_srcdir)/tests/fldiff_yaml.py $^ | tee $@

test-ccs.out.yaml: driver_ccs
	cp -f $(srcdir)/matrix1.dat $(PWD)
	cp -f $(srcdir)/matrix2.dat $(PWD)
	$(run_parallel) ./driver_ccs > test-ccs.out.yaml
