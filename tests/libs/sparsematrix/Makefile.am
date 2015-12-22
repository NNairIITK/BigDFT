## Process this file with automake to produce Makefile.in

SHORT_CHECK = test-sparsematrix.report.yaml
LONG_CHECK = $(SHORT_CHECK)

include $(srcdir)/../check.mk

check_PROGRAMS = driver driver_foe driver_ccs ice_ccs

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

driver_foe_SOURCES = driver_foe.f90
driver_foe_LDADD = -lbigdft-1 @LINALG_LIBS@

driver_ccs_SOURCES = driver_ccs.f90
driver_ccs_LDADD = -lbigdft-1 @LINALG_LIBS@

ice_ccs_SOURCES = utilities.f90 ice_ccs.f90
ice_ccs_LDADD = -lbigdft-1 @LINALG_LIBS@


CLEANFILES = test-sparsematrix.out.yaml test-sparsematrix.report.yaml test-sparsematrix.report

test-sparsematrix.report.yaml: $(srcdir)/test-sparsematrix.ref.yaml test-sparsematrix.out.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-sparsematrix.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-sparsematrix -o $@ 
##	python $(top_srcdir)/tests/fldiff_yaml.py $^ | tee $@

test-sparsematrix.out.yaml: driver
