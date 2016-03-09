## Process this file with automake to produce Makefile.in

SHORT_CHECK = test-ccs.report.yaml test-foeccs.report.yaml
LONG_CHECK = $(SHORT_CHECK)

include $(srcdir)/../check.mk

EXTRA_DIST = test-ccs.ref.yaml \
	matrix1.dat matrix2.dat \
	test-foeccs.ref.yaml \
    overlap_ccs.dat hamiltonian_ccs.dat density_kernel_ccs.dat

check_PROGRAMS = driver_ccs driver_foe_ccs driver_single

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

driver_ccs_SOURCES = utilities.f90 driver_ccs.f90
driver_ccs_LDADD = -lbigdft-1 @LINALG_LIBS@

driver_foe_ccs_SOURCES = utilities.f90 driver_foe_ccs.f90
driver_foe_ccs_LDADD = -lbigdft-1 @LINALG_LIBS@

driver_single_SOURCES = utilities.f90 driver_single.f90
driver_single_LDADD = -lbigdft-1 @LINALG_LIBS@

CLEANFILES = test-ccs.out.yaml test-ccs.report.yaml test-ccs.report \
             test-foeccs.out.yaml test-foeccs.report.yaml test-foeccs.report

test-ccs.report.yaml: $(srcdir)/test-ccs.ref.yaml test-ccs.out.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-ccs.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-ccs -o $@ 

test-foeccs.report.yaml: $(srcdir)/test-foeccs.ref.yaml test-foeccs.out.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-foeccs.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-foeccs -o $@ 

test-ccs.out.yaml: driver_ccs
	cp -f $(srcdir)/matrix1-base.dat $(PWD)/matrix1.dat
	cp -f $(srcdir)/matrix2-base.dat $(PWD)/matrix2.dat
	$(run_parallel) ./driver_ccs > test-ccs.out.yaml

test-foeccs.out.yaml: driver_foe_ccs
	cp -f $(srcdir)/overlap_ccs-base.dat $(PWD)/overlap_ccs.dat
	cp -f $(srcdir)/hamiltonian_ccs-base.dat $(PWD)/hamiltonian_ccs.dat
	cp -f $(srcdir)/density_kernel_ccs-base.dat $(PWD)/density_kernel_ccs.dat
	$(run_parallel) ./driver_foe_ccs > test-foeccs.out.yaml
