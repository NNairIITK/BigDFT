## Process this file with automake to produce Makefile.in

SHORT_CHECK = test-overlapgeneral.report.yaml
LONG_CHECK = $(SHORT_CHECK)

include $(srcdir)/../check.mk

EXTRA_DIST = test-overlapgeneral.ref.yaml \
	inputdata1.fake inputdata2.fake inputdata3.fake

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


CLEANFILES = test-overlapgeneral.out.yaml test-overlapgeneral.report.yaml test-overlapgeneral.report

test-overlapgeneral.report.yaml: $(srcdir)/test-overlapgeneral.ref.yaml test-overlapgeneral.out.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-overlapgeneral.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-overlapgeneral -o $@ 
##	python $(top_srcdir)/tests/fldiff_yaml.py $^ | tee $@

test-overlapgeneral.out.yaml: driver
	cp -f $(srcdir)/inputdata1.fake $(PWD)/inputdata.fake
	$(run_parallel) ./driver > test-overlapgeneral.out.yaml
	cp -f $(srcdir)/inputdata2.fake $(PWD)/inputdata.fake
	$(run_parallel) ./driver >> test-overlapgeneral.out.yaml
	cp -f $(srcdir)/inputdata3.fake $(PWD)/inputdata.fake
	$(run_parallel) ./driver >> test-overlapgeneral.out.yaml
