## Process this file with automake to produce Makefile.in

EXTRA_DIST = test-coordinatetransform.ref.yaml

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


CLEANFILES = test-coordinatetransform.out.yaml test-coordinatetransform.report.yaml

check: test-coordinatetransform.report.yaml

test-coordinatetransform.report.yaml: $(srcdir)/test-coordinatetransform.ref.yaml test-coordinatetransform.out.yaml
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-coordinatetransform.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-coordinatetransform -o $@ 
##	python $(top_srcdir)/tests/fldiff_yaml.py $^ | tee $@

test-coordinatetransform.out.yaml: driver
	cp $(srcdir)/posinp1.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver > test-coordinatetransform.out.yaml
	cp $(srcdir)/posinp2.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml
	cp $(srcdir)/posinp3.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml
	cp $(srcdir)/posinp4.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml
