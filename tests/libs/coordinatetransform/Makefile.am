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


CLEANFILES = *.yaml posinp.* malloc.prc report

check: test-coordinatetransform.report.yaml
check: test-bigdft-tool.report.yaml

test-coordinatetransform.report.yaml: $(srcdir)/test-coordinatetransform.ref.yaml test-coordinatetransform.out.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-coordinatetransform.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-coordinatetransform -o $@ 

test-coordinatetransform.out.yaml: driver
	cp $(srcdir)/posinp1.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver > test-coordinatetransform.out.yaml
	cp $(srcdir)/posinp2.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml
	cp $(srcdir)/posinp3.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml
	cp $(srcdir)/posinp4.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml

test-bigdft-tool.report.yaml: $(srcdir)/test-bigdft-tool.ref.yaml test-bigdft-tool.out.yaml
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-bigdft-tool.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-bigdft-tool -o $@ 

test-bigdft-tool.out.yaml: $(abs_top_builddir)/src/bigdft-tool
	cp $(srcdir)/posinp1.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml > test-bigdft-tool.out.yaml
	rm -f posout.*
	cp $(srcdir)/posinp2.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml >> test-bigdft-tool.out.yaml
	rm -f posout.*
	cp $(srcdir)/posinp3.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml >> test-bigdft-tool.out.yaml
	rm -f posout.*
	cp $(srcdir)/posinp4.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml >> test-bigdft-tool.out.yaml
	rm -f posout.*
