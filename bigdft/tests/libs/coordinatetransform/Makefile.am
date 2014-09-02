## Process this file with automake to produce Makefile.in

EXTRA_DIST = test-coordinatetransform.ref.yaml test-bigdft-tool.ref.yaml \
	posinp1.xyz posinp2.xyz posinp3.xyz posinp4.xyz posinp5.xyz posinp5.fix

check_PROGRAMS = driver

AM_FCFLAGS = -I$(top_builddir)/includes -I. @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@  @MPI_INCLUDE@
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
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-coordinatetransform.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-coordinatetransform -o $@ 

test-coordinatetransform.out.yaml: driver
	cp -f $(srcdir)/posinp1.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver > test-coordinatetransform.out.yaml
	cp -f $(srcdir)/posinp2.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml
	cp -f $(srcdir)/posinp3.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml
	cp -f $(srcdir)/posinp4.xyz $(PWD)/posinp.xyz
	$(run_parallel) ./driver >> test-coordinatetransform.out.yaml

test-bigdft-tool.report.yaml: $(srcdir)/test-bigdft-tool.ref.yaml test-bigdft-tool.out.yaml
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d test-bigdft-tool.out.yaml -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=test-bigdft-tool -o $@ 

test-bigdft-tool.out.yaml: $(abs_top_builddir)/src/bigdft-tool
	cp -f $(srcdir)/posinp1.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml > test-bigdft-tool.out.yaml
	rm -f posout.*
	cp -f $(srcdir)/posinp2.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml >> test-bigdft-tool.out.yaml
	rm -f posout.*
	cp -f $(srcdir)/posinp3.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml >> test-bigdft-tool.out.yaml
	rm -f posout.*
	cp -f $(srcdir)/posinp4.xyz $(PWD)/posinp.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml >> test-bigdft-tool.out.yaml
	rm -f posout.*
	cp -f $(srcdir)/posinp5.xyz $(PWD)/posinp.xyz
	cp -f $(srcdir)/posinp5.fix $(PWD)/posinp.fix
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posinp.xyz --outfile=posout.int
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.int --outfile=posout.xyz
	$(abs_top_builddir)/src/bigdft-tool -a transform-coordinates --infile=posout.xyz --outfile=posout.yaml
	cat posout.yaml >> test-bigdft-tool.out.yaml
	rm -f posout.*
