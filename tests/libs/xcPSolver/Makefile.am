## Process this file with automake to produce Makefile.in

EXTRA_DIST = PS_Check.ref.yaml tols.yaml

check_PROGRAMS = PS_Check PS_Integral

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

PS_Check_SOURCES = PS_Check.f90
PS_Check_LDADD = -lbigdft-1 @LINALG_LIBS@

PS_Integral_SOURCES = PS_Integral.f90 fisf.f90
PS_Integral_LDADD = -lbigdft-1 @LINALG_LIBS@


CLEANFILES = PS_Check.out.yaml fldiff.report.yaml

PS_Check.out.yaml: PS_Check
	$(run_parallel) ./PS_Check 57 48 63  0 F  > PS_Check.out.yaml
	$(run_parallel) ./PS_Check 57 48 63  1 F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 57 48 63 11 F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 57 48 63 13 F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64  0 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64  1 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64 11 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64 13 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48  0 S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48  1 S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48 11 S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48 13 S >> PS_Check.out.yaml

check: PS_Check.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/PS_Check.ref.yaml -d PS_Check.out.yaml -t $(abs_srcdir)/tols.yaml \
	        --output fldiff.report.yaml

#	python $(top_srcdir)/tests/fldiff.py --mode=psolver --discrepancy=2.e-8 \
#		PS_Check.out $(srcdir)/PS_Check.out.ref | tee fldiff.report
