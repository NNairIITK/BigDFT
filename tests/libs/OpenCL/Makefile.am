if USE_OCL
check_PROGRAMS = conv_check
my_tests = OpenCL.report.yaml
else
check_PROGRAMS =
my_tests =
endif

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif


EXTRA_DIST = conv_check.awk tols.yaml OpenCL.ref.yaml

CLEANFILES = fort.1 conv_check.out.yaml OpenCL.report.yaml report

AM_FCFLAGS = -I$(top_builddir)/src/modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = $(top_srcdir)/src/OpenCL/conv_check.f90 
conv_check_LDADD = @LIBS_DEPENDENCIES@

conv_check.out.yaml: conv_check
	echo 76 92 130 2 > fort.1
	./conv_check > $@

OpenCL.report.yaml: conv_check.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $(srcdir)/OpenCL.ref.yaml -d conv_check.out.yaml -t $(srcdir)/tols.yaml -o $@ 
#	awk -f $(srcdir)/conv_check.awk $< | tee $@

check: $(my_tests)
