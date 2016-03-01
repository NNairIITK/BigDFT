if USE_OCL
check_PROGRAMS = conv_check
my_tests = OpenCL.report.yaml
else
check_PROGRAMS =
my_tests =
endif

EXTRA_DIST = conv_check.awk tols.yaml OpenCL.ref.yaml

CLEANFILES = fort.1 conv_check.out.yaml OpenCL.report.yaml report

AM_FCFLAGS = -I$(top_builddir)/includes @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
endif

conv_check_SOURCES = $(top_srcdir)/src/OpenCL/conv_check.f90 
conv_check_LDADD = @LIB_BIGDFT_LIBS@ @LINALG_LIBS@

conv_check.out.yaml: conv_check
	echo 76 92 130 2 > fort.1
	-./conv_check > $@

OpenCL.report.yaml: conv_check.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $(srcdir)/OpenCL.ref.yaml -d conv_check.out.yaml -t $(srcdir)/tols.yaml -o $@ 
#	awk -f $(srcdir)/conv_check.awk $< | tee $@

check: $(my_tests)
