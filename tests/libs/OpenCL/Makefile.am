if USE_OCL
check_PROGRAMS = conv_check
my_tests = OpenCL.report
else
check_PROGRAMS =
my_tests =
endif

EXTRA_DIST = conv_check.awk

CLEANFILES = fort.1 conv_check.out OpenCL.report

AM_FCFLAGS = -I$(top_builddir)/src/modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = $(top_srcdir)/src/OpenCL/conv_check.f90 
conv_check_LDADD = @LIBS_DEPENDENCIES@

conv_check.out: conv_check
	echo 124 132 130 2 > fort.1
	./conv_check > $@

OpenCL.report: conv_check.out
	awk -f $(srcdir)/conv_check.awk $< | tee $@

check: $(my_tests)
