
if USE_CUDA_GPU
check_PROGRAMS = conv_check
my_tests = CUDA.report
else
check_PROGRAMS = 
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

EXTRA_DIST = conv_check.awk

AM_FCFLAGS = -I$(top_builddir)/includes @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

conv_check_SOURCES = $(top_builddir)/src/CUDA/conv_check.f90  $(top_builddir)/src/CUDA/timing.c
conv_check_LDADD = -lbigdft-1

conv_check.out: conv_check
	echo 1 64 64 2048 2048 2 > fort.1
	./conv_check > conv_check.out

CUDA.report: conv_check.out
	awk -f $(srcdir)/conv_check.awk $< | tee CUDA.report

check: $(my_tests)
