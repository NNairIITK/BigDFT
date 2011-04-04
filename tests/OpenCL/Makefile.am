EXTRA_DIST = \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

if USE_OCL
check_PROGRAMS = conv_check
else
check_PROGRAMS =
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 malloc.prc *linkinfo conv_check.out OpenCL.report

AM_FCFLAGS = -I../../src/modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = ../../src/OpenCL/conv_check.f90 
conv_check_LDADD =  ../../src/libbigdft.a ../../src/PSolver/libpoissonsolver.a ../../src/OpenCL/libbasec.a ../../src/OpenCL/libbasef.a
conv_check_LDFLAGS = -lrt


conv_check.out:
	echo 1 128 128 16384 16384 2 > fort.1
	./conv_check > conv_check.out

check: conv_check.out
	awk -f $(top_srcdir)/tests/OpenCL/conv_check.awk $< | tee OpenCL.report
