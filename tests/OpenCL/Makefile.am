EXTRA_DIST = \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

TESTS = run_check.sh
if USE_OCL
check_PROGRAMS = conv_check
else
check_PROGRAMS =
endif


CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 malloc.prc *linkinfo conv_check.out OpenCL.report run_check.sh

AM_FCFLAGS = -I../../src/modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = $(top_srcdir)/src/OpenCL/conv_check.f90 
conv_check_LDADD =  $(top_builddir)/src/libbigdft.a \
	$(top_builddir)/src/PSolver/libpoissonsolver.a \
	$(top_builddir)/src/OpenCL/libbasec.a \
	$(top_builddir)/src/OpenCL/libbasef.a
conv_check_LDFLAGS = -lrt


conv_check.out: conv_check
	echo 1 128 128 16384 16384 2 > fort.1
	./conv_check > conv_check.out

run_check.sh: conv_check.out
	echo "#!/bin/sh" > run_check.sh
	echo "awk -f $(srcdir)/conv_check.awk $< | tee OpenCL.report" >> run_check.sh
	chmod u+x run_check.sh
