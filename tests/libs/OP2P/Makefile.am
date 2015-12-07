EXTRA_DIST = \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

check_PROGRAMS =
my_tests = op2p_test

EXTRA_PROGRAMS = op2p_test op2p_check kpts_check mpi_check

CLEANFILES = op2p_check.out kpts_check.out\
	mpi_layer.@MODULE_EXT@ MPI_LAYER.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@ OVERLAP_POINT_TO_POINT.@MODULE_EXT@

AM_FCFLAGS = -I$(top_builddir)/includes @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

op2p_check_SOURCES = mpi_layer.f90 example.f90 
op2p_check_LDADD = -lbigdft-1

kpts_check_SOURCES = kpts_check.f90
kpts_check_LDADD = -lbigdft-1

mpi_check_SOURCES = mpi_check.f90
mpi_check_LDADD = -lbigdft-1

op2p_test_SOURCES = op2p_check.f90
op2p_test_LDADD = -lbigdft-1

op2p_check.o: $(top_builddir)/src/libbigdft-1.a

op2p_check.out: op2p_check
	$(run_parallel)	./op2p_check > $@

check: $(my_tests)
