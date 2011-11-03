EXTRA_DIST = \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

check_PROGRAMS =
my_tests =

EXTRA_PROGRAMS = op2p_check kpts_check

CLEANFILES = op2p_check.out kpts_check.out\
	mpi_layer.@MODULE_EXT@ MPI_LAYER.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@ OVERLAP_POINT_TO_POINT.@MODULE_EXT@

AM_FCFLAGS = -I../../src/modules -I../../src/ @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 

op2p_check_SOURCES = mpi_layer.f90 example.f90 
op2p_check_LDADD =  $(top_builddir)/src/libbigdft.a \
	$(top_builddir)/src/PSolver/libpoissonsolver.a

kpts_check_SOURCES = kpts_check.f90
kpts_check_LDADD =  $(top_builddir)/src/libbigdft.a \
	$(top_builddir)/src/PSolver/libpoissonsolver.a



op2p_check.out: op2p_check
	$(run_parallel)	./op2p_check > $@

check: $(my_tests)
