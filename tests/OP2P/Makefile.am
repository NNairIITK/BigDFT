EXTRA_DIST = \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

check_PROGRAMS =
my_tests =

EXTRA_PROGRAMS = op2p_check

if USE_DGEMMSY
libs_dgemmsy = $(top_builddir)/src/dgemmsy/libdgemmsy.a
else
libs_dgemmsy =
endif

CLEANFILES = op2p_check.out \
	mpi_layer.@MODULE_EXT@ MPI_LAYER.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@ OVERLAP_POINT_TO_POINT.@MODULE_EXT@

AM_FCFLAGS = -I../../src/modules -I../../src/ @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 

op2p_check_SOURCES = op2p_module.f90 mpi_layer.f90 example.f90 
op2p_check_LDADD =  $(top_builddir)/src/libbigdft.a \
	$(top_builddir)/src/PSolver/libpoissonsolver.a \
	$(libs_dgemmsy)

op2p_check.out: op2p_check
	$(run_parallel)	./op2p_check > $@

check: $(my_tests)
