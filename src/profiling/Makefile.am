AM_FCFLAGS = -I. -I../modules -I$(srcdir) @MPI_INCLUDE@ -I../PSolver -I../ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

noinst_LIBRARIES = libprofiling.a

libprofiling_a_SOURCES = time.f90 \
	memoryestimator.f90

CLEANFILES = timedata.@MODULE_EXT@ TIMEDATA.@MODULE_EXT@

#dependencies
