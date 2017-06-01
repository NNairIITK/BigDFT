# To be installed:
noinst_LIBRARIES = libxcPSolver.a

# Code source.
libxcPSolver_a_SOURCES = xcenergy.f90 \
	3Dgradient.f90 \
	psolvernc.f90

AM_FCFLAGS = -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ @MPI_INCLUDE@