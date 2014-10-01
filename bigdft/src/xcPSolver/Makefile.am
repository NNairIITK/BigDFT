# To be installed:
noinst_LIBRARIES = libxcPSolver.a
module_DATA = $(modules)

# Code source.
if BUILD_LIB_BIGDFT
PS_SOURCES = xcenergy.f90 \
	3Dgradient.f90 \
	psolvernc.f90
else
PS_SOURCES =
endif
libxcPSolver_a_SOURCES = $(PS_SOURCES)

AM_FCFLAGS = @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ @MPI_INCLUDE@ -I$(top_builddir)/includes