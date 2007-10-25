if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@ -I../PSolver

noinst_LIBRARIES = libprofiling.a

libprofiling_a_SOURCES = timing.f90

CLEANFILES = mpif.h

#dependencies
timing.o: $(mpi_include)

mpif.h:
	touch mpif.h
	echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h
	echo "integer :: MPI_DOUBLE_PRECISION" >> mpif.h
	echo "integer :: MPI_MIN, MPI_MAX" >> mpif.h
