## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	lazy_100.inc \
	lazy_14.inc \
	lazy_16.inc \
	lazy_20.inc \
	lazy_24.inc \
	lazy_30.inc \
	lazy_40.inc \
	lazy_50.inc \
	lazy_60.inc \
	lazy_8.inc \
	perfdata.inc

noinst_LIBRARIES = libPoisson.a

SUBDIRS =

AM_FCFLAGS = @MPI_INCLUDE@ -I../ABINIT-common

if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

libPoisson_a_SOURCES = \
	Build_Kernel.f90 \
	fft3d.f90 \
	PSolver_Kernel.f90 \
	pconv.f90 \
	pfft3d.f90 \
	scaling_function.f90 \
	gequad.f90 \
	xcenergy.f90 \
	3Dgradient.f90

CLEANFILES = $(mpi_include)

#dependencies
pconv.o: $(mpi_include)

mpif.h:
mpif.h:
	touch mpif.h
	echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h
	echo "integer :: MPI_DOUBLE_PRECISION,MPI_MAX,MPI_MIN" >> mpif.h

