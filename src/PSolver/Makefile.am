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
	perfdata.inc \
	$(SOURCES_MOD)

lib_LIBRARIES = libpoissonsolver.a

SUBDIRS =

AM_FCFLAGS = @MPI_INCLUDE@ -I. -I$(srcdir) @XC_INCLUDE@

if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

libpoissonsolver_a_SOURCES = \
	Build_Kernel_OLD.f90 \
	fft3d_OLD.f90 \
	PSolver_Kernel.f90 \
	pconv.f90 \
	pfft3d.f90 \
	scaling_function_OLD.f90 \
	gequad.f90 \
	xcenergy_OLD.f90 \
	3Dgradient_OLD.f90 \
	Poisson_Solver.f90 

CLEANFILES = $(mpi_include)

#dependencies
pconv.o: $(mpi_include)

Poisson_Solver.o: $(mpi_include)

mpif.h:
	touch mpif.h
	echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h
	echo "integer :: MPI_DOUBLE_PRECISION,MPI_MAX,MPI_MIN" >> mpif.h

SOURCES_MOD = PSolver_launch.f90 \
	Build_Kernel.f90 \
	PSolver_Base.f90 \
	xcenergy.f90 \
	3Dgradient.f90 \
	fft3d.f90 \
	scaling_function.f90

Poisson_Solver.o : $(SOURCES_MOD)