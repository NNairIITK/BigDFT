## Process this file with automake to produce Makefile.in

SOURCES_MOD = PSolver_launch.f90 \
	Build_Kernel.f90 \
	PSolver_Base.f90 \
	xcenergy.f90 \
	3Dgradient.f90 \
	fft3d.f90 \
	scaling_function.f90

OTHERS = Makefile.PSolver \
	bench.csh \
	perfdata.t \
	acc_F.20-100.ref \
	acc_S-128.20-100.ref \
	acc_P.20-100.ref \
	README \
	INSTALL \
	COPYING \
	TODO 

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
	$(SOURCES_MOD) \
	$(OTHERS)

lib_LIBRARIES = libpoissonsolver.a


if CAPITALIZE
module_DATA = POISSON_SOLVER.@MODULE_EXT@
else
module_DATA = poisson_solver.@MODULE_EXT@
endif

SUBDIRS =

check_PROGRAMS = PSolver 

AM_FCFLAGS = @MPI_INCLUDE@ -I. -I$(srcdir) @XC_INCLUDE@



if USE_MPI
mpi_source =
mpi_include =
else
mpi_source = ../MPIfake.f90
mpi_include = mpif.h
endif

libpoissonsolver_a_SOURCES = \
	Poisson_Solver.f90 

PSolver_SOURCES = PSolver.f90 timing.f90 $(mpi_source)
PSolver_LDFLAGS = $(LMPI)
PSolver_LDADD = libpoissonsolver.a \
	../ABINIT-common/libABINIT.a \
	@MPI_LIBS@


CLEANFILES = mpif.h POISSON_SOLVER.@MODULE_EXT@ poisson_solver.@MODULE_EXT@

#dependencies
mpif.h:
	touch mpif.h
	echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h
	echo "integer :: MPI_DOUBLE_PRECISION,MPI_MAX,MPI_MIN" >> mpif.h

poisson_solver.@MODULE_EXT@: Poisson_Solver.o

Poisson_Solver.o : $(mpi_include) $(SOURCES_MOD)

PSolver.o: Poisson_Solver.o
