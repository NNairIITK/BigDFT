## Process this file with automake to produce Makefile.in

if BUILD_LIB_XC
sub_dirs_xc = ABINIT-common
else
sub_dirs_xc = 
endif


SOURCES_MOD = PSolver_Main.f90 \
	Build_Kernel.f90 \
	PSolver_Base.f90 \
	xcenergy.f90 \
	3Dgradient.f90 \
	fft3d.f90 \
	scaling_function.f90

OTHERS = Makefile.PSolver \
	PS_Check.f90 \
	PS_Program.f90 \
	PS_Exercise.f90 \
	bench.csh \
	perfdata.t \
	acc_F.20-100.ref \
	acc_S-128.20-100.ref \
	acc_P.20-100.ref \
	README \
	INSTALL \
	COPYING \
	SUBDIRS \
	TODO \
	$(EXERCISE)

EXERCISE = Exercise/plot.gnuplot \
	Exercise/PS_Exercise.tex \
	Exercise/accF.pdf \
	Exercise/accS.pdf \
	Exercise/fandg.pdf \
	Exercise/PvsF.pdf \
	Exercise/rhoSpotS.pdf \
	Exercise/accP.pdf \
	Exercise/emepot.pdf \
	Exercise/rhoFpotF.pdf

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

SUBDIRS = $(sub_dirs_xc) 

check_PROGRAMS = PS_Check PS_Program PS_Exercise

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @XC_INCLUDE@ 


if USE_MPI
mpi_source =
mpi_include =
else
mpi_source = ../MPIfake.f90
mpi_include = mpif.h
endif

libpoissonsolver_a_SOURCES = \
	Poisson_Solver.f90 

libpoissonsolver_a_LIBADD = 

PS_Check_SOURCES = PS_Check.f90 timing.f90 base.f90 $(mpi_source)
PS_Check_LDFLAGS = $(LMPI)
PS_Check_LDADD = libpoissonsolver.a \
	ABINIT-common/libABINIT.a \
	@MPI_LIBS@

PS_Program_SOURCES = PS_Program.f90 timing.f90 base.f90 $(mpi_source)
PS_Program_LDFLAGS = $(LMPI)
PS_Program_LDADD = libpoissonsolver.a \
	ABINIT-common/libABINIT.a \
	@MPI_LIBS@

PS_Exercise_SOURCES = PS_Exercise.f90 timing.f90 base.f90 
PS_Exercise_LDADD = libpoissonsolver.a \
	ABINIT-common/libABINIT.a


CLEANFILES = mpif.h POISSON_SOLVER.@MODULE_EXT@ poisson_solver.@MODULE_EXT@ \
	module_base.@MODULE_EXT@

#dependencies
mpif.h:
	touch mpif.h && \
	 echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h && \
	 echo "integer :: MPI_DOUBLE_PRECISION,MPI_INTEGER" >> mpif.h && \
	 echo "integer :: MPI_MAX,MPI_MIN" >> mpif.h

poisson_solver.@MODULE_EXT@ POISSON_SOLVER.@MODULE_EXT@: Poisson_Solver.o

Poisson_Solver.o : base.o $(mpi_include) $(SOURCES_MOD)

PS_Check.o PSolver.o: Poisson_Solver.o

