## Process this file with automake to produce Makefile.in

if BUILD_LIB_ABINIT
sub_dirs_xc = ABINIT-common
else
sub_dirs_xc = 
endif


SOURCES_MOD = PSolver_Main.f90 \
	Build_Kernel.f90 \
	PSolver_Base.f90 \
	xcenergy.f90 \
	3Dgradient.f90 \
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
	PS_Check.out.ref \
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

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @XC_INCLUDE@ @LIBXC_INCLUDE@


if USE_MPI
mpi_source =
LMPI = @MPI_LDFLAGS@ 
else
mpi_source = ../MPIfake.f90
LMPI = 
endif

if USE_CUDA_GPU
libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas ../CUDA/libGPU.a @LIB_CUTIL_PATH@/libcutil.a
else
libs_cuda=
endif


libpoissonsolver_a_SOURCES = \
	Poisson_Solver.f90 \
	../lib/fft/fft3d.f90 \
	../profiling/time.f90 

libpoissonsolver_a_LIBADD = ../modules/libmodules.a

PS_Check_SOURCES = PS_Check.f90 $(mpi_source)
PS_Check_LDFLAGS =
PS_Check_LDADD = libpoissonsolver.a \
	../modules/libmodules.a \
	../profiling/libprofiling.a \
	ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)
PS_Check_DEPENDENCIES = Poisson_Solver.o

PS_Program_SOURCES = PS_Program.f90 $(mpi_source)
PS_Program_LDFLAGS = 
PS_Program_LDADD = libpoissonsolver.a \
	../modules/libmodules.a \
	../profiling/libprofiling.a \
	ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda) 

PS_Exercise_SOURCES = PS_Exercise.f90 $(mpi_source)
PS_Exercise_LDADD =	libpoissonsolver.a \
	../modules/libmodules.a \
	../profiling/libprofiling.a \
	ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)

CLEANFILES = malloc.prc time.prc PS_Check.out \
	POISSON_SOLVER.@MODULE_EXT@ poisson_solver.@MODULE_EXT@

#dependencies
poisson_solver.@MODULE_EXT@ POISSON_SOLVER.@MODULE_EXT@: Poisson_Solver.o

Poisson_Solver.o: ../modules/base.o $(SOURCES_MOD)

#PS_Check.o: Poisson_Solver.o



check:
	$(run_parallel) ./PS_Check 57 48 63  0 F  > PS_Check.out
	$(run_parallel) ./PS_Check 57 48 63  1 F >> PS_Check.out
	$(run_parallel) ./PS_Check 57 48 63 11 F >> PS_Check.out
	$(run_parallel) ./PS_Check 57 48 63 13 F >> PS_Check.out
	$(run_parallel) ./PS_Check 64 64 64  0 P >> PS_Check.out
	$(run_parallel) ./PS_Check 64 64 64  1 P >> PS_Check.out
	$(run_parallel) ./PS_Check 64 64 64 11 P >> PS_Check.out
	$(run_parallel) ./PS_Check 64 64 64 13 P >> PS_Check.out
	$(run_parallel) ./PS_Check 32 64 48  0 S >> PS_Check.out
	$(run_parallel) ./PS_Check 32 64 48  1 S >> PS_Check.out
	$(run_parallel) ./PS_Check 32 64 48 11 S >> PS_Check.out
	$(run_parallel) ./PS_Check 32 64 48 13 S >> PS_Check.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=2.e-8 PS_Check.out $(top_srcdir)/src/PSolver/PS_Check.out.ref
