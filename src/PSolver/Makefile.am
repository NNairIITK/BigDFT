## Process this file with automake to produce Makefile.in

EXTRA_DIST = gaussfit_wires.inc\
	lazy_16.inc\
	lazy_40.inc\
	lazy_ISF_100_2048.inc\
	lazy_ISF_20_2048.inc\
	lazy_ISF_50_2048.inc\
	gaussfit_Yukawa.inc\
	lazy_20.inc\
	lazy_50.inc\
	lazy_ISF_14_2048.inc\
	lazy_ISF_24_2048.inc\
	lazy_ISF_60_2048.inc\
	lazy_100.inc\
	lazy_24.inc\
	lazy_60.inc\
	lazy_ISF_16_2048.inc\
	lazy_ISF_30_2048.inc\
	lazy_ISF_8_2048.inc\
	lazy_14.inc\
	lazy_30.inc\
	lazy_8.inc\
	lazy_ISF_16_2048.preinc\
	lazy_ISF_40_2048.inc\
	perfdata.inc\
	createKernel.f90 \
	PSolver_Main.f90


SOURCES_MOD = PSolver_Main.f90 \
	createKernel.f90

noinst_LIBRARIES = libpoissonsolver.a

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

libpoissonsolver_a_SOURCES = Poisson_Solver.f90

if CAPITALIZE
module_DATA = POISSON_SOLVER.@MODULE_EXT@
else
module_DATA = poisson_solver.@MODULE_EXT@
endif

CLEANFILES = $(module_DATA)

#dependencies
poisson_solver.@MODULE_EXT@ POISSON_SOLVER.@MODULE_EXT@: Poisson_Solver.o

Poisson_Solver.o: ../modules/base.o ../modules/xc.o $(SOURCES_MOD)

Build_Kernel.o: fft3d.o
