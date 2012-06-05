## Process this file with automake to produce Makefile.in

EXTRA_DIST = *.inc \
	PSolver/createKernel.f90 \
	PSolver/PSolver_Main.f90


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
