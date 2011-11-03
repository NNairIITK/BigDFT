## Process this file with automake to produce Makefile.in

SOURCES_MOD = PSolver_Main.f90 \
	createKernel.f90 \
	xcenergy.f90 \
	3Dgradient.f90 

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

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

libpoissonsolver_a_SOURCES = \
	Poisson_Solver.f90 \
	Build_Kernel.f90 \
	xcenergy.f90 \
	3Dgradient.f90 \
	scaling_function.f90 \
	PSolver_Base_new.f90 \
	wofz.f90 \
	../lib/fft/fft3d.f90 \
	../profiling/time.f90 

#libpoissonsolver_a_LIBADD = ../modules/libmodules.a


CLEANFILES = malloc.prc time.prc \
	module_fft_sg.@MODULE_EXT@ MODULE_FFT_SG.@MODULE_EXT@ \
	POISSON_SOLVER.@MODULE_EXT@ poisson_solver.@MODULE_EXT@ \
	module_base.@MODULE_EXT@ timedata.@MODULE_EXT@ \
	timedata.@MODULE_EXT@ TIMEDATA.@MODULE_EXT@

#dependencies
poisson_solver.@MODULE_EXT@ POISSON_SOLVER.@MODULE_EXT@: Poisson_Solver.o

Poisson_Solver.o: ../modules/base.o $(SOURCES_MOD)

Build_Kernel.o: fft3d.o

xcenergy.o : Poisson_Solver.o
