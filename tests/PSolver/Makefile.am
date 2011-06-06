## Process this file with automake to produce Makefile.in

OTHERS = PS_Check.f90 \
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


check_PROGRAMS = PS_Check PS_Program PS_Exercise

AM_FCFLAGS = -I../../src/modules -I../../src/PSolver @MPI_INCLUDE@ -I. @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@


if USE_MPI
mpi_source =
LMPI = @MPI_LDFLAGS@ 
else
mpi_source = ../../src/MPIfake.f90
LMPI = 
endif

if USE_CUDA_GPU
libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas ../../src/CUDA/libGPU.a 
else
libs_cuda=
endif

if USE_OCL
libs_ocl=../../src/OpenCL/libOCL.a -lstdc++
else
libs_ocl=
endif

if USE_DGEMMSY
libs_dgemmsy=../../src/dgemmsy/libdgemmsy.a
else
libs_dgemmsy=
endif

libs_opt=$(libs_cuda)\
 $(libs_ocl)\
 $(libs_dgemmsy)


PS_Check_SOURCES = PS_Check.f90 $(mpi_source)
PS_Check_LDFLAGS =
PS_Check_LDADD = ../../src/PSolver/libpoissonsolver.a \
				 ../../src/libbigdft.a \
				 $(LMPI) \
				 @MPI_LIBS@ \
				 $(libs_opt)


PS_Program_SOURCES = PS_Program.f90 $(mpi_source)
PS_Program_LDFLAGS = 
PS_Program_LDADD = ../../src/PSolver/libpoissonsolver.a \
				   ../../src/libbigdft.a \
				   $(LMPI) \
				   @MPI_LIBS@ \
				   $(libs_opt)

PS_Exercise_SOURCES = PS_Exercise.f90 $(mpi_source)
PS_Exercise_LDADD =	../../src/PSolver/libpoissonsolver.a \
				    ../../src/libbigdft.a \
					$(LMPI) \
					@MPI_LIBS@ \
					$(libs_opt)

CLEANFILES = malloc.prc time.prc PS_Check.out fldiff.report \
	module_fft_sg.@MODULE_EXT@ MODULE_FFT_SG.@MODULE_EXT@ \
	POISSON_SOLVER.@MODULE_EXT@ poisson_solver.@MODULE_EXT@ module_base.@MODULE_EXT@

PS_Check.out: PS_Check PS_Program PS_Exercise
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

check: PS_Check.out
	python $(top_srcdir)/tests/fldiff.py --mode=psolver --discrepancy=2.e-8 \
		PS_Check.out $(top_srcdir)/tests/PSolver/PS_Check.out.ref | tee fldiff.report
