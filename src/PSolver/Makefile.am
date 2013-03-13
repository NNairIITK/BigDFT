if USE_CUDA_GPU
EXTRA_PROGRAMS = conv_check_fft

conv_check_fft_SOURCES = conv_check_fft.f90
conv_check_fft_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

conv_check_fftw_SOURCES = conv_check_fftw.f90 \
	fftw.cpp
conv_check_fftw_LDADD = @LIBS_DEPENDENCIES@ -lfftw3 -lfftw3_threads

else
EXTRA_PROGRAMS = conv_check_fftw

conv_check_fft_SOURCES =
conv_check_fft_LDADD =

conv_check_fftw_SOURCES =
conv_check_fftw_LDADD =

endif

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

CLEANFILES = $(EXTRA_PROGRAMS) $(module_DATA)

SUFFIXES = .cu

NVCC = @NVCC@
NVCC_FLAGS = -arch sm_13 @NVCC_FLAGS@

AM_CPPFLAGS = -I. -I$(top_builddir) -I@CUDA_PATH@/include/

#dependencies
poisson_solver.@MODULE_EXT@ POISSON_SOLVER.@MODULE_EXT@: Poisson_Solver.o

Poisson_Solver.o: ../modules/base.o ../modules/xc.o $(SOURCES_MOD)

Build_Kernel.o: fft3d.o

conv_check_fft.o: cufft.o

conv_check_fftw.o: fftw.o

.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@

.cpp.o:
	icc -openmp -lfftw3 -lfftw3_threads -c $< -o $@
