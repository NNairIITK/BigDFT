# To be installed:
noinst_LIBRARIES = libPSolver.a
module_DATA = $(modules)

# Code source.
if BUILD_LIB_BIGDFT
PS_SOURCES = Build_Kernel.f90 \
	xcenergy.f90 \
	3Dgradient.f90 \
	scaling_function.f90 \
	PSolver_Base_new.f90 \
	wofz.f90
else
PS_SOURCES =
endif
if USE_CUDA_GPU
CUDA_SOURCES = cufft.cu
else
CUDA_SOURCES =
endif
SOURCES_MOD = PSolver_Main.f90 createKernel.f90
libPSolver_a_SOURCES = $(PS_SOURCES) $(CUDA_SOURCES) Poisson_Solver.f90

# Modules.
if CAPITALIZE
modules = POISSON_SOLVER.@MODULE_EXT@
else
modules = poisson_solver.@MODULE_EXT@
endif

CLEANFILES = $(modules)
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
	lazy_ISF_40_2048.inc \
	perfdata.inc \
	$(SOURCES_MOD)

# Build targets.
if USE_CUDA_GPU
CONV_CHECK = conv_check_fft
else
CONV_CHECK = conv_check_fftw
endif
EXTRA_PROGRAMS = $(CONV_CHECK)
conv_check_fft_SOURCES = conv_check_fft.f90
conv_check_fft_LDADD = -lbigdft-1

conv_check_fftw_SOURCES = conv_check_fftw.f90 fftw.cpp
conv_check_fftw_LDADD = @LIBS_SHORT_DEPS@ -lfftw3 -lfftw3_threads

AM_FCFLAGS = -I../includes -I. -I$(srcdir) @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@ @MPI_INCLUDE@
AM_CPPFLAGS = -I. -I$(top_builddir) -I@CUDA_PATH@/include/ 
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

SUFFIXES = .cu

NVCC = @NVCC@
NVCC_FLAGS = -arch sm_13 @NVCC_FLAGS@ 

#dependencies
$(modules): $(libPSolver_a_OBJECTS)
Poisson_Solver.o: ../modules/base.o ../modules/xc.o $(SOURCES_MOD)
xcenergy.o: Poisson_Solver.o
conv_check_fft.o: cufft.o
conv_check_fftw.o: fftw.o

.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@

.cpp.o:
	icc -openmp -lfftw3 -lfftw3_threads -c $< -o $@

# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules)
	test -e "../includes" || $(INSTALL) -d -m 755 ../includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod ../includes ; \
	done
