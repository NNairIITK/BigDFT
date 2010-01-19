EXTRA_DIST = \
	Makefile.standalone \
	GPUparameters.h \
	kernels_compress.hcu \
	kernels_locpot.hcu \
	kernels_kinetic.hcu \
	kernels_anasyn.hcu \
	structDef_anasyn.h \
	structDef_locpot.h \
	structDef_kinetic.h \
	check_cuda.h \
	kinetic.h \
	locpot.h \
	commonDef.h \
	read_conf_exception.h \
	cpp_utils.h \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV = 	cudafct.cu \
	locham.cu 

SOURCES_CONV_CHECK = locpot.cu \
	kinetic.cu \
	anasyn.cu \
	compress.cu \
	precond.cu \
	density.cu

SOURCES_SHARED = binding_s_gpu.c \
	commonDef.cpp  \
	fortran.c

INC_CUDA = -I@LIB_CUTIL_PATH@

noinst_LIBRARIES = libGPU.a

libGPU_a_SOURCES = $(SOURCES_CU_CONV1) $(SOURCES_SHARED)  $(SOURCES_CONV)
AM_CPPFLAGS = $(INC_CUDA) -I. -I$(top_builddir) -I@CUDA_PATH@/include/ 
AM_CFLAGS = $(INC_CUDA) -I. -I@CUDA_PATH@/include/ -DCUBLAS_USE_THUNKING 

NVCC = @NVCC@
CUDA_FLAGS = @INC_GCC41CUDA@ -arch sm_13 @CUDA_FLAGS@ 
SUFFIXES = .cu

if USE_CUDA_GPU
sub_dirs_cuda=CUDA
libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas libGPU.a S_GPU/libs_gpu.a
#cutil is used for CUDA timer
#libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas libGPU.a @LIB_CUTIL_PATH@/libcutil.a

#if BUILD_BINARIES
#bin_PROGRAMS = GPUham
#else
#bin_PROGRAMS = 
#endif

check_PROGRAMS = conv_check 
#	GPUham

else
check_PROGRAMS = 
EXTRA_DIST += conv_check.f90 \
	GPUham.f90 \
	$(SOURCES_CONV)
endif

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir)

conv_check_SOURCES = conv_check.f90  ../profiling/memory.f90 ../modules/base.F90 $(mpi_source)
conv_check_LDFLAGS =
conv_check_LDADD = ../libbigdft.a \
	../PSolver/libpoissonsolver.a \
	$(libs_cuda)

#GPUham_SOURCES = GPUham.f90 ../profiling/memory.f90 ../modules/base.F90 $(mpi_source)
#GPUham_LDFLAGS =
#GPUham_LDADD = ../libbigdft.a \
#	$(libs_cuda)

.cu.o:
	$(NVCC) $(CUDA_FLAGS) $(AM_CPPFLAGS) -c $< -o $@


#dependencies

locpot.o: kernels_locpot.hcu structDef_locpot.h	commonDef.h

kinetic.o: kernels_kinetic.hcu structDef_kinetic.h commonDef.h

anasyn.o: kernels_anasyn.hcu structDef_anasyn.h commonDef.h

locham.o : kernels_anasyn.hcu \
	kernels_locpot.hcu \
	kernels_kinetic.hcu \
	kernels_compress.hcu \
	GPUparameters.h \
	precond.cu \
	density.cu \
	commonDef.h

compress.o : commonDef.h kernels_compress.hcu

read_conf_file.cpp : read_conf_file.h

check:
	echo 1 64 64 2048 2048 1 > fort.1
	$(top_builddir)/src/CUDA/conv_check
