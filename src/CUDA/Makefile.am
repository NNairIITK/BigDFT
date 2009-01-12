EXTRA_DIST = \
	structUtil.h \
	convolution.h \
	convolutionPot.h \
	convSeria.h \
	conv_shared_multi_kernel.cu \
	conv_pot_shared_kernel.cu \
	Makefile.standalone \
	read_conf_file.h \
	reduction.hcu \
	kernel_tex.cu

SOURCES_CU_CONV1 = convolution.cu \
	convolutionPot.cu \
	structUtil.cu \
	cudafct.cu

SOURCES_CONV = locpot.cu \
	kinetic.cu 

SOURCES_SHARED = convSeria.cpp \
	setCPU_GPU_aff.cpp \
	read_conf_file.cpp \
	memory_management.cu \
	reduction.cu \
	commonDef.cpp \
	fortran.c

INC_CUDA = -I@LIB_CUTIL_PATH@

noinst_LIBRARIES = libGPU.a

libGPU_a_SOURCES = $(SOURCES_CU_CONV1) $(SOURCES_SHARED) $(SOURCES_CONV)
AM_CPPFLAGS = $(INC_CUDA) -I. -I@CUDA_PATH@/include/ 
AM_CFLAGS = $(INC_CUDA) -I. -I@CUDA_PATH@/include/ -DCUBLAS_USE_THUNKING 

NVCC = @NVCC@
CUDA_FLAGS = @INC_GCC41CUDA@ @CUDA_FLAGS@
SUFFIXES = .cu

if USE_CUDA_GPU
sub_dirs_cuda=CUDA
libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas libGPU.a @LIB_CUTIL_PATH@/libcutil.a

check_PROGRAMS = conv_check GPUham

else
check_PROGRAMS = 
EXTRA_DIST += conv_check.f90 \
	GPUham.f90 \
	1Dconv.cu \
	1Dconv_new.cu \
	locpot.cu \
	kinetic.cu 
endif

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif


CLEANFILES = module_base.@MODULE_EXT@ fort.1

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @XC_INCLUDE@ 

conv_check_SOURCES = conv_check.f90 1Dconv.cu 1Dconv_new.cu ../profiling/memory.f90 ../modules/base.f90 $(mpi_source)
conv_check_LDFLAGS =
conv_check_LDADD = ../convolutions/libconvolutions.a \
	$(libs_cuda)

GPUham_SOURCES = GPUham.f90 ../profiling/memory.f90 ../modules/base.f90 $(mpi_source)
GPUham_LDFLAGS =
GPUham_LDADD = ../convolutions/libconvolutions.a \
	$(libs_cuda)

.cu.o:
	$(NVCC) $(CUDA_FLAGS) $(AM_CPPFLAGS) -c $< -o $@


#dependencies
convolution.o : conv_shared_multi_kernel.cu \
	convolution.h \
	convSeria.h \
	structUtil.h

1Dconv.o : conv_shared_multi_kernel.cu \
	conv_pot_shared_kernel.cu \
	convolution.h \
	convSeria.h \
	structUtil.h 

conv_shared_multi.cu : structUtil.h

locpot.o: reduction.hcu kernels_locpot.hcu structDef_locpot.h

convolutionPot.cu : conv_pot_shared_kernel.cu \
	convolutionPot.h \
	convSeria.h \
	structUtil.h

conv_pot_shared_kernel.cu : structUtil.h

conv_check.o : 1Dconv.cu


convSeria.cpp structUtil.cu : structUtil.h

read_conf_file.cpp : read_conf_file.h

check:
	echo 3 128 128 64 64 5 > fort.1
	$(top_builddir)/src/CUDA/conv_check
