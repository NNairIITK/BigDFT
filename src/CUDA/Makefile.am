EXTRA_DIST = \
	structUtil.h \
	convolution.h \
	convolutionPot.h \
	convSeria.h \
	conv_shared_multi_kernel.cu \
	conv_pot_shared_kernel.cu \
	Makefile.standalone \
	read_conf_file.h \
	GPUparameters.h \
	reduction.hcu \
	kernel_tex.cu \
	kernels_compress.hcu \
	kernels_locpot.hcu \
	kernels_kinetic.hcu \
	kernels_anasyn.hcu \
	structDef_anasyn.h \
	structDef_locpot.h \
	structDef_kinetic.h \
	kinetic.h \
	locpot.h \
	reduction.h \
	commonDef.h \
	$(SOURCES_CONV_CHECK)

SOURCES_CU_CONV1 = convolution.cu \
	convolutionPot.cu \
	structUtil.cu \
	cudafct.cu

SOURCES_CONV = 	memory_management.cu \
	reduction.cu \
	locham.cu 

SOURCES_CONV_CHECK = locpot.cu \
	kinetic.cu \
	anasyn.cu \
	compress.cu \
	precond.cu \
	density.cu

SOURCES_SHARED = convSeria.cpp \
	dynamicGPU/interface_stream.cpp \
	dynamicGPU/init_network.cpp \
	dynamicGPU/convolution_fct_call.cpp \
	dynamicGPU/localqueu.cpp \
	init_lib.cpp \
	read_conf_file.cpp \
	commonDef.cpp  \
	fortran.c

INC_CUDA = -I@LIB_CUTIL_PATH@

noinst_LIBRARIES = libGPU.a

libGPU_a_SOURCES = $(SOURCES_CU_CONV1) $(SOURCES_SHARED)  $(SOURCES_CONV)
AM_CPPFLAGS = $(INC_CUDA) -I. -I@CUDA_PATH@/include/ 
AM_CFLAGS = $(INC_CUDA) -I. -I@CUDA_PATH@/include/ -DCUBLAS_USE_THUNKING 

NVCC = @NVCC@
CUDA_FLAGS = @INC_GCC41CUDA@ -arch sm_13 @CUDA_FLAGS@ 
SUFFIXES = .cu

if USE_CUDA_GPU
sub_dirs_cuda=CUDA
libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas libGPU.a
#cutil is used for CUDA timer
#libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas libGPU.a @LIB_CUTIL_PATH@/libcutil.a

if BUILD_BINARIES
bin_PROGRAMS = GPUham
else
bin_PROGRAMS = 
endif


check_PROGRAMS = conv_check

else
check_PROGRAMS = 
EXTRA_DIST += conv_check.f90 \
	GPUham.f90 \
	1Dconv.cu \
	1Dconv_new.cu \
	$(SOURCES_CONV)
endif

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif


CLEANFILES = module_base.@MODULE_EXT@ fort.1

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @XC_INCLUDE@ 

conv_check_SOURCES = conv_check.f90  ../profiling/memory.f90 ../modules/base.f90 $(mpi_source)
conv_check_LDFLAGS =
conv_check_LDADD = ../libbigdft.a \
	$(libs_cuda)

GPUham_SOURCES = GPUham.f90 ../profiling/memory.f90 ../modules/base.f90 $(mpi_source)
GPUham_LDFLAGS =
GPUham_LDADD = ../libbigdft.a \
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

kinetic.o: reduction.hcu kernels_kinetic.hcu structDef_kinetic.h

anasyn.o: kernels_anasyn.hcu structDef_anasyn.h

locham.o : kernels_anasyn.hcu \
	kernels_locpot.hcu \
	kernels_kinetic.hcu \
	kernels_compress.hcu \
	reduction.hcu \
	GPUparameters.h \
	precond.cu \
	density.cu

compress.o : commonDef.h kernels_compress.hcu

convolutionPot.cu : conv_pot_shared_kernel.cu \
	convolutionPot.h \
	convSeria.h \
	structUtil.h

conv_pot_shared_kernel.cu : structUtil.h

conv_check.o : 1Dconv.cu


convSeria.cpp structUtil.cu : structUtil.h

read_conf_file.cpp : read_conf_file.h

check:
	echo 1 64 64 2048 2048 1 > fort.1
	$(top_builddir)/src/CUDA/conv_check
