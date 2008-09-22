EXTRA_DIST = \
	structUtil.h \
	convolution.h \
	convolutionPot.h \
	convSeria.h \
	conv_shared_multi_kernel.cu \
	conv_pot_shared_kernel.cu \
	newcase.cu \
	Makefile.standalone \
	read_conf_file.h

SOURCES_CU_CONV1 = convolution.cu \
	convolutionPot.cu \
	structUtil.cu \
	cudafct.cu

SOURCES_CXX_CONV1 = convSeria.cpp \
	setCPU_GPU_aff.cpp \
	read_conf_file.cpp \
	fortran.c

INC_CUDA = -I@LIB_CUTIL_PATH@

noinst_LIBRARIES = libGPU.a

libGPU_a_SOURCES = $(SOURCES_CU_CONV1) $(SOURCES_CXX_CONV1)
AM_CPPFLAGS = $(INC_CUDA) -I. -I@CUDA_PATH@/include/ 
AM_CFLAGS = $(INC_CUDA) -I. -I@CUDA_PATH@/include/ -DCUBLAS_USE_THUNKING 

NVCC = @NVCC@
CUDA_FLAGS = @INC_GCC41CUDA@ @CUDA_FLAGS@
SUFFIXES = .cu

if USE_CUDA_GPU
sub_dirs_cuda=CUDA
libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas libGPU.a @LIB_CUTIL_PATH@/libcutil.a

check_PROGRAMS = conv_check

else
check_PROGRAMS = 
EXTRA_DIST += conv_check.f90 \
	1Dconv.cu
endif

CLEANFILES = module_base.@MODULE_EXT@

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @XC_INCLUDE@ 

conv_check_SOURCES = conv_check.f90 1Dconv.cu ../profiling/memory.f90 ../modules/base.f90
conv_check_LDFLAGS =
conv_check_LDADD = ../convolutions/libconvolutions.a \
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
	structUtil.h \
	newcase.cu 


conv_shared_multi.cu : structUtil.h


convolutionPot.cu : conv_pot_shared_kernel.cu \
	convolutionPot.h \
	convSeria.h \
	structUtil.h

conv_pot_shared_kernel.cu : structUtil.h

conv_check.o : 1Dconv.cu


convSeria.cpp structUtil.cu : structUtil.h

read_conf_file.cpp : read_conf_file.h
