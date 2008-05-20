EXTRA_DIST = \
	structUtil.h \
	convolution.h \
	convSeria.h \
	conv_shared_multi_kernel.cu

SOURCES_CU_CONV1 = convolution.cu \
	structUtil.cu

SOURCES_CXX_CONV1 = convSeria.cpp

INC_CUDA = -I@LIB_CUTIL_PATH@

noinst_LIBRARIES = libGPU.a

libGPU_a_SOURCES = $(SOURCES_CU_CONV1) $(SOURCES_CXX_CONV1)
AM_CPPFLAGS = $(INC_CUDA) -I.

NVCC = @NVCC@
CUDA_FLAGS = @INC_GCC41CUDA@ @CUDA_FLAGS@
SUFFIXES = .cu
.cu.o:
	$(NVCC) $(CUDA_FLAGS) $(AM_CPPFLAGS) -c $< -o $@


#dependencies
convolution.cu : conv_shared_multi_kernel.cu \
	convolution.h \
	convSeria.h \
	structUtil.h

conv_shared_multi.cu : structUtil.h

convSeria.cpp structUtil.cu : structUtil.h
