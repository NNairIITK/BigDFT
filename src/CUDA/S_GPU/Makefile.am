S_GPU_SOURCES = src/bindings/C/s_gpu_C.cpp \
	src/fct_call_implementation.cpp  \
	src/localqueu.cpp \
	src/init_network.cpp \
	src/manage_cpu_affinity.cpp \
	src/read_conf_file.cpp \
	src/set_repartition.cpp \
	src/bindings/FORTRAN/s_gpu_F.c



S_GPU_SOURCES_CUDA = src/cudafct.cu  \
	src/memory_management.cu


INCLUDE_S_GPU = -I$(top_srcdir)/src/CUDA/S_GPU/include/  -I$(top_srcdir)/src/CUDA/S_GPU/src/

noinst_LIBRARIES = libs_gpu.a

libs_gpu_a_SOURCES = $(S_GPU_SOURCES) $(S_GPU_SOURCES_CUDA)
AM_CPPFLAGS = $(INCLUDE_S_GPU) 
AM_CFLAGS = $(INCLUDE_S_GPU) 


NVCC = @NVCC@
CUDA_FLAGS = @INC_GCC41CUDA@ -arch sm_13 @CUDA_FLAGS@ 
SUFFIXES = .cu

sub_dirs_cuda=CUDA

.cu.o:
	$(NVCC) $(CUDA_FLAGS) $(AM_CPPFLAGS) -c $< -o $@
