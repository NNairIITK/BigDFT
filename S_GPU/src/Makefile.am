## Process this file with automake to produce Makefile.in

SUBDIRS = CUDA

S_GPU_SOURCES = bindings/C/s_gpu_C.cpp \
	fct_call_implementation.cpp  \
	localqueu.cpp \
	init_network.cpp \
	manage_cpu_affinity.cpp \
	read_conf_file.cpp \
	set_repartition.cpp \
	bindings/FORTRAN/s_gpu_F.c \
	class_utils.h \
	exceptions.h \
	fct_call.h \
	fct_call_implementation.h \
	init_network.h \
	localqueu.h \
	manage_cpu_affinity.h \
	manage_global_var.h \
	message.h \
	read_conf_exception.h \
	read_conf_file.h \
	set_repartition.h \
	sg_common_def.h \
	s_gpu.h

include_HEADERS = sg_common_def.h \
	s_gpu.h

#S_GPU_OBJECTS		= $(S_GPU_SOURCES:.c*=.o) $(S_GPU_SOURCES_CUDA:.cu=.o)
#lib_LIBRARIES = libs_gpu.a
lib_LTLIBRARIES = libs_gpu.la
libs_gpu_la_LIBTOOLFLAGS= $(AM_LIBTOOLFLAGS) --tag=CUDA
libs_gpu_la_LDFLAGS= @AC_LDFLAGS_CUDA@
#libs_gpu_a_SOURCES = $(S_GPU_SOURCES) $(S_GPU_SOURCES_CUDA)
libs_gpu_la_SOURCES = $(S_GPU_SOURCES) $(S_GPU_SOURCES_CUDA)
libs_gpu_la_LIBADD  = CUDA/libcuda.la

#INCLUDE_S_GPU = -I$(top_srcdir)/include/  -I$(top_srcdir)/src/
#AM_CPPFLAGS = $(INCLUDE_S_GPU)
#AM_CFLAGS = $(INCLUDE_S_GPU)

#CXXLINK = @NVCC@ -shared
#CXXLD = @NVCC@

#NVCC = @NVCC@
#RANLIB = @NVCC@ -lib


#NVCCFLAGS=@CUDA_FLAGS@

#libs_gpu.so:
#	$(NVCC) @CUDA_FLAGS@ -shared --compiler-options '-fPIC' $(S_GPU_OBJECTS) -o libs_gpu.so -lcudart
