EXTRA_DIST = \
	bench_lib.h \
	Kinetic.h \
	Kinetic_Generator.h \
	Kinetic_k.h \
	Kinetic_k_Generator.h \
	bench_lib.h \
	MagicFilter.h \
	MagicFilter_Generator.h \
	Uncompress.h \
	Wavelet.h \
	Wavelet_Generator.h \
	Initialize.h \
	OpenCL_wrappers.h \
	Stream.h \
	Reduction.h \
	Reduction_Generator.h \
	Benchmark_Generator.h \
        Tool.h\
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =


EXTRA_PROGRAMS = conv_check conv_check_fft bench_conv bench_conv_quick bench_zgemm stream_check

noinst_LIBRARIES = libbasef.a libbasec.a libOCL.a libStream.a

libbasef_a_SOURCES = $(mpi_source)
libbasef_a_LIBADD = ../modules/defs.o ../modules/base.o \
                     ../modules/cublas_fake.o

libbasec_a_SOURCES = OpenCL_wrappers.c \
                     MagicFilter.c \
                     MagicFilter_Generator.cpp \
                     Kinetic.c \
                     Kinetic_Generator.cpp \
                     Kinetic_k_Generator.cpp \
                     Wavelet.c \
                     Wavelet_Generator.cpp \
                     Uncompress.c \
                     Initialize.c \
                     Profiling.c \
                     Hamiltonian.c \
                     Reduction.c \
                     Reduction_Generator.cpp \
                     Preconditioner.c \
                     Density.c \
                     fft.c \
                     fft_generator.cpp \
                     Benchmark.c \
                     Benchmark_Generator.cpp
libOCL_a_SOURCES = $(libbasec_a_SOURCES)
libStream_a_SOURCES = Stream.c

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = conv_check.f90 
conv_check_LDADD =  ../libbigdft.a ../PSolver/libpoissonsolver.a libbasec.a libbasef.a
conv_check_LDFLAGS = -lrt -lstdc++
conv_check_fft_SOURCES = conv_check_fft.f90 
conv_check_fft_LDADD =  ../libbigdft.a ../PSolver/libpoissonsolver.a libbasec.a libbasef.a
conv_check_fft_LDFLAGS = -lrt -lstdc++
stream_check_SOURCES = StreamCheck.c
stream_check_LDADD = libOCL.a libStream.a
stream_check_LDFLAGS = -lrt -lstdc++
bench_conv_SOURCES = bench_conv.c bench_lib.c
bench_conv_LDADD = libbasec.a
bench_conv_LDFLAGS = -lrt -lstdc++
bench_conv_quick_SOURCES = bench_conv_quick.c bench_lib.c
bench_conv_quick_LDADD = libbasec.a
bench_conv_quick_LDFLAGS = -lrt -lstdc++
bench_zgemm_SOURCES = bench_zgemm.c bench_lib.c
bench_zgemm_LDADD = libbasec.a
bench_zgemm_LDFLAGS = -lrt -lstdc++

stream: stream_check
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='$(srcdir)/profile.conf' ./stream_check

bench: bench_conv
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='$(srcdir)/profile.conf' ./bench_conv

bench_quick: bench_conv_quick
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='$(srcdir)/profile.conf' ./bench_conv_quick

bench_zg: bench_zgemm
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='$(srcdir)/profile.conf' ./bench_zgemm

check_conv: conv_check
	echo 124 132 130 100 > fort.1
	./conv_check

check_conv_fft: conv_check_fft
	echo 64 70 84 100 > fort.2
	./conv_check_fft
