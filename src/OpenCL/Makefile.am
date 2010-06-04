EXTRA_DIST = \
	Kinetic.h \
	MagicFilter.h \
	Uncompress.h \
	Wavelet.h \
	Initialize.h \
	OpenCL_wrappers.h \
	Stream.h \
	Reduction.h \
        Tool.h\
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

check_PROGRAMS = conv_check
noinst_PROGRAMS = bench_conv bench_conv_quick stream_check

noinst_LIBRARIES = libbasef.a libbasec.a libOCL.a libStream.a

libbasef_a_SOURCES = ../profiling/memory.f90 \
                     ../modules/base.F90 \
                     ../modules/cublas_fake.f90 $(mpi_source)
libbasec_a_SOURCES = OpenCL_wrappers.c \
                     MagicFilter.c \
                     Kinetic.c \
                     Wavelet.c \
                     Uncompress.c \
                     Initialize.c \
                     Profiling.c \
                     Hamiltonian.c \
                     Reduction.c \
                     Preconditioner.c \
                     Density.c \
                     Benchmark.c
libOCL_a_SOURCES = $(libbasec_a_SOURCES)
libStream_a_SOURCES = Stream.c

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir)
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = conv_check.f90 
conv_check_LDADD = libbasec.a libbasef.a ../libbigdft.a ../PSolver/libpoissonsolver.a
conv_check_LDFLAGS = -lrt
stream_check_SOURCES = StreamCheck.c
stream_check_LDADD = libOCL.a libStream.a
stream_check_LDFLAGS = -lrt
bench_conv_SOURCES = bench_conv.c bench_lib.c
bench_conv_LDADD = libbasec.a
bench_conv_LDFLAGS = -lrt
bench_conv_quick_SOURCES = bench_conv_quick.c bench_lib.c
bench_conv_quick_LDADD = libbasec.a
bench_conv_quick_LDFLAGS = -lrt

stream: stream_check

bench: bench_conv
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='profile.conf' ./bench_conv

bench_quick: bench_conv_quick
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='profile.conf' ./bench_conv_quick

check:
	echo 1 128 128 16384 16384 100 > fort.1
	./conv_check
