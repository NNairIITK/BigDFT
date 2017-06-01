# To be installed:
noinst_LIBRARIES = libOpenCL.a

# Code source.
if USE_OCL
OCL_SOURCES = \
	Stream.h Tool.h bench_lib.h \
	OpenCL_wrappers.c         OpenCL_wrappers.h \
	MagicFilter.c             MagicFilter.h \
	MagicFilter_Generator.cpp MagicFilter_Generator.h \
	Kinetic.c                 Kinetic.h \
	Kinetic_Generator.cpp     Kinetic_Generator.h \
	Kinetic_k.h \
	Kinetic_k_Generator.cpp   Kinetic_k_Generator.h \
	Wavelet.c                 Wavelet.h \
	Wavelet_Generator.cpp     Wavelet_Generator.h \
	Uncompress.c              Uncompress.h \
	Initialize.c              Initialize.h \
	Profiling.c \
	Hamiltonian.c \
	Reduction.c               Reduction.h \
	Reduction_Generator.cpp   Reduction_Generator.h \
	Preconditioner.c \
	Density.c \
	dsfft.c \
	fft_generator_noshared.cpp \
	fft_noshared.c \
	fft_generator_helper.cpp \
	dsfft_generator.cpp       fft_generator.h \
	Benchmark.c \
	Benchmark_Generator.cpp   Benchmark_Generator.h
else
OCL_SOURCES = fake_wrappers.c
endif
libOpenCL_a_SOURCES = $(OCL_SOURCES)

CLEANFILES = $(EXTRA_PROGRAMS)

EXTRA_PROGRAMS = conv_check \
	conv_check_fft \
	bench_conv \
	bench_conv_quick \
	bench_zgemm \
	stream_check

AM_FCFLAGS = -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ -I. -I$(srcdir)
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Wno-deprecated-declarations
AM_CXXFLAGS = -Wall -Wno-deprecated-declarations
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
bigdft_library = $(top_builddir)/src/libbigdft-1.so.@BIGDFT_MINOR_VERSION@
else
bigdft_library = $(top_builddir)/src/libbigdft-1.a @LIB_BIGDFT_DEPS@
endif

conv_check_SOURCES = conv_check.f90
conv_check_LDADD = $(bigdft_library) @LINALG_LIBS@

conv_check_fft_SOURCES = conv_check_fft.f90
conv_check_fft_LDADD = $(bigdft_library)

stream_check_SOURCES = Stream.c StreamCheck.c
stream_check_LDADD = $(bigdft_library)

bench_conv_SOURCES = bench_conv.c bench_lib.c
bench_conv_LDADD = $(bigdft_library)

bench_conv_quick_SOURCES = bench_conv_quick.c bench_lib.c
bench_conv_quick_LDADD = $(bigdft_library)

bench_zgemm_SOURCES = bench_zgemm.c bench_lib.c
bench_zgemm_LDADD = $(bigdft_library)

$(EXTRA_PROGRAMS): lib.stamp

lib.stamp:
	cd ..; $(MAKE) libbigdft-1.a

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
