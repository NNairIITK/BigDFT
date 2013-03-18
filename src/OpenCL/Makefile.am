EXTRA_PROGRAMS = conv_check \
	conv_check_fft \
	bench_conv \
	bench_conv_quick \
	bench_zgemm \
	stream_check

CLEANFILES = $(EXTRA_PROGRAMS)

AM_FCFLAGS = -I../modules -I../PSolver @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ -I. -I$(srcdir)
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Wno-deprecated-declarations
AM_CXXFLAGS = -Wall -Wno-deprecated-declarations

conv_check_SOURCES = conv_check.f90
conv_check_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

conv_check_fft_SOURCES = conv_check_fft.f90
conv_check_fft_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

stream_check_SOURCES = Stream.c StreamCheck.c
stream_check_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

bench_conv_SOURCES = bench_conv.c bench_lib.c
bench_conv_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

bench_conv_quick_SOURCES = bench_conv_quick.c bench_lib.c
bench_conv_quick_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

bench_zgemm_SOURCES = bench_zgemm.c bench_lib.c
bench_zgemm_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

$(EXTRA_PROGRAMS): lib.stamp

lib.stamp:
	cd ..; $(MAKE) libbigdft.a libdft_common.a

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
