EXTRA_DIST = \
	Kinetic.h \
	MagicFilter.h \
	Uncompress.h \
	Wavelet.h \
	Initialize.h \
	OpenCL_wrappers.h \
	Stream.h \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

check_PROGRAMS = conv_check
noinst_PROGRAMS = bench_conv bench_conv_quick stream_check

noinst_LIBRARIES = libbasef.a libbasec.a libOCL.a libStream.a

libbasef_a_SOURCES = ../profiling/memory.f90 ../modules/base.F90 ../modules/cublas_fake.f90 $(mpi_source)
libbasec_a_SOURCES = OpenCL_wrappers.c MagicFilter.c Kinetic.c Wavelet.c Uncompress.c Initialize.c Profiling.c Hamiltonian.c
libOCL_a_SOURCES = $(libbasec_a_SOURCES)
libStream_a_SOURCES = Stream.c

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir)

conv_check_SOURCES = conv_check.f90 
conv_check_LDFLAGS = -lOpenCL
conv_check_LDADD = libbasec.a libbasef.a ../libbigdft.a \
	../PSolver/libpoissonsolver.a
stream_check_SOURCES = StreamCheck.c
stream_check_LDFLAGS = -lOpenCL
stream_check_LDADD = libOCL.a libStream.a
bench_conv_SOURCES = bench_conv.c
bench_conv_LDFLAGS = -lOpenCL
bench_conv_LDADD = libbasec.a
bench_conv_quick_SOURCES = bench_conv_quick.c
bench_conv_quick_LDFLAGS = -lOpenCL
bench_conv_quick_LDADD = libbasec.a

stream: stream_check

bench: bench_conv bench_conv_quick
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='profile.conf' $(top_builddir)/src/OpenCL/bench_conv_quick
	OPENCL_PROFILE=1 OPENCL_PROFILE_CONFIG='profile.conf' $(top_builddir)/src/OpenCL/bench_conv

check:
	echo 1 64 64 2048 2048 1000 > fort.1
	$(top_builddir)/src/OpenCL/conv_check
