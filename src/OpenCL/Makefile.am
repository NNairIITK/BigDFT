EXTRA_PROGRAMS = conv_check conv_check_fft bench_conv bench_conv_quick bench_zgemm stream_check

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

if USE_CUDA_GPU
libs_cuda= -lcudart -lcublas ../CUDA/libGPU.a 
else
libs_cuda=
endif

libs_opt=$(libs_cuda)

CLEANFILES = $(EXTRA_PROGRAMS)

AM_FCFLAGS = -I../modules -I../PSolver @MPI_INCLUDE@ -I. -I$(srcdir)
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = conv_check.f90 $(mpi_source)
conv_check_LDADD =  ../libbigdft.a ../PSolver/libpoissonsolver.a

conv_check_fft_SOURCES = conv_check_fft.f90 $(mpi_source)
conv_check_fft_LDADD = ../PSolver/libpoissonsolver.a ../libbigdft.a $(libs_opt) 

stream_check_SOURCES = Stream.c StreamCheck.c $(mpi_source)
stream_check_LDADD = ../libbigdft.a ../PSolver/libpoissonsolver.a

bench_conv_SOURCES = bench_conv.c bench_lib.c $(mpi_source)
bench_conv_LDADD = ../libbigdft.a ../PSolver/libpoissonsolver.a 

bench_conv_quick_SOURCES = bench_conv_quick.c bench_lib.c $(mpi_source)
bench_conv_quick_LDADD = ../libbigdft.a ../PSolver/libpoissonsolver.a 

bench_zgemm_SOURCES = bench_zgemm.c bench_lib.c $(mpi_source)
bench_zgemm_LDADD = ../libbigdft.a ../PSolver/libpoissonsolver.a 

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
