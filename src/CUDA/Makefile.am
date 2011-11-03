EXTRA_DIST = Makefile.standalone

# Not used anymore, to be removed
#INC_CUDA = -I@LIB_CUTIL_PATH@
EXTRA_PROGRAMS = conv_check

AM_CPPFLAGS = -I. -I$(top_builddir) -I@CUDA_PATH@/include/ @LIBSGPU_INCLUDE@

NVCC = @NVCC@
NVCC_FLAGS = -arch sm_13 @NVCC_FLAGS@ 
SUFFIXES = .cu

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@

conv_check_SOURCES = conv_check.f90 timing.c
conv_check_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

#GPUham_SOURCES = GPUham.f90 ../profiling/memory.f90 ../modules/base.F90 $(mpi_source)
#GPUham_LDFLAGS =
#GPUham_LDADD = ../libbigdft.a \
#	$(libs_cuda)

.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@


#dependencies

locpot.o: kernels_locpot.hcu structDef_locpot.h	commonDef.h

kinetic.o: kernels_kinetic.hcu structDef_kinetic.h commonDef.h

anasyn.o: kernels_anasyn.hcu structDef_anasyn.h commonDef.h

locham.o : kernels_anasyn.hcu \
	kernels_locpot.hcu \
	kernels_kinetic.hcu \
	kernels_compress.hcu \
	GPUparameters.h \
	precond.cu \
	density.cu \
	commonDef.h

compress.o : commonDef.h kernels_compress.hcu

read_conf_file.cpp : read_conf_file.h

check_conv: conv_check
	echo 1 64 64 2048 2048 1000 > fort.1
	$(top_builddir)/src/CUDA/conv_check
