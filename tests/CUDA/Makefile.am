
if USE_CUDA_GPU
check_PROGRAMS = conv_check 
else
check_PROGRAMS = 
endif

if USE_MPI
mpi_source =
else
mpi_source = ../../src/MPIfake.f90
endif

if USE_OCL
libs_ocl=../../src/OpenCL/libOCL.a
else
libs_ocl=
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../../src/modules @MPI_INCLUDE@ -I. -I$(srcdir)

conv_check_SOURCES = ../../src/CUDA/conv_check.f90  ../../src/modules/defs.F90 ../../src/modules/base.f90 $(mpi_source)
conv_check_LDFLAGS = 
conv_check_LDADD = ../../src/libbigdft.a \
	../../src/PSolver/libpoissonsolver.a \
	../../src/CUDA/libGPU.a $(libs_ocl)

check:
	echo 1 64 64 2048 2048 2 > fort.1
	./conv_check
