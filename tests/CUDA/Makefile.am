
if USE_CUDA_GPU
check_PROGRAMS = conv_check
my_tests = CUDA.report
else
check_PROGRAMS = 
endif

if USE_MPI
mpi_source =
else
mpi_source = $(top_builddir)/src/MPIfake.f90
endif

if USE_OCL
libs_ocl = $(top_builddir)/src/OpenCL/libOCL.a
else
libs_ocl =
endif

if USE_DGEMMSY
libs_dgemmsy = $(top_builddir)/src/dgemmsy/libdgemmsy.a
else
libs_dgemmsy =
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../../src/modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@

conv_check_SOURCES = ../../src/CUDA/conv_check.f90  ../../src/CUDA/timing.c ../../src/modules/defs.F90 ../../src/modules/base.f90 $(mpi_source)
conv_check_LDFLAGS = -lrt
conv_check_LDADD = $(top_builddir)/src/libbigdft.a \
	$(top_builddir)/src/PSolver/libpoissonsolver.a \
	$(top_builddir)/src/CUDA/libGPU.a $(libs_ocl) $(libs_dgemmsy)

conv_check.out: conv_check
	echo 1 64 64 2048 2048 2 > fort.1
	./conv_check > conv_check.out

CUDA.report: conv_check.out
	awk -f $(srcdir)/conv_check.awk $< | tee CUDA.report

check: $(my_tests)
