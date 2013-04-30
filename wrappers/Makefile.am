# To be installed:
noinst_LIBRARIES = libwrappers.a
module_DATA = $(modules)

# Source code.
if USE_DGEMMSY
DGEMMSY_SOURCES = interface_dgemmsy.f90 \
	dgemmsy/dgemmsy.c       dgemmsy/dgemmsy.h \
	dgemmsy/gemm_block_c.c  dgemmsy/gemm_block_c.h \
	dgemmsy/patterns.c      dgemmsy/patterns.h \
	dgemmsy/dgemmsy_utils.c dgemmsy/utils.h \
	dgemmsy/visitors.c      dgemmsy/visitors.h
else
DGEMMSY_SOURCES = interface_dgemmsy_fake.f90
endif
if USE_CUDA_GPU
CUBLAS_SOURCES = 
else
CUBLAS_SOURCES = cublas_fake.f90 interface_cuda_fake.f90
endif
if USE_OCL
OCL_SOURCES =
else
OCL_SOURCES = interface_ocl_fake.f90
endif
libwrappers_a_SOURCES = $(DGEMMSY_SOURCES) $(CUBLAS_SOURCES) $(OCL_SOURCES) \
	fft/fft2d.f90 \
	fft/fft3d.f90 \
	linalg.f90 \
	mpi.F90

# Dependencies
$(modules): $(libwrappers_a_OBJECTS)
interface_cuda_fake.o: linalg.o
interface_dgemmsy_fake.o: linalg.o

AM_FCFLAGS = -I$(top_builddir)/includes @LIBABINIT_INCLUDE@
AM_CPPFLAGS = @DGEMMSY_CPPFLAGS@

# Module files.
if CAPITALIZE
modules = MODULE_FFT_SG.@MODULE_EXT@ \
	WRAPPER_LINALG.@MODULE_EXT@ \
	WRAPPER_MPI.@MODULE_EXT@
else
modules = module_fft_sg.@MODULE_EXT@ \
	wrapper_linalg.@MODULE_EXT@ \
	wrapper_mpi.@MODULE_EXT@
endif
CLEANFILES = $(modules)

# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done