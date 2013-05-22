# To be installed:
noinst_LIBRARIES = libwrappers.a
module_DATA = $(modules)

if USE_MPI
#mpi_source =
mpi_include =
else
#mpi_source = MPIfake.f90
mpi_include = mpif.h
endif

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "integer, parameter :: MPI_SUCCESS=0,MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1, MPI_REAL8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=10, MPI_STATUS_SIZE=1,MPI_LAND=1" >> mpif.h &&\
        echo "integer, parameter :: MPI_REQUEST_NULL=1,MPI_STATUS_IGNORE=1" >> mpif.h &&\
        echo "integer, parameter :: mpi_tag_ub=1,mpi_address_kind=8,mpi_info_null=0" >> mpif.h &&\
        echo "integer, parameter :: mpi_mode_noprecede=0" >> mpif.h &&\
        echo "real(kind=8), external :: mpi_wtime" >> mpif.h


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
CUBLAS_SOURCES = fortran.c
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
mpi.o: $(mpi_include)


AM_FCFLAGS = -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @MPI_INCLUDE@
AM_CPPFLAGS = @DGEMMSY_CPPFLAGS@
if USE_CUDA_GPU
AM_CPPFLAGS += -I@CUDA_PATH@/include
endif

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
