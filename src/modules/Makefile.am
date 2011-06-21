if USE_CUDA_GPU
GPUBLAS_INTERFACE = 
else
GPUBLAS_INTERFACE = cublas_fake.f90
endif

if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

if CAPITALIZE
mod_bigdft = MODULE_BASE.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@
else
mod_bigdft = module_base.@MODULE_EXT@ \
	module_defs.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@
endif
module_DATA = $(mod_bigdft)

EXTRA_DIST = configure.inc.in

modules_source = base.f90 \
	defs.F90 \
	types.f90 \
	interfaces.f90 \
	$(GPUBLAS_INTERFACE)

AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@ @LIBABINIT_INCLUDE@

noinst_LIBRARIES = libmodules.a

libmodules_a_SOURCES = $(modules_source)

CLEANFILES = mpif.h MODULE_BASE.@MODULE_EXT@ module_base.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ module_defs.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@ module_types.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@ module_interfaces.@MODULE_EXT@
#dependencies
interfaces.o: base.o types.o
types.o: base.o
defs.o: $(mpi_include) configure.inc 
base.o: defs.o

$(mod_bigdft) : base.o types.o interfaces.o

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "integer, parameter :: MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1, MPI_LAND=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1, MPI_REAL8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=1, MPI_STATUS_SIZE=1" >> mpif.h
