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

if BUILD_LIBXC
LIBXC_MODULE =
xc.o: base.o
else
LIBXC_MODULE = libxc_fake.f90
xc.o: base.o libxc_fake.o
endif

if CAPITALIZE
mod_bigdft = MODULE_BASE.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@ \
	MODULE_XC.@MODULE_EXT@ \
	OP2P_MODULE.@MODULE_EXT@ \
	MODULE_INPUT.@MODULE_EXT@
else
mod_bigdft = module_base.@MODULE_EXT@ \
	module_defs.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ \
	module_xc.@MODULE_EXT@ \
	op2p_module.@MODULE_EXT@ \
	module_input.@MODULE_EXT@
endif
module_DATA = $(mod_bigdft)

EXTRA_DIST = configure.inc.in

modules_source = base.f90 \
	defs.F90 \
	types.f90 \
	input.f90 \
	interfaces.f90 \
	xc.f90 \
	op2p_module.f90 \
	$(LIBXC_MODULE) \
	$(GPUBLAS_INTERFACE)

AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

noinst_LIBRARIES = libmodules.a

libmodules_a_SOURCES = $(modules_source)

CLEANFILES = mpif.h $(mod_bigdft)

# Hack for XLF strange DEFS option.
PPFCCOMPILE = $(FC) $(DEFS:-D%=@FCDEFS@%) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "integer, parameter :: MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1, MPI_REAL8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=1, MPI_STATUS_SIZE=1" >> mpif.h
	echo "integer, parameter :: MPI_LAND=1" >> mpif.h


#dependencies
defs.o: $(mpi_include) configure.inc 
base.o: defs.o
types.o: base.o
input.o: base.o
interfaces.o: base.o types.o
op2p_module.o : base.o

$(mod_bigdft) : base.o types.o interfaces.o
