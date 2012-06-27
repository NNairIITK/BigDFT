if USE_CUDA_GPU
GPUBLAS_INTERFACE = 
else
GPUBLAS_INTERFACE = cublas_fake.f90
endif

if USE_BLACS
BLACS = 
else
BLACS = blacs_fake.f90
endif

if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

if USE_LIBXC
LIBXC_MODULE =
xc.o: base.o
else
LIBXC_MODULE = libxc_fake.f90
xc.o: base.o libxc_fake.o
endif

if BUILD_BINDINGS
bindings_headers = bigdft.h bigdft_cst.h
else
bindings_headers =
endif

if CAPITALIZE
mod_bigdft = MODULE_BASE.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@ \
	MODULE_PRIVATE_API.@MODULE_EXT@ \
	MODULE_XC.@MODULE_EXT@ \
	OVERLAP_POINT_TO_POINT.@MODULE_EXT@ \
	MODULE_INPUT.@MODULE_EXT@ \
	BIGDFT_API.@MODULE_EXT@ \
	YAML_OUTPUT.@MODULE_EXT@ 
else
mod_bigdft = module_base.@MODULE_EXT@ \
	module_defs.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ \
	module_private_api.@MODULE_EXT@ \
	module_xc.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@ \
	module_input.@MODULE_EXT@ \
	bigdft_api.@MODULE_EXT@ \
	yaml_output.@MODULE_EXT@
endif
module_DATA = $(mod_bigdft)
include_HEADERS = $(bindings_headers)

EXTRA_DIST = configure.inc.in bigdft_cst.h.in

AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

noinst_LIBRARIES = libmodules.a

libmodules_a_SOURCES = base.f90 \
	defs.F90 \
	types.f90 \
	input.f90 \
	interfaces.f90 \
	private_api.f90 \
	xc.f90 \
	op2p_module.f90 \
	BigDFT_API.f90 \
	vdwcorrection.f90 \
	yaml_output.f90 \
	$(bindings_headers) \
	$(LIBXC_MODULE) \
	$(GPUBLAS_INTERFACE) \
	$(BLACS)

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
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=10, MPI_STATUS_SIZE=1,MPI_LAND=1" >> mpif.h &&\
        echo "integer, parameter :: MPI_REQUEST_NULL=1,MPI_STATUS_IGNORE=1" >> mpif.h &&\
        echo "integer, parameter :: mpi_tag_ub=1" >> mpif.h &&\
        echo "real(kind=8), external :: mpi_wtime" >> mpif.h


#dependencies
defs.o: $(mpi_include) configure.inc 
base.o: defs.o
types.o: base.o
input.o: base.o yaml_output.o
interfaces.o: base.o types.o xc.o
op2p_module.o : base.o
BigDFT_API.o: base.o types.o interfaces.o xc.o
vdwcorrection.o: base.o types.o
private_api.o: base.o types.o interfaces.o
yaml_output.o : base.o

$(mod_bigdft) : base.o types.o interfaces.o
