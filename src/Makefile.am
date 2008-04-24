## Process this file with automake to produce Makefile.in

if BUILD_LIB_PSOLVER
sub_dirs_ps = PSolver
else
sub_dirs_ps = 
endif
if BUILD_LIB_BIGDFT
lib_bigdft = libbigdft.a
if CAPITALIZE
mod_bigdft = BIGDFT_API.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@
else
mod_bigdft = bigdft_api.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ \
	module_types.@MODULE_EXT@
endif
else
lib_bigdft =
mod_bigdft =
endif
if BUILD_BINARIES
binaries_sources = cluster test_forces memguess
else
binaries_sources = 
endif


SUBDIRS = $(sub_dirs_ps) convolutions profiling

EXTRA_DIST = \
	Makefile.standalone \
	wavelib/intots.inc \
	wavelib/recs16.inc \
	wavelib/sym_16.inc 

AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@ -IPSolver -IPSolver/ABINIT-common

bin_PROGRAMS = $(binaries_sources)

lib_LIBRARIES = $(lib_bigdft)
module_DATA = $(mod_bigdft)

if USE_MPI
mpi_source =
mpi_include =
LMPI = @MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
lMPI =
endif

if USE_OPTI_CONVOLUT
convolution_source = convolutions/combined_shrink_optim.f90 \
	convolutions/combined_grow_optim.f90 \
	convolutions/combined_shrink_optim.f90 \
	convolutions/convolut_ib_optim.f90 \
	convolutions/conv_per_optim.f90 \
	convolutions/convolut_simple.f90 \
	convolutions/convolution_prev.f90 

else
convolution_source = convolutions/combined_shrink_simple.f90 \
	convolutions/combined_grow_simple.f90 \
	convolutions/combined_shrink_simple.f90 \
	convolutions/convolut_ib_simple.f90 \
	convolutions/conv_per_simple.f90 \
	convolutions/convolut_simple.f90 \
	convolutions/convolution_prev.f90 

endif

convolution_objects = $(convolution_source:.f90=.o)

#objects which depends from module_types
high_level = cluster.f90 \
	input_variables.f90 \
	init.f90 \
	sumrho.f90 \
	hpsiortho.f90 \
	forces.f90 \
	tail.f90 \
	restart.f90 \
	geometry.f90

#objects which can be compiled whithout any dependence
low_level = init/projectors.f90  \
	init/ionicpot.f90 \
	init/gridmanipulation.f90 \
	init/gauprod.f90 \
	init/denspotd.f90 \
	init/logrid.f90 \
	init/inputguess.f90 \
	init/gautowav.f90 \
	init/eleconf.f90 \
	wfn_opt/precond.f90 \
	wfn_opt/orthogonality.f90 \
	wfn_opt/diis.f90 \
	wfn_opt/applyh.f90 \
	wfn_opt/applyh_prev.f90 \
	wfn_opt/precond_per.f90 \
	wfn_opt/kernel_per.f90 \
	wfn_opt/fft_per.f90 \
	wavelib/transwaves.f90 \
	wavelib/scalprod.f90 \
	wavelib/plotting.f90 \
	wavelib/i-o.f90 \
	wavelib/un-compress.f90 \
	wavelib/razero.f90 \
	wavelib/gauss_to_daub.f90 

#common routine with abinit
abinit = ABINIT-ewald/ewald.F90 \
	 ABINIT-ewald/metric.F90 \
	 ABINIT-ewald/matr3inv.F90 \
	 ABINIT-ewald/derfc.F90 

libbigdft_a_SOURCES = \
	$(mpi_source) \
	types.f90 \
	interfaces.f90 \
	BigDFT_API.f90 \
	$(high_level) \
	$(low_level) \
	$(abinit)

libbigdft_a_LIBADD = profiling/timem.o \
	profiling/memoryestimator.o \
	$(convolution_objects)

cluster_SOURCES = BigDFT.f90 
cluster_LDFLAGS = $(LMPI)
cluster_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a \
	@MPI_LIBS@

test_forces_SOURCES = test_forces.f90
test_forces_LDFLAGS = $(LMPI)
test_forces_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a \
	@MPI_LIBS@

memguess_SOURCES = memguess.f90 
memguess_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a 


CLEANFILES = mpif.h \
	bigdft_api.@MODULE_EXT@ BIGDFT_API.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ MODULE_INTERFACES.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ MODULE_TYPES.@MODULE_EXT@

#dependencies
cluster.o geometry.o: $(mpi_include) \
	interfaces.o

init.o : interfaces.o

interfaces.o: types.o

BigDFT_API.o: types.o interfaces.o

$(high_level:.f90=.o) : types.o

mpif.h:
	touch mpif.h &&\
	 echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h &&\
	 echo "integer :: MPI_DOUBLE_PRECISION, MPI_INTEGER" >> mpif.h &&\
	 echo "integer :: MPI_MIN, MPI_MAX" >> mpif.h

test_forces.o: parameters.h \
	$(mod_bigdft)

$(mod_bigdft) : BigDFT_API.o

BigDFT.o: parameters.h \
	interfaces.o

memguess.o: types.o
