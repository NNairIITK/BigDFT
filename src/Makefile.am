## Process this file with automake to produce Makefile.in

if BUILD_LIB_PSOLVER
sub_dirs_ps = PSolver
else
sub_dirs_ps = 
endif
if BUILD_LIB_BIGDFT
lib_bigdft = libbigdft.a
if CAPITALIZE
mod_bigdft = BIGDFT_API.@MODULE_EXT@ 
else
mod_bigdft = bigdft_api.@MODULE_EXT@ 
endif
else
lib_bigdft =
mod_bigdft =
endif
if BUILD_BINARIES
binaries_sources = cluster test_forces memguess NEB frequencies $(minima_hopping_sources)
scripts_sources = NEB_driver.sh
else
binaries_sources = 
scripts_sources = 
endif

if USE_CUDA_GPU
sub_dirs_cuda=CUDA
GPU_INTERFACE = convolutions/convolut_common_interface_cuda.f90
libs_cuda=-L@CUDA_PATH@/lib/ -lcudart -lcublas CUDA/libGPU.a @LIB_CUTIL_PATH@/libcutil.a
GPU_BLAS=
else
sub_dirs_cuda=
GPU_INTERFACE = convolutions/interface_cuda_fake.f90
libs_cuda=
GPU_BLAS = modules/cublas_fake.f90
endif

if USE_MINIMA_HOPPING
minima_hopping_sources = global
else
minima_hopping_sources =
endif
global_SOURCES = # Should be global.f90, but see below.
global_LDFLAGS = 
global_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)
# Hack to avoid global.f90 inclusion into make dist.
am_global_OBJECTS = global.$(OBJEXT)
global_LINK = $(FCLD) $(AM_CFLAGS) $(CFLAGS) $(global_LDFLAGS) \
	$(LDFLAGS) -o $@


SUBDIRS = modules convolutions $(sub_dirs_ps) profiling $(sub_dirs_cuda)

EXTRA_DIST = \
	Makefile.standalone \
	wavelib/intots.inc \
	wavelib/recs16.inc \
	wavelib/sym_16.inc \
	lib \
	NEB_driver.sh

AM_FCFLAGS = -I. -Imodules -I$(srcdir) @MPI_INCLUDE@ -IPSolver -IPSolver/ABINIT-common

bin_PROGRAMS = $(binaries_sources)
bin_SCRIPTS = $(scripts_sources)

lib_LIBRARIES = $(lib_bigdft)
module_DATA = $(mod_bigdft)

if USE_MPI
mpi_source =
mpi_include =
LMPI = @MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
LMPI =
endif

if USE_OPTI_CONVOLUT
convolution_source = convolutions/combined_shrink_optim.f90 \
	convolutions/combined_grow_optim.f90 \
	convolutions/growshrink.f90 \
	convolutions/growshrink_hyb_optim.f90 \
	convolutions/growshrink_hyb_common.f90 \
	convolutions/convolut_ib_optim.f90 \
	convolutions/conv_per_optim.f90 \
	convolutions/conv_per_common.f90 \
	convolutions/convolut_common_slab.f90 \
	convolutions/convolut_common_per.f90 \
	convolutions/convolut_optim_per.f90 \
	convolutions/convolut_new_per.f90 \
	convolutions/convolut_optim_slab.f90 \
	convolutions/convolut_simple_per_k.f90 \
	convolutions/convolut_simple_slab_k.f90 \
	$(GPU_INTERFACE)

else
convolution_source = convolutions/combined_shrink_simple.f90 \
	convolutions/combined_grow_simple.f90 \
	convolutions/growshrink.f90 \
	convolutions/growshrink_hyb_simple.f90 \
	convolutions/growshrink_hyb_common.f90 \
	convolutions/convolut_ib_simple.f90 \
	convolutions/conv_per_simple.f90 \
	convolutions/conv_per_common.f90 \
	convolutions/convolut_common_slab.f90 \
	convolutions/convolut_common_per.f90 \
	convolutions/convolut_simple_per.f90 \
	convolutions/convolut_new_per.f90 \
	convolutions/convolut_simple_slab.f90 \
	convolutions/convolut_simple_per_k.f90 \
	convolutions/convolut_simple_slab_k.f90 \
	$(GPU_INTERFACE)

endif

convolution_objects = $(convolution_source:.f90=.o)

#modules for types, interfaces and base variables
modules_source = modules/base.f90 \
	modules/types.f90 \
	modules/interfaces.f90 \
	$(GPU_BLAS)

modules_objects = $(modules_source:.f90=.o)

#objects which depends from module_types
high_level = cluster.f90 \
	input_variables.f90 \
	init.f90 \
	sumrho.f90 \
	hpsiortho.f90 \
	forces.f90 \
	tail.f90 \
	restart.f90 \
	geometry.f90 \
	davidson.f90 \
	vdwcorrection.f90 \
	lanczos_interface.f90 \
	lanczos_base.f90 \
	lanczos.f90 

#objects which can be compiled whithout any dependence
low_level = init/projectors.f90  \
	init/ionicpot.f90 \
	init/gridmanipulation.f90 \
	init/gauprod.f90 \
	init/denspotd.f90 \
	init/logrid.f90 \
	init/logrid_per.f90 \
	init/inputguess.f90 \
	init/gautowav.f90 \
	init/eleconf.f90 \
	init/locreg.f90 \
	init/sysprop.f90 \
	init/gaussians.f90 \
	init/wfn_init.f90 \
	init/xabsorber.f90 \
	wfn_opt/precond.f90 \
	wfn_opt/orthogonality.f90 \
	wfn_opt/diis.f90 \
	wfn_opt/applyh.f90 \
	wfn_opt/precond_per_optim.f90 \
	wfn_opt/precond_hyb.f90 \
	wfn_opt/precond_slab.f90 \
	wfn_opt/kernel_per_optim.f90 \
	wfn_opt/kernel_slab_simple.f90 \
	wfn_opt/fft_per.f90 \
	wfn_opt/fft2d.f90 \
	wavelib/transwaves.f90 \
	wavelib/scalprod.f90 \
	wavelib/plotting.f90 \
	wavelib/i-o.f90 \
	wavelib/un-compress.f90 \
	wavelib/razero.f90 \
	wavelib/gauss_to_daub.f90 \
	wavelib/daubisf.f90

if BUILD_LIB_ABINIT
abinit = ABINIT-ewald/ewald.F90 \
	ABINIT-ewald/matr3inv.F90 \
	ABINIT-ewald/derfc.F90 \
	ABINIT-geometry/sort_dp.F90 \
	ABINIT-geometry/mati3inv.F90 \
	ABINIT-geometry/canon9.F90 \
	ABINIT-geometry/ab6_symmetry_f90.f90 \
	ABINIT-geometry/bldgrpaf.F90 \
	ABINIT-geometry/bldgrp.F90 \
	ABINIT-geometry/chkdilatmx.F90 \
	ABINIT-geometry/chkgrp.F90 \
	ABINIT-geometry/chkorthsy.F90 \
	ABINIT-geometry/chkprimit.F90 \
	ABINIT-geometry/gensymshub4.F90 \
	ABINIT-geometry/gensymshub.F90 \
	ABINIT-geometry/gensymspgr.F90 \
	ABINIT-geometry/getptgroupma.F90 \
	ABINIT-geometry/getspinrot.F90 \
	ABINIT-geometry/gridgcart.F90 \
	ABINIT-geometry/holocell.F90 \
	ABINIT-geometry/metric.F90 \
	ABINIT-geometry/mkrdim.F90 \
	ABINIT-geometry/operat.F90 \
	ABINIT-geometry/prtspgroup.F90 \
	ABINIT-geometry/ptgmadata.F90 \
	ABINIT-geometry/smallprim.F90 \
	ABINIT-geometry/spgdata.F90 \
	ABINIT-geometry/strainsym.F90 \
	ABINIT-geometry/strconv.F90 \
	ABINIT-geometry/stresssym.F90 \
	ABINIT-geometry/sym2cart.F90 \
	ABINIT-geometry/symanal.F90 \
	ABINIT-geometry/symatm.F90 \
	ABINIT-geometry/symaxes.F90 \
	ABINIT-geometry/symbrav.F90 \
	ABINIT-geometry/symchk.F90 \
	ABINIT-geometry/symdet.F90 \
	ABINIT-geometry/symfind.F90 \
	ABINIT-geometry/symlist_bcc.F90 \
	ABINIT-geometry/symlist_fcc.F90 \
	ABINIT-geometry/symlist_others.F90 \
	ABINIT-geometry/symlist_prim.F90 \
	ABINIT-geometry/symmultsg.F90 \
	ABINIT-geometry/symplanes.F90 \
	ABINIT-geometry/symredcart.F90 \
	ABINIT-geometry/symrelrot.F90 \
	ABINIT-geometry/symsgcube.F90 \
	ABINIT-geometry/symsghexa.F90 \
	ABINIT-geometry/symsgmono.F90 \
	ABINIT-geometry/symsgortho.F90 \
	ABINIT-geometry/symsgtetra.F90 \
	ABINIT-geometry/symspgr.F90 \
	ABINIT-geometry/symzat.F90 \
	ABINIT-geometry/xredxcart.F90 \
	ABINIT-geometry/smpbz.F90 \
	ABINIT-geometry/symkpt.F90 \
	ABINIT-geometry/getkgrid.F90
else
abinit = 
endif

libbigdft_a_SOURCES = \
	$(abinit) \
	$(mpi_source) \
	BigDFT_API.f90 \
	$(high_level) \
	$(low_level)

libbigdft_a_LIBADD = profiling/time.o \
	profiling/memory.o \
	profiling/memoryestimator.o \
	$(convolution_objects) \
	$(modules_objects)

cluster_SOURCES = BigDFT.f90 
cluster_LDFLAGS = 
cluster_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)

test_forces_SOURCES = test_forces.f90
test_forces_LDFLAGS = 
test_forces_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)

NEB_SOURCES = NEB.f90 \
	$(mpi_source)
NEB_LDFLAGS = 
NEB_LDADD = input_variables.o \
	$(modules_objects) \
	profiling/memory.o \
	convolutions/interface_cuda_fake.o \
	$(libs_cuda)

memguess_SOURCES = memguess.f90 \
	$(mpi_source)
memguess_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)

frequencies_SOURCES = frequencies.f90  \
	$(mpi_source)
frequencies_LDFLAGS = 
frequencies_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	PSolver/ABINIT-common/libABINIT.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)


CLEANFILES = mpif.h \
	bigdft_api.@MODULE_EXT@ BIGDFT_API.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ MODULE_INTERFACES.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ MODULE_TYPES.@MODULE_EXT@ \
	module_base.@MODULE_EXT@ MODULE_BASE.@MODULE_EXT@ \
	vdwcorrection.@MODULE_EXT@ VDWCORRECTION.@MODULE_EXT@ \
	lanczos_base.@MODULE_EXT@ LANCZOS_BASE.@MODULE_EXT@ \
	lanczos_interface.@MODULE_EXT@ LANCZOS_INTERFACE.@MODULE_EXT@ \
	formats.@MODULE_EXT@ FORMATS.@MODULE_EXT@ \
	miscellany.@MODULE_EXT@ MISCELLANY.@MODULE_EXT@ \
	neb_routines.@MODULE_EXT@ NEB_ROUTINES.@MODULE_EXT@ \
	numeric.@MODULE_EXT@ NUMERIC.@MODULE_EXT@ \
	minimization.@MODULE_EXT@ MINIMIZATION.@MODULE_EXT@ \
	minimization_routines.@MODULE_EXT@ MINIMIZATION_ROUTINES.@MODULE_EXT@ \
	neb_variables.@MODULE_EXT@ NEB_VARIABLES.@MODULE_EXT@

#dependencies
$(libbigdft_a_SOURCES:.f90=.o) BigDFT.o test_forces.o  memguess.o global.o : $(mpi_include)  $(modules_objects)

cluster.o geometry.o: vdwcorrection.o

lanczos.o : lanczos_base.o lanczos_interface.o

#init.o : base.o interfaces.o

#BigDFT_API.o: base.o types.o interfaces.o

#$(high_level:.f90=.o) : base.o types.o

mpif.h:
	touch mpif.h &&\
	 echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h &&\
	 echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1, MPI_LOGICAL=1" >> mpif.h &&\
	 echo "integer :: MPI_STATUSES_IGNORE" >> mpif.h &&\
	 echo "integer :: MPI_MIN, MPI_MAX,MPI_CHARACTER" >> mpif.h

test_forces.o: $(mod_bigdft)

$(mod_bigdft) : BigDFT_API.o

$(modules_objects) : $(mpi_include)

#BigDFT.o: \
#	base.o \
#	interfaces.o

#memguess.o: base.o types.o
