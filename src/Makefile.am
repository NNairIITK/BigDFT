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
binaries_sources = bigdft test_forces memguess splsad NEB frequencies MDanalysis rism oneatom sandbox $(minima_hopping_sources) bart abscalc
scripts_sources = NEB_driver.sh
else
binaries_sources = 
scripts_sources = 
endif

if USE_CUDA_GPU
sub_dirs_cuda=CUDA
GPU_INTERFACE = convolutions/convolut_common_interface_cuda.f90
libs_cuda=CUDA/libGPU.a
GPU_BLAS=
else
sub_dirs_cuda=
GPU_INTERFACE = convolutions/interface_cuda_fake.f90
libs_cuda=
GPU_BLAS = modules/cublas_fake.f90
endif

if USE_OCL
sub_dirs_ocl=OpenCL
OCL_INTERFACE = convolutions/convolut_common_interface_ocl.f90
libs_ocl=OpenCL/libOCL.a
else
sub_dirs_ocl=
OCL_INTERFACE = convolutions/interface_ocl_fake.f90
libs_ocl=
endif

if USE_DGEMMSY
sub_dirs_dgemmsy=dgemmsy
DGEMMSY_INTERFACE = dgemmsy/interface_dgemmsy.f90
libs_dgemmsy=dgemmsy/libdgemmsy.a
else
sub_dirs_dgemmsy=
DGEMMSY_INTERFACE = dgemmsy/interface_dgemmsy_fake.f90
libs_dgemmsy=
endif

libs_opt=$(libs_cuda)\
 $(libs_ocl)\
 $(libs_dgemmsy)

if HAVE_LIBCONFIG
sub_dirs_libconfig=convolutions-c
else
sub_dirs_libconfig=
endif

if HAVE_ETSF_IO
etsf_files = wavelib/plotting-etsf.f90 \
	wavelib/i-o-etsf.f90
else
etsf_files = wavelib/etsf_fake.f90
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
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda) \
	$(libs_ocl)
# Hack to avoid global.f90 inclusion into make dist.
am_global_OBJECTS = global.$(OBJEXT)
global_LINK = $(FCLD) $(AM_CFLAGS) $(CFLAGS) $(global_LDFLAGS) \
	$(LDFLAGS) -o $@


SUBDIRS = modules convolutions $(sub_dirs_ps) profiling $(sub_dirs_cuda) $(sub_dirs_ocl) $(sub_dirs_dgemmsy) art $(sub_dirs_libconfig)

EXTRA_DIST = \
	wavelib/intots.inc \
	wavelib/recs16.inc \
	wavelib/sym_16.inc \
	NEB_driver.sh 

AM_FCFLAGS = -I. -Imodules @MPI_INCLUDE@ -IPSolver -I$(srcdir)/PSolver \
	@LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@ @LIBETSFIO_INCLUDE@

bin_PROGRAMS = $(binaries_sources)
bin_SCRIPTS = $(scripts_sources)

lib_LIBRARIES = $(lib_bigdft)
module_DATA = $(mod_bigdft)

if USE_MPI
mpi_source =
LMPI = @MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
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
	$(GPU_INTERFACE) $(OCL_INTERFACE)

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
	$(GPU_INTERFACE) $(OCL_INTERFACE)

endif

convolution_objects = $(convolution_source:.f90=.o)

#modules for types, interfaces and base variables
modules_source = modules/base.f90 \
	modules/defs.f90 \
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
	pdos.f90 \
	vdwcorrection.f90 \
	lanczos_interface.f90 \
	lanczos_base.f90 \
	lanczos.f90 \
	scfloop_API.f90 \
	esatto.f90 \
	tddft.f90

#objects which can be compiled only with module_base
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
	wfn_opt/exctX.f90 \
	wfn_opt/precond_per_optim.f90 \
	wfn_opt/precond_hyb.f90 \
	wfn_opt/precond_slab.f90 \
	wfn_opt/kernel_per_optim.f90 \
	wfn_opt/kernel_slab_simple.f90 \
	wfn_opt/coupling_matrix.f90 \
	lib/fft/fft2d.f90 \
	$(etsf_files) \
	wavelib/transwaves.f90 \
	wavelib/scalprod.f90 \
	wavelib/plotting.f90 \
	wavelib/i-o.f90 \
	wavelib/un-compress.f90 \
	wavelib/razero.f90 \
	wavelib/gauss_to_daub.f90 \
	wavelib/daubisf.f90 \
        linear/filterModule.f90 \
        linear/getLocBasis.f90 \
        linear/applyHamiltonian.f90 \
        linear/precondLinear.f90 \
        linear/locReg.f90 \
        linear/psitolocReg.f90 \
        linear/locReg_projectors.f90 \
        linear/transposition.f90 \
        linear/linearScaling.f90 \
        linear/energyAndForces.f90 \
        linear/convolutions.f90 \
        linear/initAndUtils.f90 \
        linear/sumrhoLinear.f90 \
	$(DGEMMSY_INTERFACE)  

libbigdft_a_SOURCES = \
	$(mpi_source) \
	BigDFT_API.f90 \
	$(low_level) \
	$(high_level)

libbigdft_a_LIBADD = profiling/time.o \
	profiling/memoryestimator.o \
	$(convolution_objects) \
	$(modules_objects) 

bigdft_SOURCES = BigDFT.f90 
bigdft_LDFLAGS = 
bigdft_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt)

splsad_SOURCES = splinedsaddle.f90
splsad_LDFLAGS =
splsad_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt)

rism_SOURCES = rism.f90 rismlowlevel.f90
rism_LDFLAGS = 
rism_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda) \
	$(libs_ocl)

sandbox_SOURCES = sandbox.f90
sandbox_LDFLAGS = 
sandbox_LDADD = libbigdft.a \
        PSolver/libpoissonsolver.a \
        $(LMPI) \
        @MPI_LIBS@ \
        $(libs_opt)

 
oneatom_SOURCES = oneatom.f90
oneatom_LDFLAGS = 
oneatom_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt)


MDanalysis_SOURCES = distances.f90 
MDanalysis_LDFLAGS = 
MDanalysis_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt)  

abscalc_SOURCES = abscalc.f90 
abscalc_LDFLAGS = 
abscalc_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt) 

test_forces_SOURCES = test_forces.f90
test_forces_LDFLAGS = 
test_forces_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt) 

NEB_SOURCES = NEB.f90 \
	$(mpi_source)
NEB_LDFLAGS = 
NEB_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(libs_opt) 

memguess_SOURCES = memguess.f90 \
	$(mpi_source)
memguess_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt) 

frequencies_SOURCES = frequencies.f90  \
	$(mpi_source)
frequencies_LDFLAGS = 
frequencies_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt) 

bart_SOURCES = \
	art/bigdft_forces.f90 \
	art/art_diis.f90 \
	art/initialize_potential.f90 \
	art/art_lanczos.f90 \
	art/find_saddle.f90 \
	art/min_converge.f90 \
	art/saddle_converge.f90 \
	art/initialize.f90 \
	art/read_parameters.f90 \
	art.f90
bart_FCFLAGS = -I. -Iart -Imodules @LIBABINIT_INCLUDE@
bart_LDFLAGS = 
bart_LDADD = art/libart.a \
	libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt) 

CLEANFILES = \
	bigdft_api.@MODULE_EXT@ BIGDFT_API.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ MODULE_INTERFACES.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ MODULE_TYPES.@MODULE_EXT@ \
	module_base.@MODULE_EXT@ MODULE_BASE.@MODULE_EXT@ \
	vdwcorrection.@MODULE_EXT@ VDWCORRECTION.@MODULE_EXT@ \
	lanczos_base.@MODULE_EXT@ LANCZOS_BASE.@MODULE_EXT@ \
	esatto.@MODULE_EXT@ ESATTO.@MODULE_EXT@ \
	lanczos_interface.@MODULE_EXT@ LANCZOS_INTERFACE.@MODULE_EXT@ \
	formats.@MODULE_EXT@ FORMATS.@MODULE_EXT@ \
	miscellany.@MODULE_EXT@ MISCELLANY.@MODULE_EXT@ \
	neb_routines.@MODULE_EXT@ NEB_ROUTINES.@MODULE_EXT@ \
	numeric.@MODULE_EXT@ NUMERIC.@MODULE_EXT@ \
	geopt_minimization.@MODULE_EXT@ MINIMIZATION.@MODULE_EXT@ \
	neb_variables.@MODULE_EXT@ NEB_VARIABLES.@MODULE_EXT@ \
	scfloop_api.@MODULE_EXT@ SCFLOOP_API.@MODULE_EXT@ \
	minimization.@MODULE_EXT@ MINIMIZATION.@MODULE_EXT@ \
	bigdft_forces.@MODULE_EXT@ BIGDFT_FORCES.@MODULE_EXT@ \
	diis_defs.@MODULE_EXT@ DIIS_DEFS.@MODULE_EXT@ \
	lanczos_defs.@MODULE_EXT@ LANCZOS_DEFS.@MODULE_EXT@ \
	minimization_routines.@MODULE_EXT@ MINIMIZATION_ROUTINES.@MODULE_EXT@ \
	minpar.@MODULE_EXT@ MINPAR.@MODULE_EXT@ \
	saddles.@MODULE_EXT@ SADDLES.@MODULE_EXT@


check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

#dependencies
$(libbigdft_a_SOURCES:.f90=.o) BigDFT.o test_forces.o  memguess.o global.o : $(mpi_include)  $(modules_objects)

cluster.o geometry.o xabsorber.o: vdwcorrection.o esatto.o scfloop_API.o

lanczos.o : lanczos_base.o lanczos_interface.o

BigDFT_API.o : scfloop_API.o

distances.o global.o rism.o sandbox.o : BigDFT_API.o

xabsorber.o init/xabsorber.o : init/xabsorber.f90 esatto.o

#init.o : base.o interfaces.o

#BigDFT_API.o: base.o types.o interfaces.o

#$(high_level:.f90=.o) : base.o types.o

test_forces.o: $(mod_bigdft)

$(mod_bigdft) : BigDFT_API.o

$(modules_objects) : $(mpi_include)

modules/defs.o: modules/base.o

#BigDFT.o: \
#	base.o \
#	interfaces.o

#memguess.o: base.o types.o

# Dependencies for ART.
bart-art.o: art.f90 art/defs.o art/random.o bart-art_lanczos.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-art.o `test -f 'art.f90' || echo '$(srcdir)/'`art.f90

bart-bigdft_forces.o: art/bigdft_forces.f90 BigDFT_API.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-bigdft_forces.o `test -f 'art/bigdft_forces.f90' || echo '$(srcdir)/'`art/bigdft_forces.f90

bart-find_saddle.o: art/find_saddle.f90 art/defs.o art/random.o bart-art_lanczos.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-find_saddle.o `test -f 'art/find_saddle.f90' || echo '$(srcdir)/'`art/find_saddle.f90

bart-initialize_potential.o: art/initialize_potential.f90 art/defs.o bart-bigdft_forces.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-initialize_potential.o `test -f 'art/initialize_potential.f90' || echo '$(srcdir)/'`art/initialize_potential.f90

bart-art_lanczos.o: art/art_lanczos.f90 art/defs.o art/random.o bart-bigdft_forces.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-art_lanczos.o `test -f 'art/art_lanczos.f90' || echo '$(srcdir)/'`art/art_lanczos.f90

bart-min_converge.o: art/min_converge.f90 art/defs.o bart-bigdft_forces.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-min_converge.o `test -f 'art/min_converge.f90' || echo '$(srcdir)/'`art/min_converge.f90

bart-saddle_converge.o: art/saddle_converge.f90 art/defs.o bart-bigdft_forces.o bart-find_saddle.o art/random.o bart-art_lanczos.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-saddle_converge.o `test -f 'art/saddle_converge.f90' || echo '$(srcdir)/'`art/saddle_converge.f90

bart-art_diis.o: art/art_diis.f90 art/defs.o bart-bigdft_forces.o art/random.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-art_diis.o `test -f 'art/art_diis.f90' || echo '$(srcdir)/'`art/art_diis.f90

bart-initialize.o: art/initialize.f90 art/defs.o bart-art_lanczos.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-initialize.o `test -f 'art/initialize.f90' || echo '$(srcdir)/'`art/initialize.f90

bart-read_parameters.o: art/read_parameters.f90 art/defs.o bart-art_lanczos.o bart-find_saddle.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-read_parameters.o `test -f 'art/read_parameters.f90' || echo '$(srcdir)/'`art/read_parameters.f90

bart-utils.o: art/utils.f90 art/defs.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-utils.o `test -f 'art/utils.f90' || echo '$(srcdir)/'`art/utils.f90

bart-write_refconfig.o: art/write_refconfig.f90 art/defs.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-write_refconfig.o `test -f 'art/write_refconfig.f90' || echo '$(srcdir)/'`art/write_refconfig.f90
