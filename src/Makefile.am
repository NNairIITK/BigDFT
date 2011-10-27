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
binaries_sources = bigdft memguess splsad NEB frequencies MDanalysis rism oneatom $(minima_hopping_sources) bart abscalc
local_binaries = test_forces
scripts_sources = NEB_driver.sh
else
binaries_sources = 
local_binaries =
scripts_sources = 
endif

if HAVE_ETSF_IO
binaries_sources += BigDFT2Wannier
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
libs_ocl=OpenCL/libOCL.a -lstdc++
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

if HAVE_FC_GET_COMMAND_ARGUMENT
getarg_files =
else
getarg_files = modules/get_command_argument_fake.f90
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
	$(libs_opt)
# Hack to avoid global.f90 inclusion into make dist.
am_global_OBJECTS = global.$(OBJEXT)
global_LINK = $(FCLD) $(FCFLAGS) $(AM_CFLAGS) $(CFLAGS) $(global_LDFLAGS) \
	$(LDFLAGS) -o $@


SUBDIRS = modules convolutions $(sub_dirs_ps) profiling $(sub_dirs_cuda) $(sub_dirs_ocl) $(sub_dirs_dgemmsy) art $(sub_dirs_libconfig) tools tools/bader

EXTRA_DIST = \
	wavelib/intots.inc \
	wavelib/recs16.inc \
	wavelib/sym_16.inc \
	init/pspconf.in.f90 \
	NEB_driver.sh 

AM_FCFLAGS = -I. -Imodules @MPI_INCLUDE@ -IPSolver -I$(srcdir)/PSolver -I$(srcdir)/wavelib \
	@LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@ \
	@LIBETSFIO_INCLUDE@
AM_CPPFLAGS = @LIB_ARCHIVE_CFLAGS@

bin_PROGRAMS = $(binaries_sources)
noinst_PROGRAMS = $(local_binaries)
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
	modules/input.f90 \
	modules/xc.f90 \
	modules/op2p_module.f90 \
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
	ConstrainedDavidson.f90 \
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
	init/pspconf.f90 \
	init/locreg.f90 \
	init/sysprop.f90 \
	init/gaussians.f90 \
	init/wfn_init.f90 \
	wfn_opt/precond.f90 \
	wfn_opt/orthogonality.f90 \
	wfn_opt/diis.f90 \
	wfn_opt/applyh.f90 \
	wfn_opt/exctX.f90 \
	wfn_opt/sic.f90 \
	wfn_opt/precond_per_optim.f90 \
	wfn_opt/precond_hyb.f90 \
	wfn_opt/precond_slab.f90 \
	wfn_opt/kernel_per_optim.f90 \
	wfn_opt/kernel_slab_simple.f90 \
	wfn_opt/coupling_matrix.f90 \
	lib/fft/fft2d.f90 \
	lib/lbfgs.f90 \
	$(etsf_files) \
	wavelib/transwaves.f90 \
	wavelib/scalprod.f90 \
	wavelib/plotting.f90 \
	wavelib/i-o.f90 \
	wavelib/un-compress.f90 \
	wavelib/razero.f90 \
	wavelib/gauss_to_daub.f90 \
	wavelib/daubisf.f90 \
	sdcg.f90 \
	bfgs.f90 \
	$(DGEMMSY_INTERFACE)

# C objects
c_level = posfiles.c \
	utils.c

libbigdft_a_SOURCES = \
	$(mpi_source) \
	BigDFT_API.f90 \
	$(c_level) \
	$(low_level) \
	$(high_level) \
	$(getarg_files)

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
	$(libs_opt)

oneatom_SOURCES = oneatom.f90
oneatom_LDFLAGS =  
oneatom_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt)


BigDFT2Wannier_SOURCES = BigDFT2Wannier.f90                                                                                                                                                                                  
BigDFT2Wannier_LDFLAGS = 
BigDFT2Wannier_LDADD = libbigdft.a \
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
	MPIfake.f90
memguess_LDFLAGS = 
memguess_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
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
        art/art_step.f90 \
        art/calcforce.f90\
        art/calcfo_sw.f90\
        art/neighbour.f90\
        art/initialize_potential.f90 \
        art/art_lanczos.f90 \
        art/find_saddle.f90 \
        art/min_converge.f90 \
        art/saddle_converge.f90 \
        art/initialize.f90 \
        art/bayes.f90 \
        art/read_parameters.f90 \
        art.f90
bart_FCFLAGS = -I. -Iart -Imodules @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@
bart_LDFLAGS =  
bart_LDADD = art/libart.a \
	libbigdft.a \
	PSolver/libpoissonsolver.a \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_opt) 

CLEANFILES = *.@MODULE_EXT@


check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

#dependencies
$(low_level:.f90=.o) $(high_level:.f90=.o) BigDFT_API.o BigDFT.o test_forces.o  memguess.o global.o abscalc.o : $(mpi_include)  $(modules_objects)

cluster.o geometry.o xabsorber.o: vdwcorrection.o esatto.o scfloop_API.o

lanczos.o : lanczos_base.o lanczos_interface.o

BigDFT_API.o : scfloop_API.o

distances.o global.o rism.o oneatom.o BigDFT2Wannier.o : BigDFT_API.o

xabsorber.o init/xabsorber.o : init/xabsorber.f90 esatto.o

#init.o : base.o interfaces.o

#BigDFT_API.o: base.o types.o interfaces.o

$(high_level:.f90=.o) : $(modules_objects)

test_forces.o: $(mod_bigdft)

$(mod_bigdft) : BigDFT_API.o

$(modules_objects) : $(mpi_include)

#BigDFT.o: \
#	base.o \
#	interfaces.o

#memguess.o: base.o types.o

# Dependencies for ART.
bart-art.o: art.f90 art/defs.o art/random.o bart-art_lanczos.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-art.o `test -f 'art.f90' || echo '$(srcdir)/'`art.f90

bart-bigdft_forces.o: art/bigdft_forces.f90 BigDFT_API.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-bigdft_forces.o `test -f 'art/bigdft_forces.f90' || echo '$(srcdir)/'`art/bigdft_forces.f90

bart-bayes.o: art/bayes.f90 art/defs.o art/random.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-bayes.o `test -f 'art/bayes.f90' || echo '$(srcdir)/'`art/bayes.f90

bart-calcforce.o: art/calcforce.f90 art/defs.o bart-bigdft_forces.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-calcforce.o `test -f 'art/calcforce.f90' || echo '$(srcdir)/'`art/calcforce.f90

bart-calcfo_sw.o: art/calcfo_sw.f90 art/defs.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-calcfo_sw.o `test -f 'art/calcfo_sw.f90' || echo '$(srcdir)/'`art/calcfo_sw.f90

bart-neighbour.o: art/neighbour.f90 bart-calcfo_sw.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-neighbour.o `test -f 'art/neighbour.f90' || echo '$(srcdir)/'`art/neighbour.f90

bart-find_saddle.o: art/find_saddle.f90 art/defs.o art/random.o bart-art_lanczos.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-find_saddle.o `test -f 'art/find_saddle.f90' || echo '$(srcdir)/'`art/find_saddle.f90

bart-initialize_potential.o: art/initialize_potential.f90 art/defs.o bart-bigdft_forces.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-initialize_potential.o `test -f 'art/initialize_potential.f90' || echo '$(srcdir)/'`art/initialize_potential.f90

bart-art_lanczos.o: art/art_lanczos.f90 art/defs.o art/random.o bart-bigdft_forces.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-art_lanczos.o `test -f 'art/art_lanczos.f90' || echo '$(srcdir)/'`art/art_lanczos.f90

bart-min_converge.o: art/min_converge.f90 art/defs.o bart-bigdft_forces.o bart-art_lanczos.o bart-bayes.o bart-find_saddle.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-min_converge.o `test -f 'art/min_converge.f90' || echo '$(srcdir)/'`art/min_converge.f90

bart-saddle_converge.o: art/saddle_converge.f90 art/defs.o bart-bigdft_forces.o bart-find_saddle.o bart-art_lanczos.o  bart-art_step.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-saddle_converge.o `test -f 'art/saddle_converge.f90' || echo '$(srcdir)/'`art/saddle_converge.f90

bart-art_step.o: art/art_step.f90 art/defs.o bart-bigdft_forces.o bart-art_lanczos.o  bart-find_saddle.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-art_step.o `test -f 'art/art_step.f90' || echo '$(srcdir)/'`art/art_step.f90

bart-initialize.o: art/initialize.f90 art/defs.o bart-bigdft_forces.o bart-art_lanczos.o bart-find_saddle.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-initialize.o `test -f 'art/initialize.f90' || echo '$(srcdir)/'`art/initialize.f90

bart-read_parameters.o: art/read_parameters.f90 art/defs.o bart-art_lanczos.o bart-find_saddle.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-read_parameters.o `test -f 'art/read_parameters.f90' || echo '$(srcdir)/'`art/read_parameters.f90

bart-utils.o: art/utils.f90 art/defs.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-utils.o `test -f 'art/utils.f90' || echo '$(srcdir)/'`art/utils.f90

bart-write_refconfig.o: art/write_refconfig.f90 art/defs.o
	$(FC) $(bart_FCFLAGS) $(FCFLAGS) -c -o bart-write_refconfig.o `test -f 'art/write_refconfig.f90' || echo '$(srcdir)/'`art/write_refconfig.f90
