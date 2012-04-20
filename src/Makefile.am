## Process this file with automake to produce Makefile.in

all: all-recursive
	chmod ug+x bigdft-tool

# Main targets.
# ------------
bin_PROGRAMS    = $(binaries_sources)
noinst_PROGRAMS = $(local_binaries)
bin_SCRIPTS     = $(scripts_sources)
lib_LIBRARIES   = $(lib_bigdft)
module_DATA     = $(mod_bigdft)

# Sub-directories.
# ---------------
SUBDIRS = modules \
	PSolver \
	$(sub_dirs_libconfig) \
	art \
	$(sub_dirs_cuda) \
	$(sub_dirs_ocl) \
	$(sub_dirs_bin)

EXTRA_DIST = \
	wavelib/intots.inc \
	wavelib/recs16.inc \
	wavelib/sym_16.inc \
	init/pspconf.in.f90 \
	convolutions/v.inc \
	convolutions/v_long.inc \
	convolutions/v_17.inc \
	PSolver/lazy_100.inc \
	PSolver/lazy_14.inc \
	PSolver/lazy_16.inc \
	PSolver/lazy_20.inc \
	PSolver/lazy_24.inc \
	PSolver/lazy_30.inc \
	PSolver/lazy_40.inc \
	PSolver/lazy_50.inc \
	PSolver/lazy_60.inc \
	PSolver/lazy_8.inc \
	PSolver/createKernel.f90 \
	PSolver/PSolver_Main.f90 \
	PSolver/perfdata.inc \
	NEB_driver.sh 

# Build targets, binaries and libraries.
# -------------
if BUILD_LIB_BIGDFT
lib_bigdft = libbigdft.a libdft_common.a
PS_SOURCES = PSolver/Build_Kernel.f90 \
	PSolver/xcenergy.f90 \
	PSolver/3Dgradient.f90 \
	PSolver/scaling_function.f90 \
	PSolver/PSolver_Base_new.f90 \
	PSolver/wofz.f90 \
	lib/fft/fft3d.f90
PS_MODULES = PSolver/Poisson_Solver.f90
else
lib_bigdft =
PS_SOURCES =
PS_MODULES =
endif
if BUILD_BINARIES
binaries_sources = bigdft memguess splsad NEB frequencies MDanalysis $(minima_hopping_sources) bart abscalc BigDFT2Wannier WaCo
local_binaries = test_forces
scripts_sources = NEB_driver.sh bigdft-tool
sub_dirs_bin = tools tools/bader
else
binaries_sources = 
local_binaries =
scripts_sources = 
sub_dirs_bin =
endif

# Additional sources.
# ------------------
#if USE_BLACS
#BLACS_MODULES = 
#else
BLACS_MODULES = modules/blacs_fake.f90
#endif

EXTRA_DIST += CUDA/precond.cu \
	CUDA/density.cu \
	CUDA/anasyn.cu \
	CUDA/locpot.cu CUDA/locpot.h \
	CUDA/kinetic.cu CUDA/kinetic.h \
	CUDA/compress.cu
if USE_CUDA_GPU
sub_dirs_cuda = CUDA
CUDA_SOURCES = CUDA/binding_s_gpu.c \
	CUDA/commonDef.cpp CUDA/commonDef.h \
	CUDA/fortran.c \
	CUDA/cudafct.cu \
	CUDA/locham.cu \
	CUDA/kernels_anasyn.hcu \
	CUDA/kernels_locpot.hcu \
	CUDA/kernels_kinetic.hcu \
	CUDA/kernels_compress.hcu \
	CUDA/GPUparameters.h \
	CUDA/structDef_anasyn.h \
	CUDA/structDef_locpot.h \
	CUDA/structDef_kinetic.h \
	CUDA/check_cuda.h \
	CUDA/read_conf_exception.h \
	CUDA/cpp_utils.h
CUDA_CONV = convolutions/convolut_common_interface_cuda.f90
CUDA_MODULES =
else
sub_dirs_cuda =
CUDA_SOURCES =
CUDA_CONV = convolutions/interface_cuda_fake.f90
CUDA_MODULES = modules/cublas_fake.f90
endif

if USE_OCL
sub_dirs_ocl = OpenCL
OCL_SOURCES = \
	OpenCL/Stream.h OpenCL/Tool.h OpenCL/bench_lib.h \
	OpenCL/OpenCL_wrappers.c         OpenCL/OpenCL_wrappers.h \
	OpenCL/MagicFilter.c             OpenCL/MagicFilter.h \
	OpenCL/MagicFilter_Generator.cpp OpenCL/MagicFilter_Generator.h \
	OpenCL/Kinetic.c                 OpenCL/Kinetic.h \
	OpenCL/Kinetic_Generator.cpp     OpenCL/Kinetic_Generator.h \
	OpenCL/Kinetic_k.h \
	OpenCL/Kinetic_k_Generator.cpp   OpenCL/Kinetic_k_Generator.h \
	OpenCL/Wavelet.c                 OpenCL/Wavelet.h \
	OpenCL/Wavelet_Generator.cpp     OpenCL/Wavelet_Generator.h \
	OpenCL/Uncompress.c              OpenCL/Uncompress.h \
	OpenCL/Initialize.c              OpenCL/Initialize.h \
	OpenCL/Profiling.c \
	OpenCL/Hamiltonian.c \
	OpenCL/Reduction.c               OpenCL/Reduction.h \
	OpenCL/Reduction_Generator.cpp   OpenCL/Reduction_Generator.h \
	OpenCL/Preconditioner.c \
	OpenCL/Density.c \
	OpenCL/dsfft.c \
	OpenCL/fft_generator_noshared.cpp \
	OpenCL/fft_noshared.c \
	OpenCL/fft_generator_helper.cpp \
	OpenCL/dsfft_generator.cpp       OpenCL/fft_generator.h \
	OpenCL/Benchmark.c \
	OpenCL/Benchmark_Generator.cpp   OpenCL/Benchmark_Generator.h
OCL_CONV = convolutions/convolut_common_interface_ocl.f90
else
sub_dirs_ocl =
OCL_SOURCES =
OCL_CONV = convolutions/interface_ocl_fake.f90
endif

if USE_DGEMMSY
DGEMMSY_SOURCES = dgemmsy/interface_dgemmsy.f90 \
		  dgemmsy/dgemmsy.c       dgemmsy/dgemmsy.h \
		  dgemmsy/gemm_block_c.c  dgemmsy/gemm_block_c.h \
		  dgemmsy/patterns.c      dgemmsy/patterns.h \
		  dgemmsy/dgemmsy_utils.c dgemmsy/utils.h \
		  dgemmsy/visitors.c      dgemmsy/visitors.h
else
DGEMMSY_SOURCES = dgemmsy/interface_dgemmsy_fake.f90
endif

if USE_LIBXC
LIBXC_MODULES =
else
LIBXC_MODULES = modules/libxc_fake.f90
endif

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

if BUILD_BINDINGS
bindings_files = lib/bindings_api.h \
	lib/bindings.c lib/bindings.h \
	lib/bindings_atoms.c \
	lib/bindings_regions.c \
	lib/bindings_wavefunctions.c \
	lib/bindings_localfields.c \
	lib/bindingsf.f90
else
bindings_files = lib/bindingsf_fake.f90
endif

if HAVE_FC_FLUSH
flush_files = lib/flush.f90
else
flush_files = lib/flush_fake.f90
endif

if HAVE_FC_GET_COMMAND_ARGUMENT
getarg_files =
else
getarg_files = modules/get_command_argument_fake.f90
endif

# Specific treatment for global (not distributed).
if USE_MINIMA_HOPPING
minima_hopping_sources = global
else
minima_hopping_sources =
endif
global_SOURCES = # Should be global.f90, but see below.
global_LDADD = @LIBS_DEPENDENCIES@
# Hack to avoid global.f90 inclusion into make dist.
am_global_OBJECTS = global.$(OBJEXT)
global_LINK = $(FCLD) $(FCFLAGS) $(AM_CFLAGS) $(CFLAGS) $(global_LDFLAGS) \
	$(LDFLAGS) -o $@

AM_FCFLAGS = -I. -Imodules @MPI_INCLUDE@ -IPSolver \
	-I$(srcdir)/PSolver -I$(srcdir)/wavelib -I$(srcdir)/convolutions \
	@LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@ \
	@LIBETSFIO_INCLUDE@ @LIBSGPU_INCLUDE@
AM_CPPFLAGS = @LIB_ARCHIVE_CFLAGS@ @OCL_INCLUDE@ -I$(srcdir)/OpenCL @LIBSGPU_INCLUDE@ -I$(srcdir)/CUDA -I$(srcdir)/modules -Imodules
#AM_CFLAGS = -Wall -Werror -Wno-strict-aliasing @DGEMMSY_CPPFLAGS@ @CUDA_INCLUDE@
AM_CFLAGS = @DGEMMSY_CPPFLAGS@ @CUDA_INCLUDE@ @GLIB_CFLAGS@

if USE_MPI
mpi_source =
else
mpi_source = MPIfake.f90
endif

if USE_OPTI_CONVOLUT
CONV_SOURCES = convolutions/combined_shrink_optim.f90 \
	convolutions/combined_grow_optim.f90 \
	convolutions/growshrink_hyb_optim.f90 \
	convolutions/convolut_ib_optim.f90 \
	convolutions/conv_per_optim.f90 \
	convolutions/convolut_optim_per.f90 \
	convolutions/convolut_optim_slab.f90
else
CONV_SOURCES = convolutions/combined_shrink_simple.f90 \
	convolutions/combined_grow_simple.f90 \
	convolutions/growshrink_hyb_simple.f90 \
	convolutions/convolut_ib_simple.f90 \
	convolutions/conv_per_simple.f90 \
	convolutions/convolut_simple_per.f90 \
	convolutions/convolut_simple_slab.f90
endif
CONV_SOURCES += \
	convolutions/growshrink.f90 \
	convolutions/growshrink_hyb_common.f90 \
	convolutions/conv_per_common.f90 \
	convolutions/convolut_common_slab.f90 \
	convolutions/convolut_common_per.f90 \
	convolutions/convolut_new_per.f90 \
	convolutions/convolut_simple_per_k.f90 \
	convolutions/convolut_simple_slab_k.f90 \
	$(CUDA_CONV) $(OCL_CONV)

#modules for types, interfaces and base variables
modules_source = modules/base.f90 \
	modules/defs.f90 \
	modules/types.f90 \
	modules/interfaces.f90 \
	modules/private_api.f90 \
	modules/input.f90 \
	modules/xc.f90 \
	modules/op2p_module.f90 \
	modules/BigDFT_API.f90 \
	modules/vdwcorrection.f90 \
	modules/yaml_output.f90 \
	$(CUDA_MODULES) \
	$(LIBXC_MODULES) \
	$(PS_MODULES) \
	$(BLACS_MODULES)
modules_objects = $(modules_source:.f90=.o)

#objects which depends from module_types
high_level = cluster.f90 \
	init.f90 \
	input_variables.f90 \
	sumrho.f90 \
	hpsiortho.f90 \
	forces.f90 \
	tail.f90 \
	restart.f90 \
	geometry.f90 \
	davidson.f90 \
	ConstrainedDavidson.f90 \
	pdos.f90 \
	sdcg.f90 \
	bfgs.f90 \
	lanczos_interface.f90 \
	lanczos_base.f90 \
	lanczos.f90 \
	scfloop_API.f90 \
	esatto.f90 \
	tddft.f90 \
	output.f90

#objects which can be compiled only with module_base
low_level = \
	$(PS_SOURCES) \
	init/atoms.f90 \
	init/projectors.f90  \
	init/kswfn.f90 \
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
	init/wavefunctions.f90 \
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
	$(bindings_files) \
	$(flush_files) \
	$(etsf_files) \
	wavelib/transwaves.f90 \
	wavelib/scalprod.f90 \
	wavelib/plotting.f90 \
	wavelib/i-o.f90 \
	wavelib/un-compress.f90 \
	wavelib/razero.f90 \
	wavelib/gauss_to_daub.f90 \
	wavelib/daubisf.f90 \
	profiling/time.f90 \
	profiling/memoryestimator.f90 \
	linear/cleanup.f90 \
	linear/nullify.f90 \
	linear/filterModule.f90 \
	linear/communicatetypes.f90 \
	linear/misc.f90 \
	linear/getLocBasis.f90 \
	linear/convolutions_aux.f90 \
	linear/convolutions.f90 \
	linear/applyHamiltonian.f90 \
	linear/precondLinear.f90 \
	linear/psitolocReg.f90 \
	linear/locReg_orbitals.f90 \
	linear/transposition.f90 \
	linear/linearScaling.f90 \
	linear/energyAndForces.f90 \
	linear/initAndUtils.f90 \
	linear/sumrhoLinear.f90 \
	linear/inputguessLinear.f90 \
	linear/psitolocReg.f90 \
	linear/orthonormality.f90 \
	linear/orbitalOptimization.f90 \
	linear/mixdiis.f90 \
	linear/copyTypes.f90 \
	linear/de-allocate.f90 \
	linear/derivativeBasis.f90 \
	linear/parallelLinAlg.f90 \
	linear/locReg.f90 \
	linear/unitaryTrans.f90 \
	linear/communicatePotential.f90 \
	linear/hpsiortho_linear.f90 \
	linear/init_collective_comms.f90 \
	linear/init_collective_comms_vectors.f90 \
	$(CONV_SOURCES) \
	$(DGEMMSY_SOURCES) \
	$(OCL_SOURCES) \
	$(getarg_files)

# C objects
c_level = posfiles.c \
	$(CUDA_SOURCES)

#	convolutions/convolut_ib_sse.c \ #to be verified, not working yet

libdft_common_a_SOURCES = \
	utils.c

libbigdft_a_SOURCES = \
	$(mpi_source) \
	$(c_level) \
	$(low_level) \
	$(high_level)
libbigdft_a_LIBADD = $(modules_objects) 
$(bin_PROGRAMS) $(noinst_PROGRAMS): libbigdft.a libdft_common.a

bigdft_SOURCES = BigDFT.f90
bigdft_LDADD = @LIBS_DEPENDENCIES@

splsad_SOURCES = splinedsaddle.f90
splsad_LDADD = @LIBS_DEPENDENCIES@

#rism_SOURCES = rism.f90 rismlowlevel.f90
#rism_LDADD = @LIBS_DEPENDENCIES@

#sandbox_SOURCES = sandbox.f90
#sandbox_LDADD = @LIBS_DEPENDENCIES@

#oneatom_SOURCES = oneatom.f90
#oneatom_LDADD = @LIBS_DEPENDENCIES@

BigDFT2Wannier_SOURCES = BigDFT2Wannier.f90
BigDFT2Wannier_LDADD = @LIBS_DEPENDENCIES@

WaCo_SOURCES = WaCo.f90                                                                                                                                                                    
WaCo_LDADD = @LIBS_DEPENDENCIES@

MDanalysis_SOURCES = distances.f90
MDanalysis_LDADD = @LIBS_DEPENDENCIES@ 

abscalc_SOURCES = abscalc.f90
abscalc_LDADD = @LIBS_DEPENDENCIES@

test_forces_SOURCES = test_forces.f90
test_forces_LDADD = @LIBS_DEPENDENCIES@

NEB_SOURCES = NEB.f90
NEB_LDADD = @LIBS_DEPENDENCIES@

memguess_SOURCES = memguess.f90 MPIfake.f90
memguess_LDADD = @LIBS_DEPENDENCIES@

frequencies_SOURCES = frequencies.f90
frequencies_LDADD = @LIBS_DEPENDENCIES@

bart_SOURCES = art.f90
bart_FCFLAGS = -I. -Iart -Imodules @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@
bart_LDADD = art/libart.a @LIBS_DEPENDENCIES@

CLEANFILES = *.@MODULE_EXT@

NVCC = @NVCC@
NVCC_FLAGS = -arch sm_13 @NVCC_FLAGS@ 
SUFFIXES = .cu
.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@

check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

#dependencies
$(low_level:.f90=.o) $(high_level:.f90=.o) BigDFT.o test_forces.o  memguess.o global.o abscalc.o distances.o global.o rism.o oneatom.o BigDFT2Wannier.o WaCo.o: $(mpi_include)  $(modules_objects)

cluster.o geometry.o xabsorber.o: esatto.o scfloop_API.o

lanczos.o : lanczos_base.o lanczos_interface.o

xabsorber.o init/xabsorber.o : init/xabsorber.f90 esatto.o

$(modules_objects) : $(mpi_include)
