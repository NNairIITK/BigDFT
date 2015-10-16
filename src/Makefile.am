## Process this file with automake to produce Makefile.in

# Main targets.
# ------------
dynlibdir           = $(libdir)
bin_PROGRAMS        = $(binaries_sources)
noinst_PROGRAMS     = $(local_binaries)
bin_SCRIPTS         = $(scripts_sources)
lib_LIBRARIES       = $(lib_bigdft)
module_DATA         = $(mod_bigdft)
dynlib_DATA         = $(bigdft_dynamic_library)

# Sub-directories.
# ---------------
SUBDIRS = \
	orbitals \
	modules \
	xcPSolver \
	$(sub_dirs_libconfig) \
	CUDA \
	OpenCL \
	bindings \
	geopt \
	art \
	mhgps_files \
	global_files \
	$(sub_dirs_bin)

EXTRA_DIST = \
	wavelib/intots.inc \
	wavelib/recs16.inc \
	wavelib/sym_16.inc \
	init/pspconf.in.f90 \
	NEB_driver.sh \
	input_variables_definition.yaml

# Build targets, binaries and libraries.
# -------------
noinst_LIBRARIES = libmain.a libabscalc.a
# put all subs
subs_static_library = \
	modules/libmodules.a \
	CUDA/libCUDA.a \
	orbitals/liborbitals.a \
	geopt/libgeopt.a \
	xcPSolver/libxcPSolver.a \
	OpenCL/libOpenCL.a \
	bindings/libbindings.a
out_static_library =  \
	$(top_builddir)/PSolver/src/libPSolver-1.a \
	$(top_builddir)/wrappers/libwrappers.a \
	$(top_builddir)/flib/src/libflib-1.a
if BUILD_DYNAMIC_LIBS
bigdft_dynamic_library = libbigdft-1.so.@BIGDFT_MINOR_VERSION@.0.@BIGDFT_MICRO_VERSION@
bigdft_soname_library = libbigdft-1.so.@BIGDFT_MINOR_VERSION@
bigdft_library = $(bigdft_soname_library)
AM_LDFLAGS = -Wl,-rpath=$(DESTDIR)$(dynlibdir)
else
bigdft_dynamic_library =
bigdft_soname_library =
bigdft_library = libbigdft-1.a $(out_static_library) @LIBS_SHORT_DEPS@ $(amber_ld) $(cp2k_ld)
endif
if BUILD_LIB_BIGDFT
lib_bigdft = libbigdft-1.a
else
lib_bigdft =
endif
if BUILD_BINARIES
binaries_sources = bigdft memguess splsad NEB frequencies global globaltool abscalc mhgps mhgpstool utilities
local_binaries = test_forces
scripts_sources = NEB_driver.sh bigdft-tool
sub_dirs_bin = tools tools/bader
EXTRA_PROGRAMS = MDanalysis bart BigDFT2Wannier WaCo NEB_images
else
binaries_sources =
local_binaries =
scripts_sources =
sub_dirs_bin =
check_PROGRAMS =
endif

# Additional sources.
# ------------------
if HAVE_LIBCONFIG
sub_dirs_libconfig=convolutions-c
else
sub_dirs_libconfig=
endif

if HAVE_ETSF_IO
etsf_files = wavelib/plotting-etsf.f90
else
etsf_files = wavelib/etsf_fake.f90
endif


if HAVE_AMBERTOOLS
amber_ld= @LIB_AMBERTOOLS@  -lsff -lnab -lpbsa -lcifparse -lrism -lfftw3 -larpack
ambertools_wrap= wrap.c
###for amber supp.:
wrap.o:     wrap.c
	$(NAB) -c  -g $<
else
ambertools_wrap= wrap_fake.f90
amber_ld =
endif

if HAVE_CP2K
cp2k_ld= @CP2K_LINKLINE@ @LINALG_LIBS@
##cp2k_interfacemod = modules/module_cp2k.f90
##cp2k_interfacemodCAP = MODULES/MODULE_CP2k.F90
else
##cp2k_interfacemod = modules/cp2k_fake.f90 modules/module_cp2k.f90
##cp2k_interfacemodCAP = MODULES/CP2K_FAKE.F90 MODULES/MODULE_CP2K_FAKE.F90
cp2k_ld =
endif

AM_FCFLAGS = -I. -I$(top_builddir)/includes -I$(srcdir)/wavelib \
	@LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@ \
	@LIBETSFIO_INCLUDE@ @LIBSGPU_INCLUDE@ \
	@MPI_INCLUDE@
AM_CPPFLAGS = @LIB_ARCHIVE_CFLAGS@ @LIBSGPU_INCLUDE@ -I$(srcdir)/modules -Imodules
#AM_CFLAGS = -Wall -Werror -Wno-strict-aliasing @DGEMMSY_CPPFLAGS@ @CUDA_INCLUDE@ -Wno-deprecated-declarations -Wno-cpp
AM_CFLAGS = @CUDA_INCLUDE@

#objects which depends from module_types
high_level = cluster.f90 \
	init.f90 \
	initialization.f90 \
	input_variables.f90 \
	sumrho.f90 \
	hpsiortho.f90 \
	forces.f90 \
	tail.f90 \
	restart.f90 \
	davidson.f90 \
	ConstrainedDavidson.f90 \
	pdos.f90 \
	tddft.f90 \
	output.f90 \
	external.f90 \
	images.f90

#objects which can be compiled only with module_base
low_level = \
	init/atoms.f90 \
	init/projectors.f90  \
	init/kswfn.f90 \
	init/ionicpot.f90 \
	init/gridmanipulation.f90 \
	init/gauprod.f90 \
	init/denspotd.f90 \
	init/inputguess.f90 \
	init/rcov.f90 \
	init/gautowav.f90 \
	init/pspconf.f90 \
	init/locreg_init.f90 \
	init/sysprop.f90 \
	init/wavefunctions.f90 \
	init/gaussians_old.f90 \
	init/wfn_init.f90 \
	wfn_opt/precond.f90 \
	wfn_opt/orthogonality.f90 \
	wfn_opt/diis.f90 \
	wfn_opt/applyh.f90 \
	wfn_opt/exctX.f90 \
	wfn_opt/exctX_op2p.f90 \
	wfn_opt/sic.f90 \
	wfn_opt/precond_per_optim.f90 \
	wfn_opt/precond_hyb.f90 \
	wfn_opt/precond_slab.f90 \
	wfn_opt/kernel_per_optim.f90 \
	wfn_opt/kernel_slab_simple.f90 \
	wfn_opt/coupling_matrix.f90 \
	wavelib/transwaves.f90 \
	wavelib/scalprod.f90 \
	wavelib/plotting.f90 \
	wavelib/i-o.f90 \
	wavelib/un-compress.f90 \
	wavelib/gauss_to_daub.f90 \
	wavelib/daubisf.f90 \
	$(etsf_files) \
	profiling/memoryestimator.f90 \
	linear/cleanup.f90 \
	linear/filterModule.f90 \
	linear/misc.f90 \
	linear/getLocBasis.f90 \
	linear/convolutions_aux.f90 \
	linear/convolutions.f90 \
	linear/precondLinear.f90 \
	linear/locReg_orbitals.f90 \
	linear/linearScaling.f90 \
	linear/initAndUtils.f90 \
	linear/sumrhoLinear.f90 \
	linear/inputguessLinear.f90 \
	linear/orbitalOptimization.f90 \
	linear/mixdiis.f90 \
	linear/copyTypes.f90 \
	linear/derivativeBasis.f90 \
	linear/hpsiortho_linear.f90 \
	linear/optimize_coeffs.f90 \
	linear/pulay.f90 \
	linear/transferintegrals.f90 \
	linear/cdft.f90 \
	linear/frags.f90\
	$(ambertools_wrap)
##	$(cp2k_interfacemod)

# C objects
c_level = posfiles.c \
	input_variables_description.c

libmain_a_SOURCES = \
	$(c_level) \
	$(low_level) \
	$(high_level)

libabscalc_a_SOURCES = \
	lanczos_base.f90 \
	abscalc_module.f90 \
	lanczos_interface.f90 \
	lanczos.f90 \
	abscalc_init.f90 \
	esatto.f90

libbigdft_1_a_SOURCES =

input_variables_description.o: input_variables_definition-inc.h
#build the fortran include file with the stream of input variables into a local string
#the sed command below will transform the entire yaml document into a fortran string
input_variables_definition-inc.h: input_variables_definition.yaml
	$(SED) -e "s/^/\"/;s/$$/\\\n\"/" $< > $@

tmp-libbigdft/extract.stamp: $(subs_static_library) $(out_static_library)
	test -e "tmp-libbigdft" || $(INSTALL) -d -m 755 tmp-libbigdft
	cd tmp-libbigdft ; \
	rm -f *; \
	for lib in $^; do \
	  $(AR) x ../$$lib ; \
	done ; \
	touch extract.stamp

libbigdft-1.a: $(libmain_a_OBJECTS) tmp-libbigdft/extract.stamp
	$(AR) $(ARFLAGS) $@ $(libmain_a_OBJECTS) tmp-libbigdft/*.o
	$(RANLIB) $@

$(bigdft_dynamic_library): $(libmain_a_OBJECTS) tmp-libbigdft/extract.stamp
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) -Wl,-soname=$(bigdft_library) -o $@ $(libmain_a_OBJECTS) tmp-libbigdft/*.o @LIBS_SHORT_DEPS@ $(amber_ld) $(cp2k_ld)
	@chmod a+x $@

libbigdft-1.so.7: $(bigdft_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libbigdft-1.so

install-data-hook:
	if test -n "$(bigdft_dynamic_library)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(bigdft_dynamic_library) ; \
	  ln -fs $(bigdft_dynamic_library) $(bigdft_library) ; \
	  ln -fs $(bigdft_library) libbigdft-1.so ; \
	fi

#obsolete dependency $(bin_PROGRAMS) $(noinst_PROGRAMS): libbigdft-1.a libdft_common.a

#sources for abscalc exeutable, to be removed from libbigdft.a
#abscalc_additional= 

bigdft_SOURCES = BigDFT.f90
bigdft_LDADD = $(bigdft_library)

NEB_images_SOURCES = BigDFT.f90
NEB_images_LDADD = $(bigdft_library)

splsad_SOURCES = splinedsaddle.f90
splsad_LDADD = $(bigdft_library) @LINALG_LIBS@

#rism_SOURCES = rism.f90 rismlowlevel.f90
#rism_LDADD = $(bigdft_library)

#sandbox_SOURCES = sandbox.f90
#sandbox_LDADD = $(bigdft_library)

#oneatom_SOURCES = oneatom.f90
#oneatom_LDADD = $(bigdft_library)

BigDFT2Wannier_SOURCES = BigDFT2Wannier.f90
BigDFT2Wannier_LDADD = $(bigdft_library) @LINALG_LIBS@

WaCo_SOURCES = WaCo.f90
WaCo_LDADD = $(bigdft_library) @LINALG_LIBS@

MDanalysis_SOURCES = distances.f90
MDanalysis_LDADD = $(bigdft_library)

abscalc_SOURCES = abscalc.f90 
abscalc_FCFLAGS = -I. -I$(top_builddir)/includes -I$(srcdir)/wavelib -I$(top_srcdir)/PSolver/src \
	@LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@ \
	@LIBETSFIO_INCLUDE@ @LIBSGPU_INCLUDE@ \
	@MPI_INCLUDE@
abscalc_LDADD = libabscalc.a $(bigdft_library) @LINALG_LIBS@

test_forces_SOURCES = test_forces.f90
test_forces_LDADD = $(bigdft_library)

NEB_SOURCES = NEB.f90
NEB_LDADD = $(bigdft_library) @LINALG_LIBS@

#to be adjusted: when the compiler has MPI built-in libraries the linker finds multiple definitions
memguess_SOURCES = memguess.f90 $(top_srcdir)/flib/src/MPIfake.f90
memguess_LDADD = $(bigdft_library) @LINALG_LIBS@

utilities_SOURCES = utilities.f90
utilities_LDADD = $(bigdft_library) @LINALG_LIBS@

frequencies_SOURCES = frequencies.f90
frequencies_LDADD = $(bigdft_library) @LINALG_LIBS@

bart_SOURCES = art.f90
bart_FCFLAGS = -I. -Iart -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@
bart_LDADD = art/libart.a $(bigdft_library) @LINALG_LIBS@

mhgps_SOURCES = mhgps.f90
mhgps_FCFLAGS = -I. -Imhgps_files -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@
mhgps_LDADD = mhgps_files/libmhgps.a $(bigdft_library) @LINALG_LIBS@

mhgpstool_SOURCES = mhgpstool.f90
mhgpstool_FCFLAGS = -I. -Imhgps_files -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@
mhgpstool_LDADD = mhgps_files/libmhgps.a $(bigdft_library)  @LINALG_LIBS@

global_SOURCES = global.f90
global_LDADD = $(bigdft_library) @LINALG_LIBS@

globaltool_SOURCES = globaltool.f90
globaltool_FCFLAGS = -I. -Iglobal_files -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@
globaltool_LDADD = global_files/libglobalf.a $(bigdft_library)  @LINALG_LIBS@

CLEANFILES = *.@MODULE_EXT@ bindings_dbus.* introspect libbigdft-1.* input_variables_definition-inc.h
clean-local:
	rm -rf tmp-libbigdft

all: all-recursive
	chmod ug+x bigdft-tool
install: install-recursive
	chmod ug+x bigdft-tool

check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

#dependencies
convolutions_aux.o linear/convolutions_aux.o: linear/filterModule.f90 filterModule.o

lanczos_interface.o: abscalc_module.o lanczos_base.o

abscalc_init.o abscalc-abscalc.o: esatto.o lanczos_interface.o

lanczos.o : lanczos_interface.o

xabsorber.o init/xabsorber.o : init/xabsorber.f90

initAndUtils.o linear/initAndUtils.o: linear/initAndUtils.f90 cleanup.o

linearScaling.o transferintegrals.o: i-o.o

NEB.o: images.o

restart.o: i-o.o

# Introspection for Python bindings
# ---------------------------------
if WITH_GOBJECT_INTROSPECTION
BUILT_GIRSOURCES = BigDFT-1.7.gir

gir_DATA = $(BUILT_GIRSOURCES)
typelibs_DATA = $(BUILT_GIRSOURCES:.gir=.typelib)

if BUILD_BINDINGS
additonal_sources = $(srcdir)/bindings/bindings.c \
		$(srcdir)/bindings/bindings_regions.c \
		$(srcdir)/bindings/bindings_localfields.c \
		$(srcdir)/bindings/bindings_wavefunctions.c \
		$(srcdir)/bindings/bindings_signals.c \
		$(srcdir)/bindings/bindings_signals_inet.c

endif

$(gir_DATA): $(G_IR_SCANNER) $(bigdft_library)
	$(G_IR_SCANNER) -v --warn-all \
		--namespace BigDFT \
		--symbol-prefix=bigdft --identifier-prefix=BigDFT_ --identifier-prefix=Bigdft \
		--symbol-prefix=f90 --identifier-prefix=f90_ \
		--nsversion "1.7" \
		--add-include-path=$(girdir) \
		--include GObject-2.0 \
		--include GLib-2.0 \
		--include Gio-2.0 \
		--output=$(gir_DATA) \
		--library=bigdft-1 \
		--no-libtool \
		-DGLIB_MAJOR_VERSION=2 -L$(builddir) \
		-I$(srcdir)/modules -I$(builddir)/modules \
		$(srcdir)/modules/bigdft.h \
		$(srcdir)/modules/bigdft_input_keys.h \
		$(builddir)/modules/bigdft_cst.h \
		$(srcdir)/bindings/wrappers.c \
		$(srcdir)/bindings/wrapper_inputs.c \
		$(srcdir)/bindings/wrapper_run.c \
		$(srcdir)/bindings/wrapper_neb.c \
		$(srcdir)/bindings/wrapper_atoms.c

$(typelibs_DATA): $(gir_DATA) $(G_IR_COMPILER)
	$(G_IR_COMPILER) \
		--includedir=. \
		--verbose \
		-o $(typelibs_DATA) \
		$(gir_DATA)

CLEANFILES += $(gir_DATA) $(typelibs_DATA)
endif
