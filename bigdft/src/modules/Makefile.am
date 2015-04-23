# To be installed:
overridesdir = $(pyexecdir)/gi/overrides
noinst_LIBRARIES = libmodules.a
module_DATA = $(modules)
include_HEADERS = bigdft.h bigdft_input_keys.h
nodist_include_HEADERS = bigdft_cst.h
overrides_DATA = BigDFT.py

# Code source.
if USE_BLACS
BLACS = 
else
BLACS = blacs_fake.f90
endif

#if USE_MPI
#mpi_include =
#else
#mpi_include = mpif.h
#endif

libmodules_a_SOURCES = base.f90 \
	defs.f90 \
	public_keys.f90 \
	ab7_mixing.f90 \
	gaussians.f90 \
	locregs.f90 \
	psp_projectors.f90 \
	input_keys.f90 \
	gaussdaub.f90 \
	input_dicts.f90 \
	types.f90 \
	input.f90 \
	interfaces.f90 \
	private_api.f90 \
	xc.f90 \
	module_lj.f90 \
	module_morse_bulk.f90 \
	module_lenosky_si.f90 \
	module_tersoff.f90 \
	module_bornmayerhugginstosifumi.f90 \
	op2p_module.f90 \
	BigDFT_API.f90 \
	vdwcorrection.f90 \
	fragments.f90 \
	constrained_dft.f90 \
	diis_module.f90 \
	ao_inguess.f90 \
	reformatting.f90 \
	atoms_data.f90 \
	communications_base.f90 \
	communications_init.f90 \
	communications.f90 \
	sparsematrix_base.f90 \
	sparsematrix_init.f90 \
	sparsematrix.f90 \
	foe_base.f90 \
	foe_common.f90 \
	internal_coordinates.f90 \
	fermi_level.f90 \
	copy_utils.f90 \
	module_sqn.f90 \
	module_fingerprints.f90 \
	transposed_operations.f90 \
	module_forces.f90 \
	chebyshev.f90 \
	unitary_tests.f90 \
	matrix_operations.f90 \
	parallel_linalg.f90 \
	io.f90 \
	postprocessing_linear.f90 \
	foe.f90 \
	ice.f90 \
	locreg_operations.f90 \
	bigdft_run.f90 \
	rhopotential.f90 \
	$(bindings_headers) \
	$(BLACS) \
	configure.c

# Modules
if CAPITALIZE
modules = MODULE_BASE.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@ \
	MODULE_MIXING.@MODULE_EXT@ \
	PUBLIC_KEYS.@MODULE_EXT@ \
	PUBLIC_ENUMS.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@ \
	MODULE_PRIVATE_API.@MODULE_EXT@ \
	MODULE_XC.@MODULE_EXT@ \
	OVERLAP_POINT_TO_POINT.@MODULE_EXT@ \
	MODULE_INPUT.@MODULE_EXT@ \
	GAUSSIANS.@MODULE_EXT@ \
	BIGDFT_API.@MODULE_EXT@ \
	VDWCORRECTION.@MODULE_EXT@ \
	MODULE_FRAGMENTS.@MODULE_EXT@ \
	CONSTRAINED_DFT.@MODULE_EXT@ \
	DIIS_SD_OPTIMIZATION.@MODULE_EXT@ \
	INPUT_OLD_TEXT_FORMAT.@MODULE_EXT@ \
	MODULE_INPUT_KEYS.@MODULE_EXT@ \
	MODULE_INPUT_DICTS.@MODULE_EXT@ \
	PSP_PROJECTORS.@MODULE_EXT@ \
	LOCREGS.@MODULE_EXT@ \
	AO_INGUESS.@MODULE_EXT@ \
	REFORMATTING.@MODULE_EXT@ \
	MODULE_ATOMS.@MODULE_EXT@ \
	COMMUNICATIONS_BASE.@MODULE_EXT@ \
	COMMUNICATIONS_INIT.@MODULE_EXT@ \
	COMMUNICATIONS.@MODULE_EXT@ \
	SPARSEMATRIX_BASE.@MODULE_EXT@ \
	SPARSEMATRIX_INIT.@MODULE_EXT@ \
	SPARSEMATRIX.@MODULE_EXT@ \
	FOE_BASE.@MODULE_EXT@ \
	FOE_COMMON.@MODULE_EXT@ \
	INTERNAL_COORDINATES.@MODULE_EXT@ \
	FERMI_LEVEL.@MODULE_EXT@ \
	COPY_UTILS.@MODULE_EXT@ \
	MODULE_SQN.@MODULE_EXT@ \
	MODULE_FORCES.@MODULE_EXT@ \
	MODULE_FINGERPRINTS.@MODULE_EXT@ \
	TRANSPOSED_OPERATIONS.@MODULE_EXT@ \
	CHEBYSHEV.@MODULE_EXT@ \
	UNITARY_TESTS.@MODULE_EXT@ \
	MATRIX_OPERATIONS.@MODULE_EXT@ \
	PARALLEL_LINALG.@MODULE_EXT@ \
	IO.@MODULE_EXT@ \
	POSTPROCESSING_LINEAR.@MODULE_EXT@ \
	FOE.@MODULE_EXT@ \
	ICE.@MODULE_EXT@ \
	LOCREG_OPERATIONS.@MODULE_EXT@ \
	BIGDFT_RUN.@MODULE_EXT@ \
	RHOPOTENTIAL.@MODULE_EXT@ \
	MODULE_LJ.@MODULE_EXT@ \
	MODULE_LENOSKY_SI.@MODULE_EXT@\
	MODULE_MORSE_BULK.@MODULE_EXT@\
	MODULE_TERSOFF.@MODULE_EXT@\
	MODULE_BORNMAYERHUGGINSTOSIFUMI.@MODULE_EXT@
else
modules = module_base.@MODULE_EXT@ \
	module_defs.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ \
	module_mixing.@MODULE_EXT@ \
	public_keys.@MODULE_EXT@ \
	public_enums.@MODULE_EXT@ \
	module_interfaces.@MODULE_EXT@ \
	module_private_api.@MODULE_EXT@ \
	module_xc.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@ \
	module_input.@MODULE_EXT@ \
	gaussians.@MODULE_EXT@ \
	bigdft_api.@MODULE_EXT@ \
	vdwcorrection.@MODULE_EXT@ \
	module_fragments.@MODULE_EXT@ \
	constrained_dft.@MODULE_EXT@ \
	diis_sd_optimization.@MODULE_EXT@ \
	input_old_text_format.@MODULE_EXT@ \
	module_input_keys.@MODULE_EXT@ \
	module_input_dicts.@MODULE_EXT@ \
	psp_projectors.@MODULE_EXT@ \
	locregs.@MODULE_EXT@ \
	ao_inguess.@MODULE_EXT@ \
	reformatting.@MODULE_EXT@ \
	module_atoms.@MODULE_EXT@ \
	communications_base.@MODULE_EXT@ \
	communications_init.@MODULE_EXT@ \
	communications.@MODULE_EXT@ \
	sparsematrix_base.@MODULE_EXT@ \
	sparsematrix_init.@MODULE_EXT@ \
	sparsematrix.@MODULE_EXT@ \
	foe_base.@MODULE_EXT@ \
	foe_common.@MODULE_EXT@ \
	internal_coordinates.@MODULE_EXT@ \
	fermi_level.@MODULE_EXT@ \
	copy_utils.@MODULE_EXT@ \
	module_sqn.@MODULE_EXT@ \
	module_forces.@MODULE_EXT@ \
	module_fingerprints.@MODULE_EXT@ \
	transposed_operations.@MODULE_EXT@ \
	chebyshev.@MODULE_EXT@ \
	unitary_tests.@MODULE_EXT@ \
	matrix_operations.@MODULE_EXT@ \
	parallel_linalg.@MODULE_EXT@ \
	io.@MODULE_EXT@ \
	postprocessing_linear.@MODULE_EXT@ \
	foe.@MODULE_EXT@ \
	ice.@MODULE_EXT@ \
	locreg_operations.@MODULE_EXT@ \
	bigdft_run.@MODULE_EXT@ \
	rhopotential.@MODULE_EXT@ \
	module_lj.@MODULE_EXT@ \
	module_morse_bulk.@MODULE_EXT@ \
	module_tersoff.@MODULE_EXT@ \
	module_bornmayerhugginstosifumi.@MODULE_EXT@ \
	module_lenosky_si.@MODULE_EXT@
endif

AO_INCLUDES = eleconf-inc.f90 
ATOMS_INCLUDES = astruct-inc.f90
COMMS_INCLUDES = check_array-inc.f90

EXTRA_DIST = configure.inc.in bigdft_cst.h.in \
	$(AO_INCLUDES) $(ATOMS_INCLUDES) $(COMMS_INCLUDES) $(overrides_DATA)
CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ @LIB_FLIB_CFLAGS@

# Hack for XLF strange DEFS option.
PPFCCOMPILE = $(FC) $(DEFS:-D%=@FCDEFS@%) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

# Hack for MIC compilation.
if USE_MIC
defs.o: defs.f90 configure.inc
	$(PPFCCOMPILE) -offload-attribute-target=mic -opt-report-phase:offload -offload-option,mic,compiler,'-g -O0' -c -o $@ $<
cublas_fake.o: cublas_fake.f90
	$(FCCOMPILE) -offload-attribute-target=mic -c -o $@ $<
else
defs.o: configure.inc 
endif

#dependencies
$(modules): $(libmodules_a_OBJECTS)
input_keys.o: defs.o base.o public_keys.o
base.o: defs.o public_keys.o
types.o: base.o gaussians.o locregs.o psp_projectors.o ao_inguess.o atoms_data.o communications_base.o \
	sparsematrix_base.o xc.o foe_base.o ab7_mixing.o
input.o: base.o types.o input_keys.o
interfaces.o: base.o types.o xc.o fragments.o diis_module.o constrained_dft.o sparsematrix_base.o foe_base.o
gaussians.o: base.o
op2p_module.o : base.o
BigDFT_API.o: base.o types.o interfaces.o xc.o atoms_data.o input_dicts.o ao_inguess.o psp_projectors.o communications_base.o communications_init.o bigdft_run.o public_keys.o gaussians.o
vdwcorrection.o: base.o types.o
private_api.o: base.o types.o interfaces.o
diis_module.o: base.o types.o
fragments.o: base.o types.o locreg_operations.o 
constrained_dft.o: base.o types.o sparsematrix_base.o
input_dicts.o: defs.o ao_inguess.o interfaces.o atoms_data.o input.o public_keys.o
psp_projectors.o: locregs.o gaussians.o defs.o
locregs.o: base.o
ao_inguess.o: $(AO_INCLUDES) base.o psp_projectors.o psp_projectors.o
atoms_data.o: $(ATOMS_INCLUDES) ao_inguess.o internal_coordinates.o
reformatting.o: base.o
communications_base.o: base.o
communications_init.o: base.o types.o communications_base.o
communications.o: base.o types.o communications_base.o locregs.o locreg_operations.o xc.o
sparsematrix_base.o: base.o
sparsematrix_init.o: base.o types.o interfaces.o sparsematrix_base.o locregs.o 
sparsematrix.o: base.o types.o sparsematrix_base.o sparsematrix_init.o
foe_base.o: base.o defs.o
foe_common.o: base.o defs.o foe_base.o sparsematrix_base.o sparsematrix_init.o sparsematrix.o
fermi_level.o: base.o foe_common.o
copy_utils.o: base.o
gaussdaub.o: base.o
xc.o: base.o
module_sqn.o: base.o
bigdft_run.o: types.o input_dicts.o interfaces.o module_forces.o module_lenosky_si.o module_lj.o module_morse_bulk.o module_tersoff.o module_bornmayerhugginstosifumi.o
ab7_mixing.o: defs.o
module_lj.o: base.o defs.o
module_morse_bulk.o: base.o
module_tersoff.o: base.o atoms_data.o
module_bornmayerhugginstosifumi.o: base.o atoms_data.o
module_lenosky_si.o: base.o defs.o
module_fingerprints.o: base.o
transposed_operations.o: base.o types.o sparsematrix_base.o sparsematrix.o
module_forces.o: base.o types.o 
chebyshev.o: base.o types.o sparsematrix_base.o sparsematrix_init.o sparsematrix.o
unitary_tests.o: base.o types.o interfaces.o communications_init.o communications.o sparsematrix_base.o sparsematrix_init.o sparsematrix.o transposed_operations.o locregs.o locreg_operations.o rhopotential.o
matrix_operations.o: base.o types.o interfaces.o sparsematrix_base.o sparsematrix_init.o sparsematrix.o parallel_linalg.o ice.o
parallel_linalg.o: base.o types.o
io.o: base.o types.o interfaces.o sparsematrix_base.o sparsematrix.o
postprocessing_linear.o: base.o types.o interfaces.o communications_base.o communications.o sparsematrix_base.o sparsematrix.o transposed_operations.o matrix_operations.o rhopotential.o
foe.o: base.o types.o interfaces.o sparsematrix_base.o sparsematrix_init.o sparsematrix.o foe_base.o foe_common.o fermi_level.o chebyshev.o matrix_operations.o
ice.o: base.o types.o interfaces.o sparsematrix_base.o sparsematrix_init.o sparsematrix.o foe_base.o foe_common.o fermi_level.o chebyshev.o
locreg_operations.o: base.o types.o atoms_data.o locregs.o
rhopotential.o: base.o types.o interfaces.o xc.o communications_base.o communications.o locreg_operations.o sparsematrix_base.o

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: $(modules) 
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
#	test -z "$(mpi_include)" || $(INSTALL) -m 644 mpif.h $(top_builddir)/includes
