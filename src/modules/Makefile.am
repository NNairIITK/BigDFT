# To be installed:
noinst_LIBRARIES = libmodules.a
module_DATA = $(modules)
include_HEADERS = bigdft.h bigdft_input_keys.h
nodist_include_HEADERS = bigdft_cst.h

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

if USE_LIBXC
LIBXC_MODULE =
xc.o: base.o
else
LIBXC_MODULE = libxc_fake.f90
xc.o: base.o libxc_fake.o
endif

libmodules_a_SOURCES = base.f90 \
	defs.f90 \
	gaussians.f90 \
	locregs.f90 \
	psp_projectors.f90 \
	input_keys.f90 \
	input_dicts.f90 \
	types.f90 \
	input.f90 \
	interfaces.f90 \
	private_api.f90 \
	xc.f90 \
	op2p_module.f90 \
	BigDFT_API.f90 \
	vdwcorrection.f90 \
	fragments.f90 \
	constrained_dft.f90 \
	diis_module.f90 \
	ao_inguess.f90 \
	atoms_data.f90 \
	communications_base.f90 \
	communications_init.f90 \
	communications.f90 \
	sparsematrix_base.f90 \
	sparsematrix_init.f90 \
	sparsematrix.f90 \
	foe_base.f90 \
	internal_coordinates.f90 \
	$(bindings_headers) \
	$(LIBXC_MODULE) \
	$(BLACS) \
	configure.c

# Modules
if CAPITALIZE
modules = MODULE_BASE.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@ \
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
	MODULE_ATOMS.@MODULE_EXT@ \
	COMMUNICATIONS_BASE.@MODULE_EXT@ \
	COMMUNICATIONS_INIT.@MODULE_EXT@ \
	COMMUNICATIONS.@MODULE_EXT@ \
	SPARSEMATRIX_BASE.@MODULE_EXT@ \
	SPARSEMATRIX_INIT.@MODULE_EXT@ \
	SPARSEMATRIX.@MODULE_EXT@ \
	FOE_BASE.@MODULE_EXT@ \
	INTERNAL_COORDINATES.@MODULE_EXT@
else
modules = module_base.@MODULE_EXT@ \
	module_defs.@MODULE_EXT@ \
	module_types.@MODULE_EXT@ \
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
	module_atoms.@MODULE_EXT@ \
	communications_base.@MODULE_EXT@ \
	communications_init.@MODULE_EXT@ \
	communications.@MODULE_EXT@ \
	sparsematrix_base.@MODULE_EXT@ \
	sparsematrix_init.@MODULE_EXT@ \
	sparsematrix.@MODULE_EXT@ \
	foe_base.@MODULE_EXT@ \
	internal_coordinates.@MODULE_EXT@
else


endif

AO_INCLUDES = eleconf-inc.f90
ATOMS_INCLUDES = astruct-inc.f90
COMMS_INCLUDES = check_array-inc.f90

EXTRA_DIST = configure.inc.in bigdft_cst.h.in BigDFT.py \
	$(AO_INCLUDES) $(ATOMS_INCLUDES) $(COMMS_INCLUDES)
CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

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
input_keys.o: defs.o $(top_builddir)/flib/src/dynamic_memory.o base.o
base.o: defs.o $(top_builddir)/flib/src/dynamic_memory.o
types.o: base.o gaussians.o locregs.o psp_projectors.o ao_inguess.o atoms_data.o communications_base.o \
	sparsematrix_base.o xc.o foe_base.o
input.o: base.o types.o input_keys.o
interfaces.o: base.o types.o xc.o fragments.o diis_module.o constrained_dft.o sparsematrix_base.o foe_base.o
gaussians.o: base.o
op2p_module.o : base.o
BigDFT_API.o: base.o types.o interfaces.o xc.o
vdwcorrection.o: base.o types.o
private_api.o: base.o types.o interfaces.o
diis_module.o: base.o types.o communications.o sparsematrix_base.o sparsematrix.o
fragments.o: base.o types.o
constrained_dft.o: base.o types.o sparsematrix_base.o
input_dicts.o: $(top_builddir)/flib/src/dynamic_memory.o types.o defs.o ao_inguess.o interfaces.o atoms_data.o input_keys.o
psp_projectors.o: locregs.o
locregs.o: base.o
ao_inguess.o: $(AO_INCLUDES) base.o psp_projectors.o
atoms_data.o: $(ATOMS_INCLUDES) ao_inguess.o
communications_base.o: base.o
communications_init.o: base.o types.o communications_base.o
communications.o: base.o types.o communications_base.o
sparsematrix_base.o: base.o
sparsematrix_init.o: base.o types.o interfaces.o sparsematrix_base.o
sparsematrix.o: base.o types.o sparsematrix_base.o
foe_base.o: base.o defs.o
internat_coodinates.o:

# Copy modules in a common include directory.
all: all-am insmod
check: check-am insmod
insmod: $(modules) 
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
#	test -z "$(mpi_include)" || $(INSTALL) -m 644 mpif.h $(top_builddir)/includes
