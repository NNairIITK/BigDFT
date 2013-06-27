# To be installed:
noinst_LIBRARIES = libmodules.a
module_DATA = $(modules)
include_HEADERS = $(bindings_headers)
nodist_include_HEADERS = $(bindings_headers_nodist)

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

if BUILD_BINDINGS
bindings_headers = bigdft.h
bindings_headers_nodist = bigdft_cst.h
else
bindings_headers =
bindings_headers_nodist =
endif

libmodules_a_SOURCES = base.f90 \
	defs.f90 \
	gaussians.f90 \
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
	DIIS_SD_OPTIMIZATION.@MODULE_EXT@
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
	diis_sd_optimization.@MODULE_EXT@
endif


EXTRA_DIST = configure.inc.in bigdft_cst.h.in BigDFT.py
CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

# Hack for XLF strange DEFS option.
PPFCCOMPILE = $(FC) $(DEFS:-D%=@FCDEFS@%) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

# Hack for MIC compilation.
if USE_MIC
defs.o: defs.F90 configure.inc
	$(PPFCCOMPILE) -offload-attribute-target=mic -opt-report-phase:offload -offload-option,mic,compiler,'-g -O0' -c -o $@ $<
cublas_fake.o: cublas_fake.f90
	$(FCCOMPILE) -offload-attribute-target=mic -c -o $@ $<
else
defs.o: configure.inc 
endif

#dependencies
$(modules): $(libmodules_a_OBJECTS)
base.o: defs.o $(top_builddir)/flib/src/dynamic_memory.o
types.o: base.o gaussians.o
input.o: base.o 
interfaces.o: base.o types.o xc.o fragments.o diis_module.o constrained_dft.o
gaussians.o: base.o
op2p_module.o : base.o
BigDFT_API.o: base.o types.o interfaces.o xc.o
vdwcorrection.o: base.o types.o
private_api.o: base.o types.o interfaces.o
fragments.o diis_module.o: base.o types.o
constrained_dft.o: base.o types.o fragments.o

# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules) 
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
#	test -z "$(mpi_include)" || $(INSTALL) -m 644 mpif.h $(top_builddir)/includes
