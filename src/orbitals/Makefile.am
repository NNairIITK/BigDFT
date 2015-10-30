# Not To be installed:
noinst_LIBRARIES = liborbitals.a

if USE_OPTI_CONVOLUT
CONV_SOURCES = combined_shrink_optim.f90 \
	combined_grow_optim.f90 \
	growshrink_hyb_optim.f90 \
	convolut_ib_optim.f90 \
	conv_per_optim.f90 \
	convolut_optim_per.f90 \
	convolut_optim_slab.f90
else
CONV_SOURCES = combined_shrink_simple.f90 \
	combined_grow_simple.f90 \
	growshrink_hyb_simple.f90 \
	convolut_ib_simple.f90 \
	conv_per_simple.f90 \
	convolut_simple_per.f90 \
	convolut_simple_slab.f90
endif

mod_src = defs.f90 base.f90 \
	public_keys.f90 \
	locregs.f90 \
	gaussians.f90 \
	locreg_operations.f90

routines = $(CONV_SOURCES) \
	growshrink.f90 \
	growshrink_hyb_common.f90 \
	conv_per_common.f90 \
	convolut_common_slab.f90 \
	convolut_common_per.f90 \
	convolut_new_per.f90 \
	convolut_simple_per_k.f90 \
	convolut_simple_slab_k.f90 \
	scalar_product.f90 \
	potential.f90 \
	locregs_basics.f90 \
	precondLinear.f90 \
	daubisf.f90

liborbitals_a_SOURCES = $(mod_src) $(routines)
# convolut_ib_sse.c \ #to be verified, not working yet
EXTRA_DIST = configure.inc.in v.inc v_long.inc v_17.inc scalar_product-inc.f90 wpdot-inc.f90

# Module files.
if CAPITALIZE
modules = MODULE_BASE.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ \
	PUBLIC_KEYS.@MODULE_EXT@ \
	PUBLIC_ENUMS.@MODULE_EXT@ \
	GAUSSIANS.@MODULE_EXT@ \
	LOCREG_OPERATIONS.@MODULE_EXT@ \
	LOCREGS.@MODULE_EXT@
else
modules = module_base.@MODULE_EXT@ \
	module_defs.@MODULE_EXT@ \
	public_keys.@MODULE_EXT@ \
	public_enums.@MODULE_EXT@ \
	gaussians.@MODULE_EXT@ \
	locreg_operations.@MODULE_EXT@ \
	locregs.@MODULE_EXT@
endif

CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes

# Hack for MIC compilation.
if USE_MIC
defs.o: defs.f90 configure.inc
	$(PPFCCOMPILE) -offload-attribute-target=mic -opt-report-phase:offload -offload-option,mic,compiler,'-g -O0' -c -o $@ $<
else
defs.o: configure.inc 
endif

#dependencies
base.o: defs.o public_keys.o
locregs.o: base.o
gaussians.o: base.o
locreg_operations.o: locregs.o
potential.o: locreg_operations.o
daubisf.o precondLinear.o locregs_basics.o: locreg_operations.o
$(routines:.f90=.o): defs.o

convolut_common_per.o conv_per_common.o conv_per_optim.o: base.o

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: liborbitals.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
