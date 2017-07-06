# Not To be installed:
noinst_LIBRARIES = liborbitals.a
module_DATA = $(modules)

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

mod_src = public_keys.f90 \
	gaussians.f90 \
	ao_inguess.f90 \
	xc.f90 \
	pseudopotentials.f90 \
	communications_base.f90

routines = $(CONV_SOURCES) \
	growshrink.f90 \
	growshrink_hyb_common.f90 \
	conv_per_common.f90 \
	convolut_common_slab.f90 \
	convolut_common_per.f90 \
	convolut_new_per.f90 \
	convolut_simple_per_k.f90 \
	convolut_simple_slab_k.f90 \
	potential.f90 \
	precondLinear.f90 \
	daubisf.f90

liborbitals_a_SOURCES = $(mod_src) $(routines)

AO_INCLUDES = eleconf-inc.f90 

# convolut_ib_sse.c \ #to be verified, not working
EXTRA_DIST = v.inc v_long.inc v_17.inc $(AO_INCLUDES)

# Module files.
if CAPITALIZE
modules = PUBLIC_KEYS.@MODULE_EXT@ \
	PUBLIC_ENUMS.@MODULE_EXT@ \
	GAUSSIANS.@MODULE_EXT@ \
	AO_INGUESS.@MODULE_EXT@ \
	MODULE_XC.@MODULE_EXT@ \
	PSEUDOPOTENTIALS.@MODULE_EXT@ \
	COMMUNICATIONS_BASE.@MODULE_EXT@
else
modules = public_keys.@MODULE_EXT@ \
	public_enums.@MODULE_EXT@ \
	gaussians.@MODULE_EXT@ \
	ao_inguess.@MODULE_EXT@ \
	module_xc.@MODULE_EXT@ \
	pseudopotentials.@MODULE_EXT@ \
	communications_base.@MODULE_EXT@
endif

CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ @LIB_FUTILE_CFLAGS@

#dependencies
ao_inguess.o: $(AO_INCLUDES) public_keys.o
pseudopotentials.o: ao_inguess.o public_keys.o xc.o

$(modules): $(liborbitals_a_OBJECTS)

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: liborbitals.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
