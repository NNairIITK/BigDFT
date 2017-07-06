# Not To be installed:
noinst_LIBRARIES = liblocregs.a
module_DATA = $(modules)

sources = defs.f90 \
	base.f90 \
	compression.f90 \
	locregs.f90 \
	locreg_operations.f90 \
	locregs_basics.f90 \
	bounds.f90 \
	scalar_product.f90

liblocregs_a_SOURCES = $(sources)

EXTRA_DIST = configure.inc.in scalar_product-inc.f90 wpdot-inc.f90

# Module files.
if CAPITALIZE
modules = MODULE_BASE.@MODULE_EXT@ \
	MODULE_DEFS.@MODULE_EXT@ \
	BOUNDS.@MODULE_EXT@ \
	LOCREG_OPERATIONS.@MODULE_EXT@ \
	LOCREGS.@MODULE_EXT@ \
	COMPRESSION.@MODULE_EXT@ 
else
modules = module_base.@MODULE_EXT@ \
	module_defs.@MODULE_EXT@ \
	bounds.@MODULE_EXT@ \
	locreg_operations.@MODULE_EXT@ \
	locregs.@MODULE_EXT@ \
	compression.@MODULE_EXT@ 
endif

CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @LIB_FUTILE_CFLAGS@

check_PROGRAMS = locreg_test

locreg_test_SOURCES = locreg_test.f90
locreg_test_LDADD = -L$(builddir) -llocregs @LIB_FUTILE_LIBS@ @LINALG_LIBS@

# Hack for MIC compilation.
if USE_MIC
defs.o: defs.f90 configure.inc
	$(PPFCCOMPILE) -offload-attribute-target=mic -opt-report-phase:offload -offload-option,mic,compiler,'-g -O0' -c -o $@ $<
else
defs.o: configure.inc 
endif

#dependencies
locregs.o: base.o compression.o bounds.o
bounds.o: compression.o
locreg_operations.o: locregs.o
locregs_basics.o: locreg_operations.o
$(sources:.f90=.o): defs.o
locreg_test.o: liblocregs.a

$(modules): $(liblocregs_a_OBJECTS)

#liblocregs.a: $(liblocregs_a_OBJECTS) insmod

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: $(liblocregs_a_OBJECTS)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
