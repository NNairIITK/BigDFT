# Not To be installed:
noinst_LIBRARIES = libsparsematrix.a

# Code source.
if USE_BLACS
BLACS =
else
BLACS = blacs_fake.f90
endif
if HAVE_CP2K
BLACS =
endif
mod_src = sparsematrix_base.f90 \
          sparsematrix_init.f90 \
          sparsematrix.f90 \
		  chebyshev.f90 \
		  foe_base.f90 \
		  foe_common.f90 \
		  fermi_level.f90 \
		  foe.f90 \
		  ice.f90 \
		  matrix_operations.f90 \
		  sparsematrix_highlevel.f90 \
		  sparsematrix_io.f90 \
		  parallel_linalg.f90 \
		  $(BLACS)

libsparsematrix_a_SOURCES = $(mod_src)


# Module files.
if CAPITALIZE
modules = SPARSEMATRIX_BASE.@MODULE_EXT@ \
          SPARSEMATRIX_INIT.@MODULE_EXT@ \
          SPARSEMATRIX.@MODULE_EXT@ \
		  CHEBYSHEV.@MODULE_EXT@ \
		  FOE_BASE.@MODULE_EXT@ \
		  FOE_COMMON.@MODULE_EXT@ \
		  FERMI_LEVEL.@MODULE_EXT@ \
		  FOE.@MODULE_EXT@ \
	MODULE_FUNC.@MODULE_EXT@ \
		  ICE.@MODULE_EXT@ \
		  MATRIX_OPERATIONS.@MODULE_EXT@ \
		  PARALLEL_LINALG.@MODULE_EXT@ \
		  SPARSEMATRIX_HIGHLEVEL.@MODULE_EXT@ \
		  SPARSEMATRIX_IO.@MODULE_EXT@
else
modules = sparsematrix_base.@MODULE_EXT@ \
          sparsematrix_init.@MODULE_EXT@ \
          sparsematrix.@MODULE_EXT@ \
		  chebyshev.@MODULE_EXT@ \
		  foe_base.@MODULE_EXT@ \
		  foe_common.@MODULE_EXT@ \
		  fermi_level.@MODULE_EXT@ \
		  foe.@MODULE_EXT@ \
	module_func.@MODULE_EXT@ \
		  ice.@MODULE_EXT@ \
		  matrix_operations.@MODULE_EXT@ \
		  parallel_linalg.@MODULE_EXT@ \
		  sparsematrix_highlevel.@MODULE_EXT@ \
		  sparsematrix_io.@MODULE_EXT@
endif

CLEANFILES = *.@MODULE_EXT@

AM_FCFLAGS = -I. -I$(srcdir) @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes


#dependencies
parallel_linalg.o:
foe_base.o:
sparsematrix_base.o:
sparsematrix_init.o: sparsematrix_base.o
sparsematrix.o: sparsematrix_base.o sparsematrix_init.o
chebyshev.o: sparsematrix_base.o sparsematrix_init.o sparsematrix.o
foe_common.o: sparsematrix_base.o sparsematrix_init.o sparsematrix.o foe_base.o chebyshev.o fermi_level.o
ice.o: sparsematrix_base.o sparsematrix_init.o sparsematrix.o chebyshev.o foe_base.o foe_common.o fermi_level.o
matrix_operations.o: sparsematrix_base.o sparsematrix_init.o sparsematrix.o parallel_linalg.o ice.o
foe.o: sparsematrix_base.o sparsematrix.o chebyshev.o foe_base.o foe_common.o fermi_level.o matrix_operations.o
sparsematrix_highlevel.o: sparsematrix_base.o sparsematrix_init.o sparsematrix.o sparsematrix_io.o ice.o foe.o
sparsematrix_io.o: sparsematrix_base.o sparsematrix.o
$(routines:.f90=.o): defs.o


# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: libsparsematrix.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
