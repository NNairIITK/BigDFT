# Not To be installed:
noinst_LIBRARIES = libsparsematrix.a

# Code source.
if USE_BLACS
BLACS =
else
BLACS = blacs_fake.f90
endif

mod_src = sparsematrix_base.f90 \
          sparsematrix_init.f90 \
          sparsematrix.f90 \
	  sparsematrix_io.f90 \
	  parallel_linalg.f90 \
	  sparsematrix_errorhandling.f90 \
	  sparsematrix_timing.f90 \
	  sparsematrix_types.f90 \
	  sparsematrix_memory.f90 \
	  $(BLACS)

libsparsematrix_a_SOURCES = $(mod_src)

# Module files.
if CAPITALIZE
modules = SPARSEMATRIX_BASE.@MODULE_EXT@ \
          SPARSEMATRIX_INIT.@MODULE_EXT@ \
          SPARSEMATRIX.@MODULE_EXT@ \
	  PARALLEL_LINALG.@MODULE_EXT@ \
	  SPARSEMATRIX_IO.@MODULE_EXT@ \
	  SPARSEMATRIX_ERRORHANDLING.@MODULE_EXT@ \
	  SPARSEMATRIX_TIMING.@MODULE_EXT@ \
	  SPARSEMATRIX_TYPES.@MODULE_EXT@ \
	  SPARSEMATRIX_MEMORY.@MODULE_EXT@
else
modules = sparsematrix_base.@MODULE_EXT@ \
          sparsematrix_init.@MODULE_EXT@ \
          sparsematrix.@MODULE_EXT@ \
	  parallel_linalg.@MODULE_EXT@ \
	  sparsematrix_io.@MODULE_EXT@ \
	  sparsematrix_errorhandling.@MODULE_EXT@ \
	  sparsematrix_timing.@MODULE_EXT@ \
	  sparsematrix_types.@MODULE_EXT@ \
	  sparsematrix_memory.@MODULE_EXT@
endif

CLEANFILES = *.@MODULE_EXT@

AM_FCFLAGS = -I. -I$(srcdir) @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@

#dependencies
sparsematrix_memory.o: sparsematrix_errorhandling.o sparsematrix_types.o
sparsematrix_base.o: sparsematrix_errorhandling.o sparsematrix_timing.o sparsematrix_types.o sparsematrix_memory.o
parallel_linalg.o: sparsematrix_base.o
sparsematrix_init.o: sparsematrix_base.o
sparsematrix.o: sparsematrix_base.o sparsematrix_init.o parallel_linalg.o
sparsematrix_io.o: sparsematrix_base.o sparsematrix.o

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: libsparsematrix.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
