# Not To be installed:
noinst_LIBRARIES = libsparsematrix.a

mod_src = sparsematrix_base.f90 \
          sparsematrix_init.f90 \
          sparsematrix.f90

libsparsematrix_a_SOURCES = $(mod_src)


# Module files.
if CAPITALIZE
modules = SPARSEMATRIX_BASE.@MODULE_EXT@ \
          SPARSEMATRIX_INIT.@MODULE_EXT@ \
          SPARSEMATRIX.@MODULE_EXT@
else
modules = sparsematrix_base.@MODULE_EXT@ \
          sparsematrix_init.@MODULE_EXT@ \
          sparsematrix.@MODULE_EXT@
endif

CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes


#dependencies
sparsematrix_base.o:
sparsematrix_init.o: sparsematrix_base.o
sparsematrix.o: sparsematrix_base.o sparsematrix_init.o
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
