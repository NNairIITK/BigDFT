# Not To be installed:
noinst_LIBRARIES = libchesspexsi.a
module_DATA      = $(modules)

# Code source.
if HAVE_PEXSI
PEXSI_SRC = pexsi_base.f90 pexsi_interfaces.f90 selinv.f90 pexsi.f90
else
PEXSI_SRC = pexsi_base.f90 pexsi_interfaces_fake.f90 selinv.f90 pexsi.f90
endif

mod_src = $(PEXSI_SRC)

libchesspexsi_a_SOURCES = $(mod_src)

# Module files.
if CAPITALIZE
modules = PEXSI_BASE.@MODULE_EXT@ \
	PEXSI_INTERFACES.@MODULE_EXT@ \
	SELINV.@MODULE_EXT@ \
	PEXSI.@MODULE_EXT@
else
modules = pexsi_base.@MODULE_EXT@ \
	pexsi_interfaces.@MODULE_EXT@ \
	selinv.@MODULE_EXT@ \
	pexsi.@MODULE_EXT@

endif

CLEANFILES = *.@MODULE_EXT@

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@

#dependencies
$(modules): $(libchesspexsi_a_OBJECTS)
if HAVE_PEXSI
pexsi_interfaces.o: pexsi_base.o
selinv.o: pexsi_base.o pexsi_interfaces.o
pexsi.o: pexsi_base.o pexsi_interfaces.o
else
pexsi_interfaces_fake.o: pexsi_base.o
selinv.o: pexsi_base.o pexsi_interfaces_fake.o
pexsi.o: pexsi_base.o pexsi_interfaces_fake.o
endif

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: libchesspexsi.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
