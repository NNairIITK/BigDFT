# To be installed:
noinst_LIBRARIES = liblib.a
module_DATA = $(modules)

liblib_a_SOURCES = lbfgs.f90

# Dependencies
$(modules): $(liblib_a_OBJECTS)

AM_FCFLAGS = -I$(top_builddir)/includes @LIBABINIT_INCLUDE@

# Module files.
if CAPITALIZE
modules = MINPAR.@MODULE_EXT@
else
modules = minpar.@MODULE_EXT@
endif
CLEANFILES = $(modules)

# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done