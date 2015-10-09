# Not To be installed:
noinst_LIBRARIES = libgeopt.a

mod_src = scfloop_API.f90
routines = \
	sqnm.f90 sdcg.f90 lbfgs.f90 bfgs.f90 geometry.f90 nose_hoover_chains.f90 md_run.f90

libgeopt_a_SOURCES = $(mod_src) $(routines)
EXTRA_DIST = 

# Module files.
if CAPITALIZE
modules = MINPAR.@MODULE_EXT@ SCFLOOP_API.@MODULE_EXT@ NOSE_HOOVER_CHAINS.@MODULE_EXT@
else
modules = minpar.@MODULE_EXT@ scfloop_api.@MODULE_EXT@ nose_hoover_chains.@MODULE_EXT@
endif

CLEANFILES = $(modules)

AM_FCFLAGS = -I. -I$(srcdir) @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes

#dependencies
geometry.o bfgs.o sdcg.o sqnm.o: lbfgs.o

geometry.o: scfloop_API.o
md_run.o: nose_hoover_chains.o

# Copy modules in a common include directory.
all: all-am insmod
install: install-am insmod
check: check-am insmod
insmod: libgeopt.a
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
