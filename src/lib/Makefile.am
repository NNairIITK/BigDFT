# To be installed:
noinst_LIBRARIES = liblib.a
module_DATA = $(modules)

# Source code.
if USE_DGEMMSY
DGEMMSY_SOURCES = dgemmsy/interface_dgemmsy.f90 \
	dgemmsy/dgemmsy.c       dgemmsy/dgemmsy.h \
	dgemmsy/gemm_block_c.c  dgemmsy/gemm_block_c.h \
	dgemmsy/patterns.c      dgemmsy/patterns.h \
	dgemmsy/dgemmsy_utils.c dgemmsy/utils.h \
	dgemmsy/visitors.c      dgemmsy/visitors.h
else
DGEMMSY_SOURCES = dgemmsy/interface_dgemmsy_fake.f90
endif
liblib_a_SOURCES = $(DGEMMSY_SOURCES) \
	fft/fft2d.f90 \
	fft/fft3d.f90 \
	lbfgs.f90

# Dependencies
$(modules): $(liblib_a_OBJECTS)

AM_FCFLAGS = -I$(top_builddir)/includes @LIBABINIT_INCLUDE@
AM_CPPFLAGS = @DGEMMSY_CPPFLAGS@

# Module files.
if CAPITALIZE
modules = MODULE_FFT_SG.@MODULE_EXT@ MINPAR.@MODULE_EXT@
else
modules = module_fft_sg.@MODULE_EXT@ minpar.@MODULE_EXT@
endif
CLEANFILES = $(modules)

# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done