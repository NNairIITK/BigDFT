EXTRA_DIST = config/m4

ACLOCAL_AMFLAGS = -I config/m4

DISTCLEANFILES = config.inc bigdft.pc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bigdft.pc

SUBDIRS =

# Low level libraries provided by BigDFT team.
if BUILD_FLIB
SUBDIRS += flib
endif

if BUILD_LIB_PSOLVER
SUBDIRS += wrappers PSolver
endif

if BUILD_LIB_BIGDFT
# BigDFT itself.
SUBDIRS += src
endif

SUBDIRS += \
	doc \
	utils

# Tests subdir is include only if binaries are compiled.
if BUILD_BINARIES
if HAVE_PYYAML
SUBDIRS += tests
endif
endif

check report:
	python $(top_builddir)/tests/report.py

clean-local:
	-rm -rf includes
