EXTRA_DIST = autogen.sh config/m4 PyYAML-3.10.tar.gz

ACLOCAL_AMFLAGS = -I config/m4

DISTCLEANFILES = config.inc bigdft.pc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bigdft.pc

SUBDIRS =

# This is required for testing.
if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
SUBDIRS += yaml-0.1.4
endif

# Low level libraries provided by BigDFT team.
SUBDIRS += flib

if BUILD_LIB_PSOLVER
SUBDIRS += wrappers PSolver
endif

if BUILD_LIB_BIGDFT
# All BigDFT dependencies.
if BUILD_LIBXC
SUBDIRS += libxc-2.0.x
endif
if BUILD_LIBABINIT
SUBDIRS += libABINIT/src
endif
if BUILD_LIBSGPU
SUBDIRS += S_GPU
endif

# BigDFT itself.
SUBDIRS += src

# Additional tools around BigDFT.
if BUILD_BINARIES_PSEUDO
SUBDIRS += pseudo
endif
endif

SUBDIRS += \
	doc \
	utils

# Required by tests
if BUILD_LIBYAML
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
SUBDIRS += PyYAML-3.10
endif

# Tests subdir is include only if binaries are compiled.
if BUILD_BINARIES
SUBDIRS += tests
endif

check report:
	export PYTHONPATH=${PYTHONPATH}; export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}; python $(top_srcdir)/tests/report.py

clean-local:
	-rm -rf includes
