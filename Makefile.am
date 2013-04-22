EXTRA_DIST = autogen.sh config/m4 PyYAML-3.10.tar.gz

ACLOCAL_AMFLAGS = -I config/m4

DISTCLEANFILES = config.inc bigdft.pc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bigdft.pc

SUBDIRS =

if BUILD_LIBXC
SUBDIRS += libxc-2.0.x
endif

if BUILD_LIBABINIT
SUBDIRS += libABINIT/src
endif

if BUILD_LIBSGPU
SUBDIRS += S_GPU
endif

if BUILD_LIBYAML
SUBDIRS += yaml-0.1.4
endif

SUBDIRS += \
	flib \
	src \
	doc \
	utils

if BUILD_BINARIES_PSEUDO
SUBDIRS += pseudo
endif

if BUILD_BINARIES
SUBDIRS += tests
endif

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

check report:
	export PYTHONPATH=${PYTHONPATH}; export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ;python $(top_srcdir)/tests/report.py

clean-local:
	-rm -rf includes
