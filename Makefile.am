EXTRA_DIST = autogen.sh config/m4 PyYAML-3.10.tar.gz

ACLOCAL_AMFLAGS = -I config/m4

DISTCLEANFILES = config.inc bigdft.pc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bigdft.pc

SUBDIRS =

if BUILD_LIBXC
SUBDIRS += libxc-1.1.0
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

SUBDIRS += src \
	   doc \
	   utils

if BUILD_BINARIES_PSEUDO
SUBDIRS += pseudo
endif

if BUILD_BINARIES
SUBDIRS += tests
endif

check:
	@python $(top_srcdir)/tests/report.py

report:
	@python $(top_srcdir)/tests/report.py
