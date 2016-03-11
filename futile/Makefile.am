EXTRA_DIST = config/m4

ACLOCAL_AMFLAGS = -I config/m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = futile.pc

SUBDIRS = flib wrappers src

if HAVE_PYYAML
SUBDIRS += tests
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(top_builddir)/tests/report.py
