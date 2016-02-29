EXTRA_DIST = config/m4

ACLOCAL_AMFLAGS = -I config/m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = psolver.pc

SUBDIRS =

# Low level libraries provided by BigDFT team.
if BUILD_FLIB
SUBDIRS += flib
endif

if BUILD_LIB_PSOLVER
SUBDIRS += wrappers src
endif

if HAVE_PYYAML
SUBDIRS += tests
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py
