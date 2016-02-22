EXTRA_DIST = config/m4

ACLOCAL_AMFLAGS = -I config/m4

SUBDIRS =

# Low level libraries provided by BigDFT team.
if BUILD_FLIB
SUBDIRS += flib
endif

SUBDIRS += wrappers src

if HAVE_PYYAML
SUBDIRS += tests
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py
