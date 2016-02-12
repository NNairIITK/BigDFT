SUBDIRS = src
if HAVE_PYYAML
SUBDIRS += tests
endif

EXTRA_DIST = COPYING \
	README

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(top_builddir)/tests/report.py
