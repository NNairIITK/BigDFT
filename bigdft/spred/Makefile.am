SUBDIRS = src \
	tests

EXTRA_DIST = COPYING \
	README

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py
