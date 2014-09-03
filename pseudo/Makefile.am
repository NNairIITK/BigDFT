ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = doc

SUBDIRS = src tests

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python $(top_srcdir)/tests/report.py
