SUBDIRS = src

if HAVE_PYYAML
SUBDIRS += tests
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py
