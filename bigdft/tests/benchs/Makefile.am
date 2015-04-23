SUBDIRS = sparsealgebra

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python $(top_builddir)/tests/report.py
