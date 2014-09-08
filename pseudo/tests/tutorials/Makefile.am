## Process this file with automake to produce Makefile.in

SUBDIRS = oxygen.NLCC

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python $(top_builddir)/tests/report.py
