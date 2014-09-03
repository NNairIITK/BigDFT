## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  atoms \
	  H-pseudo \
	  C-pseudo \
	  pawpatchtest \
	  tutorials

EXTRA_DIST = \
	README \
	example.O.new \
	pawpatchtest

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python $(top_srcdir)/tests/report.py
