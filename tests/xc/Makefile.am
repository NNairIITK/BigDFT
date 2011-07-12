## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = xc.out.ref

CLEANFILES = xc.out xc.tmp.*

check_PROGRAMS = test_xc

test_xc_SOURCES = test_xc.f90 
test_xc_LDADD =  -lbigdft
test_xc_LDFLAGS = -L$(top_builddir)/src

AM_FCFLAGS = -I$(top_builddir)/src/modules

xc.out: test_xc
	./test_xc > xc.out

check: xc.out xc.out.ref
	@cut -b -70 xc.out > xc.tmp.out
	@cut -b -70 $(srcdir)/xc.out.ref > xc.tmp.ref
	python $(top_srcdir)/tests/fldiff.py xc.tmp.out xc.tmp.ref | tee fldiff.report
	@rm -f xc.tmp.*
