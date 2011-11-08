## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = xc.out.ref

CLEANFILES = xc.out xc.tmp.* xc.report

check_PROGRAMS = test_xc

test_xc_SOURCES = test_xc.f90 
test_xc_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

AM_FCFLAGS = -I$(top_builddir)/src/modules @LIBABINIT_INCLUDE@ 

xc.out: test_xc
	$(run_parallel) ./test_xc > xc.out

check: xc.out xc.out.ref
	@cut -b -70 xc.out > xc.tmp.out
	@cut -b -70 $(srcdir)/xc.out.ref > xc.tmp.ref
	python $(top_srcdir)/tests/fldiff.py --discrepancy=9e-6 xc.tmp.out xc.tmp.ref | tee xc.report
	@rm -f xc.tmp.*
