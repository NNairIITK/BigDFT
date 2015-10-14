## Process this file with automake to produce Makefile.in

SHORT_CHECK=
LONG_CHECK=xc.out xc.out.ref xc.out.report

include $(srcdir)/../check.mk


SUBDIRS = 

EXTRA_DIST = xc.out.ref

CLEANFILES = xc.out xc.tmp.* xc.report

check_PROGRAMS = test_xc

AM_FCFLAGS = @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

test_xc_SOURCES = test_xc.f90 
test_xc_LDADD = -lbigdft-1

xc.out: test_xc
	-$(run_parallel) ./test_xc > xc.out

xc.out.report: xc.out xc.out.ref
	@cut -b -70 xc.out > xc.tmp.out
	@cut -b -70 $(srcdir)/xc.out.ref > xc.tmp.ref
	python $(top_srcdir)/tests/fldiff.py --discrepancy=9e-6 xc.tmp.out xc.tmp.ref | tee xc.report
	@rm -f xc.tmp.*
