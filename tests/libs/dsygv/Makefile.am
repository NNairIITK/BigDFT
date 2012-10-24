EXTRA_DIST = \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

check_PROGRAMS =
my_tests =

EXTRA_PROGRAMS = test

CLEANFILES = test-dsygv.out\
	interfaces.@MODULE_EXT@ INTERFACES.@MODULE_EXT@

AM_FCFLAGS = -I$(top_builddir)/src/modules -I$(top_builddir)/src/ @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 

test_SOURCES = interfaces.f90 init.f90 wrappers.f90 test.f90
test_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

check: test-dsygv.report

test-dsygv.report: $(srcdir)/test-dsygv.out.ref test-dsygv.out
	python $(top_srcdir)/tests/fldiff.py $^ | tee $@

test-dsygv.out: test
	cp $(srcdir)/input.dat $(PWD)
	$(run_parallel)	./test > $@
