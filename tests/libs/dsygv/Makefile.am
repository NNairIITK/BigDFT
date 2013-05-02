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
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

test_SOURCES = interfaces.f90 init.f90 wrappers.f90 test.f90
test_LDADD = -lbigdft-1 @LINALG_LIBS@

check: test-dsygv.report

test-dsygv.report: $(srcdir)/test-dsygv.out.ref test-dsygv.out
	python $(top_srcdir)/tests/fldiff.py $^ | tee $@

test-dsygv.out: test
	cp $(srcdir)/input1.dat $(PWD)/input.dat
	$(run_parallel)	./test > $@
	cp $(srcdir)/input2.dat $(PWD)/input.dat
	$(run_parallel)	./test >> $@
	cp $(srcdir)/input3.dat $(PWD)/input.dat
	$(run_parallel)	./test >> $@
	cp $(srcdir)/input4.dat $(PWD)/input.dat
	$(run_parallel)	./test >> $@
