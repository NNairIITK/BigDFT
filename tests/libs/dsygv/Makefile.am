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

if BUILD_DYNAMIC_LIBS
bigdft_library = -L$(top_builddir)/src -lbigdft-1
bigdft_library_rpath = -Wl,-rpath=$(abs_top_builddir)/src
else
bigdft_library = $(top_builddir)/src/libbigdft.a @LIBS_SHORT_DEPS@
bigdft_library_rpath =
endif

test_SOURCES = interfaces.f90 init.f90 wrappers.f90 test.f90
test_LDADD = $(bigdft_library) @LINALG_LIBS@
test_LDFLAGS = $(bigdft_library_rpath)

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
