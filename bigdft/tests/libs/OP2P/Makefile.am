EXTRA_DIST = \
	diis_check.f90 \
	kpts_check.f90

check_PROGRAMS = diis_check

check: $(check_PROGRAMS)

report:
		@if test $(MAKELEVEL) = 0 ; then python $(top_builddir)/tests/report.py ; fi

EXTRA_PROGRAMS = kpts_check

CLEANFILES = kpts_check.out

AM_FCFLAGS = -I$(top_builddir)/includes @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
bigdft_library = $(top_builddir)/src/libbigdft-1.so.@BIGDFT_MINOR_VERSION@
else
bigdft_library = $(top_builddir)/src/libbigdft-1.a @LIB_BIGDFT_DEPS@
endif

kpts_check_SOURCES = kpts_check.f90
kpts_check_LDADD = $(bigdft_library)

diis_check_SOURCES = diis_check.f90
diis_check_LDADD = $(bigdft_library)
