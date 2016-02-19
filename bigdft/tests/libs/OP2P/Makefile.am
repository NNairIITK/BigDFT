EXTRA_DIST = \
  op2p_check.ref.yaml \
  diis_check.f90 \
	example.f90 \
	kpts_check.f90 \
	mpi_check.f90 \
	mpi_layer.f90 \
	op2p_check.f90

check_PROGRAMS = op2p_check diis_check

check: op2p_check.report.yaml $(check_PROGRAMS)

report:
		@if test $(MAKELEVEL) = 0 ; then python $(top_builddir)/tests/report.py ; fi

EXTRA_PROGRAMS = kpts_check mpi_check

CLEANFILES = op2p_check.out.yaml kpts_check.out\
	mpi_layer.@MODULE_EXT@ MPI_LAYER.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@ OVERLAP_POINT_TO_POINT.@MODULE_EXT@

AM_FCFLAGS = -I$(top_builddir)/includes @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

op2p_check_SOURCES = op2p_check.f90
op2p_check_LDADD = -lbigdft-1

kpts_check_SOURCES = kpts_check.f90
kpts_check_LDADD = -lbigdft-1

mpi_check_SOURCES = mpi_check.f90
mpi_check_LDADD = -lbigdft-1

diis_check_SOURCES = diis_check.f90
diis_check_LDADD = -lbigdft-1

op2p_check.o: $(top_builddir)/src/libbigdft-1.a

op2p_check.out.yaml: op2p_check
	$(run_parallel)	./op2p_check -o [10,10] -n 100 > $@

op2p_check.report.yaml: op2p_check.out.yaml
	python $(top_builddir)/tests/fldiff_yaml.py -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/op2p_check.ref.yaml -d $< -o $@
