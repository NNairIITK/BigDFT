EXTRA_DIST = $(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

check_PROGRAMS = conv_check 

noinst_LIBRARIES = libbasef.a libbasec.a

libbasef_a_SOURCES = ../profiling/memory.f90 ../modules/base.F90 ../modules/cublas_fake.f90 $(mpi_source)
libbasec_a_SOURCES = OpenCL_wrappers.c

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir)

conv_check_SOURCES = conv_check.f90 
conv_check_LDFLAGS = -lOpenCL
conv_check_LDADD = libbasec.a libbasef.a ../libbigdft.a \
	../PSolver/libpoissonsolver.a 

check:
	echo 1 64 64 2048 2048 1 > fort.1
	$(top_builddir)/src/OpenCL/conv_check
