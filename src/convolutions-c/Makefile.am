EXTRA_DIST = \
	tool.h \
	thread_engine.h \
	topology_muscade.cfg \
	placement_muscade.cfg \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

if HAVE_LIBCONFIG
check = conv_check
else
check =
endif

EXTRA_PROGRAMS = conv_check
noinst_LIBRARIES = libbasef.a libbasec.a

libbasef_a_SOURCES = ../profiling/memory.f90 \
                     ../modules/base.F90 \
                     ../modules/cublas_fake.f90 $(mpi_source)
libbasec_a_SOURCES = thread_engine.cpp \
                     MagicFilter_sse.c

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir)
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror -O2 -msse3 -mtune=core2 -ftree-vectorize -ftree-vectorizer-verbose=9 -ffast-math
#-funroll-loops -mtune=core2 -ftree-vectorize -floop-interchange -floop-strip-mine -floop-block

conv_check_SOURCES = conv_check.f90 
conv_check_LDADD = libbasec.a libbasef.a ../libbigdft.a ../PSolver/libpoissonsolver.a
conv_check_LDFLAGS = -lpthread -lconfig++ -lstdc++ -lrt

check_conv: conv_check
	echo 1 128 128 204900 204900 10 > fort.1
	numactl --interleave=all ./conv_check
