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

EXTRA_PROGRAMS = conv_check convolution_sse convolution_sse_t
noinst_LIBRARIES = libbasef.a libbasec.a

libbasef_a_SOURCES = ../modules/defs.F90 ../modules/base.f90 \
                     ../modules/cublas_fake.f90 $(mpi_source)
libbasec_a_SOURCES = thread_engine.cpp \
                     MagicFilter_sse.c

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

if USE_DGEMMSY
libs_dgemmsy=../dgemmsy/libdgemmsy.a
else
libs_dgemmsy=
endif


CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@
AM_CFLAGS = -Wall -Werror -O2 -msse3 -mtune=core2 -ftree-vectorize -ftree-vectorizer-verbose=9 -ffast-math
AM_LDFLAGS = -lrt -lpapi
#-funroll-loops -mtune=core2 -ftree-vectorize -floop-interchange -floop-strip-mine -floop-block

convolution_sse_SOURCES = convolution_sse.c
convolution_sse_LDFLAGS = -lrt -lpapi

convolution_sse_t_SOURCES = convolution_sse_t.c
convolution_sse_t_LDFLAGS = -lrt -lpapi


conv_check_SOURCES = conv_check.f90 
conv_check_LDADD = libbasec.a libbasef.a ../libbigdft.a ../PSolver/libpoissonsolver.a $(libs_dgemmsy)
conv_check_LDFLAGS = -lpthread -lconfig++ -lstdc++ -lrt -lpapi

check_conv: conv_check
	echo 1 128 128 204900 204900 10 > fort.1
	./conv_check
