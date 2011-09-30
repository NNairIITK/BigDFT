EXTRA_DIST = \
	tool.h \
	thread_engine.h \
	topology_muscade.cfg \
	placement_muscade.cfg \
	thread_engine.cpp \
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

libbasef_a_SOURCES = $(mpi_source)
libbasef_a_LIBADD = ../modules/defs.o ../modules/base.o \
                     ../modules/cublas_fake.o
libbasec_a_SOURCES = MagicFilter_sse.c

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


CLEANFILES = module_base.@MODULE_EXT@ module_defs.@MODULE_EXT@ \
			 fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @LIBABINIT_INCLUDE@ @MPI_INCLUDE@ -I. -I$(srcdir) 
AM_CFLAGS = -Wall -Werror -O2 -msse3 -mtune=core2 -ftree-vectorize -ftree-vectorizer-verbose=9 -ffast-math
AM_LDFLAGS = -lrt -lpapi
#-funroll-loops -mtune=core2 -ftree-vectorize -floop-interchange -floop-strip-mine -floop-block

convolution_sse_SOURCES = convolution_sse.c
convolution_sse_LDFLAGS = -lrt -lpapi

convolution_sse_t_SOURCES = convolution_sse_t.c
convolution_sse_t_LDFLAGS = -lrt -lpapi


conv_check_SOURCES = conv_check.f90 
conv_check_LDADD = ../libbigdft.a ../PSolver/libpoissonsolver.a $(libs_dgemmsy) libbasec.a libbasef.a
conv_check_LDFLAGS = -lpthread -lstdc++ -lrt 


#-lpapi

check_conv: conv_check
	echo 1 128 128 204900 204900 10 > fort.1
	./conv_check
