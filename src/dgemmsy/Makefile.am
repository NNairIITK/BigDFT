EXTRA_DIST = \
	dgemmsy.h \
	gemm_block_c.h \
	patterns.h \
	utils.h \
	visitors.h


SOURCES_CONV_CHECK =
SOURCES_SHARED =


#EXTRA_PROGRAMS = conv_check bench_conv bench_conv_quick bench_zgemm stream_check

noinst_LIBRARIES = libdgemmsy.a

#libbasef_a_SOURCES = ../profiling/memory.f90 \
#                     ../modules/base.F90 \
#                     ../modules/cublas_fake.f90 $(mpi_source)

libbasec_a_SOURCES = dgemmsy.c \
		     gemm_block_c.c \
		     patterns.c \
		     utils.c \
		     visitors.c

libdgemmsy_a_SOURCES = $(libbasec_a_SOURCES)

if USE_MPI
mpi_source =
else
mpi_source = ../MPIfake.f90
endif

CLEANFILES = module_base.@MODULE_EXT@ fort.1 fort.98 *linkinfo

#AM_FCFLAGS = -I../modules @MPI_INCLUDE@ -I. -I$(srcdir)
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror -msse3

