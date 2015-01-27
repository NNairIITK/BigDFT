vpath %.h ../modules


globalf_objects = \
               module_globaltool.f90

globalf_SOURCES = $(globalf_objects)

##globalf_LDADD = $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@

##AM_FCFLAGS = -I. -I$(top_builddir)/includes
##AM_FCFLAGS = -I. I$(top_builddir)/includes
AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@\
             @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@


CLEANFILES = *.@MODULE_EXT@ 

EXTRA_DIST = \
	hunt-inc.f90

noinst_LIBRARIES = libglobalf.a

libglobalf_a_SOURCES = $(globalf_SOURCES)
#dependencies

module_globaltool.o: ../modules/BigDFT_API.o
