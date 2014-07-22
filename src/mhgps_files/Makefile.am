vpath %.h ../modules
###bin_PROGRAMS = mhgps

mhgps_objects = \
               module_global_variables.f90 \
               module_init.f90 \
               module_io.f90 \
               module_energyandforces.f90 \
			   module_minimizers.f90 \
               module_saddle.f90 \
               module_interpol.f90 \
               module_freezingstring.f90

mhgps_SOURCES = $(mhgps_objects)

##mhgps_LDADD = $(bigdft_library) $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@ @LINALG_LIBS@
##mhgps_LDADD = $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@

##AM_FCFLAGS = -I. -I$(top_builddir)/includes
##AM_FCFLAGS = -I. I$(top_builddir)/includes
AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

CLEANFILES = *.@MODULE_EXT@ 
noinst_LIBRARIES = libmhgps.a

libmhgps_a_SOURCES = $(mhgps_SOURCES)

%.o:     %.c
	        nab -c  -g $<

#dependencies
mhgps.o : $(mhgps_objects:.f90=.o) ../modules/BigDFT_API.o wrap.o

module_energyandforces.o : module_global_variables.o module_init.o ../modules/BigDFT_API.o wrap.o
module_init.o : module_global_variables.o ../modules/BigDFT_API.o 
module_saddle.o : module_global_variables.o module_io.o module_energyandforces.o ../modules/module_sbfgs.o ../modules/BigDFT_API.o
module_minimizers.o : module_global_variables.o module_energyandforces.o ../modules/module_sbfgs.o ../modules/BigDFT_API.o
module_io.o : module_global_variables.o ../modules/BigDFT_API.o
module_freezingstring.o : module_interpol.o module_global_variables.o ../modules/BigDFT_API.o
module_interpol.o : module_global_variables.o ../modules/BigDFT_API.o
