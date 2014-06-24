vpath %.h ../modules
###bin_PROGRAMS = mhgps

mhgps_objects = \
               module_global_variables.f90 \
               module_init.f90 \
               module_energyandforces.f90 \
               module_saddle.f90

mhgps_SOURCES = $(mhgps_objects)

##mhgps_LDADD = $(bigdft_library) $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@ @LINALG_LIBS@
##mhgps_LDADD = $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@

##AM_FCFLAGS = -I. -I$(top_builddir)/includes
##AM_FCFLAGS = -I. I$(top_builddir)/includes
AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

CLEANFILES = *.@MODULE_EXT@ 
noinst_LIBRARIES = libmhgps.a

libmhgps_a_SOURCES = $(mhgps_SOURCES)

#dependencies
mhgps.o : $(mhgps_objects:.f90=.o) ../modules/BigDFT_API.o

module_energyandforces.o : module_global_variables.o module_init.o ../modules/BigDFT_API.o
module_init.o : module_global_variables.o ../modules/BigDFT_API.o 
module_saddle.o : module_global_variables.o module_energyandforces.o ../modules/BigDFT_API.o

##
##ions_mod.o matrix_mod.o options_mod.o: kind_mod.o
##
##charge_mod.o chgcar_mod.o bader_mod.o dipole_mod.o voronoi_mod.o : matrix_mod.o
##
##io_mod.o cube_mod.o : ions_mod.o charge_mod.o
##
##dipole_mod.o chgcar_mod.o : charge_mod.o
##
##bader_mod.o voronoi_mod.o dipole_mod.o: io_mod.o
##
##io_mod.o : chgcar_mod.o cube_mod.o
##
##dipole_mod.o : bader_mod.o
