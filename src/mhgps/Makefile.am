bin_PROGRAMS = mhgps

mhgps_objects = \
               globvars.f90 \
               energyandforces.f90 \
               findsad.f90

mhgps_SOURCES = mhgps.f90 $(mhgps_objects)

#Use the flib library (use only yaml_output)
mhgps_LDADD = $(bigdft_library) $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@ @LINALG_LIBS@
##mhgps_LDADD = $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@

AM_FCFLAGS = -I. -I$(top_builddir)/includes

CLEANFILES = *.@MODULE_EXT@ 

#dependencies
mhgps.o : $(mhgps_objects:.f90=.o)
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
