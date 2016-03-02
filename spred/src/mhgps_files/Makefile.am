mhgps_objects = \
               module_mhgps_state.f90 \
               module_userinput.f90 \
               module_hessian.f90 \
               module_io.f90 \
               module_energyandforces.f90 \
	       module_minimizers.f90 \
               module_saddle.f90 \
               module_interpol.f90 \
               module_ls_rmsd.f90 \
               module_connect_object.f90 \
               module_connect.f90 \
               module_misc.f90 \
               module_freezingstring.f90 \
               module_mhgpstool.f90

mhgps_SOURCES = $(mhgps_objects)

##mhgps_LDADD = $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@

##AM_FCFLAGS = -I. -I$(top_builddir)/includes
##AM_FCFLAGS = -I. I$(top_builddir)/includes
AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/src/spred_baselib @MPI_INCLUDE@ @LIB_BIGDFT_CFLAGS@

CLEANFILES = *.@MODULE_EXT@ 
noinst_LIBRARIES = libmhgps.a
EXTRA_DIST = \
	analyze.sh \
	mhgps_version_number-inc.f90

libmhgps_a_SOURCES = $(mhgps_SOURCES)
#dependencies

module_mhgps_state.o : mhgps_version_number-inc.f90

module_energyandforces.o : module_mhgps_state.o
module_userinput.o : module_mhgps_state.o
module_saddle.o : module_userinput.o module_mhgps_state.o\
                  module_energyandforces.o
module_minimizers.o : module_userinput.o module_mhgps_state.o\
                      module_energyandforces.o
module_io.o : module_userinput.o module_mhgps_state.o module_connect_object.o
module_freezingstring.o : module_energyandforces.o module_interpol.o \
                          module_userinput.o module_mhgps_state.o
module_interpol.o : module_userinput.o module_mhgps_state.o
module_connect.o: module_connect_object.o module_io.o module_userinput.o module_misc.o module_ls_rmsd.o\
                  module_mhgps_state.o module_saddle.o\
                  module_freezingstring.o module_minimizers.o

module_ls_rmsd.o: module_misc.o

module_hessian.o: module_userinput.o module_mhgps_state.o module_energyandforces.o
mhgpstool.o: module_mhgpstool.o

module_mhgpstool.o: module_io.o
