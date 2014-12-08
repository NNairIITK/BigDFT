vpath %.h ../modules
##bin_PROGRAMS = mhgpstool


mhgps_objects = \
               module_global_variables.f90 \
               module_init.f90 \
               module_hessian.f90 \
               module_io.f90 \
               module_energyandforces.f90 \
	       module_minimizers.f90 \
               module_saddle.f90 \
               module_interpol.f90 \
               module_ls_rmsd.f90 \
               module_connect.f90 \
               module_misc.f90 \
               module_fingerprints.f90 \
               module_freezingstring.f90

##mhgpstool_objects = \
##	       module_global_variables.f90 \
##	       module_fingerprints.f90 \
##	       mhgpstool.f90

mhgps_SOURCES = $(mhgps_objects)
#mhgpstool_SOURCES = $(mhgpstool_objects)

##mhgpstool_LDADD = $(bigdft_library) $(top_builddir)/flib/src/libflib-1.a\
##                  @LIB_YAML_LIBS@ @LINALG_LIBS@
##mhgps_LDADD = $(top_builddir)/flib/src/libflib.a @LIB_YAML_LIBS@

##AM_FCFLAGS = -I. -I$(top_builddir)/includes
##AM_FCFLAGS = -I. I$(top_builddir)/includes
AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@\
             @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

CLEANFILES = *.@MODULE_EXT@ 
noinst_LIBRARIES = libmhgps.a

libmhgps_a_SOURCES = $(mhgps_SOURCES)
#dependencies

module_global_variables.o : ../modules/BigDFT_API.o

module_energyandforces.o : module_global_variables.o \
                           module_init.o \
                           ../modules/BigDFT_API.o
module_init.o : module_global_variables.o ../modules/BigDFT_API.o 
module_saddle.o : module_global_variables.o module_io.o\
                  module_energyandforces.o ../modules/module_sqn.o\
                  ../modules/BigDFT_API.o
module_minimizers.o : module_global_variables.o\
                      module_energyandforces.o\
                      ../modules/module_sqn.o\
                      ../modules/BigDFT_API.o
module_io.o : module_global_variables.o ../modules/BigDFT_API.o
module_freezingstring.o : module_energyandforces.o module_interpol.o \
                          module_global_variables.o\
                          ../modules/BigDFT_API.o
module_interpol.o : module_global_variables.o ../modules/BigDFT_API.o
module_connect.o: module_misc.o module_ls_rmsd.o\
                  module_global_variables.o module_saddle.o\
                  module_freezingstring.o module_minimizers.o\
                  ../modules/BigDFT_API.o module_fingerprints.o
module_ls_rmsd.o: module_misc.o ../modules/BigDFT_API.o
module_misc.o: ../modules/BigDFT_API.o
module_fingerprints.o: module_global_variables.o\
                      ../modules/BigDFT_API.o
module_hessian.o: ../modules/BigDFT_API.o module_energyandforces.o
##mhgpstool.o: ../modules/BigDFT_API.o module_fingerprints.o
