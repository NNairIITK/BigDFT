vpath %.h ../modules

# C objects
c_level = SPRED_input_variables_fetch.c

SPRED_input_variables_fetch.o: SPRED_input_variables_definition-inc.h
#build the fortran include file with the stream of input variables into a local string
#the sed command below will transform the entire yaml document into a c string, to be passed to fortran
SPRED_input_variables_definition-inc.h: SPRED_input_variables_definition.yaml
	$(SED) -e "s/^/\"/;s/$$/\\\n\"/" $< > $@

libspredbasef_objects = $(c_level) \
                  SPREDbase.f90 \
                  SPRED_public_keys.f90 \
                  SPREDtypes.f90 \
                  module_fingerprints.f90

libspredbasef_SOURCES = $(libspredbasef_objects)


AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@\
             @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@

EXTRA_DIST = SPRED_input_variables_definition.yaml

CLEANFILES = *.@MODULE_EXT@ 

noinst_LIBRARIES = libspredbasef.a

libspredbasef_a_SOURCES = $(libspredbasef_SOURCES)
#dependencies

SPREDbase.o: SPRED_input_variables_fetch.o

SPREDtypes.o: SPREDbase.o SPRED_input_variables_fetch.o

module_fingerprints.f90: SPREDbase.o SPREDtypes.o
