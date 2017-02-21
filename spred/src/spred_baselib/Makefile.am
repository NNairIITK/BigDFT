
AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@ @LIB_FUTILE_CFLAGS@

EXTRA_DIST = 

CLEANFILES = *.@MODULE_EXT@

YAML_DATABASE = SPRED_input_variables_definition.yaml \
	spred_bibliography.yaml 

noinst_LIBRARIES = libspredbasef.a

c_level = $(DATABASE_SOURCES)

libspredbasef_a_SOURCES = $(c_level) \
                  SPREDbase.f90 \
                  SPRED_public_keys.f90 \
                  SPREDtypes.f90 \
                  module_fingerprints.f90

include $(srcdir)/import_database.mk

#dependencies
SPREDtypes.o: SPREDbase.o SPRED_public_keys.o

module_fingerprints.o: SPREDbase.o SPREDtypes.o
