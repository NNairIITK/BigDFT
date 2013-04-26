# To be installed:
noinst_LIBRARIES = libflib.a
module_DATA = $(modules)

# Source code.
if HAVE_FC_FLUSH
flush_files = flush.f90
else
flush_files = flush_fake.f90
endif
if USE_MPI
mpi_source =
else
mpi_source = MPIfake.f90
endif
if HAVE_FC_GET_COMMAND_ARGUMENT
getarg_files =
else
getarg_files = get_command_argument_fake.f90
endif
libflib_a_SOURCES = dictionaries.f90 yaml_strings.f90 utils.c utilsadd.c \
	yaml_output.f90 dynamic_memory.f90 razero.f90 getadd.f90 error_handling.f90 \
	$(flush_files) $(mpi_source) $(getarg_files)

DYN_INCLUDES =  deallocate-inc.f90 allocate-inc.f90 f_malloc-inc.f90
YAML_INCLUDES = yaml_map-inc.f90 yaml_map-arr-inc.f90

EXTRA_DIST = $(DYN_INCLUDES) $(YAML_INCLUDES)

# Module files.
if CAPITALIZE
modules = YAML_STRINGS.@MODULE_EXT@ \
	DICTIONARIES.@MODULE_EXT@ \
	YAML_OUTPUT.@MODULE_EXT@ \
	ERROR_HANDLING.@MODULE_EXT@ \
	METADATA_INTERFACES.@MODULE_EXT@ \
	DYNAMIC_MEMORY.@MODULE_EXT@
else
modules = yaml_strings.@MODULE_EXT@ \
	dictionaries.@MODULE_EXT@ \
	yaml_output.@MODULE_EXT@ \
	metadata_interfaces.@MODULE_EXT@ \
	error_handling.@MODULE_EXT@ \
	dynamic_memory.@MODULE_EXT@
endif

CLEANFILES = $(modules)

AM_FCFLAGS = -I$(srcdir) @LIBABINIT_INCLUDE@

# Dependencies
$(modules): $(libflib_a_OBJECTS)
dictionaries.o: yaml_strings.o
yaml_output.o : dictionaries.o $(YAML_INCLUDES)
error_handling.o: dictionaries.o yaml_output.o
dynamic_memory.o: yaml_strings.o yaml_output.o $(DYN_INCLUDES)


# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules)
	test -e "../includes" || $(INSTALL) -d -m 755 ../includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod ../includes ; \
	done