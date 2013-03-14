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
libflib_a_SOURCES = dictionaries.f90 yaml_strings.f90 utils.c \
	$(flush_files) $(mpi_source) $(getarg_files)

# Module files.
if CAPITALIZE
modules = YAML_STRINGS.@MODULE_EXT@ \
	DICTIONARIES.@MODULE_EXT@
else
modules = yaml_strings.@MODULE_EXT@ \
	dictionaries.@MODULE_EXT@
endif
CLEANFILES = $(modules)

# Dependencies
$(modules): $(libflib_a_OBJECTS)
dictionaries.o: yaml_strings.o

# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules)
	test -e "../includes" || $(INSTALL) -d -m 755 ../includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod ../includes ; \
	done