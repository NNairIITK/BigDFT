GPUBLAS_INTERFACE = cublas_fake.f90

modules_source = base.f90 \
	types.f90 \
	interfaces.f90 \
	$(GPUBLAS_INTERFACE)

AM_FCFLAGS = -I. -I$(srcdir)

noinst_LIBRARIES = libmodules.a

libmodules_a_SOURCES = $(modules_source)

CLEANFILES = MODULE_BASE.@MODULE_EXT@ module_base.@MODULE_EXT@ \
	MODULE_TYPES.@MODULE_EXT@ module_types.@MODULE_EXT@ \
	MODULE_INTERFACES.@MODULE_EXT@ module_interfaces.@MODULE_EXT@
#dependencies
interfaces.o: base.o types.o
types.o: base.o

