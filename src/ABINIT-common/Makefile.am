## Process this file with automake to produce Makefile.in

INCLUDES = @MPI_INCLUDE@
#LIBS =

EXTRA_DIST =

noinst_LIBRARIES = libABINIT.a

SUBDIRS =

AM_FCFLAGS = @MPI_INCLUDE@

libABINIT_a_SOURCES = \
	defs_basis.F90 \
        defs_datatypes.F90 \
        defs_xc.F90 \
        drivexc.F90 \
        xctetr.F90 \
        xcwign.F90 \
        xcxalp.F90 \
        xchelu.F90 \
        xcpzca.F90 \
        xcspol.F90 \
        xchcth.F90 \
        xclb.F90 \
        xcpbe.F90 \
        invcb.F90 \
	size_dvxc.F90

CLEANFILES = defs_basis.@MODULE_EXT@ \
	defs_datatypes.@MODULE_EXT@ \
	defs_xc.@MODULE_EXT@

#dependencies
defs_basis.@MODULE_EXT@: defs_basis.o
defs_xc.@MODULE_EXT@: defs_xc.o
defs_datatypes.@MODULE_EXT@: defs_datatypes.o

defs_datatypes.o: defs_datatypes.F90 \
	defs_basis.@MODULE_EXT@

defs_xc.o: defs_xc.F90 \
	defs_basis.@MODULE_EXT@ \
	defs_datatypes.@MODULE_EXT@

drivexc.o: drivexc.F90 \
	defs_basis.@MODULE_EXT@ \
	defs_xc.@MODULE_EXT@

invcb.o: invcb.F90 \
	defs_basis.@MODULE_EXT@

size_dvxc.o: size_dvxc.F90 \
	defs_basis.@MODULE_EXT@

xchcth.o: xchcth.F90 \
	defs_basis.@MODULE_EXT@

xchelu.o: xchelu.F90 \
	defs_basis.@MODULE_EXT@

xclb.o: xclb.F90 \
	defs_basis.@MODULE_EXT@

xcpbe.o: xcpbe.F90 \
	defs_basis.@MODULE_EXT@

xcpzca.o: xcpzca.F90 \
	defs_basis.@MODULE_EXT@

xcspol.o: xcspol.F90 \
	defs_basis.@MODULE_EXT@

xctetr.o: xctetr.F90 \
	defs_basis.@MODULE_EXT@

xcwign.o: xcwign.F90 \
	defs_basis.@MODULE_EXT@

xcxalp.o: xcxalp.F90 \
	defs_basis.@MODULE_EXT@
