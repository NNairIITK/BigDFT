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
		wrtout.f90 \
		leave_new.f90

CLEANFILES = defs_basis.@MODULE_EXT@ \
	defs_datatypes.@MODULE_EXT@ \
	defs_xc.@MODULE_EXT@

#dependencies
