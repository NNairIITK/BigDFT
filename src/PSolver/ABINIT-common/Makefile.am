## Process this file with automake to produce Makefile.in

INCLUDES = @MPI_INCLUDE@
#LIBS =

EXTRA_DIST =

noinst_LIBRARIES = libABINIT.a

SUBDIRS =

AM_FCFLAGS = @MPI_INCLUDE@ @LIBXC_INCLUDE@

libABINIT_a_SOURCES = \
	defs_basis.F90 \
	defs_parameters.F90 \
        defs_datatypes.F90 \
	defs_abitypes.F90 \
        defs_xc.F90 \
	interfaces_14_hidewrite.F90 \
	interfaces_16_hideleave.F90 \
	interfaces_56_xc.F90 \
	m_libxc_functionals.F90 \
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

CLEANFILES = \
	defs_abitypes.@MODULE_EXT@ DEFS_ABITYPES.@MODULE_EXT@ \
	defs_basis.@MODULE_EXT@ DEFS_BASIS.@MODULE_EXT@ \
	defs_datatypes.@MODULE_EXT@ DEFS_DATATYPES.@MODULE_EXT@ \
	defs_parameters.@MODULE_EXT@ DEFS_PARAMETERS.@MODULE_EXT@ \
	defs_xc.@MODULE_EXT@ DEFS_XC.@MODULE_EXT@ \
	libxc_functionals.@MODULE_EXT@ LIBXC_FUNCTIONALS.@MODULE_EXT@

#dependencies

interfaces_14_hidewrite.o interfaces_16_hideleave.o interfaces_56_xc.o m_libxc_functionals.o defs_parameters.o : defs_basis.o

m_libxc_functionals.o : interfaces_16_hideleave.o 

defs_datatypes.o: defs_datatypes.F90 \
	defs_parameters.o \
	defs_basis.o

defs_abitypes.o: defs_abitypes.F90 \
	defs_datatypes.o

defs_xc.o: defs_xc.F90 \
	defs_basis.o \
	defs_datatypes.o

drivexc.o: drivexc.F90 \
	defs_basis.o \
	defs_datatypes.o \
	defs_xc.o \
	m_libxc_functionals.o \
	interfaces_56_xc.o \
	interfaces_16_hideleave.o \
	interfaces_14_hidewrite.o

invcb.o: invcb.F90 \
	interfaces_16_hideleave.o \
	defs_basis.o

size_dvxc.o: size_dvxc.F90 \
	defs_basis.o

xchcth.o: xchcth.F90 \
	defs_basis.o \
	interfaces_14_hidewrite.o

xchelu.o: xchelu.F90 \
	interfaces_16_hideleave.o \
	defs_basis.o

xclb.o: xclb.F90 \
	defs_basis.o

xcpbe.o: xcpbe.F90 \
	defs_basis.o

xcpzca.o: xcpzca.F90 \
	interfaces_16_hideleave.o \
	defs_basis.o 

xcspol.o: xcspol.F90 \
	defs_basis.o \
	defs_xc.o

xctetr.o: xctetr.F90 \
	defs_basis.o \
	interfaces_14_hidewrite.o

xcwign.o: xcwign.F90 \
	defs_basis.o \
	interfaces_14_hidewrite.o

xcxalp.o: xcxalp.F90 \
	defs_basis.o \
	interfaces_16_hideleave.o \
	interfaces_14_hidewrite.o
