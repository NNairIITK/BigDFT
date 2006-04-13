## Process this file with automake to produce Makefile.in

#INCLUDES =
#LIBS =

EXTRA_DIST = \
	Trash \
	Deprecated

bin_PROGRAMS = cluster

#SUBDIRS = \
#	dumpModules

if USE_MPI
mpi_source =
lMPI = -lmpi -lfmpi
LMPI = -L@MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
lMPI =
LMPI =
endif

if USE_OPTI_CONVOLUT
convolution_source = convolut_optim.f90
else
convolution_source = convolut_simple.f90
endif

#AM_FCFLAGS = @DEFS@

cluster_SOURCES = $(mpi_source) \
	Build_Kernel.f90 \
	chargedens.f90 \
	cluster.f90 \
	$(convolution_source) \
	fft3d.f90 \
	gauss_to_daub.f90 \
	gridlib.f90 \
	inputwaves.f90 \
	ionicpots.f90 \
	localization.f90 \
	localpots.f90 \
	orthonorm.f90 \
	pconv.f90 \
	pfft3d.f90 \
	plotting.f90 \
	precond.f90 \
	projectors.f90 \
	PSolver_Kernel.f90 \
	scaling_function.f90 \
	waveioput.f90 \
	wavelib.f90 \
	gequad.f



#AM_CFLAGS = -Wall -Wno-unused -W -Werror -pedantic-errors -ansi $(cflagWin)

cluster_LDFLAGS = $(LMPI)

cluster_LDADD = $(lMPI)

#dependencies
MPIfake.f90: mpif.h

mpif.h:
	touch mpif.h
