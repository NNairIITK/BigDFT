## Process this file with automake to produce Makefile.in

#LIBS =

EXTRA_DIST = \
	Trash \
	Deprecated \
	daub_8.inc \
	gsint_x.inc \
	gsints.inc \
	intots.inc \
	lazy_100.inc \
	lazy_14.inc \
	lazy_16.inc \
	lazy_20.inc \
	lazy_24.inc \
	lazy_30.inc \
	lazy_40.inc \
	lazy_50.inc \
	lazy_60.inc \
	lazy_8.inc \
	recd10.inc \
	recd12.inc \
	recd14.inc \
	recd16.inc \
	recd6.inc \
	recd8.inc \
	recs16.inc \
	sym_16.inc


bin_PROGRAMS = cluster

SUBDIRS = \
	ABINIT-common

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

INCLUDES = -I. \
	@MPI_INCLUDE@ \
	-I./ABINIT-common

AM_FCFLAGS =  -I. \
	@MPI_INCLUDE@ \
	-I./ABINIT-common

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
	gequad.f90



#AM_CFLAGS = -Wall -Wno-unused -W -Werror -pedantic-errors -ansi $(cflagWin)

cluster_LDFLAGS = $(LMPI) $(lMPI)

cluster_LDADD = @LAPACK_LIB@ @BLAS_LIB@ \
	ABINIT-common/libABINIT.a

#dependencies
MPIfake.f90: mpif.h

mpif.h:
	touch mpif.h
	echo "real*8 :: MPI_double_precision" >> mpif.h
	echo "real*8 :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h

cluster.o: parameters.h
