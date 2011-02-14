SUBDIRS = 

EXTRA_DIST =

if BUILD_BINARIES_PSEUDO
binaries_sources = atom pseudo
else
binaries_sources = 
endif

AM_FCFLAGS = $(FCFLAGS) -I. @LIBXC_INCLUDE@ @MPI_INCLUDE@
AM_FFLAGS = $(FFLAGS) -I. @LIBXC_INCLUDE@ @MPI_INCLUDE@


bin_PROGRAMS = $(binaries_sources)

if USE_MPI
mpi_source =
mpi_include =
LMPI = @MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
LMPI =
endif

atom_SOURCES = atom.f atom.splines.f $(xc_level)
atom_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)


pseudo_SOURCES = pseudo.f $(low_level) $(f90_level) $(xc_level) $(mpi_source)
pseudo_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)

  
#List of sources (objects) grouped by dependencies
low_level = crtvh.f \
	penalty.f \
	amoeba.f \
	gatom.f \
	wave.f \
	wave2.f \
	wave3.f \
	detnp.f \
	resid.f \
	etot.f \
	pj2test.f \
	xpown.f \
	gamma.f \
	ppack.f \
	radgrid.f \
	zero.f \
	zbrent.f

f90_level = \
	ekin_gauss_wvlt.f90 \
	errorhandler.f90

xc_level = \
	xcfunction.f90 \
	driveXC.f90


CLEANFILES = mpif.h

# Dependencies
atom.o: xcfunction.o driveXC.o
amoeba.o: mpif.h
ekin_gauss_wvlt.cublas.o: mpif.h
ekin_gauss_wvlt.o: mpif.h
penalty.o: mpif.h
pseudo.o: mpif.h

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "      integer, parameter :: MPI_SUM=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_DOUBLE_PRECISION=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_REAL=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_INTEGER=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_STATUSES_IGNORE=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_LOGICAL=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_REAL8=1" >> mpif.h &&\
	echo "      integer, parameter :: MPI_MAX_PROCESSOR_NAME=1" >> mpif.h
	echo "      integer, parameter :: MPI_STATUS_SIZE=1" >> mpif.h
