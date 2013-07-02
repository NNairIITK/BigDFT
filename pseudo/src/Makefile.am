SUBDIRS = 

EXTRA_DIST = $(include_level)

include_level = \
	func.inc \
	intots.inc \
	recs16.inc \
	sym_16.inc

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

atom_SOURCES = atom.f90 atom.splines.f $(xc_level)
atom_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	@LIBXC_LIBS@ \
	@LINALG_LIBS@ \
	$(libs_cuda)

pseudo_SOURCES = pseudo.f90 $(low_level) $(f90_level) $(xc_level) $(cuda_level) $(mpi_source)
pseudo_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	@LIBXC_LIBS@ \
	@LINALG_LIBS@ \
	$(libs_cuda)

  
#List of sources (objects) grouped by dependencies
low_level = \
	wave.f \
	wave2.f \
	wave3.f \
	detnp.f \
	resid.f \
	pj2test.f \
	xpown.f \
	gamma.f \
	zero.f \
	zbrent.f\
	atom.splines.f

f90_level = \
	amoeba.f90 \
	crtvh.f90 \
	ekin_gauss_wvlt.f90 \
	errorhandler.f90 \
	etot.f90 \
	gatom.f90 \
	gatom_modified.f90 \
	patchgenerator.f90 \
	pawpatch.f90\
	pawpatch.dirac.f90 \
	penalty.f90 \
	ppack.f90 \
	radgrid.f90


xc_level = \
	xcfunction.f90 \
	driveXC.f90


if USE_CUDA_GPU
#cuda_level = cublas.fortran.bindings.o \
	 ekin_gauss_wvlt.cublas.o
else
cuda_level =
endif


CLEANFILES = mpif.h

# Dependencies
atom.o: xcfunction.o driveXC.o
penalty.o: $(mpi_include)
pseudo.o: $(mpi_include)
amoeba.o: $(mpi_include)
ekin_gauss_wvlt.cublas.o: $(mpi_include)
ekin_gauss_wvlt.o: $(mpi_include)

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
