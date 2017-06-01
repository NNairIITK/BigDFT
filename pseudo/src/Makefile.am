SUBDIRS = 

EXTRA_DIST = $(include_level) 

include_level = \
	func.inc \
	intots.inc \
	recs16.inc \
	sym_16.inc

AM_FCFLAGS = $(FCFLAGS) -I. @LIB_XC_CFLAGS@ @MPI_INCLUDE@
AM_FFLAGS = $(FFLAGS) -I. @LIB_XC_CFLAGS@ @MPI_INCLUDE@

bin_PROGRAMS = atom pseudo

if USE_MPI
mpi_source =
mpi_include =
LMPI = @MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
LMPI =
endif

atom_SOURCES = atom.f90 atom.splines.f90 $(xc_level)
atom_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	@LIB_XC_LIBS@ \
	@LINALG_LIBS@ \
	$(libs_cuda)

pseudo_SOURCES = pseudo.f90 $(f90_level) $(xc_level) $(cuda_level) $(mpi_source)
pseudo_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	@LIB_XC_LIBS@ \
	@LINALG_LIBS@ \
	$(libs_cuda)


#List of sources (objects) grouped by dependencies

f90_level = \
	amoeba.f90 \
	atom.splines.f90 \
	crtvh.f90 \
	detnp.f90 \
	ekin_gauss_wvlt.f90 \
	errorhandler.f90 \
	etot.f90 \
	gamma.f90 \
	gatom.f90 \
	gatom_modified.f90 \
	patchgenerator.f90 \
	pawpatch.f90\
	pawpatch.dirac.f90 \
	penalty.f90 \
	pj2test.f90 \
	ppack.f90 \
	radgrid.f90 \
	resid.f90 \
	wave.f90 \
	wave2.f90 \
	wave3.f90 \
	xpown.f90 \
	zbrent.f90 


xc_level = \
	xcfunction.f90 \
	driveXC.f90


#if USE_CUDA_GPU
#cuda_level = cublas.fortran.bindings.o \
	 ekin_gauss_wvlt.cublas.o
#else
cuda_level =
#endif


CLEANFILES = mpif.h *.@MODULE_EXT@

# Dependencies
atom.o: xcfunction.o driveXC.o
driveXC.o: xcfunction.o
penalty.o: $(mpi_include) pseudo.o
pseudo.o: xcfunction.o $(mpi_include)
amoeba.o: $(mpi_include) pseudo.o
ekin_gauss_wvlt.cublas.o: $(mpi_include)
ekin_gauss_wvlt.o: $(mpi_include)
gatom.o: pseudo.o
ppack.o: pseudo.o

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
