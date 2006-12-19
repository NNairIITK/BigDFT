## Process this file with automake to produce Makefile.in

SUBDIRS = \
	ABINIT-common \
	PSolver

EXTRA_DIST = \
	daub_8.inc \
	gsint_x.inc \
	gsints.inc \
	intots.inc \
	recd10.inc \
	recd12.inc \
	recd14.inc \
	recd16.inc \
	recd6.inc \
	recd8.inc \
	recs16.inc \
	sym_16.inc

AM_FCFLAGS = @MPI_INCLUDE@


bin_PROGRAMS = cluster test_forces
lib_LIBRARIES = libbigdft.a
module_DATA = libbigdft.@MODULE_EXT@

if USE_MPI
mpi_source =
mpi_include =
LMPI = @MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
lMPI =
endif

if USE_OPTI_CONVOLUT
convolution_source = convolut_optim_omp.f90 \
	standconv_optim.f90
else
convolution_source = convolut_simple.f90 \
	standconv_simple.f90
endif


libbigdft_a_SOURCES = \
	$(mpi_source) \
	$(convolution_source) \
	cluster_all.f90 \
	standconv_test.f90 \
	gauss_to_daub.f90 \
	precond_sg.f90 \
	forces.f90 \
	timing.f90

cluster_SOURCES = BigDFT.f90
cluster_LDFLAGS = $(LMPI)
cluster_LDADD = libbigdft.a \
	PSolver/libPoisson.a \
	ABINIT-common/libABINIT.a \
	@MPI_LIBS@

test_forces_SOURCES = test_forces.f90
test_forces_LDFLAGS = $(LMPI)
test_forces_LDADD = libbigdft.a \
	PSolver/libPoisson.a \
	ABINIT-common/libABINIT.a \
	@MPI_LIBS@

CLEANFILES = $(mpi_include) libbigdft.@MODULE_EXT@

#dependencies
cluster_all.o: $(mpi_include)

mpif.h:
	touch mpif.h
	echo "real*8 :: MPI_double_precision" >> mpif.h
	echo "real*8 :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h

cluster.o: parameters.h
test_forces.o: parameters.h
