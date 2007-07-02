## Process this file with automake to produce Makefile.in

if BUILD_LIB_XC
sub_dirs_xc = ABINIT-common
else
sub_dirs_xc = 
endif
if BUILD_LIB_PSOLVER
sub_dirs_ps = PSolver
else
sub_dirs_ps = 
endif
if BUILD_LIB_BIGDFT
lib_bigdft = libbigdft.a
mod_bigdft = libbigdft.@MODULE_EXT@
else
lib_bigdft =
mod_bigdft =
endif
if BUILD_BINARIES
binaries_sources = cluster test_forces
else
binaries_sources = 
endif


SUBDIRS = $(sub_dirs_xc) $(sub_dirs_ps)

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

AM_FCFLAGS = -I. -I$(srcdir) @MPI_INCLUDE@ -IPSolver

bin_PROGRAMS = $(binaries_sources)

lib_LIBRARIES = $(lib_bigdft)
module_DATA = $(mod_bigdft)

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
convolution_source = convolut_simple.f90 \
	conv_per_optim.f90 \
	combined_grow.f90\
	combined_shrink.f90\
	logrid.f90\
	to_cf.f90

else
convolution_source = convolut_simple.f90 \
	conv_per_optim.f90 \
	combined_grow.f90\
	combined_shrink.f90\
	logrid.f90\
	to_cf.f90
endif


libbigdft_a_SOURCES = \
	$(mpi_source) \
	$(convolution_source) \
        cluster_all.f90 \
        convolut_simple.f90 \
        conv_per_optim.f90 \
        eleconf.f90  \
        gauss_to_daub.f90 \
        precond.f90 \
        prec_diag.f90 \
        timing.f90 \
        forces.f90

cluster_SOURCES = BigDFT.f90
cluster_LDFLAGS = $(LMPI)
cluster_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	ABINIT-common/libABINIT.a \
	@MPI_LIBS@

test_forces_SOURCES = test_forces.f90
test_forces_LDFLAGS = $(LMPI)
test_forces_LDADD = libbigdft.a \
	PSolver/libpoissonsolver.a \
	ABINIT-common/libABINIT.a \
	@MPI_LIBS@

CLEANFILES = $(mpi_include) libbigdft.@MODULE_EXT@

#dependencies
cluster_all.o: $(mpi_include)
timing.o: $(mpi_include)

mpif.h:
	touch mpif.h
	echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h
	echo "integer :: MPI_DOUBLE_PRECISION" >> mpif.h
	echo "integer :: MPI_MIN, MPI_MAX" >> mpif.h

test_forces.o: parameters.h \
	libbigdft.@MODULE_EXT@

libbigdft.@MODULE_EXT@: cluster_all.o

BigDFT.o: BigDFT.f90 \
	parameters.h \
	libbigdft.@MODULE_EXT@

forces.o: forces.f90 \
	libbigdft.@MODULE_EXT@
