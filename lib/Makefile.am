## Process this file with automake to produce Makefile.in

vpath %.f90 ../src

modulelib_LIBRARIES = \
	libbigdft.a

modulemod_DATA = \
	lib_bigdft.mod

FCCOMPILE = $(FC) $(INCLUDES) $(AM_FCFLAGS) $(FCFLAGS)

INCLUDES = \
	-I../ \
	-I../src

clean-local: 
	-rm -f mpif.h lib_bigdft.mod

if USE_OPTI_CONVOLUT
convolution_source = convolut_optim_omp.f90
else
convolution_source = convolut_simple.f90
endif

if USE_MPI
mpi_source =
mpi_include =
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
endif

libbigdft_a_SOURCES = $(mpi_source)\
	module_bigdft.f90

#dependencies
module_bigdft.o: $(mpi_include) \
	$(convolution_source) \
	common_routines.f90 \
	fft3d.f90 \
	gauss_to_daub.f90 \
	gequad.f90 \
	pfft3d.f90 \
	precond_fb.f90 \
	scaling_function.f90

mpif.h:
	touch mpif.h
	echo "real*8 :: MPI_double_precision" >> mpif.h
	echo "real*8 :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h

