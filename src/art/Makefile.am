EXTRA_DIST = README 

#	bart.sh

if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

art_source = bigdft_forces.f90 \
	defs.f90 \
	diis.f90 \
	find_saddle.f90 \
	initialize.f90 \
	initialize_potential.f90 \
	lanczos.f90 \
	min_converge.f90 \
	random.f90 \
	read_parameters.f90 \
	restart.f90 \
	saddle_converge.f90 \
	utils.f90 \
	write_refconfig.f90


AM_FCFLAGS = -I. -I$(srcdir) -I../modules -I../ @MPI_INCLUDE@

noinst_LIBRARIES = libart.a

libart_a_SOURCES = $(art_source)

mpif.h:
	touch mpif.h &&\
	 echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h &&\
	 echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1" >> mpif.h &&\
	 echo "integer :: MPI_STATUSES_IGNORE, MPI_LOGICAL, MPI_REAL8" >> mpif.h &&\
	 echo "integer :: MPI_MIN, MPI_MAX,MPI_CHARACTER" >> mpif.h


#dependencies
read_parameters.o: defs.o lanczos.o find_saddle.o

bigdft_forces.o: defs.o

find_saddle.o: defs.o lanczos.o

initialize.o: defs.o lanczos.o

initialize_potential.o: defs.o bigdft_forces.o

lanczos.o: defs.o random.o  bigdft_forces.o

restart.o: defs.o

min_converge.o: defs.o bigdft_forces.o

saddle_converge.o: defs.o bigdft_forces.o find_saddle.o

utils.o: defs.o

write_refconfig.o: defs.o

diis.o: defs.o bigdft_forces.o random.o

defs.o : $(mpi_include)