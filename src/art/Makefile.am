vpath %.h ../modules

EXTRA_DIST = README \
	bart.sh

if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

art_source = \
	defs.f90 \
	random.f90 \
	restart.f90 \
	utils.f90 \
	write_refconfig.f90


AM_FCFLAGS = -I. -I$(srcdir) -I../modules @MPI_INCLUDE@

noinst_LIBRARIES = libart.a

libart_a_SOURCES = $(art_source)

#dependencies
restart.o: defs.o

utils.o: defs.o

write_refconfig.o: defs.o

defs.o : $(mpi_include)
