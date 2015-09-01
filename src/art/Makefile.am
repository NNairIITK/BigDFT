vpath %.h ../modules

EXTRA_DIST = README \
	bart.sh

#if USE_MPI
#mpi_include =
#else
#mpi_include = mpif.h
#endif

art_source = \
	defs.f90 \
	random.f90 \
	restart.f90 \
	utils.f90 \
        bigdft_forces.f90 \
        art_step.f90 \
        calcforce.f90\
        calcfo_sw.f90\
        neighbour.f90\
        initialize_potential.f90 \
        art_lanczos.f90 \
        find_saddle.f90 \
        min_converge.f90 \
        saddle_converge.f90 \
        initialize.f90 \
        bayes.f90 \
        read_parameters.f90

AM_FCFLAGS = -I. -I$(srcdir) -I$(top_builddir)/includes @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

CLEANFILES = *.@MODULE_EXT@

noinst_LIBRARIES = libart.a

libart_a_SOURCES = $(art_source)

#dependencies
restart.o: defs.o

utils.o: defs.o

#defs.o : $(mpi_include)

art.o: defs.o random.o art_lanczos.o

bigdft_forces.o: ../modules/BigDFT_API.o defs.o

bayes.o: defs.o random.o

calcforce.o: defs.o bigdft_forces.o

calcfo_sw.o: defs.o random.o

neighbour.o: calcfo_sw.o

find_saddle.o: defs.o random.o art_lanczos.o

initialize_potential.o: defs.o bigdft_forces.o

art_lanczos.o: defs.o random.o bigdft_forces.o

min_converge.o: defs.o bigdft_forces.o art_lanczos.o bayes.o find_saddle.o

saddle_converge.o: defs.o bigdft_forces.o find_saddle.o art_lanczos.o  art_step.o

art_step.o: defs.o bigdft_forces.o art_lanczos.o  find_saddle.o

initialize.o: defs.o bigdft_forces.o art_lanczos.o find_saddle.o

read_parameters.o: defs.o art_lanczos.o find_saddle.o

write_refconfig.o: defs.o

