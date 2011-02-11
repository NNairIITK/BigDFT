SUBDIRS = 

EXTRA_DIST =

if BUILD_BINARIES
binaries_sources = atom pseudo
else
binaries_sources = 
endif

AM_FCFLAGS = -I. -Imodules @MPI_INCLUDE@ @LIBXC_INCLUDE@

bin_PROGRAMS = $(binaries_sources)

if USE_MPI
mpi_source =
LMPI = @MPI_LDFLAGS@
else
mpi_source = MPIfake.f90
LMPI =
endif

atom_SOURCES = atom.f
atom_LDFLAGS =
atom_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)


pseudo_SOURCES = pseudo.f
pseudo_LDFLAGS =
pseudo_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)


#List of sources (objects) grouped by dependencies

OBJf77	=crtvh.o \
	penalty.o \
	amoeba.o \
	gatom.o \
	wave.o \
	wave2.o \
	wave3.o \
	detnp.o \
	resid.o \
	etot.o \
	pj2test.o \
	xpown.o \
	gamma.o \
	ppack.o \
	radgrid.o \
	zero.o \
	zbrent.o

OBJmain	= pseudo.o

OBJf90	= ekin_gauss_wvlt.o \
	errorhandler.o

OBJxc	= xcfunction.o \
	driveXC.o

OBJat	= atom.o \
	atom.splines.o 


# dependencies


$(OBJf77): %.o: ../src/%.f
		@echo f77 object...
		$(FC)   $(FFlags)  -c $< -o $@

$(OBJf90): %.o: ../src/%.f90
		@echo f90 objects...
		$(F90C) $(F90Flags)  -c $< -o $@

$(OBJmain): %.o: ../src/%.f
		@echo pseudo.o, needs libxc ...
		$(FC)  $(libxcInc)  $(FFlags) $(libxcLib) -c $< -o $@


$(OBJxc): %.o: ../src/%.f90
		@echo objects that need libxc ...
		$(F90C)  $(libxcInc)  $(F90Flags) $(libxcLib) -c $< -o $@


$(OBJat): %.o: ../src/%.f
		@echo objects for atom ...
		$(FC)  $(libxcInc)  $(FFlags) $(libxcLib) -c $< -o $@


clean:
		rm -f *.o

