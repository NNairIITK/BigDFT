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

atom_SOURCES = atom.f atom.splines.f $(xc_level)
atom_LDADD = \
	$(LMPI) \
	@MPI_LIBS@ \
	$(libs_cuda)


pseudo_SOURCES = pseudo.f $(low_level) $(f90_level) $(xc_level)
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


CLEANFILES = 

# Dependencies
atom.o: xcfunction.o driveXC.o
