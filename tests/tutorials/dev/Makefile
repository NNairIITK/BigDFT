FC = mpif90.openmpi
FCFLAGS = -fbounds-check

INCLUDES = -I../../../tmp/src -I../../../tmp/src/PSolver

LIBS_ABINIT = -L../../../tmp/libABINIT/src -labinit
LIBS_LIBXC = -L../../../tmp/libxc-1.1.0/src/.libs -lxc
LIBS_LAPACK = -llapack -lblas
LIBS_ARCHIVE = -larchive
LIBS_ETSF_IO = -L../../../../usr/lib -letsf_io -lnetcdf -lnetcdff

LDFLAGS = -L../../../tmp/src -L../../../tmp/src/PSolver
LIBS = -lbigdft -lpoissonsolver \
	$(LIBS_ABINIT) $(LIBS_LIBXC) $(LIBS_LAPACK) $(LIBS_ARCHIVE) $(LIBS_ETSF_IO)

test: toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(INCLUDES) $< $(LDFLAGS) $(LIBS)