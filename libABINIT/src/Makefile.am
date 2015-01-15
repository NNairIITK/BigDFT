## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libabinit.a

EXTRA_DIST = \
	72_geomoptim/isokinetic.F90 \
	72_geomoptim/isotemp.F90 \
	72_geomoptim/isothermal.F90 \
	72_geomoptim/langevin.F90 \
	72_geomoptim/moldyn.F90 \
	72_geomoptim/nose.F90 \
	72_geomoptim/others.F90 \
	72_geomoptim/quenched.F90 \
	72_geomoptim/velocity_verlet.F90 \
	72_geomoptim/xfpack.F90 \
	symm/symmetry.h \
	symm/symmetry.c \
	symm/symmetry-main.c \
	libpaw/libpaw.h

## === libXC section ===
LIBXC_DEFS =
if USE_LIBXC
if BUILD_LIB_BIGDFT
LIBXC_DEFS += -DHAVE_DFT_LIBXC
endif
endif

## === libPAW section ===
LIBPAW_DEFS =
##Define libPAW src files
include $(srcdir)/libpaw/libpaw.src
libpaw_files_1 = $(LIBPAW_SRC_FILES:%.h=libpaw/%.h)
libpaw_files = $(libpaw_files_1:%.F90=libpaw/%.F90)
##Activate libPAW BigDFT sections
LIBPAW_DEFS += -DHAVE_LIBPAW_BIGDFT
##  Deactivate libXC libPAW code sections (for testing purpose)
##  LIBPAW_DEFS += -UHAVE_DFT_LIBXC
##  LIBXC_DEFS =

## === MPI section ===
mpi_include =
MPI_DEFS =
if USE_MPI
##Force MPI level 1
mpi_include += mpif.h
MPI_DEFS += -DHAVE_MPI -DHAVE_MPI1
endif

## === Fortran compiler directives ===
DEFS = $(LIBXC_DEFS) $(LIBPAW_DEFS) $(MPI_DEFS)
DEFS_1 = $(DEFS:-D%=@FCDEFS@-D%)
DEFS_FORTRAN = $(DEFS_1:-U%=@FCDEFS@-U%)

AM_FCFLAGS = @MPI_INCLUDE@ @LIBXC_INCLUDE@ -I$(srcdir)/72_geomoptim \
	        -I$(top_builddir)/includes

## Hack for XLF strange DEFS option
PPFCCOMPILE = $(FC) $(DEFS_FORTRAN) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

## === Source files ===
libabinit_a_SOURCES = \
	abi_lowlevel/abi_defs_basis.F90 \
	abi_lowlevel/abi_interfaces_lowlevel.F90 \
	abi_lowlevel/abi_interfaces_linalg.F90 \
	abi_lowlevel/m_abi_xmpi.F90 \
	abi_lowlevel/abi_wrtout.F90 \
	abi_lowlevel/abi_leave_new.F90 \
	abi_numeric/abi_interfaces_numeric.F90 \
	abi_numeric/abi_derfc.F90 \
	abi_numeric/abi_sort_dp.F90 \
	abi_numeric/abi_sort_int.F90 \
	abi_numeric/abi_mati3inv.F90 \
	abi_numeric/abi_matr3inv.F90 \
	abi_numeric/abi_matrginv.F90 \
	abi_numeric/abi_uniformrandom.F90 \
	abi_numeric/abi_wrap2_pmhalf.F90 \
	abi_xc_lowlevel/abi_interfaces_xc_lowlevel.F90 \
	abi_xc_lowlevel/abi_drivexc.F90 \
	abi_xc_lowlevel/abi_invcb.F90 \
	abi_xc_lowlevel/abi_mkdenpos.F90 \
	abi_xc_lowlevel/abi_size_dvxc.F90 \
	abi_xc_lowlevel/abi_xchcth.F90 \
	abi_xc_lowlevel/abi_xchelu.F90 \
	abi_xc_lowlevel/abi_xclb.F90 \
	abi_xc_lowlevel/abi_xcpbe.F90 \
	abi_xc_lowlevel/abi_xcpzca.F90 \
	abi_xc_lowlevel/abi_xcspol.F90 \
	abi_xc_lowlevel/abi_xctetr.F90 \
	abi_xc_lowlevel/abi_xcwign.F90 \
	abi_xc_lowlevel/abi_xcxalp.F90 \
	abi_geometry/abi_interfaces_geometry.F90 \
	abi_geometry/abi_bldgrp.F90 \
	abi_geometry/abi_bonds_lgth_angles.F90 \
	abi_geometry/abi_chkgrp.F90 \
	abi_geometry/abi_chkprimit.F90 \
	abi_geometry/abi_getptgroupma.F90 \
	abi_geometry/abi_holocell.F90 \
	abi_geometry/abi_metric.F90 \
	abi_geometry/abi_mkrdim.F90 \
	abi_geometry/abi_ptgmadata.F90 \
	abi_geometry/abi_smallprim.F90 \
	abi_geometry/abi_spgdata.F90 \
	abi_geometry/abi_strainsym.F90 \
	abi_geometry/abi_strconv.F90 \
	abi_geometry/abi_symanal.F90 \
	abi_geometry/abi_symatm.F90 \
	abi_geometry/abi_symaxes.F90 \
	abi_geometry/abi_symbrav.F90 \
	abi_geometry/abi_symchk.F90 \
	abi_geometry/abi_symdet.F90 \
	abi_geometry/abi_symfind.F90 \
	abi_geometry/abi_symlatt.F90 \
	abi_geometry/abi_symlist_bcc.F90 \
	abi_geometry/abi_symlist_fcc.F90 \
	abi_geometry/abi_symlist_others.F90 \
	abi_geometry/abi_symlist_prim.F90 \
	abi_geometry/abi_symplanes.F90 \
	abi_geometry/abi_symptgroup.F90 \
	abi_geometry/abi_symrelrot.F90 \
	abi_geometry/abi_symsgcube.F90 \
	abi_geometry/abi_symsghexa.F90 \
	abi_geometry/abi_symsgmono.F90 \
	abi_geometry/abi_symsgortho.F90 \
	abi_geometry/abi_symsgtetra.F90 \
	abi_geometry/abi_symspgr.F90 \
	abi_geometry/abi_xredxcart.F90 \
	abi_geometry/m_ab6_symmetry.F90 \
	abi_mixing/abi_interfaces_mixing.F90 \
	abi_mixing/abi_aprxdr.F90 \
	abi_mixing/abi_dotprodm_v.F90 \
	abi_mixing/abi_dotprodm_vn.F90 \
	abi_mixing/abi_findminscf.F90 \
	abi_mixing/abi_scfcge.F90 \
	abi_mixing/abi_scfeig.F90 \
	abi_mixing/abi_scfopt.F90 \
	abi_mixing/abi_sqnormm_v.F90 \
	abi_mixing/m_ab7_mixing.F90 \
	abi_recipspace/abi_interfaces_recipspace.F90 \
	abi_recipspace/abi_irrzg.F90 \
	abi_recipspace/abi_getkgrid.F90 \
	abi_recipspace/abi_smpbz.F90 \
	abi_recipspace/abi_symkpt.F90 \
	abi_recipspace/abi_testkgrid.F90 \
	abi_recipspace/m_ab6_kpoints.F90 \
	67_common/interfaces_67_common.F90 \
	67_common/ewald.F90 \
	67_common/ewald2.F90 \
	67_common/fconv.F90 \
	67_common/prtxvf.F90 \
	72_geomoptim/ab6_moldyn.F90 \
	$(libpaw_files)

CLEANFILES = mpif.h *.@MODULE_EXT@

## === Dependencies ===
## MPI
mpif.h:
	@touch mpif.h &&\
	 echo "integer, parameter :: MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_COMM_NULL=1, MPI_COMM_SELF=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_LOR=1, MPI_COMPLEX=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_UNDEFINED=1, MPI_MAX_ERROR_STRING=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_DOUBLE_COMPLEX=1, MPI_STATUS_IGNORE=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_STATUS_SIZE=1, MPI_ERR_UNKNOWN=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_TAG_UB=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1, MPI_INTEGER8=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1,MPI_CHARACTER=1,MPI_REAL8=1" >> mpif.h
	 echo "integer, parameter :: MPI_ANY_SOURCE=1, MPI_SUCCESS=1, MPI_GROUP_NULL=1" >> mpif.h
	 echo "integer, parameter :: MPI_REQUEST_NULL=1, MPI_PACKED=1, MPI_ERRORS_RETURN=1" >> mpif.h
	 echo "integer, parameter :: MPI_ERR_COMM=1, MPI_ERR_GROUP=1" >> mpif.h
	 echo "real(kind=8), external :: MPI_WTICK" >> mpif.h
## Automatic ones...
include $(srcdir)/deps
