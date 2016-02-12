## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libabinit.a

module_DATA = $(modules)

EXTRA_DIST = \
	abi_geomoptim/abi_isokinetic.F90.inc \
	abi_geomoptim/abi_isotemp.F90.inc \
	abi_geomoptim/abi_isothermal.F90.inc \
	abi_geomoptim/abi_langevin.F90.inc \
	abi_geomoptim/abi_moldyn.F90.inc \
	abi_geomoptim/abi_nose.F90.inc \
	abi_geomoptim/abi_others.F90.inc \
	abi_geomoptim/abi_quenched.F90.inc \
	abi_geomoptim/abi_velocity_verlet.F90.inc \
	abi_geomoptim/abi_xfpack.F90.inc \
	symm/symmetry.h \
	symm/symmetry.c \
	symm/symmetry-main.c


## === libPAW section ===
LIBPAW_DEFS =
##Define libPAW src files
include $(srcdir)/libpaw/libpaw.src
libpaw_files_1 = $(LIBPAW_SRC_FILES:%.h=libpaw/%.h)
libpaw_files = $(libpaw_files_1:%.F90=libpaw/%.F90)
##Activate libPAW BigDFT sections
LIBPAW_DEFS += @LIBPAW_HOST_DEFS@ -DHAVE_DFT_LIBXC
##  Deactivate libXC libPAW code sections (for testing purpose)
##  LIBPAW_DEFS += -UHAVE_DFT_LIBXC
##  LIBXC_DEFS =

## === MPI section ===
mpi_include =
MPI_DEFS =
if USE_MPI
##Force MPI level 1
MPI_DEFS += -DHAVE_MPI -DHAVE_MPI1
else
mpi_include += mpif.h
endif

## === Fortran compiler directives ===
DEFS = $(LIBPAW_DEFS) $(MPI_DEFS)
DEFS_1 = $(DEFS:-D%=@FCDEFS@-D%)
DEFS_FORTRAN = $(DEFS_1:-U%=@FCDEFS@-U%)

AM_FCFLAGS = @MPI_INCLUDE@ -I. -I$(srcdir)/abi_geomoptim @LIB_FLIB_CFLAGS@ @LIB_XC_CFLAGS@

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
	abi_recipspace/abi_interfaces_recipspace.F90 \
	abi_recipspace/abi_irrzg.F90 \
	abi_recipspace/abi_getkgrid.F90 \
	abi_recipspace/abi_smpbz.F90 \
	abi_recipspace/abi_symkpt.F90 \
	abi_recipspace/abi_testkgrid.F90 \
	abi_recipspace/m_ab6_kpoints.F90 \
	abi_common/abi_interfaces_common.F90 \
	abi_common/abi_ewald.F90 \
	abi_common/abi_ewald2.F90 \
	abi_common/abi_fconv.F90 \
	abi_common/abi_prtxvf.F90 \
	abi_geomoptim/m_ab6_moldyn.F90 \
	add_libpaw/abi_pawinit.F90 \
	add_libpaw/abi_initrhoij.F90 \
	add_libpaw/abi_interfaces_add_libpaw.F90 \
	add_libpaw/m_abi_distribfft.F90 \
	add_libpaw/abi_wvl_nhatgrid.F90 \
	add_libpaw/abi_pawnhatfr.F90 \
	add_libpaw/abi_pawmknhat.F90 \
	add_libpaw/abi_pawaccrhoij.F90 \
	add_libpaw/abi_pawxenergy.F90 \
	add_libpaw/abi_pawdensities.F90 \
	add_libpaw/abi_pawdenpot.F90 \
	$(libpaw_files)

CLEANFILES = mpif.h *.@MODULE_EXT@

$(modules): $(libabinit_a_OBJECTS)

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
libpaw.h: $(srcdir)/libpaw/libpaw.h
	$(SED) -n '1h;1!H;$${;g;s%/\*[^*]*\*/%%g;p;}' $< > $@

# Modules
if CAPITALIZE
modules = \
	ABI_INTERFACES_LOWLEVEL.@MODULE_EXT@ \
	ABI_DEFS_BASIS.@MODULE_EXT@ \
	M_LIBPAW_DEFS.@MODULE_EXT@ \
	M_ABI_XMPI.@MODULE_EXT@ \
	ABI_INTERFACES_LINALG.@MODULE_EXT@ \
	ABI_INTERFACES_NUMERIC.@MODULE_EXT@ \
	ABI_INTERFACES_XC_LOWLEVEL.@MODULE_EXT@ \
	ABI_INTERFACES_GEOMETRY.@MODULE_EXT@ \
	M_AB6_SYMMETRY.@MODULE_EXT@ \
	ABI_INTERFACES_MIXING.@MODULE_EXT@ \
	ABI_INTERFACES_RECIPSPACE.@MODULE_EXT@ \
	ABI_INTERFACES_COMMON.@MODULE_EXT@ \
	ABI_INTERFACES_ADD_LIBPAW.@MODULE_EXT@ \
	M_LIBPAW_MPI.@MODULE_EXT@ \
	M_AB6_MOLDYN.@MODULE_EXT@ \
	M_AB6_KPOINTS.@MODULE_EXT@ \
	M_LIBPAW_TOOLS.@MODULE_EXT@ \
	M_PAW_NUMERIC.@MODULE_EXT@ \
	M_PAWRAD.@MODULE_EXT@ \
	M_PAW_SPHHARM.@MODULE_EXT@ \
	M_LIBPAW_LIBXC_FUNCS.@MODULE_EXT@ \
	M_LIBPAW_LIBXC.@MODULE_EXT@ \
	M_PAWTAB.@MODULE_EXT@ \
	M_PARAL_ATOM.@MODULE_EXT@ \
	M_PAW_IO.@MODULE_EXT@ \
	M_PAWXMLPS.@MODULE_EXT@ \
	M_PAW_GAUSSFIT.@MODULE_EXT@ \
	M_PAWANG.@MODULE_EXT@ \
	M_PAWFGRTAB.@MODULE_EXT@ \
	M_PAWXC.@MODULE_EXT@ \
	M_PAW_ATOM.@MODULE_EXT@ \
	M_PAW_IJ.@MODULE_EXT@ \
	M_PAW_AN.@MODULE_EXT@ \
	M_PAWRHOIJ.@MODULE_EXT@ \
	M_PAW_FINEGRID.@MODULE_EXT@ \
	M_PAWCPRJ.@MODULE_EXT@ \
	M_PAWPSP.@MODULE_EXT@ \
	M_PAWDIJ.@MODULE_EXT@
else
modules = \
	abi_interfaces_lowlevel.@MODULE_EXT@ \
	abi_defs_basis.@MODULE_EXT@ \
	m_libpaw_defs.@MODULE_EXT@ \
	m_abi_xmpi.@MODULE_EXT@ \
	abi_interfaces_linalg.@MODULE_EXT@ \
	abi_interfaces_numeric.@MODULE_EXT@ \
	abi_interfaces_xc_lowlevel.@MODULE_EXT@ \
	abi_interfaces_geometry.@MODULE_EXT@ \
	m_ab6_symmetry.@MODULE_EXT@ \
	abi_interfaces_mixing.@MODULE_EXT@ \
	abi_interfaces_recipspace.@MODULE_EXT@ \
	abi_interfaces_common.@MODULE_EXT@ \
	abi_interfaces_add_libpaw.@MODULE_EXT@ \
	m_libpaw_mpi.@MODULE_EXT@ \
	m_ab6_moldyn.@MODULE_EXT@ \
	m_ab6_kpoints.@MODULE_EXT@ \
	m_libpaw_tools.@MODULE_EXT@ \
	m_paw_numeric.@MODULE_EXT@ \
	m_pawrad.@MODULE_EXT@ \
	m_paw_sphharm.@MODULE_EXT@ \
	m_libpaw_libxc_funcs.@MODULE_EXT@ \
	m_libpaw_libxc.@MODULE_EXT@ \
	m_pawtab.@MODULE_EXT@ \
	m_paral_atom.@MODULE_EXT@ \
	m_paw_io.@MODULE_EXT@ \
	m_pawxmlps.@MODULE_EXT@ \
	m_paw_gaussfit.@MODULE_EXT@ \
	m_pawang.@MODULE_EXT@ \
	m_pawfgrtab.@MODULE_EXT@ \
	m_pawxc.@MODULE_EXT@ \
	m_paw_atom.@MODULE_EXT@ \
	m_paw_ij.@MODULE_EXT@ \
	m_paw_an.@MODULE_EXT@ \
	m_pawrhoij.@MODULE_EXT@ \
	m_paw_finegrid.@MODULE_EXT@ \
	m_pawcprj.@MODULE_EXT@ \
	m_pawpsp.@MODULE_EXT@ \
	m_pawdij.@MODULE_EXT@
endif

# Automatic ones...
include $(srcdir)/deps
