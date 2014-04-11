## Process this file with automake to produce Makefile.in
if USE_MPI
mpi_include =
else
mpi_include = mpif.h
endif

EXTRA_DIST = 72_geomoptim/isokinetic.F90 \
	72_geomoptim/isotemp.F90 \
	72_geomoptim/isothermal.F90 \
	72_geomoptim/langevin.F90 \
	72_geomoptim/moldyn.F90 \
	72_geomoptim/nose.F90 \
	72_geomoptim/others.F90 \
	72_geomoptim/quenched.F90 \
	72_geomoptim/velocity_verlet.F90 \
	xallgather_mpi.finc \
	xallgatherv_mpi.finc \
	xalltoall_mpi.finc \
	xalltoallv_mpi.finc \
	xcast_mpi.finc \
	xexch_mpi.finc \
	xgather_mpi.finc \
	xgatherv_mpi.finc \
	xirecv_mpi.finc \
	xisend_mpi.finc \
	xmax_mpi.finc \
	xmin_mpi.finc \
	xmpi_ialltoall.finc \
	xrecv_mpi.finc \
	xscatterv_mpi.finc \
	xsend_mpi.finc \
	xsum_master.finc \
	xsum_mpi.finc

lib_LIBRARIES = libabinit.a

if USE_LIBXC
if BUILD_LIB_BIGDFT
LIBXC_DEFS = -DHAVE_LIBXC
endif
else
LIBXC_DEFS =
endif

# Special treatment for libpaw
# Deactivate temporarily libXC in libpaw
LIBPAW_DEFS = -DHAVE_DFT_BIGDFT -UHAVE_DFT_LIBXC

# Force HAVE_MPI and HAVE_MPI1.
# Don't use the HAVE_CONFIG_H
MPI_DEFS = -DHAVE_MPI -DHAVE_MPI1

DEFS = $(MPI_DEFS) $(LIBXC_DEFS) $(EXTRA_DEFS)

AM_FCFLAGS = @MPI_INCLUDE@ @LIBXC_INCLUDE@ -I$(srcdir)/72_geomoptim -I$(srcdir)/incs

# Hack for XLF strange DEFS option.
PPFCCOMPILE = $(FC) $(DEFS:-D%=@FCDEFS@%) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

libabinit_a_SOURCES = \
	01_interfaces_ext/interfaces_linalg.F90 \
	10_defs/defs_basis.F90 \
	10_defs/defs_datatypes.F90 \
	10_defs/defs_parameters.F90 \
	12_hide_mpi/interfaces_12_hide_mpi.F90 \
	12_hide_mpi/m_xmpi.F90 \
	12_hide_mpi/m_xomp.F90 \
	14_hidewrite/interfaces_14_hidewrite.F90 \
	14_hidewrite/wrtout.F90 \
	14_hidewrite/m_fstrings.F90 \
	14_hidewrite/m_io_tools.F90 \
	16_hideleave/interfaces_16_hideleave.F90 \
	16_hideleave/leave_new.F90 \
	16_hideleave/m_errors.F90 \
	16_hideleave/m_xieee.F90 \
	18_timing/interfaces_18_timing.F90 \
	18_timing/timab.F90 \
	18_memory/m_profiling.F90 \
	28_numeric_noabirule/interfaces_28_numeric_noabirule.F90 \
        28_numeric_noabirule/m_splines.F90 \
	28_numeric_noabirule/sort_dp.F90 \
	28_numeric_noabirule/sort_int.F90 \
	28_numeric_noabirule/uniformrandom.F90 \
	32_util/interfaces_32_util.F90 \
	32_util/canon9.F90 \
	32_util/derfc.F90 \
	32_util/mati3inv.F90 \
	32_util/matr3inv.F90 \
	32_util/matrginv.F90 \
	32_util/wrap2_pmhalf.F90 \
	41_xc_lowlevel/interfaces_41_xc_lowlevel.F90 \
	41_xc_lowlevel/drivexc.F90 \
	41_xc_lowlevel/invcb.F90 \
	41_xc_lowlevel/mkdenpos.F90 \
	41_xc_lowlevel/m_libxc_functionals.F90 \
	41_xc_lowlevel/size_dvxc.F90 \
	41_xc_lowlevel/xchcth.F90 \
	41_xc_lowlevel/xchelu.F90 \
	41_xc_lowlevel/xclb.F90 \
	41_xc_lowlevel/xcpbe.F90 \
	41_xc_lowlevel/xcpzca.F90 \
	41_xc_lowlevel/xcspol.F90 \
	41_xc_lowlevel/xctetr.F90 \
	41_xc_lowlevel/xcwign.F90 \
	41_xc_lowlevel/xcxalp.F90 \
	42_geometry/bldgrpaf.F90 \
	42_geometry/bldgrp.F90 \
	42_geometry/bonds_lgth_angles.F90 \
	42_geometry/chkdilatmx.F90 \
	42_geometry/chkgrp.F90 \
	42_geometry/chkorthsy.F90 \
	42_geometry/chkprimit.F90 \
	42_geometry/gensymshub4.F90 \
	42_geometry/gensymshub.F90 \
	42_geometry/gensymspgr.F90 \
	42_geometry/getptgroupma.F90 \
	42_geometry/getspinrot.F90 \
	42_geometry/gridgcart.F90 \
	42_geometry/holocell.F90 \
	42_geometry/m_ab6_symmetry.F90 \
	42_geometry/metric.F90 \
	42_geometry/mkrdim.F90 \
	42_geometry/operat.F90 \
	42_geometry/prtspgroup.F90 \
	42_geometry/ptgmadata.F90 \
	42_geometry/smallprim.F90 \
	42_geometry/spgdata.F90 \
	42_geometry/strainsym.F90 \
	42_geometry/strconv.F90 \
	42_geometry/stresssym.F90 \
	42_geometry/sym2cart.F90 \
	42_geometry/symanal.F90 \
	42_geometry/symatm.F90 \
	42_geometry/symaxes.F90 \
	42_geometry/symbrav.F90 \
	42_geometry/symchk.F90 \
	42_geometry/symdet.F90 \
	42_geometry/symfind.F90 \
	42_geometry/symlatt.F90 \
	42_geometry/symlist_bcc.F90 \
	42_geometry/symlist_fcc.F90 \
	42_geometry/symlist_others.F90 \
	42_geometry/symlist_prim.F90 \
	42_geometry/symmultsg.F90 \
	42_geometry/symplanes.F90 \
	42_geometry/symptgroup.F90 \
	42_geometry/symredcart.F90 \
	42_geometry/symrelrot.F90 \
	42_geometry/symsgcube.F90 \
	42_geometry/symsghexa.F90 \
	42_geometry/symsgmono.F90 \
	42_geometry/symsgortho.F90 \
	42_geometry/symsgtetra.F90 \
	42_geometry/symspgr.F90 \
	42_geometry/symzat.F90 \
	42_geometry/xredxcart.F90 \
	interfaces_42_libpaw.F90 \
	m_atompaw.F90 \
	m_gaussfit.F90 \
	m_paral_atom.F90 \
	m_paw_an.F90 \
	m_pawang.F90 \
	m_pawcprj.F90 \
	m_pawdij.F90 \
	m_pawfgrtab.F90 \
	m_paw_finegrid.F90 \
	m_paw_ij.F90 \
	m_pawio.F90 \
	m_paw_numeric.F90 \
	m_pawpsp.F90 \
	m_pawrad.F90 \
	m_pawrhoij.F90 \
	m_pawtab.F90 \
	m_pawxc.F90 \
	m_pawxmlps.F90 \
	m_sphharm.F90 \
	42_geometry/interfaces_42_geometry.F90 \
	50_abitypes_defs/defs_abitypes.F90 \
	56_mixing/aprxdr.F90 \
	56_mixing/dotprodm_v.F90 \
	56_mixing/dotprodm_vn.F90 \
	56_mixing/findminscf.F90 \
	56_mixing/interfaces_56_mixing.F90 \
	56_mixing/m_ab6_mixing.F90 \
	56_mixing/scfcge.F90 \
	56_mixing/scfeig.F90 \
	56_mixing/scfopt.F90 \
	56_mixing/sqnormm_v.F90 \
	56_recipspace/interfaces_56_recipspace.F90 \
	56_recipspace/irrzg.F90 \
	56_recipspace/m_ab6_kpoints.F90 \
	56_recipspace/getkgrid.F90 \
	56_recipspace/smpbz.F90 \
	56_recipspace/symkpt.F90 \
	56_recipspace/testkgrid.F90 \
	67_common/ewald.F90 \
        67_common/ewald2.F90 \
	67_common/fconv.F90 \
	67_common/prtxvf.F90 \
	72_geomoptim/ab6_moldyn.F90 \
	72_geomoptim/xfpack.F90

CLEANFILES = mpif.h *.@MODULE_EXT@

#dependencies
mpif.h:
	touch mpif.h &&\
	 echo "integer, parameter :: MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_COMM_NULL=1, MPI_COMM_SELF=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_LOR=1, MPI_COMPLEX=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_UNDEFINED=1, MPI_MAX_ERROR_STRING=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_DOUBLE_COMPLEX=1, MPI_STATUS_IGNORE=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_STATUS_SIZE=1, MPI_ERR_UNKNOWN=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_TAG_UB=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	 echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1,MPI_CHARACTER=1,MPI_REAL8=1" >> mpif.h

# Automatic ones...
include $(srcdir)/deps
