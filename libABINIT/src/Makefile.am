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
	72_geomoptim/velocity_verlet.F90

noinst_LIBRARIES = libabinit.a

AM_FCFLAGS = @MPI_INCLUDE@ @LIBXC_INCLUDE@ -I$(srcdir)/72_geomoptim

libabinit_a_SOURCES = \
	10_defs/defs_basis.F90 \
	10_defs/defs_datatypes.F90 \
	10_defs/defs_parameters.F90 \
	28_numeric_noabirule/sort_dp.F90 \
	28_numeric_noabirule/uniformrandom.F90 \
	32_util/canon9.F90 \
	32_util/derfc.F90 \
	32_util/mati3inv.F90 \
	32_util/matr3inv.F90 \
	32_util/matrginv.F90 \
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
	42_geometry/symlist_bcc.F90 \
	42_geometry/symlist_fcc.F90 \
	42_geometry/symlist_others.F90 \
	42_geometry/symlist_prim.F90 \
	42_geometry/symmultsg.F90 \
	42_geometry/symplanes.F90 \
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
	50_abitypes_defs/defs_abitypes.F90 \
	56_recipspace/getkgrid.F90 \
	56_recipspace/smpbz.F90 \
	56_recipspace/symkpt.F90 \
	56_recipspace/testkgrid.F90 \
	56_recipspace/ab6_symmetry_f90.f90 \
	56_xc/interfaces_56_xc.F90 \
	56_xc/drivexc.F90 \
	56_xc/invcb.F90 \
	56_xc/m_libxc_functionals.F90 \
	56_xc/size_dvxc.F90 \
	56_xc/xchcth.F90 \
	56_xc/xchelu.F90 \
	56_xc/xclb.F90 \
	56_xc/xcpbe.F90 \
	56_xc/xcpzca.F90 \
	56_xc/xcspol.F90 \
	56_xc/xctetr.F90 \
	56_xc/xcwign.F90 \
	56_xc/xcxalp.F90 \
	67_common/ewald.F90 \
	67_common/fconv.F90 \
	67_common/prtxvf.F90 \
	72_geomoptim/ab6_moldyn.F90 \
	72_geomoptim/xfpack.F90

CLEANFILES = $(mpi_include) \
	defs_abitypes.@MODULE_EXT@ DEFS_ABITYPES.@MODULE_EXT@ \
	defs_basis.@MODULE_EXT@ DEFS_BASIS.@MODULE_EXT@ \
	defs_datatypes.@MODULE_EXT@ DEFS_DATATYPES.@MODULE_EXT@ \
	defs_parameters.@MODULE_EXT@ DEFS_PARAMETERS.@MODULE_EXT@ \
	defs_xc.@MODULE_EXT@ DEFS_XC.@MODULE_EXT@ \
	libxc_functionals.@MODULE_EXT@ LIBXC_FUNCTIONALS.@MODULE_EXT@ \
	ab6_moldyn.@MODULE_EXT@ AB6_MOLDYN.@MODULE_EXT@ \
	ab6_symmetry.@MODULE_EXT@ AB6_SYMMETRY.@MODULE_EXT@

#dependencies
10_defs/defs_parameters.F90: defs_basis.o
10_defs/defs_datatypes.F90: defs_parameters.o \
	defs_basis.o

32_util/derfc.F90: defs_basis.o
32_util/canon9.F90: defs_basis.o
32_util/mati3inv.F90: defs_basis.o
32_util/matr3inv.F90: defs_basis.o
32_util/matginv.F90: defs_basis.o
32_util/matrginv.F90: defs_basis.o

42_geometry/bldgrpaf.F90: defs_basis.o
42_geometry/bldgrp.F90: defs_basis.o
42_geometry/bonds_lgth_angles.F90: defs_basis.o
42_geometry/chkdilatmx.F90: defs_basis.o
42_geometry/chkgrp.F90: defs_basis.o
42_geometry/chkorthsy.F90: defs_basis.o
42_geometry/chkprimit.F90: defs_basis.o
42_geometry/gensymshub4.F90: defs_basis.o
42_geometry/gensymshub.F90: defs_basis.o
42_geometry/gensymspgr.F90: defs_basis.o
42_geometry/getptgroupma.F90: defs_basis.o
42_geometry/getspinrot.F90: defs_basis.o
42_geometry/gridgcart.F90: defs_basis.o defs_datatypes.o
42_geometry/holocell.F90: defs_basis.o
42_geometry/metric.F90: defs_basis.o
42_geometry/mkrdim.F90: defs_basis.o
42_geometry/operat.F90: defs_basis.o
42_geometry/prtspgroup.F90: defs_basis.o
42_geometry/ptgmadata.F90: defs_basis.o
42_geometry/smallprim.F90: defs_basis.o
42_geometry/spgdata.F90: defs_basis.o
42_geometry/strainsym.F90: defs_basis.o
42_geometry/strconv.F90: defs_basis.o
42_geometry/stresssym.F90: defs_basis.o
42_geometry/sym2cart.F90: defs_basis.o
42_geometry/symanal.F90: defs_basis.o
42_geometry/symatm.F90: defs_basis.o
42_geometry/symaxes.F90: defs_basis.o
42_geometry/symbrav.F90: defs_basis.o
42_geometry/symchk.F90: defs_basis.o
42_geometry/symdet.F90: defs_basis.o
42_geometry/symfind.F90: defs_basis.o
42_geometry/symlist_bcc.F90: defs_basis.o
42_geometry/symlist_fcc.F90: defs_basis.o
42_geometry/symlist_others.F90: defs_basis.o
42_geometry/symlist_prim.F90: defs_basis.o
42_geometry/symmultsg.F90: defs_basis.o
42_geometry/symplanes.F90: defs_basis.o
42_geometry/symredcart.F90: defs_basis.o
42_geometry/symrelrot.F90: defs_basis.o
42_geometry/symsgcube.F90: defs_basis.o
42_geometry/symsghexa.F90: defs_basis.o
42_geometry/symsgmono.F90: defs_basis.o
42_geometry/symsgortho.F90: defs_basis.o
42_geometry/symsgtetra.F90: defs_basis.o
42_geometry/symspgr.F90: defs_basis.o
42_geometry/symzat.F90: defs_basis.o
42_geometry/xredxcart.F90: defs_basis.o

50_abitypes_defs/defs_abitypes.F90: defs_basis.o defs_datatypes.o

56_recipspace/getkgrid.F90: defs_basis.o
56_recipspace/smpbz.F90: defs_basis.o
56_recipspace/symkpt.F90: defs_basis.o
56_recipspace/testkgrid.F90: defs_basis.o
56_recipspace/ab6_symmetry_f90.f90: defs_basis.o

56_xc/interfaces_56_xc.F90: defs_basis.o
56_xc/m_libxc_functionals.F90: defs_basis.o
56_xc/drivexc.F90: defs_basis.o \
	defs_datatypes.o \
	m_libxc_functionals.o \
	interfaces_56_xc.o
56_xc/invcb.F90: defs_basis.o
56_xc/size_dvxc.F90: defs_basis.o $(mpi_include)
56_xc/xchcth.F90: defs_basis.o
56_xc/xchelu.F90: defs_basis.o
56_xc/xclb.F90: defs_basis.o
56_xc/xcpbe.F90: defs_basis.o
56_xc/xcpzca.F90: defs_basis.o 
56_xc/xcspol.F90: defs_basis.o
56_xc/xctetr.F90: defs_basis.o
56_xc/xcwign.F90: defs_basis.o
56_xc/xcxalp.F90: defs_basis.o

67_common/ewald.F90: defs_basis.o
67_common/fconv.F90: defs_basis.o
67_common/prtxvf.F90: defs_basis.o

72_geomoptim/ab6_moldyn.F90: defs_basis.o defs_datatypes.o

mpif.h:
	touch mpif.h &&\
	 echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h &&\
	 echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1" >> mpif.h &&\
	 echo "integer :: MPI_STATUSES_IGNORE, MPI_LOGICAL" >> mpif.h &&\
	 echo "integer :: MPI_MIN, MPI_MAX,MPI_CHARACTER,MPI_REAL8" >> mpif.h
