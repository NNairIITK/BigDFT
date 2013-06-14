!> @file
!!  Module to store all dictionary keys of the input files.
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Define all static strings to store input variables
module module_input_keys
  use dictionaries, only: dictionary

  implicit none
  
  private

  character(len = *), parameter, public :: DFT_VARIABLES = "dft"
  character(len = *), parameter, public :: HGRIDS = "hgrids"
  character(len = *), parameter, public :: RMULT = "rmult"
  character(len = *), parameter, public :: IXC = "ixc"
  character(len = *), parameter, public :: NCHARGE = "ncharge"
  character(len = *), parameter, public :: ELECFIELD = "elecfield"
  character(len = *), parameter, public :: NSPIN = "nspin", MPOL = "mpol"
  character(len = *), parameter, public :: GNRM_CV = "gnrm_cv"
  character(len = *), parameter, public :: ITERMAX = "itermax", NREPMAX = "nrepmax"
  character(len = *), parameter, public :: NCONG = "ncong", IDSX = "idsx"
  character(len = *), parameter, public :: DISPERSION = "dispersion"
  character(len = *), parameter, public :: INPUTPSIID = "inputpsiid"
  character(len = *), parameter, public :: OUTPUT_WF = "output_wf"
  character(len = *), parameter, public :: OUTPUT_DENSPOT = "output_denspot"
  character(len = *), parameter, public :: RBUF = "rbuf"
  character(len = *), parameter, public :: NCONGT = "ncongt"
  character(len = *), parameter, public :: NORBV = "norbv", NVIRT = "nvirt"
  character(len = *), parameter, public :: NPLOT = "nplot"
  character(len = *), parameter, public :: DISABLE_SYM = "disablesym"

  character(len = *), parameter, public :: KPT_VARIABLES = "kpt"
  character(len = *), parameter, public :: KPT_METHOD = "method"
  character(len = *), parameter, public :: KPTRLEN = "kptrlen"
  character(len = *), parameter, public :: NGKPT = "ngkpt"
  character(len = *), parameter, public :: SHIFTK = "shiftk"
  character(len = *), parameter, public :: KPT = "kpt"
  character(len = *), parameter, public :: WKPT = "wkpt"
  character(len = *), parameter, public :: BANDS = "bands"
  character(len = *), parameter, public :: ISEG = "iseg"
  character(len = *), parameter, public :: KPTV = "kptv"
  character(len = *), parameter, public :: NGRANULARITY = "ngranularity"
  character(len = *), parameter, public :: BAND_STRUCTURE_FILENAME = "band_structure_filename"

  character(len = *), parameter, public :: GEOPT_VARIABLES = "geopt"
  character(len = *), parameter, public :: GEOPT_METHOD = "method"
  character(len = *), parameter, public :: NCOUNT_CLUSTER_X = "ncount_cluster_x"
  character(len = *), parameter, public :: FRAC_FLUCT = "frac_fluct"
  character(len = *), parameter, public :: FORCEMAX = "forcemax"
  character(len = *), parameter, public :: RANDDIS = "randdis"
  character(len = *), parameter, public :: IONMOV = "ionmov"
  character(len = *), parameter, public :: DTION = "dtion"
  character(len = *), parameter, public :: MDITEMP = "mditemp"
  character(len = *), parameter, public :: MDFTEMP = "mdftemp"
  character(len = *), parameter, public :: NOSEINERT = "noseinert"
  character(len = *), parameter, public :: FRICTION = "friction"
  character(len = *), parameter, public :: MDWALL = "mdwall"
  character(len = *), parameter, public :: NNOS = "nnos"
  character(len = *), parameter, public :: QMASS = "qmass"
  character(len = *), parameter, public :: BMASS = "bmass"
  character(len = *), parameter, public :: VMASS = "vmass"
  character(len = *), parameter, public :: BETAX = "betax"
  character(len = *), parameter, public :: HISTORY = "history"
  character(len = *), parameter, public :: DTINIT = "dtinit"
  character(len = *), parameter, public :: DTMAX = "dtmax"

  !> Error ids for this module.
  integer, public :: INPUT_VAR_NOT_IN_LIST
  integer, public :: INPUT_VAR_NOT_IN_RANGE
  integer, public :: INPUT_VAR_ILLEGAL
  type(dictionary), pointer :: failed_exclusive

  character(len = *), parameter :: RANGE = "__range__", EXCLUSIVE = "__exclusive__"
  character(len = *), parameter :: DEFAULT = "default", COMMENT = "__comment__"
  character(len = *), parameter :: COND = "__condition__", WHEN = "__when__"
  character(len = *), parameter :: MASTER_KEY = "master_key", USER_DEFINED = "__user__"
  character(len = *), parameter :: PROF_KEY = "__profile__", ATTRS = "_attributes"

  public :: input_keys_init, input_keys_finalize
  public :: input_keys_set, input_keys_fill, input_keys_dump
  public :: input_keys_equal, input_keys_get_source
  public :: input_keys_get_profiles

  type(dictionary), pointer :: parameters

contains

  subroutine abort_excl()
    use yaml_output
    use dictionaries
    implicit none
    
    call yaml_open_map("allowed values")
    call yaml_dict_dump(failed_exclusive)
    call yaml_close_map()
    call f_err_severe()
  end subroutine abort_excl

  subroutine warn_illegal()
    use yaml_output
    use dictionaries
    implicit none
    
    call yaml_warning("variable not allowed in context")
  end subroutine warn_illegal

  subroutine input_keys_init()
    use yaml_output
    use dictionaries
    implicit none

    call dict_init(parameters)
    call set(parameters // DFT_VARIABLES, get_dft_parameters())
    call set(parameters // KPT_VARIABLES, get_kpt_parameters())
    call set(parameters // GEOPT_VARIABLES, get_geopt_parameters())

    !call yaml_dict_dump(parameters, comment_key = COMMENT)
    call f_err_define(err_name='INPUT_VAR_NOT_IN_LIST',&
         err_msg='given value not in allowed list.',&
         err_action='choose a value from the list below.',&
         err_id=INPUT_VAR_NOT_IN_LIST,callback=abort_excl)
    call f_err_define(err_name='INPUT_VAR_NOT_IN_RANGE',&
         err_msg='given value not in allowed range.',&
         err_action='adjust the given value.',&
         err_id=INPUT_VAR_NOT_IN_RANGE)
    call f_err_define(err_name='INPUT_VAR_ILLEGAL',&
         err_msg='provided variable is not allowed in this context.',&
         err_action='remove the input variable.',&
         err_id=INPUT_VAR_ILLEGAL,callback=warn_illegal)
  END SUBROUTINE input_keys_init
  
  subroutine input_keys_finalize()
    use dictionaries
    implicit none

    call dict_free(parameters)
  END SUBROUTINE input_keys_finalize
  
  function get_dft_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_dft_parameters, excl

    call dict_init(p)

    ! Settings
    call set(p // HGRIDS, dict_new( (/ &
         & COMMENT   .is. "grid spacing in the three directions (bohr)", &
         & RANGE     .is. list_new((/ .item."0.", .item."2." /)), &
         & DEFAULT   .is. list_new((/ .item."0.45", .item."0.45", .item."0.45" /)), &
         & "fast"    .is. list_new((/ .item."0.55", .item."0.55", .item."0.55" /)), &
         & "accurate".is. list_new((/ .item."0.30", .item."0.30", .item."0.30" /)) /) ))

    call set(p // RMULT, dict_new( (/ &
         & COMMENT   .is. "c(f)rmult*radii_cf(:,1(2))=coarse(fine) atom-basec radius", &
         & RANGE     .is. list_new((/ .item."0.", .item."100." /)), &
         & DEFAULT   .is. list_new((/ .item."5.", .item."8." /)) /)))

    call set(p // IXC, dict_new((/ &
         & COMMENT   .is. "exchange-correlation parameter (LDA=1,PBE=11)", &
         & DEFAULT   .is. "1", &
         & "PBE"     .is. "11", &
         & "Hybrid"  .is. "-406" /)))

    call set(p // NCHARGE, dict_new((/ &
         & COMMENT   .is. "charge of the system", &
         & RANGE     .is. list_new((/ .item."-500.", .item."500." /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // ELECFIELD, dict_new((/ &
         & COMMENT   .is. "electric fielc (Ex,Ey,Ez)", &
         & DEFAULT   .is. list_new((/ .item."0.", .item."0.", .item."0." /)) /)))

    call set(p // NSPIN, dict_new((/ &
         & COMMENT   .is. "spin polarization", &
         & EXCLUSIVE .is. dict_new((/ "1".is."no spin", "2".is."collinear", "4".is."non-collinear" /)), &
         & DEFAULT   .is. "1" /)))

    call set(p // MPOL, dict_new((/ &
         & COMMENT .is. "total magnetic moment", &
         & DEFAULT  . is. "0" /)))

    call set(p // GNRM_CV, dict_new((/ &
         & COMMENT   .is. "convergence criterion gradient", &
         & RANGE     .is. list_new((/ .item."1.e-20", .item."1." /)), &
         & DEFAULT   .is. "1e-4", &
         & "fast"    .is. "1e-3", &
         & "accurate".is. "1e-5"/)))

    call set(p // ITERMAX, dict_new((/ &
         & COMMENT   .is. "max. # of wfn. opt. steps", &
         & RANGE     .is. list_new((/ .item."0", .item."10000" /)), &
         & DEFAULT   .is. "50" /)))

    call set(p // NREPMAX, dict_new((/ &
         & COMMENT   .is. "max. # of re-diag. runs", &
         & RANGE     .is. list_new((/ .item."0", .item."1000" /)), &
         & DEFAULT   .is. "1", &
         & "accurate".is. "10" /)))

    call set(p // NCONG, dict_new((/ &
         & COMMENT   .is. "# of CG it. for preconditioning eq.", &
         & RANGE     .is. list_new((/ .item."0", .item."20" /)), &
         & DEFAULT   .is. "6" /)))

    call set(p // IDSX, dict_new((/ &
         & COMMENT     .is. "wfn. diis history", &
         & RANGE       .is. list_new((/ .item."0", .item."15" /)), &
         & DEFAULT     .is. "6", &
         & "low memory".is. "2" /)))

    call set(p // DISPERSION, dict_new((/ &
         & COMMENT   .is. "dispersion correction potential (values 1,2,3,4,5), 0=none", &
         & RANGE     .is. list_new((/ .item."0", .item."5" /)), &
         & DEFAULT   .is. "0" /)))

    call dict_init(excl)
    call set(excl // "-1000", "empty")
    call set(excl // "-2", "random")
    call set(excl // "-1", "CP2K")
    call set(excl // "0", "LCAO")
    call set(excl // "1", "wvl. in mem.")
    call set(excl // "2", "wvl. on disk")
    call set(excl // "10", "LCAO + gauss.")
    call set(excl // "11", "gauss. in mem.")
    call set(excl // "12", "gauss. on disk")
    call set(excl // "100", "Linear AO")
    call set(excl // "101", "Linear restart")
    call set(excl // "102", "Linear on disk")
    call set(p // INPUTPSIID, dict_new((/ &
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. excl /)))

    call dict_init(excl)
    call set(excl // "0", "none")
    call set(excl // "1", "plain text")
    call set(excl // "2", "Fortran binary")
    call set(excl // "3", "ETSF file format")
    call set(p // OUTPUT_WF, dict_new((/ &
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. excl /)))

    call dict_init(excl)
    call set(excl // "0", "none")
    call set(excl // "1", "density in plain text")
    call set(excl // "2", "density and potentials in plain text")
    call set(excl // "11", "density in ETSF file format")
    call set(excl // "12", "density and potentials in ETSF file format")
    call set(excl // "21", "density in cube file format")
    call set(excl // "22", "density and potentials in cube file format")
    call set(p // OUTPUT_DENSPOT, dict_new((/ &
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. excl /)))

    call set(p // RBUF, dict_new((/ &
         & COMMENT     .is. "length of the tail (AU)", &
         & RANGE       .is. list_new((/ .item."0.", .item."10." /)), &
         & DEFAULT     .is. "0.", &
         & "with tails".is. "5." /)))

    call set(p // NCONGT, dict_new((/ &
         & COMMENT   .is. "# tail CG iterations", &
         & RANGE     .is. list_new((/ .item."0", .item."50" /)), &
         & DEFAULT   .is. "30" /)))

    call set(p // NORBV, dict_new((/ &
         & COMMENT   .is. "Davidson subspace dim.", &
         & RANGE     .is. list_new((/ .item."-9999", .item."9999" /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // NVIRT, dict_new((/ &
         & COMMENT   .is. "# of opt. orbs", &
         & RANGE     .is. list_new((/ .item."0", .item."9999" /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // NPLOT, dict_new((/ &
         & COMMENT   .is. "# of plottec orbs", &
         & RANGE     .is. list_new((/ .item."0", .item."9999" /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // DISABLE_SYM, dict_new((/ &
         & COMMENT   .is. "disable the symmetry detection", &
         & DEFAULT   .is. "No" /)))

    get_dft_parameters => p
  END FUNCTION get_dft_parameters

  function get_kpt_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_kpt_parameters

    call dict_init(p)

    ! Settings
    call set(p // KPT_METHOD, dict_new( (/ &
         & COMMENT   .is. 'K-point sampling method', &
         & EXCLUSIVE .is. dict_new((/ "auto".is."based on kptrlen", &
         & "MPGrid".is."Monkhorst-Pack", "manual".is."based on raw coordinates" /)), &
         & DEFAULT   .is. "manual" /) ))

    call set(p // KPTRLEN, dict_new( (/ &
         & COMMENT   .is. 'Equivalent length of K-space resolution (Bohr)', &
         & COND      .is. dict_new((/ MASTER_KEY .is. KPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "auto" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e4" /)), &
         & DEFAULT   .is. "0." /) ))

    call set(p // NGKPT, dict_new( (/ &
         & COMMENT   .is. 'No. of Monkhorst-Pack grid points', &
         & COND      .is. dict_new((/ MASTER_KEY .is. KPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "MPGrid" /)) /)), &
         & RANGE     .is. list_new((/ .item."1", .item."10000" /)), &
         & DEFAULT   .is. list_new((/ .item."1", .item."1", .item."1" /)) /) ))

    call set(p // SHIFTK, dict_new( (/ &
         & COMMENT   .is. 'Grid shifts', &
         & COND      .is. dict_new((/ MASTER_KEY .is. KPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "MPGrid" /)) /)), &
         & DEFAULT   .is. list_new((/ .item."0.", .item."0.", .item."0." /)) /) ))

    call set(p // KPT, dict_new( (/ &
         & COMMENT   .is. 'Kpt coordinates', &
         & COND      .is. dict_new((/ MASTER_KEY .is. KPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "manual" /)) /)), &
         & DEFAULT   .is. list_new((/ .item.list_new((/ .item."0.", .item."0.", .item."0." /)) /)) /) ))

    call set(p // WKPT, dict_new( (/ &
         & COMMENT   .is. 'Kpt weights', &
         & COND      .is. dict_new((/ MASTER_KEY .is. KPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "manual" /)) /)), &
         & RANGE     .is. list_new((/ .item."0.", .item."1e4" /)), &
         & DEFAULT   .is. list_new((/ .item."1." /)) /) ))

    call set(p // BANDS, dict_new( (/ &
         & COMMENT   .is. 'For doing band structure calculation', &
         & DEFAULT   .is. "No" /) ))

    call set(p // ISEG, dict_new( (/ &
         & COMMENT   .is. 'points for each segment of the BZ path', &
         & COND      .is. dict_new((/ MASTER_KEY .is. BANDS, &
         &                            WHEN .is. list_new((/ .item. "Yes" /)) /)), &
         & DEFAULT   .is. list_new((/ .item."1" /)) /) ))

    call set(p // KPTV, dict_new( (/ &
         & COND      .is. dict_new((/ MASTER_KEY .is. BANDS, &
         &                            WHEN .is. list_new((/ .item. "Yes" /)) /)), &
         & DEFAULT   .is. list_new((/ &
         & .item. list_new((/ .item."0. ", .item."0. ", .item."0. " /)), &
         & .item. list_new((/ .item."0.5", .item."0.5", .item."0.5" /)) /)) /) ))

    call set(p // NGRANULARITY, dict_new( (/ &
         & COMMENT   .is. '# of points done for each group', &
         & COND      .is. dict_new((/ MASTER_KEY .is. BANDS, &
         &                            WHEN .is. list_new((/ .item. "Yes" /)) /)), &
         & RANGE     .is. list_new((/ .item."1", .item."1000" /)), &
         & DEFAULT   .is. "1" /) ))

    call set(p // BAND_STRUCTURE_FILENAME, dict_new( (/ &
         & COND      .is. dict_new((/ MASTER_KEY .is. BANDS, &
         &                            WHEN .is. list_new((/ .item. "Yes" /)) /)), &
         & DEFAULT   .is. "''" /) ))

    get_kpt_parameters => p
  END FUNCTION get_kpt_parameters

  function get_geopt_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_geopt_parameters

    call dict_init(p)

    ! Settings
    call set(p // GEOPT_METHOD, dict_new( (/ &
         & COMMENT   .is. 'Geometry optimisation method', &
         & EXCLUSIVE .is. dict_new((/ &
         & "SDCG"    .is. "a combination of Steepest Descent and Conjugate Gradient", &
         & "VSSD"    .is. "Variable Stepsize Steepest Descent method", &
         & "LBFGS"   .is. "limited-memory BFGS", &
         & "BFGS"    .is. "Broyden–Fletcher–Goldfarb–Shanno", &
         & "PBFGS"   .is. "Same as BFGS with an initial Hessian obtained from a force field", &
         & "AB6MD"   .is. "molecular dynamics from ABINIT", &
         & "DIIS"    .is. "direct inversion of iterative subspace", &
         & "FIRE"    .is. "Fast Inertial Relaxation Engine as described by Bitzek et al." /)), &
         & DEFAULT   .is. "none" /) ))

    call set(p // NCOUNT_CLUSTER_X, dict_new( (/ &
         & COMMENT   .is. 'Maximum number of force evaluations', &
         & RANGE     .is. list_new((/ .item."0", .item."2000" /)), &
         & PROF_KEY  .is. GEOPT_METHOD, &
         & DEFAULT   .is. "50", &
         & "none"    .is. "1"/) ))

    call set(p // FRAC_FLUCT, dict_new( (/ &
         & RANGE     .is. list_new((/ .item."0.", .item."10." /)), &
         & DEFAULT   .is. "1." /) ))

    call set(p // FORCEMAX, dict_new( (/ &
         & RANGE     .is. list_new((/ .item."0.", .item."10." /)), &
         & DEFAULT   .is. "0." /) ))

    call set(p // RANDDIS, dict_new( (/ &
         & COMMENT   .is. 'random displacement amplitude', &
         & RANGE     .is. list_new((/ .item."0.", .item."10." /)), &
         & DEFAULT   .is. "0." /) ))

    call set(p // BETAX, dict_new( (/ &
         & COMMENT   .is. 'Stepsize for the geometry optimisation', &
         & COND      .is. dict_new((/ MASTER_KEY .is. GEOPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "SDCG", &
         &                .item."VSSD", .item."LBFGS", .item."BFGS", .item."PBFGS", &
         &                .item."DIIS", .item."FIRE", .item."none" /)) /)), &
         & PROF_KEY  .is. GEOPT_METHOD, &
         & RANGE     .is. list_new((/ .item."0.", .item."100." /)), &
         & DEFAULT   .is. "4.", &
         & "DIIS"    .is. "2." /) ))

    call set(p // HISTORY, dict_new( (/ &
         & COMMENT   .is. 'history for DIIS method', &
         & COND      .is. dict_new((/ MASTER_KEY .is. GEOPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "DIIS" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."100" /)), &
         & DEFAULT   .is. "4" /) ))

    call set(p // DTINIT, dict_new( (/ &
         & COMMENT   .is. 'initial time step for the FIRE method', &
         & COND      .is. dict_new((/ MASTER_KEY .is. GEOPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "FIRE" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e4" /)), &
         & DEFAULT   .is. "0.75" /) ))

    call set(p // DTMAX, dict_new( (/ &
         & COMMENT   .is. 'maximal time step for the FIRE method', &
         & COND      .is. dict_new((/ MASTER_KEY .is. GEOPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "FIRE" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e4" /)), &
         & DEFAULT   .is. "1.5" /) ))

    call set(p // IONMOV, dict_new( (/ &
         & COMMENT   .is. 'movement ion method', &
         & COND      .is. dict_new((/ MASTER_KEY .is. GEOPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "AB6MD" /)) /)), &
         & EXCLUSIVE .is. dict_new((/ "6".is."simple velocity-Verlet molecular dynamic", &
         & "7".is."quenched molecular dynamic", "8".is."Nose-Hoover thermostat", &
         & "9".is."Langevin dynamic", "12".is."Isokinetic ensemble molecular dynamics", &
         & "13".is."Iosthermal/isenthalpic ensemble" /)), &
         & DEFAULT   .is. "6" /) ))

    call set(p // DTION, dict_new( (/ &
         & COMMENT   .is. 'time step - Atomic Units (20.670689 AU=0.5 fs)', &
         & COND      .is. dict_new((/ MASTER_KEY .is. GEOPT_METHOD, &
         &                            WHEN .is. list_new((/ .item. "AB6MD" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e3" /)), &
         & DEFAULT   .is. "20.670689" /) ))

    call set(p // MDITEMP, dict_new( (/ &
         & COMMENT   .is. 'temperature of molecular dynamics', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "6" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e9" /)), &
         & DEFAULT   .is. "300." /) ))

    call set(p // MDFTEMP, dict_new( (/ &
         & COMMENT   .is. 'final temperature of molecular dynamics', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "8", .item."9", &
         &                                                  .item."12", .item."13" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e9" /)), &
         & DEFAULT   .is. "300." /) ))

    call set(p // NOSEINERT, dict_new( (/ &
         & COMMENT   .is. 'thermostat inertia coefficient for Nose_Hoover dynamics', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "8" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e9" /)), &
         & DEFAULT   .is. "1e5" /) ))

    call set(p // FRICTION, dict_new( (/ &
         & COMMENT   .is. 'Friction coefficient for Langevin dynamics', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "9" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e5" /)), &
         & DEFAULT   .is. "1e-3" /) ))

    call set(p // MDWALL, dict_new( (/ &
         & COMMENT   .is. 'Distance in bohr where atoms can bounce for Langevin dynamics', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "9" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e5" /)), &
         & DEFAULT   .is. "1e4" /) ))

    call set(p // QMASS, dict_new( (/ &
         & COMMENT   .is. 'Mass of each thermostat (isothermal/isenthalpic ensemble)', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "13" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e9" /)), &
         & DEFAULT   .is. list_new((/ .item."0." /)) /) ))

    call set(p // BMASS, dict_new( (/ &
         & COMMENT   .is. 'Barostat masses (isothermal/isenthalpic ensemble)', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "13" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e9" /)), &
         & DEFAULT   .is. "10." /) ))

    call set(p // VMASS, dict_new( (/ &
         & COMMENT   .is. 'Barostat masses (isothermal/isenthalpic ensemble)', &
         & COND      .is. dict_new((/ MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new((/ .item. "13" /)) /)), &
         & RANGE     .is. list_new((/ .item."0", .item."1e9" /)), &
         & DEFAULT   .is. "1." /) ))

    get_geopt_parameters => p
  END FUNCTION get_geopt_parameters

  function input_keys_get_profiles(file)
    use dictionaries
    implicit none
    character(len = *), intent(in) :: file
    type(dictionary), pointer :: input_keys_get_profiles

    type(dictionary), pointer :: p
    integer :: i
    character(max_field_length), dimension(:), allocatable :: keys

    call dict_init(p)
    
    if (has_key(parameters, file)) then
       call vars(p, parameters // file)
    else
       allocate(keys(dict_size(parameters)))
       keys = dict_keys(parameters)
       do i = 1, size(keys), 1
          call vars(p // keys(i), parameters // keys(i))
       end do
       deallocate(keys)
    end if

    input_keys_get_profiles => p
  contains
    subroutine vars(dict, ref)
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict, ref

      integer :: i
      character(max_field_length), dimension(:), allocatable :: var

      allocate(var(dict_size(ref)))
      var = dict_keys(ref)
      do i = 1, size(var), 1
         call generate(dict,  var(i), ref // var(i))
      end do
      deallocate(var)
      if (dict_size(dict) == 0) then
         call set(dict, "no profile")
      end if
    end subroutine vars

    subroutine generate(dict, key, ref)
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict, ref
      character(len = *), intent(in) :: key

      integer :: i
      character(max_field_length), dimension(:), allocatable :: keys

      allocate(keys(dict_size(ref)))
      keys = dict_keys(ref)
      do i = 1, size(keys), 1
         if (trim(keys(i)) /= COMMENT .and. &
              & trim(keys(i)) /= COND .and. &
              & trim(keys(i)) /= RANGE .and. &
              & trim(keys(i)) /= PROF_KEY .and. &
              & trim(keys(i)) /= EXCLUSIVE .and. &
              & trim(keys(i)) /= DEFAULT) then
            call add(dict // key, keys(i))
         end if
      end do
      deallocate(keys)
    end subroutine generate
  END FUNCTION input_keys_get_profiles

  !> Compare two strings (case-insensitive). Blanks are relevant!
  function input_keys_equal(stra,strb)
    implicit none
    character(len=*), intent(in) :: stra,strb
    logical :: input_keys_equal
    !Local variables
    integer :: i,ica,icb,ila,ilb,ilength
    ila=len(stra)
    ilb=len(strb)
    ilength=min(ila,ilb)
    ica=ichar(stra(1:1))
    icb=ichar(strb(1:1))
    input_keys_equal=(modulo(ica-icb,32) == 0) .and. (ila==ilb)
    do i=2,ilength
       ica=ichar(stra(i:i))
       icb=ichar(strb(i:i))
       input_keys_equal=input_keys_equal .and. &
            &   (modulo(ica-icb,32) == 0)
       if (.not. input_keys_equal) exit
    end do
  END FUNCTION input_keys_equal

  function input_keys_get_source(dict, key)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: key

    character(len = max_field_length) :: input_keys_get_source

    input_keys_get_source(1:max_field_length) = DEFAULT
    if (has_key(dict, trim(key) // ATTRS)) then
       if (has_key(dict // (trim(key) // ATTRS), PROF_KEY)) then
          input_keys_get_source = dict // (trim(key) // ATTRS) // PROF_KEY
       end if
    end if
  end function input_keys_get_source

  subroutine input_keys_fill(dict, file, profile)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file
    character(len = *), intent(in), optional :: profile

    integer :: i
    type(dictionary), pointer :: ref
    character(max_field_length), dimension(:), allocatable :: keys

    ref => parameters // file
    allocate(keys(dict_size(ref)))
    keys = dict_keys(ref)

    if (present(profile)) then
       do i = 1, size(keys), 1
          call input_keys_set(dict, file, keys(i), profile)
       end do
    else
       do i = 1, size(keys), 1
          call input_keys_set(dict, file, keys(i))
       end do
    end if

    deallocate(keys)
  END SUBROUTINE input_keys_fill

  subroutine input_keys_set(dict, file, key, profile)
    use dictionaries
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file, key
    character(len = *), intent(in), optional :: profile

    integer :: i
    type(dictionary), pointer :: ref
    character(len = max_field_length) :: val, profile_
    character(len = max_field_length), dimension(:), allocatable :: keys
    double precision, dimension(2) :: rg
    logical :: found

    ref => parameters // file // key

    profile_(1:max_field_length) = " "
    if (present(profile)) then
       ! User defined profile.
       profile_(1:max_field_length) = profile
    else
       ! Hard-coded profile from key.
       if (has_key(ref, PROF_KEY)) then
          val = ref // PROF_KEY
          if (has_key(dict, val)) then
             profile_ = dict // val
          end if
       end if
    end if
    if (trim(profile_) == "") profile_(1:max_field_length) = DEFAULT
    
    if (has_key(dict, key)) then
       ! Key should be present only for some unmet conditions.
       if (f_err_raise(.not.set_(dict, ref), err_id = INPUT_VAR_ILLEGAL, &
            & err_msg = trim(key) // " is not allowed in this context.")) return

       val = dict // key
       ! There is already a value in dict.
       if (has_key(ref, val)) then
          ! The key was asking for a profile, we copy it.
          call copy(dict // key, ref // val)
          profile_ = val
       else
          ! The value is user-defined, we validate it with info from dict from.
          if (has_key(ref, RANGE)) then
             rg(1) = ref // RANGE // 0
             rg(2) = ref // RANGE // 1
             call validate(dict // key, key, rg)
          else if (has_key(ref, EXCLUSIVE)) then
             failed_exclusive => ref // EXCLUSIVE
             allocate(keys(dict_size(failed_exclusive)))
             keys = dict_keys(failed_exclusive)
             found = .false.
             do i = 1, size(keys), 1
                found = input_keys_equal(trim(val), trim(keys(i)))
                if (found) exit
             end do
             deallocate(keys)
             if (f_err_raise(.not. found, err_id = INPUT_VAR_NOT_IN_LIST, &
                  & err_msg = trim(key) // " = '" // trim(val) // "' is not allowed.")) return
             nullify(failed_exclusive)
          end if
          profile_(1:max_field_length) = USER_DEFINED
       end if
    else
       ! Key should be present only for some unmet conditions.
       if (.not.set_(dict, ref)) return

       ! There is no value in dict, we take it from ref.
       if (.not. has_key(ref, profile_)) profile_ = DEFAULT
       call copy(dict // key, ref // profile_)
    end if

    ! Copy the comment.
    if (has_key(ref, COMMENT)) &
         & call copy(dict // (trim(key) // ATTRS) // COMMENT, ref // COMMENT)
    if (trim(profile_) /= DEFAULT) &
         & call set(dict // (trim(key) // ATTRS) // PROF_KEY, profile_)

  contains

    function set_(dict, ref)
      implicit none
      type(dictionary), pointer :: dict, ref
      logical :: set_

      integer :: j
      character(max_field_length) :: mkey, val_master, val_when

      set_ = .true.
      if (has_key(ref, COND)) then
         mkey = ref // COND // MASTER_KEY
         if (.not. has_key(dict, mkey)) then
            set_ = .false.
            return
         end if
         val_master = dict // mkey
         set_ = .false.
         do j = 0, dict_len(ref // COND // WHEN) - 1, 1
            val_when = ref // COND // WHEN // j
            set_ = set_ .or. &
                 & (input_keys_equal(trim(val_master), trim(val_when)))
         end do
      end if
    end function set_

    recursive subroutine validate(dict, key, rg)
      implicit none
      type(dictionary), pointer :: dict
      character(len = *), intent(in) :: key
      double precision, dimension(2), intent(in) :: rg
      
      character(len = max_field_length) :: val
      character(max_field_length), dimension(:), allocatable :: keys
      double precision :: var

      if (associated(dict%child)) then
         if (dict_len(dict) >= 1) then
            ! List case.
            do i = 0, dict_len(dict) - 1, 1
               call validate(dict // i, key, rg)
            end do            
         else
            ! Dictionary case
            allocate(keys(dict_size(dict)))
            keys = dict_keys(dict)
            do i = 1, size(keys), 1
               call validate(dict // keys(i), key, rg)
            end do
            deallocate(keys)
         end if
      else
         var = dict
         if (f_err_raise(var < rg(1) .or. var > rg(2), err_id = INPUT_VAR_NOT_IN_RANGE, &
              & err_msg = trim(key) // " = '" // trim(val) // &
              & "' not in range.")) return
      end if
    end subroutine validate

    recursive subroutine copy(dict, ref)
      implicit none
      type(dictionary), pointer :: dict, ref

      integer :: i
      character(max_field_length), dimension(:), allocatable :: keys
      character(len = max_field_length) :: val
      
      if (associated(ref%child)) then
         if (dict_len(ref) >= 1) then
            ! List case.
            do i = 0, dict_len(ref) - 1, 1
               call copy(dict // i, ref // i)
            end do            
         else
            ! Dictionary case
            allocate(keys(dict_size(ref)))
            keys = dict_keys(ref)
            do i = 1, size(keys), 1
               call copy(dict // keys(i), ref // keys(i))
            end do
            deallocate(keys)
         end if
      else
         ! Leaf case.
         val = ref
         call set(dict, val)
      end if
    end subroutine copy
  END SUBROUTINE input_keys_set

  !> Dump a dictionary
  subroutine input_keys_dump(dict)
    use yaml_output
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict   !< Dictionary to dump

    !local variables
    integer :: i
    character(max_field_length), dimension(:), allocatable :: keys

    call yaml_comment("Input parameters", hfill = "-")
    !TEST (the first dictionary has no key)
    !if (.not. associated(dict%parent)) then
    if (associated(dict%child)) then
       if (dict_len(dict) >= 1) then
          ! List case.
          do i = 0, dict_len(dict) - 1, 1
             call dict_dump_(dict // i)
          end do
       else
          ! Dictionary case
          allocate(keys(dict_size(dict)))
          keys = dict_keys(dict)
          do i = 1, size(keys), 1
             call dict_dump_(dict // keys(i))
          end do
          deallocate(keys)
       end if
    else
       call yaml_scalar(dict%data%value)
    end if

  contains

    recursive subroutine dict_dump_(dict)
      use yaml_output
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict

      logical :: flow
      integer :: i
      type(dictionary), pointer :: parent, attr, iter
      character(max_field_length) :: descr
      character(max_field_length), dimension(:), allocatable :: keys

      descr = " "
      if (index(dict%data%key, ATTRS) > 0) return
      parent => dict%parent
      if (associated(parent)) then
         if (has_key(parent, trim(dict%data%key) // ATTRS)) then
            attr => parent // (trim(dict%data%key) // ATTRS)
            if (has_key(attr, COMMENT)) descr = attr // COMMENT
         end if
      end if

      if (associated(dict%child)) then
         if (dict_len(dict) >= 1) then
            ! List case.
            flow = (.not.associated(dict%child%child))
            if (.not.flow .and. trim(descr) /= "") then
               call yaml_open_sequence(trim(dict%data%key), advance = "no")
               call yaml_comment(trim(descr), tabbing = 50)
            else
               call yaml_open_sequence(trim(dict%data%key), flow=flow)
            end if
            do i = 0, dict_len(dict) - 1, 1
               call yaml_sequence("", advance = "no")
               call dict_dump_(dict // i)
            end do
            if (flow .and. trim(descr) /= "") then
               call yaml_close_sequence(advance = "no")
               call yaml_comment(trim(descr), tabbing = 50)
            else
               call yaml_close_sequence()
            end if
         else
            ! Dictionary case
            call yaml_open_map(trim(dict%data%key),flow=.false.)
            iter => dict_next(dict)
            allocate(keys(dict_size(dict)))
            keys = dict_keys(dict)
            do i = 1, size(keys), 1
               call dict_dump_(dict // keys(i))
            end do
            deallocate(keys)
            call yaml_close_map()
         end if
      else
         ! Leaf case.
         if (dict%data%item >= 0) then
            call yaml_sequence(trim(dict%data%value))
         else
            if (trim(descr) /= "") then
               call yaml_map(trim(dict%data%key),trim(dict%data%value),advance = "no")
               call yaml_comment(trim(descr), tabbing = 50)
            else
               call yaml_map(trim(dict%data%key),trim(dict%data%value))
            end if
         end if
      end if

    end subroutine dict_dump_
  end subroutine input_keys_dump

end module module_input_keys
