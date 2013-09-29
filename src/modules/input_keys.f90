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

  character(len = *), parameter, public :: MIX_VARIABLES = "mix"
  character(len = *), parameter, public :: ISCF = "iscf"
  character(len = *), parameter, public :: ITRPMAX = "itrpmax"
  character(len = *), parameter, public :: RPNRM_CV = "rpnrm_cv"
  character(len = *), parameter, public :: NORBSEMPTY = "norbsempty"
  character(len = *), parameter, public :: TEL = "Tel"
  character(len = *), parameter, public :: OCCOPT = "occopt"
  character(len = *), parameter, public :: ALPHAMIX = "alphamix"
  character(len = *), parameter, public :: ALPHADIIS = "alphadiis"

  character(len = *), parameter, public :: SIC_VARIABLES = "sic"
  character(len = *), parameter, public :: SIC_APPROACH = "sic_approach"
  character(len = *), parameter, public :: SIC_ALPHA = "sic_alpha"
  character(len = *), parameter, public :: SIC_FREF = "sic_fref"

  character(len = *), parameter, public :: TDDFT_VARIABLES = "tddft"
  character(len = *), parameter, public :: TDDFT_APPROACH = "tddft_approach"

  character(len = *), parameter, public :: PERF_VARIABLES = "perf"
  character(len = *), parameter, public :: DEBUG = "debug"
  character(len = *), parameter, public :: FFTCACHE = "fftcache"
  character(len = *), parameter, public :: ACCEL = "accel"
  character(len = *), parameter, public :: OCL_PLATFORM = "ocl_platform"
  character(len = *), parameter, public :: OCL_DEVICES = "ocl_devices"
  character(len = *), parameter, public :: BLAS = "blas"
  character(len = *), parameter, public :: PROJRAD = "projrad"
  character(len = *), parameter, public :: EXCTXPAR = "exctxpar"
  character(len = *), parameter, public :: IG_DIAG = "ig_diag"
  character(len = *), parameter, public :: IG_NORBP = "ig_norbp"
  character(len = *), parameter, public :: IG_BLOCKS = "ig_blocks"
  character(len = *), parameter, public :: IG_TOL = "ig_tol"
  character(len = *), parameter, public :: METHORTHO = "methortho"
  character(len = *), parameter, public :: RHO_COMMUN = "rho_commun"
  character(len = *), parameter, public :: PSOLVER_GROUPSIZE = "psolver_groupsize"
  character(len = *), parameter, public :: PSOLVER_ACCEL = "psolver_accel"
  character(len = *), parameter, public :: UNBLOCK_COMMS = "unblock_comms"
  character(len = *), parameter, public :: LINEAR = "linear"
  character(len = *), parameter, public :: TOLSYM = "tolsym"
  character(len = *), parameter, public :: SIGNALING = "signaling"
  character(len = *), parameter, public :: SIGNALTIMEOUT = "signaltimeout"
  character(len = *), parameter, public :: DOMAIN = "domain"
  character(len = *), parameter, public :: INGUESS_GEOPT = "inguess_geopt"
  character(len = *), parameter, public :: STORE_INDEX = "store_index"
  character(len = *), parameter, public :: VERBOSITY = "verbosity"
  character(len = *), parameter, public :: OUTDIR = "outdir"
  character(len = *), parameter, public :: PSP_ONFLY = "psp_onfly"
  character(len = *), parameter, public :: PDSYEV_BLOCKSIZE = "pdsyev_blocksize"
  character(len = *), parameter, public :: PDGEMM_BLOCKSIZE = "pdgemm_blocksize"
  character(len = *), parameter, public :: MAXPROC_PDSYEV = "maxproc_pdsyev"
  character(len = *), parameter, public :: MAXPROC_PDGEMM = "maxproc_pdgemm"
  character(len = *), parameter, public :: EF_INTERPOL_DET = "ef_interpol_det"
  character(len = *), parameter, public :: EF_INTERPOL_CHARGEDIFF = "ef_interpol_chargediff"
  character(len = *), parameter, public :: MIXING_AFTER_INPUTGUESS = "mixing_after_inputguess"
  character(len = *), parameter, public :: ITERATIVE_ORTHOGONALIZATION = "iterative_orthogonalization"
  character(len = *), parameter, public :: CHECK_SUMRHO = "check_sumrho"

  !> Error ids for this module.
  integer, public :: INPUT_VAR_NOT_IN_LIST = 0
  integer, public :: INPUT_VAR_NOT_IN_RANGE = 0
  integer, public :: INPUT_VAR_ILLEGAL = 0
  type(dictionary), pointer :: failed_exclusive

  character(len = *), parameter :: RANGE = "__range__", EXCLUSIVE = "__exclusive__"
  character(len = *), parameter :: DEFAULT = "default", COMMENT = "__comment__"
  character(len = *), parameter :: COND = "__condition__", WHEN = "__when__"
  character(len = *), parameter :: MASTER_KEY = "__master_key__", USER_DEFINED = "__user__"
  character(len = *), parameter :: PROF_KEY = "__profile__", ATTRS = "_attributes"

  public :: input_keys_init, input_keys_finalize
  public :: input_keys_set, input_keys_fill, input_keys_fill_all, input_keys_dump
  public :: input_keys_equal, input_keys_get_source, input_keys_dump_def
  public :: input_keys_get_profiles

  type(dictionary), pointer :: parameters

contains

  subroutine abort_excl()
    use yaml_output
    use dictionaries
    implicit none

    call f_dump_last_error()
    call yaml_open_map("allowed values")
    call yaml_dict_dump(failed_exclusive)
    call yaml_close_map()
    call f_err_severe()
  end subroutine abort_excl

  subroutine warn_illegal()
    use yaml_output
    use dictionaries
    implicit none
    
    integer :: ierr
    character(len = max_field_length) :: val

    ierr = f_get_last_error(val)
    call yaml_warning(trim(val))
  end subroutine warn_illegal

  subroutine input_keys_init()
    use yaml_output
    use dictionaries
    use dynamic_memory
    use yaml_parse
    implicit none
    !local variables
    integer :: params_size
    integer(kind = 8) :: cbuf_add !< address of c buffer
    character, dimension(:), allocatable :: params

!!$    !alternative filling of parameters from hard-coded source file
!!$    call getstaticinputdef(cbuf_add,params_size)
!!$    !allocate array
!!$    params=f_malloc(params_size,id='params')
!!$    !fill it and parse dictionary
!!$    call copycbuffer(params,cbuf_add,params_size)
!!$    call yaml_parse_from_char_array(parameters,params)
!!$    call f_free(params)
!!$    call yaml_dict_dump_all(parameters,unit=17)
!!$    call dict_free(parameters)

    call dict_init(parameters)
    call set(parameters // DFT_VARIABLES, get_dft_parameters())
    call set(parameters // KPT_VARIABLES, get_kpt_parameters())
    call set(parameters // GEOPT_VARIABLES, get_geopt_parameters())
    call set(parameters // MIX_VARIABLES, get_mix_parameters())
    call set(parameters // SIC_VARIABLES, get_sic_parameters())
    call set(parameters // TDDFT_VARIABLES, get_tddft_parameters())
    call set(parameters // PERF_VARIABLES, get_perf_parameters())

    !call yaml_dict_dump(parameters, comment_key = COMMENT)
    if (INPUT_VAR_NOT_IN_LIST == 0) then
       call f_err_define(err_name='INPUT_VAR_NOT_IN_LIST',&
            err_msg='given value not in allowed list.',&
            err_action='choose a value from the list below.',&
            err_id=INPUT_VAR_NOT_IN_LIST,callback=abort_excl)
    end if
    if (INPUT_VAR_NOT_IN_RANGE == 0) then
       call f_err_define(err_name='INPUT_VAR_NOT_IN_RANGE',&
            err_msg='given value not in allowed range.',&
            err_action='adjust the given value.',&
            err_id=INPUT_VAR_NOT_IN_RANGE)
    end if
    if (INPUT_VAR_ILLEGAL == 0) then
       call f_err_define(err_name='INPUT_VAR_ILLEGAL',&
            err_msg='provided variable is not allowed in this context.',&
            err_action='remove the input variable.',&
            err_id=INPUT_VAR_ILLEGAL,callback=warn_illegal)
    end if
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
    call set(p // HGRIDS, dict_new( &
         & COMMENT   .is. "grid spacing in the three directions (bohr)", &
         & RANGE     .is. list_new(.item."0.", .item."2."), &
         & DEFAULT   .is. list_new(.item."0.45", .item."0.45", .item."0.45"), &
         & "fast"    .is. list_new(.item."0.55", .item."0.55", .item."0.55"), &
         & "accurate".is. list_new(.item."0.30", .item."0.30", .item."0.30") ))

    call set(p // RMULT, dict_new( &
         & COMMENT   .is. "c(f)rmult*radii_cf(:,1(2))=coarse(fine) atom-basec radius", &
         & RANGE     .is. list_new(.item."0.", .item."100."), &
         & DEFAULT   .is. list_new(.item."5.", .item."8.")))

    call set(p // IXC, dict_new(&
         & COMMENT   .is. "exchange-correlation parameter (LDA=1,PBE=11)", &
         & DEFAULT   .is. "1", &
         & "PBE"     .is. "11", &
         & "Hybrid"  .is. "-406"))

    call set(p // NCHARGE, dict_new(&
         & COMMENT   .is. "charge of the system", &
         & RANGE     .is. list_new(.item."-500.", .item."500."), &
         & DEFAULT   .is. "0"))

    call set(p // ELECFIELD, dict_new(&
         & COMMENT   .is. "electric field (Ex,Ey,Ez)", &
         & DEFAULT   .is. list_new(.item."0.", .item."0.", .item."0.")))

    call set(p // NSPIN, dict_new(&
         & COMMENT   .is. "spin polarization", &
         & EXCLUSIVE .is. dict_new("1".is."no spin", "2".is."collinear", "4".is."non-collinear"), &
         & DEFAULT   .is. "1"))

    call set(p // MPOL, dict_new(&
         & COMMENT .is. "total magnetic moment", &
         & DEFAULT .is. "0"))

    call set(p // GNRM_CV, dict_new(&
         & COMMENT   .is. "convergence criterion gradient", &
         & RANGE     .is. list_new(.item."1.e-20", .item."1."), &
         & DEFAULT   .is. "1e-4", &
         & "fast"    .is. "1e-3", &
         & "accurate".is. "1e-5"))

    call set(p // ITERMAX, dict_new(&
         & COMMENT   .is. "max. # of wfn. opt. steps", &
         & RANGE     .is. list_new(.item."0", .item."10000"), &
         & DEFAULT   .is. "50"))

    call set(p // NREPMAX, dict_new(&
         & COMMENT   .is. "max. # of re-diag. runs", &
         & RANGE     .is. list_new(.item."0", .item."1000"), &
         & DEFAULT   .is. "1", &
         & "accurate".is. "10"))

    call set(p // NCONG, dict_new(&
         & COMMENT   .is. "# of CG it. for preconditioning eq.", &
         & RANGE     .is. list_new(.item."0", .item."20"), &
         & DEFAULT   .is. "6"))

    call set(p // IDSX, dict_new(&
         & COMMENT     .is. "wfn. diis history", &
         & RANGE       .is. list_new(.item."0", .item."15"), &
         & DEFAULT     .is. "6", &
         & "low memory".is. "2"))

    call set(p // DISPERSION, dict_new(&
         & COMMENT   .is. "dispersion correction potential (values 1,2,3,4,5), 0=none", &
         & RANGE     .is. list_new(.item."0", .item."5"), &
         & DEFAULT   .is. "0"))

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
    call set(p // INPUTPSIID, dict_new(&
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. excl))

    call dict_init(excl)
    call set(excl // "0", "none")
    call set(excl // "1", "plain text")
    call set(excl // "2", "Fortran binary")
    call set(excl // "3", "ETSF file format")
    call set(p // OUTPUT_WF, dict_new(&
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. excl))

    call dict_init(excl)
    call set(excl // "0", "none")
    call set(excl // "1", "density in plain text")
    call set(excl // "2", "density and potentials in plain text")
    call set(excl // "11", "density in ETSF file format")
    call set(excl // "12", "density and potentials in ETSF file format")
    call set(excl // "21", "density in cube file format")
    call set(excl // "22", "density and potentials in cube file format")
    call set(p // OUTPUT_DENSPOT, dict_new(&
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. excl))

    call set(p // RBUF, dict_new(&
         & COMMENT     .is. "length of the tail (AU)", &
         & RANGE       .is. list_new(.item."0.", .item."10."), &
         & DEFAULT     .is. "0.", &
         & "with tails".is. "5."))

    call set(p // NCONGT, dict_new(&
         & COMMENT   .is. "# tail CG iterations", &
         & RANGE     .is. list_new(.item."0", .item."50"), &
         & DEFAULT   .is. "30"))

    call set(p // NORBV, dict_new(&
         & COMMENT   .is. "Davidson subspace dim.", &
         & RANGE     .is. list_new(.item."-9999", .item."9999"), &
         & DEFAULT   .is. "0"))

    call set(p // NVIRT, dict_new(&
         & COMMENT   .is. "# of opt. orbs", &
         & RANGE     .is. list_new(.item."0", .item."9999"), &
         & DEFAULT   .is. "0"))

    call set(p // NPLOT, dict_new(&
         & COMMENT   .is. "# of plottec orbs", &
         & RANGE     .is. list_new(.item."0", .item."9999"), &
         & DEFAULT   .is. "0"))

    call set(p // DISABLE_SYM, dict_new(&
         & COMMENT   .is. "disable the symmetry detection", &
         & DEFAULT   .is. "No"))

    get_dft_parameters => p
  END FUNCTION get_dft_parameters

  function get_kpt_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_kpt_parameters

    call dict_init(p)

    ! Settings
    call set(p // KPT_METHOD, dict_new( &
         & COMMENT   .is. 'K-point sampling method', &
         & EXCLUSIVE .is. dict_new("auto".is."based on kptrlen", &
         & "MPGrid"  .is."Monkhorst-Pack", "manual".is."based on raw coordinates"), &
         & DEFAULT   .is. "manual" ))

    call set(p // KPTRLEN, dict_new( &
         & COMMENT   .is. 'Equivalent length of K-space resolution (Bohr)', &
         & COND      .is. dict_new(MASTER_KEY .is. KPT_METHOD, &
         &                         WHEN .is. list_new(.item. "auto")), &
         & RANGE     .is. list_new(.item."0", .item."1e4"), &
         & DEFAULT   .is. "0." ))

    call set(p // NGKPT, dict_new( &
         & COMMENT   .is. 'No. of Monkhorst-Pack grid points', &
         & COND      .is. dict_new(MASTER_KEY .is. KPT_METHOD, &
         &                         WHEN .is. list_new(.item. "MPGrid")), &
         & RANGE     .is. list_new(.item."1", .item."10000"), &
         & DEFAULT   .is. list_new(.item."1", .item."1", .item."1") ))

    call set(p // SHIFTK, dict_new( &
         & COMMENT   .is. 'Grid shifts', &
         & COND      .is. dict_new(MASTER_KEY .is. KPT_METHOD, &
         &                         WHEN .is. list_new(.item. "MPGrid")), &
         & DEFAULT   .is. list_new(.item."0.", .item."0.", .item."0.") ))

    call set(p // KPT, dict_new( &
         & COMMENT   .is. 'Kpt coordinates', &
         & COND      .is. dict_new(MASTER_KEY .is. KPT_METHOD, &
         &                         WHEN .is. list_new(.item. "manual")), &
         & DEFAULT   .is. list_new(.item.list_new(.item."0.", .item."0.", .item."0.")) ))

    call set(p // WKPT, dict_new( &
         & COMMENT   .is. 'Kpt weights', &
         & COND      .is. dict_new(MASTER_KEY .is. KPT_METHOD, &
         &                         WHEN .is. list_new(.item. "manual")), &
         & RANGE     .is. list_new(.item."0.", .item."1e4"), &
         & DEFAULT   .is. list_new(.item."1.") ))

    call set(p // BANDS, dict_new( &
         & COMMENT   .is. 'For doing band structure calculation', &
         & DEFAULT   .is. "No" ))

    call set(p // ISEG, dict_new( &
         & COMMENT   .is. 'points for each segment of the BZ path', &
         & COND      .is. dict_new(MASTER_KEY .is. BANDS, &
         &                            WHEN .is. list_new(.item. "Yes")), &
         & DEFAULT   .is. list_new(.item."1") ))

    call set(p // KPTV, dict_new( &
         & COND      .is. dict_new(MASTER_KEY .is. BANDS, &
         &                         WHEN .is. list_new(.item. "Yes")), &
         & DEFAULT   .is. list_new(&
         & .item. list_new(.item."0. ", .item."0. ", .item."0. "), &
         & .item. list_new(.item."0.5", .item."0.5", .item."0.5")) ))

    call set(p // NGRANULARITY, dict_new( &
         & COMMENT   .is. '# of points done for each group', &
         & COND      .is. dict_new(MASTER_KEY .is. BANDS, &
         &                         WHEN .is. list_new(.item. "Yes")), &
         & RANGE     .is. list_new(.item."1", .item."1000"), &
         & DEFAULT   .is. "1" ))

    call set(p // BAND_STRUCTURE_FILENAME, dict_new( &
         & COND      .is. dict_new(MASTER_KEY .is. BANDS, &
         &                         WHEN .is. list_new(.item. "Yes")), &
         & DEFAULT   .is. "''" ))

    get_kpt_parameters => p
  END FUNCTION get_kpt_parameters

  function get_geopt_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_geopt_parameters

    call dict_init(p)

    ! Settings
    call set(p // GEOPT_METHOD, dict_new( &
         & COMMENT   .is. 'Geometry optimisation method', &
         & EXCLUSIVE .is. dict_new(&
         & "none"    .is. "no geometry optimization", &
         & "SDCG"    .is. "a combination of Steepest Descent and Conjugate Gradient", &
         & "VSSD"    .is. "Variable Stepsize Steepest Descent method", &
         & "LBFGS"   .is. "limited-memory BFGS", &
         & "BFGS"    .is. "Broyden-Fletcher-Goldfarb-Shanno", &
         & "PBFGS"   .is. "Same as BFGS with an initial Hessian obtained from a force field", &
         & "AB6MD"   .is. "molecular dynamics from ABINIT", &
         & "DIIS"    .is. "direct inversion of iterative subspace", &
         & "FIRE"    .is. "Fast Inertial Relaxation Engine as described by Bitzek et al."), &
         & DEFAULT   .is. "none" ))

    call set(p // NCOUNT_CLUSTER_X, dict_new( &
         & COMMENT   .is. 'Maximum number of force evaluations', &
         & RANGE     .is. list_new(.item."0", .item."2000"), &
         & PROF_KEY  .is. GEOPT_METHOD, &
         & DEFAULT   .is. "50", &
         & "none"    .is. "1" ))

    call set(p // FRAC_FLUCT, dict_new( &
         & RANGE     .is. list_new(.item."0.", .item."10."), &
         & DEFAULT   .is. "1." ))

    call set(p // FORCEMAX, dict_new( &
         & RANGE     .is. list_new(.item."0.", .item."10."), &
         & DEFAULT   .is. "0." ))

    call set(p // RANDDIS, dict_new( &
         & COMMENT   .is. 'random displacement amplitude', &
         & RANGE     .is. list_new(.item."0.", .item."10."), &
         & DEFAULT   .is. "0." ))

    call set(p // BETAX, dict_new( &
         & COMMENT   .is. 'Stepsize for the geometry optimisation', &
         & COND      .is. dict_new(MASTER_KEY .is. GEOPT_METHOD, &
         &                         WHEN .is. list_new(.item. "SDCG", &
         &                .item."VSSD", .item."LBFGS", .item."BFGS", .item."PBFGS", &
         &                .item."DIIS", .item."FIRE", .item."none")), &
         & PROF_KEY  .is. GEOPT_METHOD, &
         & RANGE     .is. list_new(.item."0.", .item."100."), &
         & DEFAULT   .is. "4.", &
         & "DIIS"    .is. "2." ))

    call set(p // HISTORY, dict_new( &
         & COMMENT   .is. 'history for DIIS method', &
         & COND      .is. dict_new(MASTER_KEY .is. GEOPT_METHOD, &
         &                         WHEN .is. list_new(.item. "DIIS")), &
         & RANGE     .is. list_new(.item."0", .item."100"), &
         & DEFAULT   .is. "4" ))

    call set(p // DTINIT, dict_new( &
         & COMMENT   .is. 'initial time step for the FIRE method', &
         & COND      .is. dict_new(MASTER_KEY .is. GEOPT_METHOD, &
         &                         WHEN .is. list_new(.item. "FIRE")), &
         & RANGE     .is. list_new(.item."0", .item."1e4"), &
         & DEFAULT   .is. "0.75" ))

    call set(p // DTMAX, dict_new( &
         & COMMENT   .is. 'maximal time step for the FIRE method', &
         & COND      .is. dict_new(MASTER_KEY .is. GEOPT_METHOD, &
         &                         WHEN .is. list_new(.item. "FIRE")), &
         & RANGE     .is. list_new(.item."0", .item."1e4"), &
         & DEFAULT   .is. "1.5" ))

    call set(p // IONMOV, dict_new( &
         & COMMENT   .is. 'movement ion method', &
         & COND      .is. dict_new(MASTER_KEY .is. GEOPT_METHOD, &
         &                            WHEN .is. list_new(.item. "AB6MD")), &
         & EXCLUSIVE .is. dict_new("6".is."simple velocity-Verlet molecular dynamic", &
         & "7".is."quenched molecular dynamic", "8".is."Nose-Hoover thermostat", &
         & "9".is."Langevin dynamic", "12".is."Isokinetic ensemble molecular dynamics", &
         & "13".is."Iosthermal/isenthalpic ensemble"), &
         & DEFAULT   .is. "6" ))

    call set(p // DTION, dict_new( &
         & COMMENT   .is. 'time step - Atomic Units (20.670689 AU=0.5 fs)', &
         & COND      .is. dict_new(MASTER_KEY .is. GEOPT_METHOD, &
         &                         WHEN .is. list_new(.item. "AB6MD")), &
         & RANGE     .is. list_new(.item."0", .item."1e3"), &
         & DEFAULT   .is. "20.670689" ))

    call set(p // MDITEMP, dict_new( &
         & COMMENT   .is. 'temperature of molecular dynamics', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                         WHEN .is. list_new(.item. "6")), &
         & RANGE     .is. list_new(.item."0", .item."1e9"), &
         & DEFAULT   .is. "300." ))

    call set(p // MDFTEMP, dict_new( &
         & COMMENT   .is. 'final temperature of molecular dynamics', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                         WHEN .is. list_new(.item. "8", .item."9", &
         &                                            .item."12", .item."13")), &
         & RANGE     .is. list_new(.item."0", .item."1e9"), &
         & DEFAULT   .is. "300." ))

    call set(p // NOSEINERT, dict_new( &
         & COMMENT   .is. 'thermostat inertia coefficient for Nose_Hoover dynamics', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                         WHEN .is. list_new(.item. "8")), &
         & RANGE     .is. list_new(.item."0", .item."1e9"), &
         & DEFAULT   .is. "1e5" ))

    call set(p // FRICTION, dict_new( &
         & COMMENT   .is. 'Friction coefficient for Langevin dynamics', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                         WHEN .is. list_new(.item. "9")), &
         & RANGE     .is. list_new(.item."0", .item."1e5"), &
         & DEFAULT   .is. "1e-3" ))

    call set(p // MDWALL, dict_new( &
         & COMMENT   .is. 'Distance in bohr where atoms can bounce for Langevin dynamics', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                         WHEN .is. list_new(.item. "9")), &
         & RANGE     .is. list_new(.item."0", .item."1e5"), &
         & DEFAULT   .is. "1e4" ))

    call set(p // QMASS, dict_new( &
         & COMMENT   .is. 'Mass of each thermostat (isothermal/isenthalpic ensemble)', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                            WHEN .is. list_new(.item. "13")), &
         & RANGE     .is. list_new(.item."0", .item."1e9"), &
         & DEFAULT   .is. list_new(.item."0.") ))

    call set(p // BMASS, dict_new( &
         & COMMENT   .is. 'Barostat masses (isothermal/isenthalpic ensemble)', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                         WHEN .is. list_new(.item. "13")), &
         & RANGE     .is. list_new(.item."0", .item."1e9"), &
         & DEFAULT   .is. "10." ))

    call set(p // VMASS, dict_new( &
         & COMMENT   .is. 'Barostat masses (isothermal/isenthalpic ensemble)', &
         & COND      .is. dict_new(MASTER_KEY .is. IONMOV, &
         &                         WHEN .is. list_new(.item. "13")), &
         & RANGE     .is. list_new(.item."0", .item."1e9"), &
         & DEFAULT   .is. "1." ))

    get_geopt_parameters => p
  END FUNCTION get_geopt_parameters

  function get_mix_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_mix_parameters

    call dict_init(p)

    ! Settings
    call set(p // ISCF, dict_new( &
         & COMMENT   .is. 'mixing parameters', &
         & EXCLUSIVE .is. dict_new(&
         & "-1" .is. "reserved, don't use it.", &
         & "0"  .is. "direct minimization", &
         & "1"  .is. "get the largest eigenvalue of the SCF cycle", &
         & "2"  .is. "simple mixing of the potential", &
         & "3"  .is. "Anderson mixing of the potential", &
         & "4"  .is. "Anderson mixing of the potential based on the two previous iterations", &
         & "5"  .is. "CG based on the minim. of the energy with respect to the potential", &
         & "7"  .is. "Pulay mixing of the potential based on npulayit previous iterations", &
         & "12" .is. "simple mixing of the density", &
         & "13" .is. "Anderson mixing of the density", &
         & "14" .is. "Anderson mixing of the density based on the two previous iterations", &
         & "15" .is. "CG based on the minim. of the energy with respect to the density", &
         & "17" .is. "Pulay mixing of the density"), &
         & DEFAULT   .is. "0" ))

    call set(p // ITRPMAX, dict_new( &
         & COMMENT   .is. 'maximum number of diagonalisation iterations', &
         & RANGE     .is. list_new(.item."0", .item."10000"), &
         & DEFAULT   .is. "1" ))

    call set(p // RPNRM_CV, dict_new( &
         & COMMENT   .is. 'stop criterion on the residue of potential or density', &
         & RANGE     .is. list_new(.item."0.", .item."10."), &
         & DEFAULT   .is. "1e-4" ))

    call set(p // NORBSEMPTY, dict_new( &
         & COMMENT   .is. 'No. of additional bands', &
         & RANGE     .is. list_new(.item."0", .item."10000"), &
         & DEFAULT   .is. "0" ))

    call set(p // TEL, dict_new( &
         & COMMENT   .is. 'electronic temperature', &
         & RANGE     .is. list_new(.item."0.", .item."1e6"), &
         & DEFAULT   .is. "0." ))

    call set(p // OCCOPT, dict_new( &
         & COMMENT   .is. 'smearing method', &
         & EXCLUSIVE .is. dict_new(&
         &       "1" .is. "error function smearing", &
         &       "2" .is. "normal Fermi distribution", &
         &       "3" .is. "Marzari's cold smearing a=-.5634 (bumb minimization)", &
         &       "4" .is. "Marzari's cold smearing with a=-.8165 (monotonic tail)", &
         &       "5" .is. "Methfessel and Paxton (same as cold with a=0)"), &
         & DEFAULT   .is. "1" ))

    call set(p // ALPHAMIX, dict_new( &
         & COMMENT   .is. 'Multiplying factors for the mixing', &
         & RANGE     .is. list_new(.item."0", .item."1."), &
         & DEFAULT   .is. "0." ))

    call set(p // ALPHADIIS, dict_new( &
         & COMMENT   .is. 'Multiplying factors for the electronic DIIS', &
         & RANGE     .is. list_new(.item."0.", .item."10."), &
         & DEFAULT   .is. "2." ))

    get_mix_parameters => p
  END FUNCTION get_mix_parameters

  function get_sic_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_sic_parameters

    call dict_init(p)

    ! Settings
    call set(p // SIC_APPROACH, dict_new( &
         & COMMENT   .is. 'SIC method', &
         & EXCLUSIVE .is. dict_new(&
         &    "none" .is. "no self-interaction correction", &
         &    "PZ"   .is. "Perdew-Zunger SIC scheme", &
         &    "NK"   .is. "Non-Koopmans correction (Experimental)"), &
         & DEFAULT.is. "none" ))

    call set(p // SIC_ALPHA, dict_new( &
         & COMMENT   .is. 'SIC downscaling parameter', &
         & RANGE     .is. list_new(.item."0", .item."1."), &
         & DEFAULT   .is. "0." ))

    call set(p // SIC_FREF, dict_new( &
         & COMMENT   .is. 'reference occupation fref (for NK case only)', &
         & COND      .is. dict_new(MASTER_KEY .is. SIC_APPROACH, &
         &                            WHEN .is. list_new(.item. "NK")), &
         & RANGE     .is. list_new(.item."0.", .item."1."), &
         & DEFAULT   .is. "0." ))

    get_sic_parameters => p
  END FUNCTION get_sic_parameters

  function get_tddft_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_tddft_parameters

    call dict_init(p)

    ! Settings
    call set(p // TDDFT_APPROACH, dict_new( &
         & COMMENT   .is. 'TDDFT method', &
         & EXCLUSIVE .is. dict_new(&
         &    "none" .is. "no tddft post-processing", &
         &    "TDA"  .is. "???"), &
         & DEFAULT.is. "none" ))

    get_tddft_parameters => p
  END FUNCTION get_tddft_parameters

  function get_perf_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_perf_parameters

    call dict_init(p)

    ! Settings
    call set(p // DEBUG, dict_new( &
         & COMMENT .is. 'debug option', &
         & DEFAULT .is. "No" ))

    call set(p // FFTCACHE, dict_new( &
         & COMMENT .is. 'cache size for the FFT', &
         & DEFAULT .is. "8192" ))

    call set(p // ACCEL, dict_new( &
         & COMMENT   .is. 'acceleration', &
         & EXCLUSIVE .is. dict_new(&
         & "no"      .is. "no material acceleration", &
         & "CUDAGPU" .is. "CUDA", &
         & "OCLGPU"  .is. "OpenCL on GPU", &
         & "OCLCPU"  .is. "OpenCL on CPU", &
         & "OCLACC"  .is. "???"), &
         & DEFAULT   .is. "No" ))

    call set(p // OCL_PLATFORM, dict_new( &
         & COMMENT   .is. 'Chosen OCL platform', &
         & DEFAULT.is. " " ))

    call set(p // OCL_DEVICES, dict_new( &
         & COMMENT   .is. 'Chosen OCL devices', &
         & DEFAULT.is. " " ))

    call set(p // BLAS, dict_new( &
         & COMMENT   .is. 'CUBLAS acceleration', &
         & DEFAULT.is. "No" ))

    call set(p // PROJRAD, dict_new( &
         & COMMENT   .is. 'Radius of the projector as a function of the maxrad', &
         & RANGE     .is. list_new(.item."0.", .item."100."), &
         & DEFAULT.is. "15." ))

    call set(p // EXCTXPAR, dict_new( &
         & COMMENT   .is. 'Exact exchange parallelisation scheme', &
         & EXCLUSIVE .is. dict_new(&
         & "BC"      .is. "???", &
         & "OP2P"    .is. "???"), &
         & DEFAULT.is. "OP2P" ))

    call set(p // IG_DIAG, dict_new( &
         & COMMENT .is. 'Input guess: (T:Direct, F:Iterative) diag. of Ham.', &
         & DEFAULT .is. "Yes" ))

    call set(p // IG_NORBP, dict_new( &
         & COMMENT   .is. 'Input guess: Orbitals per process for iterative diag.', &
         & RANGE     .is. list_new(.item."1", .item."1000"), &
         & DEFAULT.is. "5" ))

    call set(p // IG_BLOCKS, dict_new( &
         & COMMENT   .is. 'Input guess: Block sizes for orthonormalisation', &
         & RANGE     .is. list_new(.item."1", .item."100000"), &
         & DEFAULT   .is. list_new(.item."300", .item."800") ))

    call set(p // IG_TOL, dict_new( &
         & COMMENT   .is. 'Input guess: Tolerance criterion', &
         & RANGE     .is. list_new(.item."0.", .item."1."), &
         & DEFAULT.is. "1e-4" ))

    call set(p // METHORTHO, dict_new( &
         & COMMENT   .is. 'Orthogonalisation', &
         & EXCLUSIVE .is. dict_new(&
         & "0"    .is. "Cholesky", &
         & "1"    .is. "GS/Chol", &
         & "2"    .is. "Loewdin"), &
         & DEFAULT.is. "0" ))

    call set(p // RHO_COMMUN, dict_new( &
         & COMMENT   .is. 'Density communication scheme (DBL, RSC, MIX)', &
         & DEFAULT.is. "DEF" ))

    call set(p // PSOLVER_GROUPSIZE, dict_new( &
         & COMMENT   .is. 'Size of Poisson Solver taskgroups (0=nproc)', &
         & RANGE     .is. list_new(.item."0", .item."1000000"), &
         & DEFAULT.is. "0" ))

    call set(p // PSOLVER_ACCEL, dict_new( &
         & COMMENT   .is. 'Acceleration of the Poisson Solver (0=none, 1=CUDA)', &
         & EXCLUSIVE .is. dict_new(&
         & "0"    .is. "none", &
         & "1"    .is. "CUDA"), &
         & DEFAULT.is. "0" ))

    call set(p // UNBLOCK_COMMS, dict_new( &
         & COMMENT   .is. 'Overlap Communications of fields (OFF,DEN,POT)', &
         & EXCLUSIVE .is. dict_new(&
         & "OFF"    .is. "synchronous run", &
         & "DEN"    .is. "???", &
         & "POT"    .is. "???"), &
         & DEFAULT.is. "OFF" ))

    call set(p // LINEAR, dict_new( &
         & COMMENT   .is. 'Linear Input Guess approach', &
         & EXCLUSIVE .is. dict_new(&
         & "OFF"    .is. "???", &
         & "LIG"    .is. "???", &
         & "FUL"    .is. "???", &
         & "TMO"    .is. "???"), &
         & DEFAULT.is. "OFF" ))

    call set(p // TOLSYM, dict_new( &
         & COMMENT   .is. 'Tolerance for symmetry detection', &
         & RANGE     .is. list_new(.item."0.", .item."1."), &
         & DEFAULT.is. "1e-8" ))

    call set(p // SIGNALING, dict_new( &
         & COMMENT   .is. 'Expose calculation results on Network', &
         & DEFAULT.is. "No" ))

    call set(p // SIGNALTIMEOUT, dict_new( &
         & COMMENT   .is. 'Time out on startup for signal connection (in seconds)', &
         & RANGE     .is. list_new(.item."-1", .item."3600"), &
         & DEFAULT.is. "0" ))

    call set(p // DOMAIN, dict_new( &
         & COMMENT   .is. 'Domain to add to the hostname to find the IP', &
         & DEFAULT.is. " " ))

    call set(p // INGUESS_GEOPT, dict_new( &
         & COMMENT   .is. '0= wavlet input guess, 1= real space input guess', &
         & EXCLUSIVE .is. dict_new(&
         & "0"    .is. "wavlet input guess", &
         & "1"    .is. "real space input guess"), &
         & DEFAULT.is. "0" ))

    call set(p // STORE_INDEX, dict_new( &
         & COMMENT   .is. 'linear scaling: store indices or recalculate them', &
         & DEFAULT.is. "Yes" ))

    call set(p // VERBOSITY, dict_new( &
         & COMMENT   .is. 'verbosity of the output', &
         & EXCLUSIVE .is. dict_new(&
         & "0"    .is. "???", &
         & "1"    .is. "???", &
         & "2"    .is. "???", &
         & "3"    .is. "???"), &
         & DEFAULT.is. "2" ))

    call set(p // OUTDIR, dict_new( &
         & COMMENT   .is. 'Writing directory', &
         & DEFAULT.is. "." ))

    call set(p // PSP_ONFLY, dict_new( &
         & COMMENT   .is. 'Calculate pseudopotential projectors on the fly', &
         & DEFAULT.is. "Yes" ))

    call set(p // PDSYEV_BLOCKSIZE, dict_new( &
         & COMMENT   .is. 'SCALAPACK linear scaling blocksize', &
         & DEFAULT.is. "-8" ))

    call set(p // PDGEMM_BLOCKSIZE, dict_new( &
         & COMMENT   .is. 'SCALAPACK linear scaling blocksize', &
         & DEFAULT.is. "-8" ))

    call set(p // MAXPROC_PDSYEV, dict_new( &
         & COMMENT   .is. 'SCALAPACK linear scaling max num procs', &
         & DEFAULT.is. "4" ))

    call set(p // MAXPROC_PDGEMM, dict_new( &
         & COMMENT   .is. 'SCALAPACK linear scaling max num procs', &
         & DEFAULT.is. "4" ))

    call set(p // EF_INTERPOL_DET, dict_new( &
         & COMMENT   .is. 'FOE: max determinant of cubic interpolation matrix', &
         & RANGE     .is. list_new(.item."0.", .item."1."), &
         & DEFAULT.is. "1e-20" ))

    call set(p // EF_INTERPOL_CHARGEDIFF, dict_new( &
         & COMMENT   .is. 'FOE: max charge difference for interpolation', &
         & RANGE     .is. list_new(.item."0.", .item."1000."), &
         & DEFAULT.is. "10." ))

    call set(p // MIXING_AFTER_INPUTGUESS, dict_new( &
         & COMMENT   .is. 'mixing step after linear input guess (T/F)', &
         & DEFAULT.is. "Yes" ))

    call set(p // ITERATIVE_ORTHOGONALIZATION, dict_new( &
         & COMMENT   .is. 'iterative_orthogonalization for input guess orbitals', &
         & DEFAULT.is. "No" ))

    call set(p // CHECK_SUMRHO,dict_new( &
         COMMENT .is. 'enables linear sumrho check', &
         DEFAULT .is. '2', &
         EXCLUSIVE .is. dict_new('0' .is. 'no check',&
                                 '1' .is. 'light check',&
                                 '2' .is. 'full check')))
    !the opportunity of entering this dictionary already in yaml format should be discussed
    !for example the above variable becomes:
    !check_sumrho:
    !   COMMENT: enables linear sumrho check
    !   DEFAULT: 2
    !   EXCLUSIVE: {0: no check,1: light check, 2: full check}
    !which is *way* simpler


    get_perf_parameters => p
  END FUNCTION get_perf_parameters

  subroutine input_keys_dump_def(fname, file)
    use dictionaries
    use yaml_output
    implicit none
    character(len = *), intent(in) :: fname
    character(len = *), intent(in), optional :: file
    !local variables
    integer, parameter :: unt=789159 !to be sure is not opened
    integer :: iunit_def, ierr

    ! Switch YAML output stream
    call yaml_get_default_stream(iunit_def)
    call yaml_set_stream(unit=unt,filename=trim(fname),tabbing=0,record_length=100)

    if (ierr == 0) then
       call input_keys_init()
       if (present(file)) then
          if (has_key(parameters, file)) then
             call yaml_dict_dump(parameters // file)
          else
             call yaml_dict_dump(parameters)
          end if
       else
          call yaml_dict_dump(parameters)
       end if
       call input_keys_finalize()
    end if
    call yaml_close_stream(unit=unt)
    ! Set back normal YAML output
    call yaml_set_default_stream(iunit_def,ierr)
  end subroutine input_keys_dump_def

  !> get for each keys available profiles.
  function input_keys_get_profiles(file)
    use dictionaries
    implicit none
    character(len = *), intent(in) :: file
    type(dictionary), pointer :: input_keys_get_profiles

    type(dictionary), pointer :: p
    integer :: i
    character(max_field_length), dimension(:), allocatable :: keys

    call input_keys_init()
    
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

    call input_keys_finalize()

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
         call generate(dict, var(i), ref // var(i))
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
      character(max_field_length), intent(in) :: key

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
            call add(dict // key // "profiles", "'" // trim(keys(i)) // "'")
            !call dict_copy(dict // key // "profiles" // keys(i), ref // keys(i))
!!$         else if (trim(keys(i)) == EXCLUSIVE) then
!!$            ! Add exclusive values.
!!$            call dict_copy(dict // key // "allowed values", ref // EXCLUSIVE)
!!$         else if (trim(keys(i)) == RANGE) then
!!$            ! Add range definition
!!$            call dict_copy(dict // key // "within range", ref // RANGE)
!!$         else if (trim(keys(i)) == DEFAULT) then
!!$            ! Add range definition
!!$            call dict_copy(dict // key // "default value", ref // DEFAULT)
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

  subroutine input_keys_fill_all(dict)
    use dictionaries
    !use yaml_output
    implicit none
    type(dictionary), pointer :: dict

    !call yaml_dict_dump(dict)

    ! Check and complete dictionary.
    call input_keys_init()

    call input_keys_fill(dict//PERF_VARIABLES, PERF_VARIABLES)
    call input_keys_fill(dict//DFT_VARIABLES, DFT_VARIABLES)
    call input_keys_fill(dict//KPT_VARIABLES, KPT_VARIABLES)
    call input_keys_fill(dict//GEOPT_VARIABLES, GEOPT_VARIABLES)
    call input_keys_fill(dict//MIX_VARIABLES, MIX_VARIABLES)
    call input_keys_fill(dict//SIC_VARIABLES, SIC_VARIABLES)
    call input_keys_fill(dict//TDDFT_VARIABLES, TDDFT_VARIABLES)

    call input_keys_finalize()
  end subroutine input_keys_fill_all

  subroutine input_keys_fill(dict, file, profile)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file
    character(len = *), intent(in), optional :: profile

    integer :: i
    type(dictionary), pointer :: ref
    character(len=max_field_length), dimension(:), allocatable :: keys

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
            & err_msg = trim(file) // "/" // trim(key) // " is not allowed in this context.")) return

       val = dict // key
       ! There is already a value in dict.
       if (has_key(ref, val)) then
          ! The key was asking for a profile, we copy it.
          call dict_copy(dict // key, ref // val)
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
       call dict_copy(dict // key, ref // profile_)
    end if

    ! Copy the comment.
    if (has_key(ref, COMMENT)) &
         & call dict_copy(dict // (trim(key) // ATTRS) // COMMENT, ref // COMMENT)
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
         if (var < rg(1) .or. var > rg(2)) then
            val = dict
            call f_err_throw(err_id = INPUT_VAR_NOT_IN_RANGE, &
              & err_msg = trim(key) // " = '" // trim(val) // &
              & "' not in range.")
            return
         end if
      end if
    end subroutine validate
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
      character(max_field_length) :: descr, tag
      character(max_field_length), dimension(:), allocatable :: keys

      if (index(dict%data%key, ATTRS) > 0) return

      descr = " "
      tag = " "
      parent => dict%parent
      if (associated(parent)) then
         if (has_key(parent, trim(dict%data%key) // ATTRS)) then
            attr => parent // (trim(dict%data%key) // ATTRS)
            if (has_key(attr, COMMENT)) descr = attr // COMMENT
            !if (has_key(attr, PROF_KEY)) tag = attr // PROF_KEY
         end if
      end if

      if (dict_len(dict) > 0) then
         ! List case.
         flow = (.not.associated(dict%child%child))
         if (.not.flow .and. trim(descr) /= "") then
            call yaml_open_sequence(trim(dict%data%key), tag = tag, advance = "no")
            call yaml_comment(trim(descr), tabbing = 50)
         else
            call yaml_open_sequence(trim(dict%data%key), tag = tag, flow=flow)
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
      else if (dict_size(dict) > 0) then
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
      else if (associated(dict)) then
         ! Leaf case.
         if (dict%data%item >= 0) then
            call yaml_sequence(trim(dict%data%value))
         else
            if (trim(descr) /= "") then
               call yaml_map(trim(dict%data%key), trim(dict%data%value), tag = tag, advance = "no")
               call yaml_comment(trim(descr), tabbing = 50)
            else
               call yaml_map(trim(dict%data%key), trim(dict%data%value), tag = tag)
            end if
         end if
      end if

    end subroutine dict_dump_
  end subroutine input_keys_dump

end module module_input_keys
