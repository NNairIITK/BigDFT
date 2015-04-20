!> @file
!! Contains public parameters that have to be used here and there in the code
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module public_keys
  implicit none

  public ! guess why?

  character(len = *), parameter :: MODE_VARIABLES = "mode"
  character(len = *), parameter :: METHOD_KEY = "method"
  character(len = *), parameter :: RUN_NAME_KEY = "name"

  character(len = *), parameter :: POSINP = "posinp"
  character(len = *), parameter :: OCCUPATION = "occupation"
  character(len = *), parameter :: IG_OCCUPATION = "ig_occupation"
  character(len = *), parameter :: DFT_VARIABLES = "dft"
  character(len = *), parameter :: HGRIDS = "hgrids"
  character(len = *), parameter :: RMULT = "rmult"
  character(len = *), parameter :: IXC = "ixc"
  character(len = *), parameter :: NCHARGE = "qcharge"
  character(len = *), parameter :: ELECFIELD = "elecfield"
  character(len = *), parameter :: NSPIN = "nspin", MPOL = "mpol"
  character(len = *), parameter :: GNRM_CV = "gnrm_cv"
  character(len = *), parameter :: GNRM_IG = "gnrm_ig"
  character(len = *), parameter :: NIT_IG = "nit_ig"
  character(len = *), parameter :: ITERMAX = "itermax",ITERMIN = "itermin", NREPMAX = "nrepmax"
  character(len = *), parameter :: NCONG = "ncong", IDSX = "idsx"
  character(len = *), parameter :: DISPERSION = "dispersion"
  character(len = *), parameter :: INPUTPSIID = "inputpsiid"
  character(len = *), parameter :: OUTPUT_WF = "output_wf"
  character(len = *), parameter :: OUTPUT_DENSPOT = "output_denspot"
  character(len = *), parameter :: RBUF = "rbuf"
  character(len = *), parameter :: NCONGT = "ncongt"
  character(len = *), parameter :: NORBV = "norbv", NVIRT = "nvirt"
  character(len = *), parameter :: NPLOT = "nplot"
  character(len = *), parameter :: DISABLE_SYM = "disablesym"
  character(len = *), parameter :: SOLVENT = "solvent"

  character(len = *), parameter :: KPT_VARIABLES = "kpt"
  character(len = *), parameter :: KPT_METHOD = "method"
  character(len = *), parameter :: KPTRLEN = "kptrlen"
  character(len = *), parameter :: NGKPT = "ngkpt"
  character(len = *), parameter :: SHIFTK = "shiftk"
  character(len = *), parameter :: KPT = "kpt"
  character(len = *), parameter :: WKPT = "wkpt"
  character(len = *), parameter :: BANDS = "bands"
  character(len = *), parameter :: ISEG = "iseg"
  character(len = *), parameter :: KPTV = "kptv"
  character(len = *), parameter :: NGRANULARITY = "ngranularity"
  character(len = *), parameter :: BAND_STRUCTURE_FILENAME = "band_structure_filename"

  character(len = *), parameter :: GEOPT_VARIABLES = "geopt"
  character(len = *), parameter :: GEOPT_METHOD = "method"
  character(len = *), parameter :: NCOUNT_CLUSTER_X = "ncount_cluster_x"
  character(len = *), parameter :: FRAC_FLUCT = "frac_fluct"
  character(len = *), parameter :: FORCEMAX = "forcemax"
  character(len = *), parameter :: RANDDIS = "randdis"
  character(len = *), parameter :: IONMOV = "ionmov"
  character(len = *), parameter :: DTION = "dtion"
  character(len = *), parameter :: MDITEMP = "mditemp"
  character(len = *), parameter :: MDFTEMP = "mdftemp"
  character(len = *), parameter :: NOSEINERT = "noseinert"
  character(len = *), parameter :: FRICTION = "friction"
  character(len = *), parameter :: MDWALL = "mdwall"
  character(len = *), parameter :: NNOS = "nnos"
  character(len = *), parameter :: QMASS = "qmass"
  character(len = *), parameter :: BMASS = "bmass"
  character(len = *), parameter :: VMASS = "vmass"
  character(len = *), parameter :: BETAX = "betax"
  character(len = *), parameter :: HISTORY = "history"
  character(len = *), parameter :: DTINIT = "dtinit"
  character(len = *), parameter :: DTMAX = "dtmax"
  character(len = *), parameter :: NEB_RESTART = "restart"
  character(len = *), parameter :: NEB_CLIMBING = "climbing"
  character(len = *), parameter :: EXTREMA_OPT = "extrema_opt"
  character(len = *), parameter :: NEB_METHOD = "neb_method"
  character(len = *), parameter :: TEMP = "temp"
  character(len = *), parameter :: NEB_DAMP = "damp"
  character(len = *), parameter :: SPRINGS_K = "springs_k"
  character(len = *), parameter :: FIX_TOL = "fix_tol"
  character(len = *), parameter :: NIMG = "nimg"
  !SQNM parameters:
  character(len = *), parameter :: NHISTX = "nhistx"
  character(len = *), parameter :: BIOMODE = "biomode"
  character(len = *), parameter :: BETA_STRETCHX = "beta_stretchx"
  character(len = *), parameter :: MAXRISE = "maxrise"
  character(len = *), parameter :: CUTOFFRATIO = "cutoffratio"
  character(len = *), parameter :: STEEPTHRESH = "steepthresh"
  character(len = *), parameter :: TRUSTR = "trustr"

  !Force field parameter keyword
  character(len = *), parameter :: MM_PARAMSET = "mm_paramset" !for hard-coded parameter sets
  character(len = *), parameter :: MM_PARAMFILE = "mm_paramfile" !for parameter sets given by file


  character(len = *), parameter :: MIX_VARIABLES = "mix"
  character(len = *), parameter :: ISCF = "iscf"
  character(len = *), parameter :: ITRPMAX = "itrpmax"
  character(len = *), parameter :: RPNRM_CV = "rpnrm_cv"
  character(len = *), parameter :: NORBSEMPTY = "norbsempty"
  character(len = *), parameter :: TEL = "tel"
  character(len = *), parameter :: OCCOPT = "occopt"
  character(len = *), parameter :: ALPHAMIX = "alphamix"
  character(len = *), parameter :: ALPHADIIS = "alphadiis"

  character(len = *), parameter :: SIC_VARIABLES = "sic"
  character(len = *), parameter :: SIC_APPROACH = "sic_approach"
  character(len = *), parameter :: SIC_ALPHA = "sic_alpha"
  character(len = *), parameter :: SIC_FREF = "sic_fref"

  character(len = *), parameter :: TDDFT_VARIABLES = "tddft"
  character(len = *), parameter :: TDDFT_APPROACH = "tddft_approach"

  character(len = *), parameter :: PERF_VARIABLES = "perf"
  character(len = *), parameter :: DEBUG = "debug"
  character(len = *), parameter :: FFTCACHE = "fftcache"
  character(len = *), parameter :: ACCEL = "accel"
  character(len = *), parameter :: OCL_PLATFORM = "ocl_platform"
  character(len = *), parameter :: OCL_DEVICES = "ocl_devices"
  character(len = *), parameter :: BLAS = "blas"
  character(len = *), parameter :: PROJRAD = "projrad"
  character(len = *), parameter :: EXCTXPAR = "exctxpar"
  character(len = *), parameter :: IG_DIAG = "ig_diag"
  character(len = *), parameter :: IG_NORBP = "ig_norbp"
  character(len = *), parameter :: IG_BLOCKS = "ig_blocks"
  character(len = *), parameter :: IG_TOL = "ig_tol"
  character(len = *), parameter :: METHORTHO = "methortho"
  character(len = *), parameter :: RHO_COMMUN = "rho_commun"
  character(len = *), parameter :: PSOLVER_GROUPSIZE = "psolver_groupsize"
  character(len = *), parameter :: PSOLVER_ACCEL = "psolver_accel"
  character(len = *), parameter :: UNBLOCK_COMMS = "unblock_comms"
  character(len = *), parameter :: LINEAR = "linear"
  character(len = *), parameter :: TOLSYM = "tolsym"
  character(len = *), parameter :: SIGNALING = "signaling"
  character(len = *), parameter :: SIGNALTIMEOUT = "signaltimeout"
  character(len = *), parameter :: DOMAIN = "domain"
  character(len = *), parameter :: INGUESS_GEOPT = "inguess_geopt"
  character(len = *), parameter :: STORE_INDEX = "store_index"
  character(len = *), parameter :: VERBOSITY = "verbosity"
  character(len = *), parameter :: PSP_ONFLY = "psp_onfly"
  character(len = *), parameter :: PDSYEV_BLOCKSIZE = "pdsyev_blocksize"
  character(len = *), parameter :: PDGEMM_BLOCKSIZE = "pdgemm_blocksize"
  character(len = *), parameter :: MAXPROC_PDSYEV = "maxproc_pdsyev"
  character(len = *), parameter :: MAXPROC_PDGEMM = "maxproc_pdgemm"
  character(len = *), parameter :: EF_INTERPOL_DET = "ef_interpol_det"
  character(len = *), parameter :: EF_INTERPOL_CHARGEDIFF = "ef_interpol_chargediff"
  character(len = *), parameter :: MIXING_AFTER_INPUTGUESS = "mixing_after_inputguess"
  character(len = *), parameter :: ITERATIVE_ORTHOGONALIZATION = "iterative_orthogonalization"
  character(len = *), parameter :: MULTIPOLE_PRESERVING = "multipole_preserving"
  character(len = *), parameter :: MP_ISF = "mp_isf"
  character(len = *), parameter :: CHECK_SUMRHO = "check_sumrho"
  character(len = *), parameter :: CHECK_OVERLAP = "check_overlap"
  character(len = *), parameter :: EXPERIMENTAL_MODE = "experimental_mode"
  character(len = *), parameter :: WRITE_ORBITALS = "write_orbitals"
  character(len = *), parameter :: EXPLICIT_LOCREGCENTERS = "explicit_locregcenters"
  character(len = *), parameter :: CALCULATE_KS_RESIDUE = "calculate_KS_residue"
  character(len = *), parameter :: INTERMEDIATE_FORCES = "intermediate_forces"
  character(len = *), parameter :: KAPPA_CONV = "kappa_conv"
  character(len = *), parameter :: EVBOUNDS_NSATUR = "evbounds_nsatur"
  character(len = *), parameter :: EVBOUNDSSHRINK_NSATUR = "evboundsshrink_nsatur"
  character(len = *), parameter :: METHOD_UPDATEKERNEL = "method_updatekernel"
  character(len = *), parameter :: PURIFICATION_QUICKRETURN = "purification_quickreturn"
  character(len = *), parameter :: ADJUST_FOE_TEMPERATURE = "adjust_FOE_temperature"
  character(len = *), parameter :: CALCULATE_GAP = "calculate_gap"
  character(len = *), parameter :: LOEWDIN_CHARGE_ANALYSIS = "loewdin_charge_analysis"
  character(len = *), parameter :: CHECK_MATRIX_COMPRESSION = "check_matrix_compression"
  character(len = *), parameter :: CORRECTION_CO_CONTRA = "correction_co_contra"
  character(len = *), parameter :: GPS_METHOD = "gps_method"
  character(len = *), parameter :: FOE_GAP = "foe_gap"

  !keys for linear input variables
  !level keys
  character(len=*), parameter :: LIN_GENERAL     ='lin_general'
  character(len=*), parameter :: LIN_BASIS       ='lin_basis'
  character(len=*), parameter :: LIN_KERNEL      ='lin_kernel'
  character(len=*), parameter :: LIN_BASIS_PARAMS='lin_basis_params'
  character(len=*), parameter :: FRAG_VARIABLES  ='frag'
  character(len=*), parameter :: HYBRID          ='hybrid'
  character(len=*), parameter :: LINEAR_METHOD   ='linear_method'
  character(len=*), parameter :: MIXING_METHOD   ='mixing_method'
  character(len=*), parameter :: NIT             ='nit'
  character(len=*), parameter :: NSTEP           ='nstep'
  character(len=*), parameter :: IDSX_COEFF      ='idsx_coeff'
  character(len=*), parameter :: GNRM_CV_COEFF   ='gnrm_cv_coeff'
  character(len=*), parameter :: DELTAE_CV       ='deltae_cv'
  character(len=*), parameter :: GNRM_DYN        ='gnrm_dyn'
  character(len=*), parameter :: MIN_GNRM_FOR_DYNAMIC = 'min_gnrm_for_dynamic'
  character(len=*), parameter :: CONF_DAMPING    ='conf_damping'
  character(len=*), parameter :: TAYLOR_ORDER    ='taylor_order'
  character(len=*), parameter :: CALC_DIPOLE     ='calc_dipole'
  character(len=*), parameter :: CALC_PULAY      ='calc_pulay'
  character(len=*), parameter :: SUBSPACE_DIAG   ='subspace_diag'
  character(len=*), parameter :: ALPHA_DIIS      ='alpha_diis'
  character(len=*), parameter :: ALPHA_SD        ='alpha_sd'
  character(len=*), parameter :: ALPHA_SD_COEFF  ='alpha_sd_coeff'
  character(len=*), parameter :: ALPHA_FIT_COEFF ='alpha_fit_coeff'
  character(len=*), parameter :: NSTEP_PREC      ='nstep_prec'
  character(len=*), parameter :: EVAL_RANGE_FOE  ='eval_range_foe'
  character(len=*), parameter :: FSCALE_FOE      ='fscale_foe'
  character(len=*), parameter :: AO_CONFINEMENT  ='ao_confinement'
  character(len=*), parameter :: CONFINEMENT     ='confinement'
  character(len=*), parameter :: RLOC            ='rloc'
  character(len=*), parameter :: RLOC_KERNEL     ='rloc_kernel'
  character(len=*), parameter :: RLOC_KERNEL_FOE ='rloc_kernel_foe'
  character(len=*), parameter :: NBASIS          ='nbasis'
  character(len=*), parameter :: TRANSFER_INTEGRALS='transfer_integrals'
  character(len=*), parameter :: CONSTRAINED_DFT  ='constrained_dft'
  character(len=*), parameter :: FIX_BASIS       ='fix_basis' 
  character(len=*), parameter :: CORRECTION_ORTHOCONSTRAINT='correction_orthoconstraint'
  character(len=*), parameter :: FSCALE_LOWERBOUND="fscale_lowerbound"
  character(len=*), parameter :: FSCALE_UPPERBOUND="fscale_upperbound"
  character(len=*), parameter :: EXTRA_STATES="extra_states"
  character(len=*), parameter :: FRAGMENT_NO="Fragment No. "
  character(len=*), parameter :: MAX_INVERSION_ERROR = "max_inversion_error"
  character(len=*), parameter :: FOE_RESTART="FOE_restart"
  character(len=*), parameter :: IMETHOD_OVERLAP = "imethod_overlap"
  character(len=*), parameter :: EXTRA_SHELLS_KEY='empty_shells'
  character(len=*), parameter :: ENABLE_MATRIX_TASKGROUPS='enable_matrix_taskgroups'
  character(len=*), parameter :: HAMAPP_RADIUS_INCR='hamapp_radius_incr'
  character(len=*), parameter :: ADJUST_KERNEL_ITERATIONS='adjust_kernel_iterations'
  character(len=*), parameter :: WF_EXTENT_ANALYSIS='wf_extent_analysis'
  character(len=*), parameter :: CALCULATE_ONSITE_OVERLAP='calculate_onsite_overlap'

  !> Parameters to avoid typos in dictionary keys
  character(len=*), parameter :: ASTRUCT_UNITS = 'units' 
  character(len=*), parameter :: ASTRUCT_CELL = 'cell' 
  character(len=*), parameter :: ASTRUCT_POSITIONS = 'positions' 
  character(len=*), parameter :: ASTRUCT_PROPERTIES = 'properties' 
  character(len=*), parameter :: ASTRUCT_ATT_FROZEN = 'Frozen' 
  character(len=*), parameter :: ASTRUCT_ATT_IGSPIN = 'IGSpin' 
  character(len=*), parameter :: ASTRUCT_ATT_IGCHRG = 'IGChg' 
  character(len=*), parameter :: ASTRUCT_ATT_IXYZ_1 = 'int_ref_atoms_1' 
  character(len=*), parameter :: ASTRUCT_ATT_IXYZ_2 = 'int_ref_atoms_2' 
  character(len=*), parameter :: ASTRUCT_ATT_IXYZ_3 = 'int_ref_atoms_3' 
  character(len=*), parameter :: ASTRUCT_ATT_RXYZ_INT_1 = 'rxyz_int_atoms_1' 
  character(len=*), parameter :: ASTRUCT_ATT_RXYZ_INT_2 = 'rxyz_int_atoms_2' 
  character(len=*), parameter :: ASTRUCT_ATT_RXYZ_INT_3 = 'rxyz_int_atoms_3' 

  character(len=*), parameter :: GOUT_ENERGY = 'energy (Ha)' 
  character(len=*), parameter :: GOUT_FORCES = 'forces (Ha/Bohr)' 
  character(len=*), parameter :: FORMAT_KEY = 'format' 
  character(len=*), parameter :: FORMAT_YAML = 'yaml' 
  character(len=*), parameter :: RADII_KEY = 'Radii of active regions (AU)' 
  character(len=*), parameter :: LPSP_KEY = 'Local Pseudo Potential (HGH convention)' 
  character(len=*), parameter :: NLPSP_KEY = 'NonLocal PSP Parameters'
  character(len=*), parameter :: PSPXC_KEY = 'Pseudopotential XC'
  character(len=*), parameter :: PSP_TYPE = 'Pseudopotential type'
  character(len=*), parameter :: COARSE = 'Coarse'
  character(len=*), parameter :: COARSE_PSP = 'Coarse PSP'
  character(len=*), parameter :: FINE = 'Fine'
  character(len=*), parameter :: SOURCE_KEY = 'Source'
  character(len=*), parameter :: ATOMIC_NUMBER = 'Atomic number'
  character(len=*), parameter :: ELECTRON_NUMBER = 'No. of Electrons'
  character(len=*), parameter :: POSINP_SOURCE = 'source'

end module public_keys

!>module identifying constants that have to be used as enumerators
!! they can be used to define f_enumerator types or directly as integers
module public_enums
  use f_enums
  implicit none
  
  public
  
  type(f_enumerator), parameter :: LENNARD_JONES_RUN_MODE      =f_enumerator('LENNARD_JONES_RUN_MODE',-1000,null())
  type(f_enumerator), parameter :: LENOSKY_SI_CLUSTERS_RUN_MODE=f_enumerator('LENOSKY_SI_CLUSTERS_RUN_MODE',-999,null())
  type(f_enumerator), parameter :: LENOSKY_SI_BULK_RUN_MODE    =f_enumerator('LENOSKY_SI_BULK_RUN_MODE',-998,null())
  type(f_enumerator), parameter :: AMBER_RUN_MODE              =f_enumerator('AMBER_RUN_MODE',-997,null())
  type(f_enumerator), parameter :: MORSE_BULK_RUN_MODE         =f_enumerator('MORSE_BULK_RUN_MODE',-996,null())
  type(f_enumerator), parameter :: MORSE_SLAB_RUN_MODE         =f_enumerator('MORSE_SLAB_RUN_MODE',-995,null())
  type(f_enumerator), parameter :: QM_RUN_MODE                 =f_enumerator('QM_RUN_MODE',-994,null())
  type(f_enumerator), parameter :: TERSOFF_RUN_MODE         =f_enumerator('TERSOFF_RUN_MODE',-993,null())
  type(f_enumerator), parameter :: BMHTF_RUN_MODE         =f_enumerator('BMHTF_RUN_MODE',-992,null())
  
end module public_enums





