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
  
  character(len = *), parameter :: PY_HOOKS = "py_hooks"
  character(len = *), parameter :: PLUGINS = "plugins"
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
  character(len = *), parameter :: GNRM_CV_VIRT = "gnrm_cv_virt"
  character(len = *), parameter :: GNRM_IG = "gnrm_ig"
  character(len = *), parameter :: NIT_IG = "nit_ig"
  character(len = *), parameter :: ITERMAX = "itermax",ITERMIN = "itermin", NREPMAX = "nrepmax"
  character(len = *), parameter :: ITERMAX_VIRT = "itermax_virt"
  character(len = *), parameter :: NCONG = "ncong", IDSX = "idsx"
  character(len = *), parameter :: DISPERSION = "dispersion"
  character(len = *), parameter :: INPUTPSIID = "inputpsiid"
  character(len = *), parameter :: OUTPUT_WF = "output_wf"
  character(len = *), parameter :: OUTPUT_MAT = "output_mat"
  character(len = *), parameter :: OUTPUT_COEFF = "output_coeff"
  character(len = *), parameter :: OUTPUT_DENSPOT = "output_denspot"
  character(len = *), parameter :: OUTPUT_FRAGMENTS = "output_fragments"
  character(len = *), parameter :: KERNEL_RESTART_MODE = "kernel_restart_mode"
  character(len = *), parameter :: KERNEL_RESTART_NOISE = "kernel_restart_noise"
  character(len = *), parameter :: FRAG_NUM_NEIGHBOURS = "frag_num_neighbours"
  character(len = *), parameter :: FRAG_NEIGHBOUR_CUTOFF = "frag_neighbour_cutoff"
  character(len = *), parameter :: RBUF = "rbuf"
  character(len = *), parameter :: NCONGT = "ncongt"
  character(len = *), parameter :: NORBV = "norbv", NVIRT = "nvirt"
  character(len = *), parameter :: NPLOT = "nplot"
  character(len = *), parameter :: DISABLE_SYM = "disablesym"
  character(len = *), parameter :: SOLVENT = "solvent"
  character(len = *), parameter :: EXTERNAL_POTENTIAL = "external_potential"
  character(len = *), parameter :: CHARGE_MULTIPOLES = "charge_multipoles"
  character(len = *), parameter :: CALCULATE_STRTEN = "calculate_strten"
  character(len = *), parameter :: OCCUPANCY_CONTROL = "occupancy_control"
  character(len = *), parameter :: OCCUPANCY_CONTROL_ITERMAX= "itermax_occ_ctrl"
  character(len = *), parameter :: RESET_DIIS_HISTORY = "reset_DIIS_history"

  character(len = *), parameter :: PSOLVER = "psolver"

  character(len = *), parameter :: OUTPUT_VARIABLES = "output"
  character(len = *), parameter :: ATOMIC_DENSITY_MATRIX = "atomic_density_matrix" 

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
  character(len = *), parameter :: SOCKINET = "sockinet"
  character(len = *), parameter :: SOCKPORT = "sockport"
  character(len = *), parameter :: SOCKHOST = "sockhost"
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

!MD keywords
  character(len = *), parameter :: MD_VARIABLES = "md"
  character(len = *), parameter :: MDSTEPS = "mdsteps"
  character(len = *), parameter :: PRINT_FREQUENCY = "print_frequency"
  character(len = *), parameter :: TEMPERATURE = "temperature"
  character(len = *), parameter :: TIMESTEP = "timestep"
  character(len = *), parameter :: NO_TRANSLATION = "no_translation"
  character(len = *), parameter :: THERMOSTAT = "thermostat"
  character(len = *), parameter :: NOSE_CHAIN_LENGTH = "nose_chain_length"
  character(len = *), parameter :: NOSE_MTS_SIZE = "nose_mts_size"
  character(len = *), parameter :: NOSE_YOSHIDA_FACTOR = "nose_yoshida_factor"
  character(len = *), parameter :: NOSE_FREQUENCY = "nose_frequency"
  character(len = *), parameter :: WAVEFUNCTION_EXTRAPOLATION="wavefunction_extrapolation"
  character(len = *), parameter :: RESTART_POS="restart_pos"
  character(len = *), parameter :: RESTART_VEL="restart_vel"
  character(len = *), parameter :: RESTART_NOSE="restart_nose"

  !mode parameter keywords
  character(len = *), parameter :: ADD_COULOMB_FORCE_KEY = "add_coulomb_force"
  character(len = *), parameter :: MM_PARAMSET = "mm_paramset" !for hard-coded parameter sets
  character(len = *), parameter :: MM_PARAMFILE = "mm_paramfile" !for parameter sets given by file
  character(len = *), parameter :: SW_EQFACTOR = "sw_eqfactor"
  character(len = *), parameter :: SECTIONS = "sections"
  character(len = *), parameter :: SECTION_BUFFER = "section_buffer"
  character(len = *), parameter :: SECTION_PASSIVATION = "section_passivation"
  character(len = *), parameter :: NAB_OPTIONS = "nab_options"

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
  character(len = *), parameter :: DECOMPOSE_PERTURBATION = 'decompose_perturbation'

  character(len = *), parameter :: PERF_VARIABLES = "perf"
  character(len = *), parameter :: DEBUG = "debug"
  character(len = *), parameter :: PROFILING_DEPTH = "profiling_depth"
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
!!$  character(len = *), parameter :: PSOLVER_GROUPSIZE = "psolver_groupsize"
!!$  character(len = *), parameter :: PSOLVER_ACCEL = "psolver_accel"
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
  character(len = *), parameter :: CALCULATE_GAP = "calculate_gap"
  character(len = *), parameter :: LOEWDIN_CHARGE_ANALYSIS = "loewdin_charge_analysis"
  character(len = *), parameter :: COEFF_WEIGHT_ANALYSIS = "coeff_weight_analysis"
  character(len = *), parameter :: CHECK_MATRIX_COMPRESSION = "check_matrix_compression"
  character(len = *), parameter :: CORRECTION_CO_CONTRA = "correction_co_contra"
  character(len = *), parameter :: GPS_METHOD = "gps_method"
  character(len = *), parameter :: FOE_GAP = "foe_gap"
  character(len = *), parameter :: SUPPORT_FUNCTION_MULTIPOLES = "support_function_multipoles"
  character(len = *), parameter :: PLOT_MPPOT_AXES = "plot_mppot_axes"
  character(len = *), parameter :: PLOT_POT_AXES = "plot_pot_axes"
  character(len = *), parameter :: PLOT_LOCREG_GRIDS = "plot_locreg_grids"
  character(len = *), parameter :: CALCULATE_FOE_EIGENVALUES = "calculate_FOE_eigenvalues"
  character(len = *), parameter :: PRECISION_FOE_EIGENVALUES = "precision_FOE_eigenvalues"

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
  character(len=*), parameter :: CALC_QUADRUPOLE ='calc_quadrupole'
  character(len=*), parameter :: CDFT_LAG_MULT_INIT='cdft_lag_mult_init'
  character(len=*), parameter :: CDFT_CONV_CRIT  ='cdft_conv_crit'
  character(len=*), parameter :: SUBSPACE_DIAG   ='subspace_diag'
  character(len=*), parameter :: ALPHA_DIIS      ='alpha_diis'
  character(len=*), parameter :: ALPHA_SD        ='alpha_sd'
  character(len=*), parameter :: ALPHA_SD_COEFF  ='alpha_sd_coeff'
  character(len=*), parameter :: ALPHA_FIT_COEFF ='alpha_fit_coeff'
  character(len=*), parameter :: NSTEP_PREC      ='nstep_prec'
  character(len=*), parameter :: EVAL_RANGE_FOE  ='eval_range_foe'
  character(len=*), parameter :: FSCALE_FOE      ='fscale_foe'
  character(len=*), parameter :: COEFF_SCALING_FACTOR='coeff_scaling_factor'
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
  character(len=*), parameter :: ORTHOGONALIZE_AO = 'orthogonalize_ao'
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
  character(len=*), parameter :: ADJUST_KERNEL_THRESHOLD='adjust_kernel_threshold'
  character(len=*), parameter :: WF_EXTENT_ANALYSIS='wf_extent_analysis'
  character(len=*), parameter :: CALCULATE_ONSITE_OVERLAP='calculate_onsite_overlap'
  character(len=*), parameter :: DELTA_PNRM='delta_pnrm'
  character(len=*), parameter :: PEXSI_NPOLES='pexsi_npoles'
  character(len=*), parameter :: PEXSI_MUMIN='pexsi_mumin'
  character(len=*), parameter :: PEXSI_MUMAX='pexsi_mumax'
  character(len=*), parameter :: PEXSI_MU='pexsi_mu'
  character(len=*), parameter :: PEXSI_TEMPERATURE='pexsi_temperature'
  character(len=*), parameter :: PEXSI_TOL_CHARGE='pexsi_tol_charge'

  !> Parameters to avoid typos in dictionary keys
  character(len=*), parameter :: ASTRUCT_UNITS = 'units'
  character(len=*), parameter :: ASTRUCT_CELL = 'cell'
  character(len=*), parameter :: ASTRUCT_POSITIONS = 'positions'
  character(len=*), parameter :: ASTRUCT_PROPERTIES = 'properties'
  character(len=*), parameter :: ASTRUCT_ATT_FROZEN = 'Frozen'
  character(len=*), parameter :: ASTRUCT_ATT_IGSPIN = 'IGSpin'
  character(len=*), parameter :: ASTRUCT_ATT_IGCHRG = 'IGChg'
  character(len=*), parameter :: ASTRUCT_ATT_CAVRAD = 'rcav' !< custom radius of the cavity
  character(len=*), parameter :: ASTRUCT_ATT_IXYZ_1 = 'int_ref_atoms_1'
  character(len=*), parameter :: ASTRUCT_ATT_IXYZ_2 = 'int_ref_atoms_2'
  character(len=*), parameter :: ASTRUCT_ATT_IXYZ_3 = 'int_ref_atoms_3'
  character(len=*), parameter :: ASTRUCT_ATT_MODE = 'mode'
  character(len=*), parameter :: ASTRUCT_ATT_RXYZ_INT_1 = 'rxyz_int_atoms_1'
  character(len=*), parameter :: ASTRUCT_ATT_RXYZ_INT_2 = 'rxyz_int_atoms_2'
  character(len=*), parameter :: ASTRUCT_ATT_RXYZ_INT_3 = 'rxyz_int_atoms_3'
  character(len=*), parameter :: ASTRUCT_ATT_ORIG_ID = 'fromNode'

  character(len=*), parameter :: GOUT_ENERGY = 'energy (Ha)'
  character(len=*), parameter :: GOUT_FORCES = 'forces (Ha/Bohr)'
  character(len=*), parameter :: FORMAT_KEY = 'format'
  character(len=*), parameter :: FORMAT_YAML = 'yaml'
  character(len=*), parameter :: RADII_KEY = 'Radii of active regions (AU)'
  character(len=*), parameter :: LPSP_KEY = 'Local Pseudo Potential (HGH convention)'
  character(len=*), parameter :: NLPSP_KEY = 'NonLocal PSP Parameters'
  character(len=*), parameter :: NLCC_KEY = 'Non Linear Core Correction term'
  character(len=*), parameter :: PSPXC_KEY = 'Pseudopotential XC'
  character(len=*), parameter :: PSP_TYPE = 'Pseudopotential type'
  character(len=*), parameter :: COARSE = 'Coarse'
  character(len=*), parameter :: COARSE_PSP = 'Coarse PSP'
  character(len=*), parameter :: FINE = 'Fine'
  character(len=*), parameter :: SOURCE_KEY = 'Source'
  character(len=*), parameter :: COEFF_KEY = 'Coefficients (c1 .. c4)'
  character(len=*), parameter :: ATOMIC_NUMBER = 'Atomic number'
  character(len=*), parameter :: ELECTRON_NUMBER = 'No. of Electrons'
  character(len=*), parameter :: POSINP_SOURCE = 'source'
end module public_keys

!>module identifying constants that have to be used as enumerators
!! they can be used to define f_enumerator types or directly as integers
module public_enums
  use f_enums
  implicit none

  private !as we use at least one module, private become compulsory

  !> Error codes, to be documented little by little
  integer, parameter, public :: BIGDFT_SUCCESS        = 0   !< No errors
  integer, parameter, public :: BIGDFT_UNINITIALIZED  = -10 !< The quantities we want to access seem not yet defined
  integer, parameter, public :: BIGDFT_INCONSISTENCY  = -11 !< Some of the quantities is not correct
  integer, parameter, public :: BIGDFT_INVALID        = -12 !< Invalid entry

  !> How to do the partition of the support functions (linear scaling version)
  integer,parameter,public :: LINEAR_PARTITION_SIMPLE = 61
  integer,parameter,public :: LINEAR_PARTITION_OPTIMAL = 62
  integer,parameter,public :: LINEAR_PARTITION_NONE = 63

  !> Flags for rhov status
  integer, parameter, public :: EMPTY              = -1980
  integer, parameter, public :: ELECTRONIC_DENSITY = -1979
  integer, parameter, public :: CHARGE_DENSITY     = -1978
  integer, parameter, public :: KS_POTENTIAL       = -1977
  integer, parameter, public :: HARTREE_POTENTIAL  = -1976

  !> Flags for the restart (linear scaling only)
  integer, parameter, public :: LINEAR_LOWACCURACY  = 101 !< Low accuracy after restart
  integer, parameter, public :: LINEAR_HIGHACCURACY = 102 !< High accuracy after restart

  !> Flags for optimization loop id
  integer, parameter, public :: OPTLOOP_HAMILTONIAN   = 0
  integer, parameter, public :: OPTLOOP_SUBSPACE      = 1
  integer, parameter, public :: OPTLOOP_WAVEFUNCTIONS = 2
  integer, parameter, public :: OPTLOOP_N_LOOPS       = 3

  !> Constants to determine between cubic version and linear version
  integer, parameter, public :: CUBIC_VERSION =  0
  integer, parameter, public :: LINEAR_VERSION = 100
  integer, parameter  :: GAUSSIAN_ENHANCEMENT = 102

  integer, parameter :: NONE=0
  integer, parameter :: TEXT=1
  integer, parameter :: BINARY=2
  integer, parameter :: ETSF=3
  integer, parameter :: CUBE=22
  integer, parameter :: MPI_NATIVE=4 !<native (i.e. non-portable) MPI format


  !> Output wf parameters.
  integer, parameter, public :: WF_FORMAT_NONE       = NONE
  integer, parameter, public :: WF_FORMAT_PLAIN      = TEXT
  integer, parameter, public :: WF_FORMAT_BINARY     = BINARY
  integer, parameter, public :: WF_FORMAT_ETSF       = ETSF

  !> Output matrix parameters.
  integer, parameter, public :: MATRIX_FORMAT_NONE       = NONE
  integer, parameter, public :: MATRIX_FORMAT_PLAIN      = TEXT
  integer, parameter, public :: MATRIX_FORMAT_BINARY     = BINARY
  integer, parameter, public :: MATRIX_FORMAT_ETSF       = ETSF
  integer, parameter, public :: MATRIX_FORMAT_MPI_NATIVE = MPI_NATIVE

  !> Output grid parameters.
  ! with these options we would have
  ! 0 : none
  ! 1 : density, plain
  ! 2 : denspot, plain
  !
  integer, parameter, public :: OUTPUT_DENSPOT_NONE    = 0
  integer, parameter, public :: OUTPUT_DENSPOT_DENSITY = 1
  integer, parameter, public :: OUTPUT_DENSPOT_DENSPOT = 2

  integer, parameter, public :: OUTPUT_DENSPOT_FORMAT_TEXT = TEXT
  integer, parameter, public :: OUTPUT_DENSPOT_FORMAT_ETSF = ETSF
  integer, parameter, public :: OUTPUT_DENSPOT_FORMAT_CUBE = CUBE


  !>empty enumerator
  type(f_enumerator), public :: ENUM_EMPTY =f_enumerator('NONE',0,null())

  !> enumerators defining the psi policy
  type(f_enumerator), public :: ENUM_SCRATCH =f_enumerator('SCRATCH',10000,null())
  type(f_enumerator), public :: ENUM_MEMORY =f_enumerator('MEMORY',10001,null())
  type(f_enumerator), public :: ENUM_FILE =f_enumerator('FILE',10002,null())

  !> enumerators defining the mode
  type(f_enumerator), public :: ENUM_CUBIC =f_enumerator('CUBIC',CUBIC_VERSION,null())
  type(f_enumerator), public :: ENUM_LINEAR =f_enumerator('LINEAR',LINEAR_VERSION,null())
  type(f_enumerator), public :: ENUM_GAUSSIAN =f_enumerator('GAUSSIAN',20002,null())

  !> enumerators defining the i/o format
  type(f_enumerator), public :: ENUM_TEXT =f_enumerator('TEXT',TEXT,null())
  type(f_enumerator), public :: ENUM_ETSF =f_enumerator('ETSF',ETSF,null())
  type(f_enumerator), public :: ENUM_CUBE =f_enumerator('CUBE',CUBE,null())
  type(f_enumerator), public :: ENUM_BINARY =f_enumerator('BINARY',BINARY,null())


  !> Input wf parameters. @relates module_types::input_variables::inputpsiid @relates inputpsiid
  !! used to define the inputpsiid enumerator and the corresponding attributes
  integer, parameter, public :: INPUT_PSI_EMPTY        = -1000  !< Input PSI to 0
  integer, parameter, public :: INPUT_PSI_RANDOM       = -2     !< Input Random PSI
  integer, parameter, public :: INPUT_PSI_CP2K         = -1     !< Input PSI coming from cp2k
  integer, parameter, public :: INPUT_PSI_LCAO         = 0      !< Input PSI coming from Localised ATomic Orbtials
  integer, parameter, public :: INPUT_PSI_MEMORY_WVL   = 1      !< Input PSI from memory
  integer, parameter, public :: INPUT_PSI_DISK_WVL     = 2      !< Input PSI from disk (wavelet coefficients)
  integer, parameter, public :: INPUT_PSI_DISK_PW      = 3
  integer, parameter, public :: INPUT_PSI_LCAO_GAUSS   = 10
  integer, parameter, public :: INPUT_PSI_MEMORY_GAUSS = 11
  integer, parameter, public :: INPUT_PSI_DISK_GAUSS   = 12
  integer, parameter, public :: INPUT_PSI_LINEAR_AO    = 100    !< Input PSI for linear from Atomic Orbital
  integer, parameter, public :: INPUT_PSI_MEMORY_LINEAR= 101    !< Input PSI for linear in memory
  integer, parameter, public :: INPUT_PSI_DISK_LINEAR  = 102    !< Input PSI for linear from disk

  !> Variables for kernel restart in linear
  integer, parameter, public :: LIN_RESTART_KERNEL       = 0  !< Use kernel
  integer, parameter, public :: LIN_RESTART_COEFF        = 1  !< Use coefficients
  integer, parameter, public :: LIN_RESTART_RANDOM       = 2  !< Use randomized kernel
  integer, parameter, public :: LIN_RESTART_DIAG_KERNEL  = 3  !< Use diagonal kernel scaled by number of electrons
  integer, parameter, public :: LIN_RESTART_TMB_WEIGHT   = 4  !< Use support function weights

  !> Variables for fragment output
  integer, parameter, public :: OUTPUT_FRAGMENTS_AND_FULL = 0  !< Output both
  integer, parameter, public :: OUTPUT_FRAGMENTS_ONLY     = 1  !< Output fragments only
  integer, parameter, public :: OUTPUT_FULL_ONLY          = 2  !< Output full system onlu

  !> Source of the radii coefficients
  integer, parameter, public :: RADII_SOURCE_HARD_CODED = 1
  integer, parameter, public :: RADII_SOURCE_FILE = 2
  integer, parameter, public :: RADII_SOURCE_USER = 3
  integer, parameter, public :: RADII_SOURCE_UNKNOWN = 4
  character(len=*), dimension(4), parameter, public :: RADII_SOURCE = (/ &
       "Hard-Coded  ", &
       "PSP File    ", &
       "User defined", &
       "Unknown     " /)



  !> SCF mixing parameters (@todo mixing parameters to be added).
  integer, parameter, public :: SCF_KIND_GENERALIZED_DIRMIN = -1
  integer, parameter, public :: SCF_KIND_DIRECT_MINIMIZATION = 0

  !> Function to determine the occupation numbers
  integer, parameter, public :: SMEARING_DIST_ERF   = 1  !< Tends to 0 and 1 faster \f$1/2\left[1-erf\left(\frac{E-\mu}{\delta E}\right)\right]\f$
  integer, parameter, public :: SMEARING_DIST_FERMI = 2  !< Normal Fermi distribution i.e.\f$\frac{1}{1+e^{E-\mu}/k_BT}\f$
  integer, parameter, public :: SMEARING_DIST_COLD1 = 3  !< Marzari's cold smearing with a=-.5634 (bumb minimization)
  integer, parameter, public :: SMEARING_DIST_COLD2 = 4  !< Marzari's cold smearing with a=-.8165 (monotonic tail)
  integer, parameter, public :: SMEARING_DIST_METPX = 5  !< Methfessel and Paxton (same as COLD with a=0)

  !> Target function for the optimization of the basis functions (linear scaling version)
  integer, parameter, public :: TARGET_FUNCTION_IS_TRACE=0
  integer, parameter, public :: TARGET_FUNCTION_IS_ENERGY=1
  integer, parameter, public :: TARGET_FUNCTION_IS_HYBRID=2
  !!integer, parameter, public :: DECREASE_LINEAR=0
  !!integer, parameter, public :: DECREASE_ABRUPT=1
  !!integer, parameter, public :: COMMUNICATION_COLLECTIVE=0
  !!integer, parameter, public :: COMMUNICATION_P2P=1
  integer, parameter, public :: LINEAR_DIRECT_MINIMIZATION=100
  integer, parameter, public :: LINEAR_MIXDENS_SIMPLE=101
  integer, parameter, public :: LINEAR_MIXPOT_SIMPLE=102
  integer, parameter, public :: LINEAR_FOE=103
  integer, parameter, public :: LINEAR_PEXSI=104
  integer, parameter, public :: KERNELMODE_DIRMIN = 10
  integer, parameter, public :: KERNELMODE_DIAG = 11
  integer, parameter, public :: KERNELMODE_FOE = 12
  integer, parameter, public :: KERNELMODE_PEXSI = 13
  integer, parameter, public :: MIXINGMODE_DENS = 20
  integer, parameter, public :: MIXINGMODE_POT = 21
  integer, parameter, public :: FOE_ACCURATE = 30
  integer, parameter, public :: FOE_FAST = 31

  integer, parameter :: MIXING= 1000
  integer, parameter :: DIRMIN= 0

  integer, parameter, public :: AB7_MIXING_NONE        = 0
  integer, parameter, public :: AB7_MIXING_EIG         = 1
  integer, parameter, public :: AB7_MIXING_SIMPLE      = 2
  integer, parameter, public :: AB7_MIXING_ANDERSON    = 3
  integer, parameter, public :: AB7_MIXING_ANDERSON_2  = 4
  integer, parameter, public :: AB7_MIXING_CG_ENERGY   = 5
  integer, parameter, public :: AB7_MIXING_CG_ENERGY_2 = 6
  integer, parameter, public :: AB7_MIXING_PULAY       = 7

  integer, parameter, public :: AB7_MIXING_POTENTIAL  = 0
  integer, parameter, public :: AB7_MIXING_DENSITY    = 1

  integer, parameter, public :: AB7_MIXING_REAL_SPACE    = 1
  integer, parameter, public :: AB7_MIXING_FOURIER_SPACE = 2

  type(f_enumerator), public :: MIX_ENUM=f_enumerator('MIXING',MIXING,null())

  type(f_enumerator), public :: POT_MIX_ENUM=f_enumerator('MIXING_ON',AB7_MIXING_POTENTIAL,null())
  type(f_enumerator), public :: DEN_MIX_ENUM=f_enumerator('MIXING_ON',AB7_MIXING_DENSITY,null())

  type(f_enumerator), public :: RSPACE_MIX_ENUM=f_enumerator('MIXING_SPACE',AB7_MIXING_REAL_SPACE,null())
  type(f_enumerator), public :: FSPACE_MIX_ENUM=f_enumerator('MIXING_SPACE',AB7_MIXING_FOURIER_SPACE,null())



  !> How to update the density kernel during teh support function optimization
  integer, parameter, public :: UPDATE_BY_PURIFICATION = 0
  integer, parameter, public :: UPDATE_BY_FOE = 1
  integer, parameter, public :: UPDATE_BY_RENORMALIZATION = 2

  integer, parameter, public :: INPUT_IG_OFF  = 0
  integer, parameter, public :: INPUT_IG_LIG  = 1
  integer, parameter, public :: INPUT_IG_FULL = 2
  integer, parameter, public :: INPUT_IG_TMO  = 3

  !> Flags for the input files.
  integer, parameter, public :: INPUTS_NONE  =   0
  integer, parameter, public :: INPUTS_DFT   =   1
  integer, parameter, public :: INPUTS_GEOPT =   2
  integer, parameter, public :: INPUTS_PERF  =   4
  integer, parameter, public :: INPUTS_KPT   =   8
  integer, parameter, public :: INPUTS_MIX   =  16
  integer, parameter, public :: INPUTS_TDDFT =  32
  integer, parameter, public :: INPUTS_SIC   =  64
  integer, parameter, public :: INPUTS_FREQ  = 128
  integer, parameter, public :: INPUTS_LIN   = 256
  integer, parameter, public :: INPUTS_FRAG  = 512

  !> Type of pseudopotential
  integer, parameter, public :: PSPCODE_UNINITIALIZED = 1
  integer, parameter, public :: PSPCODE_GTH = 2
  integer, parameter, public :: PSPCODE_HGH = 3
  integer, parameter, public :: PSPCODE_PAW = 7
  integer, parameter, public :: PSPCODE_HGH_K = 10
  integer, parameter, public :: PSPCODE_HGH_K_NLCC = 12

  !> Output for run modes
  type(f_enumerator), public :: RUN_MODE_CREATE_DOCUMENT = &
       & f_enumerator('RUN_MODE_CREATE_DOCUMENT',1,null())
  type(f_enumerator), parameter, public :: TDPOT_RUN_MODE      =f_enumerator('TDPOT_RUN_MODE',-989,null())

  !run modes
  type(f_enumerator), parameter, public :: LENNARD_JONES_RUN_MODE      = &
       & f_enumerator('LENNARD_JONES_RUN_MODE',-1000,null())
  type(f_enumerator), parameter, public :: LENOSKY_SI_CLUSTERS_RUN_MODE= &
       & f_enumerator('LENOSKY_SI_CLUSTERS_RUN_MODE',-999,null())
  type(f_enumerator), parameter, public :: LENOSKY_SI_BULK_RUN_MODE    = &
       & f_enumerator('LENOSKY_SI_BULK_RUN_MODE',-998,null())
  type(f_enumerator), parameter, public :: AMBER_RUN_MODE              = &
       & f_enumerator('AMBER_RUN_MODE',-997,null())
  type(f_enumerator), parameter, public :: MORSE_BULK_RUN_MODE         = &
       & f_enumerator('MORSE_BULK_RUN_MODE',-996,null())
  type(f_enumerator), parameter, public :: MORSE_SLAB_RUN_MODE         = &
       & f_enumerator('MORSE_SLAB_RUN_MODE',-995,null())
  type(f_enumerator), parameter, public :: QM_RUN_MODE                 = &
       & f_enumerator('QM_RUN_MODE',-994,null())
  type(f_enumerator), parameter, public :: TERSOFF_RUN_MODE            = &
       & f_enumerator('TERSOFF_RUN_MODE',-993,null())
  type(f_enumerator), parameter, public :: BMHTF_RUN_MODE              = &
       & f_enumerator('BMHTF_RUN_MODE',-992,null())
  type(f_enumerator), parameter, public :: CP2K_RUN_MODE               = &
       & f_enumerator('CP2K_RUN_MODE',-991,null())
  type(f_enumerator), parameter, public :: DFTBP_RUN_MODE              = &
       & f_enumerator('DFTBP_RUN_MODE',-990,null())
  type(f_enumerator), parameter, public :: MULTI_RUN_MODE              = &
       & f_enumerator('MULTI_RUN_MODE',-989,null())
  type(f_enumerator), parameter, public :: SW_RUN_MODE                 = &
       & f_enumerator('SW_RUN_MODE',-988,null())
  type(f_enumerator), parameter, public :: BAZANT_RUN_MODE             = &
       & f_enumerator('BAZANT_RUN_MODE',-987,null()) 

end module public_enums
