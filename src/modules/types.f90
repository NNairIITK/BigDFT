!> @file
!!  Define the fortran types
!! @author
!!    Copyright (C) 2008-2013 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Module which contains the Fortran data structures
!! and the routines of allocations and de-allocations
module module_types

  use m_ab7_mixing, only : ab7_mixing_object
  use module_base, only : gp,wp,dp,tp,uninitialized,mpi_environment,mpi_environment_null,&
       bigdft_mpi,ndebug,memocc!,vcopy
  use gaussians, only: gaussian_basis
  use Poisson_Solver, only: coulomb_operator
  use dictionaries, only: dictionary
  use locregs
  use psp_projectors
  use module_atoms, only: atoms_data,symmetry_data,atomic_structure
  use communications_base, only: comms_linear, comms_cubic
  use sparsematrix_base, only: sparse_matrix

  implicit none

  !> Constants to determine between cubic version and linear version
  integer, parameter :: CUBIC_VERSION =  0
  integer, parameter :: LINEAR_VERSION = 100

  !> Error codes, to be documented little by little
  integer, parameter :: BIGDFT_SUCCESS        = 0   !< No errors
  integer, parameter :: BIGDFT_UNINITIALIZED  = -10 !< The quantities we want to access seem not yet defined
  integer, parameter :: BIGDFT_INCONSISTENCY  = -11 !< Some of the quantities is not correct
  integer, parameter :: BIGDFT_INVALID        = -12 !< Invalid entry
  integer :: BIGDFT_RUNTIME_ERROR                   !< error during runtime
  integer :: BIGDFT_MPI_ERROR                       !< see error definitions below
  integer :: BIGDFT_LINALG_ERROR                    !< to be moved to linalg wrappers
  integer :: BIGDFT_INPUT_VARIABLES_ERROR           !< problems in parsing or in consistency of input variables

  !> Input wf parameters.
  integer, parameter :: INPUT_PSI_EMPTY        = -1000  !< Input PSI to 0
  integer, parameter :: INPUT_PSI_RANDOM       = -2     !< Input Random PSI
  integer, parameter :: INPUT_PSI_CP2K         = -1     !< Input PSI coming from cp2k
  integer, parameter :: INPUT_PSI_LCAO         = 0      !< Input PSI coming from Localised ATomic Orbtials
  integer, parameter :: INPUT_PSI_MEMORY_WVL   = 1
  integer, parameter :: INPUT_PSI_DISK_WVL     = 2
  integer, parameter :: INPUT_PSI_LCAO_GAUSS   = 10
  integer, parameter :: INPUT_PSI_MEMORY_GAUSS = 11
  integer, parameter :: INPUT_PSI_DISK_GAUSS   = 12
  integer, parameter :: INPUT_PSI_LINEAR_AO    = 100
  integer, parameter :: INPUT_PSI_MEMORY_LINEAR= 101
  integer, parameter :: INPUT_PSI_DISK_LINEAR  = 102

  !> All possible values of input psi (determnitation of the input guess)
  integer, dimension(12), parameter :: input_psi_values = &
       (/ INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_CP2K, &
       INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, INPUT_PSI_DISK_WVL, &
       INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS, &
       INPUT_PSI_LINEAR_AO, INPUT_PSI_DISK_LINEAR, INPUT_PSI_MEMORY_LINEAR /)

  !> Output wf parameters.
  integer, parameter :: WF_FORMAT_NONE   = 0
  integer, parameter :: WF_FORMAT_PLAIN  = 1
  integer, parameter :: WF_FORMAT_BINARY = 2
  integer, parameter :: WF_FORMAT_ETSF   = 3
  integer, parameter :: WF_N_FORMAT      = 4
  character(len = 12), dimension(0:WF_N_FORMAT-1), parameter :: wf_format_names = &
       (/ "none        ", &
          "plain text  ", &
          "Fortran bin.", &
          "ETSF        " /)

  !> Output grid parameters.
  integer, parameter :: OUTPUT_DENSPOT_NONE    = 0
  integer, parameter :: OUTPUT_DENSPOT_DENSITY = 1
  integer, parameter :: OUTPUT_DENSPOT_DENSPOT = 2
  character(len = 12), dimension(0:2), parameter :: OUTPUT_DENSPOT_names = &
       (/ "none        ", &
          "density     ", &
          "dens. + pot." /)
  integer, parameter :: OUTPUT_DENSPOT_FORMAT_TEXT = 0
  integer, parameter :: OUTPUT_DENSPOT_FORMAT_ETSF = 1
  integer, parameter :: OUTPUT_DENSPOT_FORMAT_CUBE = 2
  character(len = 4), dimension(0:2), parameter :: OUTPUT_DENSPOT_format_names = &
       (/ "text", &
          "ETSF", &
          "cube" /)

  !> SCF mixing parameters. (mixing parameters to be added)
  integer, parameter :: SCF_KIND_GENERALIZED_DIRMIN = -1
  integer, parameter :: SCF_KIND_DIRECT_MINIMIZATION = 0

  !> Function to determine the occupation numbers
  integer, parameter :: SMEARING_DIST_ERF   = 1  !< tends to 0 and 1 faster \f$1/2\left[1-erf\left(\frac{E-\mu}{\delta E}\right)\right]\f$
  integer, parameter :: SMEARING_DIST_FERMI = 2  !< Normal Fermi distribution i.e.\f$\frac{1}{1+e^{E-\mu}/k_BT}\f$
  integer, parameter :: SMEARING_DIST_COLD1 = 3  !< Marzari's cold smearing with a=-.5634 (bumb minimization)
  integer, parameter :: SMEARING_DIST_COLD2 = 4  !< Marzari's cold smearing with a=-.8165 (monotonic tail)
  integer, parameter :: SMEARING_DIST_METPX = 5  !< Methfessel and Paxton (same as COLD with a=0)
  character(len = 11), dimension(5), parameter :: smearing_names = &
       (/ "Error func.", &
          "Fermi      ", &
          "Cold (bumb)", &
          "Cold (mono)", &
          "Meth.-Pax. " /) !< Name of the smearing methods 

  !> Target function for the optimization of the basis functions (linear scaling version)
  integer, parameter :: TARGET_FUNCTION_IS_TRACE=0
  integer, parameter :: TARGET_FUNCTION_IS_ENERGY=1
  integer, parameter :: TARGET_FUNCTION_IS_HYBRID=2
  !!integer, parameter :: DECREASE_LINEAR=0
  !!integer, parameter :: DECREASE_ABRUPT=1
  !!integer, parameter :: COMMUNICATION_COLLECTIVE=0
  !!integer, parameter :: COMMUNICATION_P2P=1
  integer, parameter :: LINEAR_DIRECT_MINIMIZATION=100
  integer, parameter :: LINEAR_MIXDENS_SIMPLE=101
  integer, parameter :: LINEAR_MIXPOT_SIMPLE=102
  integer, parameter :: LINEAR_FOE=103

  !> How to update the density kernel during teh support function optimization
  integer, parameter :: UPDATE_BY_PURIFICATION = 0
  integer, parameter :: UPDATE_BY_FOE = 1
  
  !> Type used for the orthogonalisation parameters
  type, public :: orthon_data
     !> directDiag decides which input guess is chosen:
     !!   if .true. -> as usual direct diagonalization of the Hamiltonian with dsyev (suitable for small systems)
     !!   if .false. -> iterative diagonalization (suitable for large systems)
     logical :: directDiag
     !> norbpInguess indicates how many orbitals shall be treated by each process during the input guess
     !! if directDiag=.false.
     integer :: norbpInguess
     !> You have to choose two numbers for the block size, bsLow and bsUp:
     !!   if bsLow<bsUp, then the program will choose an appropriate block size in between these two numbers
     !!   if bsLow==bsUp, then the program will take exactly this blocksize
     integer :: bsLow
     !> Block size up value (see bsLow)
     integer :: bsUp
     !> The variable methOrtho indicates which orthonormalization procedure is used:
     !!   methOrtho==0 -> Gram-Schmidt with Cholesky decomposition
     !!   methOrtho==1 -> combined block wise classical Gram-Schmidt and Cholesky
     !!   methOrtho==2 -> Loewdin
     integer :: methOrtho
     real(gp) :: iguessTol           !< Gives the tolerance to which the input guess will converged (maximal residue of all orbitals).
     integer :: methTransformOverlap, nItOrtho, blocksize_pdsyev, blocksize_pdgemm, nproc_pdsyev
  end type orthon_data

  type, public :: SIC_data
     character(len=4) :: approach !< Approach for the Self-Interaction-Correction (PZ, NK)
     integer :: ixc               !< Base for the SIC correction
     real(gp) :: alpha            !< Downscaling coefficient
     real(gp) :: fref             !< Reference occupation (for alphaNK case)
  end type SIC_data

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

  !> Contains all parameters related to the linear scaling version.
  type,public:: linearInputParameters 
    integer :: DIIS_hist_lowaccur, DIIS_hist_highaccur, nItPrecond
    integer :: nItSCCWhenOptimizing, nItBasis_lowaccuracy, nItBasis_highaccuracy
    integer :: mixHist_lowaccuracy, mixHist_highaccuracy
    integer :: dmin_hist_lowaccuracy, dmin_hist_highaccuracy
    integer :: methTransformOverlap, blocksize_pdgemm, blocksize_pdsyev
    integer :: correctionOrthoconstraint, nproc_pdsyev, nproc_pdgemm
    integer :: nit_lowaccuracy, nit_highaccuracy, nItdmin_lowaccuracy, nItdmin_highaccuracy
    integer :: nItSCCWhenFixed_lowaccuracy, nItSCCWhenFixed_highaccuracy
    real(kind=8) :: convCrit_lowaccuracy, convCrit_highaccuracy, alphaSD, alphaDIIS, evlow, evhigh, ef_interpol_chargediff
    real(kind=8) :: alpha_mix_lowaccuracy, alpha_mix_highaccuracy, reduce_confinement_factor, ef_interpol_det
    integer :: plotBasisFunctions
    real(kind=8) ::  fscale, deltaenergy_multiplier_TMBexit, deltaenergy_multiplier_TMBfix
    real(kind=8) :: lowaccuracy_conv_crit, convCritMix_lowaccuracy, convCritMix_highaccuracy
    real(kind=8) :: highaccuracy_conv_crit, support_functions_converged, alphaSD_coeff
    real(kind=8) :: convCritDmin_lowaccuracy, convCritDmin_highaccuracy
    real(kind=8), dimension(:), pointer :: locrad, locrad_lowaccuracy, locrad_highaccuracy, kernel_cutoff_FOE
    real(kind=8), dimension(:,:), pointer :: locrad_type
    real(kind=8), dimension(:), pointer :: potentialPrefac_lowaccuracy, potentialPrefac_highaccuracy, potentialPrefac_ao
    real(kind=8), dimension(:),pointer :: kernel_cutoff, locrad_kernel
    real(kind=8) :: early_stop, gnrm_dynamic
    integer, dimension(:), pointer :: norbsPerType
    integer :: scf_mode, nlevel_accuracy
    logical :: calc_dipole, pulay_correction, mixing_after_inputguess, iterative_orthogonalization, new_pulay_correction
    logical :: fragment_calculation, calc_transfer_integrals, constrained_dft, curvefit_dmin, diag_end, diag_start
    integer :: extra_states, order_taylor
  end type linearInputParameters

  type,public:: fragmentInputParameters
    integer :: nfrag_ref, nfrag
    integer, dimension(:), pointer :: frag_index ! array matching system fragments to reference fragments
    real(kind=8), dimension(:), pointer :: charge ! array giving the charge on each fragment for constrained DFT calculations
    !integer, dimension(:,:), pointer :: frag_info !array giving number of atoms in fragment and environment for reference fragments
    character(len=100), dimension(:), pointer :: label ! array of fragment names
    character(len=100), dimension(:), pointer :: dirname ! array of fragment directories, blank if not a fragment calculation
  end type fragmentInputParameters


  integer, parameter, public :: INPUT_IG_OFF  = 0
  integer, parameter, public :: INPUT_IG_LIG  = 1
  integer, parameter, public :: INPUT_IG_FULL = 2
  integer, parameter, public :: INPUT_IG_TMO  = 3

  !> Structure controlling the nature of the accelerations (Convolutions, Poisson Solver)
  type, public :: material_acceleration
     !> variable for material acceleration
     !! values 0: traditional CPU calculation
     !!        1: CUDA acceleration with CUBLAS
     !!        2: OpenCL acceleration (with CUBLAS one day)
     integer :: iacceleration
     integer :: Psolver_igpu !< acceleration of the Poisson solver
     character(len=11) :: OCL_platform
     character(len=11) :: OCL_devices
  end type material_acceleration


  !> Structure of the variables read by input.* files (*.dft, *.geopt...)
  type, public :: input_variables
     !strings of the input files
     character(len=100) :: file_occnum,file_igpop,file_lin,file_frag
     character(len=100) :: dir_output !< Strings of the directory which contains all data output files
     character(len=100) :: run_name   !< Contains the prefix (by default input) used for input files as input.dft
     integer :: files                 !< Existing files.
     !miscellaneous variables
     logical :: gaussian_help
     integer :: itrpmax
     integer :: iscf,norbsempty,norbsuempty,norbsdempty, occopt
     integer :: last_run
     real(gp) :: frac_fluct,gnrm_sw,alphamix,Tel, alphadiis
     real(gp) :: rpnrm_cv,gnrm_startmix
     integer :: verbosity
     ! DFT basic parameters.
     integer :: ixc,ncharge,itermax,nrepmax,ncong,idsx,ncongt,inputPsiId,nspin,mpol
     integer :: norbv,nvirt,nplot
     integer :: output_denspot,dispersion,output_wf_format,output_denspot_format
     real(gp) :: hx,hy,hz,crmult,frmult,gnrm_cv,rbuf
     real(gp) :: elecfield(3)
     logical :: disableSym

     ! For absorption calculations
     integer :: iabscalc_type   !< 0 non calc, 1 cheb ,  2 lanc
     !! integer :: iat_absorber, L_absorber, N_absorber, rpower_absorber, Linit_absorber
     integer :: iat_absorber,  L_absorber
     real(gp), pointer :: Gabs_coeffs(:)
     real(gp) :: abscalc_bottomshift
     logical ::  c_absorbtion , abscalc_alterpot, abscalc_eqdiff,abscalc_S_do_cg, abscalc_Sinv_do_cg
     integer ::  potshortcut
     integer ::  nsteps
     character(len=100) :: extraOrbital
     character(len=1000) :: xabs_res_prefix
   
     ! Frequencies calculations (finite difference)
     real(gp) :: freq_alpha  !< Factor for the finite difference step (step = alpha * hgrid)
     integer :: freq_order   !< Order of the finite difference scheme
     integer :: freq_method  !< Method to calculate the frequencies

     ! kpoints related input variables
     ! generated results
     integer :: gen_nkpt
     real(gp), pointer :: gen_kpt(:,:), gen_wkpt(:)
     ! Band structure path
     integer :: nkptv,ngroups_kptv
     integer, dimension(:), pointer :: nkptsv_group
     real(gp), pointer :: kptv(:,:)
     character(len=100) :: band_structure_filename

     ! Orbitals and input occupation.
     integer :: gen_norb, gen_norbu, gen_norbd
     real(gp), dimension(:), pointer :: gen_occup

     ! Geometry variables from *.geopt
     character(len=10) :: geopt_approach !<id of geopt driver
     integer :: ncount_cluster_x !< Maximum number of geopt steps 
     integer :: wfn_history !< number of previous steps saved for wfn reformatting
     integer :: history !< History of DIIS method
     real(gp) :: betax,forcemax,randdis
     integer :: optcell, ionmov, nnos
     real(gp) :: dtion, mditemp, mdftemp, noseinert, friction, mdwall
     real(gp) :: bmass, vmass, strprecon, strfact
     real(gp) :: strtarget(6)
     real(gp), pointer :: qmass(:)
     real(gp) :: dtinit,dtmax !for FIRE
     ! tddft variables from *.tddft
     character(len=10) :: tddft_approach
     !variables for SIC
     type(SIC_data) :: SIC !<parameters for the SIC methods

     ! Performance variables from input.perf
     logical :: debug      !< Debug option (used by memocc)
     integer :: ncache_fft !< Cache size for FFT
     real(gp) :: projrad   !< Coarse radius of the projectors in units of the maxrad
     real(gp) :: symTol    !< Tolerance for symmetry detection.
     integer :: linear
     logical :: signaling  !< Expose results on DBus or Inet.
     integer :: signalTimeout !< Timeout for inet connection.
     character(len = 64) :: domain !< Domain to get the IP from hostname.
     character(len=500) :: writing_directory !< absolute path of the local directory to write the data on
     double precision :: gmainloop !< Internal C pointer on the signaling structure.
     integer :: inguess_geopt !< 0= Wavelet input guess, 1 = real space input guess 

     !orthogonalisation data
     type(orthon_data) :: orthpar
  
     !linear scaling data
     type(linearInputParameters) :: lin

     !fragment data
     type(fragmentInputParameters) :: frag

     !acceleration parameters
     type(material_acceleration) :: matacc

     !> parallelisation scheme of the exact exchange operator
     !!   BC (Blocking Collective)
     !!   OP2P (Overlap Point-to-Point)
     character(len=4) :: exctxpar

     !> paradigm for unblocking global communications via OMP_NESTING
     character(len=3) :: unblock_comms

     !> communication scheme for the density
     !!  DBL traditional scheme with double precision
     !!  MIX mixed single-double precision scheme (requires rho_descriptors)
     character(len=3) :: rho_commun
     !> number of taskgroups for the poisson solver
     !! works only if the number of MPI processes is a multiple of it
     integer :: PSolver_groupsize
     
     !> Global MPI group size (will be written in the mpi_environment)
     ! integer :: mpi_groupsize 

     !> linear scaling: store indices of the sparse matrices or recalculate them 
     logical :: store_index

     !> linear scaling: perform a check of sumrho (no check, light check or full check)
     integer :: check_sumrho

     !>linear scaling: activate the experimental mode
     logical :: experimental_mode

     !> linear scaling: write KS orbitals for cubic restart
     logical :: write_orbitals

     !> linear scaling: explicitely specify localization centers
     logical :: explicit_locregcenters

     !> linear scaling: calculate Kohn-Sham residue
     logical :: calculate_KS_residue
     
     !> linear scaling: calculate intermediate forces
     logical :: intermediate_forces

     !> linear scaling: exit kappa for extended input guess (experimental mode)
     real(kind=8) :: kappa_conv

     !> linear scaling: number of FOE cycles before the eigenvalue bounds are shrinked
     integer :: evbounds_nsatur

     !> linear scaling: maximal number of unsuccessful eigenvalue bounds shrinkings
     integer :: evboundsshrink_nsatur

     !> linear scaling: how to update the density kernel during the support function optimization (0: purification, 1: FOE)
     integer :: method_updatekernel

     !> linear scaling: quick return in purification
     logical :: purification_quickreturn
  end type input_variables

  !> Contains all energy terms
  type, public :: energy_terms
     real(gp) :: eh      =0.0_gp !< Hartree energy
     real(gp) :: exc     =0.0_gp !< Exchange-correlation
     real(gp) :: evxc    =0.0_gp
     real(gp) :: eion    =0.0_gp !< Ion-Ion interaction
     real(gp) :: edisp   =0.0_gp !< Dispersion force
     real(gp) :: ekin    =0.0_gp !< Kinetic term
     real(gp) :: epot    =0.0_gp
     real(gp) :: eproj   =0.0_gp
     real(gp) :: eexctX  =0.0_gp
     real(gp) :: ebs     =0.0_gp
     real(gp) :: eKS     =0.0_gp
     real(gp) :: trH     =0.0_gp
     real(gp) :: evsum   =0.0_gp
     real(gp) :: evsic   =0.0_gp 
     real(gp) :: excrhoc =0.0_gp 
     real(gp) :: eTS     =0.0_gp
     real(gp) :: ePV     =0.0_gp !< pressure term
     real(gp) :: energy  =0.0_gp !< the functional which is minimized
     real(gp) :: e_prev  =0.0_gp !< the previous value, to show the delta
     real(gp) :: trH_prev=0.0_gp !< the previous value, to show the delta
     !real(gp), dimension(:,:), pointer :: fion,f

     integer(kind = 8) :: c_obj = 0  !< Storage of the C wrapper object.
  end type energy_terms

  !> Memory estimation requirements
  type, public :: memory_estimation
     double precision :: submat
     integer :: ncomponents, norb, norbp
     double precision :: oneorb, allpsi_mpi, psistorage
     double precision :: projarr
     double precision :: grid
     double precision :: workarr

     double precision :: kernel, density, psolver, ham

     double precision :: peak
  end type memory_estimation

  !> Used to split between points to be treated in simple or in double precision
  type, public :: rho_descriptors
     character(len=1) :: geocode !< @copydoc poisson_solver::doc::geocode
     integer :: icomm !< method for communicating the density
     integer :: nrhotot !< dimension of the partial density array before communication
     integer :: n_csegs,n_fsegs,dp_size,sp_size
     integer, dimension(:,:), pointer :: spkey,dpkey
     integer, dimension(:), pointer :: cseg_b,fseg_b
  end type rho_descriptors

!> Contains arguments needed for rho_local for WVL+PAW

  type, public :: rholoc_objects
    integer ,pointer,dimension(:)    :: msz ! mesh size for local rho
    real(gp),pointer,dimension(:,:,:) :: d! local rho and derivatives
    real(gp),pointer,dimension(:,:)  :: rad!radial mesh for local rho
    real(gp),pointer,dimension(:) :: radius !after this radius, rholoc is zero
  end type rholoc_objects

  !> Structure to store the density / potential distribution among processors.
  type, public :: denspot_distribution
     integer :: n3d,n3p,n3pi,i3xcsh,i3s,nrhodim,i3rho_add
     integer :: ndimpot,ndimgrid,ndimrhopot 
     integer, dimension(3) :: ndims !< box containing the grid dimensions in ISF basis
     real(gp), dimension(3) :: hgrids !< grid spacings of the box (half of wavelet ones)
     integer, dimension(:,:), pointer :: nscatterarr, ngatherarr
     type(mpi_environment) :: mpi_env
  end type denspot_distribution

!>   Structures of basis of gaussian functions of the form exp(-a*r2)cos/sin(b*r2)
  type, public :: gaussian_basis_c
     integer :: nat,ncoeff,nshltot,nexpo
     integer, dimension(:), pointer :: nshell,ndoc,nam
     complex(gp), dimension(:), pointer :: expof,psiat
     real(gp), dimension(:,:), pointer :: rxyz
  end type gaussian_basis_c

  !> All the parameters which are important for describing the orbitals
  !! Add also the objects related to k-points sampling, after symmetries applications
  type, public :: orbitals_data 
     integer :: norb          !< Total number of orbitals per k point
     integer :: norbp         !< Total number of orbitals for the given processors
     integer :: norbu,norbd,nspin,nspinor,isorb
     integer :: nkpts,nkptsp,iskpts
     real(gp) :: efermi,HLgap,eTS
     integer, dimension(:), pointer :: iokpt,ikptproc,isorb_par,ispot
     integer, dimension(:), pointer :: inwhichlocreg,onwhichatom
     integer, dimension(:,:), pointer :: norb_par
     real(wp), dimension(:), pointer :: eval
     real(gp), dimension(:), pointer :: occup,spinsgn,kwgts
     real(gp), dimension(:,:), pointer :: kpts
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
  end type orbitals_data


  !> Contains the pointers to be handled to control GPU information
  !! Given that they are pointers on GPU address, they are C pointers
  !! which take 8 bytes
  !! So they are declared as kind=8 variables either if the GPU works in simple precision
  !! Also other information concerning the GPU runs can be stored in this structure
  type, public :: GPU_pointers
     logical :: useDynamic,full_locham
     integer :: id_proc,ndevices
     real(kind=8) :: keys,work1,work2,work3,rhopot,r,d
     real(kind=8) :: rhopot_down, rhopot_up
     real(kind=8) :: work1_i,work2_i,work3_i,d_i
     real(kind=8) :: pinned_in,pinned_out
     real(kind=8), dimension(:), pointer :: psi
     real(kind=8) :: psi_c,psi_f
     real(kind=8) :: psi_c_i,psi_f_i
     real(kind=8) :: psi_c_r,psi_f_r,psi_c_b,psi_f_b,psi_c_d,psi_f_d
     real(kind=8) :: psi_c_r_i,psi_f_r_i,psi_c_b_i,psi_f_b_i,psi_c_d_i,psi_f_d_i
     real(kind=8) :: keyg_c,keyg_f,keyv_c,keyv_f
     real(kind=8) :: keyg_c_host,keyg_f_host,keyv_c_host,keyv_f_host
     real(kind=8) :: context,queue
     !host pointers to be freed
     real(kind=8) :: rhopot_down_host, rhopot_up_host
     real(kind=8), dimension(:,:,:), pointer :: ekinpot_host
     real(kind=8), dimension(:,:), pointer :: psicf_host
     real(kind=8), dimension(:,:), pointer :: hpsicf_host
     real(kind=8), dimension(:), pointer :: bprecond_host

     real(gp), dimension(:,:), pointer :: ekin, epot !< values of the kinetic and potential energies to be passed to local_hamiltonian
     real(wp), dimension(:), pointer :: hpsi_ASYNC !<pointer to the wavefunction allocated in the case of asyncronous local_hamiltonian
  end type GPU_pointers

  !> Contains all the descriptors necessary for splitting the calculation in different locregs 
  type,public:: local_zone_descriptors
     logical :: linear                         !< if true, use linear part of the code
     integer :: nlr                            !< Number of localization regions 
     integer :: lintyp                         !< if 0 cubic, 1 locreg and 2 TMB
     integer :: ndimpotisf                      !< total dimension of potential in isf (including exctX)
     real(gp), dimension(3) :: hgrids          !<grid spacings of wavelet grid
     type(locreg_descriptors) :: Glr           !< Global region descriptors
     type(locreg_descriptors),dimension(:),pointer :: Llr                !< Local region descriptors (dimension = nlr)
  end type local_zone_descriptors

  !> Contains the work arrays needed for expressing wavefunction in real space
  !! with all the BC
  type, public :: workarr_sumrho
     integer :: nw1,nw2,nxc,nxf
     real(wp), dimension(:), pointer :: x_c,x_f,w1,w2
  end type workarr_sumrho


  !> Contains the work arrays needed for hamiltonian application with all the BC
  type, public :: workarr_locham
     integer :: nw1,nw2,nxc,nyc,nxf1,nxf2,nxf3,nxf,nyf
     real(wp), dimension(:), pointer :: w1,w2
     !for the periodic BC case, these arrays substitute 
     !psifscf,psifscfk,psig,ww respectively
     real(wp), dimension(:,:), pointer :: x_c,y_c,x_f1,x_f2,x_f3,x_f,y_f
  end type workarr_locham


  !> Contains the work arrays needed for th preconditioner with all the BC
  !! Take different pointers depending on the boundary conditions
  type, public :: workarr_precond
     integer, dimension(:), pointer :: modul1,modul2,modul3
     real(wp), dimension(:), pointer :: psifscf,ww,x_f1,x_f2,x_f3,kern_k1,kern_k2,kern_k3
     real(wp), dimension(:,:), pointer :: af,bf,cf,ef
     real(wp), dimension(:,:,:), pointer :: xpsig_c,ypsig_c,x_c
     real(wp), dimension(:,:,:,:), pointer :: xpsig_f,ypsig_f,x_f,y_f
     real(wp), dimension(:,:,:,:,:), pointer :: z1,z3 ! work array for FFT
  end type workarr_precond

  !> Contains all parameters needed for point to point communication
  type,public:: p2pComms
    integer,dimension(:),pointer :: noverlaps
    real(kind=8),dimension(:),pointer :: recvBuf
    integer,dimension(:,:,:),pointer :: comarr
    integer :: nrecvBuf, window
    integer,dimension(:,:),pointer :: ise ! starting / ending index of recvBuf in x,y,z dimension after communication (glocal coordinates)
    integer,dimension(:,:),pointer :: mpi_datatypes
    logical :: communication_complete
  end type p2pComms

  type,public :: foe_data
    integer,dimension(:),pointer :: kernel_nseg
    integer,dimension(:,:,:),pointer :: kernel_segkeyg
    real(kind=8) :: ef !< Fermi energy for FOE
    real(kind=8) :: evlow, evhigh !< eigenvalue bounds for FOE 
    real(kind=8) :: bisection_shift !< bisection shift to find Fermi energy (FOE)
    real(kind=8) :: fscale !< length scale for complementary error function (FOE)
    real(kind=8) :: ef_interpol_det !<FOE: max determinant of cubic interpolation matrix
    real(kind=8) :: ef_interpol_chargediff !<FOE: max charge difference for interpolation
    real(kind=8) :: charge !total charge of the system
    integer :: evbounds_isatur, evboundsshrink_isatur, evbounds_nsatur, evboundsshrink_nsatur !< variables to check whether the eigenvalue bounds might be too big
  end type foe_data

!!$  type, public :: sparse_matrix_metadata
!!$     integer :: nvctr, nseg, full_dim1, full_dim2
!!$     integer,dimension(:),pointer :: noverlaps
!!$     integer,dimension(:,:),pointer :: overlaps
!!$     integer,dimension(:),pointer :: keyv, nsegline, istsegline
!!$     integer,dimension(:,:),pointer :: keyg
!!$     integer,dimension(:,:),pointer :: matrixindex_in_compressed, orb_from_index
!!$  end type sparse_matrix_metadata

  !!type,public :: sparse_matrix
  !!    integer :: nvctr, nseg, nvctrp, isvctr, parallel_compression, nfvctr, nfvctrp, isfvctr
  !!    integer,dimension(:),pointer :: keyv, nsegline, istsegline, isvctr_par, nvctr_par, isfvctr_par, nfvctr_par
  !!    integer,dimension(:,:),pointer :: keyg
  !!    !type(sparse_matrix_metadata), pointer :: pattern
  !!    real(kind=8),dimension(:),pointer :: matrix_compr,matrix_comprp
  !!    real(kind=8),dimension(:,:),pointer :: matrix,matrixp
  !!    !integer,dimension(:,:),pointer :: matrixindex_in_compressed, orb_from_index
  !!    integer,dimension(:,:),pointer :: matrixindex_in_compressed_arr, orb_from_index
  !!    integer,dimension(:,:),pointer :: matrixindex_in_compressed_fortransposed
  !!    logical :: store_index, can_use_dense
  !!    !!contains
  !!    !!  procedure,pass :: matrixindex_in_compressed
  !!end type sparse_matrix

  type,public :: linear_matrices !may not keep
      type(sparse_matrix) :: ham, ovrlp, denskern_large, inv_ovrlp_large
  end type linear_matrices

  type,public:: workarrays_quartic_convolutions
    real(wp),dimension(:,:,:),pointer :: xx_c, xy_c, xz_c
    real(wp),dimension(:,:,:),pointer :: xx_f1
    real(wp),dimension(:,:,:),pointer :: xy_f2
    real(wp),dimension(:,:,:),pointer :: xz_f4
    real(wp),dimension(:,:,:,:),pointer :: xx_f, xy_f, xz_f
    real(wp),dimension(:,:,:),pointer :: y_c
    real(wp),dimension(:,:,:,:),pointer :: y_f
    ! The following arrays are work arrays within the subroutine
    real(wp),dimension(:,:),pointer :: aeff0array, beff0array, ceff0array, eeff0array
    real(wp),dimension(:,:),pointer :: aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array
    real(wp),dimension(:,:),pointer :: aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray
    real(wp),dimension(:,:,:),pointer :: xya_c, xyc_c
    real(wp),dimension(:,:,:),pointer :: xza_c, xzc_c
    real(wp),dimension(:,:,:),pointer :: yza_c, yzb_c, yzc_c, yze_c
    real(wp),dimension(:,:,:,:),pointer :: xya_f, xyb_f, xyc_f, xye_f
    real(wp),dimension(:,:,:,:),pointer :: xza_f, xzb_f, xzc_f, xze_f
    real(wp),dimension(:,:,:,:),pointer :: yza_f, yzb_f, yzc_f, yze_f
  end type workarrays_quartic_convolutions

  type,public:: localizedDIISParameters
    integer :: is, isx, mis, DIISHistMax, DIISHistMin
    integer :: icountSDSatur, icountDIISFailureCons, icountSwitch, icountDIISFailureTot, itBest
    real(kind=8),dimension(:),pointer :: phiHist, hphiHist, energy_hist
    real(kind=8) :: alpha_coeff !step size for optimization of coefficients
    real(kind=8),dimension(:,:,:),pointer :: mat
    real(kind=8) :: trmin, trold, alphaSD, alphaDIIS
    logical :: switchSD, immediateSwitchToSD, resetDIIS
  end type localizedDIISParameters


  type,public:: mixrhopotDIISParameters
    integer :: is, isx, mis
    real(kind=8),dimension(:),pointer :: rhopotHist, rhopotresHist
    real(kind=8),dimension(:,:),pointer :: mat
  end type mixrhopotDIISParameters

  !> Contains the arguments needed for the diis procedure
  type, public :: diis_objects
     logical :: switchSD
     integer :: idiistol,mids,ids,idsx
     real(gp) :: energy_min,energy_old,energy,alpha,alpha_max
     real(tp), dimension(:), pointer :: psidst
     real(tp), dimension(:), pointer :: hpsidst
     real(tp), dimension(:,:,:,:,:,:), pointer :: ads
  end type diis_objects

  !> Contains the information needed for the preconditioner
  type, public :: precond_data
    integer :: confPotOrder                           !< The order of the algebraic expression for Confinement potential
    integer :: ncong                                  !< Number of CG iterations for the preconditioning equation
    logical, dimension(:), pointer :: withConfPot     !< Use confinement potentials
    real(kind=8), dimension(:), pointer :: potentialPrefac !< Prefactor for the potentiar : Prefac * f(r) 
  end type precond_data

  !> Information for the confining potential to be used in TMB scheme
  !! The potential is supposed to be defined as prefac*(r-rC)**potorder
  type, public :: confpot_data
     integer :: potorder                !< order of the confining potential
     integer, dimension(3) :: ioffset   !< offset for the coordinates of potential lr in global region
     real(gp) :: prefac                 !< prefactor
     real(gp), dimension(3) :: hh       !< grid spacings in ISF grid
     real(gp), dimension(3) :: rxyzConf !< confining potential center in global coordinates
  end type confpot_data

  !> Defines the important information needed to reformat a old wavefunctions
  type, public :: old_wavefunction
     type(local_zone_descriptors) :: Lzd !< local zone descriptors of the corresponding run
     real(wp), dimension(:), pointer :: psi !<wavelets coefficients in compressed form
     real(gp), dimension(:,:), pointer :: rxyz !<atomic positions of the step
  end type old_wavefunction

  !> Densities and potentials, and related metadata, needed for their creation/application
  !! Not all these quantities are available, some of them may point to the same memory space
  type, public :: DFT_local_fields
     real(dp), dimension(:), pointer :: rhov !< generic workspace. What is there is indicated by rhov_is
     
     type(ab7_mixing_object), pointer :: mix          !< History of rhov, allocated only when using diagonalisation
     !local fields which are associated to their name
     !normally given in parallel distribution
     real(dp), dimension(:,:), pointer :: rho_psi !< density as given by square of el. WFN
     real(dp), dimension(:,:,:,:), pointer :: rho_C   !< core density
     real(wp), dimension(:,:,:,:), pointer :: V_ext   !< local part of pseudopotientials
     real(wp), dimension(:,:,:,:), pointer :: V_XC    !< eXchange and Correlation potential (local)
     real(wp), dimension(:,:,:,:), pointer :: Vloc_KS !< complete local potential of KS Hamiltonian (might point on rho_psi)
     real(wp), dimension(:,:,:,:), pointer :: f_XC !< dV_XC[rho]/d_rho
     !temporary arrays
     real(wp), dimension(:), pointer :: rho_work,pot_work !<full grid arrays
     !metadata
     integer :: rhov_is
     real(gp) :: psoffset !< offset of the Poisson Solver in the case of Periodic BC
     type(rho_descriptors) :: rhod !< descriptors of the density for parallel communication
     type(denspot_distribution) :: dpbox !< distribution of density and potential box
     character(len=3) :: PSquiet
     !real(gp), dimension(3) :: hgrids !<grid spacings of denspot grid (half of the wvl grid)
     type(coulomb_operator) :: pkernel !< kernel of the Poisson Solver used for V_H[rho]
     type(coulomb_operator) :: pkernelseq !<for monoproc PS (useful for exactX, SIC,...)

     integer(kind = 8) :: c_obj = 0                !< Storage of the C wrapper object.
  end type DFT_local_fields

  !> Flags for rhov status
  integer, parameter, public :: EMPTY              = -1980
  integer, parameter, public :: ELECTRONIC_DENSITY = -1979
  integer, parameter, public :: CHARGE_DENSITY     = -1978
  integer, parameter, public :: KS_POTENTIAL       = -1977
  integer, parameter, public :: HARTREE_POTENTIAL  = -1976

  !> Flags for the restart (linear scaling only)
  integer, parameter, public :: LINEAR_LOWACCURACY  = 101 !low accuracy after restart
  integer, parameter, public :: LINEAR_HIGHACCURACY = 102 !high accuracy after restart

  !check if all comms are necessary here
  type, public :: hamiltonian_descriptors
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     type(local_zone_descriptors) :: Lzd !< data on the localisation regions, if associated
     type(comms_linear) :: collcom ! describes collective communication
     type(p2pComms) :: comgp           !<describing p2p communications for distributing the potential
     real(wp), dimension(:), pointer :: psi,psit_c,psit_f !< these should eventually be eliminated
     logical :: can_use_transposed
  end type hamiltonian_descriptors

  !> The wavefunction which have to be considered at the DFT level
  type, public :: DFT_wavefunction
     !coefficients
     real(wp), dimension(:), pointer :: psi,hpsi,psit,psit_c,psit_f !< orbitals, or support functions, in wavelet basis
     real(wp), dimension(:), pointer :: spsi !< Metric operator applied to psi (To be used for PAW)
     real(wp), dimension(:,:), pointer :: gaucoeffs !orbitals in gbd basis
     !basis sets
     type(gaussian_basis) :: gbd !<gaussian basis description, if associated
     type(local_zone_descriptors) :: Lzd !< data on the localisation regions, if associated
     !restart objects (consider to move them in rst structure)
     type(old_wavefunction), dimension(:), pointer :: oldpsis !< previously calculated wfns
     integer :: istep_history !< present step of wfn history
     !data properties
     logical :: can_use_transposed !< true if the transposed quantities are allocated and can be used
     type(orbitals_data) :: orbs !<wavefunction specification in terms of orbitals
     type(comms_cubic) :: comms !< communication objects for the cubic approach
     type(diis_objects) :: diis
     type(confpot_data), dimension(:), pointer :: confdatarr !<data for the confinement potential
     type(SIC_data) :: SIC !<control the activation of SIC scheme in the wavefunction
     type(orthon_data) :: orthpar !< control the application of the orthogonality scheme for cubic DFT wavefunction
     character(len=4) :: exctxpar !< Method for exact exchange parallelisation for the wavefunctions, in case
     type(p2pComms) :: comgp !<describing p2p communications for distributing the potential
     type(comms_linear) :: collcom ! describes collective communication
     type(comms_linear) :: collcom_sr ! describes collective communication for the calculation of the charge density
     integer(kind = 8) :: c_obj !< Storage of the C wrapper object. it has to be initialized to zero
     type(foe_data) :: foe_obj        !<describes the structure of the matrices for the linear method foe
     type(linear_matrices) :: linmat
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     type(hamiltonian_descriptors) :: ham_descr
     real(kind=8),dimension(:,:),pointer :: coeff !<expansion coefficients
  end type DFT_wavefunction

  !> Flags for optimization loop id
  integer, parameter, public :: OPTLOOP_HAMILTONIAN   = 0
  integer, parameter, public :: OPTLOOP_SUBSPACE      = 1
  integer, parameter, public :: OPTLOOP_WAVEFUNCTIONS = 2
  integer, parameter, public :: OPTLOOP_N_LOOPS       = 3

  !> Used to control the optimization of wavefunctions
  type, public :: DFT_optimization_loop
     integer :: iscf !< Kind of optimization scheme.

     integer :: itrpmax !< specify the maximum number of mixing cycle on potential or density
     integer :: nrepmax !< specify the maximum number of restart after re-diagonalization
     integer :: itermax !< specify the maximum number of minimization iterations, self-consistent or not

     integer :: itrp    !< actual number of mixing cycle.
     integer :: itrep   !< actual number of re-diagonalisation runs.
     integer :: iter    !< actual number of minimization iterations.

     integer :: infocode !< return value after optimization loop.

     real(gp) :: gnrm   !< actual value of cv criterion of the minimization loop.
     real(gp) :: rpnrm  !< actual value of cv criterion of the mixing loop.

     real(gp) :: gnrm_cv       !< convergence criterion of the minimization loop.
     real(gp) :: rpnrm_cv      !< convergence criterion of the mixing loop.
     real(gp) :: gnrm_startmix !< gnrm value to start mixing after.

     integer(kind = 8) :: c_obj = 0 !< Storage of the C wrapper object.
  end type DFT_optimization_loop

  !>  Used to restart a new DFT calculation or to save information 
  !!  for post-treatment
  type, public :: restart_objects
     integer :: version !< 0=cubic, 100=linear
     integer :: n1,n2,n3,nat
     real(gp) :: hx_old,hy_old,hz_old
     real(gp), dimension(:,:), pointer :: rxyz_old,rxyz_new
     type(DFT_wavefunction) :: KSwfn !< Kohn-Sham wavefunctions
     type(DFT_wavefunction) :: tmb !<support functions for linear scaling
     type(GPU_pointers) :: GPU 
  end type restart_objects

  !> Public container to be used with call_bigdft().
  type, public :: run_objects
     type(dictionary), pointer :: user_inputs

     type(input_variables), pointer    :: inputs
     type(atoms_data), pointer         :: atoms
     type(restart_objects), pointer    :: rst
     real(gp), dimension(:,:), pointer :: radii_cf
  end type run_objects

  !> Used to store results of a DFT calculation.
  type, public :: DFT_global_output
     real(gp) :: energy, fnoise, pressure
     type(energy_terms) :: energs
     integer :: fdim                           !< Dimension of allocated forces (second dimension)
     real(gp), dimension(:,:), pointer :: fxyz
     real(gp), dimension(6) :: strten
  end type DFT_global_output

!> type paw_ij_objects

 type paw_ij_objects

!Integer scalars

  integer :: cplex
   ! cplex=1 if all on-site PAW quantities are real, 2 if they are complex
   ! cplex=2 is useful for RF calculations

  integer :: cplex_dij
   ! cplex=1 if dij are real, 2 if they are complex

  !!!!$integer :: has_dijexxcore !> does this makes sense?
   ! 1 if dijexxcore is allocated
   ! 2 if dijexxcore is already computed

  integer :: has_dij
   ! 1 if dij is allocated
   ! 2 if dij is already computed

  integer :: has_dijfr
   ! 1 if dijfr is allocated
   ! 2 if dijfr is already computed

  integer :: has_dijhartree
   ! 1 if dijhartree is allocated
   ! 2 if dijhartree is already computed

  integer :: has_dijhat
   ! 1 if dijhat is allocated
   ! 2 if dijhat is already computed

  integer :: has_dijso
   ! 1 if dijso is associated and used, 0 otherwise
   ! 2 if dijso is already computed

  integer :: has_dijU
   ! 1 if dijU is associated and used, 0 otherwise
   ! 2 if dijU is already computed

  integer :: has_dijxc
   ! 1 if dijxc is associated and used, 0 otherwise
   ! 2 if dijxc is already computed

  integer :: has_dijxc_val
   ! 1 if dijxc_val is associated and used, 0 otherwise
   ! 2 if dijxc_val is already computed

  integer :: has_exexch_pot
   ! 1 if PAW+(local exact exchange) potential is allocated

  integer :: has_pawu_occ
   ! 1 if PAW+U occupations are allocated

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: ndij
   ! Number of components of dij
   ! Usually ndij=nspden, except for nspinor==2 (where ndij=nspinor**2)

  integer :: nspden
   ! Number of spin-density components (may be different from dtset%nspden if spin-orbit)

  integer :: nsppol
   ! Number of independant spin-components

!Real (real(dp)) arrays

  real(dp), pointer :: dij(:,:)
   ! dij(cplex_dij*lmn2_size,ndij)
   ! Dij term (non-local operator)
   ! May be complex if cplex_dij=2
   !  dij(:,:,1) contains Dij^up-up
   !  dij(:,:,2) contains Dij^dn-dn
   !  dij(:,:,3) contains Dij^up-dn (only if nspinor=2)
   !  dij(:,:,4) contains Dij^dn-up (only if nspinor=2)

  !real(dp),pointer :: dijexxcore(:,:)
  ! dijexxcore(cplex_dij*lmn2_size,ndij)
  ! Onsite matrix elements of the Fock operator generated by core electrons

!  real(dp), pointer :: dijfr(:,:)
!   ! dijhat(cplex_dij*lmn2_size,ndij)
!   ! For response function calculation only
!   ! RF Frozen part of Dij (depends on q vector but not on 1st-order wave function)
!   ! Same storage as Dij (see above)
!
!  real(dp), pointer :: dijhartree(:)
!   ! dijhartree(cplex*lmn2_size)
!   ! Dij_hartree term
!   ! Contains all contributions to Dij from hartree
!   ! Warning: Dimensioned by cplex, not cplex_dij
!   ! Same storage as Dij (see above)
!
!  real(dp), pointer :: dijhat(:,:)
!   ! dijhat(cplex_dij*lmn2_size,ndij)
!   ! Dij_hat term (non-local operator) i.e \sum_LM \int_FFT Q_{ij}^{LM} vtrial
!   ! Same storage as Dij (see above)
!
!  real(dp), pointer :: dijU(:,:)
!   ! dijU(cplex_dij*lmn2_size,ndij)
!   ! Onsite matrix elements of the U part of the PAW Hamiltonian.
!   ! Same storage as Dij (see above)
!
!  real(dp), pointer :: dijso(:,:)
!   ! dijso(cplex_dij*lmn2_size,ndij)
!   ! Onsite matrix elements of L.S i.e <phi_i|L.S|phi_j>
!   ! Same storage as Dij (see above)
!
!  real(dp), pointer :: dijxc(:,:)
!   ! dijxc(cplex_dij*lmn2_size,ndij)
!   ! Onsite matrix elements of vxc i.e
!   ! <phi_i|vxc[n1+nc]|phi_j> - <tphi_i|vxc(tn1+nhat+tnc]|tphi_j>
!   ! Same storage as Dij (see above)
!
!  real(dp), pointer :: dijxc_val(:,:)
!   ! dijxc_val(cplex_dij*lmn2_size,ndij)
!   ! Onsite matrix elements of valence-only vxc i.e
!   ! <phi_i|vxc[n1]|phi_j> - <tphi_i|vxc(tn1+nhat]|tphi_j>
!   ! Same storage as Dij (see above)
!
!  real(dp), pointer :: noccmmp(:,:,:,:)
!   ! noccmmp(cplex_dij,2*lpawu+1,2*lpawu+1,nocc_nspden)
!   ! cplex_dij=1 if collinear
!   ! cplex_dij=2 if spin orbit is used
!   ! cplex_dij=2 is used if non-collinear (for coherence, it is not necessary in this case, however)
!   ! gives occupation matrix for lda+u (computed in setnoccmmp)
!   ! Stored as: noccmmp(:,:,1)=   n^{up,up}_{m,mp}
!   !            noccmmp(:,:,2)=   n^{dn,dn}_{m,mp}
!   !            noccmmp(:,:,3)=   n^{up,dn}_{m,mp}
!   !            noccmmp(:,:,4)=   n^{dn,up}_{m,mp}
!   ! noccmmp(m,mp,:) is computed from rhoij(klmn) with  m=klmntomn(2)>mp=klmntomn(1)
!
!  real(dp), pointer :: nocctot(:)
!   ! nocctot(nspden)
!   ! gives trace of occupation matrix for lda+u (computed in pawdenpot)
!   ! for each value of ispden (1 or 2)
!
!  real(dp), pointer :: vpawx(:,:,:)
!   ! vpawx(2*lexexch+1,2*lexexch+1,nspden)
!   ! exact exchange potential

 end type paw_ij_objects

!>This is cprj_type in ABINIT,
!!this will be obsolete with the PAW Library
 type cprj_objects

!Integer scalars

  integer :: ncpgr
   ! Number of gradients of cp=<p_lmn|Cnk>

  integer :: nlmn
   ! Number of (l,m,n) non-local projectors

!Real (real(dp)) arrays

  real(wp), pointer :: cp (:,:)
   ! cp(2,nlmn)
   ! <p_lmn|Cnk> projected scalars for a given atom and wave function

  real(wp), pointer :: dcp (:,:,:)
   ! dcp(2,ncpgr,nlmn)
   ! derivatives of <p_lmn|Cnk> projected scalars for a given atom and wave function

 end type cprj_objects

!> Contains the arguments needed for the PAW implementation:
  type, public :: paw_objects
    integer :: lmnmax
    integer :: ntypes
    integer :: natom
    integer :: usepaw
    integer,dimension(:,:,:),pointer :: indlmn
    type(paw_ij_objects),dimension(:),allocatable :: paw_ij
    type(cprj_objects),dimension(:,:),allocatable :: cprj
    real(wp),dimension(:),pointer :: spsi
    real(wp),dimension(:,:),pointer :: sij
    real(gp),dimension(:),pointer :: rpaw
  end type paw_objects

  interface input_set
     module procedure input_set_char, input_set_int, input_set_dbl, input_set_bool, &
          & input_set_int_array, input_set_dbl_array, input_set_bool_array, &
          & input_set_dict
  end interface input_set

  !>timing categories
  character(len=*), parameter, private :: tgrp_pot='Potential'
  integer, save, public :: TCAT_EXCHANGECORR=TIMING_UNINITIALIZED
  integer, parameter, private :: ncls_max=6,ncat_bigdft=138   ! define timimg categories and classes
  character(len=14), dimension(ncls_max), parameter, private :: clss = (/ &
       'Communications'    ,  &
       'Convolutions  '    ,  &
       'Linear Algebra'    ,  &
       'Other         '    ,  &
!       'Potential     '    ,  &
       'Initialization'    ,  &
       'Finalization  '    /)
  character(len=14), dimension(3,ncat_bigdft), parameter, private :: cats = reshape((/ &
                                !       Name           Class       Operation Kind
       'ReformatWaves ','Initialization' ,'Small Convol  ' ,  &  !< Reformatting of input waves
       'CrtDescriptors','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of descriptor arrays
       'CrtLocPot     ','Initialization' ,'Miscellaneous ' ,  &  !< Calculation of local potential
       'CrtProjectors ','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of projectors
       'CrtPcProjects ','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of preconditioning projectors
       'CrtPawProjects','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of abscalc-pawprojectors
       'ApplyLocPotKin','Convolutions  ' ,'OpenCL ported ' ,  &  !< Application of PSP, kinetic energy
       'ApplyProj     ','Other         ' ,'RMA pattern   ' ,  &  !< Application of nonlocal PSP
       'Precondition  ','Convolutions  ' ,'OpenCL ported ' ,  &  !< Precondtioning
       'Rho_comput    ','Convolutions  ' ,'OpenCL ported ' ,  &  !< Calculation of charge density (sumrho) computation
       'Rho_commun    ','Communications' ,'AllReduce grid' ,  &  !< Calculation of charge density (sumrho) communication
       'Pot_commun    ','Communications' ,'AllGathrv grid' ,  &  !< Communication of potential
       'Pot_comm start','Communications' ,'MPI_types/_get' ,  &  !< Communication of potential
       'Un-TransSwitch','Other         ' ,'RMA pattern   ' ,  &  !< Transposition of wavefunction, computation
       'Un-TransComm  ','Communications' ,'ALLtoALLV     ' ,  &  !< Transposition of wavefunction, communication
       'GramS_comput  ','Linear Algebra' ,'DPOTRF        ' ,  &  !< Gram Schmidt computation        
       'GramS_commun  ','Communications' ,'ALLReduce orbs' ,  &  !< Gram Schmidt communication
       'LagrM_comput  ','Linear Algebra' ,'DGEMM         ' ,  &  !< Lagrange Multipliers computation
       'LagrM_commun  ','Communications' ,'ALLReduce orbs' ,  &  !< Lagrange Multipliers communication
       'Diis          ','Other         ' ,'Other         ' ,  &  
       !       'PSolv_comput  ','Potential     ' ,'3D FFT        ' ,  &  
       !       'PSolv_commun  ','Communications' ,'ALLtoALL      ' ,  &  
       !       'PSolvKernel   ','Initialization' ,'Miscellaneous ' ,  &  
!       'Exchangecorr  ','Potential     ' ,'Miscellaneous ' ,  &  
       'Forces        ','Finalization  ' ,'Miscellaneous ' ,  &  
       'Tail          ','Finalization  ' ,'Miscellaneous ' ,  &
       'Loewdin_comput','Linear Algebra' ,'              ' ,  &
       'Loewdin_commun','Communications' ,'ALLReduce orbs' ,  &
       'Chol_commun   ','Communications' ,'              ' ,  &
       'Chol_comput   ','Linear Algebra' ,'ALLReduce orbs' ,  &
       'GS/Chol_comput','Linear Algebra' ,'              ' ,  &
       'GS/Chol_commun','Communications' ,'ALLReduce orbs' ,  &
       'Input_comput  ','Initialization' ,'Miscellaneous ' ,  &
       'Input_commun  ','Communications' ,'ALLtoALL+Reduc' ,  &
       'Davidson      ','Finalization  ' ,'Complete SCF  ' ,  &
       'check_IG      ','Initialization' ,'Linear Scaling' ,  &
       'constrc_locreg','Initialization' ,'Miscellaneous ' ,  &
       'wavefunction  ','Initialization' ,'Miscellaneous ' ,  &
       'create_nlpspd ','Initialization' ,'RMA pattern   ' ,  &
       'p2pOrtho_post ','Communications' ,'irecv / irsend' ,  &
       'p2pOrtho_wait ','Communications' ,'mpi_waitany   ' ,  &
       'lovrlp_comm   ','Communications' ,'mpi_allgatherv' ,  &
       'lovrlp_comp   ','Linear Algebra' ,'many ddots    ' ,  &
       'lovrlp_compr  ','Other         ' ,'cut out zeros ' ,  &
       'lovrlp_uncompr','Other         ' ,'insert zeros  ' ,  &
       'extract_orbs  ','Other         ' ,'copy to sendb ' ,  &
       'lovrlp^-1/2   ','Linear Algebra' ,'exact or appr ' ,  &
       'lovrlp^-1/2old','Linear Algebra' ,'exact or appr ' ,  &
       'lovrlp^-1/2com','Linear Algebra' ,'exact or appr ' ,  &
       'lovrlp^-1/2par','Linear Algebra' ,'exact or appr ' ,  &
       'build_lincomb ','Linear Algebra' ,'many daxpy    ' ,  &
       'convolQuartic ','Convolutions  ' ,'No OpenCL     ' ,  &
       'p2pSumrho_wait','Communications' ,'mpi_test/wait ' ,  &
       'sumrho_TMB    ','Other         ' ,'port to GPU?  ' ,  &
       'TMB_kernel    ','Linear Algebra' ,'dgemm         ' ,  &
       'diagonal_seq  ','Linear Algebra' ,'dsygv         ' ,  &
       'diagonal_par  ','Linear Algebra' ,'pdsygvx       ' ,  &
       'lovrlp^-1     ','Linear Algebra' ,'exact or appr ' ,  &
       'lagmat_orthoco','Linear Algebra' ,'dgemm seq/par ' ,  &
       'optimize_DIIS ','Other         ' ,'Other         ' ,  &
       'optimize_SD   ','Other         ' ,'Other         ' ,  &
       'mix_linear    ','Other         ' ,'Other         ' ,  &
       'mix_DIIS      ','Other         ' ,'Other         ' ,  &
       'ig_matric_comm','Communications' ,'mpi p2p       ' ,  &
       'wf_signals    ','Communications' ,'Socket transf.' ,  &
       'energs_signals','Communications' ,'Socket transf.' ,  &
       'rhov_signals  ','Communications' ,'Socket transf.' ,  &
       'init_locregs  ','Initialization' ,'Miscellaneous ' ,  &
       'init_commSumro','Initialization' ,'Miscellaneous ' ,  &
       'init_commPot  ','Initialization' ,'Miscellaneous ' ,  &
       'init_commOrtho','Initialization' ,'Miscellaneous ' ,  &
       'init_inguess  ','Initialization' ,'Miscellaneous ' ,  &
       'init_matrCompr','Initialization' ,'Miscellaneous ' ,  &
       'init_collcomm ','Initialization' ,'Miscellaneous ' ,  &
       'init_collco_sr','Initialization' ,'Miscellaneous ' ,  &
       'init_orbs_lin ','Initialization' ,'Miscellaneous ' ,  &
       'init_repart   ','Initialization' ,'Miscellaneous ' ,  &
       'initMatmulComp','Initialization' ,'Miscellaneous ' ,  &
       'Pot_after_comm','Other         ' ,'global_to_loca' ,  & 
       !       'Init to Zero  ','Other         ' ,'Memset        ' ,  &
       'calc_kernel   ','Other         ' ,'Miscellaneous ' ,  &
       'commun_kernel ','Communications' ,'mpi_allgatherv' ,  &
       'getlocbasinit ','Other         ' ,'Miscellaneous ' ,  &
       'updatelocreg1 ','Other         ' ,'Miscellaneous ' ,  &
       'linscalinit   ','Other         ' ,'Miscellaneous ' ,  &
       'commbasis4dens','Communications' ,'Miscellaneous ' ,  &
       'eglincomms    ','Communications' ,'Miscellaneous ' ,  &
       'allocommsumrho','Communications' ,'Miscellaneous ' ,  &
       'ovrlptransComp','Other         ' ,'Miscellaneous ' ,  &
       'ovrlptransComm','Communications' ,'mpi_allreduce ' ,  &
       'lincombtrans  ','Other         ' ,'Miscellaneous ' ,  &
       'glsynchham1   ','Other         ' ,'Miscellaneous ' ,  &
       'glsynchham2   ','Other         ' ,'Miscellaneous ' ,  &
       'gauss_proj    ','Other         ' ,'Miscellaneous ' ,  &
       'sumrho_allred ','Communications' ,'mpiallred     ' ,  &
       'deallocprec   ','Other         ' ,'Miscellaneous ' ,  &
       'large2small   ','Other         ' ,'Miscellaneous ' ,  &
       'small2large   ','Other         ' ,'Miscellaneous ' ,  &
       'renormCoefCom1','Linear Algebra' ,'Miscellaneous ' ,  &
       'renormCoefCom2','Linear Algebra' ,'Miscellaneous ' ,  &
       'renormCoefComm','Communications' ,'Miscellaneous ' ,  &
       'waitAllgatKern','Other         ' ,'Miscellaneous ' ,  &
       'UnBlockPot    ','Other         ' ,'Overlap comms ' ,  &
       'UnBlockDen    ','Other         ' ,'Overlap comms ' ,  &
       'global_local  ','Initialization' ,'Unknown       ' ,  &
       'wfd_creation  ','Other         ' ,'Miscellaneous ' ,  & 
       'comm_llr      ','Communications' ,'Miscellaneous ' ,  &
       !       'AllocationProf','Other         ' ,'Allocate arrs ' ,  &
       'dirmin_lagmat1','Linear Algebra' ,'grad calc     ' ,  &
       'dirmin_lagmat2','Linear Algebra' ,'allgatherv    ' ,  &
       'dirmin_dgesv  ','Linear Algebra' ,'dgesv/pdgesv  ' ,  &
       'dirmin_sddiis ','Linear Algebra' ,'Miscellaneous ' ,  &
       'dirmin_allgat ','Linear Algebra' ,'allgatherv    ' ,  &
       'dirmin_sdfit  ','Linear Algebra' ,'allgatherv etc' ,  &
       'chebyshev_comp','Linear Algebra' ,'matmul/matadd ' ,  &
       'chebyshev_comm','Communications' ,'allreduce     ' ,  &
       'chebyshev_coef','Other         ' ,'Miscellaneous ' ,  &
       'FOE_auxiliary ','Other         ' ,'Miscellaneous ' ,  &
       'FOE_init      ','Other         ' ,'Miscellaneous ' ,  &
       'compress_uncom','Other         ' ,'Miscellaneous ' ,  &
       'norm_trans    ','Other         ' ,'Miscellaneous ' ,  &
       'misc          ','Other         ' ,'Miscellaneous ' ,  &
       'sparse_copy   ','Other         ' ,'Miscellaneous ' ,  &
       'constraineddft','Other         ' ,'Miscellaneous ' ,  &
       'transfer_int  ','Other         ' ,'Miscellaneous ' ,  &
       'Reformatting  ','Initialization' ,'Interpolation ' ,  &
       'restart_wvl   ','Initialization' ,'inguess    rst' ,  &
       'restart_rsp   ','Initialization' ,'inguess    rst' ,  &
       'check_sumrho  ','Initialization' ,'unitary check ' ,  &
       'check_pot     ','Initialization' ,'unitary check ' ,  &
       'ApplyLocPot   ','Convolutions  ' ,'OpenCL ported ' ,  &
       'ApplyLocKin   ','Convolutions  ' ,'OpenCL ported ' ,  &
       'kernel_init   ','Other         ' ,'Fragment calc ' ,  &
       'calc_energy   ','Linear Algebra' ,'allred etc    ' ,  &
       'new_pulay_corr','Other         ' ,'Pulay forces  ' ,  &
       'dev_from_unity','Other         ' ,'Miscellaneous ' ,  &
       'ks_residue    ','Linear Algebra' ,'Miscellaneous ' ,  &
       'weightanalysis','Linear Algebra' ,'Fragment calc ' ,  &
       'tmbrestart    ','Initialization' ,'Miscellaneous ' ,  &
       'readtmbfiles  ','Initialization' ,'Miscellaneous ' ,  &
       'readisffiles  ','Initialization' ,'Miscellaneous ' ,  &
       'purify_kernel ','Linear Algebra' ,'dgemm         ' ,  &
       'potential_dims','Other         ' ,'auxiliary     ' ,  &
       'calc_bounds   ','Other         ' ,'Miscellaneous ' /),(/3,ncat_bigdft/))
  integer, dimension(ncat_bigdft), private, save :: cat_ids !< id of the categories to be converted

contains

  function old_wavefunction_null() result(wfn)
    implicit none
    type(old_wavefunction) :: wfn
    wfn%Lzd=default_lzd()
    nullify(wfn%psi)
    nullify(wfn%rxyz)
  end function old_wavefunction_null

  function dpbox_null() result(dd)
    implicit none
    type(denspot_distribution) :: dd
    dd%n3d=0
    dd%n3p=0
    dd%n3pi=0
    dd%i3xcsh=0
    dd%i3s=0
    dd%nrhodim=0
    dd%i3rho_add=0
    dd%ndimpot=0
    dd%ndimgrid=0
    dd%ndimrhopot=0
    dd%ndims=(/0,0,0/)
    dd%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    nullify(dd%nscatterarr)
    nullify(dd%ngatherarr)
    dd%mpi_env=mpi_environment_null()
  end function dpbox_null

  function material_acceleration_null() result(ma)
    type(material_acceleration) :: ma
    ma%iacceleration=0
    ma%Psolver_igpu=0
    ma%OCL_platform=repeat(' ',len(ma%OCL_platform))
    ma%OCL_platform=repeat(' ',len(ma%OCL_devices))
  end function material_acceleration_null

  function default_lzd() result(lzd)
    type(local_zone_descriptors) :: lzd
    lzd%linear=.false.
    lzd%nlr=0
    lzd%lintyp=0
    lzd%ndimpotisf=0
    lzd%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    lzd%Glr=locreg_null()
    nullify(lzd%Llr)
  end function default_lzd
 
  function bigdft_run_id_toa()
    use yaml_output
    implicit none
    character(len=20) :: bigdft_run_id_toa

    bigdft_run_id_toa=repeat(' ',len(bigdft_run_id_toa))

    if (bigdft_mpi%ngroup>1) then
       bigdft_run_id_toa=adjustl(trim(yaml_toa(bigdft_mpi%igroup,fmt='(i15)')))
    end if

  end function bigdft_run_id_toa

  !> Fills the old_wavefunction structure with corresponding data
  !! Deallocate previous workspaces if already existing
  subroutine old_wavefunction_set(wfn,nat,norbp,Lzd,rxyz,psi)
    
    implicit none
    integer, intent(in) :: nat,norbp
    type(local_zone_descriptors), intent(in) :: Lzd
    real(gp), dimension(3,nat), intent(in) :: rxyz
    real(wp), dimension((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*norbp), intent(in) :: psi
    type(old_wavefunction), intent(inout) :: wfn
    !local variables
    character(len=*), parameter :: subname='old_wavefunction_set'
    integer :: i_stat

    !first, free the workspace if not already done
    call old_wavefunction_free(wfn,subname)
    !then allocate the workspaces and fill them
    allocate(wfn%psi((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*norbp+ndebug),stat=i_stat)
    call memocc(i_stat,wfn%psi,'psi',subname)
    
    if (norbp>0) call vcopy((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*norbp,&
         psi(1),1,wfn%psi(1),1)

    allocate(wfn%rxyz(3,nat+ndebug),stat=i_stat)
    call memocc(i_stat,wfn%rxyz,'rxyz',subname)
    if (nat>0) call vcopy(3*nat,rxyz(1,1),1,wfn%rxyz(1,1),1)
    call copy_local_zone_descriptors(Lzd,wfn%Lzd,subname)

  end subroutine old_wavefunction_set

  subroutine old_wavefunction_free(wfn,subname)
    implicit none
    character(len=*), intent(in) :: subname
    type(old_wavefunction), intent(inout) :: wfn
    !local variables
    integer :: i_all,i_stat

    if (associated(wfn%psi)) then
       i_all=-product(shape(wfn%psi))*kind(wfn%psi)
       deallocate(wfn%psi,stat=i_stat)
       call memocc(i_stat,i_all,'psi',subname)
    end if
    if (associated(wfn%rxyz)) then
       i_all=-product(shape(wfn%rxyz))*kind(wfn%rxyz)
       deallocate(wfn%rxyz,stat=i_stat)
       call memocc(i_stat,i_all,'rxyz',subname)
    end if
    !lzd should be deallocated also (to be checked again)
    call deallocate_local_zone_descriptors(wfn%Lzd, subname)

  end subroutine old_wavefunction_free
   

!> De-Allocate comms_cubic
  subroutine deallocate_comms(comms,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(comms_cubic), intent(inout) :: comms
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(comms%nvctr_par))*kind(comms%nvctr_par)
    deallocate(comms%nvctr_par,stat=i_stat)
    call memocc(i_stat,i_all,'nvctr_par',subname)
    i_all=-product(shape(comms%ncntd))*kind(comms%ncntd)
    deallocate(comms%ncntd,stat=i_stat)
    call memocc(i_stat,i_all,'ncntd',subname)
    i_all=-product(shape(comms%ncntt))*kind(comms%ncntt)
    deallocate(comms%ncntt,stat=i_stat)
    call memocc(i_stat,i_all,'ncntt',subname)
    i_all=-product(shape(comms%ndspld))*kind(comms%ndspld)
    deallocate(comms%ndspld,stat=i_stat)
    call memocc(i_stat,i_all,'ndspld',subname)
    i_all=-product(shape(comms%ndsplt))*kind(comms%ndsplt)
    deallocate(comms%ndsplt,stat=i_stat)
    call memocc(i_stat,i_all,'ndsplt',subname)
  END SUBROUTINE deallocate_comms


  subroutine deallocate_abscalc_input(in, subname)
    use module_base
    implicit none
    type(input_variables) :: in
    character(len=*), intent(in) :: subname

    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(in%Gabs_coeffs))*kind(in%Gabs_coeffs)
    deallocate(in%Gabs_coeffs, stat=i_stat)
    call memocc(i_stat,i_all,'in%Gabs_coeffs',subname)

  END SUBROUTINE deallocate_abscalc_input


!> De-Allocate orbitals data structure, except eval pointer
!! which is not allocated in the orbitals_descriptor routine
subroutine deallocate_orbs(orbs,subname)
  use module_base
  implicit none
    character(len=*), intent(in) :: subname
    type(orbitals_data), intent(inout) :: orbs
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(orbs%norb_par))*kind(orbs%norb_par)
    deallocate(orbs%norb_par,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%norb_par',subname)
    i_all=-product(shape(orbs%occup))*kind(orbs%occup)
    deallocate(orbs%occup,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%occup',subname)
    i_all=-product(shape(orbs%spinsgn))*kind(orbs%spinsgn)
    deallocate(orbs%spinsgn,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%spinsgn',subname)
    i_all=-product(shape(orbs%kpts))*kind(orbs%kpts)
    deallocate(orbs%kpts,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%kpts',subname)
    i_all=-product(shape(orbs%kwgts))*kind(orbs%kwgts)
    deallocate(orbs%kwgts,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%kwgts',subname)
    i_all=-product(shape(orbs%iokpt))*kind(orbs%iokpt)
    deallocate(orbs%iokpt,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%iokpt',subname)
    i_all=-product(shape(orbs%ikptproc))*kind(orbs%ikptproc)
    deallocate(orbs%ikptproc,stat=i_stat)
    call memocc(i_stat,i_all,'ikptproc',subname)
    i_all=-product(shape(orbs%inwhichlocreg))*kind(orbs%inwhichlocreg)
    deallocate(orbs%inwhichlocreg,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%inwhichlocreg',subname)
    i_all=-product(shape(orbs%onwhichatom))*kind(orbs%onwhichatom)
    deallocate(orbs%onwhichatom,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%onwhichatom',subname)
    i_all=-product(shape(orbs%isorb_par))*kind(orbs%isorb_par)
    deallocate(orbs%isorb_par,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%isorb_par',subname)
    if (associated(orbs%ispot)) then
       i_all=-product(shape(orbs%ispot))*kind(orbs%ispot)
       deallocate(orbs%ispot,stat=i_stat)
       call memocc(i_stat,i_all,'orbs%ispot',subname)
    end if

  END SUBROUTINE deallocate_orbs

  !> All in one routine to initialise and set-up restart objects.
  subroutine init_restart_objects(iproc,inputs,atoms,rst,subname)
    use module_base
    implicit none
    !Arguments
    character(len=*), intent(in) :: subname
    integer, intent(in) :: iproc
    type(input_variables), intent(in) :: inputs
    type(atoms_data), intent(in) :: atoms
    type(restart_objects), intent(out) :: rst

    call restart_objects_new(rst)
    call restart_objects_set_mode(rst, inputs%inputpsiid)
    call restart_objects_set_nat(rst, atoms%astruct%nat, subname)
    call restart_objects_set_mat_acc(rst, iproc, inputs%matacc)
  END SUBROUTINE init_restart_objects

  !> Allocate and nullify restart objects
  subroutine restart_objects_new(rst)
    use module_base
    implicit none
    !Arguments
    type(restart_objects), intent(out) :: rst

    ! Decide whether we use the cubic or the linear version
    rst%version = UNINITIALIZED(CUBIC_VERSION)

    !allocate pointers
    rst%nat = 0
    nullify(rst%rxyz_new)
    nullify(rst%rxyz_old)

    !nullify unallocated pointers
    rst%KSwfn%c_obj = 0
    nullify(rst%KSwfn%psi)
    nullify(rst%KSwfn%orbs%eval)

    nullify(rst%KSwfn%gaucoeffs)
    nullify(rst%KSwfn%oldpsis)

    rst%KSwfn%Lzd%Glr = locreg_null()
    nullify(rst%KSwfn%Lzd%Glr%wfd%keyglob)
    nullify(rst%KSwfn%Lzd%Glr%wfd%keygloc)
    nullify(rst%KSwfn%Lzd%Glr%wfd%keyvloc)
    nullify(rst%KSwfn%Lzd%Glr%wfd%keyvglob)
                
    nullify(rst%KSwfn%gbd%nshell)
    nullify(rst%KSwfn%gbd%ndoc)
    nullify(rst%KSwfn%gbd%nam)
    nullify(rst%KSwfn%gbd%xp)
    nullify(rst%KSwfn%gbd%psiat)
    nullify(rst%KSwfn%gbd%rxyz)

    !Nullify LZD for cubic version (new input guess)
    call nullify_local_zone_descriptors(rst%tmb%lzd)
  END SUBROUTINE restart_objects_new

  subroutine restart_objects_set_mode(rst, inputpsiid)
    implicit none
    type(restart_objects), intent(inout) :: rst
    integer, intent(in) :: inputpsiid

    select case (inputpsiid)
    case (INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_CP2K, INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, &
         INPUT_PSI_DISK_WVL, INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS)
       rst%version = CUBIC_VERSION
    case (INPUT_PSI_LINEAR_AO, INPUT_PSI_MEMORY_LINEAR, INPUT_PSI_DISK_LINEAR)
       rst%version = LINEAR_VERSION
    end select
  END SUBROUTINE restart_objects_set_mode
  subroutine restart_objects_set_nat(rst, nat, subname)
    use module_base
    implicit none
    !Arguments
    character(len=*), intent(in) :: subname
    integer, intent(in) :: nat
    type(restart_objects), intent(inout) :: rst
    !local variables
    integer :: i_all,i_stat
    if (associated(rst%rxyz_old)) then
       i_all=-product(shape(rst%rxyz_old))*kind(rst%rxyz_old)
       deallocate(rst%rxyz_old,stat=i_stat)
       call memocc(i_stat,i_all,'rxyz_old',subname)
    end if
    if (associated(rst%rxyz_new)) then
       i_all=-product(shape(rst%rxyz_new))*kind(rst%rxyz_new)
       deallocate(rst%rxyz_new,stat=i_stat)
       call memocc(i_stat,i_all,'rxyz_new',subname)
    end if

    rst%nat = nat
    allocate(rst%rxyz_new(3,nat+ndebug),stat=i_stat)
    call memocc(i_stat,rst%rxyz_new,'rxyz_new',subname)
    allocate(rst%rxyz_old(3,nat+ndebug),stat=i_stat)
    call memocc(i_stat,rst%rxyz_old,'rxyz_old',subname)
  END SUBROUTINE restart_objects_set_nat

  subroutine restart_objects_set_mat_acc(rst, iproc, matacc)
    implicit none
    !Arguments
    type(restart_objects), intent(inout) :: rst
    integer, intent(in) :: iproc
    type(material_acceleration), intent(in) :: matacc
    !initialise the acceleration strategy if required
    call init_material_acceleration(iproc,matacc,rst%GPU)
  END SUBROUTINE restart_objects_set_mat_acc

!>  De-Allocate restart_objects
  subroutine free_restart_objects(rst,subname)
    use module_base
    use locregs
    implicit none
    character(len=*), intent(in) :: subname
    type(restart_objects) :: rst
    !local variables
    integer :: i_all,i_stat,istep

    if (rst%version == LINEAR_VERSION) then
       call destroy_DFT_wavefunction(rst%tmb)
    end if
    !always deallocate lzd for new input guess
    !call deallocate_lzd(rst%tmb%lzd, subname)
    ! Modified by SM
    call deallocate_local_zone_descriptors(rst%tmb%lzd, subname)

    call deallocate_locreg_descriptors(rst%KSwfn%Lzd%Glr)

    if (associated(rst%KSwfn%psi)) then
       i_all=-product(shape(rst%KSwfn%psi))*kind(rst%KSwfn%psi)
       deallocate(rst%KSwfn%psi,stat=i_stat)
       call memocc(i_stat,i_all,'psi',subname)
    end if

    if (associated(rst%KSwfn%orbs%eval)) then
       i_all=-product(shape(rst%KSwfn%orbs%eval))*kind(rst%KSwfn%orbs%eval)
       deallocate(rst%KSwfn%orbs%eval,stat=i_stat)
       call memocc(i_stat,i_all,'eval',subname)
    end if

    if (associated(rst%KSwfn%oldpsis)) then
       do istep=0,product(shape(rst%KSwfn%oldpsis))-1
          call old_wavefunction_free(rst%KSwfn%oldpsis(istep),subname)
       end do
       deallocate(rst%KSwfn%oldpsis)
    end if


    if (associated(rst%rxyz_old)) then
       i_all=-product(shape(rst%rxyz_old))*kind(rst%rxyz_old)
       deallocate(rst%rxyz_old,stat=i_stat)
       call memocc(i_stat,i_all,'rxyz_old',subname)
    end if
    if (associated(rst%rxyz_new)) then
       i_all=-product(shape(rst%rxyz_new))*kind(rst%rxyz_new)
       deallocate(rst%rxyz_new,stat=i_stat)
       call memocc(i_stat,i_all,'rxyz_new',subname)
    end if

    !The gaussian basis descriptors are always allocated together
    !with the gaussian coefficients
    if (associated(rst%KSwfn%gbd%rxyz)) then
       nullify(rst%KSwfn%gbd%rxyz)
       call deallocate_gwf(rst%KSwfn%gbd,subname)
    end if

    if (associated(rst%KSwfn%gaucoeffs)) then
       i_all=-product(shape(rst%KSwfn%gaucoeffs))*kind(rst%KSwfn%gaucoeffs)
       deallocate(rst%KSwfn%gaucoeffs,stat=i_stat)
       call memocc(i_stat,i_all,'gaucoeffs',subname)
    end if

    !finalise the material accelearion usage
    call release_material_acceleration(rst%GPU)

  END SUBROUTINE free_restart_objects


!> Allocate wavefunctions_descriptors
  subroutine deallocate_rho_descriptors(rhodsc,subname)
    use module_base
    implicit none
    type(rho_descriptors) :: rhodsc
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    if (associated(rhodsc%spkey))then
       i_all=-product(shape(rhodsc%spkey))*kind(rhodsc%spkey)
       deallocate(rhodsc%spkey,stat=i_stat)
       call memocc(i_stat,i_all,'spkey',subname)
    end if
    if (associated(rhodsc%dpkey))then
       i_all=-product(shape(rhodsc%dpkey))*kind(rhodsc%dpkey)
       deallocate(rhodsc%dpkey,stat=i_stat)
       call memocc(i_stat,i_all,'dpkey',subname)
    end if
    if (associated(rhodsc%cseg_b))then
       i_all=-product(shape(rhodsc%cseg_b))*kind(rhodsc%cseg_b)
       deallocate(rhodsc%cseg_b,stat=i_stat)
       call memocc(i_stat,i_all,'csegb',subname)
    end if
    if (associated(rhodsc%fseg_b))then
       i_all=-product(shape(rhodsc%fseg_b))*kind(rhodsc%fseg_b)
       deallocate(rhodsc%fseg_b,stat=i_stat)
       call memocc(i_stat,i_all,'fsegb',subname)
    end if

  end subroutine deallocate_rho_descriptors


!> De-Allocate gaussian_basis type
  subroutine deallocate_gwf(G,subname)
    use module_base
    implicit none
    type(gaussian_basis) :: G
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    !normally positions should be deallocated outside
    
    i_all=-product(shape(G%ndoc))*kind(G%ndoc)
    deallocate(G%ndoc,stat=i_stat)
    call memocc(i_stat,i_all,'ndoc',subname)
    i_all=-product(shape(G%nam))*kind(G%nam)
    deallocate(G%nam,stat=i_stat)
    call memocc(i_stat,i_all,'nam',subname)
    i_all=-product(shape(G%nshell))*kind(G%nshell)
    deallocate(G%nshell,stat=i_stat)
    call memocc(i_stat,i_all,'nshell',subname)
    i_all=-product(shape(G%psiat))*kind(G%psiat)
    deallocate(G%psiat,stat=i_stat)
    call memocc(i_stat,i_all,'psiat',subname)
    i_all=-product(shape(G%xp))*kind(G%xp)
    deallocate(G%xp,stat=i_stat)
    call memocc(i_stat,i_all,'xp',subname)

  END SUBROUTINE deallocate_gwf


!> De-Allocate gaussian_basis type
  subroutine deallocate_gwf_c(G,subname)
    use module_base
    implicit none
    type(gaussian_basis_c) :: G
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    !normally positions should be deallocated outside
    
    i_all=-product(shape(G%ndoc))*kind(G%ndoc)
    deallocate(G%ndoc,stat=i_stat)
    call memocc(i_stat,i_all,'G%ndoc',subname)
    i_all=-product(shape(G%nam))*kind(G%nam)
    deallocate(G%nam,stat=i_stat)
    call memocc(i_stat,i_all,'nam',subname)
    i_all=-product(shape(G%nshell))*kind(G%nshell)
    deallocate(G%nshell,stat=i_stat)
    call memocc(i_stat,i_all,'G%nshell',subname)
    i_all=-product(shape(G%psiat))*kind(G%psiat)
    deallocate(G%psiat,stat=i_stat)
    call memocc(i_stat,i_all,'G%psiat',subname)

    i_all=-product(shape(G%expof))*kind(G%expof)
    deallocate(G%expof,stat=i_stat)
    call memocc(i_stat,i_all,'G%expof',subname)

    i_all=-product(shape(G%rxyz))*kind(G%rxyz)
    deallocate(G%rxyz,stat=i_stat)
    call memocc(i_stat,i_all,'G%rxyz',subname)

  END SUBROUTINE 

  !> Deallocate lr (obsolete)
  !! @todo Remove this function.
  subroutine deallocate_lr(lr,subname)
    use module_base
    character(len=*), intent(in) :: subname
    type(locreg_descriptors) :: lr
!    integer :: i_all,i_stat

    write(0,*) "deallocate_lr : TODO, remove me"
    
    call deallocate_wfd(lr%wfd)

    call deallocate_bounds(lr%geocode,lr%hybrid_on,lr%bounds,subname)

!    if (associated(lr%projflg)) then
!       i_all=-product(shape(lr%projflg)*kind(lr%projflg))
!       deallocate(lr%projflg,stat=i_stat)
!       call memocc(i_stat,i_all,'lr%projflg',subname)
!    end if
  END SUBROUTINE deallocate_lr

  subroutine deallocate_Lzd(Lzd,subname)
    use module_base
    character(len=*), intent(in) :: subname
    type(local_zone_descriptors) :: Lzd
    integer :: ilr

!   nullify the bounds of Glr
    if ((Lzd%Glr%geocode == 'P' .and. Lzd%Glr%hybrid_on) .or. Lzd%Glr%geocode == 'F') then
       nullify(Lzd%Glr%bounds%kb%ibyz_f)
       nullify(Lzd%Glr%bounds%kb%ibxz_f)
       nullify(Lzd%Glr%bounds%kb%ibxy_f)
       nullify(Lzd%Glr%bounds%sb%ibxy_ff)
       nullify(Lzd%Glr%bounds%sb%ibzzx_f)
       nullify(Lzd%Glr%bounds%sb%ibyyzz_f)
       nullify(Lzd%Glr%bounds%gb%ibyz_ff)
       nullify(Lzd%Glr%bounds%gb%ibzxx_f)
       nullify(Lzd%Glr%bounds%gb%ibxxyy_f)
    end if
    !the arrays which are needed only for free BC
    if (Lzd%Glr%geocode == 'F') then
       nullify(Lzd%Glr%bounds%kb%ibyz_c)
       nullify(Lzd%Glr%bounds%kb%ibxz_c)
       nullify(Lzd%Glr%bounds%kb%ibxy_c)
       nullify(Lzd%Glr%bounds%sb%ibzzx_c)
       nullify(Lzd%Glr%bounds%sb%ibyyzz_c)
       nullify(Lzd%Glr%bounds%gb%ibzxx_c)
       nullify(Lzd%Glr%bounds%gb%ibxxyy_c)
       nullify(Lzd%Glr%bounds%ibyyzz_r)
    end if

! nullify the wfd of Glr
   nullify(Lzd%Glr%wfd%keyglob)
   nullify(Lzd%Glr%wfd%keygloc)
!   nullify(Lzd%Glr%wfd%keyv)
   nullify(Lzd%Glr%wfd%keyvloc)
   nullify(Lzd%Glr%wfd%keyvglob)

! nullify the Gnlpspd
!   call deallocate_proj_descr(Lzd%Gnlpspd,subname)
!!$   nullify(Lzd%Gnlpspd%nvctr_p)
!!$   nullify(Lzd%Gnlpspd%nseg_p)
!!$   nullify(Lzd%Gnlpspd%keyv_p)
!!$   nullify(Lzd%Gnlpspd%keyg_p)
!!$   nullify(Lzd%Gnlpspd%nboxp_c)
!!$   nullify(Lzd%Gnlpspd%nboxp_f)
 
!Now destroy the Llr
    do ilr = 1, Lzd%nlr 
       call deallocate_lr(Lzd%Llr(ilr),subname)
!       call deallocate_Lnlpspd(Lzd%Lnlpspd(ilr),subname)
    end do
     nullify(Lzd%Llr)
!     nullify(Lzd%Lnlpspd)

  END SUBROUTINE deallocate_Lzd


  function input_psi_names(id)
    integer, intent(in) :: id
    character(len = 14) :: input_psi_names

    select case(id)
    case(INPUT_PSI_EMPTY)
       write(input_psi_names, "(A)") "empty"
    case(INPUT_PSI_RANDOM)
       write(input_psi_names, "(A)") "random"
    case(INPUT_PSI_CP2K)
       write(input_psi_names, "(A)") "CP2K"
    case(INPUT_PSI_LCAO)
       write(input_psi_names, "(A)") "LCAO"
    case(INPUT_PSI_MEMORY_WVL)
       write(input_psi_names, "(A)") "wvl. in mem."
    case(INPUT_PSI_DISK_WVL)
       write(input_psi_names, "(A)") "wvl. on disk"
    case(INPUT_PSI_LCAO_GAUSS)
       write(input_psi_names, "(A)") "LCAO + gauss."
    case(INPUT_PSI_MEMORY_GAUSS)
       write(input_psi_names, "(A)") "gauss. in mem."
    case(INPUT_PSI_DISK_GAUSS)
       write(input_psi_names, "(A)") "gauss. on disk"
    case(INPUT_PSI_LINEAR_AO)
       write(input_psi_names, "(A)") "Linear AO"
    case(INPUT_PSI_MEMORY_LINEAR)
       write(input_psi_names, "(A)") "Linear restart"
    case(INPUT_PSI_DISK_LINEAR)
       write(input_psi_names, "(A)") "Linear on disk"
    case default
       write(input_psi_names, "(A)") "Error"
    end select
  end function input_psi_names

  subroutine input_psi_help()
    integer :: i

    write(*, "(1x,A)") "Available values of inputPsiId are:"
    do i = 1, size(input_psi_values)
       write(*, "(1x,A,I5,A,A)") " | ", input_psi_values(i), &
            & " - ", input_psi_names(input_psi_values(i))
    end do
  end subroutine input_psi_help

  function input_psi_validate(id)
    integer, intent(in) :: id
    logical :: input_psi_validate

    integer :: i

    input_psi_validate = .false.
    do i = 1, size(input_psi_values)
       if (id == input_psi_values(i)) then
          input_psi_validate = .true.
          return
       end if
    end do
  end function input_psi_validate

  subroutine output_wf_format_help()
    integer :: i

    write(*, "(1x,A)") "Available values of output_wf are:"
    do i = 0, size(wf_format_names) - 1
       write(*, "(1x,A,I5,A,A)") " | ", i, &
            & " - ", wf_format_names(i)
    end do
  end subroutine output_wf_format_help

  function output_wf_format_validate(id)
    integer, intent(in) :: id
    logical :: output_wf_format_validate

    output_wf_format_validate = (id >= 0 .and. id < size(wf_format_names))
  end function output_wf_format_validate

  subroutine output_denspot_help()
    integer :: i, j

    write(*, "(1x,A)") "Available values of output_denspot are:"
    do i = 0, size(output_denspot_format_names) - 1
       do j = 0, size(output_denspot_names) - 1
          if (j == 0 .and. i == 0) then
             write(*, "(1x,A,I5,A,A,A)") " | ", i * 10 + j, &
                  & " - ", trim(output_denspot_names(j)), "."
          else if (j /= 0) then
             write(*, "(1x,A,I5,A,A,A,A,A)") " | ", i * 10 + j, &
                  & " - ", trim(output_denspot_names(j)), &
                  & " in ", trim(output_denspot_format_names(i)), " format."
          end if
       end do
    end do
  end subroutine output_denspot_help

  function output_denspot_validate(id, fid)
    integer, intent(in) :: id, fid
    logical :: output_denspot_validate

    output_denspot_validate = (id >= 0 .and. id < size(output_denspot_names)) .and. &
         & (fid >= 0 .and. fid < size(output_denspot_format_names))
  end function output_denspot_validate

subroutine nullify_DFT_local_fields(denspot)
  implicit none
  type(DFT_local_fields),intent(out) :: denspot

  nullify(denspot%rhov)
  nullify(denspot%mix)
  nullify(denspot%rho_psi)
  nullify(denspot%rho_C)
  nullify(denspot%V_ext)
  nullify(denspot%V_XC)
  nullify(denspot%Vloc_KS)
  nullify(denspot%f_XC)
  nullify(denspot%rho_work)
  nullify(denspot%pot_work)
  call nullify_rho_descriptors(denspot%rhod)
  call nullify_denspot_distribution(denspot%dpbox)
  call nullify_coulomb_operator(denspot%pkernel)
  call nullify_coulomb_operator(denspot%pkernelseq)
  
end subroutine nullify_DFT_local_fields


subroutine deallocate_denspot_distribution(dpbox,subname)
  implicit none
  character(len=*), intent(in) :: subname
  type(denspot_distribution),intent(inout)::dpbox
  !local variables
  integer :: i_all,i_stat
  
  if(associated(dpbox%nscatterarr)) then
    i_all=-product(shape(dpbox%nscatterarr))*kind(dpbox%nscatterarr)
    deallocate(dpbox%nscatterarr,stat=i_stat)
    call memocc(i_stat,i_all,'nscatterarr',subname)
  end if
  if(associated(dpbox%ngatherarr)) then
    i_all=-product(shape(dpbox%ngatherarr))*kind(dpbox%ngatherarr)
    deallocate(dpbox%ngatherarr,stat=i_stat)
    call memocc(i_stat,i_all,'ngatherarr',subname)
  end if

end subroutine deallocate_denspot_distribution

subroutine nullify_coulomb_operator(coul_op)
  implicit none
  type(coulomb_operator),intent(out) :: coul_op
  nullify(coul_op%kernel)
end subroutine nullify_coulomb_operator

subroutine copy_coulomb_operator(coul1,coul2,subname)
  implicit none
  type(coulomb_operator),intent(in)::coul1
  type(coulomb_operator),intent(out)::coul2
  character(len=*), intent(in) :: subname
  !local variables
  integer :: i_all,i_stat

  if(associated(coul2%kernel)) then
    i_all=-product(shape(coul2%kernel))*kind(coul2%kernel)
    deallocate(coul2%kernel,stat=i_stat)
    call memocc(i_stat,i_all,'coul%kernel',subname)
  end if
  coul2%kernel   =>coul1%kernel
  coul2%itype_scf= coul1%itype_scf
  coul2%mu       = coul1%mu
  coul2%geocode  = coul1%geocode
  coul2%ndims    = coul1%ndims
  coul2%hgrids   = coul1%hgrids
  coul2%angrad   = coul1%angrad
  coul2%work1_GPU= coul1%work1_GPU
  coul2%work2_GPU= coul1%work2_GPU
  coul2%k_GPU    = coul1%k_GPU
  coul2%plan     = coul1%plan
  coul2%geo      = coul1%geo
  coul2%igpu     = coul1%igpu
  coul2%initCufftPlan=coul1%initCufftPlan
  coul2%keepGPUmemory=coul1%keepGPUmemory
! mpi_env:
  coul2%mpi_env%mpi_comm = coul1%mpi_env%mpi_comm
  coul2%mpi_env%iproc    = coul1%mpi_env%iproc
  coul2%mpi_env%nproc    = coul1%mpi_env%nproc
  coul2%mpi_env%igroup   = coul1%mpi_env%igroup
  coul2%mpi_env%ngroup   = coul1%mpi_env%ngroup

end subroutine copy_coulomb_operator

subroutine deallocate_coulomb_operator(coul_op,subname)
  implicit none
  type(coulomb_operator),intent(out)::coul_op
  character(len=*), intent(in) :: subname
  !local variables
  integer :: i_all,i_stat

  if(associated(coul_op%kernel)) then
    i_all=-product(shape(coul_op%kernel))*kind(coul_op%kernel)
    deallocate(coul_op%kernel,stat=i_stat)
    call memocc(i_stat,i_all,'coul_op%kernel',subname)
  end if
  call nullify_coulomb_operator(coul_op)
end subroutine deallocate_coulomb_operator


subroutine nullify_denspot_distribution(dpbox)
  implicit none
  type(denspot_distribution),intent(out) :: dpbox
  
  nullify(dpbox%nscatterarr)
  nullify(dpbox%ngatherarr)
end subroutine nullify_denspot_distribution

subroutine nullify_rho_descriptors(rhod)
  implicit none
  type(rho_descriptors),intent(out) :: rhod

  nullify(rhod%spkey)
  nullify(rhod%dpkey)
  nullify(rhod%cseg_b)
  nullify(rhod%fseg_b)
end subroutine nullify_rho_descriptors

subroutine nullify_GPU_pointers(gpup)
  implicit none
  type(GPU_pointers), intent(out) :: gpup

  nullify(gpup%psi)
  nullify(gpup%ekinpot_host)
  nullify(gpup%psicf_host)
  nullify(gpup%hpsicf_host)
  nullify(gpup%bprecond_host)
  nullify(gpup%ekin)
  nullify(gpup%epot)

end subroutine nullify_GPU_pointers

!subroutine nullify_wfn_metadata(wfnmd)
!  implicit none
!  type(wfn_metadata),intent(inout) :: wfnmd
!
!  nullify(wfnmd%coeff)
!  nullify(wfnmd%coeff_proj)
!  nullify(wfnmd%alpha_coeff)
!  nullify(wfnmd%grad_coeff_old)
!
!end subroutine nullify_wfn_metadata

subroutine nullify_diis_objects(diis)
  implicit none
  type(diis_objects),intent(inout) :: diis

  nullify(diis%psidst)
  nullify(diis%hpsidst)
  nullify(diis%ads)

end subroutine nullify_diis_objects

subroutine nullify_rholoc_objects(rholoc)
  implicit none
  type(rholoc_objects),intent(inout) :: rholoc
  
  nullify(rholoc%msz)
  nullify(rholoc%d)
  nullify(rholoc%rad)
  nullify(rholoc%radius) 
end subroutine nullify_rholoc_objects

subroutine nullify_paw_objects(paw,rholoc)
  implicit none
  type(paw_objects),intent(inout) :: paw
  type(rholoc_objects),optional :: rholoc
  
  nullify(paw%indlmn) 
  nullify(paw%spsi) 
  nullify(paw%sij) 
  nullify(paw%rpaw)

  if(present(rholoc)) then
   nullify(rholoc%msz)
   nullify(rholoc%d)
   nullify(rholoc%rad)
   nullify(rholoc%radius) 
  end if
end subroutine nullify_paw_objects

subroutine nullify_paw_ij_objects(paw_ij)
  implicit none
  !Arguments
  type(paw_ij_objects), intent(inout) :: paw_ij

  nullify(paw_ij%dij) 
end subroutine nullify_paw_ij_objects

subroutine nullify_cprj_objects(cprj)
  implicit none
  type(cprj_objects),intent(inout) :: cprj

  nullify(cprj%cp)
  nullify(cprj%dcp)
end subroutine nullify_cprj_objects

subroutine nullify_gaussian_basis(G)

  implicit none
  !Arguments
  type(gaussian_basis),intent(inout) :: G 

  G%ncplx=1
  nullify(G%nshell)
  nullify(G%ndoc)
  nullify(G%nam)
  nullify(G%psiat)
  nullify(G%xp)
  nullify(G%rxyz)

END SUBROUTINE nullify_gaussian_basis

subroutine nullify_global_output(outs)
  implicit none
  type(DFT_global_output), intent(out) :: outs

  outs%fdim      = 0
  nullify(outs%fxyz)
  outs%energy    = UNINITIALIZED(1.0_gp)
  outs%fnoise    = UNINITIALIZED(1.0_gp)
  outs%pressure  = UNINITIALIZED(1.0_gp)
  outs%strten(:) = UNINITIALIZED(1.0_gp)
END SUBROUTINE nullify_global_output

subroutine init_global_output(outs, nat)
  use module_base
  implicit none
  type(DFT_global_output), intent(out) :: outs
  integer, intent(in) :: nat

  call nullify_global_output(outs)
  outs%fdim = nat
  allocate(outs%fxyz(3, outs%fdim))
  outs%fxyz(:,:) = UNINITIALIZED(1.0_gp)
END SUBROUTINE init_global_output

subroutine deallocate_global_output(outs, fxyz)
  use module_base
  implicit none
  type(DFT_global_output), intent(inout) :: outs
  real(gp), intent(out), optional :: fxyz

  if (associated(outs%fxyz)) then
     if (present(fxyz)) then
        call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz, 1)
     end if
     deallocate(outs%fxyz)
  end if
END SUBROUTINE deallocate_global_output

!> cprj_clean will be obsolete with the PAW library
!! this is cprj_free in abinit.
 subroutine cprj_clean(cprj)

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(cprj_objects),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)
!write(std_out,*) "cprj_free ndim = ", n1dim, n2dim
 do jj=1,n2dim
   do ii=1,n1dim
     if (associated(cprj(ii,jj)%cp))  then
       deallocate(cprj(ii,jj)%cp)
     end if
     if (associated(cprj(ii,jj)%dcp))  then
       deallocate(cprj(ii,jj)%dcp)
     end if
   end do
 end do
end subroutine cprj_clean

!> this routine is cprj_alloc in abinit
!! with the PAW library this will be obsolet.
 subroutine cprj_paw_alloc(cprj,ncpgr,nlmn)

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncpgr
!arrays
 integer,intent(in) :: nlmn(:)
 type(cprj_objects),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim,nn

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2);nn=size(nlmn,dim=1)
 if (nn/=n1dim) then
   write(*,*)"Error in cprj_alloc: wrong sizes !",nn,n1dim
   stop
 end if
!write(std_out,*) "cprj_alloc ndim = ", n1dim, n2dim
 do jj=1,n2dim
   do ii=1,n1dim
     nullify (cprj(ii,jj)%cp)
     nullify (cprj(ii,jj)%dcp)

     nn=nlmn(ii)
     cprj(ii,jj)%nlmn=nn
     ALLOCATE(cprj(ii,jj)%cp(2,nn))
!    XG 080820 Was needed to get rid of problems with test paral#R with four procs
     cprj(ii,jj)%cp=0.0_dp
!    END XG 080820

     cprj(ii,jj)%ncpgr=ncpgr
     if (ncpgr>0) then
       ALLOCATE(cprj(ii,jj)%dcp(2,ncpgr,nn))
       cprj(ii,jj)%dcp=0.0_dp
     end if
   end do
 end do
end subroutine cprj_paw_alloc

subroutine cprj_to_array(cprj,array,norb,nspinor,shift,option)
  implicit none
  integer,intent(in) :: option,norb,nspinor,shift
  real(kind=8),intent(inout) :: array(:,:)
  type(cprj_objects),intent(inout) :: cprj(:)
  !
  integer :: ii,jj,ilmn,iorb
  !
  if(option==1) then
    do iorb=1,norb*nspinor
      ii=0
      do ilmn=1,cprj(iorb+shift)%nlmn
        do jj=1,2
          ii=ii+1
          array(ii,iorb)=cprj(iorb+shift)%cp(jj,ilmn)
        end do
      end do
    end do
  elseif(option==2) then
    do iorb=1,norb*nspinor
      ii=0
      do ilmn=1,cprj(iorb+shift)%nlmn
        do jj=1,2
          ii=ii+1
          cprj(iorb+shift)%cp(jj,ilmn)=array(ii,iorb)
        end do
      end do
    end do
  end if
end subroutine cprj_to_array

!> create a null Lzd. Note: this is the correct way of defining 
!! association through prure procedures.
!! A pure subroutine has to be defined to create a null structure.
!! this is important when using the nullification inside other
!! nullification routines since the usage of a pure function is forbidden
!! otherwise the routine cannot be pure
pure function local_zone_descriptors_null() result(lzd)
  implicit none
  type(local_zone_descriptors) :: lzd
  call nullify_local_zone_descriptors(lzd)
end function local_zone_descriptors_null
pure subroutine nullify_local_zone_descriptors(lzd)
  implicit none
  type(local_zone_descriptors), intent(out) :: lzd

  lzd%linear=.false.
  lzd%nlr=0
  lzd%lintyp=0
  lzd%ndimpotisf=0
  lzd%hgrids=0.0_gp
  call nullify_locreg_descriptors(lzd%glr)
  nullify(lzd%llr) 
end subroutine nullify_local_zone_descriptors

subroutine bigdft_init_errors()
  use dictionaries
  implicit none
  external :: bigdft_severe_abort

  call f_err_define('BIGDFT_RUNTIME_ERROR',&
       'An invalid operation has been done during runtime',&
       BIGDFT_RUNTIME_ERROR,&
       err_action='Check the exact unrolling of runtime operations,'//&
       ' likely something has been initialized/finalized twice')

  call f_err_define('BIGDFT_MPI_ERROR',&
       'An error of MPI library occurred',&
       BIGDFT_MPI_ERROR,&
       err_action='Check if the error is related to MPI library or runtime condtions')

    call f_err_define('BIGDFT_LINALG_ERROR',&
       'An error of linear algebra occurred',&
       BIGDFT_LINALG_ERROR,&
       err_action='Check if the matrix is correct at input, also look at the info value')

    call f_err_define('BIGDFT_INPUT_VARIABLES_ERROR',&
       'An error while parsing the input variables occured',&
       BIGDFT_INPUT_VARIABLES_ERROR,&
       err_action='Check above which input variable has been not correctly parsed, or check their values')

  !define the severe operation via MPI_ABORT
  call f_err_severe_override(bigdft_severe_abort)
end subroutine bigdft_init_errors

!> initialize the timing categories for BigDFT runs.
!! It is of course assumed that f_lib_initialize has already been called
subroutine bigdft_init_timing_categories()
  use Poisson_Solver, only: PS_initialize_timing_categories
  implicit none
  !local variables
  integer :: icls,icat

  !initialize categories for the Poisson Solver
  call PS_initialize_timing_categories()

  !initialize groups
  call f_timing_category_group(tgrp_pot,'Operations for local potential construction (mainly XC)')

  do icls=2,ncls_max
     call f_timing_category_group(trim(clss(icls)),'Empty description for the moment')
  end do

  !define the timing categories for exchange and correlation
  call f_timing_category('Exchange-Correlation',tgrp_pot,&
       'Operations needed to construct local XC potential',&
       TCAT_EXCHANGECORR)


  !! little by little, these categories should be transformed in the 
  !! new scheme dictated by f_timing API in time_profiling module of f_lib.

  !initialize categories
  do icat=1,ncat_bigdft
     call f_timing_category(trim(cats(1,icat)),trim(cats(2,icat)),trim(cats(3,icat)),&
          cat_ids(icat))
  end do

end subroutine bigdft_init_timing_categories

!> routine to convert timing categories from the old scheme to the API of f_lib
!! as soon as the timing categories are identified with their ID, this routine should disappear
subroutine find_category(category,cat_id)
  use yaml_output, only: yaml_warning
  implicit none
  character(len=*), intent(in) :: category
  integer, intent(out) :: cat_id !< id of the found category
  !local variables
  integer :: i
  !controls if the category exists
  cat_id=0
  do i=1,ncat_bigdft
     if (trim(category) == trim(cats(1,i))) then
        cat_id=cat_ids(i)
        exit
     endif
  enddo
  if (cat_id==0) then
     call f_err_throw('Timing routine error,'//&
          ' requested category '//trim(category)//' has not been found',&
          err_id=TIMING_INVALID)
!!$     if (bigdft_mpi%iproc==0) &
!!$          call yaml_warning('Requested timing category '//trim(category)//&
!!$          ' has not been found')
     cat_id=TIMING_UNINITIALIZED
  end if
end subroutine find_category


!!integer function matrixindex_in_compressed(this, iorb, jorb)
!!  implicit none
!!
!!  ! Calling arguments
!!  class(sparse_matrix),intent(in) :: this
!!  integer,intent(in) :: iorb, jorb
!!
!!  ! Local variables
!!  integer :: compressed_index
!!
!!  if (this%store_index) then
!!      ! Take the value from the array
!!      matrixindex_in_compressed = this%matrixindex_in_compressed_arr(iorb,jorb)
!!  else
!!      ! Recalculate the value
!!      matrixindex_in_compressed = compressed_index_fn(iorb, jorb, this%full_dim1, this)
!!  end if
!!
!!  contains
!!    ! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
!!    integer function compressed_index_fn(irow, jcol, norb, sparsemat)
!!      !use module_base
!!      !use module_types
!!      implicit none
!!    
!!      ! Calling arguments
!!      integer,intent(in) :: irow, jcol, norb
!!      type(sparse_matrix),intent(in) :: sparsemat
!!    
!!      ! Local variables
!!      integer :: ii, iseg
!!    
!!      ii=(jcol-1)*norb+irow
!!    
!!      iseg=sparsemat%istsegline(jcol)
!!      do
!!          if (ii>=sparsemat%keyg(1,iseg) .and. ii<=sparsemat%keyg(2,iseg)) then
!!              ! The matrix element is in this segment
!!               compressed_index_fn = sparsemat%keyv(iseg) + ii - sparsemat%keyg(1,iseg)
!!              return
!!          end if
!!          iseg=iseg+1
!!          if (iseg>sparsemat%nseg) exit
!!          if (ii<sparsemat%keyg(1,iseg)) then
!!              compressed_index_fn=0
!!              return
!!          end if
!!      end do
!!    
!!      ! Not found
!!      compressed_index_fn=0
!!    
!!    end function compressed_index_fn
!!
!!end function matrixindex_in_compressed

  subroutine input_set_int(in, key, val)
    use dictionaries, only: dictionary, operator(//), set, dict_init, dict_free
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    integer, intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_int

  subroutine input_set_dbl(in, key, val)
    use dictionaries, only: dictionary, operator(//), set, dict_init, dict_free
    use module_defs, only: gp
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    real(gp), intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_dbl

  subroutine input_set_bool(in, key, val)
    use dictionaries, only: dictionary, operator(//), set, dict_init, dict_free
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    logical, intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_bool

  subroutine input_set_char(in, key, val)
    use dictionaries, only: dictionary, operator(//), set, dict_init, dict_free
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    character(len = *), intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_char

  subroutine input_set_int_array(in, key, val)
    use dictionaries, only: dictionary, operator(//), set, dict_init, dict_free
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    integer :: i
    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    do i = 0, size(val) - 1, 1
       call set(dict // key(sep + 1:) // i, val(i + 1))
    end do
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_int_array

  subroutine input_set_dbl_array(in, key, val)
    use dictionaries, only: dictionary, operator(//), set, dict_init, dict_free
    use module_defs, only: gp
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    real(gp), dimension(:), intent(in) :: val

    integer :: i
    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    do i = 0, size(val) - 1, 1
       call set(dict // key(sep + 1:) // i, val(i + 1))
    end do
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_dbl_array

  subroutine input_set_bool_array(in, key, val)
    use dictionaries, only: dictionary, operator(//), set, dict_init, dict_free
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    logical, dimension(:), intent(in) :: val

    integer :: i
    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    do i = 0, size(val) - 1, 1
       call set(dict // key(sep + 1:) // i, val(i + 1))
    end do
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_bool_array

  subroutine input_set_dict(in, level, val)
    use dictionaries, only: dictionary, operator(//), assignment(=)
    use dictionaries, only: dict_key, max_field_length, dict_value, dict_len
    use module_defs, only: DistProjApply, GPUblas, gp
    use module_input_keys
    use dynamic_memory
    use yaml_output, only: yaml_warning
    implicit none
    type(input_variables), intent(inout) :: in
    type(dictionary), pointer :: val
    character(len = *), intent(in) :: level
    integer, dimension(2) :: dummy_int !<to use as filling for input variables
    real(gp), dimension(3) :: dummy_gp !< to fill the input variables
    character(len = max_field_length) :: str
    integer :: i, ipos

    if (index(dict_key(val), "_attributes") > 0) return

    select case(trim(level))
    case (DFT_VARIABLES)
       ! the DFT variables ------------------------------------------------------
       select case (trim(dict_key(val)))
       case (HGRIDS)
          !grid spacings (profiles can be used if we already read PSPs)
          dummy_gp(1:3)=val
          in%hx = dummy_gp(1)
          in%hy = dummy_gp(2)
          in%hz = dummy_gp(3)
       case (RMULT)
          !coarse and fine radii around atoms
          dummy_gp(1:2)=val
          in%crmult = dummy_gp(1)
          in%frmult = dummy_gp(2)
       case (IXC)
          in%ixc = val !XC functional (ABINIT XC codes)
       case (NCHARGE)
          in%ncharge = val !charge and electric field
       case (ELECFIELD)
          in%elecfield = val
       case (NSPIN)
          in%nspin = val !spin and polarization
       case (MPOL)
          in%mpol = val
       case (GNRM_CV)
          in%gnrm_cv = val !convergence parameters
       case (ITERMAX)
          in%itermax = val
       case (NREPMAX)
          in%nrepmax = val
       case (NCONG)
          in%ncong = val !convergence parameters
       case (IDSX)
          in%idsx = val
       case (DISPERSION)
          in%dispersion = val !dispersion parameter
          ! Now the variables which are to be used only for the last run
       case (INPUTPSIID)
          in%inputPsiId = val
       case (OUTPUT_WF)
          in%output_wf_format = val
       case (OUTPUT_DENSPOT)
          in%output_denspot = val
       case (RBUF)
          in%rbuf = val ! Tail treatment.
       case (NCONGT)
          in%ncongt = val
       case (NORBV)
          in%norbv = val !davidson treatment
       case (NVIRT)
          in%nvirt = val
       case (NPLOT)
          in%nplot = val
       case (DISABLE_SYM)
          in%disableSym = val ! Line to disable symmetries.
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
       ! the KPT variables ------------------------------------------------------
    case (KPT_VARIABLES)
       stop "kpt set_input not implemented"
    case (PERF_VARIABLES)
       ! the PERF variables -----------------------------------------------------
       select case (trim(dict_key(val)))       
       case (DEBUG)
          in%debug = val
       case (FFTCACHE)
          in%ncache_fft = val
       case (VERBOSITY)
          in%verbosity = val
       case (OUTDIR)
          in%writing_directory = val
       case (TOLSYM)
          in%symTol = val
       case (PROJRAD)
          in%projrad = val
       case (EXCTXPAR)
          in%exctxpar = val
       case (INGUESS_GEOPT)
          in%inguess_geopt = val
       case (ACCEL)
          str = dict_value(val)
          if (input_keys_equal(trim(str), "CUDAGPU")) then
             in%matacc%iacceleration = 1
          else if (input_keys_equal(trim(str), "OCLGPU")) then
             in%matacc%iacceleration = 2
          else if (input_keys_equal(trim(str), "OCLCPU")) then
             in%matacc%iacceleration = 3
          else if (input_keys_equal(trim(str), "OCLACC")) then
             in%matacc%iacceleration = 4
          else 
             in%matacc%iacceleration = 0
          end if
       case (OCL_PLATFORM)
          !determine desired OCL platform which is used for acceleration
          in%matacc%OCL_platform = val
          ipos=min(len(in%matacc%OCL_platform),len(trim(in%matacc%OCL_platform))+1)
          do i=ipos,len(in%matacc%OCL_platform)
             in%matacc%OCL_platform(i:i)=achar(0)
          end do
       case (OCL_DEVICES)
          in%matacc%OCL_devices = val
          ipos=min(len(in%matacc%OCL_devices),len(trim(in%matacc%OCL_devices))+1)
          do i=ipos,len(in%matacc%OCL_devices)
             in%matacc%OCL_devices(i:i)=achar(0)
          end do
       case (PSOLVER_ACCEL)
          in%matacc%PSolver_igpu = val
       case (SIGNALING)
          in%signaling = val ! Signaling parameters
       case (SIGNALTIMEOUT)
          in%signalTimeout = val
       case (DOMAIN)
          in%domain = val
       case (BLAS)
          GPUblas = val !!@TODO to relocate
       case (PSP_ONFLY)
          DistProjApply = val
       case (IG_DIAG)
          in%orthpar%directDiag = val
       case (IG_NORBP)
          in%orthpar%norbpInguess = val
       case (IG_TOL)
          in%orthpar%iguessTol = val
       case (METHORTHO)
          in%orthpar%methOrtho = val
       case (IG_BLOCKS)
          !Block size used for the orthonormalization
          dummy_int(1:2)=val
          in%orthpar%bsLow = dummy_int(1)
          in%orthpar%bsUp  = dummy_int(2)
       case (RHO_COMMUN)
          in%rho_commun = val
       case (PSOLVER_GROUPSIZE)
          in%PSolver_groupsize = val
       case (UNBLOCK_COMMS)
          in%unblock_comms = val
       case (LINEAR)
          !Use Linear scaling methods
          str = dict_value(val)
          if (input_keys_equal(trim(str), "LIG")) then
             in%linear = INPUT_IG_LIG
          else if (input_keys_equal(trim(str), "FUL")) then
             in%linear = INPUT_IG_FULL
          else if (input_keys_equal(trim(str), "TMO")) then
             in%linear = INPUT_IG_TMO
          else
             in%linear = INPUT_IG_OFF
          end if
       case (STORE_INDEX)
          in%store_index = val
       case (PDSYEV_BLOCKSIZE)
          !block size for pdsyev/pdsygv, pdgemm (negative -> sequential)
          in%lin%blocksize_pdsyev = val
       case (PDGEMM_BLOCKSIZE)
          in%lin%blocksize_pdgemm = val
       case (MAXPROC_PDSYEV)
          !max number of process uses for pdsyev/pdsygv, pdgemm
          in%lin%nproc_pdsyev = val
       case (MAXPROC_PDGEMM)
          in%lin%nproc_pdgemm = val
       case (EF_INTERPOL_DET)
          !FOE: if the determinant of the interpolation matrix to find the Fermi energy
          !is smaller than this value, switch from cubic to linear interpolation.
          in%lin%ef_interpol_det = val
       case (EF_INTERPOL_CHARGEDIFF)
          in%lin%ef_interpol_chargediff = val
          !determines whether a mixing step shall be preformed after the input guess !(linear version)
       case (MIXING_AFTER_INPUTGUESS)
          in%lin%mixing_after_inputguess = val
          !determines whether the input guess support functions are orthogonalized iteratively (T) or in the standard way (F)
       case (ITERATIVE_ORTHOGONALIZATION)
          in%lin%iterative_orthogonalization = val
       case (CHECK_SUMRHO)
          in%check_sumrho = val
          !  call input_var("mpi_groupsize",0, "number of MPI processes for BigDFT run (0=nproc)", in%mpi_groupsize)
       case (EXPERIMENTAL_MODE)
          in%experimental_mode = val
       case (WRITE_ORBITALS)
          ! linear scaling: write KS orbitals for cubic restart
          in%write_orbitals = val
       case (EXPLICIT_LOCREGCENTERS)
          ! linear scaling: explicitely specify localization centers
          in%explicit_locregcenters = val
       case (CALCULATE_KS_RESIDUE)
          ! linear scaling: calculate Kohn-Sham residue
          in%calculate_KS_residue = val
       case (INTERMEDIATE_FORCES)
          ! linear scaling: calculate intermediate forces
          in%intermediate_forces = val
       case (KAPPA_CONV)
          ! linear scaling: exit kappa for extended input guess (experimental mode)
          in%kappa_conv = val
       case (EVBOUNDS_NSATUR)
           ! linear scaling: number of FOE cycles before the eigenvalue bounds are shrinked
           in%evbounds_nsatur = val
       case(EVBOUNDSSHRINK_NSATUR)
           ! linear scaling: maximal number of unsuccessful eigenvalue bounds shrinkings
           in%evboundsshrink_nsatur = val
       case (METHOD_UPDATEKERNEL)
           ! linear scaling: how to update the density kernel during the support function optimization (0: purification, 1: FOE)
           in%method_updatekernel = val
       case (PURIFICATION_QUICKRETURN)
           ! linear scaling: quick return in purification
           in%purification_quickreturn = val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case ("geopt")

       ! the GEOPT variables ----------------------------------------------------
       select case (trim(dict_key(val)))
       case (GEOPT_METHOD)
          in%geopt_approach = val !geometry input parameters
       case (NCOUNT_CLUSTER_X)
          in%ncount_cluster_x = val
       case (FRAC_FLUCT)
          in%frac_fluct = val
       case (FORCEMAX)
          in%forcemax = val
       case (RANDDIS)
          in%randdis = val
       case (IONMOV)
          in%ionmov = val
       case (DTION)
          in%dtion = val
       case (MDITEMP)
          in%mditemp = val
       case (MDFTEMP)
          in%mdftemp = val
       case (NOSEINERT)
          in%noseinert = val
       case (FRICTION)
          in%friction = val
       case (MDWALL)
          in%mdwall = val
       case (QMASS)
          in%nnos = dict_len(val)
          call f_free_ptr(in%qmass)
          in%qmass = f_malloc_ptr(in%nnos, id = "in%qmass")
          do i=1,in%nnos-1
             in%qmass(i) = dict_len(val // (i-1))
          end do
       case (BMASS)
          in%bmass = val
       case (VMASS)
          in%vmass = val
       case (BETAX)
          in%betax = val
       case (HISTORY)
          in%history = val
       case (DTINIT)
          in%dtinit = val
       case (DTMAX)
          in%dtmax = val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (MIX_VARIABLES)
       ! the MIX variables ------------------------------------------------------
       select case (trim(dict_key(val)))
       case (ISCF)
          in%iscf = val
          !put the startmix if the mixing has to be done
          if (in%iscf >  SCF_KIND_DIRECT_MINIMIZATION) in%gnrm_startmix=1.e300_gp
       case (ITRPMAX)
          in%itrpmax = val
       case (RPNRM_CV)
          in%rpnrm_cv = val
       case (NORBSEMPTY)
          in%norbsempty = val
       case (TEL)
          in%Tel = val
       case (OCCOPT)
          in%occopt = val
       case (ALPHAMIX)
          in%alphamix = val
       case (ALPHADIIS)
          in%alphadiis = val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (SIC_VARIABLES)
       ! the SIC variables ------------------------------------------------------
       select case (trim(dict_key(val)))
       case (SIC_APPROACH)
          in%SIC%approach = val
       case (SIC_ALPHA)
          in%SIC%alpha = val
       case (SIC_FREF)
          in%SIC%fref = val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (TDDFT_VARIABLES)
       ! the TDDFT variables ----------------------------------------------------
       select case (trim(dict_key(val)))
       case (TDDFT_APPROACH)
          in%tddft_approach = val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select      
    case DEFAULT
       call yaml_warning("unknown level '" // trim(level) //"'")
    end select
  END SUBROUTINE input_set_dict

  subroutine basis_params_set_dict(dict_basis,lin,jtype)
    use module_input_keys
    use dictionaries
    implicit none
    integer, intent(in) :: jtype !< local type of which we are filling the values
    type(dictionary), pointer :: dict_basis
    type(linearInputParameters),intent(inout) :: lin
    !local variables
    real(gp), dimension(2) :: dummy_darr
    !-- default parameters for the basis set of linear scaling

    !then update the values of each parameter if present
    select case(trim(dict_key(dict_basis)))
    case(NBASIS)
       lin%norbsPerType(jtype)=dict_basis !npt
    case(AO_CONFINEMENT)
       lin%potentialPrefac_ao(jtype)=dict_basis !ppao
    case(CONFINEMENT)
       dummy_darr=dict_basis
       lin%potentialPrefac_lowaccuracy(jtype)=dummy_darr(1)!ppl
       lin%potentialPrefac_highaccuracy(jtype)=dummy_darr(2)!pph
    case(RLOC)
       dummy_darr=dict_basis
       !locradType(jtype)=dummy_darr(1) !lrl
       lin%locrad_type(jtype,1)=dummy_darr(1) !lrl
       lin%locrad_type(jtype,2)=dummy_darr(2) !lrh
       !locradType_lowaccur(jtype)=dummy_darr(1) !lrl
       !locradType_highaccur(jtype)=dummy_darr(2) !lrh
       !atoms%rloc(jtype,:)=locradType(jtype)
    case(RLOC_KERNEL) 
         lin%kernel_cutoff(jtype)=dict_basis !kco
    case(RLOC_KERNEL_FOE) 
       lin%kernel_cutoff_FOE(jtype)=dict_basis !kco_FOE
    end select
    
  end subroutine basis_params_set_dict

end module module_types
