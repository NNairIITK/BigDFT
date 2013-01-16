!> @file
!!  Define the fortran types
!! @author
!!    Copyright (C) 2008-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Modules which contains the Fortran data structures
!! and the routines of allocations and de-allocations
module module_types

  use m_ab6_mixing, only : ab6_mixing_object
  use module_base, only : gp,wp,dp,tp,uninitialized,mpi_environment,mpi_environment_null,&
       bigdft_mpi,ndebug,memocc,vcopy
  use gaussians, only: gaussian_basis
  implicit none

  !> Constants to determine between cubic version and linear version
  integer,parameter :: CUBIC_VERSION =  0
  integer,parameter :: LINEAR_VERSION = 100

  !> Error codes, to be documented little by little
  integer, parameter :: BIGDFT_SUCCESS        = 0   !< No errors
  integer, parameter :: BIGDFT_UNINITIALIZED  = -10 !< The quantities we want to access seem not yet defined
  integer, parameter :: BIGDFT_INCONSISTENCY  = -11 !< Some of the quantities is not correct


  !> Input wf parameters.
  integer, parameter :: INPUT_PSI_EMPTY        = -1000
  integer, parameter :: INPUT_PSI_RANDOM       = -2
  integer, parameter :: INPUT_PSI_CP2K         = -1
  integer, parameter :: INPUT_PSI_LCAO         = 0
  integer, parameter :: INPUT_PSI_MEMORY_WVL   = 1
  integer, parameter :: INPUT_PSI_DISK_WVL     = 2
  integer, parameter :: INPUT_PSI_LCAO_GAUSS   = 10
  integer, parameter :: INPUT_PSI_MEMORY_GAUSS = 11
  integer, parameter :: INPUT_PSI_DISK_GAUSS   = 12
  integer, parameter :: INPUT_PSI_LINEAR_AO    = 100
  integer, parameter :: INPUT_PSI_MEMORY_LINEAR= 101
  integer, parameter :: INPUT_PSI_DISK_LINEAR  = 102

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
       (/ "none        ", "plain text  ", "Fortran bin.", "ETSF        " /)

  !> Output grid parameters.
  integer, parameter :: OUTPUT_DENSPOT_NONE    = 0
  integer, parameter :: OUTPUT_DENSPOT_DENSITY = 1
  integer, parameter :: OUTPUT_DENSPOT_DENSPOT = 2
  character(len = 12), dimension(0:2), parameter :: OUTPUT_DENSPOT_names = &
       (/ "none        ", "density     ", "dens. + pot." /)
  integer, parameter :: OUTPUT_DENSPOT_FORMAT_TEXT = 0
  integer, parameter :: OUTPUT_DENSPOT_FORMAT_ETSF = 1
  integer, parameter :: OUTPUT_DENSPOT_FORMAT_CUBE = 2
  character(len = 4), dimension(0:2), parameter :: OUTPUT_DENSPOT_format_names = &
       (/ "text", "ETSF", "cube" /)

  !> SCF mixing parameters. (mixing parameters to be added)
  integer, parameter :: SCF_KIND_GENERALIZED_DIRMIN = -1
  integer, parameter :: SCF_KIND_DIRECT_MINIMIZATION = 0

  !> Occupation parameters.
  integer, parameter :: SMEARING_DIST_ERF   = 1
  integer, parameter :: SMEARING_DIST_FERMI = 2
  integer, parameter :: SMEARING_DIST_COLD1 = 3  !< Marzari's cold smearing  with a=-.5634 (bumb minimization)
  integer, parameter :: SMEARING_DIST_COLD2 = 4  !< Marzari's cold smearing  with a=-.8165 (monotonic tail)
  integer, parameter :: SMEARING_DIST_METPX = 5  !< Methfessel and Paxton (same as COLD with a=0)
  character(len = 11), dimension(5), parameter :: smearing_names = &
       (/ "Error func.", &
          "Fermi      ", &
          "Cold (bumb)", &
          "Cold (mono)", &
          "Meth.-Pax. " /)

  !> Target function for the optimization of the basis functions (linear scaling version)
  integer,parameter:: TARGET_FUNCTION_IS_TRACE=0
  integer,parameter:: TARGET_FUNCTION_IS_ENERGY=1
  integer,parameter:: TARGET_FUNCTION_IS_HYBRID=2
  !!integer,parameter:: DECREASE_LINEAR=0
  !!integer,parameter:: DECREASE_ABRUPT=1
  !!integer,parameter:: COMMUNICATION_COLLECTIVE=0
  !!integer,parameter:: COMMUNICATION_P2P=1
  integer,parameter:: LINEAR_DIRECT_MINIMIZATION=100
  integer,parameter:: LINEAR_MIXDENS_SIMPLE=101
  integer,parameter:: LINEAR_MIXPOT_SIMPLE=102
  integer,parameter:: LINEAR_FOE=103
  

  !> Type used for the orthogonalisation parameter
  type, public :: orthon_data
     !> directDiag decides which input guess is chosen:
     !!   if .true. -> as usual direct diagonalization of the Hamiltonian with dsyev (suitable for small systems)
     !!   if .false. -> iterative diagonalization (suitable for large systems)
     logical:: directDiag
     !> norbpInguess indicates how many orbitals shall be treated by each process during the input guess
     !! if directDiag=.false.
     integer:: norbpInguess
     !> You have to choose two numbers for the block size, bsLow and bsUp:
     !!   if bsLow<bsUp, then the program will choose an appropriate block size in between these two numbers
     !!   if bsLow==bsUp, then the program will take exactly this blocksize
     integer:: bsLow, bsUp
     !> the variable methOrtho indicates which orthonormalization procedure is used:
     !!   methOrtho==0 -> Gram-Schmidt with Cholesky decomposition
     !!   methOrtho==1 -> combined block wise classical Gram-Schmidt and Cholesky
     !!   methOrtho==2 -> Loewdin
     integer:: methOrtho
     !> iguessTol gives the tolerance to which the input guess will converged (maximal
     !! residue of all orbitals).
     real(gp):: iguessTol
     integer:: methTransformOverlap, nItOrtho, blocksize_pdsyev, blocksize_pdgemm, nproc_pdsyev
  end type orthon_data

  type, public :: SIC_data
     character(len=4) :: approach !< approach for the Self-Interaction-Correction (PZ, NK)
     integer :: ixc !< base for the SIC correction
     real(gp) :: alpha !<downscaling coefficient
     real(gp) :: fref !< reference occupation (for alphaNK case)
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

  !> Contains all parameters related to the linear scaling version.
  type,public:: linearInputParameters 
    integer :: DIIS_hist_lowaccur, DIIS_hist_highaccur, nItPrecond
    integer :: nItSCCWhenOptimizing, nItBasis_lowaccuracy, nItBasis_highaccuracy
    integer :: mixHist_lowaccuracy, mixHist_highaccuracy
    integer :: methTransformOverlap, blocksize_pdgemm, blocksize_pdsyev
    integer :: correctionOrthoconstraint, nproc_pdsyev, nproc_pdgemm
    integer :: nit_lowaccuracy, nit_highaccuracy
    integer :: nItSCCWhenFixed_lowaccuracy, nItSCCWhenFixed_highaccuracy
    real(kind=8) :: convCrit_lowaccuracy, convCrit_highaccuracy, alphaSD, alphaDIIS
    real(kind=8) :: alpha_mix_lowaccuracy, alpha_mix_highaccuracy, reduce_confinement_factor
    integer :: plotBasisFunctions
    real(kind=8) ::  fscale, deltaenergy_multiplier_TMBexit, deltaenergy_multiplier_TMBfix
    real(kind=8) :: lowaccuracy_conv_crit, convCritMix_lowaccuracy, convCritMix_highaccuracy
    real(kind=8) :: highaccuracy_conv_crit, support_functions_converged
    real(kind=8), dimension(:), pointer :: locrad, locrad_lowaccuracy, locrad_highaccuracy, locrad_type, kernel_cutoff
    real(kind=8), dimension(:), pointer :: potentialPrefac_lowaccuracy, potentialPrefac_highaccuracy
    integer, dimension(:), pointer :: norbsPerType
    integer :: scf_mode, nlevel_accuracy
    logical :: calc_dipole, pulay_correction
  end type linearInputParameters

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
  end type material_acceleration


  !> Structure of the variables read by input.* files (*.dft, *.geopt...)
  type, public :: input_variables
     !strings of the input files
     character(len=100) :: file_dft,file_geopt,file_kpt,file_perf,file_tddft, &
                           file_mix,file_sic,file_occnum,file_igpop,file_lin
     character(len=100) :: dir_output !< Strings of the directory which contains all data output files
     character(len=100) :: run_name   !< Contains the prefix (by default input) used for input files as input.dft
     integer :: files                 !< Existing files.
     !miscellaneous variables
     logical :: gaussian_help
     integer :: ixc,ncharge,itermax,nrepmax,ncong,idsx,ncongt,inputPsiId,nspin,mpol,itrpmax
     integer :: norbv,nvirt,nplot,iscf,norbsempty,norbsuempty,norbsdempty, occopt
     integer :: OUTPUT_DENSPOT,dispersion,last_run,output_wf_format,OUTPUT_DENSPOT_format
     real(gp) :: frac_fluct,gnrm_sw,alphamix,Tel, alphadiis
     real(gp) :: hx,hy,hz,crmult,frmult,gnrm_cv,rbuf,rpnrm_cv,gnrm_startmix
     integer :: verbosity
     real(gp) :: elecfield(3)
     logical :: disableSym

     ! For absorption calculations
     integer :: iabscalc_type   !< 0 non calc, 1 cheb ,  2 lanc
     !! integer :: iat_absorber, L_absorber, N_absorber, rpower_absorber, Linit_absorber
     integer :: iat_absorber,  L_absorber
     real(gp), pointer:: Gabs_coeffs(:)
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
     integer :: nkpt, nkptv,ngroups_kptv
     integer, dimension(:), pointer :: nkptsv_group
     real(gp), pointer :: kpt(:,:), wkpt(:), kptv(:,:)
     character(len=100) :: band_structure_filename

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

     !orthogonalisation data
     type(orthon_data) :: orthpar
  
     !linear scaling data
     type(linearInputParameters) :: lin

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

  !> Bounds for coarse and fine grids for kinetic operations
  !! Useful only for isolated systems AND in CPU
  type, public :: kinetic_bounds
     integer, dimension(:,:,:), pointer :: ibyz_c,ibxz_c,ibxy_c
     integer, dimension(:,:,:), pointer :: ibyz_f,ibxz_f,ibxy_f
  end type kinetic_bounds


  !> Bounds to compress the wavefunctions
  !! Useful only for isolated systems AND in CPU
  type, public :: shrink_bounds
     integer, dimension(:,:,:), pointer :: ibzzx_c,ibyyzz_c
     integer, dimension(:,:,:), pointer :: ibxy_ff,ibzzx_f,ibyyzz_f
  end type shrink_bounds


  !> Bounds to uncompress the wavefunctions
  !! Useful only for isolated systems AND in CPU
  type, public :: grow_bounds
     integer, dimension(:,:,:), pointer :: ibzxx_c,ibxxyy_c
     integer, dimension(:,:,:), pointer :: ibyz_ff,ibzxx_f,ibxxyy_f
  end type grow_bounds


  !> Bounds for convolutions operations
  !! Useful only for isolated systems AND in CPU
  type, public :: convolutions_bounds
     type(kinetic_bounds) :: kb
     type(shrink_bounds) :: sb
     type(grow_bounds) :: gb
     integer, dimension(:,:,:), pointer :: ibyyzz_r !< real space border
  end type convolutions_bounds

  !> Used for lookup table for compressed wavefunctions
  type, public :: wavefunctions_descriptors
     integer :: nvctr_c,nvctr_f,nseg_c,nseg_f
     integer, dimension(:,:), pointer :: keyglob
     integer, dimension(:,:), pointer :: keygloc
     integer, dimension(:), pointer :: keyvloc,keyvglob
  end type wavefunctions_descriptors

  !> Grid dimensions in old different wavelet basis
  type, public :: grid_dimensions
     integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i
  end type grid_dimensions

  !> Contains the information needed for describing completely a
  !! wavefunction localisation region
  type, public :: locreg_descriptors
     character(len=1) :: geocode
     logical :: hybrid_on   !< interesting for global, periodic, localisation regions
     integer :: ns1,ns2,ns3 !< starting point of the localisation region in global coordinates
     integer :: nsi1,nsi2,nsi3  !< starting point of locreg for interpolating grid
     integer :: Localnorb              !< number of orbitals contained in locreg
     integer,dimension(3) :: outofzone  !< vector of points outside of the zone outside Glr for periodic systems
     real(kind=8),dimension(3):: locregCenter !< center of the locreg 
     real(kind=8):: locrad !< cutoff radius of the localization region
     type(grid_dimensions) :: d
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
  end type locreg_descriptors

  !> Non local pseudopotential descriptors
  type, public :: nonlocal_psp_descriptors
     integer :: nproj,nprojel,natoms                  !< Number of projectors and number of elements
     type(locreg_descriptors), dimension(:), pointer :: plr !< pointer which indicates the different localization region per processor
  end type nonlocal_psp_descriptors


  !> Used to split between points to be treated in simple or in double precision
  type, public :: rho_descriptors
     character(len=1) :: geocode
     integer :: icomm !< method for communicating the density
     integer :: nrhotot !< dimension of the partial density array before communication
     integer :: n_csegs,n_fsegs,dp_size,sp_size
     integer, dimension(:,:), pointer :: spkey,dpkey
     integer, dimension(:), pointer :: cseg_b,fseg_b
  end type rho_descriptors

  !> Quantities used for the symmetry operators.
  type, public :: symmetry_data
     integer :: symObj    !< The symmetry object from ABINIT
     integer, dimension(:,:,:), pointer :: irrzon
     real(dp), dimension(:,:,:), pointer :: phnons
  end type symmetry_data

  !> Atomic data (name, polarisation, ...)
  type, public :: atoms_data
     character(len=1) :: geocode
     character(len=5) :: format
     character(len=20) :: units
     integer :: nat                                        !< nat            Number of atoms
     integer :: ntypes                                     !< ntypes         Number of type of atoms
     integer :: natsc
     character(len=20), dimension(:), pointer :: atomnames !< atomnames(ntypes) Name of type of atoms
     real(gp) :: alat1,alat2,alat3                         !< dimension of the periodic supercell
     integer, dimension(:), pointer :: iatype              !< iatype(nat)    Type of the atoms
     integer, dimension(:), pointer :: iasctype
     integer, dimension(:), pointer :: natpol
     integer, dimension(:), pointer :: nelpsp
     integer, dimension(:), pointer :: npspcode
     integer, dimension(:), pointer :: ixcpsp
     integer, dimension(:), pointer :: nzatom
     real(gp), dimension(:,:), pointer :: radii_cf         !< user defined radii_cf, overridden in sysprop.f90
     integer, dimension(:), pointer :: ifrztyp             !< ifrztyp(nat) Frozen atoms
     real(gp), dimension(:), pointer :: amu                !< amu(ntypes)  Atomic Mass Unit for each type of atoms
     real(gp), dimension(:,:), pointer :: aocc,rloc
     real(gp), dimension(:,:,:), pointer :: psppar         !< pseudopotential parameters (HGH SR section)
     logical :: donlcc                                     !< activate non-linear core correction treatment
     integer, dimension(:), pointer :: nlcc_ngv,nlcc_ngc   !<number of valence and core gaussians describing NLCC 
     real(gp), dimension(:,:), pointer :: nlccpar    !< parameters for the non-linear core correction, if present
     real(gp), dimension(:,:), pointer :: ig_nlccpar !< parameters for the input NLCC
     type(symmetry_data) :: sym                      !< The symmetry operators

     !! for abscalc with pawpatch
     integer, dimension(:), pointer ::  paw_NofL, paw_l, paw_nofchannels
     integer, dimension(:), pointer ::  paw_nofgaussians
     real(gp), dimension(:), pointer :: paw_Greal, paw_Gimag, paw_Gcoeffs
     real(gp), dimension(:), pointer :: paw_H_matrices, paw_S_matrices, paw_Sm1_matrices
     integer :: iat_absorber 

  end type atoms_data

  !> Structure to store the density / potential distribution among processors.
  type, public :: denspot_distribution
     integer :: n3d,n3p,n3pi,i3xcsh,i3s,nrhodim,i3rho_add
     integer :: ndimpot,ndimgrid,ndimrhopot 
     integer, dimension(3) :: ndims !< box containing the grid dimensions in ISF basis
     real(gp), dimension(3) :: hgrids !< grid spacings of the box (half of wavelet ones)
     integer, dimension(:,:), pointer :: nscatterarr, ngatherarr
     type(mpi_environment) :: mpi_env
  end type denspot_distribution

  !> Structures of basis of gaussian functions of the form exp(-a*r2)cos/sin(b*r2)
  type, public :: gaussian_basis_c
     integer :: nat,ncoeff,nshltot,nexpo
     integer, dimension(:), pointer :: nshell,ndoc,nam
     complex(gp), dimension(:), pointer :: expof,psiat
     real(gp), dimension(:,:), pointer :: rxyz
  end type gaussian_basis_c

  !> Contains all array necessary to apply preconditioning projectors 
  type, public :: pcproj_data_type
     type(nonlocal_psp_descriptors) :: pc_nlpspd
     real(gp), pointer :: pc_proj(:)
     integer , pointer , dimension(:):: ilr_to_mproj, iproj_to_l
     real(gp) , pointer ::  iproj_to_ene(:)
     real(gp) , pointer ::  iproj_to_factor(:)
     integer, pointer :: iorbtolr(:)
     integer :: mprojtot
     type(gaussian_basis)  :: G          
     real(gp), pointer :: gaenes(:)
     real(gp) :: ecut_pc
     logical :: DistProjApply
  end type pcproj_data_type

  !> Contains all array necessary to apply preconditioning projectors 
  type, public :: PAWproj_data_type
     type(nonlocal_psp_descriptors) :: paw_nlpspd

     integer , pointer , dimension(:):: iproj_to_paw_nchannels
     integer , pointer , dimension(:):: ilr_to_mproj, iproj_to_l
     integer , pointer , dimension(:):: iprojto_imatrixbeg

     integer :: mprojtot
     real(gp), pointer :: paw_proj(:)
     type(gaussian_basis_c)  :: G          
     integer, pointer :: iorbtolr(:)
     logical :: DistProjApply
  end type PAWproj_data_type


  !> All the parameters which are important for describing the orbitals
  !! Add also the objects related to k-points sampling, after symmetries applications
  type, public :: orbitals_data 
     integer :: norb          !< Total number of orbitals per k point
     integer :: norbp         !< Total number of orbitals for the given processors
     integer :: norbu,norbd,nspin,nspinor,isorb
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: nkpts,nkptsp,iskpts
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     real(gp) :: efermi,HLgap, eTS
     integer, dimension(:), pointer :: iokpt,ikptproc,isorb_par,ispot
     integer, dimension(:), pointer :: inwhichlocreg,onWhichMPI,onwhichatom
     integer, dimension(:,:), pointer :: norb_par
     real(wp), dimension(:), pointer :: eval
     real(gp), dimension(:), pointer :: occup,spinsgn,kwgts
     real(gp), dimension(:,:), pointer :: kpts
  end type orbitals_data

  !> Contains the information needed for communicating the wavefunctions
  !! between processors for the transposition
  type, public :: communications_arrays
     integer, dimension(:), pointer :: ncntd,ncntt,ndspld,ndsplt
     integer, dimension(:,:), pointer :: nvctr_par
  end type communications_arrays


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
     integer:: ndimpotisf                      !< total dimension of potential in isf (including exctX)
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


  !> Contains the arguments needed for the application of the hamiltonian
  type, public :: lanczos_args
     !arguments for the hamiltonian
     integer :: iproc,nproc,ndimpot,nspin, in_iat_absorber, Labsorber
     real(gp) :: hx,hy,hz
     type(energy_terms) :: energs
     !real(gp) :: ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC
     type(atoms_data), pointer :: at
     type(orbitals_data), pointer :: orbs
     type(communications_arrays) :: comms
     type(nonlocal_psp_descriptors), pointer :: nlpspd
     type(local_zone_descriptors), pointer :: Lzd
     type(gaussian_basis), pointer :: Gabsorber    
     type(SIC_data), pointer :: SIC
     integer, dimension(:,:), pointer :: ngatherarr 
     real(gp), dimension(:,:),  pointer :: rxyz,radii_cf
     real(wp), dimension(:), pointer :: proj
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), pointer :: psi
     real(wp), dimension(:), pointer :: potential
     real(wp), dimension(:), pointer :: Gabs_coeffs
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp) :: hpsi
     type(GPU_pointers), pointer :: GPU
     type(pcproj_data_type), pointer :: PPD
     type(pawproj_data_type), pointer :: PAWD
  end type lanczos_args


  !> Contains all parameters needed for point to point communication
  type,public:: p2pComms
    integer,dimension(:),pointer:: noverlaps
    real(kind=8),dimension(:),pointer:: sendBuf, recvBuf
    integer,dimension(:,:,:),pointer:: comarr
    integer:: nsendBuf, nrecvBuf, window
    integer,dimension(:,:),pointer:: ise ! starting / ending index of recvBuf in x,y,z dimension after communication (glocal coordinates)
    integer,dimension(:,:),pointer:: mpi_datatypes
    logical:: communication_complete
  end type p2pComms

  !! Contains the parameters for calculating the overlap matrix for the orthonormalization etc...
  type,public:: overlapParameters
      integer,dimension(:),pointer:: noverlaps
      integer,dimension(:,:),pointer:: overlaps
  end type overlapParameters

  type,public:: matrixDescriptors
      integer:: nvctr, nseg
      integer,dimension(:),pointer:: keyv, nsegline, istsegline
      integer,dimension(:,:),pointer:: keyg
      integer,dimension(:),pointer :: kernel_nseg
      integer,dimension(:,:,:),pointer :: kernel_segkeyg
  end type matrixDescriptors


  type:: collective_comms
    integer:: nptsp_c, ndimpsi_c, ndimind_c, ndimind_f, nptsp_f, ndimpsi_f
    integer,dimension(:),pointer:: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
    integer,dimension(:),pointer:: isendbuf_c, iextract_c, iexpand_c, irecvbuf_c
    integer,dimension(:),pointer:: norb_per_gridpoint_c, indexrecvorbital_c
    integer,dimension(:),pointer:: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
    integer,dimension(:),pointer:: isendbuf_f, iextract_f, iexpand_f, irecvbuf_f
    integer,dimension(:),pointer:: norb_per_gridpoint_f, indexrecvorbital_f
    integer,dimension(:),pointer:: isptsp_c, isptsp_f !<starting index of a given gridpoint (basically summation of norb_per_gridpoint_*)
    real(kind=8),dimension(:),pointer :: psit_c, psit_f
    integer,dimension(:),pointer :: nsendcounts_repartitionrho, nrecvcounts_repartitionrho
    integer,dimension(:),pointer :: nsenddspls_repartitionrho, nrecvdspls_repartitionrho
    integer,dimension(:,:),pointer :: matrixindex_in_compressed
  end type collective_comms


  type,public:: workarrays_quartic_convolutions
    real(wp),dimension(:,:,:),pointer:: xx_c, xy_c, xz_c
    real(wp),dimension(:,:,:),pointer:: xx_f1
    real(wp),dimension(:,:,:),pointer:: xy_f2
    real(wp),dimension(:,:,:),pointer:: xz_f4
    real(wp),dimension(:,:,:,:),pointer:: xx_f, xy_f, xz_f
    real(wp),dimension(:,:,:),pointer:: y_c
    real(wp),dimension(:,:,:,:),pointer:: y_f
    ! The following arrays are work arrays within the subroutine
    real(wp),dimension(:,:),pointer:: aeff0array, beff0array, ceff0array, eeff0array
    real(wp),dimension(:,:),pointer:: aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array
    real(wp),dimension(:,:),pointer:: aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray
    real(wp),dimension(:,:,:),pointer:: xya_c, xyc_c
    real(wp),dimension(:,:,:),pointer:: xza_c, xzc_c
    real(wp),dimension(:,:,:),pointer:: yza_c, yzb_c, yzc_c, yze_c
    real(wp),dimension(:,:,:,:),pointer:: xya_f, xyb_f, xyc_f, xye_f
    real(wp),dimension(:,:,:,:),pointer:: xza_f, xzb_f, xzc_f, xze_f
    real(wp),dimension(:,:,:,:),pointer:: yza_f, yzb_f, yzc_f, yze_f
    real(wp),dimension(-17:17) :: aeff0, aeff1, aeff2, aeff3
    real(wp),dimension(-17:17) :: beff0, beff1, beff2, beff3
    real(wp),dimension(-17:17) :: ceff0, ceff1, ceff2, ceff3
    real(wp),dimension(-14:14) :: eeff0, eeff1, eeff2, eeff3
    real(wp),dimension(-17:17) :: aeff0_2, aeff1_2, aeff2_2, aeff3_2
    real(wp),dimension(-17:17) :: beff0_2, beff1_2, beff2_2, beff3_2
    real(wp),dimension(-17:17) :: ceff0_2, ceff1_2, ceff2_2, ceff3_2
    real(wp),dimension(-14:14) :: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  end type workarrays_quartic_convolutions
  

  type,public:: localizedDIISParameters
    integer:: is, isx, mis, DIISHistMax, DIISHistMin
    integer:: icountSDSatur, icountDIISFailureCons, icountSwitch, icountDIISFailureTot, itBest
    real(kind=8),dimension(:),pointer:: phiHist, hphiHist
    real(kind=8),dimension(:),pointer:: alpha_coeff !step size for optimization of coefficients
    real(kind=8),dimension(:,:),pointer:: grad_coeff_old !coefficients gradient of previous iteration
    real(kind=8),dimension(:,:,:),pointer:: mat
    real(kind=8):: trmin, trold, alphaSD, alphaDIIS
    logical:: switchSD, immediateSwitchToSD, resetDIIS
  end type localizedDIISParameters


  type,public:: mixrhopotDIISParameters
    integer:: is, isx, mis
    real(kind=8),dimension(:),pointer:: rhopotHist, rhopotresHist
    real(kind=8),dimension(:,:),pointer:: mat
  end type mixrhopotDIISParameters


  type,public:: wfn_metadata
    real(kind=8),dimension(:,:),pointer:: coeff !<expansion coefficients
    real(kind=8),dimension(:,:),pointer:: coeffp !<coefficients distributed over processes
    real(kind=8),dimension(:,:),pointer:: density_kernel !<density kernel
    real(8),dimension(:),pointer :: density_kernel_compr !<compressed density kernel
    real(kind=8) :: ef !< Fermi energy for FOE
    real(kind=8) :: evlow, evhigh !< eigenvalue bounds for FOE 
    real(kind=8) :: bisection_shift !< bisection shift to find Fermi energy (FOE)
    real(kind=8) :: fscale !< length scale for complementary error function (FOE)
  end type wfn_metadata


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
    real(kind=8), dimension(:), pointer :: potentialPrefac !< Prefactor for the potential: Prefac * f(r) 
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

    !> Defines the fundamental structure for the kernel
  type, public :: coulomb_operator
     !variables with physical meaning
     integer :: itype_scf !< order of the ISF family to be used
     real(gp) :: mu !< inverse screening length for the Helmholtz Eq. (Poisson Eq. -> mu=0)
     character(len=1) :: geocode !< Code for the boundary conditions
     integer, dimension(3) :: ndims !< dimension of the box of the density
     real(gp), dimension(3) :: hgrids !<grid spacings in each direction
     real(gp), dimension(3) :: angrad !< angles in radiants between each of the axis
     real(dp), dimension(:), pointer :: kernel !< kernel of the Poisson Solver
     real(dp) :: work1_GPU,work2_GPU,k_GPU
     integer, dimension(5) :: plan
     integer, dimension(3) :: geo
     !variables with computational meaning
!     integer :: iproc_world !iproc in the general communicator
!     integer :: iproc,nproc
!     integer :: mpi_comm
     type(mpi_environment) :: mpi_env
     integer :: igpu !< control the usage of the GPU
  end type coulomb_operator

  !> Densities and potentials, and related metadata, needed for their creation/application
  !! Not all these quantities are available, some of them may point to the same memory space
  type, public :: DFT_local_fields
     real(dp), dimension(:), pointer :: rhov !< generic workspace. What is there is indicated by rhov_is
     
     type(ab6_mixing_object), pointer :: mix          !< History of rhov, allocated only when using diagonalisation
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
  integer,parameter,public :: LINEAR_LOWACCURACY  = 101 !low accuracy after restart
  integer,parameter,public :: LINEAR_HIGHACCURACY = 102 !high accuracy after restart

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
     logical:: can_use_transposed !< true if the transposed quantities are allocated and can be used
     type(orbitals_data) :: orbs !<wavefunction specification in terms of orbitals
     type(communications_arrays) :: comms !< communication objects for the cubic approach
     	type(diis_objects) :: diis
     type(confpot_data), dimension(:), pointer :: confdatarr !<data for the confinement potential
     type(SIC_data) :: SIC !<control the activation of SIC scheme in the wavefunction
     type(orthon_data) :: orthpar !< control the application of the orthogonality scheme for cubic DFT wavefunction
     character(len=4) :: exctxpar !< Method for exact exchange parallelisation for the wavefunctions, in case
     	type(wfn_metadata) :: wfnmd !<specifications of the kind of wavefunction
     type(p2pComms):: comon !<describing p2p communications for orthonormality
     type(overlapParameters):: op !<describing the overlaps
     	type(p2pComms):: comgp !<describing p2p communications for distributing the potential
     type(p2pComms):: comrp !<describing the repartition of the orbitals (for derivatives)
     type(p2pComms):: comsr !<describing the p2p communications for sumrho
     	type(matrixDescriptors):: mad !<describes the structure of the matrices
     type(collective_comms):: collcom ! describes collective communication
     type(collective_comms):: collcom_sr ! describes collective communication for the calculation of the charge density
     integer(kind = 8) :: c_obj !< Storage of the C wrapper object. it has to be initialized to zero
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
     integer :: n1,n2,n3
     real(gp) :: hx_old,hy_old,hz_old
     real(gp), dimension(:,:), pointer :: rxyz_old,rxyz_new
     type(DFT_wavefunction) :: KSwfn !< Kohn-Sham wavefunctions
     type(DFT_wavefunction) :: tmb !<support functions for linear scaling
     type(GPU_pointers) :: GPU 
  end type restart_objects

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
  end function material_acceleration_null

  function pkernel_null() result(k)
    type(coulomb_operator) :: k
    k%itype_scf=0
    k%geocode='F'
    k%mu=0.0_gp
    k%ndims=(/0,0,0/)
    k%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    k%angrad=(/0.0_gp,0.0_gp,0.0_gp/)
    nullify(k%kernel)
    k%work1_GPU=0.d0
    k%work2_GPU=0.d0
    k%k_GPU=0.d0
    k%plan=(/0,0,0,0,0/)
    k%geo=(/0,0,0/)
    k%mpi_env=mpi_environment_null()
    k%igpu=0
  end function pkernel_null

  function default_grid() result(g)
    type(grid_dimensions) :: g
    g%n1   =0
    g%n2   =0
    g%n3   =0
    g%nfl1 =0
    g%nfu1 =0
    g%nfl2 =0
    g%nfu2 =0
    g%nfl3 =0
    g%nfu3 =0
    g%n1i  =0
    g%n2i  =0
    g%n3i  =0
  end function default_grid

  function default_wfd() result(wfd)
    type(wavefunctions_descriptors) :: wfd
    wfd%nvctr_c=0
    wfd%nvctr_f=0
    wfd%nseg_c=0
    wfd%nseg_f=0
    nullify(wfd%keyglob)
    nullify(wfd%keygloc)
    nullify(wfd%keyvglob)
    nullify(wfd%keyvloc)
  end function default_wfd

  function default_kb() result(kb)
    type(kinetic_bounds) :: kb
    nullify(&
         kb%ibyz_c,kb%ibxz_c,kb%ibxy_c, &
         kb%ibyz_f,kb%ibxz_f,kb%ibxy_f)
  end function default_kb

  function default_sb() result(sb)
    type(shrink_bounds) :: sb
    nullify(sb%ibzzx_c,sb%ibyyzz_c,sb%ibxy_ff,sb%ibzzx_f,sb%ibyyzz_f)
  end function default_sb
  
  function default_gb() result(gb)
    type(grow_bounds) :: gb
    nullify(gb%ibzxx_c,gb%ibxxyy_c,gb%ibyz_ff,gb%ibzxx_f,gb%ibxxyy_f)
  end function default_gb

  function default_bounds() result(b)
    type(convolutions_bounds) :: b
    b%kb=default_kb()
    b%sb=default_sb()
    b%gb=default_gb()
    nullify(b%ibyyzz_r)
  end function default_bounds

  function default_locreg() result(lr)
    type(locreg_descriptors) :: lr

    lr%geocode='F'
    lr%hybrid_on=.false.   
    lr%ns1=0
    lr%ns2=0
    lr%ns3=0 
    lr%nsi1=0
    lr%nsi2=0
    lr%nsi3=0  
    lr%Localnorb=0  
    lr%outofzone=(/0,0,0/) 
    lr%d=default_grid()
    lr%wfd=default_wfd()
    lr%bounds=default_bounds()
    lr%locregCenter=(/0.0_gp,0.0_gp,0.0_gp/) 
    lr%locrad=0 

  end function default_locreg

  function default_lzd() result(lzd)
    type(local_zone_descriptors) :: lzd
    lzd%linear=.false.
    lzd%nlr=0
    lzd%lintyp=0
    lzd%ndimpotisf=0
    lzd%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    lzd%Glr=default_locreg()
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

  !accessors for external programs
  !> Get the number of orbitals of the run in rst
  function bigdft_get_number_of_orbitals(rst,istat) result(norb)
    use module_base
    implicit none
    type(restart_objects), intent(in) :: rst !> BigDFT restart variables. call_bigdft already called
    integer :: norb !> Number of orbitals of run in rst
    integer, intent(out) :: istat
    
    istat=BIGDFT_SUCCESS

    norb=rst%KSwfn%orbs%norb
    if (norb==0) istat = BIGDFT_UNINITIALIZED
    
  end function bigdft_get_number_of_orbitals
  
  !> Fill the array eval with the number of orbitals of the last run
  subroutine bigdft_get_eigenvalues(rst,eval,istat)
    use module_base
    implicit none
    type(restart_objects), intent(in) :: rst !> BigDFT restart variables. call_bigdft already called
    real(gp), dimension(*), intent(out) :: eval !> Buffer for eigenvectors. Should have at least dimension equal to bigdft_get_number_of_orbitals(rst,istat)
    integer, intent(out) :: istat !> Error code
    !local variables
    integer :: norb
    
    norb=bigdft_get_number_of_orbitals(rst,istat)

    if (istat /= BIGDFT_SUCCESS) return

    if (.not. associated(rst%KSwfn%orbs%eval)) then
       istat = BIGDFT_UNINITIALIZED
       return
    end if

    if (product(shape(rst%KSwfn%orbs%eval)) < norb) then
       istat = BIGDFT_INCONSISTENCY
       return
    end if

    call vcopy(norb,rst%KSwfn%orbs%eval(1),1,eval(1),1)

  end subroutine bigdft_get_eigenvalues

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
   

!> De-Allocate communications_arrays
  subroutine deallocate_comms(comms,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(communications_arrays), intent(inout) :: comms
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
    i_all=-product(shape(orbs%onWhichMPI))*kind(orbs%onWhichMPI)
    deallocate(orbs%onWhichMPI,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%onWhichMPI',subname)
    if (associated(orbs%ispot)) then
       i_all=-product(shape(orbs%ispot))*kind(orbs%ispot)
       deallocate(orbs%ispot,stat=i_stat)
       call memocc(i_stat,i_all,'orbs%ispot',subname)
    end if

END SUBROUTINE deallocate_orbs


!> Allocate and nullify restart objects
  subroutine init_restart_objects(iproc,matacc,atoms,rst,subname)
    use module_base
    implicit none
    !Arguments
    character(len=*), intent(in) :: subname
    integer, intent(in) :: iproc
    type(material_acceleration), intent(in) :: matacc
    type(atoms_data), intent(in) :: atoms
    type(restart_objects), intent(out) :: rst
    !local variables
    integer :: i_stat

    !allocate pointers
    allocate(rst%rxyz_new(3,atoms%nat+ndebug),stat=i_stat)
    call memocc(i_stat,rst%rxyz_new,'rxyz_new',subname)
    allocate(rst%rxyz_old(3,atoms%nat+ndebug),stat=i_stat)
    call memocc(i_stat,rst%rxyz_old,'rxyz_old',subname)

    !nullify unallocated pointers
    rst%KSwfn%c_obj = 0
    nullify(rst%KSwfn%psi)
    nullify(rst%KSwfn%orbs%eval)

    nullify(rst%KSwfn%gaucoeffs)
    nullify(rst%KSwfn%oldpsis)

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

    !initialise the acceleration strategy if required
    call init_material_acceleration(iproc,matacc,rst%GPU)

  END SUBROUTINE init_restart_objects


!>  De-Allocate restart_objects
  subroutine free_restart_objects(rst,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(restart_objects) :: rst
    !local variables
    integer :: i_all,i_stat,istep

    call deallocate_locreg_descriptors(rst%KSwfn%Lzd%Glr,subname)

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

    if (associated(rst%rxyz_old)) then
       i_all=-product(shape(rst%rxyz_old))*kind(rst%rxyz_old)
       deallocate(rst%rxyz_old,stat=i_stat)
       call memocc(i_stat,i_all,'rxyz_old',subname)
    end if

    if (associated(rst%KSwfn%oldpsis)) then
       do istep=0,product(shape(rst%KSwfn%oldpsis))-1
          call old_wavefunction_free(rst%KSwfn%oldpsis(istep),subname)
       end do
       deallocate(rst%KSwfn%oldpsis)
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
  subroutine allocate_wfd(wfd,subname)
    use module_base
    implicit none
    type(wavefunctions_descriptors), intent(inout) :: wfd
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_stat

    allocate(wfd%keyglob(2,max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
    call memocc(i_stat,wfd%keyglob,'keyglob',subname)
    allocate(wfd%keygloc(2,max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
    call memocc(i_stat,wfd%keygloc,'keygloc',subname)
    allocate(wfd%keyvloc(max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
    call memocc(i_stat,wfd%keyvloc,'keyvloc',subname)
    allocate(wfd%keyvglob(max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
    call memocc(i_stat,wfd%keyvglob,'keyvglob',subname)

  END SUBROUTINE allocate_wfd


!> De-Allocate wavefunctions_descriptors
  subroutine deallocate_wfd(wfd,subname)
    use module_base
    implicit none
    type(wavefunctions_descriptors) :: wfd
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    if (associated(wfd%keyglob, target = wfd%keygloc)) then
       i_all=-product(shape(wfd%keyglob))*kind(wfd%keyglob)
       deallocate(wfd%keyglob,stat=i_stat)
       call memocc(i_stat,i_all,'wfd%keyglob',subname)
       nullify(wfd%keyglob)
    else
       if(associated(wfd%keyglob)) then
          i_all=-product(shape(wfd%keyglob))*kind(wfd%keyglob)
          deallocate(wfd%keyglob,stat=i_stat)
          call memocc(i_stat,i_all,'wfd%keyglob',subname)
          nullify(wfd%keyglob)
       end if
       if(associated(wfd%keygloc)) then 
          i_all=-product(shape(wfd%keygloc))*kind(wfd%keygloc)
          deallocate(wfd%keygloc,stat=i_stat)
          call memocc(i_stat,i_all,'wfd%keygloc',subname)
          nullify(wfd%keygloc)
       end if
    end if
    if (associated(wfd%keyvloc, target= wfd%keyvglob)) then
       i_all=-product(shape(wfd%keyvloc))*kind(wfd%keyvloc)
       deallocate(wfd%keyvloc,stat=i_stat)
       call memocc(i_stat,i_all,'wfd%keyvloc',subname)
       nullify(wfd%keyvloc)
    else
       if (associated(wfd%keyvloc)) then
          i_all=-product(shape(wfd%keyvloc))*kind(wfd%keyvloc)
          deallocate(wfd%keyvloc,stat=i_stat)
          call memocc(i_stat,i_all,'wfd%keyvloc',subname)
          nullify(wfd%keyvloc)
       end if
       if (associated(wfd%keyvglob)) then
          i_all=-product(shape(wfd%keyvglob))*kind(wfd%keyvglob)
          deallocate(wfd%keyvglob,stat=i_stat)
          call memocc(i_stat,i_all,'wfd%keyvglob',subname)
          nullify(wfd%keyvglob)
       end if
    end if
  END SUBROUTINE deallocate_wfd

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




!>   De-Allocate gaussian_basis type

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






!> De-Allocate convolutions_bounds type, depending of the geocode and the hybrid_on
  subroutine deallocate_bounds(geocode,hybrid_on,bounds,subname)
    use module_base
    implicit none
    character(len=1), intent(in) :: geocode
    logical, intent(in) :: hybrid_on 
    type(convolutions_bounds) :: bounds
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    if ((geocode == 'P' .and. hybrid_on) .or. geocode == 'F') then
       ! Just test the first one...
       if (associated(bounds%kb%ibyz_f)) then
          i_all=-product(shape(bounds%kb%ibyz_f))*kind(bounds%kb%ibyz_f)
          deallocate(bounds%kb%ibyz_f,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%kb%ibyz_f',subname)
          i_all=-product(shape(bounds%kb%ibxz_f))*kind(bounds%kb%ibxz_f)
          deallocate(bounds%kb%ibxz_f,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%kb%ibxz_f',subname)
          i_all=-product(shape(bounds%kb%ibxy_f))*kind(bounds%kb%ibxy_f)
          deallocate(bounds%kb%ibxy_f,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%kb%ibxy_f',subname)

          i_all=-product(shape(bounds%sb%ibxy_ff))*kind(bounds%sb%ibxy_ff)
          deallocate(bounds%sb%ibxy_ff,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%sb%ibxy_ff',subname)
          i_all=-product(shape(bounds%sb%ibzzx_f))*kind(bounds%sb%ibzzx_f)
          deallocate(bounds%sb%ibzzx_f,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%sb%ibzzx_f',subname)
          i_all=-product(shape(bounds%sb%ibyyzz_f))*kind(bounds%sb%ibyyzz_f)
          deallocate(bounds%sb%ibyyzz_f,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%sb%ibyyzz_f',subname)
          i_all=-product(shape(bounds%gb%ibyz_ff))*kind(bounds%gb%ibyz_ff)
          deallocate(bounds%gb%ibyz_ff,stat=i_stat)

          call memocc(i_stat,i_all,'bounds%gb%ibyz_ff',subname)
          i_all=-product(shape(bounds%gb%ibzxx_f))*kind(bounds%gb%ibzxx_f)
          deallocate(bounds%gb%ibzxx_f,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%gb%ibzxx_f',subname)
          i_all=-product(shape(bounds%gb%ibxxyy_f))*kind(bounds%gb%ibxxyy_f)
          deallocate(bounds%gb%ibxxyy_f,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%gb%ibxxyy_f',subname)

          nullify(bounds%kb%ibyz_f)
          nullify(bounds%kb%ibxz_f)
          nullify(bounds%kb%ibxy_f)
          nullify(bounds%sb%ibxy_ff)
          nullify(bounds%sb%ibzzx_f)
          nullify(bounds%sb%ibyyzz_f)
          nullify(bounds%gb%ibyz_ff)
          nullify(bounds%gb%ibzxx_f)
          nullify(bounds%gb%ibxxyy_f)
       end if
    end if

    !the arrays which are needed only for free BC
    if (geocode == 'F') then
       ! Just test the first one...
       if (associated(bounds%kb%ibyz_c)) then
          i_all=-product(shape(bounds%kb%ibyz_c))*kind(bounds%kb%ibyz_c)
          deallocate(bounds%kb%ibyz_c,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%kb%ibyz_c',subname)
          i_all=-product(shape(bounds%kb%ibxz_c))*kind(bounds%kb%ibxz_c)
          deallocate(bounds%kb%ibxz_c,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%kb%ibxz_c',subname)
          i_all=-product(shape(bounds%kb%ibxy_c))*kind(bounds%kb%ibxy_c)
          deallocate(bounds%kb%ibxy_c,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%kb%ibxy_c',subname)
          i_all=-product(shape(bounds%sb%ibzzx_c))*kind(bounds%sb%ibzzx_c)
          deallocate(bounds%sb%ibzzx_c,stat=i_stat)

          call memocc(i_stat,i_all,'bounds%sb%ibzzx_c',subname)
          i_all=-product(shape(bounds%sb%ibyyzz_c))*kind(bounds%sb%ibyyzz_c)
          deallocate(bounds%sb%ibyyzz_c,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%sb%ibyyzz_c',subname)
          i_all=-product(shape(bounds%gb%ibzxx_c))*kind(bounds%gb%ibzxx_c)
          deallocate(bounds%gb%ibzxx_c,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%gb%ibzxx_c',subname)
          i_all=-product(shape(bounds%gb%ibxxyy_c))*kind(bounds%gb%ibxxyy_c)
          deallocate(bounds%gb%ibxxyy_c,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%gb%ibxxyy_c',subname)

          i_all=-product(shape(bounds%ibyyzz_r))*kind(bounds%ibyyzz_r)
          deallocate(bounds%ibyyzz_r,stat=i_stat)
          call memocc(i_stat,i_all,'bounds%ibyyzz_r',subname)

          nullify(bounds%kb%ibyz_c)
          nullify(bounds%kb%ibxz_c)
          nullify(bounds%kb%ibxy_c)
          nullify(bounds%sb%ibzzx_c)
          nullify(bounds%sb%ibyyzz_c)
          nullify(bounds%gb%ibzxx_c)
          nullify(bounds%gb%ibxxyy_c)
          nullify(bounds%ibyyzz_r)
       end if
    end if

  END SUBROUTINE deallocate_bounds


  !> todo: remove this function.
  subroutine deallocate_lr(lr,subname)
    use module_base
    character(len=*), intent(in) :: subname
    type(locreg_descriptors) :: lr
!    integer :: i_all,i_stat

    write(0,*) "TODO, remove me"
    
    call deallocate_wfd(lr%wfd,subname)

    call deallocate_bounds(lr%geocode,lr%hybrid_on,lr%bounds,subname)

!    if (associated(lr%projflg)) then
!       i_all=-product(shape(lr%projflg)*kind(lr%projflg))
!       deallocate(lr%projflg,stat=i_stat)
!       call memocc(i_stat,i_all,'lr%projflg',subname)
!    end if
  END SUBROUTINE deallocate_lr

  subroutine deallocate_symmetry(sym, subname)
    use module_base
    use m_ab6_symmetry
    implicit none
    type(symmetry_data), intent(inout) :: sym
    character(len = *), intent(in) :: subname

    integer :: i_stat, i_all

    if (sym%symObj >= 0) then
       call symmetry_free(sym%symObj)
    end if

    if (associated(sym%irrzon)) then
       i_all=-product(shape(sym%irrzon))*kind(sym%irrzon)
       deallocate(sym%irrzon,stat=i_stat)
       call memocc(i_stat,i_all,'irrzon',subname)
       nullify(sym%irrzon)
    end if

    if (associated(sym%phnons)) then
       i_all=-product(shape(sym%phnons))*kind(sym%phnons)
       deallocate(sym%phnons,stat=i_stat)
       call memocc(i_stat,i_all,'phnons',subname)
       nullify(sym%phnons)
    end if
  end subroutine deallocate_symmetry

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
!!
  subroutine deallocate_pawproj_data(pawproj_data,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(pawproj_data_type), intent(inout) :: pawproj_data
    !local variables
    integer :: i_all,i_stat
    if(associated(pawproj_data%paw_proj)) then

       i_all=-product(shape(  pawproj_data% paw_proj ))*kind( pawproj_data% paw_proj  )
       deallocate(pawproj_data%  paw_proj  ,stat=i_stat)
       call memocc(i_stat,i_all,'paw_proj',subname)

       i_all=-product(shape( pawproj_data%ilr_to_mproj   ))*kind(pawproj_data% ilr_to_mproj   )
       deallocate( pawproj_data% ilr_to_mproj  ,stat=i_stat)
       call memocc(i_stat,i_all,'ilr_to_mproj',subname)

       i_all=-product(shape( pawproj_data% iproj_to_l  ))*kind( pawproj_data% iproj_to_l  )
       deallocate(pawproj_data%  iproj_to_l  ,stat=i_stat)
       call memocc(i_stat,i_all,'iproj_to_l',subname)

       i_all=-product(shape( pawproj_data% iproj_to_paw_nchannels  ))*kind( pawproj_data% iproj_to_paw_nchannels  )
       deallocate(pawproj_data%  iproj_to_paw_nchannels  ,stat=i_stat)
       call memocc(i_stat,i_all,'iproj_to_paw_nchannels',subname)

       i_all=-product(shape( pawproj_data% iprojto_imatrixbeg  ))*kind( pawproj_data% iprojto_imatrixbeg  )
       deallocate(pawproj_data%  iprojto_imatrixbeg  ,stat=i_stat)
       call memocc(i_stat,i_all,'iorbto_imatrixbeg',subname)

       i_all=-product(shape( pawproj_data% iorbtolr   ))*kind( pawproj_data% iorbtolr  )
       deallocate(pawproj_data%  iorbtolr  ,stat=i_stat)
       call memocc(i_stat,i_all,'iorbtolr',subname)

       call deallocate_proj_descr(pawproj_data%paw_nlpspd,subname)
!!$       i_all=-product(shape(pawproj_data%paw_nlpspd%nboxp_c))*kind(pawproj_data%paw_nlpspd%nboxp_c)
!!$       deallocate(pawproj_data%paw_nlpspd%nboxp_c,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nboxp_c',subname)
!!$       i_all=-product(shape(pawproj_data%paw_nlpspd%nboxp_f))*kind(pawproj_data%paw_nlpspd%nboxp_f)
!!$       deallocate(pawproj_data%paw_nlpspd%nboxp_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nboxp_f',subname)
!!$       i_all=-product(shape(pawproj_data%paw_nlpspd%keyg_p))*kind(pawproj_data%paw_nlpspd%keyg_p)
!!$       deallocate(pawproj_data%paw_nlpspd%keyg_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'keyg_p',subname)
!!$       i_all=-product(shape(pawproj_data%paw_nlpspd%keyv_p))*kind(pawproj_data%paw_nlpspd%keyv_p)
!!$       deallocate(pawproj_data%paw_nlpspd%keyv_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'keyv_p',subname)
!!$       i_all=-product(shape(pawproj_data%paw_nlpspd%nvctr_p))*kind(pawproj_data%paw_nlpspd%nvctr_p)
!!$       deallocate(pawproj_data%paw_nlpspd%nvctr_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nvctr_p',subname)
!!$       i_all=-product(shape(pawproj_data%paw_nlpspd%nseg_p))*kind(pawproj_data%paw_nlpspd%nseg_p)
!!$       deallocate(pawproj_data%paw_nlpspd%nseg_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nseg_p',subname)

       if(pawproj_data%DistProjApply) then
          call deallocate_gwf_c(pawproj_data%G,subname)
       endif
       nullify(pawproj_data%paw_proj)
    end if
  END SUBROUTINE deallocate_pawproj_data


  !> deallocate_pcproj_data
  subroutine deallocate_pcproj_data(pcproj_data,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(pcproj_data_type), intent(inout) :: pcproj_data
    !local variables
    integer :: i_all,i_stat
    if(associated(pcproj_data%pc_proj)) then

       i_all=-product(shape(  pcproj_data% pc_proj ))*kind( pcproj_data% pc_proj  )
       deallocate(pcproj_data%  pc_proj  ,stat=i_stat)
       call memocc(i_stat,i_all,'pc_proj',subname)
       
       i_all=-product(shape( pcproj_data%ilr_to_mproj   ))*kind(pcproj_data% ilr_to_mproj   )
       deallocate( pcproj_data% ilr_to_mproj  ,stat=i_stat)
       call memocc(i_stat,i_all,'ilr_to_mproj',subname)
       
       i_all=-product(shape( pcproj_data% iproj_to_ene  ))*kind( pcproj_data% iproj_to_ene  )
       deallocate(  pcproj_data% iproj_to_ene ,stat=i_stat)
       call memocc(i_stat,i_all,'iproj_to_ene',subname)

       i_all=-product(shape( pcproj_data% iproj_to_factor  ))*kind( pcproj_data% iproj_to_factor  )
       deallocate(  pcproj_data% iproj_to_factor ,stat=i_stat)
       call memocc(i_stat,i_all,'iproj_to_factor',subname)
       
       i_all=-product(shape( pcproj_data% iproj_to_l  ))*kind( pcproj_data% iproj_to_l  )
       deallocate(pcproj_data%  iproj_to_l  ,stat=i_stat)
       call memocc(i_stat,i_all,'iproj_to_l',subname)
       
       i_all=-product(shape( pcproj_data% iorbtolr   ))*kind( pcproj_data% iorbtolr  )
       deallocate(pcproj_data%  iorbtolr  ,stat=i_stat)
       call memocc(i_stat,i_all,'iorbtolr',subname)
       
       i_all=-product(shape( pcproj_data% gaenes   ))*kind( pcproj_data% gaenes  )
       deallocate(pcproj_data%  gaenes  ,stat=i_stat)
       call memocc(i_stat,i_all,'gaenes',subname)
       

       call deallocate_proj_descr(pcproj_data%pc_nlpspd,subname)

!!$       i_all=-product(shape(pcproj_data%pc_nlpspd%nboxp_c))*kind(pcproj_data%pc_nlpspd%nboxp_c)
!!$       deallocate(pcproj_data%pc_nlpspd%nboxp_c,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nboxp_c',subname)
!!$       i_all=-product(shape(pcproj_data%pc_nlpspd%nboxp_f))*kind(pcproj_data%pc_nlpspd%nboxp_f)
!!$       deallocate(pcproj_data%pc_nlpspd%nboxp_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nboxp_f',subname)
!!$       i_all=-product(shape(pcproj_data%pc_nlpspd%keyg_p))*kind(pcproj_data%pc_nlpspd%keyg_p)
!!$       deallocate(pcproj_data%pc_nlpspd%keyg_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'keyg_p',subname)
!!$       i_all=-product(shape(pcproj_data%pc_nlpspd%keyv_p))*kind(pcproj_data%pc_nlpspd%keyv_p)
!!$       deallocate(pcproj_data%pc_nlpspd%keyv_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'keyv_p',subname)
!!$       i_all=-product(shape(pcproj_data%pc_nlpspd%nvctr_p))*kind(pcproj_data%pc_nlpspd%nvctr_p)
!!$       deallocate(pcproj_data%pc_nlpspd%nvctr_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nvctr_p',subname)
!!$       i_all=-product(shape(pcproj_data%pc_nlpspd%nseg_p))*kind(pcproj_data%pc_nlpspd%nseg_p)
!!$       deallocate(pcproj_data%pc_nlpspd%nseg_p,stat=i_stat)
!!$       call memocc(i_stat,i_all,'nseg_p',subname)


       if(pcproj_data%DistProjApply) then
          call deallocate_gwf(pcproj_data%G,subname)
       endif


    end if
  END SUBROUTINE deallocate_pcproj_data

end module module_types
