!> @file
!!  Define the fortran types
!! @author
!!    Copyright (C) 2008-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!>  Modules which contains the Fortran data structures
!!  and the routines of allocations and de-allocations
module module_types

  use module_base, only : gp,wp,dp,tp
  implicit none

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
  integer, parameter :: INPUT_PSI_LINEAR       = 100
  integer, dimension(10), parameter :: input_psi_values = &
       & (/ INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_CP2K, &
       & INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, INPUT_PSI_DISK_WVL, &
       & INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS, &
       & INPUT_PSI_LINEAR /)

  !> Output wf parameters.
  integer, parameter :: WF_FORMAT_NONE   = 0
  integer, parameter :: WF_FORMAT_PLAIN  = 1
  integer, parameter :: WF_FORMAT_BINARY = 2
  integer, parameter :: WF_FORMAT_ETSF   = 3
  integer, parameter :: WF_N_FORMAT      = 4
  character(len = 12), dimension(0:WF_N_FORMAT-1), parameter :: wf_format_names = &
       & (/ "none        ", "plain text  ", "Fortran bin.", "ETSF        " /)

  !> Output grid parameters.
  integer, parameter :: OUTPUT_GRID_NONE    = 0
  integer, parameter :: OUTPUT_GRID_DENSITY = 1
  integer, parameter :: OUTPUT_GRID_DENSPOT = 2
  character(len = 12), dimension(0:2), parameter :: output_grid_names = &
       & (/ "none        ", "density     ", "dens. + pot." /)
  integer, parameter :: OUTPUT_GRID_FORMAT_TEXT = 0
  integer, parameter :: OUTPUT_GRID_FORMAT_ETSF = 1
  integer, parameter :: OUTPUT_GRID_FORMAT_CUBE = 2
  character(len = 4), dimension(0:2), parameter :: output_grid_format_names = &
       & (/ "text", "ETSF", "cube" /)

  !> Occupation parameters.
  integer, parameter :: SMEARING_DIST_ERF   = 1
  integer, parameter :: SMEARING_DIST_FERMI = 2
  character(len = 11), dimension(2), parameter :: smearing_names = &
       & (/ "Error func.", "Fermi      " /)
  ! To be moved as an input parameter later
  integer, parameter :: occopt = SMEARING_DIST_ERF


!> Structure of the variables read by input.* files (*.dft, *.geopt...)
  type, public :: input_variables
     !strings of the input files
     character(len=100) :: file_dft,file_geopt,file_kpt,file_perf,file_tddft,file_mix
     !miscellaneous variables
     logical :: output_wf,calc_tail,gaussian_help,read_ref_den,correct_offset
     integer :: ixc,ncharge,itermax,nrepmax,ncong,idsx,ncongt,inputPsiId,nspin,mpol,itrpmax
     integer :: norbv,nvirt,nplot,iscf,norbsempty,norbsuempty,norbsdempty
     integer :: output_grid, dispersion,last_run,output_wf_format,output_grid_format
     real(gp) :: frac_fluct,gnrm_sw,alphamix,Tel,alphadiis
     real(gp) :: hx,hy,hz,crmult,frmult,gnrm_cv,rbuf,rpnrm_cv,gnrm_startmix
     integer :: nvacancy,verbosity
     real(gp) :: elecfield
     logical :: disableSym

     ! For absorption calculations
     integer :: iabscalc_type   ! 0 non calc, 1 cheb ,  2 lanc
     integer :: iat_absorber, L_absorber
     real(gp), pointer:: Gabs_coeffs(:)
     logical ::  c_absorbtion , abscalc_alterpot, abscalc_eqdiff 
     integer ::  potshortcut
     integer ::  nsteps
     character(len=100) :: extraOrbital

     ! Frequencies calculations (finite difference)
     real(gp) :: freq_alpha
     integer :: freq_order
     integer :: freq_method

     ! kpoints related input variables
     integer :: nkpt, nkptv,ngroups_kptv
     integer, dimension(:), pointer :: nkptsv_group
     real(gp), pointer :: kpt(:,:), wkpt(:), kptv(:,:)
     character(len=100) :: band_structure_filename

     ! Geometry variables from *.geopt
     character(len=10) :: geopt_approach
     integer :: ncount_cluster_x, history
     real(gp) :: betax,forcemax,randdis
     integer :: optcell, ionmov, nnos
     real(gp) :: dtion, mditemp, mdftemp, noseinert, friction, mdwall
     real(gp) :: bmass, vmass, strprecon, strfact
     real(gp) :: strtarget(6)
     real(gp), pointer :: qmass(:)
     real(gp) :: dtinit,dtmax !for FIRE
     ! tddft vaiables from *.tddft
     character(len=10) :: tddft_approach

     !> variable for material acceleration
     !! values 0: traditional CPU calculation
     !!        1: CUDA acceleration with CUBLAS
     !!        2: OpenCL acceleration (with CUBLAS one day)
     integer :: iacceleration

     ! Performance variables from input.perf
     logical :: debug      !< Debug option (used by memocc)
     integer :: ncache_fft !< Cache size for FFT
     real(gp) :: projrad   !<coarse radius of the projectors in units of the maxrad

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

     !> parallelisation scheme of the exact exchange operator
     !!   BC (Blocking Collective)
     !!   OP2P (Overlap Point-to-Point)
     character(len=4) :: exctxpar
  end type input_variables


!>  Bounds for coarse and fine grids for kinetic operations
!!  Useful only for isolated systems AND in CPU
  type, public :: kinetic_bounds
     integer, dimension(:,:,:), pointer :: ibyz_c,ibxz_c,ibxy_c
     integer, dimension(:,:,:), pointer :: ibyz_f,ibxz_f,ibxy_f
  end type kinetic_bounds


!>  Bounds to compress the wavefunctions
!!  Useful only for isolated systems AND in CPU
  type, public :: shrink_bounds
     integer, dimension(:,:,:), pointer :: ibzzx_c,ibyyzz_c
     integer, dimension(:,:,:), pointer :: ibxy_ff,ibzzx_f,ibyyzz_f
  end type shrink_bounds


!>  Bounds to uncompress the wavefunctions
!!  Useful only for isolated systems AND in CPU
  type, public :: grow_bounds
     integer, dimension(:,:,:), pointer :: ibzxx_c,ibxxyy_c
     integer, dimension(:,:,:), pointer :: ibyz_ff,ibzxx_f,ibxxyy_f
  end type grow_bounds


!>  Bounds for convolutions operations
!!  Useful only for isolated systems AND in CPU
  type, public :: convolutions_bounds
     type(kinetic_bounds) :: kb
     type(shrink_bounds) :: sb
     type(grow_bounds) :: gb
     integer, dimension(:,:,:), pointer :: ibyyzz_r !< real space border
  end type convolutions_bounds


!>  Used for lookup table for compressed wavefunctions
  type, public :: wavefunctions_descriptors
     integer :: nvctr_c,nvctr_f,nseg_c,nseg_f
     integer, dimension(:,:), pointer :: keyg
     integer, dimension(:), pointer :: keyv
  end type wavefunctions_descriptors


!>  Non local pseudopotential descriptors
  type, public :: nonlocal_psp_descriptors
     integer :: nproj,nprojel                     !< Number of projectors and number of elements
     !> Projector segments on real space grid
     integer, dimension(:), pointer :: nvctr_p,nseg_p,keyv_p
     integer, dimension(:,:), pointer :: keyg_p 
     !> Parameters for the boxes containing the projectors
     integer, dimension(:,:,:), pointer :: nboxp_c,nboxp_f
  end type nonlocal_psp_descriptors


!>  Atomic data (name, polarisation, ...)
  type, public :: atoms_data
     character(len=1) :: geocode
     character(len=5) :: format
     character(len=20) :: units
     integer :: nat                                !< nat          Number of atoms
     integer :: ntypes                             !< ntypes       Number of type of atoms
     integer :: natsc
     character(len=20), dimension(:), pointer :: atomnames
     real(gp) :: alat1,alat2,alat3
     integer, dimension(:), pointer :: iatype      !< iatype(nat)  Type of the atoms
     integer, dimension(:), pointer :: iasctype,natpol,nelpsp,npspcode,nzatom
     integer, dimension(:), pointer :: ifrztyp     !< ifrztyp(nat) Frozen atoms
     real(gp), dimension(:), pointer :: amu        !< amu(ntypes)  Atomic Mass Unit for each type of atoms
     real(gp), dimension(:,:), pointer :: aocc
     real(gp), dimension(:,:,:), pointer :: psppar
     integer :: symObj                             !< The symmetry object from ABINIT
     integer :: iat_absorber 
  end type atoms_data


!>  Grid dimensions in old different wavelet basis
  type, public :: grid_dimensions
     integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i
  end type grid_dimensions


!>  Structures of basis of gaussian functions
  type, public :: gaussian_basis
     integer :: nat,ncoeff,nshltot,nexpo
     integer, dimension(:), pointer :: nshell,ndoc,nam
     real(gp), dimension(:), pointer :: xp,psiat
     real(gp), dimension(:,:), pointer :: rxyz
  end type gaussian_basis


!> All the parameters which are important for describing the orbitals
!! Add also the objects related to k-points sampling, after symmetries applications
  type, public :: orbitals_data
     integer :: norb,norbp,norbu,norbd,nspin,nspinor,isorb,npsidim,nkpts,nkptsp,iskpts
     real(gp) :: efermi
     integer, dimension(:), pointer :: norb_par,iokpt,ikptproc,inwhichlocreg, inWhichLocregP !,ikptsp
     integer,dimension(:),pointer:: onWhichMPI, isorb_par
     real(wp), dimension(:), pointer :: eval
     real(gp), dimension(:), pointer :: occup,spinsgn,kwgts
     real(gp), dimension(:,:), pointer :: kpts
  end type orbitals_data


!> Contains the information needed for describing completely a
!! wavefunction localisation region
  type, public :: locreg_descriptors
     character(len=1) :: geocode
     logical :: hybrid_on               !< interesting for global, periodic, localisation regions
     integer :: ns1,ns2,ns3             !< starting point of the localisation region in global coordinates
     integer :: nsi1,nsi2,nsi3          !< starting point of locreg for interpolating grid
     integer :: Localnorb               !< number of orbitals contained in locreg
     integer,dimension(3) :: outofzone  !< vector of points outside of the zone outside Glr for periodic systems
     integer,dimension(:),pointer :: projflg    !< atoms contributing nlpsp projectors to locreg
     type(grid_dimensions) :: d
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
  end type locreg_descriptors


!> Contains the information needed for communicating the wavefunctions
!! between processors for the transposition
  type, public :: communications_arrays
     integer, dimension(:), pointer :: ncntd,ncntt,ndspld,ndsplt
     integer, dimension(:,:), pointer :: nvctr_par
  !integer,dimension(:,:,:,:),pointer:: nvctr_parLIN
  !integer, dimension(:), pointer :: ncntdLIN,ncnttLIN,ndspldLIN,ndspltLIN

  end type communications_arrays


!> Contains the pointers to be handled to control GPU information
!! Given that they are pointers on GPU address, they are C pointers
!! which take 8 bytes
!! So they are declared as kind=8 variables either if the GPU works in simple precision
!! Also other information concerning the GPU runs can be stored in this structure
  type, public :: GPU_pointers
     logical :: useDynamic,full_locham
     integer :: id_proc
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
  end type GPU_pointers


!>  Used to restart a new DFT calculation or to save information 
!!  for post-treatment
  type, public :: restart_objects
     integer :: n1,n2,n3
     real(gp) :: hx_old,hy_old,hz_old
     real(wp), dimension(:), pointer :: psi 
     real(wp), dimension(:,:), pointer :: gaucoeffs
     real(gp), dimension(:,:), pointer :: rxyz_old,rxyz_new
     type(locreg_descriptors) :: Glr
     type(gaussian_basis) :: gbd
     type(orbitals_data) :: orbs
     type(GPU_pointers) :: GPU
  end type restart_objects


!> Contains the work arrays needed for expressing wavefunction in real space
!!  with all the BC
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
     integer :: iproc,nproc,ndimpot,nspin
     real(gp) :: hx,hy,hz
     real(gp) :: ekin_sum,epot_sum,eexctX,eproj_sum
     type(atoms_data), pointer :: at
     type(orbitals_data) :: orbs
     type(communications_arrays) :: comms
     type(nonlocal_psp_descriptors), pointer :: nlpspd
     type(locreg_descriptors), pointer :: lr 
     type(gaussian_basis), pointer :: Gabsorber    
     integer, dimension(:,:), pointer :: ngatherarr 
     real(gp), dimension(:,:),  pointer :: rxyz,radii_cf
     real(wp), dimension(:), pointer :: proj
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), pointer :: psi
     real(wp), dimension(:), pointer :: potential
     real(wp), dimension(:), pointer :: Gabs_coeffs
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp) :: hpsi
     type(GPU_pointers), pointer :: GPU
  end type lanczos_args


!> Contains the dimensions of some arrays.
  type,public:: arraySizes
      integer:: size_rhopot
      integer,dimension(4):: size_potxc
      integer:: size_rhocore
      integer:: size_pot_ion
      integer,dimension(3):: size_phnons
      integer,dimension(3):: size_irrzon
      integer:: size_pkernel
      integer:: size_pkernelseq
  end type

!> Contains the parameters needed for the point to point communications
!! for sumrho in the linear scaling version.
  type,public:: p2pCommsSumrho
    integer,dimension(:),pointer:: noverlaps, overlaps, istarr, istrarr
    real(8),dimension(:),pointer:: sendBuf, recvBuf
    integer,dimension(:,:,:),pointer:: comarr
    integer:: nsendBuf, nrecvBuf
    logical,dimension(:,:),pointer:: communComplete, computComplete
  end type p2pCommsSumrho

!> Contains the parameters neeed for the point to point communications
!! for gathering the potential (for the application of the Hamiltonian)
   type,public:: p2pCommsGatherPot
       integer,dimension(:),pointer:: noverlaps, overlaps
       integer,dimension(:,:),pointer:: ise3 ! starting / ending index of recvBuf in z dimension after communication (glocal coordinates)
       integer,dimension(:,:,:),pointer:: comarr
       real(8),dimension(:),pointer:: recvBuf
       integer:: nrecvBuf
       logical,dimension(:,:),pointer:: communComplete
   end type p2pCommsGatherPot

!> Contains the parameter needed for the point to point communication for
!! the orthonormlization.
   type,public:: p2pCommsOrthonormality
       integer:: nsendBuf, nrecvBuf
       integer,dimension(:),pointer:: noverlaps
       integer,dimension(:,:),pointer:: overlaps
       integer,dimension(:,:,:),pointer:: comarr
       real(8),dimension(:),pointer:: sendBuf, recvBuf
       logical,dimension(:,:),pointer:: communComplete
   end type p2pCommsOrthonormality


!> Contains the parameters for the communications of the derivative orbitals
!! to mathc their partition.
  type,public:: p2pCommsRepartition
      integer,dimension(:,:,:),pointer:: comarr
       logical,dimension(:,:),pointer:: communComplete
  end type p2pCommsRepartition

!! Contains the parameters for calculating the overlap matrix for the orthonormalization etc...
  type,public:: overlapParameters
      integer:: ndim_lphiovrlp
      integer,dimension(:),pointer:: noverlaps, indexExpand, indexExtract
      integer,dimension(:,:),pointer:: overlaps
      integer,dimension(:,:),pointer:: indexInRecvBuf
      integer,dimension(:,:),pointer:: indexInSendBuf
      type(locreg_descriptors),dimension(:,:),pointer:: olr
  end type overlapParameters


  type,public:: matrixLocalizationRegion
      integer:: norbinlr
      integer,dimension(:),pointer:: indexInGlobal
  end type matrixLocalizationRegion


  type,public:: p2pCommsOrthonormalityMatrix
      integer:: nrecvBuf, nsendBuf
      integer,dimension(:),pointer:: noverlap, noverlapProc
      integer,dimension(:,:),pointer:: overlaps, indexInRecvBuf, overlapsProc
      integer,dimension(:,:,:),pointer:: comarr, olrForExpansion
      real(8),dimension(:),pointer:: recvBuf, sendBuf
      logical,dimension(:,:),pointer:: communComplete
      type(matrixLocalizationRegion),dimension(:,:),pointer:: olr
  end type p2pCommsOrthonormalityMatrix


  type,public:: matrixMinimization
    type(matrixLocalizationRegion),dimension(:),pointer:: mlr
    integer:: norbmax ! maximal matrix size handled by a given process
    integer:: nlrp ! number of localization regions handled by a given process
    integer,dimension(:),pointer:: inWhichLocregExtracted
    integer,dimension(:),pointer:: inWhichLocregOnMPI
    integer,dimension(:),pointer:: indexInLocreg
  end type matrixMinimization


!!!> Contains all the descriptors necessary for splitting the calculation in different locregs 
  type,public:: linear_zone_descriptors
    integer :: nlr                                              !> Number of localization regions 
    integer :: Lpsidimtot                                       !> Total dimension of the wavefunctions in the locregs
    type(orbitals_data) :: orbs                                !> Global orbitals descriptors
    type(orbitals_data),dimension(:),pointer:: Lorbs            !> Orbitals descriptors for each locreg
    type(communications_arrays) :: comms                        !> Global communication descriptors
    type(locreg_descriptors) :: Glr                             !> Global region descriptors
    type(nonlocal_psp_descriptors) :: Gnlpspd                   !> Global nonlocal pseudopotential descriptors
    type(locreg_descriptors),dimension(:),pointer :: Llr                !> Local region descriptors (dimension = nlr)
    type(nonlocal_psp_descriptors),dimension(:),pointer :: Lnlpspd      !> Nonlocal pseudopotential descriptors for locreg (dimension = nlr)
  end type


!> Contains all parameters for the basis with which we calculate the properties
!! like energy and forces. Since we may also use the derivative of the trace
!! minimizing orbitals, this basis may be larger than only the trace minimizing
!! orbitals. In case we don't use the derivatives, these parameters are identical
!! from those in lin%orbs etc.
type,public:: largeBasis
    type(communications_arrays):: comms, gcomms
    type(orbitals_data):: orbs, gorbs
    !type(linear_zone_descriptors):: lzd
    type(p2pCommsRepartition):: comrp
    type(p2pCommsOrthonormality):: comon
    type(overlapParameters):: op
    type(p2pCommsGatherPot):: comgp
end type largeBasis

  type,public:: matrixDescriptors
      integer:: nvctr, nseg, nvctrmatmul, nsegmatmul
      integer,dimension(:),pointer:: keyv, keyvmatmul, nsegline
      integer,dimension(:,:),pointer:: keyg, keygmatmul
      integer,dimension(:,:,:),pointer:: keygline
  end type matrixDescriptors


  !> Contains the parameters for the parallel input guess for the O(N) version.
  type,public:: inguessParameters
    integer:: nproc, norb, norbtot, norbtotPad, sizeWork, nvctrp, isorb
    integer,dimension(:),pointer:: norb_par, onWhichMPI, isorb_par, nvctrp_nz, sendcounts, senddispls, recvcounts, recvdispls
    !!type(matrixLocalizationRegion),dimension(:),pointer:: mlr
  end type inguessParameters

  type,public:: localizedDIISParameters
    integer:: is, isx, mis
    real(8),dimension(:),pointer:: phiHist, hphiHist
    real(8),dimension(:,:,:),pointer:: mat
    real(8):: trmin, trold
    logical:: switchSD
  end type localizedDIISParameters

  type,public:: mixrhopotDIISParameters
    integer:: is, isx, mis
    real(8),dimension(:),pointer:: rhopotHist, rhopotresHist
    real(8),dimension(:,:),pointer:: mat
  end type mixrhopotDIISParameters

  type,public:: linearInputGuess
      type(linear_zone_descriptors):: lzdig, lzdGauss
      type(orbitals_data):: orbsig, orbsGauss
      type(p2pCommsOrthonormality):: comon
      type(overlapParameters):: op
      type(p2pCommsGatherPot):: comgp
      type(matrixDescriptors):: mad
  end type linearInputGuess



!> Contains all parameters related to the linear scaling version.
  type,public:: linearParameters
    integer:: DIISHistMin, DIISHistMax, nItBasisFirst, nItBasis, nItPrecond, nItCoeff, nItSCC, confPotOrder, norbsPerProcIG
    integer:: nItInguess, nlr, nLocregOverlap, nItOrtho, mixHist, methTransformOverlap, blocksize_pdgemm, blocksize_pdsyev
    integer:: correctionOrthoconstraint, nproc_pdsyev, nproc_pdgemm
    real(8):: convCrit, alphaSD, alphaDIIS, startDIIS, convCritCoeff, alphaMix, convCritMix, convCritOrtho
    real(8),dimension(:),pointer:: potentialPrefac, locrad, lphiRestart, lphiold, lhphiold
    real(8),dimension(:,:),pointer:: hamold
    type(orbitals_data):: orbs, gorbs
    type(communications_arrays):: comms, gcomms
    integer,dimension(:),pointer:: norbsPerType
    type(arraySizes):: as
    logical:: plotBasisFunctions, startWithSD, useDerivativeBasisFunctions
    character(len=4):: getCoeff, mixingMethod
    type(p2pCommsSumrho):: comsr
    type(p2pCommsGatherPot):: comgp
    type(largeBasis):: lb
    type(linear_zone_descriptors):: lzd
    type(p2pCommsOrthonormality):: comon
    type(overlapParameters):: op
    type(linearInputGuess):: lig
    type(matrixDescriptors):: mad
  end type linearParameters

!> Contains the arguments needed for the diis procedure
  type, public :: diis_objects
     logical :: switchSD
     integer :: idiistol,mids,ids,idsx
     real(gp) :: energy_min,energy_old,energy,alpha,alpha_max
     real(wp), dimension(:), pointer :: psidst
     real(tp), dimension(:), pointer :: hpsidst
     real(wp), dimension(:,:,:,:), pointer :: ads
  end type diis_objects


contains


!> Allocate diis objects
  subroutine allocate_diis_objects(idsx,npsidim,nkptsp,diis,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    integer, intent(in) :: idsx,npsidim,nkptsp
    type(diis_objects), intent(inout) :: diis
    !local variables
    integer :: i_stat
    allocate(diis%psidst(npsidim*idsx+ndebug),stat=i_stat)
    call memocc(i_stat,diis%psidst,'psidst',subname)
    allocate(diis%hpsidst(npsidim*idsx+ndebug),stat=i_stat)
    call memocc(i_stat,diis%hpsidst,'hpsidst',subname)
    allocate(diis%ads(idsx+1,idsx+1,nkptsp,3+ndebug),stat=i_stat)
    call memocc(i_stat,diis%ads,'ads',subname)
    call razero(nkptsp*3*(idsx+1)**2,diis%ads)
  END SUBROUTINE allocate_diis_objects


!> De-Allocate diis objects
  subroutine deallocate_diis_objects(diis,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(diis_objects), intent(inout) :: diis
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(diis%psidst))*kind(diis%psidst)
    deallocate(diis%psidst,stat=i_stat)
    call memocc(i_stat,i_all,'psidst',subname)
    i_all=-product(shape(diis%hpsidst))*kind(diis%hpsidst)
    deallocate(diis%hpsidst,stat=i_stat)
    call memocc(i_stat,i_all,'hpsidst',subname)
    i_all=-product(shape(diis%ads))*kind(diis%ads)
    deallocate(diis%ads,stat=i_stat)
    call memocc(i_stat,i_all,'ads',subname)

  END SUBROUTINE deallocate_diis_objects


!> Allocate communications_arrays
  subroutine allocate_comms(nproc,orbs,comms,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    integer, intent(in) :: nproc
    type(orbitals_data), intent(in) :: orbs
    type(communications_arrays), intent(out) :: comms
    !local variables
    integer :: i_stat

    allocate(comms%nvctr_par(0:nproc-1,orbs%nkptsp+ndebug),stat=i_stat)
    call memocc(i_stat,comms%nvctr_par,'nvctr_par',subname)
    allocate(comms%ncntd(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ncntd,'ncntd',subname)
    allocate(comms%ncntt(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ncntt,'ncntt',subname)
    allocate(comms%ndspld(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ndspld,'ndspld',subname)
    allocate(comms%ndsplt(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ndsplt,'ndsplt',subname)
  END SUBROUTINE allocate_comms


!> De-Allocate communications_arrays
  subroutine deallocate_comms(comms,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(communications_arrays), intent(out) :: comms
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
    i_all=-product(shape(orbs%isorb_par))*kind(orbs%isorb_par)
    deallocate(orbs%isorb_par,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%isorb_par',subname)
    i_all=-product(shape(orbs%onWhichMPI))*kind(orbs%onWhichMPI)
    deallocate(orbs%onWhichMPI,stat=i_stat)
    call memocc(i_stat,i_all,'orbs%onWhichMPI',subname)

    !contradictory: needed for component distribution and allocated for
    !               orbital distribution. Better to deal with scalars
    !i_all=-product(shape(orbs%ikptsp))*kind(orbs%ikptsp)
    !deallocate(orbs%ikptsp,stat=i_stat)
    !call memocc(i_stat,i_all,'orbs%ikptsp',subname)

END SUBROUTINE deallocate_orbs


!> Allocate and nullify restart objects
  subroutine init_restart_objects(iproc,iacceleration,atoms,rst,subname)
    use module_base
    implicit none
    !Arguments
    character(len=*), intent(in) :: subname
    integer, intent(in) :: iproc,iacceleration
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
    nullify(rst%psi)
    nullify(rst%orbs%eval)

    nullify(rst%gaucoeffs)

    nullify(rst%Glr%wfd%keyg)
    nullify(rst%Glr%wfd%keyv)

    nullify(rst%gbd%nshell)
    nullify(rst%gbd%ndoc)
    nullify(rst%gbd%nam)
    nullify(rst%gbd%xp)
    nullify(rst%gbd%psiat)
    nullify(rst%gbd%rxyz)

    !initialise the acceleration strategy if required
    call init_material_acceleration(iproc,iacceleration,rst%GPU)

  END SUBROUTINE init_restart_objects


!>  De-Allocate restart_objects
  subroutine free_restart_objects(rst,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(restart_objects) :: rst
    !local variables
    integer :: i_all,i_stat

    call deallocate_wfd(rst%Glr%wfd,subname)

    if (associated(rst%psi)) then
       i_all=-product(shape(rst%psi))*kind(rst%psi)
       deallocate(rst%psi,stat=i_stat)
       call memocc(i_stat,i_all,'psi',subname)
    end if

    if (associated(rst%orbs%eval)) then
       i_all=-product(shape(rst%orbs%eval))*kind(rst%orbs%eval)
       deallocate(rst%orbs%eval,stat=i_stat)
       call memocc(i_stat,i_all,'eval',subname)
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
    if (associated(rst%gbd%rxyz)) then
       nullify(rst%gbd%rxyz)
       call deallocate_gwf(rst%gbd,subname)

       i_all=-product(shape(rst%gaucoeffs))*kind(rst%gaucoeffs)
       deallocate(rst%gaucoeffs,stat=i_stat)
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

    allocate(wfd%keyg(2,wfd%nseg_c+wfd%nseg_f+ndebug),stat=i_stat)
    call memocc(i_stat,wfd%keyg,'keyg',subname)
    allocate(wfd%keyv(wfd%nseg_c+wfd%nseg_f+ndebug),stat=i_stat)
    call memocc(i_stat,wfd%keyv,'keyv',subname)
  END SUBROUTINE allocate_wfd


!> De-Allocate wavefunctions_descriptors
  subroutine deallocate_wfd(wfd,subname)
    use module_base
    implicit none
    type(wavefunctions_descriptors) :: wfd
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    if (associated(wfd%keyg)) then
       i_all=-product(shape(wfd%keyg))*kind(wfd%keyg)
       deallocate(wfd%keyg,stat=i_stat)
       call memocc(i_stat,i_all,'wfd%keyg',subname)
    end if
    if (associated(wfd%keyv)) then
       i_all=-product(shape(wfd%keyv))*kind(wfd%keyv)
       deallocate(wfd%keyv,stat=i_stat)
       call memocc(i_stat,i_all,'wfd%keyv',subname)
    end if
  END SUBROUTINE deallocate_wfd


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
    end if

    !the arrays which are needed only for free BC
    if (geocode == 'F') then
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
    end if

  END SUBROUTINE deallocate_bounds


  subroutine deallocate_lr(lr,subname)
    use module_base
    character(len=*), intent(in) :: subname
    type(locreg_descriptors) :: lr
    integer :: i_all,i_stat

    call deallocate_wfd(lr%wfd,subname)

    call deallocate_bounds(lr%geocode,lr%hybrid_on,lr%bounds,subname)
    if(associated(lr%projflg)) then
       nullify(lr%projflg)
!       i_all=-product(shape(lr%projflg)*kind(lr%projflg))
!       deallocate(lr%projflg,stat=i_stat)
!       call memocc(i_stat,i_all,'lr%projflg',subname)
    end if

  END SUBROUTINE deallocate_lr

  subroutine deallocate_Lzd(Lzd,subname)
    use module_base
    character(len=*), intent(in) :: subname
    type(linear_zone_descriptors) :: Lzd
    integer :: i_all,i_stat,ilr

!    call deallocate_comms(Lzd%comms,subname)

!    call deallocate_lr(Lzd%Glr,subname)

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
   nullify(Lzd%Glr%wfd%keyg)
   nullify(Lzd%Glr%wfd%keyv)
 
!Now destroy the Llr
    do ilr = 1, Lzd%nlr 
       call deallocate_lr(Lzd%Llr(ilr),subname)
       call deallocate_Lnlpspd(Lzd%Lnlpspd(ilr),subname)
    end do
     nullify(Lzd%Llr)
     nullify(Lzd%Lnlpspd)

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

  subroutine output_grid_help()
    integer :: i, j

    write(*, "(1x,A)") "Available values of output_grid are:"
    do i = 0, size(output_grid_format_names) - 1
       do j = 0, size(output_grid_names) - 1
          if (j == 0 .and. i == 0) then
             write(*, "(1x,A,I5,A,A,A)") " | ", i * 10 + j, &
                  & " - ", trim(output_grid_names(j)), "."
          else if (j /= 0) then
             write(*, "(1x,A,I5,A,A,A,A,A)") " | ", i * 10 + j, &
                  & " - ", trim(output_grid_names(j)), &
                  & " in ", trim(output_grid_format_names(i)), " format."
          end if
       end do
    end do
  end subroutine output_grid_help

  function output_grid_validate(id, fid)
    integer, intent(in) :: id, fid
    logical :: output_grid_validate

    output_grid_validate = (id >= 0 .and. id < size(output_grid_names)) .and. &
         & (fid >= 0 .and. fid < size(output_grid_format_names))
  end function output_grid_validate

end module module_types
