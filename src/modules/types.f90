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
  integer, dimension(9), parameter :: input_psi_values = &
       & (/ INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_CP2K, &
       & INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, INPUT_PSI_DISK_WVL, &
       & INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS /)

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
  end type orthon_data

  type, public :: SIC_data
     character(len=4) :: approach !< approach for the Self-Interaction-Correction (PZ, NK)
     integer :: ixc !< base for the SIC correction
     real(gp) :: alpha !<downscaling coefficient
     real(gp) :: fref !< reference occupation (for alphaNK case)
  end type SIC_data

!> Structure of the variables read by input.* files (*.dft, *.geopt...)
  type, public :: input_variables
     !strings of the input files
     character(len=100) :: file_dft,file_geopt,file_kpt,file_perf,file_tddft, &
          & file_mix,file_sic, file_output
     !miscellaneous variables
     logical :: gaussian_help
     integer :: ixc,ncharge,itermax,nrepmax,ncong,idsx,ncongt,inputPsiId,nspin,mpol,itrpmax
     integer :: norbv,nvirt,nplot,iscf,norbsempty,norbsuempty,norbsdempty
     integer :: output_grid, dispersion,last_run,output_wf_format,output_grid_format
     real(gp) :: frac_fluct,gnrm_sw,alphamix,Tel,alphadiis
     real(gp) :: hx,hy,hz,crmult,frmult,gnrm_cv,rbuf,rpnrm_cv,gnrm_startmix
     integer :: verbosity
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
     ! tddft variables from *.tddft
     character(len=10) :: tddft_approach
     !variables for SIC
     type(SIC_data) :: SIC !<parameters for the SIC methods

     !> variable for material acceleration
     !! values 0: traditional CPU calculation
     !!        1: CUDA acceleration with CUBLAS
     !!        2: OpenCL acceleration (with CUBLAS one day)
     integer :: iacceleration

     ! Performance variables from input.perf
     logical :: debug      !< Debug option (used by memocc)
     integer :: ncache_fft !< Cache size for FFT
     real(gp) :: projrad   !<coarse radius of the projectors in units of the maxrad

     !orthogonalisation data
     type(orthon_data) :: orthpar

     !> parallelisation scheme of the exact exchange operator
     !!   BC (Blocking Collective)
     !!   OP2P (Overlap Point-to-Point)
     character(len=4) :: exctxpar

     !> communication scheme for the density
     !!  DBL traditional scheme with double precision
     !!  MIX mixed single-double precision scheme (requires rho_descriptors)
     character(len=3) :: rho_commun
  end type input_variables

  type, public :: energy_terms
     real(gp) :: eh,exc,vxc,eion,edisp,ekin,epot,eproj,eexctX
     real(gp) :: ebs,eKS,trH
  end type energy_terms

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

!> Used to split between points to be treated in simple or in double precision
  type, public :: rho_descriptors
     character(len=1) :: geocode
     integer :: icomm !< method for communicating the density
     integer :: n_csegs,n_fsegs,dp_size,sp_size
     integer, dimension(:,:), pointer :: spkey,dpkey
     integer, dimension(:), pointer :: cseg_b,fseg_b
  end type rho_descriptors

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
     real(gp), dimension(:,:), pointer :: aocc
     real(gp), dimension(:,:,:), pointer :: psppar         !< pseudopotential parameters (HGH SR section)
     logical :: donlcc                                     !< activate non-linear core correction treatment
     integer, dimension(:), pointer :: nlcc_ngv,nlcc_ngc   !<number of valence and core gaussians describing NLCC 
     real(gp), dimension(:,:), pointer :: nlccpar          !< parameters for the non-linear core correction, if present
     integer :: symObj                                     !< The symmetry object from ABINIT
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
     real(gp) :: efermi,HLgap
     integer, dimension(:), pointer :: norb_par,iokpt,ikptproc!,ikptsp
     real(wp), dimension(:), pointer :: eval
     real(gp), dimension(:), pointer :: occup,spinsgn,kwgts
     real(gp), dimension(:,:), pointer :: kpts
  end type orbitals_data


!> Contains the information needed for describing completely a
!! wavefunction localisation region
  type, public :: locreg_descriptors
     character(len=1) :: geocode
     logical :: hybrid_on             !<interesting for global, periodic, localisation regions
     integer :: ns1,ns2,ns3           !<starting points of the localisation region in global coordinates
     type(grid_dimensions) :: d
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
  end type locreg_descriptors


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
     real(gp) :: ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC
     type(atoms_data), pointer :: at
     type(orbitals_data) :: orbs
     type(communications_arrays) :: comms
     type(nonlocal_psp_descriptors), pointer :: nlpspd
     type(locreg_descriptors), pointer :: lr 
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
  end type lanczos_args


!> Contains the arguments needed for the diis procedure
  type, public :: diis_objects
     logical :: switchSD
     integer :: idiistol,mids,ids,idsx
     real(gp) :: energy_min,energy_old,energy,alpha,alpha_max
     real(tp), dimension(:), pointer :: psidst
     real(tp), dimension(:), pointer :: hpsidst
     real(tp), dimension(:,:,:,:,:,:), pointer :: ads
  end type diis_objects


contains


!> Allocate diis objects
  subroutine allocate_diis_objects(idsx,alphadiis,npsidim,nkptsp,nspinor,norbd,diis,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    integer, intent(in) :: idsx,npsidim,nkptsp,nspinor,norbd
    real(gp), intent(in) :: alphadiis
    type(diis_objects), intent(inout) :: diis
    !local variables
    integer :: i_stat,ncplx,ngroup

    !calculate the number of complex components
    if (nspinor > 1) then
       ncplx=2
    else
       ncplx=1
    end if

    !always better to allow real combination of the wavefunctions
    ncplx=1

    !add the possibility of more than one diis group
    ngroup=1

    allocate(diis%psidst(npsidim*idsx+ndebug),stat=i_stat)
    call memocc(i_stat,diis%psidst,'psidst',subname)
    allocate(diis%hpsidst(npsidim*idsx+ndebug),stat=i_stat)
    call memocc(i_stat,diis%hpsidst,'hpsidst',subname)
    allocate(diis%ads(ncplx,idsx+1,idsx+1,ngroup,nkptsp,1+ndebug),stat=i_stat)
    call memocc(i_stat,diis%ads,'ads',subname)
    call to_zero(nkptsp*ncplx*ngroup*(idsx+1)**2,diis%ads(1,1,1,1,1,1))

    !initialize scalar variables
    !diis initialisation variables
    diis%alpha=alphadiis
    diis%alpha_max=alphadiis
    diis%energy=1.d10
    !minimum value of the energy during the minimisation procedure
    diis%energy_min=1.d10
    !previous value already fulfilled
    diis%energy_old=diis%energy
    !local variable for the diis history
    diis%idsx=idsx
    !logical control variable for switch DIIS-SD
    diis%switchSD=.false.
    

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
    call memocc(i_stat,i_all,'orbs%ikptproc',subname)

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
    end if

    if (associated(rst%gaucoeffs)) then
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

    call deallocate_wfd(lr%wfd,subname)

    call deallocate_bounds(lr%geocode,lr%hybrid_on,lr%bounds,subname)

  END SUBROUTINE deallocate_lr

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
