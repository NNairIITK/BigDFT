!!****m* BigDFT/module_types
!! FUNCTION
!!  Modules which contains the Fortran data structures
!!  and the routines of allocations and de-allocations
!! AUTHOR
!!    Luigi Genovese
!! COPYRIGHT
!!    Copyright (C) 2008 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!! 
module module_types

  use module_base, only : gp,wp,dp
  implicit none
!!***

!!****t* module_types/input_variables
!! DESCRIPTION
!!   Input variable structure
!!   Structure of the variables read by input.* files (*.dft, *.geopt...)
!! SOURCE
!!
  type, public :: input_variables
     logical :: output_wf,calc_tail,gaussian_help,read_ref_den,correct_offset
     integer :: ixc,ncharge,itermax,nrepmax,ncong,idsx,ncongt,inputPsiId,nspin,mpol,nvirt,nplot
     integer :: output_grid, dispersion
     real(gp) :: frac_fluct,gnrm_sw
     real(gp) :: hx,hy,hz,crmult,frmult,gnrm_cv,rbuf 
     integer :: nvacancy,verbosity
     real(gp) :: elecfield
     real(gp):: absorber_gnrm
     integer :: iabscalc_type   ! 0 non calc, 1 cheb ,  2 lanc
     integer :: iat_absorber, L_absorber
     real(gp), pointer:: Gabs_coeffs(:)
     logical ::  c_absorbtion , abscalc_alterpot, abscalc_eqdiff 
     integer ::  potshortcut
     integer ::  nsteps
     
     ! kpoints related input variables
     integer :: nkpt
     real(gp), pointer :: kpt(:,:), wkpt(:)
     ! Geometry variables from *.geopt
     character(len=10) :: geopt_approach
     integer :: ncount_cluster_x
     real(gp) :: betax,forcemax,randdis
     integer :: optcell, ionmov, nnos
     real(gp) :: dtion, mditemp, mdftemp, noseinert, friction, mdwall
     real(gp) :: bmass, vmass, strprecon, strfact
     real(gp) :: strtarget(6)
     real(gp), pointer :: qmass(:)
  end type input_variables
!!***

!!****t* convolution_bounds/kinetic_bounds
!! DESCRIPTION
!!   Bounds for coarse and fine grids for kinetic operations
!!   Useful only for isolated systems AND in CPU
!! SOURCE
!!
  type, public :: kinetic_bounds
     integer, dimension(:,:,:), pointer :: ibyz_c,ibxz_c,ibxy_c
     integer, dimension(:,:,:), pointer :: ibyz_f,ibxz_f,ibxy_f
  end type kinetic_bounds
!!***

!!****t* convolution_bounds/shrink_bounds
!! DESCRIPTION
!!   Bounds to compress the wavefunctions
!!   Useful only for isolated systems AND in CPU
!! SOURCE
!!
  type, public :: shrink_bounds
     integer, dimension(:,:,:), pointer :: ibzzx_c,ibyyzz_c
     integer, dimension(:,:,:), pointer :: ibxy_ff,ibzzx_f,ibyyzz_f
  end type shrink_bounds
!!***

!!****t* convolution_bounds/grow_bounds
!! DESCRIPTION
!!   Bounds to uncompress the wavefunctions
!!   Useful only for isolated systems AND in CPU
!! SOURCE
!!
  type, public :: grow_bounds
     integer, dimension(:,:,:), pointer :: ibzxx_c,ibxxyy_c
     integer, dimension(:,:,:), pointer :: ibyz_ff,ibzxx_f,ibxxyy_f
  end type grow_bounds
!!***

!!****t* module_types/convolutions_bounds
!! DESCRIPTION
!!   Bounds for convolutions operations
!!   Useful only for isolated systems AND in CPU
!! SOURCE
!!
  type, public :: convolutions_bounds
     type(kinetic_bounds) :: kb
     type(shrink_bounds) :: sb
     type(grow_bounds) :: gb
     integer, dimension(:,:,:), pointer :: ibyyzz_r ! real space border
  end type convolutions_bounds
!!***
  
!!****t* module_types/wavefunctions_descriptors
!! DESCRIPTION
!!   Used for lookup table for compressed wavefunctions
!! SOURCE
!!
  type, public :: wavefunctions_descriptors
     integer :: nvctr_c,nvctr_f,nseg_c,nseg_f
     integer, dimension(:,:), pointer :: keyg
     integer, dimension(:), pointer :: keyv
  end type wavefunctions_descriptors
!!***

!!****t* module_types/nonlocal_psp_descriptors
!! DESCRIPTION
!!   Non local pseudopotential descriptors
!! SOURCE
!!
  type, public :: nonlocal_psp_descriptors
     !number of projectors and number of elements
     integer :: nproj,nprojel
     ! projector segments on real space grid
     integer, dimension(:), pointer :: nvctr_p,nseg_p,keyv_p
     integer, dimension(:,:), pointer :: keyg_p 
     ! Parameters for the boxes containing the projectors
     integer, dimension(:,:,:), pointer :: nboxp_c,nboxp_f
  end type nonlocal_psp_descriptors
!!***

!!****t* module_types/atoms_data
!! DESCRIPTION
!!   Atomic data (name, polarisation, ...)
!! nat          Number of atoms
!! ntypes       Number of type of atoms
!! iatype(nat)  Type of the atoms
!! lfrztyp(nat) Frozen atoms
!! amu(ntypes)  Atomic Mass Unit for each type of atoms
!! SOURCE
!!
  type, public :: atoms_data
     character(len=1) :: geocode
     character(len=5) :: format
     character(len=20) :: units
     integer :: nat,ntypes,natsc
     character(len=20), dimension(:), pointer :: atomnames
     real(gp) :: alat1,alat2,alat3
     integer, dimension(:), pointer :: iatype,iasctype,natpol,nelpsp,npspcode,nzatom,ifrztyp
     real(gp), dimension(:,:,:), pointer :: psppar
     real(gp), dimension(:), pointer :: amu
     ! The symmetry object from ABINIT
     integer :: symObj
     ! AMmodif
     integer :: iat_absorber 
     ! AMmodif end
  end type atoms_data
!!***

!!****t* module_types/grid_dimensions
!! DESCRIPTION
!!   Grid dimensions in old different wavelet basis
!! SOURCE
!!
  type, public :: grid_dimensions
     integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i
  end type grid_dimensions
!!***

!!****t* module_types/gaussian_basis
!! DESCRIPTION
!!   Structures of basis of gaussian functions
!! SOURCE
!!
  type, public :: gaussian_basis
     integer :: nat,ncoeff,nshltot,nexpo
     integer, dimension(:), pointer :: nshell,ndoc,nam
     real(gp), dimension(:), pointer :: xp,psiat
     real(gp), dimension(:,:), pointer :: rxyz
  end type gaussian_basis
!!***

!!****t* module_types/orbitals_data
!! DESCRIPTION
!! All the parameters which are important for describing the orbitals
!! Add also the objects related to k-points sampling, after symmetries applications
!!
!! SOURCE
!!
  type, public :: orbitals_data
     integer :: norb,norbp,norbu,norbd,nspinor,isorb,npsidim,nkpts,nkptsp,iskpts
     integer, dimension(:), pointer :: norb_par,iokpt,ikptproc
     real(wp), dimension(:), pointer :: eval
     real(gp), dimension(:), pointer :: occup,spinsgn,kwgts
     real(gp), dimension(:,:), pointer :: kpts
  end type orbitals_data
!!***

!!****t* module_types/locreg_descriptors
!! DESCRIPTION
!! Contains the information needed for describing completely a
!! wavefunction localisation region
!! SOURCE
!!
  type, public :: locreg_descriptors
     character(len=1) :: geocode
     logical :: hybrid_on !interesting for global, periodic, localisation regions
     integer :: ns1,ns2,ns3 !starting points of the localisation region in global coordinates
     type(grid_dimensions) :: d
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
  end type locreg_descriptors
!!***

!!****t* module_types/restart_objects
!! DESCRIPTION
!!  Used to restart a new DFT calculation or to save information 
!!  for post-treatment
!! SOURCE
!!
  type, public :: restart_objects
     integer :: n1,n2,n3
     real(gp) :: hx_old,hy_old,hz_old
     real(wp), dimension(:), pointer :: psi 
     real(wp), dimension(:,:), pointer :: gaucoeffs
     real(gp), dimension(:,:), pointer :: rxyz_old
     type(locreg_descriptors) :: Glr
     type(gaussian_basis) :: gbd
     type(orbitals_data) :: orbs
  end type restart_objects
!!***

!!****t* module_types/communications_arrays
!! DESCRIPTION
!! Contains the information needed for communicating the wavefunctions
!! between processors for the transposition
!!
!! SOURCE
!!
  type, public :: communications_arrays
     integer, dimension(:), pointer :: ncntd,ncntt,ndspld,ndsplt
     integer, dimension(:,:), pointer :: nvctr_par
  end type communications_arrays
!!***

!!****t* module_types/GPU_pointers
!! DESCRIPTION
!! Contains the pointers to be handled to control GPU information
!! Given that they are pointers on GPU address, they are C pointers
!! which take 8 bytes
!! So they are declared as kind=8 variables either if the GPU works in simple precision
!! Also other information concerning the GPU runs can be stored in this structure
!!
!! SOURCE
!!
  type, public :: GPU_pointers
     logical :: useDynamic,full_locham
     real(kind=8) :: keys,work1,work2,work3,rhopot,r,d
     real(kind=8) :: pinned_in,pinned_out
     real(kind=8), dimension(:), pointer :: psi
  end type GPU_pointers
!!***


!!****t* module_types/workarr_sumrho
!! DESCRIPTION
!! Contains the work arrays needed for expressing wavefunction in real space
!!  with all the BC
!!
!! SOURCE
!!
  type, public :: workarr_sumrho
     integer :: nw1,nw2,nxc,nxf
     real(wp), dimension(:), pointer :: x_c,x_f,w1,w2
  end type workarr_sumrho
!!***


!!****t* module_types/workarr_locham
!! DESCRIPTION
!! Contains the work arrays needed for hamiltonian application with all the BC
!!
!! SOURCE
!!
  type, public :: workarr_locham
     integer :: nw1,nw2,nxc,nyc,nxf1,nxf2,nxf3,nxf,nyf
     real(wp), dimension(:), pointer :: w1,w2
     !for the periodic BC case, these arrays substitute 
     !psifscf,psifscfk,psig,ww respectively
     real(wp), dimension(:,:), pointer :: x_c,y_c,x_f1,x_f2,x_f3,x_f,y_f
  end type workarr_locham
!!***

!!****t* module_types/workarr_precond
!! DESCRIPTION
!! Contains the work arrays needed for th preconditioner with all the BC
!! Take different pointers depending on the boundary conditions
!!
!! SOURCE
!!
  type, public :: workarr_precond
     integer, dimension(:), pointer :: modul1,modul2,modul3
     real(wp), dimension(:), pointer :: psifscf,ww,x_f1,x_f2,x_f3,kern_k1,kern_k2,kern_k3
     real(wp), dimension(:,:), pointer :: af,bf,cf,ef
     real(wp), dimension(:,:,:), pointer :: xpsig_c,ypsig_c,x_c
     real(wp), dimension(:,:,:,:), pointer :: xpsig_f,ypsig_f,x_f,y_f
     real(wp), dimension(:,:,:,:,:), pointer :: z1,z3 ! work array for FFT
  end type workarr_precond
!!***

!!****t* module_types/lanczos_args
!! DESCRIPTION
!! Contains the arguments needed for the application of the hamiltonian
!!
!! SOURCE
!!
  type, public :: lanczos_args
     !arguments for the hamiltonian
     integer :: iproc,nproc,ndimpot,nspin
     real(gp) :: hx,hy,hz,cpmult,fpmult
     real(gp) :: ekin_sum,epot_sum,eproj_sum
     type(atoms_data), pointer :: at
     type(orbitals_data) :: orbs
     type(communications_arrays) :: comms
     type(nonlocal_psp_descriptors), pointer :: nlpspd
     type(locreg_descriptors), pointer :: lr 
     type(gaussian_basis), pointer :: Gabsorber    
     integer, dimension(:,:), pointer :: ngatherarr 
     real(gp), dimension(:,:),  pointer :: rxyz
     real(gp), dimension(:,:), pointer :: radii_cf  
     real(wp), dimension(:), pointer :: proj
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), pointer :: psi
     real(wp), dimension(:,:), pointer :: potential
     real(wp), dimension(:), pointer :: Gabs_coeffs
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp) :: hpsi
     type(GPU_pointers), pointer :: GPU
  end type lanczos_args
!!***


contains


!!****f* module_types/allocate_comms
!! FUNCTION
!!   Allocate communications_arrays
!! SOURCE
!!
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
  end subroutine allocate_comms
!!***

!!****f* module_types/deallocate_comms
!! FUNCTION
!!   De-Allocate communications_arrays
!! SOURCE
!!
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
  end subroutine deallocate_comms
!!***

!!****f* module_types/deallocate_orbs
!! FUNCTION
!!   De-Allocate orbitals data structure, except eval pointer
!!   which is not allocated in the orbitals_descriptor routine
!! SOURCE
!!
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


end subroutine deallocate_orbs

!!****f* module_types/init_restart_objects
!! FUNCTION
!!   Allocate and nullify restart objects
!! SOURCE
!!
  subroutine init_restart_objects(atoms,rst,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(atoms_data) :: atoms
    type(restart_objects) :: rst
    !local variables
    integer :: i_stat

    !allocate pointers
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

  end subroutine init_restart_objects
!!***

!!****f* module_types/free_restart_objects
!! FUNCTION
!!   De-Allocate restart_objects
!! SOURCE
!!
  subroutine free_restart_objects(rst,subname)
    use module_base
    implicit none
    character(len=*), intent(in) :: subname
    type(restart_objects) :: rst
    !local variables
    integer :: i_all,i_stat

    call deallocate_wfd(rst%Glr%wfd,subname)

    i_all=-product(shape(rst%psi))*kind(rst%psi)
    deallocate(rst%psi,stat=i_stat)
    call memocc(i_stat,i_all,'psi',subname)
    i_all=-product(shape(rst%orbs%eval))*kind(rst%orbs%eval)
    deallocate(rst%orbs%eval,stat=i_stat)
    call memocc(i_stat,i_all,'eval',subname)
    i_all=-product(shape(rst%rxyz_old))*kind(rst%rxyz_old)
    deallocate(rst%rxyz_old,stat=i_stat)
    call memocc(i_stat,i_all,'rxyz_old',subname)

    !the gaussian basis descriptors are always allocated together
    !with the gaussian coefficients
    if (associated(rst%gbd%rxyz)) then
       nullify(rst%gbd%rxyz)
       call deallocate_gwf(rst%gbd,subname)

       i_all=-product(shape(rst%gaucoeffs))*kind(rst%gaucoeffs)
       deallocate(rst%gaucoeffs,stat=i_stat)
       call memocc(i_stat,i_all,'gaucoeffs',subname)

    end if

  end subroutine free_restart_objects
!!***

!!****f* module_types/allocate_wfd
!! FUNCTION
!!   Allocate wavefunctions_descriptors
!! SOURCE
!!
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
  end subroutine allocate_wfd
!!***

!!****f* module_types/deallocate_wfd
!! FUNCTION
!!   De-Allocate wavefunctions_descriptors
!! SOURCE
!!
  subroutine deallocate_wfd(wfd,subname)
    use module_base
    implicit none
    type(wavefunctions_descriptors) :: wfd
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(wfd%keyg))*kind(wfd%keyg)
    deallocate(wfd%keyg,stat=i_stat)
    call memocc(i_stat,i_all,'wfd%keyg',subname)
    i_all=-product(shape(wfd%keyv))*kind(wfd%keyv)
    deallocate(wfd%keyv,stat=i_stat)
    call memocc(i_stat,i_all,'wfd%keyv',subname)

  end subroutine deallocate_wfd
!!***

!!****f* module_types/deallocate_gwf
!! FUNCTION
!!   De-Allocate gaussian_basis type
!! SOURCE
!!
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

  end subroutine deallocate_gwf
!!***


!!****f* module_types/deallocate_bounds
!! FUNCTION
!!   De-Allocate convolutions_bounds type
!! SOURCE
!!
  subroutine deallocate_bounds(bounds,subname)
    use module_base
    implicit none
    type(convolutions_bounds) :: bounds
    character(len=*), intent(in) :: subname
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(bounds%kb%ibyz_c))*kind(bounds%kb%ibyz_c)
    deallocate(bounds%kb%ibyz_c,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%kb%ibyz_c',subname)
    i_all=-product(shape(bounds%kb%ibxz_c))*kind(bounds%kb%ibxz_c)
    deallocate(bounds%kb%ibxz_c,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%kb%ibxz_c',subname)
    i_all=-product(shape(bounds%kb%ibxy_c))*kind(bounds%kb%ibxy_c)
    deallocate(bounds%kb%ibxy_c,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%kb%ibxy_c',subname)
    i_all=-product(shape(bounds%kb%ibyz_f))*kind(bounds%kb%ibyz_f)
    deallocate(bounds%kb%ibyz_f,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%kb%ibyz_f',subname)
    i_all=-product(shape(bounds%kb%ibxz_f))*kind(bounds%kb%ibxz_f)
    deallocate(bounds%kb%ibxz_f,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%kb%ibxz_f',subname)
    i_all=-product(shape(bounds%kb%ibxy_f))*kind(bounds%kb%ibxy_f)
    deallocate(bounds%kb%ibxy_f,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%kb%ibxy_f',subname)

    i_all=-product(shape(bounds%sb%ibzzx_c))*kind(bounds%sb%ibzzx_c)
    deallocate(bounds%sb%ibzzx_c,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%sb%ibzzx_c',subname)
    i_all=-product(shape(bounds%sb%ibyyzz_c))*kind(bounds%sb%ibyyzz_c)
    deallocate(bounds%sb%ibyyzz_c,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%sb%ibyyzz_c',subname)
    i_all=-product(shape(bounds%sb%ibxy_ff))*kind(bounds%sb%ibxy_ff)
    deallocate(bounds%sb%ibxy_ff,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%sb%ibxy_ff',subname)
    i_all=-product(shape(bounds%sb%ibzzx_f))*kind(bounds%sb%ibzzx_f)
    deallocate(bounds%sb%ibzzx_f,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%sb%ibzzx_f',subname)
    i_all=-product(shape(bounds%sb%ibyyzz_f))*kind(bounds%sb%ibyyzz_f)
    deallocate(bounds%sb%ibyyzz_f,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%sb%ibyyzz_f',subname)

    i_all=-product(shape(bounds%gb%ibzxx_c))*kind(bounds%gb%ibzxx_c)
    deallocate(bounds%gb%ibzxx_c,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%gb%ibzxx_c',subname)
    i_all=-product(shape(bounds%gb%ibxxyy_c))*kind(bounds%gb%ibxxyy_c)
    deallocate(bounds%gb%ibxxyy_c,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%gb%ibxxyy_c',subname)
    i_all=-product(shape(bounds%gb%ibyz_ff))*kind(bounds%gb%ibyz_ff)
    deallocate(bounds%gb%ibyz_ff,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%gb%ibyz_ff',subname)
    i_all=-product(shape(bounds%gb%ibzxx_f))*kind(bounds%gb%ibzxx_f)
    deallocate(bounds%gb%ibzxx_f,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%gb%ibzxx_f',subname)
    i_all=-product(shape(bounds%gb%ibxxyy_f))*kind(bounds%gb%ibxxyy_f)
    deallocate(bounds%gb%ibxxyy_f,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%gb%ibxxyy_f',subname)

    i_all=-product(shape(bounds%ibyyzz_r))*kind(bounds%ibyyzz_r)
    deallocate(bounds%ibyyzz_r,stat=i_stat)
    call memocc(i_stat,i_all,'bounds%ibyyzz_r',subname)

  end subroutine deallocate_bounds

end module module_types
!!***
