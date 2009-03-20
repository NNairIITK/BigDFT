!!****m* BigDFT/module_types
!! FUNCTION
!!  Modules which contains the Fortran data structures
!! AUTHOR
!!    Luigi Genovese
!! COPYRIGHT
!!    Copyright (C) 2008 CEA
!! SOURCE
!! 
module module_types
  use module_base, only : gp,wp,dp
  implicit none
!!***

!!****t* module_types/input_variables
!! DESCRIPTION
!!   Input variable structure
!!   Structure of the variables read by input.dat file
!! SOURCE
!!
  type, public :: input_variables
     logical :: output_wf,calc_tail,gaussian_help
     integer :: ncount_cluster_x
     integer :: ixc,ncharge,itermax,nrepmax,ncong,idsx,ncongt,inputPsiId,nspin,mpol,nvirt,nplot
     integer :: output_grid
     real(gp) :: frac_fluct,randdis,betax,forcemax
     real(gp) :: hgrid,crmult,frmult,gnrm_cv,rbuf
     real(gp), dimension(3) :: ef
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
!! SOURCE
!!
  type, public :: atoms_data
     character(len=1) :: geocode
     character(len=20) :: units
     integer :: nat,ntypes,natsc
     character(len=20), dimension(:), pointer :: atomnames
     real(gp) :: alat1,alat2,alat3
     logical, dimension(:), pointer :: lfrztyp
     integer, dimension(:), pointer :: iatype,iasctype,natpol,nelpsp,npspcode,nzatom
     real(gp), dimension(:,:,:), pointer :: psppar
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
!!
!! SOURCE
!!
  type, public :: orbitals_data
     integer :: norb,norbp,norbu,norbd,nspinor,isorb,npsidim
     integer, dimension(:), pointer :: norb_par
     real(wp), dimension(:), pointer :: eval
     real(gp), dimension(:), pointer :: occup,spinsgn
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
     real(gp) :: hgrid_old
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
  end type communications_arrays
!!***

!!****t* module_types/GPU_pointers
!! DESCRIPTION
!! Contains the pointers to be handled to control GPU information
!! Given that they are pointers on GPU address, they are C pointers
!! which takes 8 bytes
!! So they are declared as kind=8 variables either if the GPU works in simple precision
!! Also other information concerning the GPU runs can be stored in this structure
!!
!! SOURCE
!!
  type, public :: GPU_pointers
     real(kind=8) :: keys,work1,work2,work3,rhopot,r,d
     real(kind=8), dimension(:), pointer :: psi
  end type GPU_pointers
!!***

contains

  subroutine allocate_comms(nproc,comms,routine)
    use module_base
    implicit none
    character(len=*), intent(in) :: routine
    integer, intent(in) :: nproc
    type(communications_arrays), intent(out) :: comms
    !local variables
    integer :: i_all,i_stat

    allocate(comms%ncntd(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ncntd,'ncntd',routine)
    allocate(comms%ncntt(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ncntt,'ncntt',routine)
    allocate(comms%ndspld(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ndspld,'ndspld',routine)
    allocate(comms%ndsplt(0:nproc-1+ndebug),stat=i_stat)
    call memocc(i_stat,comms%ndsplt,'ndsplt',routine)
  end subroutine allocate_comms

  subroutine deallocate_comms(comms,routine)
    use module_base
    implicit none
    character(len=*), intent(in) :: routine
    type(communications_arrays), intent(out) :: comms
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(comms%ncntd))*kind(comms%ncntd)
    deallocate(comms%ncntd,stat=i_stat)
    call memocc(i_stat,i_all,'ncntd',routine)
    i_all=-product(shape(comms%ncntt))*kind(comms%ncntt)
    deallocate(comms%ncntt,stat=i_stat)
    call memocc(i_stat,i_all,'ncntt',routine)
    i_all=-product(shape(comms%ndspld))*kind(comms%ndspld)
    deallocate(comms%ndspld,stat=i_stat)
    call memocc(i_stat,i_all,'ndspld',routine)
    i_all=-product(shape(comms%ndsplt))*kind(comms%ndsplt)
    deallocate(comms%ndsplt,stat=i_stat)
    call memocc(i_stat,i_all,'ndsplt',routine)
  end subroutine deallocate_comms

  subroutine init_restart_objects(atoms,rst,routine)
    use module_base
    implicit none
    character(len=*), intent(in) :: routine
    type(atoms_data) :: atoms
    type(restart_objects) :: rst
    !local variables
    integer :: i_all,i_stat

    !allocate pointers
    allocate(rst%rxyz_old(3,atoms%nat+ndebug),stat=i_stat)
    call memocc(i_stat,rst%rxyz_old,'rxyz_old',routine)

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

  subroutine free_restart_objects(rst,routine)
    use module_base
    implicit none
    character(len=*), intent(in) :: routine
    type(restart_objects) :: rst
    !local variables
    integer :: i_all,i_stat

    call deallocate_wfd(rst%Glr%wfd,routine)

    i_all=-product(shape(rst%psi))*kind(rst%psi)
    deallocate(rst%psi,stat=i_stat)
    call memocc(i_stat,i_all,'psi',routine)
    i_all=-product(shape(rst%orbs%eval))*kind(rst%orbs%eval)
    deallocate(rst%orbs%eval,stat=i_stat)
    call memocc(i_stat,i_all,'eval',routine)
    i_all=-product(shape(rst%rxyz_old))*kind(rst%rxyz_old)
    deallocate(rst%rxyz_old,stat=i_stat)
    call memocc(i_stat,i_all,'rxyz_old',routine)

    !the gaussian basis descriptors are always allocated together
    !with the gaussian coefficients
    if (associated(rst%gbd%rxyz)) then
       nullify(rst%gbd%rxyz)
       call deallocate_gwf(rst%gbd,routine)

       i_all=-product(shape(rst%gaucoeffs))*kind(rst%gaucoeffs)
       deallocate(rst%gaucoeffs,stat=i_stat)
       call memocc(i_stat,i_all,'gaucoeffs',routine)

    end if
       

  end subroutine free_restart_objects

  subroutine allocate_wfd(wfd,routine)
    use module_base
    implicit none
    type(wavefunctions_descriptors), intent(inout) :: wfd
    character(len=*), intent(in) :: routine
    !local variables
    integer :: i_all,i_stat

    allocate(wfd%keyg(2,wfd%nseg_c+wfd%nseg_f+ndebug),stat=i_stat)
    call memocc(i_stat,wfd%keyg,'keyg',routine)
    allocate(wfd%keyv(wfd%nseg_c+wfd%nseg_f+ndebug),stat=i_stat)
    call memocc(i_stat,wfd%keyv,'keyv',routine)
  end subroutine allocate_wfd

  subroutine deallocate_wfd(wfd,routine)
    use module_base
    implicit none
    type(wavefunctions_descriptors) :: wfd
    character(len=*), intent(in) :: routine
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(wfd%keyg))*kind(wfd%keyg)
    deallocate(wfd%keyg,stat=i_stat)
    call memocc(i_stat,i_all,'keyg',routine)
    i_all=-product(shape(wfd%keyv))*kind(wfd%keyv)
    deallocate(wfd%keyv,stat=i_stat)
    call memocc(i_stat,i_all,'keyv',routine)

  end subroutine deallocate_wfd

  subroutine deallocate_gwf(G,routine)
    use module_base
    implicit none
    type(gaussian_basis) :: G
    character(len=*), intent(in) :: routine
    !local variables
    integer :: i_all,i_stat

    !normally positions should be deallocated outside
    
    i_all=-product(shape(G%ndoc))*kind(G%ndoc)
    deallocate(G%ndoc,stat=i_stat)
    call memocc(i_stat,i_all,'ndoc',routine)
    i_all=-product(shape(G%nam))*kind(G%nam)
    deallocate(G%nam,stat=i_stat)
    call memocc(i_stat,i_all,'nam',routine)
    i_all=-product(shape(G%nshell))*kind(G%nshell)
    deallocate(G%nshell,stat=i_stat)
    call memocc(i_stat,i_all,'nshell',routine)
    i_all=-product(shape(G%psiat))*kind(G%psiat)
    deallocate(G%psiat,stat=i_stat)
    call memocc(i_stat,i_all,'psiat',routine)
    i_all=-product(shape(G%xp))*kind(G%xp)
    deallocate(G%xp,stat=i_stat)
    call memocc(i_stat,i_all,'xp',routine)

  end subroutine deallocate_gwf


  subroutine deallocate_bounds(bounds,routine)
    use module_base
    implicit none
    type(convolutions_bounds) :: bounds
    character(len=*), intent(in) :: routine
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(bounds%kb%ibyz_c))*kind(bounds%kb%ibyz_c)
    deallocate(bounds%kb%ibyz_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibyz_c',routine)
    i_all=-product(shape(bounds%kb%ibxz_c))*kind(bounds%kb%ibxz_c)
    deallocate(bounds%kb%ibxz_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibxz_c',routine)
    i_all=-product(shape(bounds%kb%ibxy_c))*kind(bounds%kb%ibxy_c)
    deallocate(bounds%kb%ibxy_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibxy_c',routine)
    i_all=-product(shape(bounds%kb%ibyz_f))*kind(bounds%kb%ibyz_f)
    deallocate(bounds%kb%ibyz_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibyz_f',routine)
    i_all=-product(shape(bounds%kb%ibxz_f))*kind(bounds%kb%ibxz_f)
    deallocate(bounds%kb%ibxz_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibxz_f',routine)
    i_all=-product(shape(bounds%kb%ibxy_f))*kind(bounds%kb%ibxy_f)
    deallocate(bounds%kb%ibxy_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibxy_f',routine)

    i_all=-product(shape(bounds%sb%ibzzx_c))*kind(bounds%sb%ibzzx_c)
    deallocate(bounds%sb%ibzzx_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibzzx_c',routine)
    i_all=-product(shape(bounds%sb%ibyyzz_c))*kind(bounds%sb%ibyyzz_c)
    deallocate(bounds%sb%ibyyzz_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_c',routine)
    i_all=-product(shape(bounds%sb%ibxy_ff))*kind(bounds%sb%ibxy_ff)
    deallocate(bounds%sb%ibxy_ff,stat=i_stat)
    call memocc(i_stat,i_all,'ibxy_ff',routine)
    i_all=-product(shape(bounds%sb%ibzzx_f))*kind(bounds%sb%ibzzx_f)
    deallocate(bounds%sb%ibzzx_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibzzx_f',routine)
    i_all=-product(shape(bounds%sb%ibyyzz_f))*kind(bounds%sb%ibyyzz_f)
    deallocate(bounds%sb%ibyyzz_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_f',routine)

    i_all=-product(shape(bounds%gb%ibzxx_c))*kind(bounds%gb%ibzxx_c)
    deallocate(bounds%gb%ibzxx_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibzxx_c',routine)
    i_all=-product(shape(bounds%gb%ibxxyy_c))*kind(bounds%gb%ibxxyy_c)
    deallocate(bounds%gb%ibxxyy_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibxxyy_c',routine)
    i_all=-product(shape(bounds%gb%ibyz_ff))*kind(bounds%gb%ibyz_ff)
    deallocate(bounds%gb%ibyz_ff,stat=i_stat)
    call memocc(i_stat,i_all,'ibyz_ff',routine)
    i_all=-product(shape(bounds%gb%ibzxx_f))*kind(bounds%gb%ibzxx_f)
    deallocate(bounds%gb%ibzxx_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibzxx_f',routine)
    i_all=-product(shape(bounds%gb%ibxxyy_f))*kind(bounds%gb%ibxxyy_f)
    deallocate(bounds%gb%ibxxyy_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibxxyy_f',routine)

    i_all=-product(shape(bounds%ibyyzz_r))*kind(bounds%ibyyzz_r)
    deallocate(bounds%ibyyzz_r,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_r',routine)

  end subroutine deallocate_bounds

end module module_types
