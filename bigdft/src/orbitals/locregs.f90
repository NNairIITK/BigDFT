!> @file
!! Datatypes and associated methods relative to the localization regions (mesh grid)
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Datatypes for localization regions descriptors
module locregs
  use module_base
  implicit none

  private 

  !> Bounds for coarse and fine grids for kinetic operations
  !! Useful only for isolated systems AND in CPU
  type, public :: kinetic_bounds
     integer, dimension(:,:,:), pointer :: ibyz_c !< coarse (2,0:n2,0:n3)
     integer, dimension(:,:,:), pointer :: ibxz_c !< coarse (2,0:n1,0:n3)
     integer, dimension(:,:,:), pointer :: ibxy_c !< coarse (2,0:n1,0:n2)
     integer, dimension(:,:,:), pointer :: ibyz_f !< fine (2,0:n2,0:n3)
     integer, dimension(:,:,:), pointer :: ibxz_f !< fine (2,0:n1,0:n3)
     integer, dimension(:,:,:), pointer :: ibxy_f !< fine (2,0:n1,0:n2)
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
     integer :: n1,n2,n3                      !< Coarse grid dimensions
     integer :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3 !< Lower and upper indices of fine grid in 3D 
     integer :: n1i,n2i,n3i                   !< ISF grid dimension (roughly 2*n+buffer)
  end type grid_dimensions


  !> Contains the information needed for describing completely a wavefunction localisation region
  type, public :: locreg_descriptors
     character(len=1) :: geocode            !< @copydoc poisson_solver::doc::geocode
     logical :: hybrid_on                   !< Interesting for global, periodic, localisation regions
     integer :: ns1,ns2,ns3                 !< Starting point of the localisation region in global coordinates
     integer :: nsi1,nsi2,nsi3              !< Starting point of locreg for interpolating grid
     integer :: Localnorb                   !< Number of orbitals contained in locreg
     integer, dimension(3) :: outofzone     !< Vector of points outside of the zone outside Glr for periodic systems
     real(gp), dimension(3) :: locregCenter !< Center of the locreg 
     real(gp) :: locrad                     !< Cutoff radius of the localization region
     real(gp) :: locrad_kernel              !< Cutoff radius of the localization region (kernel)
     real(gp) :: locrad_mult                !< Cutoff radius of the localization region for the sparse matrix multiplications
     type(grid_dimensions) :: d             !< Grid dimensions in old different wavelet basis
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
  end type locreg_descriptors


  public :: nullify_locreg_descriptors,locreg_null
  public :: deallocate_locreg_descriptors,deallocate_wfd
  public :: allocate_wfd,copy_locreg_descriptors,copy_grid_dimensions,nullify_wfd
  public :: check_overlap,check_overlap_cubic_periodic,check_overlap_from_descriptors_periodic
  public :: check_whether_bounds_overlap
  public :: get_extent_of_overlap


  interface check_whether_bounds_overlap
    module procedure check_whether_bounds_overlap_int
    module procedure check_whether_bounds_overlap_long
  end interface check_whether_bounds_overlap
  
  interface get_extent_of_overlap
    module procedure get_extent_of_overlap_int
    module procedure get_extent_of_overlap_long
  end interface get_extent_of_overlap


contains
  

  !constructors
  pure function convolutions_bounds_null() result(bounds)
    implicit none
    type(convolutions_bounds) :: bounds
    call nullify_convolutions_bounds(bounds)
  end function convolutions_bounds_null

  pure subroutine nullify_convolutions_bounds(bounds)
    implicit none
    type(convolutions_bounds), intent(out) :: bounds
    call nullify_kinetic_bounds(bounds%kb)
    call nullify_shrink_bounds(bounds%sb)
    call nullify_grow_bounds(bounds%gb)
    nullify(bounds%ibyyzz_r)
  end subroutine nullify_convolutions_bounds

  pure function kinetic_bounds_null() result(kb)
    implicit none
    type(kinetic_bounds) :: kb
    call nullify_kinetic_bounds(kb)
  end function kinetic_bounds_null

  pure subroutine nullify_kinetic_bounds(kb)
    implicit none
    type(kinetic_bounds), intent(out) :: kb
    nullify(kb%ibyz_c)
    nullify(kb%ibxz_c)
    nullify(kb%ibxy_c)
    nullify(kb%ibyz_f)
    nullify(kb%ibxz_f)
    nullify(kb%ibxy_f)
  end subroutine nullify_kinetic_bounds

  pure function shrink_bounds_null() result(sb)
    implicit none
    type(shrink_bounds) :: sb
    call nullify_shrink_bounds(sb)
  end function shrink_bounds_null

  pure subroutine nullify_shrink_bounds(sb)
    implicit none
    type(shrink_bounds), intent(out) :: sb
    nullify(sb%ibzzx_c)
    nullify(sb%ibyyzz_c)
    nullify(sb%ibxy_ff)
    nullify(sb%ibzzx_f)
    nullify(sb%ibyyzz_f)
  end subroutine nullify_shrink_bounds

  pure function grow_bounds_null() result(gb)
    implicit none
    type(grow_bounds) :: gb
    call nullify_grow_bounds(gb)
  end function grow_bounds_null

  pure subroutine nullify_grow_bounds(gb)
    implicit none
    type(grow_bounds), intent(out) :: gb
    nullify(gb%ibzxx_c)
    nullify(gb%ibxxyy_c)
    nullify(gb%ibyz_ff)
    nullify(gb%ibzxx_f)
    nullify(gb%ibxxyy_f)
  end subroutine nullify_grow_bounds

  pure function grid_null() result(g)
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
  end function grid_null

  pure function wfd_null() result(wfd)
    implicit none
    type(wavefunctions_descriptors) :: wfd
    call nullify_wfd(wfd)
  end function wfd_null

  pure subroutine nullify_wfd(wfd)
    implicit none
    type(wavefunctions_descriptors), intent(out) :: wfd
    wfd%nvctr_c=0
    wfd%nvctr_f=0
    wfd%nseg_c=0
    wfd%nseg_f=0
    nullify(wfd%keyglob)
    nullify(wfd%keygloc)
    nullify(wfd%keyvglob)
    nullify(wfd%keyvloc)
  end subroutine nullify_wfd

  pure function locreg_null() result(lr)
    implicit none
    type(locreg_descriptors) :: lr
    call nullify_locreg_descriptors(lr)
  end function locreg_null

  pure subroutine nullify_locreg_descriptors(lr)
    implicit none
    type(locreg_descriptors), intent(out) :: lr
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
    lr%d=grid_null()
    call nullify_wfd(lr%wfd)
    call nullify_convolutions_bounds(lr%bounds)
    lr%locregCenter=(/0.0_gp,0.0_gp,0.0_gp/) 
    lr%locrad=0 
  end subroutine nullify_locreg_descriptors

  !> Initializations
  subroutine allocate_wfd(wfd)
    use module_base
    implicit none
    type(wavefunctions_descriptors), intent(inout) :: wfd
    !local variables
    integer :: nsegs

    nsegs=max(1,wfd%nseg_c+wfd%nseg_f)
    wfd%keyvloc=f_malloc_ptr(nsegs,id='wfd%keyvloc')
    wfd%keyvglob=f_malloc_ptr(nsegs,id='wfd%keyvglob')
    wfd%keyglob=f_malloc_ptr((/2,nsegs/),id='wfd%keyglob')
    wfd%keygloc=f_malloc_ptr((/2,nsegs/),id='wfd%keygloc')
  END SUBROUTINE allocate_wfd

  !> De-Allocate wavefunctions_descriptors
  subroutine deallocate_wfd(wfd)
    use module_base
    implicit none
    type(wavefunctions_descriptors), intent(inout) :: wfd

    !in case the two objects points to the same target
    if (associated(wfd%keyglob, target = wfd%keygloc)) then
       !assuming that globals has been created afterwards
       nullify(wfd%keygloc)
       call f_free_ptr(wfd%keyglob)
    else
       call f_free_ptr(wfd%keygloc)
       call f_free_ptr(wfd%keyglob)
    end if
    if (associated(wfd%keyvloc, target= wfd%keyvglob)) then
       nullify(wfd%keyvloc)
       call f_free_ptr(wfd%keyvglob)
    else
       call f_free_ptr(wfd%keyvloc)
       call f_free_ptr(wfd%keyvglob)
    end if
  END SUBROUTINE deallocate_wfd

  !> Destructors
  subroutine deallocate_locreg_descriptors(lr)
    implicit none
    ! Calling arguments
    type(locreg_descriptors),intent(inout):: lr

    call deallocate_wfd(lr%wfd)
    call deallocate_convolutions_bounds(lr%bounds)
    
  end subroutine deallocate_locreg_descriptors

  subroutine deallocate_convolutions_bounds(bounds)
    implicit none
    type(convolutions_bounds),intent(inout):: bounds

    call f_free_ptr(bounds%ibyyzz_r)

    call deallocate_kinetic_bounds(bounds%kb)
    call deallocate_shrink_bounds(bounds%sb)
    call deallocate_grow_bounds(bounds%gb)

  end subroutine deallocate_convolutions_bounds

  subroutine deallocate_kinetic_bounds(kb)
    implicit none
    ! Calling arguments
    type(kinetic_bounds),intent(inout):: kb

    call f_free_ptr(kb%ibyz_c)
    call f_free_ptr(kb%ibxz_c)
    call f_free_ptr(kb%ibxy_c)
    call f_free_ptr(kb%ibyz_f)
    call f_free_ptr(kb%ibxz_f)
    call f_free_ptr(kb%ibxy_f)

  end subroutine deallocate_kinetic_bounds


  subroutine deallocate_shrink_bounds(sb)
    implicit none
    ! Calling arguments
    type(shrink_bounds),intent(inout):: sb

    call f_free_ptr(sb%ibzzx_c)
    call f_free_ptr(sb%ibyyzz_c)
    call f_free_ptr(sb%ibxy_ff)
    call f_free_ptr(sb%ibzzx_f)
    call f_free_ptr(sb%ibyyzz_f)

  end subroutine deallocate_shrink_bounds


  subroutine deallocate_grow_bounds(gb)
    implicit none
    ! Calling arguments
    type(grow_bounds),intent(inout):: gb

    call f_free_ptr(gb%ibzxx_c)
    call f_free_ptr(gb%ibxxyy_c)
    call f_free_ptr(gb%ibyz_ff)
    call f_free_ptr(gb%ibzxx_f)
    call f_free_ptr(gb%ibxxyy_f)

  end subroutine deallocate_grow_bounds

  !> Methods for copying the structures, can be needed to avoid recalculating them
  !! should be better by defining a f_malloc inheriting the shapes and the structure from other array
  !! of the type dest=f_malloc(src=source,id='dest')
  subroutine copy_locreg_descriptors(glrin, glrout)
    implicit none
    ! Calling arguments
    type(locreg_descriptors), intent(in) :: glrin !<input locreg. Unchanged on exit.
    type(locreg_descriptors), intent(out):: glrout !<output locreg. Must be freed on input.

    glrout%geocode = glrin%geocode
    glrout%hybrid_on = glrin%hybrid_on
    glrout%ns1 = glrin%ns1
    glrout%ns2 = glrin%ns2
    glrout%ns3 = glrin%ns3
    glrout%nsi1 = glrin%nsi1
    glrout%nsi2 = glrin%nsi2
    glrout%nsi3 = glrin%nsi3
    glrout%Localnorb = glrin%Localnorb
    glrout%locrad=glrin%locrad
    glrout%locrad_kernel=glrin%locrad_kernel
    glrout%locrad_mult=glrin%locrad_mult
    glrout%locregCenter=glrin%locregCenter
    glrout%outofzone= glrin%outofzone

    call copy_grid_dimensions(glrin%d, glrout%d)
    call copy_wavefunctions_descriptors(glrin%wfd, glrout%wfd)
    call copy_convolutions_bounds(glrin%bounds, glrout%bounds)

!!$    !copy bound when needed
!!$    if(glrin%geocode == 'F' .or. (glrin%geocode == 'P' .and. glrin%hybrid_on)) then
!!$       call copy_convolutions_bounds(glrin%geocode, glrin%bounds, glrout%bounds,&
!!$            'copy_locreg_descriptors') !to be removed when bounds are allocated properly
!!$    else
!!$       call nullify_convolutions_bounds(glrout%bounds)
!!$    end if

  end subroutine copy_locreg_descriptors
  pure subroutine copy_grid_dimensions(din, dout)
    implicit none
    ! Calling arguments
    type(grid_dimensions),intent(in):: din
    type(grid_dimensions),intent(out):: dout

    dout%n1 = din%n1
    dout%n2 = din%n2
    dout%n3 = din%n3
    dout%nfl1 = din%nfl1
    dout%nfu1 = din%nfu1
    dout%nfl2 = din%nfl2
    dout%nfu2 = din%nfu2
    dout%nfl3 = din%nfl3
    dout%nfu3 = din%nfu3
    dout%n1i = din%n1i
    dout%n2i = din%n2i
    dout%n3i = din%n3i

  end subroutine copy_grid_dimensions

  subroutine copy_wavefunctions_descriptors(wfdin, wfdout)
    implicit none
    ! Calling arguments
    type(wavefunctions_descriptors), intent(in) :: wfdin
    type(wavefunctions_descriptors), intent(out) :: wfdout

    ! Local variables
    !integer:: istat,iis1, iie1, iis2, iie2,i1, i2, iall

    !nullify all pointers first
    call nullify_wfd(wfdout)

    wfdout%nvctr_c = wfdin%nvctr_c
    wfdout%nvctr_f = wfdin%nvctr_f
    wfdout%nseg_c = wfdin%nseg_c
    wfdout%nseg_f = wfdin%nseg_f

    !new method
    wfdout%keygloc=f_malloc_ptr(src_ptr=wfdin%keygloc,id='wfdout%keygloc')
    wfdout%keyglob=f_malloc_ptr(src_ptr=wfdin%keyglob,id='wfdout%keyglob')
    wfdout%keyvloc=f_malloc_ptr(src_ptr=wfdin%keyvloc,id='wfdout%keyvloc')
    wfdout%keyvglob=f_malloc_ptr(src_ptr=wfdin%keyvglob,id='wfdout%keyvglob')

!!$    !no need to insert lbounds as the allocation start from 1
!!$    if (associated(wfdin%keygloc)) wfdout%keygloc=f_malloc_ptr(src=wfdin%keygloc,id='wfdout%keygloc')
!!$    if (associated(wfdin%keyglob)) wfdout%keyglob=f_malloc_ptr(src=wfdin%keyglob,id='wfdout%keyglob')
!!$    if (associated(wfdin%keyvloc)) wfdout%keyvloc=f_malloc_ptr(src=wfdin%keyvloc,id='wfdout%keyvloc')
!!$    if (associated(wfdin%keyvglob))wfdout%keyvglob=f_malloc_ptr(src=wfdin%keyvglob,id='wfdout%keyvglob')


  end subroutine copy_wavefunctions_descriptors


  subroutine copy_convolutions_bounds(boundsin, boundsout)
    implicit none
    type(convolutions_bounds),intent(in):: boundsin
    type(convolutions_bounds),intent(inout):: boundsout

    call copy_kinetic_bounds(boundsin%kb, boundsout%kb)
    call copy_shrink_bounds(boundsin%sb, boundsout%sb)
    call copy_grow_bounds(boundsin%gb, boundsout%gb)
    boundsout%ibyyzz_r = f_malloc_ptr(src_ptr=boundsin%ibyyzz_r,id='boundsout%ibyyzz_r')

  end subroutine copy_convolutions_bounds

  subroutine copy_kinetic_bounds(kbin, kbout)
    implicit none
    ! Calling arguments
    type(kinetic_bounds),intent(in):: kbin
    type(kinetic_bounds),intent(inout):: kbout

    kbout%ibyz_c = f_malloc_ptr(src_ptr=kbin%ibyz_c,id='kbout%ibyz_c')
    kbout%ibxz_c = f_malloc_ptr(src_ptr=kbin%ibxz_c,id='kbout%ibxz_c')
    kbout%ibxy_c = f_malloc_ptr(src_ptr=kbin%ibxy_c,id='kbout%ibxy_c')
    kbout%ibyz_f = f_malloc_ptr(src_ptr=kbin%ibyz_f,id='kbout%ibyz_f')
    kbout%ibxz_f = f_malloc_ptr(src_ptr=kbin%ibxz_f,id='kbout%ibxz_f')
    kbout%ibxy_f = f_malloc_ptr(src_ptr=kbin%ibxy_f,id='kbout%ibxy_f')

  end subroutine copy_kinetic_bounds



  subroutine copy_shrink_bounds(sbin, sbout)
    implicit none

    ! Calling arguments
    type(shrink_bounds),intent(in):: sbin
    type(shrink_bounds),intent(inout):: sbout

    sbout%ibzzx_c = f_malloc_ptr(src_ptr=sbin%ibzzx_c,id='sbout%ibzzx_c')
    sbout%ibyyzz_c= f_malloc_ptr(src_ptr=sbin%ibyyzz_c,id='sbout%ibyyzz_c')
    sbout%ibxy_ff = f_malloc_ptr(src_ptr=sbin%ibxy_ff,id='sbout%ibxy_ff')
    sbout%ibzzx_f = f_malloc_ptr(src_ptr=sbin%ibzzx_f,id='sbout%ibzzx_f')
    sbout%ibyyzz_f= f_malloc_ptr(src_ptr=sbin%ibyyzz_f,id='sbout%ibyyzz_f')

  end subroutine copy_shrink_bounds


  subroutine copy_grow_bounds(gbin, gbout)
    implicit none

    ! Calling arguments
    type(grow_bounds),intent(in):: gbin
    type(grow_bounds),intent(inout):: gbout

    gbout%ibzxx_c = f_malloc_ptr(src_ptr=gbin%ibzxx_c,id='gbout%ibzxx_c')
    gbout%ibxxyy_c= f_malloc_ptr(src_ptr=gbin%ibxxyy_c,id='gbout%ibxxyy_c')
    gbout%ibyz_ff = f_malloc_ptr(src_ptr=gbin%ibyz_ff,id='gbout%ibyz_ff')
    gbout%ibzxx_f = f_malloc_ptr(src_ptr=gbin%ibzxx_f,id='gbout%ibzxx_f')
    gbout%ibxxyy_f= f_malloc_ptr(src_ptr=gbin%ibxxyy_f,id='gbout%ibxxyy_f')

  end subroutine copy_grow_bounds


  !> Almost degenerate with get_number_of_overlap_region
  !! should merge the two... prefering this one since argument list is better 
  subroutine check_overlap_cubic_periodic(Glr,Ilr,Jlr,isoverlap)
    use module_base
    !use communications_init, only: check_whether_bounds_overlap
    implicit none
    type(locreg_descriptors), intent(in) :: Glr
    type(locreg_descriptors), intent(in) :: Ilr
    type(locreg_descriptors), intent(in) :: Jlr
    logical, intent(out) :: isoverlap
    !Local variables
    integer :: is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
    logical :: overlap1, overlap2, overlap3
  !!  integer :: azones,bzones,ii,izones,jzones !, i_stat, i_all
  !!  logical :: go1, go2, go3
  !!  integer,dimension(3,8) :: astart,bstart,aend,bend
  
  !!  azones = 1
  !!  bzones = 1
  !!! Calculate the number of regions to cut alr and blr
  !!  do ii=1,3
  !!     if(Ilr%outofzone(ii) > 0) azones = azones * 2
  !!     if(Jlr%outofzone(ii) > 0) bzones = bzones * 2
  !!  end do
  !!
  !!!FRACTURE THE FIRST LOCALIZATION REGION
  !!  call fracture_periodic_zone(azones,Glr,Ilr,Ilr%outofzone,astart,aend)
  !!
  !!!FRACTURE SECOND LOCREG
  !!  call fracture_periodic_zone(bzones,Glr,Jlr,Jlr%outofzone,bstart,bend)
  !!
  !!! Now check if they overlap
  !!  isoverlap = .false.
  !!  loop_izones: do izones=1,azones
  !!    do jzones=1,bzones
  !!      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones))
  !!      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones))
  !!      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones))
  !!      if(go1 .and. go2 .and. go3) then
  !!        isoverlap = .true.
  !!        exit loop_izones
  !!      end if
  !!    end do
  !!  end do loop_izones
  
  
    !@ NEW VERSION #########################################
    ! Shift all the indices into the periodic cell. This can result is starting
    ! indices being larger than ending indices
    isoverlap = .false.
    is3 = modulo(ilr%ns3,glr%d%n3+1)
    ie3 = modulo(ilr%ns3+ilr%d%n3,glr%d%n3+1)
    js3 = modulo(jlr%ns3,glr%d%n3+1)
    je3 = modulo(jlr%ns3+jlr%d%n3,glr%d%n3+1)
    overlap3 = check_whether_bounds_overlap(is3, ie3, js3, je3)
    if (overlap3) then
        is2 = modulo(ilr%ns2,glr%d%n2+1)
        ie2 = modulo(ilr%ns2+ilr%d%n2,glr%d%n2+1)
        js2 = modulo(jlr%ns2,glr%d%n2+1)
        je2 = modulo(jlr%ns2+jlr%d%n2,glr%d%n2+1)
        overlap2 = check_whether_bounds_overlap(is2, ie2, js2, je2)
        if (overlap2) then
            is1 = modulo(ilr%ns1,glr%d%n1+1)
            ie1 = modulo(ilr%ns1+ilr%d%n1,glr%d%n1+1)
            js1 = modulo(jlr%ns1,glr%d%n1+1)
            je1 = modulo(jlr%ns1+jlr%d%n1,glr%d%n1+1)
            overlap1 = check_whether_bounds_overlap(is1, ie1, js1, je1)
            if (overlap1) then
                ! If we are here, all three overlaps are true
                isoverlap = .true.
            end if
        end if
    end if
  
    !!if (overlap1 .and. overlap2 .and. overlap3) then
    !!    isoverlap = .true.
    !!else
    !!    isoverlap = .false.
    !!end if
        
    !@ END NEW VERSION #####################################
  
    !!!debug
    !!isoverlap=.true.
  
  end subroutine check_overlap_cubic_periodic



    !!!> Checks whether a segment with bounds i1,i2 (where i2 might be smaller
    !!!! than i1 due to periodic boundary conditions) overlaps with a segment with
    !!!! bounds j1,2 (where j1<=j2)
    !> Checks whether a segment with bounds i1,i2 (where i2 might be smaller
    !! than i1 due to periodic boundary conditions) overlaps with a segment with
    !! bounds j1,2 (where j2 might be smaller than j1)
    function check_whether_bounds_overlap_int(i1, i2, j1, j2) result(overlap)
      implicit none
      ! Calling arguments
      integer(kind=4),intent(in) :: i1, i2, j1, j2
      logical :: overlap
      ! Local variables
      integer :: periodic

      ! If the end is smaller than the start, we have a periodic wrap around
      periodic = 0
      if (i2<i1) then
          periodic = periodic + 1
      end if
      if (j2<j1) then
          periodic = periodic + 1
      end if

      ! Check whether there is an overlap
      select case(periodic)
      case(2)
          ! If both segments have a wrap around, they necessarily overlap
          overlap = .true.
      case(1)
          overlap = (i1<=j2 & !i2>=j1 due to periodic wrap around 
               .or. i2>=j1)   !i1<=j2 due to periodic wrap around
      case(0)
          overlap = (i2>=j1 .and. i1<=j2)
      case default
          stop 'wrong value of periodic'
      end select

    end function check_whether_bounds_overlap_int


    function check_whether_bounds_overlap_long(i1, i2, j1, j2) result(overlap)
      implicit none
      ! Calling arguments
      integer(kind=8),intent(in) :: i1, i2, j1, j2
      logical :: overlap
      ! Local variables
      integer :: periodic

      ! If the end is smaller than the start, we have a periodic wrap around
      periodic = 0
      if (i2<i1) then
          periodic = periodic + 1
      end if
      if (j2<j1) then
          periodic = periodic + 1
      end if

      ! Check whether there is an overlap
      select case(periodic)
      case(2)
          ! If both segments have a wrap around, they necessarily overlap
          overlap = .true.
      case(1)
          overlap = (i1<=j2 & !i2>=j1 due to periodic wrap around 
               .or. i2>=j1)   !i1<=j2 due to periodic wrap around
      case(0)
          overlap = (i2>=j1 .and. i1<=j2)
      case default
          stop 'wrong value of periodic'
      end select

    end function check_whether_bounds_overlap_long


    !> Checks whether a segment with bounds i1,i2 (where i2 might be smaller
    !! than i1 due to periodic boundary conditions) overlaps with a segment with
    !! bounds j1,2 (where j1<=j2). Is so, it gives the starting point, ending
    !! point and the extent of the (possibly two) overlaps.
    subroutine get_extent_of_overlap_int(i1, i2, j1, j2, n, ks, ke, nlen)
      use dictionaries, only: f_err_throw
      use yaml_strings, only: operator(//)
      implicit none
      ! Calling arguments
      integer(kind=4),intent(in) :: i1, i2, j1, j2
      integer(kind=4),intent(out) :: n !<number of overlaps
      integer(kind=4),dimension(2),intent(out) :: ks, ke, nlen
      ! Local variables
      integer :: ks1, ke1, ks2, ke2
      logical :: periodic, case1, case2, found_case

      ks(:) = 0
      ke(:) = 0
      nlen(:) = 0

      if (j2<j1) then
          call f_err_throw('j2<j1: '//&
               'i1='//i1//', i2='//i2//&
               ', j1='//j1//', j2='//j2,&
               err_name='BIGDFT_RUNTIME_ERROR')
!!$          call f_err_throw('j2<j1: '//&
!!$               &'i1='//trim(yaml_toa(i1,fmt='(i0)'))//&
!!$               &', i2='//trim(yaml_toa(i2,fmt='(i0)'))//&
!!$               &', j1='//trim(yaml_toa(j1,fmt='(i0)'))//&
!!$               &', j2='//trim(yaml_toa(j2,fmt='(i0)'))&
!!$               ,err_name='BIGDFT_RUNTIME_ERROR')
      end if

      ! Check whether there is an overlap
      if (check_whether_bounds_overlap(i1, i2, j1, j2)) then
          ! If the end is smaller than the start, we have a periodic wrap around
          periodic = (i2<i1)
          if (periodic) then
              found_case = .false.
              if (i2>=j1) then
                  ks1 = j1 !don't need to check i1 due to periodic wrap around
                  ke1 = min(i2,j2)
                  found_case = .true.
                  case1 = .true.
              else
                  ks1=huge(i2)
                  ke1=-huge(i2)
                  case1 = .false.
              end if
              if (i1<=j2) then
                  ks2 = max(i1,j1)
                  ke2 = j2 !don't need to check i2 due to periodic wrap around
                  found_case = .true.
                  case2 = .true.
              else
                  ks2=huge(i1)
                  ke2=-huge(i1)
                  case2 = .false.
              end if
              if (.not. found_case) then
                  call f_err_throw('Cannot determine overlap',err_name='BIGDFT_RUNTIME_ERROR')
              end if
              if (case1 .and. case2) then
                  ! There are two overlaps
                  n = 2
                  ks(1) = ks1
                  ke(1) = ke1
                  nlen(1) = ke(1) - ks(1) + 1
                  ks(2) = ks2
                  ke(2) = ke2
                  nlen(2) = ke(2) - ks(2) + 1
              else
                  n = 1
                  ks = min(ks1,ks2)
                  ke = max(ke1,ke2)
                  nlen = ke(1) - ks(1) + 1
              end if
          else
              n = 1
              ks(1) = max(i1,j1)
              ke(1) = min(i2,j2)
              nlen(1) = ke(1) - ks(1) + 1
          end if
          !write(*,'(a,7i8)') 'i1, i2, j1, j2, is, ie, n', i1, i2, j1, j2, is, ie, n
      else
          n = 0
          ks(1) = -1
          ke(1) = -1
          nlen(1) = 0
      end if

      if (nlen(1)<0) then
          call f_err_throw('nlen(1)<0: '//&
               &'i1='//  i1//&
               &', i2='//i2//&
               &', j1='//j1//&
               &', j2='//j2//&
               &', ks='//ks(1)//&
               &', ke='//ke(1)&
               ,err_name='BIGDFT_RUNTIME_ERROR')
      end if

      if (nlen(2)<0) then
         call f_err_throw('nlen(2)<0: i1='//  i1//&
              ', i2='//i2//&
              ', j1='//j1//&
              ', j2='//j2//&
              ', ks='//ks(2)//&
              ', ke='//ke(2)&
              ,err_name='BIGDFT_RUNTIME_ERROR')
      end if

    end subroutine get_extent_of_overlap_int


    subroutine get_extent_of_overlap_long(i1, i2, j1, j2, n, ks, ke, nlen)
      use dictionaries, only: f_err_throw
      use yaml_strings, only: operator(//)
      implicit none
      ! Calling arguments
      integer(kind=8),intent(in) :: i1, i2, j1, j2
      integer(kind=8),intent(out) :: n
      integer(kind=8),dimension(2),intent(out) :: ks, ke, nlen
      ! Local variables
      integer(kind=8) :: ks1, ke1, ks2, ke2
      logical :: periodic, case1, case2, found_case

      ks(:) = 0
      ke(:) = 0
      nlen(:) = 0

      if (j2<j1) then
          call f_err_throw('j2<j1: '//&
               &'i1='  //i1//&
               &', i2='//i2//&
               &', j1='//j1//&
               &', j2='//j2&
               ,err_name='BIGDFT_RUNTIME_ERROR')
      end if

      ! Check whether there is an overlap
      if (check_whether_bounds_overlap(i1, i2, j1, j2)) then
          ! If the end is smaller than the start, we have a periodic wrap around
          periodic = (i2<i1)
          if (periodic) then
              found_case = .false.
              if (i2>=j1) then
                  ks1 = j1 !don't need to check i1 due to periodic wrap around
                  ke1 = min(i2,j2)
                  found_case = .true.
                  case1 = .true.
              else
                  ks1=huge(i2)
                  ke1=-huge(i2)
                  case1 = .false.
              end if
              if (i1<=j2) then
                  ks2 = max(i1,j1)
                  ke2 = j2 !don't need to check i2 due to periodic wrap around
                  found_case = .true.
                  case2 = .true.
              else
                  ks2=huge(i1)
                  ke2=-huge(i1)
                  case2 = .false.
              end if
              if (.not. found_case) then
                  call f_err_throw('Cannot determine overlap',err_name='BIGDFT_RUNTIME_ERROR')
              end if
              if (case1 .and. case2) then
                  ! There are two overlaps
                  n = 2
                  ks(1) = ks1
                  ke(1) = ke1
                  nlen(1) = ke(1) - ks(1) + 1
                  ks(2) = ks2
                  ke(2) = ke2
                  nlen(2) = ke(2) - ks(2) + 1
              else
                  n = 1
                  ks = min(ks1,ks2)
                  ke = max(ke1,ke2)
                  nlen = ke(1) - ks(1) + 1
              end if
          else
              n = 1
              ks(1) = max(i1,j1)
              ke(1) = min(i2,j2)
              nlen(1) = ke(1) - ks(1) + 1
          end if
          !write(*,'(a,7i8)') 'i1, i2, j1, j2, is, ie, n', i1, i2, j1, j2, is, ie, n
      else
          n = 0
          ks(1) = -1
          ke(1) = -1
          nlen(1) = 0
      end if

      if (nlen(1)<0) then
          call f_err_throw('nlen(1)<0: i1='//  i1//&
               &', i2='//i2//&
               &', j1='//j1//&
               &', j2='//j2//&
               &', ks='//ks(1)//&
               &', ke='//ke(1)&
               ,err_name='BIGDFT_RUNTIME_ERROR')
      end if

      if (nlen(2)<0) then
          call f_err_throw('nlen(2)<0: i1='//  i1//&
               &', i2='//i2//&
               &', j1='//j1//&
               &', j2='//j2//&
               &', ks='//ks(2)//&
               &', ke='//ke(2)&
               ,err_name='BIGDFT_RUNTIME_ERROR')
      end if

    end subroutine get_extent_of_overlap_long

    subroutine check_overlap(Llr_i, Llr_j, Glr, overlap)
      implicit none

      ! Calling arguments
      type(locreg_descriptors),intent(in) :: Llr_i, Llr_j, Glr
      logical, intent(out) :: overlap

      ! Local variables
      integer :: onseg

      call check_overlap_cubic_periodic(Glr,Llr_i,Llr_j,overlap)
      if(overlap) then
         call check_overlap_from_descriptors_periodic(Llr_i%wfd%nseg_c, Llr_j%wfd%nseg_c,&
              Llr_i%wfd%keyglob, Llr_j%wfd%keyglob, overlap, onseg)
      end if

    end subroutine check_overlap

    ! check if Llrs overlap from there descriptors
    ! The periodicity is hidden in the fact that we are using the keyglobs
    ! which are correctly defined. 
    subroutine check_overlap_from_descriptors_periodic(nseg_i, nseg_j, keyg_i, keyg_j,  &
         isoverlap, onseg)
      implicit none
      ! Calling arguments
      integer :: nseg_i, nseg_j
      integer,dimension(2,nseg_i),intent(in) :: keyg_i
      integer,dimension(2,nseg_j),intent(in) :: keyg_j
      logical,intent(out) :: isoverlap
      integer, intent(out) :: onseg
      ! Local variables
      integer :: iseg, jseg, istart, jstart, kstartg
      integer :: iend, jend, kendg, nseg_k


      ! Initialize some counters
      iseg=1
      jseg=1
      nseg_k=0
      isoverlap = .false.
      onseg = 0  ! in case they don't overlap
      ! Check whether all segments of both localization regions have been processed.
      if ((iseg>=nseg_i .and. jseg>=nseg_j) .or. nseg_i==0 .or. nseg_j==0) return

      segment_loop: do

         ! Starting point already in global coordinates
         istart=keyg_i(1,iseg)
         jstart=keyg_j(1,jseg)

         ! Ending point already in global coordinates
         iend=keyg_i(2,iseg)
         jend=keyg_j(2,jseg)
         ! Determine starting and ending point of the common segment in global coordinates.
         kstartg=max(istart,jstart)
         kendg=min(iend,jend)

         ! Check whether this common segment has a non-zero length
         if(kendg-kstartg+1>0) then
            isoverlap = .true.
            nseg_k=nseg_k+1
         end if

         ! Check whether all segments of both localization regions have been processed.
         if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

         ! Increase the segment index
         if((iend<=jend .and. iseg<nseg_i) .or. jseg==nseg_j) then
            iseg=iseg+1
         else if(jseg<nseg_j) then
            jseg=jseg+1
         end if

      end do segment_loop

      if(isoverlap) then
         onseg = nseg_k
      end if

    end subroutine check_overlap_from_descriptors_periodic


end module locregs
