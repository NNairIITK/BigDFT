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


!!$  !> De-Allocate convolutions_bounds type, depending of the geocode and the hybrid_on
!!$  subroutine deallocate_bounds(geocode,hybrid_on,bounds)
!!$    use module_base
!!$    implicit none
!!$    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
!!$    logical, intent(in) :: hybrid_on 
!!$    type(convolutions_bounds) :: bounds
!!$
!!$    !if ((geocode == 'P' .and. hybrid_on) .or. geocode == 'F') then
!!$    !   ! Just test the first one...
!!$    !   if (associated(bounds%kb%ibyz_f)) then
!!$    call f_free_ptr(bounds%kb%ibyz_f)
!!$    call f_free_ptr(bounds%kb%ibxz_f)
!!$    call f_free_ptr(bounds%kb%ibxy_f)
!!$
!!$    call f_free_ptr(bounds%sb%ibxy_ff)
!!$    call f_free_ptr(bounds%sb%ibzzx_f)
!!$    call f_free_ptr(bounds%sb%ibyyzz_f)
!!$
!!$    call f_free_ptr(bounds%gb%ibyz_ff)
!!$
!!$    call f_free_ptr(bounds%gb%ibzxx_f)
!!$    call f_free_ptr(bounds%gb%ibxxyy_f)
!!$    !   end if
!!$    !end if
!!$
!!$    !the arrays which are needed only for free BC
!!$    !if (geocode == 'F') then
!!$    ! Just test the first one...
!!$    !   if (associated(bounds%kb%ibyz_c)) then
!!$    call f_free_ptr(bounds%kb%ibyz_c)
!!$    call f_free_ptr(bounds%kb%ibxz_c)
!!$    call f_free_ptr(bounds%kb%ibxy_c)
!!$
!!$
!!$    call f_free_ptr(bounds%sb%ibzzx_c)
!!$    call f_free_ptr(bounds%sb%ibyyzz_c)
!!$
!!$    call f_free_ptr(bounds%gb%ibzxx_c)
!!$    call f_free_ptr(bounds%gb%ibxxyy_c)
!!$
!!$    call f_free_ptr(bounds%ibyyzz_r)
!!$
!!$    !  end if
!!$    !end if
!!$
!!$  END SUBROUTINE deallocate_bounds


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
!!$    character(len=*),intent(in):: subname
!!$    ! Calling arguments
!!$    character(len=1),intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
!!$
!!$    ! Local variables
!!$    integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3

    call copy_kinetic_bounds(boundsin%kb, boundsout%kb)
    call copy_shrink_bounds(boundsin%sb, boundsout%sb)
    call copy_grow_bounds(boundsin%gb, boundsout%gb)
    boundsout%ibyyzz_r = f_malloc_ptr(src_ptr=boundsin%ibyyzz_r,id='boundsout%ibyyzz_r')

!!$
!!$    if(geocode == 'F') then
!!$       if(associated(boundsout%ibyyzz_r)) then
!!$          call f_free_ptr(boundsout%ibyyzz_r)
!!$       end if
!!$
!!$       if(associated(boundsin%ibyyzz_r)) then
!!$          iis1=lbound(boundsin%ibyyzz_r,1)
!!$          iie1=ubound(boundsin%ibyyzz_r,1)
!!$          iis2=lbound(boundsin%ibyyzz_r,2)
!!$          iie2=ubound(boundsin%ibyyzz_r,2)
!!$          iis3=lbound(boundsin%ibyyzz_r,3)
!!$          iie3=ubound(boundsin%ibyyzz_r,3)
!!$          boundsout%ibyyzz_r = f_malloc_ptr((/ iis1.to.iie1,iis2.to.iie2,iis3.to.iie3 /),id='boundsout%ibyyzz_r')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   boundsout%ibyyzz_r(i1,i2,i3) = boundsin%ibyyzz_r(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if
!!$    end if
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

!!$    character(len=1),intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode 
!!$    character(len=*),intent(in):: subname
!!$
!!$    ! Local variables
!!$    integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3
!!$
!!$    if(geocode == 'F') then
!!$       if(associated(kbout%ibyz_c)) then
!!$          call f_free_ptr(kbout%ibyz_c)
!!$       end if
!!$       if(associated(kbin%ibyz_c)) then
!!$          iis1=lbound(kbin%ibyz_c,1)
!!$          iie1=ubound(kbin%ibyz_c,1)
!!$          iis2=lbound(kbin%ibyz_c,2)
!!$          iie2=ubound(kbin%ibyz_c,2)
!!$          iis3=lbound(kbin%ibyz_c,3)
!!$          iie3=ubound(kbin%ibyz_c,3)
!!$          kbout%ibyz_c = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 , iis3.to.iie3 /),id='kbout%ibyz_c')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   kbout%ibyz_c(i1,i2,i3) = kbin%ibyz_c(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if

!!$       if(associated(kbout%ibxz_c)) then
!!$          call f_free_ptr(kbout%ibxz_c)
!!$       end if
!!$       if(associated(kbin%ibxz_c)) then
!!$          iis1=lbound(kbin%ibxz_c,1)
!!$          iie1=ubound(kbin%ibxz_c,1)
!!$          iis2=lbound(kbin%ibxz_c,2)
!!$          iie2=ubound(kbin%ibxz_c,2)
!!$          iis3=lbound(kbin%ibxz_c,3)
!!$          iie3=ubound(kbin%ibxz_c,3)
!!$          kbout%ibxz_c = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 , iis3.to.iie3 /),id='kbout%ibxz_c')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   kbout%ibxz_c(i1,i2,i3) = kbin%ibxz_c(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if

!!$       if(associated(kbout%ibxy_c)) then
!!$          call f_free_ptr(kbout%ibxy_c)
!!$       end if
!!$       if(associated(kbin%ibxy_c)) then
!!$          iis1=lbound(kbin%ibxy_c,1)
!!$          iie1=ubound(kbin%ibxy_c,1)
!!$          iis2=lbound(kbin%ibxy_c,2)
!!$          iie2=ubound(kbin%ibxy_c,2)
!!$          iis3=lbound(kbin%ibxy_c,3)
!!$          iie3=ubound(kbin%ibxy_c,3)
!!$          kbout%ibxy_c = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 , iis3.to.iie3 /),id='kbout%ibxy_c')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   kbout%ibxy_c(i1,i2,i3) = kbin%ibxy_c(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if
!!$    end if

!!$    if(associated(kbout%ibyz_f)) then
!!$       call f_free_ptr(kbout%ibyz_f)
!!$    end if
!!$    if(associated(kbin%ibyz_f)) then
!!$       iis1=lbound(kbin%ibyz_f,1)
!!$       iie1=ubound(kbin%ibyz_f,1)
!!$       iis2=lbound(kbin%ibyz_f,2)
!!$       iie2=ubound(kbin%ibyz_f,2)
!!$       iis3=lbound(kbin%ibyz_f,3)
!!$       iie3=ubound(kbin%ibyz_f,3)
!!$       kbout%ibyz_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='kbout%ibyz_f')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                kbout%ibyz_f(i1,i2,i3) = kbin%ibyz_f(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

!!$    if(associated(kbout%ibxz_f)) then
!!$       call f_free_ptr(kbout%ibxz_f)
!!$    end if
!!$    if(associated(kbin%ibxz_f)) then
!!$       iis1=lbound(kbin%ibxz_f,1)
!!$       iie1=ubound(kbin%ibxz_f,1)
!!$       iis2=lbound(kbin%ibxz_f,2)
!!$       iie2=ubound(kbin%ibxz_f,2)
!!$       iis3=lbound(kbin%ibxz_f,3)
!!$       iie3=ubound(kbin%ibxz_f,3)
!!$       kbout%ibxz_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='kbout%ibxz_f')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                kbout%ibxz_f(i1,i2,i3) = kbin%ibxz_f(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

!!$    if(associated(kbout%ibxy_f)) then
!!$       call f_free_ptr(kbout%ibxy_f)
!!$    end if
!!$    if(associated(kbin%ibxy_f)) then
!!$       iis1=lbound(kbin%ibxy_f,1)
!!$       iie1=ubound(kbin%ibxy_f,1)
!!$       iis2=lbound(kbin%ibxy_f,2)
!!$       iie2=ubound(kbin%ibxy_f,2)
!!$       iis3=lbound(kbin%ibxy_f,3)
!!$       iie3=ubound(kbin%ibxy_f,3)
!!$       kbout%ibxy_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='kbout%ibxy_f')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                kbout%ibxy_f(i1,i2,i3) = kbin%ibxy_f(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

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

!!$    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
!!$    character(len=*),intent(in):: subname
!!$    ! Local variables
!!$    integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3
!!$
!!$    if(geocode == 'F') then
!!$       if(associated(sbout%ibzzx_c)) then
!!$          call f_free_ptr(sbout%ibzzx_c)
!!$       end if
!!$       if(associated(sbin%ibzzx_c)) then
!!$          iis1=lbound(sbin%ibzzx_c,1)
!!$          iie1=ubound(sbin%ibzzx_c,1)
!!$          iis2=lbound(sbin%ibzzx_c,2)
!!$          iie2=ubound(sbin%ibzzx_c,2)
!!$          iis3=lbound(sbin%ibzzx_c,3)
!!$          iie3=ubound(sbin%ibzzx_c,3)
!!$          sbout%ibzzx_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibzzx_c')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   sbout%ibzzx_c(i1,i2,i3) = sbin%ibzzx_c(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if

!!$       if(associated(sbout%ibyyzz_c)) then
!!$          call f_free_ptr(sbout%ibyyzz_c)
!!$       end if
!!$       if(associated(sbin%ibyyzz_c)) then
!!$          iis1=lbound(sbin%ibyyzz_c,1)
!!$          iie1=ubound(sbin%ibyyzz_c,1)
!!$          iis2=lbound(sbin%ibyyzz_c,2)
!!$          iie2=ubound(sbin%ibyyzz_c,2)
!!$          iis3=lbound(sbin%ibyyzz_c,3)
!!$          iie3=ubound(sbin%ibyyzz_c,3)
!!$          sbout%ibyyzz_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibyyzz_c')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   sbout%ibyyzz_c(i1,i2,i3) = sbin%ibyyzz_c(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if
!!$    end if

!!$    if(associated(sbout%ibxy_ff)) then
!!$       call f_free_ptr(sbout%ibxy_ff)
!!$    end if
!!$    if(associated(sbin%ibxy_ff)) then
!!$       iis1=lbound(sbin%ibxy_ff,1)
!!$       iie1=ubound(sbin%ibxy_ff,1)
!!$       iis2=lbound(sbin%ibxy_ff,2)
!!$       iie2=ubound(sbin%ibxy_ff,2)
!!$       iis3=lbound(sbin%ibxy_ff,3)
!!$       iie3=ubound(sbin%ibxy_ff,3)
!!$       sbout%ibxy_ff = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibxy_ff')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                sbout%ibxy_ff(i1,i2,i3) = sbin%ibxy_ff(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if


!!$    if(associated(sbout%ibzzx_f)) then
!!$       call f_free_ptr(sbout%ibzzx_f)
!!$    end if
!!$    if(associated(sbin%ibzzx_f)) then
!!$       iis1=lbound(sbin%ibzzx_f,1)
!!$       iie1=ubound(sbin%ibzzx_f,1)
!!$       iis2=lbound(sbin%ibzzx_f,2)
!!$       iie2=ubound(sbin%ibzzx_f,2)
!!$       iis3=lbound(sbin%ibzzx_f,3)
!!$       iie3=ubound(sbin%ibzzx_f,3)
!!$       sbout%ibzzx_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibzzx_f')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                sbout%ibzzx_f(i1,i2,i3) = sbin%ibzzx_f(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

!!$    if(associated(sbout%ibyyzz_f)) then
!!$       call f_free_ptr(sbout%ibyyzz_f)
!!$    end if
!!$    if(associated(sbin%ibyyzz_f)) then
!!$       iis1=lbound(sbin%ibyyzz_f,1)
!!$       iie1=ubound(sbin%ibyyzz_f,1)
!!$       iis2=lbound(sbin%ibyyzz_f,2)
!!$       iie2=ubound(sbin%ibyyzz_f,2)
!!$       iis3=lbound(sbin%ibyyzz_f,3)
!!$       iie3=ubound(sbin%ibyyzz_f,3)
!!$       sbout%ibyyzz_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibyyzz_f')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                sbout%ibyyzz_f(i1,i2,i3) = sbin%ibyyzz_f(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

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

!!$    character(len=1),intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
!!$    character(len=*),intent(in):: subname
!!$
!!$    ! Local variables
!!$    integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3
!!$
!!$
!!$
!!$    if(geocode == 'F')then
!!$       if(associated(gbout%ibzxx_c)) then
!!$          call f_free_ptr(gbout%ibzxx_c)
!!$       end if
!!$       if(associated(gbin%ibzxx_c)) then
!!$          iis1=lbound(gbin%ibzxx_c,1)
!!$          iie1=ubound(gbin%ibzxx_c,1)
!!$          iis2=lbound(gbin%ibzxx_c,2)
!!$          iie2=ubound(gbin%ibzxx_c,2)
!!$          iis3=lbound(gbin%ibzxx_c,3)
!!$          iie3=ubound(gbin%ibzxx_c,3)
!!$          gbout%ibzxx_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibzxx_c')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   gbout%ibzxx_c(i1,i2,i3) = gbin%ibzxx_c(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if
!!$
!!$
!!$       if(associated(gbout%ibxxyy_c)) then
!!$          call f_free_ptr(gbout%ibxxyy_c)
!!$       end if
!!$       if(associated(gbin%ibxxyy_c)) then
!!$          iis1=lbound(gbin%ibxxyy_c,1)
!!$          iie1=ubound(gbin%ibxxyy_c,1)
!!$          iis2=lbound(gbin%ibxxyy_c,2)
!!$          iie2=ubound(gbin%ibxxyy_c,2)
!!$          iis3=lbound(gbin%ibxxyy_c,3)
!!$          iie3=ubound(gbin%ibxxyy_c,3)
!!$          gbout%ibxxyy_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibxxyy_c')
!!$          do i3=iis3,iie3
!!$             do i2=iis2,iie2
!!$                do i1=iis1,iie1
!!$                   gbout%ibxxyy_c(i1,i2,i3) = gbin%ibxxyy_c(i1,i2,i3)
!!$                end do
!!$             end do
!!$          end do
!!$       end if
!!$    end if
!!$
!!$    if(associated(gbout%ibyz_ff)) then
!!$       call f_free_ptr(gbout%ibyz_ff)
!!$    end if
!!$    if(associated(gbin%ibyz_ff)) then
!!$       iis1=lbound(gbin%ibyz_ff,1)
!!$       iie1=ubound(gbin%ibyz_ff,1)
!!$       iis2=lbound(gbin%ibyz_ff,2)
!!$       iie2=ubound(gbin%ibyz_ff,2)
!!$       iis3=lbound(gbin%ibyz_ff,3)
!!$       iie3=ubound(gbin%ibyz_ff,3)
!!$       gbout%ibyz_ff = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibyz_ff')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                gbout%ibyz_ff(i1,i2,i3) = gbin%ibyz_ff(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

!!$    if(associated(gbout%ibzxx_f)) then
!!$       call f_free_ptr(gbout%ibzxx_f)
!!$    end if
!!$    if(associated(gbin%ibzxx_f)) then
!!$       iis1=lbound(gbin%ibzxx_f,1)
!!$       iie1=ubound(gbin%ibzxx_f,1)
!!$       iis2=lbound(gbin%ibzxx_f,2)
!!$       iie2=ubound(gbin%ibzxx_f,2)
!!$       iis3=lbound(gbin%ibzxx_f,3)
!!$       iie3=ubound(gbin%ibzxx_f,3)
!!$       gbout%ibzxx_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibzxx_f')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                gbout%ibzxx_f(i1,i2,i3) = gbin%ibzxx_f(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

!!$    if(associated(gbout%ibxxyy_f)) then
!!$       call f_free_ptr(gbout%ibxxyy_f)
!!$    end if
!!$    if(associated(gbin%ibxxyy_f)) then
!!$       iis1=lbound(gbin%ibxxyy_f,1)
!!$       iie1=ubound(gbin%ibxxyy_f,1)
!!$       iis2=lbound(gbin%ibxxyy_f,2)
!!$       iie2=ubound(gbin%ibxxyy_f,2)
!!$       iis3=lbound(gbin%ibxxyy_f,3)
!!$       iie3=ubound(gbin%ibxxyy_f,3)
!!$       gbout%ibxxyy_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibxxyy_f')
!!$       do i3=iis3,iie3
!!$          do i2=iis2,iie2
!!$             do i1=iis1,iie1
!!$                gbout%ibxxyy_f(i1,i2,i3) = gbin%ibxxyy_f(i1,i2,i3)
!!$             end do
!!$          end do
!!$       end do
!!$    end if


  end subroutine copy_grow_bounds


end module locregs
