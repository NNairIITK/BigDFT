!> @file
!! Datatypes and associated methods relative to the localization regions
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!> Datatypes for localization regions descriptors
module locregs
  use module_base
  implicit none

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
     character(len=1) :: geocode                !< @copydoc poisson_solver::doc::geocode
     logical :: hybrid_on                       !< Interesting for global, periodic, localisation regions
     integer :: ns1,ns2,ns3                     !< Starting point of the localisation region in global coordinates
     integer :: nsi1,nsi2,nsi3                  !< Starting point of locreg for interpolating grid
     integer :: Localnorb                       !< Number of orbitals contained in locreg
     integer, dimension(3) :: outofzone         !< Vector of points outside of the zone outside Glr for periodic systems
     real(gp), dimension(3) :: locregCenter !< Center of the locreg 
     real(gp) :: locrad                     !< Cutoff radius of the localization region
     real(gp) :: locrad_kernel              !< Cutoff radius of the localization region (kernel)
     type(grid_dimensions) :: d
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
  end type locreg_descriptors

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

  !initializations
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
    
!!$    allocate(wfd%keyglob(2,max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
!!$    call memocc(i_stat,wfd%keyglob,'keyglob',subname)
!!$    allocate(wfd%keygloc(2,max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
!!$    call memocc(i_stat,wfd%keygloc,'keygloc',subname)
!!$    allocate(wfd%keyvloc(max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
!!$    call memocc(i_stat,wfd%keyvloc,'keyvloc',subname)
!!$    allocate(wfd%keyvglob(max(1,wfd%nseg_c+wfd%nseg_f+ndebug)),stat=i_stat)
!!$    call memocc(i_stat,wfd%keyvglob,'keyvglob',subname)

  END SUBROUTINE allocate_wfd

  !> De-Allocate wavefunctions_descriptors
  subroutine deallocate_wfd(wfd)
    use module_base
    implicit none
    type(wavefunctions_descriptors) :: wfd

    !in case the two objects points to the same target
    !pay attention that in this case odd behaviour of f_mallo may occur as the
    !pointers have not been associated by the f_associate routine (to be implemented to date)
    if (associated(wfd%keyglob, target = wfd%keygloc)) then
       !assuming that globals has been created afterwards
       call f_free_ptr(wfd%keyglob)
       nullify(wfd%keygloc)
!!$       i_all=-product(shape(wfd%keyglob))*kind(wfd%keyglob)
!!$       deallocate(wfd%keyglob,stat=i_stat)
!!$       call memocc(i_stat,i_all,'wfd%keyglob',subname)
!!$       nullify(wfd%keyglob)
    else
       call f_free_ptr(wfd%keygloc)
       call f_free_ptr(wfd%keyglob)
!!$       if(associated(wfd%keyglob)) then
!!$          i_all=-product(shape(wfd%keyglob))*kind(wfd%keyglob)
!!$          deallocate(wfd%keyglob,stat=i_stat)
!!$          call memocc(i_stat,i_all,'wfd%keyglob',subname)
!!$          nullify(wfd%keyglob)
!!$       end if
!!$       if(associated(wfd%keygloc)) then 
!!$          i_all=-product(shape(wfd%keygloc))*kind(wfd%keygloc)
!!$          deallocate(wfd%keygloc,stat=i_stat)
!!$          call memocc(i_stat,i_all,'wfd%keygloc',subname)
!!$          nullify(wfd%keygloc)
!!$       end if
    end if
    if (associated(wfd%keyvloc, target= wfd%keyvglob)) then
!!$       i_all=-product(shape(wfd%keyvloc))*kind(wfd%keyvloc)
!!$       deallocate(wfd%keyvloc,stat=i_stat)
!!$       call memocc(i_stat,i_all,'wfd%keyvloc',subname)
!!$       nullify(wfd%keyvloc)
       call f_free_ptr(wfd%keyvglob)
       nullify(wfd%keyvloc)
    else
       call f_free_ptr(wfd%keyvloc)
       call f_free_ptr(wfd%keyvglob)
!!$       if (associated(wfd%keyvloc)) then
!!$          i_all=-product(shape(wfd%keyvloc))*kind(wfd%keyvloc)
!!$          deallocate(wfd%keyvloc,stat=i_stat)
!!$          call memocc(i_stat,i_all,'wfd%keyvloc',subname)
!!$          nullify(wfd%keyvloc)
!!$       end if
!!$       if (associated(wfd%keyvglob)) then
!!$          i_all=-product(shape(wfd%keyvglob))*kind(wfd%keyvglob)
!!$          deallocate(wfd%keyvglob,stat=i_stat)
!!$          call memocc(i_stat,i_all,'wfd%keyvglob',subname)
!!$          nullify(wfd%keyvglob)
!!$       end if
    end if
  END SUBROUTINE deallocate_wfd

  !>desctructors
  subroutine deallocate_locreg_descriptors(lr)
    implicit none
    ! Calling arguments
    type(locreg_descriptors),intent(inout):: lr

    call deallocate_wfd(lr%wfd)
    call deallocate_convolutions_bounds(lr%bounds,'deallocate_locreg_descriptors')

  end subroutine deallocate_locreg_descriptors

  !> De-Allocate convolutions_bounds type, depending of the geocode and the hybrid_on
  subroutine deallocate_bounds(geocode,hybrid_on,bounds,subname)
    use module_base
    implicit none
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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

  !methods for copying the structures, can be needed to avoid recalculating them
  !should be better by defining a f_malloc inheriting the shapes and the structure from other array
  !of the type dest=f_malloc(src=source,id='dest')
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
    glrout%locregCenter=glrin%locregCenter
    glrout%outofzone= glrin%outofzone

    call copy_grid_dimensions(glrin%d, glrout%d)
    call copy_wavefunctions_descriptors(glrin%wfd, glrout%wfd)
    !copy bound when needed
    if(glrin%geocode == 'F' .or. (glrin%geocode == 'P' .and. glrin%hybrid_on)) then
       call copy_convolutions_bounds(glrin%geocode, glrin%bounds, glrout%bounds,&
            'copy_locreg_descriptors') !to be removed when bounds are allocated properly
    else
       call nullify_convolutions_bounds(glrout%bounds)
    end if

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
!    integer:: i1, i2, iis1, iie1, iis2, iie2, istat, iall

    !nullify all pointers first
    call nullify_wfd(wfdout)

    wfdout%nvctr_c = wfdin%nvctr_c
    wfdout%nvctr_f = wfdin%nvctr_f
    wfdout%nseg_c = wfdin%nseg_c
    wfdout%nseg_f = wfdin%nseg_f

    if (associated(wfdin%keygloc)) wfdout%keygloc=f_malloc_ptr(src=wfdin%keygloc,id='wfdout%keygloc')
    if (associated(wfdin%keyglob)) wfdout%keyglob=f_malloc_ptr(src=wfdin%keyglob,id='wfdout%keyglob')
    if (associated(wfdin%keyvloc)) wfdout%keyvloc=f_malloc_ptr(src=wfdin%keyvloc,id='wfdout%keyvloc')
    if (associated(wfdin%keyvglob))wfdout%keyvglob=f_malloc_ptr(src=wfdin%keyvglob,id='wfdout%keyvglob')

!!$    if(associated(wfdout%keygloc)) then
!!$       iall=-product(shape(wfdout%keygloc))*kind(wfdout%keygloc)
!!$       deallocate(wfdout%keygloc, stat=istat)
!!$       call memocc(istat, iall, 'wfdout%keygloc', subname)
!!$    end if
!!$    if(associated(wfdin%keygloc)) then
!!$       iis1=lbound(wfdin%keygloc,1)
!!$       iie1=ubound(wfdin%keygloc,1)
!!$       iis2=lbound(wfdin%keygloc,2)
!!$       iie2=ubound(wfdin%keygloc,2)
!!$
!!$       allocate(wfdout%keygloc(iis1:iie1,iis2:iie2), stat=istat)
!!$       call memocc(istat, wfdout%keygloc, 'wfdout%keygloc', subname)
!!$       do i2=iis2,iie2
!!$          do i1=iis1,iie1
!!$             wfdout%keygloc(i1,i2) = wfdin%keygloc(i1,i2)
!!$          end do
!!$       end do
!!$    end if
!!$
!!$    if(associated(wfdout%keyglob)) then
!!$       iall=-product(shape(wfdout%keyglob))*kind(wfdout%keygloc)
!!$       deallocate(wfdout%keyglob, stat=istat)
!!$       call memocc(istat, iall, 'wfdout%keyglob', subname)
!!$    end if
!!$    if(associated(wfdin%keyglob)) then
!!$       iis1=lbound(wfdin%keyglob,1)
!!$       iie1=ubound(wfdin%keyglob,1)
!!$       iis2=lbound(wfdin%keyglob,2)
!!$       iie2=ubound(wfdin%keyglob,2)
!!$       allocate(wfdout%keyglob(iis1:iie1,iis2:iie2), stat=istat)
!!$       call memocc(istat, wfdout%keyglob, 'wfdout%keyglob', subname)
!!$       do i2=iis2,iie2
!!$          do i1=iis1,iie1
!!$             wfdout%keyglob(i1,i2) = wfdin%keyglob(i1,i2)
!!$          end do
!!$       end do
!!$    end if
!!$
!!$    if(associated(wfdout%keyvloc)) then
!!$       iall=-product(shape(wfdout%keyvloc))*kind(wfdout%keyvloc)
!!$       deallocate(wfdout%keyvloc, stat=istat)
!!$       call memocc(istat, iall, 'wfdout%keyvloc', subname)
!!$    end if
!!$    if(associated(wfdin%keyvloc)) then
!!$       iis1=lbound(wfdin%keyvloc,1)
!!$       iie1=ubound(wfdin%keyvloc,1)
!!$       allocate(wfdout%keyvloc(iis1:iie1), stat=istat)
!!$       call memocc(istat, wfdout%keyvloc, 'wfdout%keyvloc', subname)
!!$       do i1=iis1,iie1
!!$          wfdout%keyvloc(i1) = wfdin%keyvloc(i1)
!!$       end do
!!$    end if
!!$
!!$    if(associated(wfdout%keyvglob)) then
!!$       iall=-product(shape(wfdout%keyvglob))*kind(wfdout%keyvglob)
!!$       deallocate(wfdout%keyvglob, stat=istat)
!!$       call memocc(istat, iall, 'wfdout%keyvglob', subname)
!!$    end if
!!$    if(associated(wfdin%keyvglob)) then
!!$       iis1=lbound(wfdin%keyvglob,1)
!!$       iie1=ubound(wfdin%keyvglob,1)
!!$       allocate(wfdout%keyvglob(iis1:iie1), stat=istat)
!!$       call memocc(istat, wfdout%keyvglob, 'wfdout%keyvglob', subname)
!!$       do i1=iis1,iie1
!!$          wfdout%keyvglob(i1) = wfdin%keyvglob(i1)
!!$       end do
!!$    end if

  end subroutine copy_wavefunctions_descriptors


end module locregs
