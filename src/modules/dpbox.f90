!> @file
!!  Define the density and potential grid datastructure 
!! @author
!!    Copyright (C) 2015-2015 BigDFT group (TD)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Module which contains the data structure associated to the dnesity and potential grid
module module_dpbox

  use module_base, only: gp,mpi_environment,mpi_environment_null,f_err_throw

  implicit none

  integer, parameter, public :: DPB_POT_ION = 0 !< Use n3pi for the iterations over z planes
  integer, parameter, public :: DPB_RHO     = 1 !< Use n3d  for the iterations over z planes
  integer, parameter, public :: DPB_POT     = 2 !< Use n3p  for the iterations over z planes

  private

  !> Structure to store the density / potential distribution among processors.
  type, public :: denspot_distribution
     integer :: n3d                  !< Number of z planes distributed for density
     integer :: n3p                  !< Number of z planes distirbuted for potential (except pot_ion)
     integer :: n3pi                 !< Number of distributed planes in z dimension for pot_ion AND to calculate charges
                                     !! BECAUSE n3d has an overlap!
                                     !! ONLY FOR POT_ION
     integer :: i3xcsh               !< GGA XC shift between density and potential
     integer :: i3s                  !< Index of the first z plane (offset) for the mpi process i.e. from i3s:i3s+n3pi-1 
     integer :: nrhodim              !< nspin !
     !> Integer which controls the presence of a density after the potential array
     !! if different than zero, at the address ndimpot*nspin+i3rho_add starts the spin up component of the density
     !! the spin down component can be found at the ndimpot*nspin+i3rho_add+ndimpot, contiguously
     !! the same holds for non-collinear calculations
     integer :: i3rho_add             !< dpbox%ndims(1)*dpbox%ndims(2)*dpbox%i3xcsh+1
     integer :: ndimpot               !< n1i*n2i*n3p = dpbox%ndims(1)dpbox%ndims(2)*dpbox%n3p
     integer :: ndimgrid              !< n1i*n2i*n3i = dpbox%ndims(1)*dpbox%ndims(2)*dpbox%ndims(3)
     integer :: ndimrhopot            !< dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3d*dpbox%nrhodim
     integer, dimension(3) :: ndims   !< Box containing the grid dimensions in ISF basis in x,y and z direction (n1i,n2i,n3i)
     real(gp), dimension(3) :: hgrids !< Grid spacings of the box (half of wavelet ones)
     character(len=1) :: geocode      !< @copydoc poisson_solver::doc::geocode
     integer, dimension(:,:), pointer :: nscatterarr !< dim(nproc,4) for each proc (n3d,n3p,i3s+i3xcsh-1,i3xcsh) see @link dpbox_repartition @endlink
     integer, dimension(:,:), pointer :: ngatherarr  !< dim(nproc,3) (dpbox%ndimpot,n1i*n2i*nscatteradd(:,3),n1i*n2i*n3d) see @link dpbox_repartition @endlink
     type(mpi_environment) :: mpi_env
  end type denspot_distribution


  !> Define an iterator over the points of the grid which should be also inside a given box (for instance centered on an atom)
  type, public :: dpbox_iterator
    integer :: ix,iy,iz                  !< Indices of the three-dimensional arrays in distributed PSolver data scheme
    integer :: ind                       !< One dimensional index (for pot_ion)
    integer, dimension(3)  :: ibox       !< 3D indices in absolute coordinates in the given box specified by boxat
    real(gp) :: x,y,z                    !< 3D absolute coordinates inside the given box
    !private
    integer, dimension(2,3) :: nbox      !< Specify a sub-box to iterate over the points (ex. around atoms)
    integer :: n1i,n2i,n3i               !< 3D dimension of the whole grid
    integer :: i3s                       !< Index of the first z plane for the mpi process i.e. from i3s:i3s+n3pi-1 
    integer :: n3_iter                   !< Indicate Z dimension when iter depending on should be n3pi,n3p,n3d
    integer :: nbl1,nbr1                 !< Size of left and right buffers in x direction
    integer :: nbl2,nbr2                 !< Size of left and right buffers in y direction
    integer :: nbl3,nbr3                 !< Size of left and right buffers in z direction
    logical :: perx,pery,perz            !< Conditions for periodicity in the three directions
    logical :: whole                     !< Iterate over the whole box or not
    type(denspot_distribution), pointer :: dpbox_ptr !< Private pointer to the original dpbox on which we are iterating
  end type dpbox_iterator


  !Public routines
  public :: dpbox_null,deallocate_denspot_distribution,dpbox_free
  public :: dpbox_iter,dpbox_iter_next
  public :: dpbox_iterator_null, nullify_dpbox_iterator
  

contains


  !> Nullify the denspot_distribution structure
  function dpbox_null() result(dpbox)
    implicit none
    type(denspot_distribution) :: dpbox
    call nullify_denspot_distribution(dpbox)
  end function dpbox_null


  !> Nullify the denspot_distribution structure
  subroutine nullify_denspot_distribution(dpbox)
    implicit none
    type(denspot_distribution),intent(out) :: dpbox
    dpbox%n3d=0
    dpbox%n3p=0
    dpbox%i3xcsh=0
    dpbox%i3s=0
    dpbox%nrhodim=0
    dpbox%i3rho_add=0
    dpbox%ndimpot=0
    dpbox%ndimgrid=0
    dpbox%ndimrhopot=0
    dpbox%ndims=(/0,0,0/)
    dpbox%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    dpbox%geocode = "F"
    nullify(dpbox%nscatterarr)
    nullify(dpbox%ngatherarr)
    dpbox%mpi_env=mpi_environment_null()
  end subroutine nullify_denspot_distribution


  !> Deallocate the denspot_distribution structure
  subroutine deallocate_denspot_distribution(dpbox)
    use module_base, only: f_free_ptr
    implicit none
    type(denspot_distribution),intent(inout)::dpbox
    
    if(associated(dpbox%nscatterarr)) then
      call f_free_ptr(dpbox%nscatterarr)
    end if
    if(associated(dpbox%ngatherarr)) then
      call f_free_ptr(dpbox%ngatherarr)
    end if

  end subroutine deallocate_denspot_distribution


  !> Free the denspot_distribution structure
  subroutine dpbox_free(dpbox)
    use module_base, only: f_free_ptr,mpi_environment_free,bigdft_mpi
    implicit none
    type(denspot_distribution), intent(inout) :: dpbox

    if (associated(dpbox%nscatterarr)) then
       call f_free_ptr(dpbox%nscatterarr)
    end if

    if (associated(dpbox%ngatherarr)) then
       call f_free_ptr(dpbox%ngatherarr)
    end if
    
    if (dpbox%mpi_env%mpi_comm /= bigdft_mpi%mpi_comm) then
       call mpi_environment_free(dpbox%mpi_env)
    end if

    dpbox=dpbox_null()

  END SUBROUTINE dpbox_free


  !> Function nullify an iterator over dpbox
  pure function dpbox_iterator_null() result (boxit)
    implicit none
    type(dpbox_iterator) :: boxit
    call  nullify_dpbox_iterator(boxit)
  end function dpbox_iterator_null


  !> Nullify the iterator dpbox type
  pure subroutine nullify_dpbox_iterator(boxit)
    implicit none
    type(dpbox_iterator), intent(out) :: boxit
    boxit%ix = -1
    boxit%iy = -1
    boxit%iz = -1
    boxit%ind = -1
    boxit%ibox(:) = -1
    boxit%nbox(:,:) = -1
    boxit%n1i = -1
    boxit%n2i = -1
    boxit%n3i = -1
    boxit%i3s = -1
    boxit%nbl1 = -1
    boxit%nbr1 = -1
    boxit%nbl2 = -1
    boxit%nbr2 = -1
    boxit%nbl3 = -1
    boxit%nbr3 = -1
    boxit%x = 0.0_gp
    boxit%y = 0.0_gp
    boxit%z = 0.0_gp
    nullify(boxit%dpbox_ptr)
  end subroutine nullify_dpbox_iterator


  !> Create an iterator dpbox to iterate over points of the (potential) grid 
  function dpbox_iter(dpbox,idpbox,nbox) result(boxit)
    implicit none
    !Arguments
    type(denspot_distribution), intent(in), target :: dpbox !< Density-potential descriptors for the box
    !> Indicate if we iterate over pot_ion (n3pi), over rho (n3d) or over rhov (n3p)
    integer, intent(in) :: idpbox 
    !> Box of start and end points which have to be considered
    integer, dimension(2,3), intent(in), optional :: nbox
    type(dpbox_iterator) :: boxit

    call nullify_dpbox_iterator(boxit)

    !Associate the original objects
    boxit%dpbox_ptr => dpbox
    !Distributed dimension over dpbox%ndims(3) in parallel
    boxit%n1i = boxit%dpbox_ptr%ndims(1)
    boxit%n2i = boxit%dpbox_ptr%ndims(2)
    boxit%n3i = boxit%dpbox_ptr%ndims(3)
    !This is correct for a potential not a density
    !Index of the first z plane between 1:n3_iter
    boxit%i3s = boxit%dpbox_ptr%i3s + boxit%dpbox_ptr%i3xcsh

    select case(idpbox)
    case(DPB_POT_ION)
      boxit%n3_iter = boxit%dpbox_ptr%n3pi
    case(DPB_RHO)
      boxit%n3_iter = boxit%dpbox_ptr%n3d
    case(DPB_POT)
      boxit%n3_iter = boxit%dpbox_ptr%n3p
    case default
      call f_err_throw('dpbox_iterator: Wrong choice for the iterations over the z dimension', &
           err_name='BIGDFT_RUNTIME_ERROR')
    end select

    if (boxit%n3_iter== 0) then
      !No iteration, the iterator is destroyed and we leave!
      call nullify_dpbox_iterator(boxit)
      return
    end if

    !Conditions for periodicity in the three directions
    boxit%perx=(boxit%dpbox_ptr%geocode /= 'F')
    boxit%pery=(boxit%dpbox_ptr%geocode == 'P')
    boxit%perz=(boxit%dpbox_ptr%geocode /= 'F')

    !Calculate external buffers for each direction
    call ext_buffers(boxit%perx,boxit%nbl1,boxit%nbr1)
    call ext_buffers(boxit%pery,boxit%nbl2,boxit%nbr2)
    call ext_buffers(boxit%perz,boxit%nbl3,boxit%nbr3)

    if (present(nbox)) then
      !We iterate in a box around an atom
      boxit%whole = .false.
      boxit%nbox = nbox
    else
      !We iterate over the whole box
      boxit%whole = .true.
      boxit%nbox(1,3) = -boxit%nbl3
      boxit%nbox(2,3) = dpbox%ndims(3) - boxit%nbl3-1
      !ndims(2) contains nbr2
      boxit%nbox(1,2) = -boxit%nbl2
      boxit%nbox(2,2) = dpbox%ndims(2) - boxit%nbl2-1
      !ndims(1) contains nbr1
      boxit%nbox(1,1) = -boxit%nbl1
      boxit%nbox(2,1) = dpbox%ndims(1) - boxit%nbl1-1
    end if

    ! Start counting
    boxit%ix=0
    boxit%iy=0
    boxit%iz=0
    boxit%x=0.0_gp
    boxit%y=0.0_gp
    boxit%z=0.0_gp
    boxit%ind=0
    ! Iterate
    !First indices to change
    boxit%ibox(1) = boxit%nbox(1,1) - 1
    boxit%ibox(2) = boxit%nbox(1,2)
    boxit%ibox(3) = boxit%nbox(1,3)

  end function dpbox_iter


  !> Increment a valid iterator
  !! the control for validity has to be done outside
  !pure subroutine dpbox_refresh_iterator(boxit)
  !  implicit none
  !  type(dpbox_iterator), intent(inout) :: boxit
  !end subroutine dpbox_refresh_iterator


  !> Increment, and nullify if ended
  !! if the iterator is nullified, it does nothing
   pure subroutine increment_dpbox_iter(boxit)
     implicit none
     !Arguments
     type(dpbox_iterator), intent(inout) :: boxit
     !Local variables
     logical :: gox,goy,goz
    
    if (associated(boxit%dpbox_ptr)) then
      !There are distributed z planes in this proc: we start a loop
      loop_ind: do
        if (boxit%ibox(1) < boxit%nbox(2,1)) then
          boxit%ibox(1) = boxit%ibox(1) + 1
        else if (boxit%ibox(2) < boxit%nbox(2,2)) then
          !First index finished, increment the second one
          boxit%ibox(1) = boxit%nbox(1,1)
          boxit%ibox(2) = boxit%ibox(2) + 1
        else if (boxit%ibox(3) < boxit%nbox(2,3)) then
          !First and second indices finished, increment the last one
          boxit%ibox(1) = boxit%nbox(1,1)
          boxit%ibox(2) = boxit%nbox(1,2)
          boxit%ibox(3) = boxit%ibox(3) + 1
        else
          !End iteration, the iterator is destroyed and we leave!
          call nullify_dpbox_iterator(boxit)
          return
        end if
        !Check if this point is inside the box
        call ind_positions_new(boxit%perz,boxit%ibox(3),boxit%n3i,boxit%iz,goz) 
        boxit%iz = boxit%iz + boxit%nbl3 + 1
        if ( .not.(boxit%iz >= boxit%i3s .and. boxit%iz <= boxit%i3s+boxit%n3_iter-1) ) cycle
        call ind_positions_new(boxit%pery,boxit%ibox(2),boxit%n2i,boxit%iy,goy)
        if (.not.goy) cycle
        call ind_positions_new(boxit%perx,boxit%ibox(1),boxit%n1i,boxit%ix,gox)
        !Check if in the box
        !if (boxit%iz >= boxit%i3s .and. boxit%iz <= boxit%i3s+boxit%n3_iter-1 .and. goy .and. gox ) then
        if (gox) then
          !This point is valid: we calculate ind (index for pot_ion) and leave!
          boxit%ind = boxit%ix+1 + boxit%nbl1 &
                  & + (boxit%iy+boxit%nbl2)*boxit%n1i &
                  & + (boxit%iz-boxit%i3s)*boxit%n1i*boxit%n2i
          boxit%x = real(boxit%ibox(1),gp)*boxit%dpbox_ptr%hgrids(1)
          boxit%y = real(boxit%ibox(2),gp)*boxit%dpbox_ptr%hgrids(2)
          boxit%z = real(boxit%ibox(3),gp)*boxit%dpbox_ptr%hgrids(3)
          return
        end if
      end do loop_ind
    end if
   end subroutine increment_dpbox_iter


  !> Logical function, returns .true. if the iterator is still valid
  pure function dpbox_iter_is_valid(boxit)
    implicit none
    type(dpbox_iterator), intent(in) :: boxit
    logical :: dpbox_iter_is_valid
    
    dpbox_iter_is_valid=associated(boxit%dpbox_ptr)
  end function dpbox_iter_is_valid


  !> Logical function for iterating above atoms
  function dpbox_iter_next(boxit)
    implicit none
    type(dpbox_iterator), intent(inout) :: boxit
    logical :: dpbox_iter_next

    call increment_dpbox_iter(boxit)
    dpbox_iter_next = dpbox_iter_is_valid(boxit)
  end function dpbox_iter_next


  !> Calculate the size of the buffers in each direction
  subroutine ext_buffers(periodic,nl,nr)
    implicit none
    logical, intent(in) :: periodic !< Periodic or not
    integer, intent(out) :: nl,nr   !< Size of left and right buffer

    if (periodic) then
       nl=0
       nr=0
    else
       nl=14
       nr=15
    end if
  END SUBROUTINE ext_buffers


  !> Determine the index in which the potential must be inserted, following the BC
  !! Determine also whether the index is inside or outside the box for free BC
  !!!pure subroutine ind_positions(periodic,i,n,j,go)
  !!!  implicit none
  !!!  logical, intent(in) :: periodic !< Periodic or not
  !!!  integer, intent(in) :: i
  !!!  integer, intent(in) :: n      
  !!!  logical, intent(out) :: go      !< True if in the box
  !!!  integer, intent(out) :: j
  !!!
  !!!  if (periodic) then
  !!!     go=.true.
  !!!     j=modulo(i,2*n+2)
  !!!  else
  !!!     j=i
  !!!     if (i >= -14 .and. i <= 2*n+16) then
  !!!        go=.true.
  !!!     else
  !!!        go=.false.
  !!!     end if
  !!!  end if
  !!!
  !!!END SUBROUTINE ind_positions


  !> Determine the index in which the potential must be inserted, following the BC
  !! Determine also whether the index is inside or outside the box for free BC
  pure subroutine ind_positions_new(periodic,i,ni,j,go)
    implicit none
    logical, intent(in) :: periodic
    integer, intent(in) :: i,ni
    logical, intent(out) :: go
    integer, intent(out) :: j

    if (periodic) then
       go=.true.
       j=modulo(i,ni)
    else
       j=i
       if (i >= -14 .and. i <= ni-15) then
          go=.true.
       else
          go=.false.
       end if
    end if

  END SUBROUTINE ind_positions_new


end module module_dpbox
