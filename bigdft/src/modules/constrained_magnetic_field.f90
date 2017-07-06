!> @file
!! Handling of the constrained magnetic field of the system
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module module_cfd
  use module_base
  implicit none
  
  private

  !>associated to an instance of the cfd calculation
  type, public :: cfd_data
     !> number of magnetic centers
     integer :: nat=0
     !> position of the centers in the simulation domain
     real(gp), dimension(:,:), pointer :: rxyz => null()
     !> radius of each of the magnetic atoms
     real(gp), dimension(:), pointer :: radii => null()
     !> value of the magnetic field close to each of the centers
     real(gp), dimension(:,:), pointer :: B_at => null()
     !> local magnetization of each of the centers
     real(gp), dimension(:,:), pointer :: m_at => null()
     !> electronic charge inside each of the center
     real(gp), dimension(:), pointer :: rho_at => null()
  end type cfd_data

  public :: cfd_allocate,cfd_free,cfd_set_radius,cfd_dump_info
  public :: cfd_set_centers

contains
  
  pure function cfd_data_null() result(cfd)
    implicit none
    type(cfd_data) :: cfd
    call nullify_cfd_data(cfd)
  end function cfd_data_null
  pure subroutine nullify_cfd_data(cfd)
    implicit none
    type(cfd_data), intent(out) :: cfd
    cfd%nat=0
    nullify(cfd%rxyz)
    nullify(cfd%radii)
    nullify(cfd%B_at)
    nullify(cfd%m_at)
    nullify(cfd%rho_at)
  end subroutine nullify_cfd_data

  subroutine cfd_free(cfd)
    implicit none
    type(cfd_data), intent(inout) :: cfd
    
    call f_free_ptr(cfd%rxyz)
    call f_free_ptr(cfd%radii)
    call f_free_ptr(cfd%B_at)
    call f_free_ptr(cfd%m_at)
    call f_free_ptr(cfd%rho_at)
    call nullify_cfd_data(cfd)
  end subroutine cfd_free

  subroutine cfd_allocate(cfd,nat)
    implicit none
    integer, intent(in) :: nat
    type(cfd_data), intent(inout) :: cfd

    call cfd_free(cfd) !we can as the initial status of the data is defined
    
    cfd%nat=nat
    cfd%rxyz=f_malloc_ptr([3,nat],id='cfd%rxyz')
    cfd%radii=f_malloc_ptr(nat,id='cfd%radii')
    cfd%B_at=f_malloc_ptr([3,nat],id='cfd%B_at')
    cfd%m_at=f_malloc_ptr([3,nat],id='cfd%m_at')
    cfd%rho_at=f_malloc_ptr(nat,id='cfd%rho_at')

  end subroutine cfd_allocate

!!$  function cfd_get_centers_ptr(cfd) result(ptr)
!!$    implicit none
!!$    type(cfd_data), intent(inout) :: cfd
!!$    real(gp), dimension(:,:), pointer :: ptr
!!$
!!$    ptr => cfd%rxyz
!!$
!!$  end function cfd_get_centers_ptr

  subroutine cfd_set_centers(cfd,rxyz)
    implicit none
    type(cfd_data), intent(inout) :: cfd
    real(gp), dimension(3,cfd%nat), intent(in) :: rxyz
    
    call f_memcpy(src=rxyz,dest=cfd%rxyz)

  end subroutine cfd_set_centers
    

  pure subroutine cfd_set_radius(cfd,iat,radius)
    implicit none
    integer, intent(in) :: iat
    real(gp), intent(in) :: radius
    type(cfd_data), intent(inout) :: cfd

    cfd%radii(iat)=radius
  end subroutine cfd_set_radius

  subroutine cfd_dump_info(cfd)
    use yaml_output
    implicit none
    type(cfd_data), intent(in) :: cfd
    !local variables
    integer :: iat

    call yaml_newline()
    call yaml_sequence_open('Local information on the magnetic centers')
    do iat=1,cfd%nat
       call yaml_newline()
       call yaml_sequence(advance='no')
       call yaml_mapping_open(flow=.true.)
       call yaml_map('R',cfd%rxyz(:,iat)) !position
       call yaml_map('D',cfd%radii(iat)) !radius
       call yaml_newline()
       call yaml_map('M',cfd%m_at(:,iat),fmt='(1pe12.5)') !mag mom
       call yaml_map('C',cfd%rho_at(iat)) !charge
       call yaml_mapping_close()
    end do
    call yaml_sequence_close()
    call yaml_newline()

  end subroutine cfd_dump_info

end module module_cfd
