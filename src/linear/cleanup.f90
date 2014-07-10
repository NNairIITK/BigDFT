!> @file
!! Linear version: deallocations
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine deallocate_local_zone_descriptors(lzd)
  use locregs, only: deallocate_locreg_descriptors
  use module_types, only: local_zone_descriptors
  implicit none
  
  ! Calling arguments
  type(local_zone_descriptors),intent(inout):: lzd
  ! Local variables
  integer:: iis1, iie1, i1

  call deallocate_locreg_descriptors(lzd%Glr)

  if(associated(lzd%llr)) then  
     iis1=lbound(lzd%llr,1)
     iie1=ubound(lzd%llr,1)
     do i1=iis1,iie1
         call deallocate_locreg_descriptors(lzd%llr(i1))
     end do
     deallocate(lzd%llr)
     nullify(lzd%llr)
  end if

end subroutine deallocate_local_zone_descriptors


subroutine deallocate_Lzd_except_Glr(lzd)
  use locregs, only: deallocate_locreg_descriptors
  use module_types, only: local_zone_descriptors
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(inout):: lzd
  ! Local variables
  integer:: iis1, iie1, i1

  if(associated(lzd%llr)) then
     iis1=lbound(lzd%llr,1)
     iie1=ubound(lzd%llr,1)
     do i1=iis1,iie1
         call deallocate_locreg_descriptors(lzd%llr(i1))
     end do
     deallocate(lzd%llr)
     nullify(lzd%llr)
  end if

end subroutine deallocate_Lzd_except_Glr


subroutine deallocate_orbitals_data(orbs)
  use module_types
  implicit none
  
  ! Calling arguments
  type(orbitals_data),intent(inout):: orbs
  
  call f_free_ptr(orbs%norb_par)
  call f_free_ptr(orbs%iokpt)
  call f_free_ptr(orbs%ikptproc)
  call f_free_ptr(orbs%inwhichlocreg)
  call f_free_ptr(orbs%onwhichatom)
  call f_free_ptr(orbs%isorb_par)
  call f_free_ptr(orbs%eval)
  call f_free_ptr(orbs%occup)
  call f_free_ptr(orbs%spinsgn)
  call f_free_ptr(orbs%kwgts)
  call f_free_ptr(orbs%kpts)
  call f_free_ptr(orbs%ispot)
  
end subroutine deallocate_orbitals_data


subroutine deallocate_comms_cubic(comms)
  use module_base
  use communications_base, only: comms_cubic
  implicit none
  
  ! Calling arguments
  type(comms_cubic),intent(inout):: comms
  
  call f_free_ptr(comms%ncntd)
  call f_free_ptr(comms%ncntt)
  call f_free_ptr(comms%ndspld)
  call f_free_ptr(comms%ndsplt)
  call f_free_ptr(comms%nvctr_par)
  
end subroutine deallocate_comms_cubic


subroutine deallocate_convolutions_bounds(bounds)
  
  use module_types
  implicit none
  
  ! Calling arguments
  type(convolutions_bounds),intent(inout):: bounds

  call f_free_ptr(bounds%ibyyzz_r)

  call deallocate_kinetic_bounds(bounds%kb)
  call deallocate_shrink_bounds(bounds%sb)
  call deallocate_grow_bounds(bounds%gb)
  
end subroutine deallocate_convolutions_bounds


subroutine deallocate_kinetic_bounds(kb)
  
  use module_types
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
  
  use module_types
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

  use module_types
  implicit none
 
  ! Calling arguments
  type(grow_bounds),intent(inout):: gb

  call f_free_ptr(gb%ibzxx_c)
  call f_free_ptr(gb%ibxxyy_c)
  call f_free_ptr(gb%ibyz_ff)
  call f_free_ptr(gb%ibzxx_f)
  call f_free_ptr(gb%ibxxyy_f)

end subroutine deallocate_grow_bounds
