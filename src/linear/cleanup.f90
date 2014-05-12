!> @file
!! Linear version: deallocations
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module for the deallocation pointers used by linear version
module deallocatePointers

  implicit none

  interface checkAndDeallocatePointer
     module procedure checkAndDeallocatePointer_int_1, checkAndDeallocatePointer_sgl_1, &
                      checkAndDeallocatePointer_dbl_1, checkAndDeallocatePointer_log_1
     module procedure checkAndDeallocatePointer_int_2, checkAndDeallocatePointer_sgl_2, &
                      checkAndDeallocatePointer_dbl_2, checkAndDeallocatePointer_log_2
     module procedure checkAndDeallocatePointer_int_3, checkAndDeallocatePointer_sgl_3, &
                      checkAndDeallocatePointer_dbl_3, checkAndDeallocatePointer_log_3
     module procedure checkAndDeallocatePointer_int_4, checkAndDeallocatePointer_sgl_4, &
                      checkAndDeallocatePointer_dbl_4, checkAndDeallocatePointer_log_4
     module procedure checkAndDeallocatePointer_int_5, checkAndDeallocatePointer_sgl_5, &
                      checkAndDeallocatePointer_dbl_5, checkAndDeallocatePointer_log_5
     module procedure checkAndDeallocatePointer_int_6, checkAndDeallocatePointer_sgl_6, &
                      checkAndDeallocatePointer_dbl_6, checkAndDeallocatePointer_log_6
     module procedure checkAndDeallocatePointer_int_7, checkAndDeallocatePointer_sgl_7, &
                      checkAndDeallocatePointer_dbl_7, checkAndDeallocatePointer_log_7
  end interface

  contains


    subroutine checkAndDeallocatePointer_int_1(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      integer,dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatepointer_int_1


    subroutine checkAndDeallocatePointer_sgl_1(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(4),dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_sgl_1


    subroutine checkAndDeallocatePointer_dbl_1(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(8),dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_dbl_1


    subroutine checkAndDeallocatePointer_log_1(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      logical,dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_log_1


    subroutine checkAndDeallocatePointer_int_2(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      integer,dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatepointer_int_2


    subroutine checkAndDeallocatePointer_sgl_2(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(4),dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_sgl_2


    subroutine checkAndDeallocatePointer_dbl_2(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(8),dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_dbl_2


    subroutine checkAndDeallocatePointer_log_2(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      logical,dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_log_2


    subroutine checkAndDeallocatePointer_int_3(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatepointer_int_3


    subroutine checkAndDeallocatePointer_sgl_3(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_sgl_3


    subroutine checkAndDeallocatePointer_dbl_3(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_dbl_3


    subroutine checkAndDeallocatePointer_log_3(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_log_3


    subroutine checkAndDeallocatePointer_int_4(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatepointer_int_4


    subroutine checkAndDeallocatePointer_sgl_4(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_sgl_4


    subroutine checkAndDeallocatePointer_dbl_4(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_dbl_4


    subroutine checkAndDeallocatePointer_log_4(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_log_4



    subroutine checkAndDeallocatePointer_int_5(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatepointer_int_5


    subroutine checkAndDeallocatePointer_sgl_5(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_sgl_5


    subroutine checkAndDeallocatePointer_dbl_5(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_dbl_5


    subroutine checkAndDeallocatePointer_log_5(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_log_5


    subroutine checkAndDeallocatePointer_int_6(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatepointer_int_6


    subroutine checkAndDeallocatePointer_sgl_6(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_sgl_6


    subroutine checkAndDeallocatePointer_dbl_6(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_dbl_6


    subroutine checkAndDeallocatePointer_log_6(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_log_6


    subroutine checkAndDeallocatePointer_int_7(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatepointer_int_7


    subroutine checkAndDeallocatePointer_sgl_7(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_sgl_7


    subroutine checkAndDeallocatePointer_dbl_7(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_dbl_7


    subroutine checkAndDeallocatePointer_log_7(array, arrayname, subname)
      use module_base, only: memocc
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: arrayname, subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          if(size(array)>0) then
              iall=-product(shape(array))*kind(array)
              deallocate(array, stat=istat)
              call memocc(istat, iall, arrayname, subname)
          else
              nullify(array)
          end if
      end if

    end subroutine checkAndDeallocatePointer_log_7


end module deallocatePointers


!!$!> Module to deallocate the interfaces for the linear scaling version
!!$module deallocationInterfaces
!!$  implicit none
!!$
!!$  interface
!!$    subroutine deallocate_local_zone_descriptors(lzd, subname)
!!$      use module_types
!!$      implicit none
!!$      type(local_zone_descriptors),intent(inout):: lzd
!!$      character(len=*),intent(in):: subname
!!$    end subroutine deallocate_local_zone_descriptors
!!$    
!!$    subroutine deallocate_orbitals_data(orbs, subname)
!!$      use module_types
!!$      use deallocatePointers
!!$      implicit none
!!$      type(orbitals_data),intent(inout):: orbs
!!$      character(len=*),intent(in):: subname
!!$    end subroutine deallocate_orbitals_data
!!$    
!!$    subroutine deallocate_comms_cubic(comms, subname)
!!$      
!!$      use module_types
!!$      use deallocatePointers
!!$      implicit none
!!$      type(comms_cubic),intent(inout):: comms
!!$      character(len=*),intent(in):: subname
!!$    end subroutine deallocate_comms_cubic
!!$    
!!$  end interface
!!$
!!$end module deallocationInterfaces

subroutine deallocate_local_zone_descriptors(lzd, subname)
  use locregs, only: deallocate_locreg_descriptors
  use module_types, only: local_zone_descriptors
  implicit none
  
  ! Calling arguments
  type(local_zone_descriptors),intent(inout):: lzd
  character(len=*),intent(in):: subname
  
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


subroutine deallocate_Lzd_except_Glr(lzd, subname)
  use locregs, only: deallocate_locreg_descriptors
  use module_types, only: local_zone_descriptors
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(inout):: lzd
  character(len=*),intent(in):: subname

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


subroutine deallocate_orbitals_data(orbs, subname)
  use module_types
  implicit none
  
  ! Calling arguments
  type(orbitals_data),intent(inout):: orbs
  character(len=*),intent(in):: subname
  
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


subroutine deallocate_comms_cubic(comms, subname)
  use communications_base, only: comms_cubic
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(comms_cubic),intent(inout):: comms
  character(len=*),intent(in):: subname
  
  call checkAndDeallocatePointer(comms%ncntd, 'comms%ncntd', subname)
  call checkAndDeallocatePointer(comms%ncntt, 'comms%ncntt', subname)
  call checkAndDeallocatePointer(comms%ndspld, 'comms%ndspld', subname)
  call checkAndDeallocatePointer(comms%ndsplt, 'comms%ndsplt', subname)
  call checkAndDeallocatePointer(comms%nvctr_par, 'comms%nvctr_par', subname)
  
end subroutine deallocate_comms_cubic


subroutine deallocate_convolutions_bounds(bounds, subname)
  
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(convolutions_bounds),intent(inout):: bounds
  character(len=*),intent(in):: subname

  call f_free_ptr(bounds%ibyyzz_r)

  call deallocate_kinetic_bounds(bounds%kb, subname)
  call deallocate_shrink_bounds(bounds%sb, subname)
  call deallocate_grow_bounds(bounds%gb, subname)
  
end subroutine deallocate_convolutions_bounds


subroutine deallocate_kinetic_bounds(kb, subname)
  
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(kinetic_bounds),intent(inout):: kb
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(kb%ibyz_c, 'kb%ibyz_c', subname)
  call checkAndDeallocatePointer(kb%ibxz_c, 'kb%ibxz_c', subname)
  call checkAndDeallocatePointer(kb%ibxy_c, 'kb%ibxy_c', subname)
  call checkAndDeallocatePointer(kb%ibyz_f, 'kb%ibyz_f', subname)
  call checkAndDeallocatePointer(kb%ibxz_f, 'kb%ibxz_f', subname)
  call checkAndDeallocatePointer(kb%ibxy_f, 'kb%ibxy_f', subname)

end subroutine deallocate_kinetic_bounds


subroutine deallocate_shrink_bounds(sb, subname)
  
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(shrink_bounds),intent(inout):: sb
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(sb%ibzzx_c, 'sb%ibzzx_c', subname)
  call checkAndDeallocatePointer(sb%ibyyzz_c, 'sb%ibyyzz_c', subname)
  call checkAndDeallocatePointer(sb%ibxy_ff, 'sb%ibxy_ff,', subname)
  call checkAndDeallocatePointer(sb%ibzzx_f, 'sb%ibzzx_f,', subname)
  call checkAndDeallocatePointer(sb%ibyyzz_f, 'sb%ibyyzz_f,', subname)

end subroutine deallocate_shrink_bounds


subroutine deallocate_grow_bounds(gb, subname)
  
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(grow_bounds),intent(inout):: gb
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(gb%ibzxx_c, 'gb%ibzxx_c', subname)
  call checkAndDeallocatePointer(gb%ibxxyy_c, 'gb%ibxxyy_c', subname)
  call checkAndDeallocatePointer(gb%ibyz_ff, 'gb%ibyz_ff', subname)
  call checkAndDeallocatePointer(gb%ibzxx_f, 'gb%ibzxx_f', subname)
  call checkAndDeallocatePointer(gb%ibxxyy_f, 'gb%ibxxyy_f', subname)

end subroutine deallocate_grow_bounds






