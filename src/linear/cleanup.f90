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
!!$    subroutine deallocate_communications_arrays(comms, subname)
!!$      
!!$      use module_types
!!$      use deallocatePointers
!!$      implicit none
!!$      type(communications_arrays),intent(inout):: comms
!!$      character(len=*),intent(in):: subname
!!$    end subroutine deallocate_communications_arrays
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
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(orbitals_data),intent(inout):: orbs
  character(len=*),intent(in):: subname
  
  call checkAndDeallocatePointer(orbs%norb_par, 'orbs%norb_par', subname)
  call checkAndDeallocatePointer(orbs%iokpt, 'orbs%iokpt', subname)
  call checkAndDeallocatePointer(orbs%ikptproc, 'orbs%ikptproc', subname)
  call checkAndDeallocatePointer(orbs%inwhichlocreg, 'orbs%inwhichlocreg', subname)
  call checkAndDeallocatePointer(orbs%onwhichatom, 'orbs%onwhichatom', subname)
  call checkAndDeallocatePointer(orbs%onwhichatom, 'orbs%onwhichfragment', subname)
  call checkAndDeallocatePointer(orbs%isorb_par, 'orbs%isorb_par', subname)
  call checkAndDeallocatePointer(orbs%eval, 'orbs%eval', subname)
  call checkAndDeallocatePointer(orbs%occup, 'orbs%occup', subname)
  call checkAndDeallocatePointer(orbs%spinsgn, 'orbs%spinsgn', subname)
  call checkAndDeallocatePointer(orbs%kwgts, 'orbs%kwgts', subname)
  call checkAndDeallocatePointer(orbs%kpts, 'orbs%kpts', subname)
  call checkAndDeallocatePointer(orbs%ispot, 'orbs%ispot', subname)
  
end subroutine deallocate_orbitals_data


subroutine deallocate_communications_arrays(comms, subname)
  
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(communications_arrays),intent(inout):: comms
  character(len=*),intent(in):: subname
  
  call checkAndDeallocatePointer(comms%ncntd, 'comms%ncntd', subname)
  call checkAndDeallocatePointer(comms%ncntt, 'comms%ncntt', subname)
  call checkAndDeallocatePointer(comms%ndspld, 'comms%ndspld', subname)
  call checkAndDeallocatePointer(comms%ndsplt, 'comms%ndsplt', subname)
  call checkAndDeallocatePointer(comms%nvctr_par, 'comms%nvctr_par', subname)
  
end subroutine deallocate_communications_arrays


subroutine deallocate_convolutions_bounds(bounds, subname)
  
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(convolutions_bounds),intent(inout):: bounds
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(bounds%ibyyzz_r, 'bounds%ibyyzz_r', subname)

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






subroutine deallocate_p2pComms(p2pcomm, subname)
  
  use module_types
  use deallocatePointers
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: p2pcomm
  character(len=*),intent(in):: subname
  
  ! Local variables
  integer :: is, ie, i, ierr

  call checkAndDeallocatePointer(p2pcomm%noverlaps, 'p2pcomm%noverlaps', subname)
  call checkAndDeallocatePointer(p2pcomm%recvBuf, 'p2pcomm%recvBuf', subname)
  call checkAndDeallocatePointer(p2pcomm%comarr, 'p2pcomm%comarr', subname)
  call checkAndDeallocatePointer(p2pcomm%ise, 'p2pcomm%ise', subname)

  if (.not.p2pcomm%communication_complete) then
      stop 'cannot deallocate mpi data types if communication has not completed'
  end if
  !!if (associated(p2pcomm%mpi_datatypes)) then
  !!    is=lbound(p2pcomm%mpi_datatypes,2)
  !!    ie=ubound(p2pcomm%mpi_datatypes,2)
  !!    do i=is,ie
  !!        if (p2pcomm%mpi_datatypes(2,i)==1) then
  !!             call mpi_type_free(p2pcomm%mpi_datatypes(1,i), ierr)
  !!         end if
  !!    end do
  !!end if
  call checkAndDeallocatePointer(p2pcomm%mpi_datatypes, 'p2pcomm%mpi_datatypes', subname)

end subroutine deallocate_p2pComms

subroutine deallocate_foe(foe_obj, subname)
  
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_foe
  implicit none
  
  ! Calling arguments
  type(foe_data),intent(inout):: foe_obj
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(foe_obj%kernel_nseg, 'foe_obj%kernel_nseg', subname)
  call checkAndDeallocatePointer(foe_obj%kernel_segkeyg, 'foe_obj%kernel_segkeyg', subname)

end subroutine deallocate_foe

subroutine deallocate_sparseMatrix(sparsemat, subname)
  
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_sparseMatrix
  implicit none
  
  ! Calling arguments
  type(sparseMatrix),intent(inout):: sparsemat
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(sparseMat%keyg, 'sparseMat%keyg', subname)
  call checkAndDeallocatePointer(sparseMat%keyv, 'sparseMat%keyv', subname)
  call checkAndDeallocatePointer(sparseMat%nsegline, 'sparseMat%nsegline', subname)
  call checkAndDeallocatePointer(sparseMat%istsegline, 'sparseMat%istsegline', subname)
  call checkAndDeallocatePointer(sparseMat%matrix_compr, 'sparseMat%matrix_compr', subname)
  call checkAndDeallocatePointer(sparseMat%matrix, 'sparseMat%matrix', subname)
  !call checkAndDeallocatePointer(sparseMat%matrixindex_in_compressed, 'sparseMat%matrixindex_in_compressed', subname)
  call checkAndDeallocatePointer(sparseMat%matrixindex_in_compressed_arr, 'sparseMat%matrixindex_in_compressed_arr', subname)
  call checkAndDeallocatePointer(sparseMat%orb_from_index, 'sparseMat%orb_from_index', subname)
  call checkAndDeallocatePointer(sparseMat%matrixindex_in_compressed_fortransposed, &
       'sparseMat%matrixindex_in_compressed_fortransposed', subname)
  call f_free_ptr(sparseMat%isvctr_par)
  call f_free_ptr(sparseMat%nvctr_par)
  call f_free_ptr(sparseMat%isfvctr_par)
  call f_free_ptr(sparseMat%nfvctr_par)
  call f_free_ptr(sparseMat%matrixp)
  call f_free_ptr(sparseMat%matrix_comprp)

end subroutine deallocate_sparseMatrix



subroutine deallocate_collective_comms(collcom, subname)
  
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_collective_comms
  implicit none
  
  ! Calling arguments
  type(collective_comms),intent(inout):: collcom
  character(len=*),intent(in):: subname

  ! Local variables

  call checkAndDeallocatePointer(collcom%nsendcounts_c, 'collcom%nsendcounts_c', subname)
  call checkAndDeallocatePointer(collcom%nsenddspls_c, 'collcom%nsenddspls_c', subname)
  call checkAndDeallocatePointer(collcom%nrecvcounts_c, 'collcom%nrecvcounts_c', subname)
  call checkAndDeallocatePointer(collcom%nrecvdspls_c, 'collcom%nrecvdspls_c', subname)
  call checkAndDeallocatePointer(collcom%isendbuf_c, 'collcom%isendbuf_c', subname)
  call checkAndDeallocatePointer(collcom%iextract_c, 'collcom%iextract_c', subname)
  call checkAndDeallocatePointer(collcom%iexpand_c, 'collcom%iexpand_c', subname)
  call checkAndDeallocatePointer(collcom%irecvbuf_c, 'collcom%irecvbuf_c', subname)
  call checkAndDeallocatePointer(collcom%norb_per_gridpoint_c, 'collcom%norb_per_gridpoint_c', subname)
  call checkAndDeallocatePointer(collcom%indexrecvorbital_c, 'collcom%indexrecvorbital_c', subname)
  call checkAndDeallocatePointer(collcom%isptsp_c, 'collcom%isptsp_c', subname)
  call checkAndDeallocatePointer(collcom%psit_c, 'collcom%psit_c', subname)
  call checkAndDeallocatePointer(collcom%nsendcounts_f, 'collcom%nsendcounts_f', subname)
  call checkAndDeallocatePointer(collcom%nsenddspls_f, 'collcom%nsenddspls_f', subname)
  call checkAndDeallocatePointer(collcom%nrecvcounts_f, 'collcom%nrecvcounts_f', subname)
  call checkAndDeallocatePointer(collcom%nrecvdspls_f, 'collcom%nrecvdspls_f', subname)
  call checkAndDeallocatePointer(collcom%isendbuf_f, 'collcom%isendbuf_f', subname)
  call checkAndDeallocatePointer(collcom%iextract_f, 'collcom%iextract_f', subname)
  call checkAndDeallocatePointer(collcom%iexpand_f, 'collcom%iexpand_f', subname)
  call checkAndDeallocatePointer(collcom%irecvbuf_f, 'collcom%irecvbuf_f', subname)
  call checkAndDeallocatePointer(collcom%norb_per_gridpoint_f, 'collcom%norb_per_gridpoint_f', subname)
  call checkAndDeallocatePointer(collcom%indexrecvorbital_f, 'collcom%indexrecvorbital_f', subname)
  call checkAndDeallocatePointer(collcom%isptsp_f, 'collcom%isptsp_f', subname)
  call checkAndDeallocatePointer(collcom%psit_f, 'collcom%psit_f', subname)
  call checkAndDeallocatePointer(collcom%nsendcounts_repartitionrho, 'collcom%nsendcounts_repartitionrho', subname)
  call checkAndDeallocatePointer(collcom%nrecvcounts_repartitionrho, 'collcom%nrecvcounts_repartitionrho', subname)
  call checkAndDeallocatePointer(collcom%nsenddspls_repartitionrho, 'collcom%nsenddspls_repartitionrho', subname)
  call checkAndDeallocatePointer(collcom%nrecvdspls_repartitionrho, 'collcom%nrecvdspls_repartitionrho', subname)
  call checkAndDeallocatePointer(collcom%commarr_repartitionrho, 'collcom%commarr_repartitionrho', subname)

end subroutine deallocate_collective_comms

