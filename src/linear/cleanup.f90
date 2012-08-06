!> @file
!! Linear version: deallocations
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Deallocation pointers used by linear version
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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
      use module_base
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


module deallocationInterfaces
  implicit none

  interface


    subroutine deallocate_local_zone_descriptors(lzd, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(local_zone_descriptors),intent(inout):: lzd
      character(len=*),intent(in):: subname
    end subroutine deallocate_local_zone_descriptors
    
    subroutine deallocate_orbitals_data(orbs, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(orbitals_data),intent(inout):: orbs
      character(len=*),intent(in):: subname
    end subroutine deallocate_orbitals_data
    
    subroutine deallocate_communications_arrays(comms, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(communications_arrays),intent(inout):: comms
      character(len=*),intent(in):: subname
    end subroutine deallocate_communications_arrays
    
    subroutine deallocate_locreg_descriptors(lr, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(locreg_descriptors),intent(inout):: lr
      character(len=*),intent(in):: subname
    end subroutine deallocate_locreg_descriptors
    
    subroutine deallocate_wavefunctions_descriptors(wfd, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(wavefunctions_descriptors),intent(inout):: wfd
      character(len=*),intent(in):: subname
    end subroutine deallocate_wavefunctions_descriptors
    
    subroutine deallocate_convolutions_bounds(bounds, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(convolutions_bounds),intent(inout):: bounds
      character(len=*),intent(in):: subname
    end subroutine deallocate_convolutions_bounds
    
    subroutine deallocate_kinetic_bounds(kb, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(kinetic_bounds),intent(inout):: kb
      character(len=*),intent(in):: subname
    end subroutine deallocate_kinetic_bounds
    
    subroutine deallocate_shrink_bounds(sb, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(shrink_bounds),intent(inout):: sb
      character(len=*),intent(in):: subname
    end subroutine deallocate_shrink_bounds
    
    subroutine deallocate_grow_bounds(gb, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(grow_bounds),intent(inout):: gb
      character(len=*),intent(in):: subname
    end subroutine deallocate_grow_bounds
    
    subroutine deallocate_nonlocal_psp_descriptors(nlpspd, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(nonlocal_psp_descriptors),intent(inout):: nlpspd
      character(len=*),intent(in):: subname
    end subroutine deallocate_nonlocal_psp_descriptors
    
    subroutine deallocate_matrixMinimization(matmin, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(matrixMinimization),intent(inout):: matmin
      character(len=*),intent(in):: subname
    end subroutine deallocate_matrixMinimization
    
    subroutine deallocate_matrixLocalizationRegion(mlr, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(matrixLocalizationRegion),intent(inout):: mlr
      character(len=*),intent(in):: subname
    end subroutine deallocate_matrixLocalizationRegion


  end interface

end module deallocationInterfaces

subroutine deallocate_local_zone_descriptors(lzd, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_local_zone_descriptors
  implicit none
  
  ! Calling arguments
  type(local_zone_descriptors),intent(inout):: lzd
  character(len=*),intent(in):: subname
  
  ! Local variables
  integer:: iis1, iie1, i1

!  call checkAndDeallocatePointer(lzd%Glr%projflg, 'lzd%Glr%projflg', subname)
  call checkAndDeallocatePointer(lzd%doHamAppl, 'lzd%doHamAppl', subname)
  call deallocate_locreg_descriptors(lzd%Glr, subname)

  if(associated(lzd%llr)) then  
     iis1=lbound(lzd%llr,1)
     iie1=ubound(lzd%llr,1)
     !write(*,*) 'iis1,iie1',iis1,iie1
     do i1=iis1,iie1
!         call checkAndDeallocatePointer(lzd%llr(i1)%projflg, 'lzd%llr(i1)%projflg', subname)
         !write(*,*) 'i1',i1
         call deallocate_locreg_descriptors(lzd%llr(i1), subname)
     end do
     deallocate(lzd%llr)
     nullify(lzd%llr)
  end if

  !!call checkAndDeallocatePointer(lzd%cutoffweight, 'cutoffweight', subname)

end subroutine deallocate_local_zone_descriptors


subroutine deallocate_Lzd_except_Glr(lzd, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_Lzd_except_Glr
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(inout):: lzd
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: iis1, iie1, i1


!  call deallocate_orbitals_data(lzd%orbs, subname)

!  iis1=lbound(lzd%lorbs,1)
!  iie1=ubound(lzd%lorbs,1)
!  do i1=iis1,iie1
!      call deallocate_orbitals_data(lzd%lorbs(i1), subname)
!  end do

!  call deallocate_communications_arrays(lzd%comms, subname)

!  call checkAndDeallocatePointer(lzd%Glr%projflg, 'lzd%Glr%projflg', subname)
  call checkAndDeallocatePointer(lzd%doHamAppl, 'lzd%doHamAppl', subname)
!  call deallocate_locreg_descriptors(lzd%Glr, subname)

!  call deallocate_nonlocal_psp_descriptors(lzd%Gnlpspd, subname)

  if(associated(lzd%llr)) then
     iis1=lbound(lzd%llr,1)
     iie1=ubound(lzd%llr,1)
     !write(*,*) 'iis1,iie1',iis1,iie1
     do i1=iis1,iie1
         !if(associated(lzd%llr(i1)%projflg)) then
         !    nullify(lzd%llr(i1)%projflg)
         !end if
!         call checkAndDeallocatePointer(lzd%llr(i1)%projflg, 'lzd%llr(i1)%projflg', subname)
         !write(*,*) 'i1',i1
         call deallocate_locreg_descriptors(lzd%llr(i1), subname)
     end do
     deallocate(lzd%llr)
     nullify(lzd%llr)
  end if

end subroutine deallocate_Lzd_except_Glr


subroutine deallocate_orbitals_data(orbs, subname)
  use module_base
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
  !call checkAndDeallocatePointer(orbs%inWhichLocregP, 'orbs%inWhichLocregP', subname)
  call checkAndDeallocatePointer(orbs%onWhichMPI, 'orbs%onWhichMPI', subname)
  call checkAndDeallocatePointer(orbs%isorb_par, 'orbs%isorb_par', subname)
  call checkAndDeallocatePointer(orbs%eval, 'orbs%eval', subname)
  call checkAndDeallocatePointer(orbs%occup, 'orbs%occup', subname)
  call checkAndDeallocatePointer(orbs%spinsgn, 'orbs%spinsgn', subname)
  call checkAndDeallocatePointer(orbs%kwgts, 'orbs%kwgts', subname)
  call checkAndDeallocatePointer(orbs%kpts, 'orbs%kpts', subname)
  call checkAndDeallocatePointer(orbs%ispot, 'orbs%ispot', subname)
  
end subroutine deallocate_orbitals_data


subroutine deallocate_communications_arrays(comms, subname)
  use module_base
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



subroutine deallocate_locreg_descriptors(lr, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_locreg_descriptors
  implicit none
  
  ! Calling arguments
  type(locreg_descriptors),intent(inout):: lr
  character(len=*),intent(in):: subname

!  call checkAndDeallocatePointer(lr%projflg,'lr%projflg', subname)

  call deallocate_wavefunctions_descriptors(lr%wfd, subname)

  call deallocate_convolutions_bounds(lr%bounds, subname)
  
end subroutine deallocate_locreg_descriptors

subroutine deallocate_locreg_descriptors2(lr,subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_locreg_descriptors
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(inout):: lr
  character(len=*),intent(in):: subname

!  call checkAndDeallocatePointer(lr%projflg, 'lr%projflg', subname)

  call deallocate_wavefunctions_descriptors(lr%wfd, subname)
!  call deallocate_convolutions_bounds(lr%bounds, subname)
! Don't need to deallocate the bounds, since they are not associated for
! overlap regions

end subroutine deallocate_locreg_descriptors2



subroutine deallocate_wavefunctions_descriptors(wfd, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(wavefunctions_descriptors),intent(inout):: wfd
  character(len=*),intent(in):: subname

  if(associated(wfd%keyglob, target = wfd%keygloc)) then
     call checkAndDeallocatePointer(wfd%keyglob, 'wfd%keyglob', subname)
     nullify(wfd%keyglob)
  else
     call checkAndDeallocatePointer(wfd%keyglob, 'wfd%keyglob', subname)
     nullify(wfd%keyglob)
     call checkAndDeallocatePointer(wfd%keygloc, 'wfd%keygloc', subname)
     nullify(wfd%keygloc)
  end if
  if(associated(wfd%keyvglob, target = wfd%keyvloc)) then
     call checkAndDeallocatePointer(wfd%keyvglob, 'wfd%keyvglob', subname)
     nullify(wfd%keyvglob)
  else
     call checkAndDeallocatePointer(wfd%keyvglob, 'wfd%keyvglob', subname)
     nullify(wfd%keyvglob)
     call checkAndDeallocatePointer(wfd%keyvloc, 'wfd%keyvloc', subname)
     nullify(wfd%keyvloc)
  end if
end subroutine deallocate_wavefunctions_descriptors


subroutine deallocate_convolutions_bounds(bounds, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_convolutions_bounds
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
  use module_base
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(kinetic_bounds),intent(inout):: kb
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(kb%ibyz_c, 'kb%ibyz_c', subname)
  call checkAndDeallocatePointer(kb%ibxz_c, 'kb%ibxz_c', subname)
  !!write(*,*) 'check, size',associated(kb%ibxz_c), size(kb%ibxz_c)
  !!if(associated(kb%ibxz_c)) then
  !!  deallocate(kb%ibxz_c)
  !!  !nullify(kb%ibxz_c)
  !!end if
  call checkAndDeallocatePointer(kb%ibxy_c, 'kb%ibxy_c', subname)
  !if(associated(kb%ibxy_c)) then
  !    if(size(kb%ibxy_c)>0) then
  !        deallocate(kb%ibxy_c)
  !    else
  !        nullify(kb%ibxy_c)
  !    end if
  !end if
  !write(*,*) 'check', associated(kb%ibyz_f)
  !if(associated(kb%ibyz_f)) then
  !    write(*,*) 'size',size(kb%ibyz_f)
  !end if
  call checkAndDeallocatePointer(kb%ibyz_f, 'kb%ibyz_f', subname)
  call checkAndDeallocatePointer(kb%ibxz_f, 'kb%ibxz_f', subname)
  call checkAndDeallocatePointer(kb%ibxy_f, 'kb%ibxy_f', subname)

end subroutine deallocate_kinetic_bounds


subroutine deallocate_shrink_bounds(sb, subname)
  use module_base
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
  use module_base
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


subroutine deallocate_nonlocal_psp_descriptors(nlpspd, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(nonlocal_psp_descriptors),intent(inout):: nlpspd
  character(len=*),intent(in):: subname
  integer :: i_stat,iat

  do iat=1,nlpspd%natoms
     call deallocate_wfd(nlpspd%plr(iat)%wfd,subname)
  end do
  if (nlpspd%natoms /=0) then
     deallocate(nlpspd%plr,stat=i_stat)
     if (i_stat /= 0) stop 'plr deallocation error'
     nlpspd%natoms=0
  end if
  nullify(nlpspd%plr)

!!$  call checkAndDeallocatePointer(nlpspd%nvctr_p, 'nlpspd%nvctr_p', subname)
!!$  call checkAndDeallocatePointer(nlpspd%nseg_p, 'nlpspd%nseg_p', subname)
!!$  call checkAndDeallocatePointer(nlpspd%keyv_p, 'nlpspd%keyv_p', subname)
!!$  call checkAndDeallocatePointer(nlpspd%keyg_p, 'nlpspd%keyg_p', subname)
!!$  call checkAndDeallocatePointer(nlpspd%nboxp_c, 'nlpspd%nboxp_c', subname)
!!$  call checkAndDeallocatePointer(nlpspd%nboxp_f, 'nlpspd%nboxp_f', subname)

end subroutine deallocate_nonlocal_psp_descriptors



subroutine deallocate_matrixMinimization(matmin, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_matrixMinimization
  implicit none
  
  ! Calling arguments
  type(matrixMinimization),intent(inout):: matmin
  character(len=*),intent(in):: subname
  
  ! Local variables
  integer:: iis1, iie1, i1
  
  iis1=lbound(matmin%mlr,1)
  iie1=ubound(matmin%mlr,1)
  do i1=iis1,iie1
      call deallocate_matrixLocalizationRegion(matmin%mlr(i1), subname)
  end do
  deallocate(matmin%mlr)
  nullify(matmin%mlr)
  
  call checkAndDeallocatePointer(matmin%inWhichLocregExtracted, 'matmin%inWhichLocregExtracted', subname)
  
  call checkAndDeallocatePointer(matmin%inWhichLocregOnMPI, 'matmin%inWhichLocregOnMPI', subname)
  
  call checkAndDeallocatePointer(matmin%indexInLocreg, 'matmin%indexInLocreg', subname)

end subroutine deallocate_matrixMinimization



subroutine deallocate_matrixLocalizationRegion(mlr, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(matrixLocalizationRegion),intent(inout):: mlr
  character(len=*),intent(in):: subname
  
  call checkAndDeallocatePointer(mlr%indexInGlobal, 'mlr%indexInGlobal', subname)
  
end subroutine deallocate_matrixLocalizationRegion

subroutine deallocate_p2pComms(p2pcomm, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: p2pcomm
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(p2pcomm%noverlaps, 'p2pcomm%noverlaps', subname)
  call checkAndDeallocatePointer(p2pcomm%overlaps, 'p2pcomm%overlaps', subname)
  !!call checkAndDeallocatePointer(p2pcomm%istarr, 'p2pcomm%istarr', subname)
  !!call checkAndDeallocatePointer(p2pcomm%istrarr, 'p2pcomm%istrarr', subname)
  call checkAndDeallocatePointer(p2pcomm%sendBuf, 'p2pcomm%sendBuf', subname)
  call checkAndDeallocatePointer(p2pcomm%recvBuf, 'p2pcomm%recvBuf', subname)
  call checkAndDeallocatePointer(p2pcomm%comarr, 'p2pcomm%comarr', subname)
  !!call checkAndDeallocatePointer(p2pcomm%communComplete, 'p2pcomm%communComplete', subname)
  !!call checkAndDeallocatePointer(p2pcomm%computComplete, 'p2pcomm%computComplete', subname)
  !!call checkAndDeallocatePointer(p2pcomm%startingindex, 'p2pcomm%startingindex', subname)
  call checkAndDeallocatePointer(p2pcomm%ise3, 'p2pcomm%ise3', subname)
  call checkAndDeallocatePointer(p2pcomm%requests, 'p2pcomm%requests', subname)

end subroutine deallocate_p2pComms


!!!subroutine deallocate_largeBasis(lb, subname)
!!!  use module_base
!!!  use module_types
!!!  use deallocatePointers
!!!  use module_interfaces, exceptThisOne => deallocate_largeBasis
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  type(largeBasis),intent(inout):: lb
!!!  character(len=*),intent(in):: subname
!!!
!!!  call deallocate_communications_arrays(lb%comms, subname)
!!!  call deallocate_communications_arrays(lb%gcomms, subname)
!!!  call deallocate_orbitals_data(lb%orbs, subname)
!!!  call deallocate_orbitals_data(lb%gorbs, subname)
!!!  !call deallocate_local_zone_descriptors(lb%lzd, subname)
!!!  !call deallocate_p2pCommsRepartition(lb%comrp, subname)
!!!  call deallocate_p2pComms(lb%comrp, subname)
!!!  !call deallocate_p2pCommsOrthonormality(lb%comon, subname)
!!!  call deallocate_p2pComms(lb%comon, subname)
!!!  call deallocate_overlapParameters(lb%op, subname)
!!!  !call deallocate_p2pCommsGatherPot(lb%comgp, subname)
!!!  call deallocate_p2pComms(lb%comgp, subname)
!!!  call deallocate_matrixDescriptors(lb%mad, subname)
!!!  call deallocate_collectiveComms(lb%collComms, subname)
!!!  !call deallocate_p2pCommsSumrho(lb%comsr, subname)
!!!  call deallocate_p2pComms(lb%comsr, subname)
!!!
!!!
!!!end subroutine deallocate_largeBasis


subroutine deallocate_overlapParameters(op, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_overlapParameters
  implicit none
  
  ! Calling arguments
  type(overlapParameters),intent(inout):: op
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: iis1, iie1, iis2, iie2, i1, i2

  call checkAndDeallocatePointer(op%noverlaps, 'op%noverlaps', subname)
  call checkAndDeallocatePointer(op%overlaps, 'op%overlaps', subname)
  call checkAndDeallocatePointer(op%indexInRecvBuf, 'op%indexInRecvBuf', subname)
  call checkAndDeallocatePointer(op%indexInSendBuf, 'op%indexInSendBuf', subname)


if(associated(op%wfd_overlap)) then
   iis1=lbound(op%wfd_overlap,1)
   iie1=ubound(op%wfd_overlap,1)
   iis2=lbound(op%wfd_overlap,2)
   iie2=ubound(op%wfd_overlap,2)
   do i2=iis2,iie2
       do i1=iis1,iie1
              call deallocate_wavefunctions_descriptors(op%wfd_overlap(i1,i2), subname)
          end do
      end do
      deallocate(op%wfd_overlap)
      nullify(op%wfd_overlap)
  end if



end subroutine deallocate_overlapParameters

!!subroutine deallocate_p2pCommsOrthonormalityMatrix(comom, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  use module_interfaces, exceptThisOne => deallocate_p2pCommsOrthonormalityMatrix
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
!!  character(len=*),intent(in):: subname
!!
!!  ! Local variables
!!  integer:: iis1, iis2, iie1, iie2, i1, i2
!!
!!  call checkAndDeallocatePointer(comom%noverlap, 'comom%noverlap', subname)
!!  !!call checkAndDeallocatePointer(comom%noverlapProc, 'comom%noverlapProc', subname)
!!  call checkAndDeallocatePointer(comom%overlaps, 'comom%overlaps', subname)
!!  !!call checkAndDeallocatePointer(comom%indexInRecvBuf, 'comom%indexInRecvBuf', subname)
!!  !!call checkAndDeallocatePointer(comom%overlapsProc, 'comom%overlapsProc', subname)
!!  call checkAndDeallocatePointer(comom%comarr, 'comom%comarr', subname)
!!  !call checkAndDeallocatePointer(comom%olrForExpansion, 'comom%olrForExpansion', subname)
!!  call checkAndDeallocatePointer(comom%recvBuf, 'comom%recvBuf', subname)
!!  call checkAndDeallocatePointer(comom%sendBuf, 'comom%sendBuf', subname)
!!  !!call checkAndDeallocatePointer(comom%communComplete, 'comom%communComplete', subname)
!!  call checkAndDeallocatePointer(comom%requests, 'comom%requests', subname)
!!
!!  !if(associated(comom%olr)) then
!!  !    iis1=lbound(comom%olr,1)
!!  !    iie1=ubound(comom%olr,1)
!!  !    iis2=lbound(comom%olr,2)
!!  !    iie2=ubound(comom%olr,2)
!!  !    do i2=iis2,iie2
!!  !        do i1=iis1,iie1
!!  !            call deallocate_matrixLocalizationRegion(comom%olr(i1,i2), subname)
!!  !        end do
!!  !    end do
!!  !    deallocate(comom%olr)
!!  !    nullify(comom%olr)
!!  !end if
!!
!!end subroutine deallocate_p2pCommsOrthonormalityMatrix


subroutine deallocate_matrixDescriptors(mad, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_matrixDescriptors
  implicit none
  
  ! Calling arguments
  type(matrixDescriptors),intent(inout):: mad
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(mad%keyg, 'mad%keyg', subname)
  call checkAndDeallocatePointer(mad%keyv, 'mad%keyv', subname)
  call checkAndDeallocatePointer(mad%keygmatmul, 'mad%keygmatmul', subname)
  call checkAndDeallocatePointer(mad%keyvmatmul, 'mad%keyvmatmul', subname)
  call checkAndDeallocatePointer(mad%nsegline, 'mad%nsegline', subname)
  call checkAndDeallocatePointer(mad%keygline, 'mad%keygline', subname)

end subroutine deallocate_matrixDescriptors



subroutine deallocate_collectiveComms(collComms, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_collectiveComms
  implicit none
  
  ! Calling arguments
  type(collectiveComms),intent(inout):: collComms
  character(len=*),intent(in):: subname

  ! Local variables

  call checkAndDeallocatePointer(collComms%nvctr_par, 'collComms%nvctr_par', subname)
  call checkAndDeallocatePointer(collComms%sendcnts, 'collComms%sendcnts', subname)
  call checkAndDeallocatePointer(collComms%senddspls, 'collComms%senddspls', subname)
  call checkAndDeallocatePointer(collComms%recvcnts, 'collComms%recvcnts', subname)
  call checkAndDeallocatePointer(collComms%recvdspls, 'collComms%recvdspls', subname)
  call checkAndDeallocatePointer(collComms%indexarray, 'collComms%indexarray', subname)

end subroutine deallocate_collectiveComms



subroutine destroy_wfn_metadata(wfnmd)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(wfn_metadata),intent(inout):: wfnmd

  ! Local variables
  integer:: istat, iall
  character(len=*),parameter:: subname='destroy_wfn_metadata'

  iall=-product(shape(wfnmd%coeff))*kind(wfnmd%coeff)
  deallocate(wfnmd%coeff, stat=istat)
  call memocc(istat, iall, 'wfnmd%coeff', subname)

  iall=-product(shape(wfnmd%coeff_proj))*kind(wfnmd%coeff_proj)
  deallocate(wfnmd%coeff_proj, stat=istat)
  call memocc(istat, iall, 'wfnmd%coeff_proj', subname)

  iall=-product(shape(wfnmd%coeffp))*kind(wfnmd%coeffp)
  deallocate(wfnmd%coeffp, stat=istat)
  call memocc(istat, iall, 'wfnmd%coeffp', subname)

  iall=-product(shape(wfnmd%density_kernel))*kind(wfnmd%density_kernel)
  deallocate(wfnmd%density_kernel, stat=istat)
  call memocc(istat, iall, 'wfnmd%density_kernel', subname)

  iall=-product(shape(wfnmd%alpha_coeff))*kind(wfnmd%alpha_coeff)
  deallocate(wfnmd%alpha_coeff, stat=istat)
  call memocc(istat, iall, 'wfnmd%alpha_coeff', subname)

  iall=-product(shape(wfnmd%grad_coeff_old))*kind(wfnmd%grad_coeff_old)
  deallocate(wfnmd%grad_coeff_old, stat=istat)
  call memocc(istat, iall, 'wfnmd%grad_coeff_old', subname)

end subroutine destroy_wfn_metadata




subroutine deallocate_collective_comms(collcom, subname)
  use module_base
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

end subroutine deallocate_collective_comms


subroutine deallocate_overlap_parameters_matrix(opm, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, except_this_one => deallocate_overlap_parameters_matrix
  implicit none
  
  ! Calling arguments
  type(overlap_parameters_matrix),intent(inout):: opm
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: i1, i2, iis1, iie1, iis2, iie2

  call checkAndDeallocatePointer(opm%noverlap, 'opm%noverlap', subname)
  call checkAndDeallocatePointer(opm%overlaps, 'opm%overlaps', subname)
  call checkAndDeallocatePointer(opm%olrForExpansion, 'opm%olrForExpansion', subname)

  if(associated(opm%olr)) then
      iis1=lbound(opm%olr,1)
      iie1=ubound(opm%olr,1)
      iis2=lbound(opm%olr,2)
      iie2=ubound(opm%olr,2)
      do i2=iis2,iie2
          do i1=iis1,iie1
              call deallocate_matrixLocalizationRegion(opm%olr(i1,i2), subname)
          end do
      end do
      deallocate(opm%olr)
      nullify(opm%olr)
  end if

end subroutine deallocate_overlap_parameters_matrix
