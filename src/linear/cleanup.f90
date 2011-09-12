module deallocatePointers

  interface checkAndDeallocatePointer
     module procedure checkAndDeallocatePointer_int_1, checkAndDeallocatePointer_sgl_1, checkAndDeallocatePointer_dbl_1,&
                      checkAndDeallocatePointer_log_1
     module procedure checkAndDeallocatePointer_int_2, checkAndDeallocatePointer_sgl_2, checkAndDeallocatePointer_dbl_2,&
                      checkAndDeallocatePointer_log_2
     module procedure checkAndDeallocatePointer_int_3, checkAndDeallocatePointer_sgl_3, checkAndDeallocatePointer_dbl_3,&
                      checkAndDeallocatePointer_log_3
     module procedure checkAndDeallocatePointer_int_4, checkAndDeallocatePointer_sgl_4, checkAndDeallocatePointer_dbl_4,&
                      checkAndDeallocatePointer_log_4
     module procedure checkAndDeallocatePointer_int_5, checkAndDeallocatePointer_sgl_5, checkAndDeallocatePointer_dbl_5,&
                      checkAndDeallocatePointer_log_5
     module procedure checkAndDeallocatePointer_int_6, checkAndDeallocatePointer_sgl_6, checkAndDeallocatePointer_dbl_6,&
                      checkAndDeallocatePointer_log_6
     module procedure checkAndDeallocatePointer_int_7, checkAndDeallocatePointer_sgl_7, checkAndDeallocatePointer_dbl_7,&
                      checkAndDeallocatePointer_log_7
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


!!module deallocationInterfaces
!!  implicit none
!!
!!  interface
!!
!!
!!    subroutine deallocate_local_zone_descriptors(lzd, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(local_zone_descriptors),intent(inout):: lzd
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_local_zone_descriptors
!!    
!!    subroutine deallocate_orbitals_data(orbs, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(orbitals_data),intent(inout):: orbs
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_orbitals_data
!!    
!!    subroutine deallocate_communications_arrays(comms, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(communications_arrays),intent(inout):: comms
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_communications_arrays
!!    
!!    subroutine deallocate_locreg_descriptors(lr, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(locreg_descriptors),intent(inout):: lr
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_locreg_descriptors
!!    
!!    subroutine deallocate_wavefunctions_descriptors(wfd, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(wavefunctions_descriptors),intent(inout):: wfd
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_wavefunctions_descriptors
!!    
!!    subroutine deallocate_convolutions_bounds(bounds, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(convolutions_bounds),intent(inout):: bounds
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_convolutions_bounds
!!    
!!    subroutine deallocate_kinetic_bounds(kb, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(kinetic_bounds),intent(inout):: kb
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_kinetic_bounds
!!    
!!    subroutine deallocate_shrink_bounds(sb, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(shrink_bounds),intent(inout):: sb
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_shrink_bounds
!!    
!!    subroutine deallocate_grow_bounds(gb, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(grow_bounds),intent(inout):: gb
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_grow_bounds
!!    
!!    subroutine deallocate_nonlocal_psp_descriptors(nlpspd, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(nonlocal_psp_descriptors),intent(inout):: nlpspd
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_nonlocal_psp_descriptors
!!    
!!    subroutine deallocate_matrixMinimization(matmin, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(matrixMinimization),intent(inout):: matmin
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_matrixMinimization
!!    
!!    subroutine deallocate_matrixLocalizationRegion(mlr, subname)
!!      use module_base
!!      use module_types
!!      use deallocatePointers
!!      implicit none
!!      type(matrixLocalizationRegion),intent(inout):: mlr
!!      character(len=*),intent(in):: subname
!!    end subroutine deallocate_matrixLocalizationRegion
!!
!!
!!  end interface
!!
!!end module deallocationInterfaces



subroutine deallocate_linearParameters(lin, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_linearParameters
  implicit none

  ! Calling arguments
  type(linearParameters),intent(inout):: lin
  character(len=*),intent(in):: subname

  integer:: ierr

  call checkAndDeallocatePointer(lin%potentialPrefac, 'lin%potentialPrefac', subname)
  call checkAndDeallocatePointer(lin%locrad, 'lin%locrad', subname)
  call checkAndDeallocatePointer(lin%lphiRestart, 'lphi%Restart', subname)
  !call checkAndDeallocatePointer(lin%lphiold, 'lphiold', subname)
  !call checkAndDeallocatePointer(lin%lhphiold, 'lhphiold', subname)
  !call checkAndDeallocatePointer(lin%hamold, 'lin%hamold', subname)
  call deallocate_orbitals_data(lin%orbs, subname)
  call deallocate_orbitals_data(lin%gorbs, subname)
  call deallocate_communications_arrays(lin%comms, subname)
  call deallocate_communications_arrays(lin%gcomms, subname)
  call checkAndDeallocatePointer(lin%norbsPerType, 'lin%norbsPerType', subname)
  call deallocate_p2pCommsSumrho(lin%comsr, subname)
  call deallocate_p2pCommsGatherPot(lin%comgp, subname)
  call deallocate_largeBasis(lin%lb, subname)
  call deallocate_local_zone_descriptors(lin%lzd, subname)
  call deallocate_p2pCommsOrthonormality(lin%comon, subname)
  call deallocate_overlapParameters(lin%op, subname)
  call deallocate_matrixDescriptors(lin%mad, subname)

end subroutine deallocate_linearParameters



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
  integer:: istat, iall, iis1, iie1, i1
  
  
!  call deallocate_orbitals_data(lzd%orbs, subname)
  
!  iis1=lbound(lzd%lorbs,1)
!  iie1=ubound(lzd%lorbs,1)
!  do i1=iis1,iie1
!      call deallocate_orbitals_data(lzd%lorbs(i1), subname)
!  end do
  
!  call deallocate_communications_arrays(lzd%comms, subname)
  
  call checkAndDeallocatePointer(lzd%Glr%projflg, 'lzd%Glr%projflg', subname)
  call checkAndDeallocatePointer(lzd%doHamAppl, 'lzd%doHamAppl', subname)
  call deallocate_locreg_descriptors(lzd%Glr, subname)

  call deallocate_nonlocal_psp_descriptors(lzd%Gnlpspd, subname)

  if(associated(lzd%llr)) then  
     iis1=lbound(lzd%llr,1)
     iie1=ubound(lzd%llr,1)
     !write(*,*) 'iis1,iie1',iis1,iie1
     do i1=iis1,iie1
         !if(associated(lzd%llr(i1)%projflg)) then
         !    nullify(lzd%llr(i1)%projflg)
         !end if
         call checkAndDeallocatePointer(lzd%llr(i1)%projflg, 'lzd%llr(i1)%projflg', subname)
         !write(*,*) 'i1',i1
         call deallocate_locreg_descriptors(lzd%llr(i1), subname)
     end do
  end if

  if(associated(lzd%lnlpspd)) then 
     iis1=lbound(lzd%lnlpspd,1)
     iie1=ubound(lzd%lnlpspd,1)
     do i1=iis1,iie1
         call deallocate_nonlocal_psp_descriptors(lzd%lnlpspd(i1), subname)
     end do
  end if
end subroutine deallocate_local_zone_descriptors


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
  call checkAndDeallocatePointer(orbs%inWhichLocregP, 'orbs%inWhichLocreg', subname)
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

  call checkAndDeallocatePointer(lr%projflg, 'lr%projflg', subname)

  call deallocate_wavefunctions_descriptors(lr%wfd, subname)
  call deallocate_convolutions_bounds(lr%bounds, subname)
  
  
end subroutine deallocate_locreg_descriptors

subroutine deallocate_locreg_descriptors2(lr, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_locreg_descriptors
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(inout):: lr
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(lr%projflg, 'lr%projflg', subname)

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

  call checkAndDeallocatePointer(wfd%keyg, 'wfd%keyg', subname)
  call checkAndDeallocatePointer(wfd%keyv, 'wfd%keyv', subname)

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

  call checkAndDeallocatePointer(nlpspd%nvctr_p, 'nlpspd%nvctr_p', subname)
  call checkAndDeallocatePointer(nlpspd%nseg_p, 'nlpspd%nseg_p', subname)
  call checkAndDeallocatePointer(nlpspd%keyv_p, 'nlpspd%keyv_p', subname)
  call checkAndDeallocatePointer(nlpspd%keyg_p, 'nlpspd%keyg_p', subname)
  call checkAndDeallocatePointer(nlpspd%nboxp_c, 'nlpspd%nboxp_c', subname)
  call checkAndDeallocatePointer(nlpspd%nboxp_f, 'nlpspd%nboxp_f', subname)

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


subroutine deallocate_p2pCommsSumrho(comsr, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(p2pCommsSumrho),intent(inout):: comsr
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(comsr%noverlaps, 'comsr%noverlaps', subname)
  call checkAndDeallocatePointer(comsr%overlaps, 'comsr%overlaps', subname)
  call checkAndDeallocatePointer(comsr%istarr, 'comsr%istarr', subname)
  call checkAndDeallocatePointer(comsr%istrarr, 'comsr%istrarr', subname)
  call checkAndDeallocatePointer(comsr%sendBuf, 'comsr%sendBuf', subname)
  call checkAndDeallocatePointer(comsr%recvBuf, 'comsr%recvBuf', subname)
  call checkAndDeallocatePointer(comsr%comarr, 'comsr%comarr', subname)
  call checkAndDeallocatePointer(comsr%communComplete, 'comsr%communComplete', subname)
  call checkAndDeallocatePointer(comsr%computComplete, 'comsr%computComplete', subname)

end subroutine deallocate_p2pCommsSumrho


subroutine deallocate_p2pCommsGatherPot(comgp, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(p2pCommsGatherPot),intent(inout):: comgp
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(comgp%noverlaps, 'comgp%noverlaps', subname)
  call checkAndDeallocatePointer(comgp%overlaps, 'comgp%overlaps', subname)
  call checkAndDeallocatePointer(comgp%ise3, 'comgp%ise3', subname)
  call checkAndDeallocatePointer(comgp%comarr, 'comgp%comarr', subname)
  call checkAndDeallocatePointer(comgp%recvBuf, 'comgp%recvBuf', subname)
  call checkAndDeallocatePointer(comgp%communComplete, 'comgp%communComplete', subname)

end subroutine deallocate_p2pCommsGatherPot


subroutine deallocate_largeBasis(lb, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_largeBasis
  implicit none
  
  ! Calling arguments
  type(largeBasis),intent(inout):: lb
  character(len=*),intent(in):: subname

  call deallocate_communications_arrays(lb%comms, subname)
  call deallocate_communications_arrays(lb%gcomms, subname)
  call deallocate_orbitals_data(lb%orbs, subname)
  call deallocate_orbitals_data(lb%gorbs, subname)
  !call deallocate_local_zone_descriptors(lb%lzd, subname)
  call dealloctae_p2pCommsRepartition(lb%comrp, subname)
  call deallocate_p2pCommsOrthonormality(lb%comon, subname)
  call deallocate_overlapParameters(lb%op, subname)
  call deallocate_p2pCommsGatherPot(lb%comgp, subname)


end subroutine deallocate_largeBasis


subroutine dealloctae_p2pCommsRepartition(comrp, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(p2pCommsRepartition),intent(inout):: comrp
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(comrp%comarr, 'comrp%comarr', subname)
  call checkAndDeallocatePointer(comrp%communComplete, 'comrp%communComplete', subname)

end subroutine dealloctae_p2pCommsRepartition



subroutine deallocate_p2pCommsOrthonormality(comon, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(p2pCommsOrthonormality),intent(inout):: comon
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(comon%noverlaps, 'comon%noverlaps', subname)
  call checkAndDeallocatePointer(comon%overlaps, 'comon%overlaps', subname)
  call checkAndDeallocatePointer(comon%comarr, 'comon%comarr', subname)
  call checkAndDeallocatePointer(comon%sendBuf, 'comon%sendBuf', subname)
  call checkAndDeallocatePointer(comon%recvBuf, 'comon%recvBuf', subname)
  call checkAndDeallocatePointer(comon%communComplete, 'comon%communComplete', subname)

end subroutine deallocate_p2pCommsOrthonormality


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
  call checkAndDeallocatePointer(op%indexExpand, 'op%indexExpand', subname)
  call checkAndDeallocatePointer(op%indexExtract, 'op%indexExtract', subname)
  call checkAndDeallocatePointer(op%overlaps, 'op%overlaps', subname)
  call checkAndDeallocatePointer(op%indexInRecvBuf, 'op%indexInRecvBuf', subname)
  call checkAndDeallocatePointer(op%indexInSendBuf, 'op%indexInSendBuf', subname)


  iis1=lbound(op%olr,1)
  iie1=ubound(op%olr,1)
  iis2=lbound(op%olr,2)
  iie2=ubound(op%olr,2)
  do i2=iis2,iie2
      do i1=iis1,iie1
          call deallocate_locreg_descriptors2(op%olr(i1,i2), subname)
      end do
  end do

end subroutine deallocate_overlapParameters



subroutine deallocate_inguessParameters(ip, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_inguessParameters
  implicit none
  
  ! Calling arguments
  type(inguessParameters),intent(inout):: ip
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: iis1, iie1, i1

  call checkAndDeallocatePointer(ip%norb_par, 'ip%norb_par', subname)
  call checkAndDeallocatePointer(ip%onWhichMPI, 'ip%onWhichMPI', subname)
  call checkAndDeallocatePointer(ip%isorb_par, 'ip%isorb_par', subname)
  call checkAndDeallocatePointer(ip%nvctrp_nz, 'ip%nvctrp_nz', subname)
  call checkAndDeallocatePointer(ip%sendcounts, 'ip%sendcounts', subname)
  call checkAndDeallocatePointer(ip%senddispls, 'ip%senddispls', subname)
  call checkAndDeallocatePointer(ip%recvcounts, 'ip%recvcounts',  subname)
  call checkAndDeallocatePointer(ip%recvdispls, 'ip%recvdispls', subname)

end subroutine deallocate_inguessParameters


subroutine deallocate_p2pCommsOrthonormalityMatrix(comom, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_p2pCommsOrthonormalityMatrix
  implicit none
  
  ! Calling arguments
  type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: iis1, iis2, iie1, iie2, i1, i2

  call checkAndDeallocatePointer(comom%noverlap, 'comom%noverlap', subname)
  call checkAndDeallocatePointer(comom%noverlapProc, 'comom%noverlapProc', subname)
  call checkAndDeallocatePointer(comom%overlaps, 'comom%overlaps', subname)
  call checkAndDeallocatePointer(comom%indexInRecvBuf, 'comom%indexInRecvBuf', subname)
  call checkAndDeallocatePointer(comom%overlapsProc, 'comom%overlapsProc', subname)
  call checkAndDeallocatePointer(comom%comarr, 'comom%comarr', subname)
  call checkAndDeallocatePointer(comom%olrForExpansion, 'comom%olrForExpansion', subname)
  call checkAndDeallocatePointer(comom%recvBuf, 'comom%recvBuf', subname)
  call checkAndDeallocatePointer(comom%sendBuf, 'comom%sendBuf', subname)
  call checkAndDeallocatePointer(comom%communComplete, 'comom%communComplete', subname)

  iis1=lbound(comom%olr,1)
  iie1=ubound(comom%olr,1)
  iis2=lbound(comom%olr,2)
  iie2=ubound(comom%olr,2)
  do i2=iis2,iie2
      do i1=iis1,iie1
          call deallocate_matrixLocalizationRegion(comom%olr(i1,i2), subname)
      end do
  end do

end subroutine deallocate_p2pCommsOrthonormalityMatrix


subroutine deallocate_matrixDescriptors(mad, subname)
  use module_base
  use module_types
  use deallocatePointers
  use module_interfaces, exceptThisOne => deallocate_matrixDescriptors
  implicit none
  
  ! Calling arguments
  type(matrixDescriptors),intent(inout):: mad
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: iis1, iis2, iie1, iie2, i1, i2

  call checkAndDeallocatePointer(mad%keyg, 'mad%keyg', subname)
  call checkAndDeallocatePointer(mad%keyv, 'mad%keyv', subname)
  call checkAndDeallocatePointer(mad%keygmatmul, 'mad%keygmatmul', subname)
  call checkAndDeallocatePointer(mad%keyvmatmul, 'mad%keyvmatmul', subname)
  call checkAndDeallocatePointer(mad%nsegline, 'mad%nsegline', subname)
  call checkAndDeallocatePointer(mad%keygline, 'mad%keygline', subname)

end subroutine deallocate_matrixDescriptors
