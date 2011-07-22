module deallocatePointers

  interface checkAndDeallocatePointer
     module procedure checkAndDeallocatePointer_int_1, checkAndDeallocatePointer_sgl_1, checkAndDeallocatePointer_dbl_1, checkAndDeallocatePointer_log_1
     module procedure checkAndDeallocatePointer_int_2, checkAndDeallocatePointer_sgl_2, checkAndDeallocatePointer_dbl_2, checkAndDeallocatePointer_log_2
     module procedure checkAndDeallocatePointer_int_3, checkAndDeallocatePointer_sgl_3, checkAndDeallocatePointer_dbl_3, checkAndDeallocatePointer_log_3
     module procedure checkAndDeallocatePointer_int_4, checkAndDeallocatePointer_sgl_4, checkAndDeallocatePointer_dbl_4, checkAndDeallocatePointer_log_4
     module procedure checkAndDeallocatePointer_int_5, checkAndDeallocatePointer_sgl_5, checkAndDeallocatePointer_dbl_5, checkAndDeallocatePointer_log_5
     module procedure checkAndDeallocatePointer_int_6, checkAndDeallocatePointer_sgl_6, checkAndDeallocatePointer_dbl_6, checkAndDeallocatePointer_log_6
     module procedure checkAndDeallocatePointer_int_7, checkAndDeallocatePointer_sgl_7, checkAndDeallocatePointer_dbl_7, checkAndDeallocatePointer_log_7
  end interface

  contains


    subroutine checkAndDeallocatePointer_int_1(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      integer,dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatepointer_int_1


    subroutine checkAndDeallocatePointer_sgl_1(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(4),dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_sgl_1


    subroutine checkAndDeallocatePointer_dbl_1(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(8),dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_dbl_1


    subroutine checkAndDeallocatePointer_log_1(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      logical,dimension(:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_log_1


    subroutine checkAndDeallocatePointer_int_2(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      integer,dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatepointer_int_2


    subroutine checkAndDeallocatePointer_sgl_2(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(4),dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_sgl_2


    subroutine checkAndDeallocatePointer_dbl_2(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(8),dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_dbl_2


    subroutine checkAndDeallocatePointer_log_2(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      logical,dimension(:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_log_2


    subroutine checkAndDeallocatePointer_int_3(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatepointer_int_3


    subroutine checkAndDeallocatePointer_sgl_3(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_sgl_3


    subroutine checkAndDeallocatePointer_dbl_3(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_dbl_3


    subroutine checkAndDeallocatePointer_log_3(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_log_3


    subroutine checkAndDeallocatePointer_int_4(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatepointer_int_4


    subroutine checkAndDeallocatePointer_sgl_4(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_sgl_4


    subroutine checkAndDeallocatePointer_dbl_4(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_dbl_4


    subroutine checkAndDeallocatePointer_log_4(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_log_4



    subroutine checkAndDeallocatePointer_int_5(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatepointer_int_5


    subroutine checkAndDeallocatePointer_sgl_5(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_sgl_5


    subroutine checkAndDeallocatePointer_dbl_5(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_dbl_5


    subroutine checkAndDeallocatePointer_log_5(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_log_5


    subroutine checkAndDeallocatePointer_int_6(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatepointer_int_6


    subroutine checkAndDeallocatePointer_sgl_6(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_sgl_6


    subroutine checkAndDeallocatePointer_dbl_6(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_dbl_6


    subroutine checkAndDeallocatePointer_log_6(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_log_6


    subroutine checkAndDeallocatePointer_int_7(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      integer,dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatepointer_int_7


    subroutine checkAndDeallocatePointer_sgl_7(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(4),dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_sgl_7


    subroutine checkAndDeallocatePointer_dbl_7(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      real(8),dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_dbl_7


    subroutine checkAndDeallocatePointer_log_7(array, subname)
      use module_base
      implicit none

      ! Calling arguments
      logical,dimension(:,:,:,:,:,:,:),pointer,intent(inout):: array
      character(len=*),intent(in):: subname
      
      ! Local variables
      integer:: istat, iall

      if(associated(array)) then
          iall=-product(shape(array))*kind(array)
          deallocate(array, stat=istat)
          !call memocc(istat, iall, 'array', subname)
      end if

    end subroutine checkAndDeallocatePointer_log_7


end module deallocatePointers


module deallocationInterfaces
  implicit none

  interface


    subroutine deallocate_linear_zone_descriptors(lzd, subname)
      use module_base
      use module_types
      use deallocatePointers
      implicit none
      type(linear_zone_descriptors),intent(inout):: lzd
      character(len=*),intent(in):: subname
    end subroutine deallocate_linear_zone_descriptors
    
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



subroutine deallocate_linear_zone_descriptors(lzd, subname)
  use module_base
  use module_types
  use deallocatePointers
  use deallocationInterfaces, exceptThisOne => deallocate_linear_zone_descriptors
  implicit none
  
  ! Calling arguments
  type(linear_zone_descriptors),intent(inout):: lzd
  character(len=*),intent(in):: subname
  
  ! Local variables
  integer:: istat, iall, iis1, iie1, i1
  
  
  call deallocate_orbitals_data(lzd%orbs, subname)
  
  iis1=lbound(lzd%lorbs,1)
  iie1=ubound(lzd%lorbs,1)
  do i1=iis1,iie1
      call deallocate_orbitals_data(lzd%lorbs(i1), subname)
  end do
  
  call deallocate_communications_arrays(lzd%comms, subname)
  
  call checkAndDeallocatePointer(lzd%Glr%projflg, subname)
  call deallocate_locreg_descriptors(lzd%Glr, subname)
  
  call deallocate_nonlocal_psp_descriptors(lzd%Gnlpspd, subname)
  
  iis1=lbound(lzd%llr,1)
  iie1=ubound(lzd%llr,1)
  do i1=iis1,iie1
      call checkAndDeallocatePointer(lzd%llr(i1)%projflg, subname)
      call deallocate_locreg_descriptors(lzd%llr(i1), subname)
  end do
  
  iis1=lbound(lzd%lnlpspd,1)
  iie1=ubound(lzd%lnlpspd,1)
  do i1=iis1,iie1
      call deallocate_nonlocal_psp_descriptors(lzd%lnlpspd(i1), subname)
  end do
  
end subroutine deallocate_linear_zone_descriptors


subroutine deallocate_orbitals_data(orbs, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(orbitals_data),intent(inout):: orbs
  character(len=*),intent(in):: subname
  
  call checkAndDeallocatePointer(orbs%norb_par, subname)
  call checkAndDeallocatePointer(orbs%iokpt, subname)
  call checkAndDeallocatePointer(orbs%ikptproc, subname)
  call checkAndDeallocatePointer(orbs%inwhichlocreg, subname)
  call checkAndDeallocatePointer(orbs%inWhichLocregP, subname)
  call checkAndDeallocatePointer(orbs%onWhichMPI, subname)
  call checkAndDeallocatePointer(orbs%isorb_par, subname)
  call checkAndDeallocatePointer(orbs%eval, subname)
  call checkAndDeallocatePointer(orbs%occup, subname)
  call checkAndDeallocatePointer(orbs%spinsgn, subname)
  call checkAndDeallocatePointer(orbs%kwgts, subname)
  call checkAndDeallocatePointer(orbs%kpts, subname)
  
end subroutine deallocate_orbitals_data


subroutine deallocate_communications_arrays(comms, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(communications_arrays),intent(inout):: comms
  character(len=*),intent(in):: subname
  
  call checkAndDeallocatePointer(comms%ncntd, subname)
  call checkAndDeallocatePointer(comms%ncntt, subname)
  call checkAndDeallocatePointer(comms%ndspld, subname)
  call checkAndDeallocatePointer(comms%ndsplt, subname)
  call checkAndDeallocatePointer(comms%nvctr_par, subname)
  
end subroutine deallocate_communications_arrays



subroutine deallocate_locreg_descriptors(lr, subname)
  use module_base
  use module_types
  use deallocatePointers
  use deallocationInterfaces, exceptThisOne => deallocate_locreg_descriptors
  implicit none
  
  ! Calling arguments
  type(locreg_descriptors),intent(inout):: lr
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(lr%projflg, subname)

  call deallocate_wavefunctions_descriptors(lr%wfd, subname)
  call deallocate_convolutions_bounds(lr%bounds, subname)
  
  
end subroutine deallocate_locreg_descriptors


subroutine deallocate_wavefunctions_descriptors(wfd, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(wavefunctions_descriptors),intent(inout):: wfd
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(wfd%keyg, subname)
  call checkAndDeallocatePointer(wfd%keyv, subname)

end subroutine deallocate_wavefunctions_descriptors


subroutine deallocate_convolutions_bounds(bounds, subname)
  use module_base
  use module_types
  use deallocatePointers
  use deallocationInterfaces, exceptThisOne => deallocate_convolutions_bounds
  implicit none
  
  ! Calling arguments
  type(convolutions_bounds),intent(inout):: bounds
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(bounds%ibyyzz_r, subname)

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

  call checkAndDeallocatePointer(kb%ibyz_c, subname)
  call checkAndDeallocatePointer(kb%ibxz_c, subname)
  call checkAndDeallocatePointer(kb%ibxy_c, subname)
  call checkAndDeallocatePointer(kb%ibyz_f, subname)
  call checkAndDeallocatePointer(kb%ibxz_f, subname)
  call checkAndDeallocatePointer(kb%ibxy_f, subname)

end subroutine deallocate_kinetic_bounds


subroutine deallocate_shrink_bounds(sb, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(shrink_bounds),intent(inout):: sb
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(sb%ibzzx_c, subname)
  call checkAndDeallocatePointer(sb%ibyyzz_c, subname)
  call checkAndDeallocatePointer(sb%ibxy_ff, subname)
  call checkAndDeallocatePointer(sb%ibzzx_f, subname)
  call checkAndDeallocatePointer(sb%ibyyzz_f, subname)

end subroutine deallocate_shrink_bounds


subroutine deallocate_grow_bounds(gb, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(grow_bounds),intent(inout):: gb
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(gb%ibzxx_c, subname)
  call checkAndDeallocatePointer(gb%ibxxyy_c, subname)
  call checkAndDeallocatePointer(gb%ibyz_ff, subname)
  call checkAndDeallocatePointer(gb%ibzxx_f, subname)
  call checkAndDeallocatePointer(gb%ibxxyy_f, subname)

end subroutine deallocate_grow_bounds


subroutine deallocate_nonlocal_psp_descriptors(nlpspd, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
 
  ! Calling arguments
  type(nonlocal_psp_descriptors),intent(inout):: nlpspd
  character(len=*),intent(in):: subname

  call checkAndDeallocatePointer(nlpspd%nvctr_p, subname)
  call checkAndDeallocatePointer(nlpspd%nseg_p, subname)
  call checkAndDeallocatePointer(nlpspd%keyv_p, subname)
  call checkAndDeallocatePointer(nlpspd%keyg_p, subname)
  call checkAndDeallocatePointer(nlpspd%nboxp_c, subname)
  call checkAndDeallocatePointer(nlpspd%nboxp_f, subname)

end subroutine deallocate_nonlocal_psp_descriptors



subroutine deallocate_matrixMinimization(matmin, subname)
  use module_base
  use module_types
  use deallocatePointers
  use deallocationInterfaces, exceptThisOne => deallocate_matrixMinimization
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
  
  call checkAndDeallocatePointer(matmin%inWhichLocregExtracted, subname)
  
  call checkAndDeallocatePointer(matmin%inWhichLocregOnMPI, subname)
  
  call checkAndDeallocatePointer(matmin%indexInLocreg, subname)

end subroutine deallocate_matrixMinimization



subroutine deallocate_matrixLocalizationRegion(mlr, subname)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(matrixLocalizationRegion),intent(inout):: mlr
  character(len=*),intent(in):: subname
  
  call checkAndDeallocatePointer(mlr%indexInGlobal, subname)
  
end subroutine deallocate_matrixLocalizationRegion



