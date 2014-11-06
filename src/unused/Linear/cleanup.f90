!!subroutine deallocate_expansionSegments(expseg, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(expansionSegments),intent(inout):: expseg
!!  character(len=*),intent(in):: subname
!!  
!!  call checkAndDeallocatePointer(expseg%segborders, 'expseg%segborders', subname)
!!
!!end subroutine deallocate_expansionSegments

!!subroutine deallocate_p2pCommsSumrho(comsr, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(p2pCommsSumrho),intent(inout):: comsr
!!  character(len=*),intent(in):: subname
!!
!!  call checkAndDeallocatePointer(comsr%noverlaps, 'comsr%noverlaps', subname)
!!  call checkAndDeallocatePointer(comsr%overlaps, 'comsr%overlaps', subname)
!!  call checkAndDeallocatePointer(comsr%istarr, 'comsr%istarr', subname)
!!  call checkAndDeallocatePointer(comsr%istrarr, 'comsr%istrarr', subname)
!!  call checkAndDeallocatePointer(comsr%sendBuf, 'comsr%sendBuf', subname)
!!  call checkAndDeallocatePointer(comsr%recvBuf, 'comsr%recvBuf', subname)
!!  call checkAndDeallocatePointer(comsr%comarr, 'comsr%comarr', subname)
!!  call checkAndDeallocatePointer(comsr%communComplete, 'comsr%communComplete', subname)
!!  call checkAndDeallocatePointer(comsr%computComplete, 'comsr%computComplete', subname)
!!  call checkAndDeallocatePointer(comsr%auxarray, 'comsr%auxarray', subname)
!!  call checkAndDeallocatePointer(comsr%startingindex, 'comsr%startingindex', subname)
!!
!!end subroutine deallocate_p2pCommsSumrho

!!subroutine deallocate_p2pCommsGatherPot(comgp, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(p2pCommsGatherPot),intent(inout):: comgp
!!  character(len=*),intent(in):: subname
!!
!!  call checkAndDeallocatePointer(comgp%noverlaps, 'comgp%noverlaps', subname)
!!  call checkAndDeallocatePointer(comgp%overlaps, 'comgp%overlaps', subname)
!!  call checkAndDeallocatePointer(comgp%ise3, 'comgp%ise3', subname)
!!  call checkAndDeallocatePointer(comgp%comarr, 'comgp%comarr', subname)
!!  call checkAndDeallocatePointer(comgp%recvBuf, 'comgp%recvBuf', subname)
!!  call checkAndDeallocatePointer(comgp%communComplete, 'comgp%communComplete', subname)
!!
!!end subroutine deallocate_p2pCommsGatherPot

!!subroutine deallocate_p2pCommsRepartition(comrp, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  implicit none
!!
!!  ! Calling arguments
!!  type(p2pCommsRepartition),intent(inout):: comrp
!!  character(len=*),intent(in):: subname
!!
!!  call checkAndDeallocatePointer(comrp%comarr, 'comrp%comarr', subname)
!!  call checkAndDeallocatePointer(comrp%communComplete, 'comrp%communComplete', subname)
!!  call checkAndDeallocatePointer(comrp%requests, 'comrp%requests', subname)
!!
!!end subroutine deallocate_p2pCommsRepartition



!!subroutine deallocate_p2pCommsOrthonormality(comon, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(p2pCommsOrthonormality),intent(inout):: comon
!!  character(len=*),intent(in):: subname
!!
!!  call checkAndDeallocatePointer(comon%noverlaps, 'comon%noverlaps', subname)
!!  !call checkAndDeallocatePointer(comon%overlaps, 'comon%overlaps', subname)
!!  call checkAndDeallocatePointer(comon%comarr, 'comon%comarr', subname)
!!  call checkAndDeallocatePointer(comon%sendBuf, 'comon%sendBuf', subname)
!!  call checkAndDeallocatePointer(comon%recvBuf, 'comon%recvBuf', subname)
!!  call checkAndDeallocatePointer(comon%communComplete, 'comon%communComplete', subname)
!!  call checkAndDeallocatePointer(comon%requests, 'comon%requests', subname)
!!
!!end subroutine deallocate_p2pCommsOrthonormality

!!!subroutine deallocate_linearParameters(lin, subname)
!!!  use module_base
!!!  use module_types
!!!  use deallocatePointers
!!!  use module_interfaces, exceptThisOne => deallocate_linearParameters
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  type(linearParameters),intent(inout):: lin
!!!  character(len=*),intent(in):: subname
!!!
!!!  integer:: ierr
!!!
!!!  call checkAndDeallocatePointer(lin%potentialPrefac, 'lin%potentialPrefac', subname)
!!!  call checkAndDeallocatePointer(lin%potentialPrefac_lowaccuracy, 'lin%potentialPrefac_lowaccuracy', subname)
!!!  call checkAndDeallocatePointer(lin%potentialPrefac_highaccuracy, 'lin%potentialPrefac_highaccuracy', subname)
!!!  call checkAndDeallocatePointer(lin%locrad, 'lin%locrad', subname)
!!!  !call checkAndDeallocatePointer(lin%lphiRestart, 'lin%lphiRestart', subname)
!!!  call checkAndDeallocatePointer(lin%locrad_lowaccuracy, 'lin%locrad_lowaccuracy', subname)
!!!  call checkAndDeallocatePointer(lin%locrad_highaccuracy, 'lin%locrad_highaccuracy', subname)
!!!  !call checkAndDeallocatePointer(lin%lphiold, 'lin%lphiold', subname)
!!!  !call checkAndDeallocatePointer(lin%hamold, 'lin%hamold', subname)
!!!  !call checkAndDeallocatePointer(lin%lphiold, 'lphiold', subname)
!!!  !call checkAndDeallocatePointer(lin%lhphiold, 'lhphiold', subname)
!!!  !call checkAndDeallocatePointer(lin%hamold, 'lin%hamold', subname)
!!!  call deallocate_orbitals_data(lin%orbs, subname)
!!!  call deallocate_orbitals_data(lin%gorbs, subname)
!!!  call deallocate_comms_cubic(lin%comms, subname)
!!!  call deallocate_comms_cubic(lin%gcomms, subname)
!!!  call checkAndDeallocatePointer(lin%norbsPerType, 'lin%norbsPerType', subname)
!!!  !call deallocate_p2pCommsSumrho(lin%comsr, subname)
!!!  call deallocate_p2pComms(lin%comsr, subname)
!!!  !call deallocate_p2pCommsGatherPot(lin%comgp, subname)
!!!  call deallocate_p2pComms(lin%comgp, subname)
!!!  call deallocate_largeBasis(lin%lb, subname)
!!!!  call deallocate_nonlocal_psp_descriptors(lin%lzd%Gnlpspd, subname)
!!!  call deallocate_local_zone_descriptors(lin%lzd, subname)
!!!  !call deallocate_p2pCommsOrthonormality(lin%comon, subname)
!!!  call deallocate_p2pComms(lin%comon, subname)
!!!  call deallocate_overlapParameters(lin%op, subname)
!!!  call deallocate_matrixDescriptors(lin%mad, subname)
!!!  call deallocate_collectiveComms(lin%collComms, subname)
!!!
!!!end subroutine deallocate_linearParameters

!!subroutine deallocate_inguessParameters(ip, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  use module_interfaces, exceptThisOne => deallocate_inguessParameters
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(inguessParameters),intent(inout):: ip
!!  character(len=*),intent(in):: subname
!!
!!  ! Local variables
!!  integer:: iis1, iie1, i1
!!
!!  !!call checkAndDeallocatePointer(ip%norb_par, 'ip%norb_par', subname)
!!  !!call checkAndDeallocatePointer(ip%onWhichMPI, 'ip%onWhichMPI', subname)
!!  !!call checkAndDeallocatePointer(ip%isorb_par, 'ip%isorb_par', subname)
!!  !!call checkAndDeallocatePointer(ip%nvctrp_nz, 'ip%nvctrp_nz', subname)
!!  call checkAndDeallocatePointer(ip%sendcounts, 'ip%sendcounts', subname)
!!  call checkAndDeallocatePointer(ip%senddispls, 'ip%senddispls', subname)
!!  call checkAndDeallocatePointer(ip%recvcounts, 'ip%recvcounts',  subname)
!!  call checkAndDeallocatePointer(ip%recvdispls, 'ip%recvdispls', subname)
!!
!!end subroutine deallocate_inguessParameters



!!subroutine deallocate_collectiveComms(collComms, subname)
!!  use module_base
!!  use module_types
!!  use deallocatePointers
!!  use module_interfaces, exceptThisOne => deallocate_collectiveComms
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(collectiveComms),intent(inout):: collComms
!!  character(len=*),intent(in):: subname
!!
!!  ! Local variables
!!
!!  call checkAndDeallocatePointer(collComms%nvctr_par, 'collComms%nvctr_par', subname)
!!  call checkAndDeallocatePointer(collComms%sendcnts, 'collComms%sendcnts', subname)
!!  call checkAndDeallocatePointer(collComms%senddspls, 'collComms%senddspls', subname)
!!  call checkAndDeallocatePointer(collComms%recvcnts, 'collComms%recvcnts', subname)
!!  call checkAndDeallocatePointer(collComms%recvdspls, 'collComms%recvdspls', subname)
!!  call checkAndDeallocatePointer(collComms%indexarray, 'collComms%indexarray', subname)
!!
!!end subroutine deallocate_collectiveComms


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

end subroutine deallocate_nonlocal_psp_descriptors



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
