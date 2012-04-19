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


