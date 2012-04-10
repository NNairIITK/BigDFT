!!subroutine nullify_p2pCommsSumrho(comsr)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_p2pCommsSumrho
!!  implicit none
!!
!!  ! Calling argument
!!  type(p2pCommsSumrho),intent(out):: comsr
!!
!!  nullify(comsr%noverlaps)
!!  nullify(comsr%overlaps)
!!  nullify(comsr%istarr)
!!  nullify(comsr%istrarr)
!!  nullify(comsr%sendBuf)
!!  nullify(comsr%recvBuf)
!!  nullify(comsr%comarr)
!!  nullify(comsr%communComplete)
!!  nullify(comsr%computComplete)
!!  nullify(comsr%startingindex)
!!end subroutine nullify_p2pCommsSumrho
!!
!!
!!subroutine nullify_p2pCommsGatherPot(comgp)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_p2pCommsGatherPot
!!  implicit none
!!
!!  ! Calling argument
!!  type(p2pCommsGatherPot),intent(out):: comgp
!!
!!  nullify(comgp%noverlaps)
!!  nullify(comgp%overlaps)
!!  nullify(comgp%ise3)
!!  nullify(comgp%comarr)
!!  nullify(comgp%recvBuf)
!!  nullify(comgp%communComplete)
!!
!!end subroutine nullify_p2pCommsGatherPot

!!!subroutine nullify_p2pCommsRepartition(comrp)
!!!  use module_base
!!!  use module_types
!!!  use module_interfaces, exceptThisOne => nullify_p2pCommsRepartition
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  type(p2pCommsRepartition),intent(out):: comrp
!!!
!!!  nullify(comrp%comarr)
!!!  nullify(comrp%communComplete)
!!!  nullify(comrp%requests)
!!!
!!!end subroutine nullify_p2pCommsRepartition


!!subroutine nullify_p2pCommsOrthonormality(comon)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_p2pCommsOrthonormality
!!  implicit none
!!
!!  ! Calling argument
!!  type(p2pCommsOrthonormality),intent(out):: comon
!!
!!  nullify(comon%noverlaps)
!!  !!nullify(comon%overlaps)
!!  nullify(comon%comarr)
!!  nullify(comon%sendBuf)
!!  nullify(comon%recvBuf)
!!  nullify(comon%communComplete)
!!  nullify(comon%requests)
!!
!!end subroutine nullify_p2pCommsOrthonormality

!!subroutine nullify_expansionSegments(expseg)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_expansionSegments
!!  implicit none
!!
!!  ! Calling argument
!!  type(expansionSegments),intent(out):: expseg
!!
!!  nullify(expseg%segborders)
!!
!!end subroutine nullify_expansionSegments



!!subroutine nullify_linearInputGuess(lig)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => nullify_linearInputGuess
!!  implicit none
!!
!!  ! Calling argument
!!  type(linearInputGuess),intent(out):: lig
!!
!!  call nullify_local_zone_descriptors(lig%lzdig)
!!  call nullify_local_zone_descriptors(lig%lzdGauss)
!!  call nullify_orbitals_data(lig%orbsig)
!!  call nullify_orbitals_data(lig%orbsGauss)
!!  !call nullify_p2pCommsOrthonormality(lig%comon)
!!  call nullify_p2pComms(lig%comon)
!!  call nullify_overlapParameters(lig%op)
!!  !call nullify_p2pCommsGatherPot(lig%comgp)
!!  call nullify_p2pComms(lig%comgp)
!!  call nullify_matrixDescriptors(lig%mad)
!!
!!end subroutine nullify_linearInputGuess


