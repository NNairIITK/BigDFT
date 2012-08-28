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

!!!subroutine nullify_linearParameters(lin)
!!!  use module_base
!!!  use module_types
!!!  use module_interfaces, exceptThisOne => nullify_linearParameters
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  type(linearParameters),intent(out):: lin
!!!
!!!  nullify(lin%potentialPrefac)
!!!  nullify(lin%locrad)
!!!  !nullify(lin%lphiRestart)
!!!  !nullify(lin%lphiold)
!!!  !nullify(lin%lhphiold)
!!!  !nullify(lin%hamold)
!!!  call nullify_orbitals_data(lin%orbs)
!!!  call nullify_orbitals_data(lin%gorbs)
!!!  call nullify_communications_arrays(lin%comms)
!!!  call nullify_communications_arrays(lin%gcomms)
!!!  nullify(lin%norbsPerType)
!!!  !call nullify_p2pCommsSumrho(lin%comsr)
!!!  call nullify_p2pComms(lin%comsr)
!!!  !call nullify_p2pCommsGatherPot(lin%comgp)
!!!  call nullify_p2pComms(lin%comgp)
!!!  call nullify_largeBasis(lin%lb)
!!!  call nullify_local_zone_descriptors(lin%lzd)
!!!  !call nullify_p2pCommsOrthonormality(lin%comon)
!!!  call nullify_p2pComms(lin%comon)
!!!  call nullify_overlapParameters(lin%op)
!!!  call nullify_matrixDescriptors(lin%mad)
!!!
!!!end subroutine nullify_linearParameters


subroutine nullify_overlap_parameters_matrix(opm)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(overlap_parameters_matrix),intent(inout):: opm

  nullify(opm%noverlap)
  nullify(opm%overlaps)
  nullify(opm%olrForExpansion)
  nullify(opm%olr)

end subroutine nullify_overlap_parameters_matrix



subroutine nullify_matrixLocalizationRegion(mlr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(matrixLocalizationRegion),intent(out):: mlr

  nullify(mlr%indexInGlobal)

end subroutine nullify_matrixLocalizationRegion

