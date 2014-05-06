!> @file
!! Routines to nullify the structures used in the linear version
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Nullify the linear Input parameters
subroutine nullifyInputLinparameters(lin)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  type(linearInputParameters),intent(inout):: lin

  nullify(lin%locrad)
  nullify(lin%potentialPrefac_lowaccuracy)
  nullify(lin%potentialPrefac_highaccuracy)
  nullify(lin%norbsPerType)
  nullify(lin%potentialPrefac_ao)
  nullify(lin%locrad_type)
  nullify(lin%kernel_cutoff_FOE)
  nullify(lin%kernel_cutoff)

end subroutine nullifyInputLinparameters


subroutine nullify_p2pComms(p2pcomm)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_p2pComms
  implicit none

  ! Calling argument
  type(p2pComms),intent(inout):: p2pcomm

  nullify(p2pcomm%noverlaps)
  nullify(p2pcomm%recvBuf)
  nullify(p2pcomm%comarr)
  nullify(p2pcomm%ise)
  nullify(p2pcomm%mpi_datatypes)

end subroutine nullify_p2pComms



subroutine nullify_orbitals_data(orbs)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(out):: orbs
  
  nullify(orbs%norb_par)
  nullify(orbs%iokpt)
  nullify(orbs%ikptproc)
  nullify(orbs%inwhichlocreg)
  nullify(orbs%onwhichatom)
  nullify(orbs%isorb_par)
  nullify(orbs%eval)
  nullify(orbs%occup)
  nullify(orbs%spinsgn)
  nullify(orbs%kwgts)
  nullify(orbs%kpts)
  nullify(orbs%ispot)
  orbs%npsidim_orbs=1
  orbs%npsidim_comp=1

end subroutine nullify_orbitals_data
