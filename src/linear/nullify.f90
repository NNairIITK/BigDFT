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



subroutine nullify_foe(foe_obj)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_foe
  implicit none

  ! Calling argument
  type(foe_data),intent(out):: foe_obj

  nullify(foe_obj%kernel_nseg)
  nullify(foe_obj%kernel_segkeyg)

end subroutine nullify_foe


pure subroutine nullify_sparsematrix(sparsemat)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => nullify_sparseMatrix
  implicit none

  ! Calling argument
  type(sparseMatrix),intent(out):: sparsemat

  nullify(sparsemat%keyv)
  nullify(sparsemat%nsegline)
  nullify(sparsemat%keyg)
  nullify(sparsemat%istsegline)
  nullify(sparsemat%matrix)
  nullify(sparsemat%matrix_compr)
  nullify(sparsemat%matrixp)
  nullify(sparsemat%matrix_comprp)
  nullify(sparsemat%matrixindex_in_compressed_arr)
  nullify(sparsemat%orb_from_index)
  nullify(sparsemat%matrixindex_in_compressed_fortransposed)
  nullify(sparsemat%nvctr_par)
  nullify(sparsemat%isvctr_par)
  nullify(sparsemat%nfvctr_par)
  nullify(sparsemat%isfvctr_par)

end subroutine nullify_sparsematrix

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


subroutine nullify_communications_arrays(comms)
  use module_base
  use communications_base, only: communications_arrays
  implicit none

  ! Calling arguments
  type(communications_arrays),intent(out):: comms

  nullify(comms%ncntd)
  nullify(comms%ncntt)
  nullify(comms%ndspld)
  nullify(comms%ndsplt)
  nullify(comms%nvctr_par)
  
end subroutine nullify_communications_arrays




!!! replace calls to these routines with functions defined above
!subroutine nullify_local_zone_descriptors(lzd)
!  use module_base
!  use module_types
!  use module_interfaces, exceptThisOne => nullify_local_zone_descriptors
!  implicit none
!
!  ! Calling arguments
!  type(local_zone_descriptors),intent(out):: lzd
! 
!  call nullify_locreg_descriptors(lzd%glr)
!  nullify(lzd%llr)
! 
!end subroutine nullify_local_zone_descriptors


!!! replace calls to these routines with functions defined above



!subroutine nullify_collective_comms(collcom)
!  use module_base
!  use module_types
!  implicit none
!  
!  ! Calling arguments
!  type(collective_comms),intent(inout):: collcom
!
!  ! Local variables
!  nullify(collcom%nsendcounts_c)
!  nullify(collcom%nsenddspls_c)
!  nullify(collcom%nrecvcounts_c)
!  nullify(collcom%nrecvdspls_c)
!  nullify(collcom%isendbuf_c)
!  nullify(collcom%iextract_c)
!  nullify(collcom%iexpand_c)
!  nullify(collcom%irecvbuf_c)
!  nullify(collcom%norb_per_gridpoint_c)
!  nullify(collcom%indexrecvorbital_c)
!  nullify(collcom%isptsp_c)
!  nullify(collcom%psit_c)
!  nullify(collcom%nsendcounts_f)
!  nullify(collcom%nsenddspls_f)
!  nullify(collcom%nrecvcounts_f)
!  nullify(collcom%nrecvdspls_f)
!  nullify(collcom%isendbuf_f)
!  nullify(collcom%iextract_f)
!  nullify(collcom%iexpand_f)
!  nullify(collcom%irecvbuf_f)
!  nullify(collcom%norb_per_gridpoint_f)
!  nullify(collcom%indexrecvorbital_f)
!  nullify(collcom%isptsp_f)
!  nullify(collcom%psit_f)
!  nullify(collcom%nsendcounts_repartitionrho)
!  nullify(collcom%nrecvcounts_repartitionrho)
!  nullify(collcom%nsenddspls_repartitionrho)
!  nullify(collcom%nrecvdspls_repartitionrho)
!  nullify(collcom%commarr_repartitionrho)
!
!end subroutine nullify_collective_comms
