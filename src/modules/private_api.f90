!> @file
!!  Define interface for private API
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining the inetrfaces for a private API
module module_private_api

  implicit none

  interface

     subroutine atoms_new(atoms)
       use module_types
       implicit none
       type(atoms_data), pointer :: atoms
     END SUBROUTINE atoms_new

     subroutine atoms_new_from_file(lstat, atoms, rxyz, filename, ln)
       use module_base
       use module_types
       use module_interfaces
       implicit none
       logical, intent(out) :: lstat
       type(atoms_data), pointer :: atoms
       integer, intent(in) :: ln
       character, intent(in) :: filename(ln)
       real(gp), dimension(:,:), pointer :: rxyz
     END SUBROUTINE atoms_new_from_file

     subroutine atoms_free(atoms)
       use module_types
       implicit none
       type(atoms_data), pointer :: atoms
     END SUBROUTINE atoms_free

     subroutine atoms_set_n_atoms(atoms, rxyz, nat)
       use module_types
       use memory_profiling
       implicit none
       type(atoms_data), intent(inout) :: atoms
       real(gp), dimension(:,:), pointer :: rxyz
       integer, intent(in) :: nat
     END SUBROUTINE atoms_set_n_atoms

     subroutine atoms_set_n_types(atoms, ntypes)
       use module_types
       implicit none
       type(atoms_data), intent(inout) :: atoms
       integer, intent(in) :: ntypes
     END SUBROUTINE atoms_set_n_types

     subroutine atoms_set_name(atoms, ityp, name)
       use module_types
       implicit none
       type(atoms_data), intent(inout) :: atoms
       integer, intent(in) :: ityp
       character, intent(in) :: name(20)
     END SUBROUTINE atoms_set_name

     subroutine atoms_sync(atoms, alat1, alat2, alat3, geocode, format, units)
       use module_types
       implicit none
       type(atoms_data), intent(inout) :: atoms
       real(gp), intent(in) :: alat1, alat2, alat3
       character(len=1), intent(in) :: geocode
       character, intent(in) :: format(5)
       character, intent(in) :: units(20)
     END SUBROUTINE atoms_sync

     subroutine atoms_copy_nat(atoms, nat)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, intent(out) :: nat
     END SUBROUTINE atoms_copy_nat

     subroutine atoms_copy_ntypes(atoms, ntypes)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, intent(out) :: ntypes
     END SUBROUTINE atoms_copy_ntypes

     subroutine atoms_get_iatype(atoms, iatype)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: iatype
     END SUBROUTINE atoms_get_iatype

     subroutine atoms_get_iasctype(atoms, iasctype)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: iasctype
     END SUBROUTINE atoms_get_iasctype

     subroutine atoms_get_natpol(atoms, natpol)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: natpol
     END SUBROUTINE atoms_get_natpol

     subroutine atoms_get_ifrztyp(atoms, ifrztyp)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: ifrztyp
     END SUBROUTINE atoms_get_ifrztyp

     subroutine atoms_get_nelpsp(atoms, nelpsp)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: nelpsp
     END SUBROUTINE atoms_get_nelpsp

     subroutine atoms_get_npspcode(atoms, npspcode)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: npspcode
     END SUBROUTINE atoms_get_npspcode

     subroutine atoms_get_nzatom(atoms, nzatom)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: nzatom
     END SUBROUTINE atoms_get_nzatom

     subroutine atoms_get_nlcc_ngv(atoms, nlcc_ngv)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: nlcc_ngv
     END SUBROUTINE atoms_get_nlcc_ngv

     subroutine atoms_get_nlcc_ngc(atoms, nlcc_ngc)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: nlcc_ngc
     END SUBROUTINE atoms_get_nlcc_ngc

     subroutine atoms_get_ixcpsp(atoms, ixcpsp)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, dimension(:), pointer :: ixcpsp
     END SUBROUTINE atoms_get_ixcpsp

     subroutine atoms_get_amu(atoms, amu)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       real(gp), dimension(:), pointer :: amu
     END SUBROUTINE atoms_get_amu

     subroutine atoms_get_aocc(atoms, aocc)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       real(gp), dimension(:,:), pointer :: aocc
     END SUBROUTINE atoms_get_aocc

     subroutine atoms_get_radii_cf(atoms, radii_cf)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       real(gp), dimension(:,:), pointer :: radii_cf
     END SUBROUTINE atoms_get_radii_cf

     subroutine atoms_get_psppar(atoms, psppar)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       real(gp), dimension(:,:,:), pointer :: psppar
     END SUBROUTINE atoms_get_psppar

     subroutine atoms_get_nlccpar(atoms, nlccpar)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       real(gp), dimension(:,:), pointer :: nlccpar
     END SUBROUTINE atoms_get_nlccpar

     subroutine atoms_get_ig_nlccpar(atoms, ig_nlccpar)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       real(gp), dimension(:,:), pointer :: ig_nlccpar
     END SUBROUTINE atoms_get_ig_nlccpar

     subroutine atoms_copy_geometry_data(atoms, geocode, format, units)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       character(len=1), intent(out) :: geocode !< @copydoc poisson_solver::doc::geocode
       character, intent(out) :: format(5)
       character, intent(out) :: units(20)
     END SUBROUTINE atoms_copy_geometry_data

     subroutine atoms_copy_psp_data(atoms, natsc, donlcc)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, intent(out) :: natsc
       logical, intent(out) :: donlcc
     END SUBROUTINE atoms_copy_psp_data

     subroutine atoms_copy_name(atoms, ityp, name, ln)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer, intent(in) :: ityp
       character(len=*), intent(out) :: name
       integer, intent(out) :: ln
     END SUBROUTINE atoms_copy_name

     subroutine atoms_copy_alat(atoms, alat1, alat2, alat3)
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       real(gp), intent(out) :: alat1, alat2, alat3
     END SUBROUTINE atoms_copy_alat

     subroutine atoms_write(atoms, filename, filelen, rxyz, forces, energy, comment, ln)
       use module_types
       implicit none
       integer, intent(in) :: ln, filelen
       character, intent(in) :: comment(ln)
       character, intent(in) :: filename(filelen)
       type(atoms_data), intent(in) :: atoms
       real(gp), intent(in) :: energy
       real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
       real(gp), dimension(:,:), pointer :: forces
     END SUBROUTINE atoms_write

     subroutine localfields_copy_metadata(denspot, rhov_is, hgrid, psoffset)
       use module_types
       implicit none
       type(DFT_local_fields), intent(in) :: denspot
       integer, intent(out) :: rhov_is
       real(gp), intent(out) :: hgrid(3)
       real(dp), intent(out) :: psoffset
     END SUBROUTINE localfields_copy_metadata
  END INTERFACE
END MODULE module_private_api
