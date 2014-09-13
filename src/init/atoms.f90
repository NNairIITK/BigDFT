!> @file
!! Routines for handling the structure atoms_data 
!! @author
!!    Copyright (C) 2011-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine atoms_write(atoms, filename, forces, energy, comment)
  use module_defs, only: gp
  use module_types
  use module_atoms, only: astruct_dump_to_file
  implicit none
  character(len = *), intent(in) :: comment
  character(len = *), intent(in) :: filename
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(:,:), pointer :: forces

  if (associated(forces)) then
     call astruct_dump_to_file(atoms%astruct,filename,comment,&
          energy,forces=forces)
  else
     call astruct_dump_to_file(atoms%astruct,filename,comment,&
          energy)
     !call write_atomic_file(filename,energy,atoms%astruct%rxyz,atoms%astruct,comment)
  end if
END SUBROUTINE atoms_write


!> Deallocate a new atoms_data type, for bindings.
subroutine atoms_empty(atoms)
  use module_atoms, only: atoms_data, deallocate_atoms_data
  implicit none
  type(atoms_data), intent(inout) :: atoms

  call deallocate_atoms_data(atoms)
END SUBROUTINE atoms_empty


subroutine atoms_set_name(atoms, ityp, name)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: ityp
  character(len=1), dimension(20), intent(in) :: name

  write(atoms%astruct%atomnames(ityp), "(20A1)") name
END SUBROUTINE atoms_set_name


subroutine astruct_set_geometry(astruct, alat, geocode, format, units)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  real(gp), intent(in) :: alat(3)
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character, intent(in) :: format(5)
  character, intent(in) :: units(20)

  astruct%cell_dim(:) = alat(:)
  astruct%geocode = geocode
  write(astruct%inputfile_format, "(5A1)") format
  write(astruct%units, "(20A1)") units
END SUBROUTINE astruct_set_geometry

!> Accessors for bindings.
subroutine atoms_get(atoms, astruct, symObj)
  use module_types
  implicit none
  type(atoms_data), intent(in), target :: atoms
  type(atomic_structure), pointer :: astruct
  type(symmetry_data), pointer :: symObj

  astruct => atoms%astruct
  symObj => atoms%astruct%sym
END SUBROUTINE atoms_get


subroutine astruct_copy_nat(astruct, nat)
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  integer, intent(out) :: nat

  nat = astruct%nat
END SUBROUTINE astruct_copy_nat


subroutine astruct_copy_ntypes(astruct, ntypes)
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  integer, intent(out) :: ntypes

  ntypes = astruct%ntypes
END SUBROUTINE astruct_copy_ntypes


subroutine atoms_get_iatype(atoms, iatype)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: iatype

  iatype => atoms%astruct%iatype
END SUBROUTINE atoms_get_iatype


!!$subroutine atoms_get_iasctype(atoms, iasctype)
!!$  use module_types
!!$  implicit none
!!$  type(atoms_data), intent(in) :: atoms
!!$  integer, dimension(:), pointer :: iasctype
!!$  !local variables
!!$  integer
!!$
!!$  iasctype => atoms%iasctype
!!$END SUBROUTINE atoms_get_iasctype


subroutine atoms_get_natpol(atoms, natpol)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: natpol

  natpol => atoms%astruct%input_polarization
END SUBROUTINE atoms_get_natpol


subroutine atoms_get_ifrztyp(atoms, ifrztyp)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: ifrztyp

  ifrztyp => atoms%astruct%ifrztyp
END SUBROUTINE atoms_get_ifrztyp


subroutine atoms_get_rxyz(atoms, rxyz)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: rxyz

  rxyz => atoms%astruct%rxyz
END SUBROUTINE atoms_get_rxyz


subroutine atoms_get_nelpsp(atoms, nelpsp)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nelpsp

  nelpsp => atoms%nelpsp
END SUBROUTINE atoms_get_nelpsp


subroutine atoms_get_npspcode(atoms, npspcode)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: npspcode

  npspcode => atoms%npspcode
END SUBROUTINE atoms_get_npspcode


subroutine atoms_get_nzatom(atoms, nzatom)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nzatom

  nzatom => atoms%nzatom
END SUBROUTINE atoms_get_nzatom


subroutine atoms_get_nlcc_ngv(atoms, nlcc_ngv)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nlcc_ngv

  nlcc_ngv => atoms%nlcc_ngv
END SUBROUTINE atoms_get_nlcc_ngv
subroutine atoms_get_nlcc_ngc(atoms, nlcc_ngc)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nlcc_ngc

  nlcc_ngc => atoms%nlcc_ngc
END SUBROUTINE atoms_get_nlcc_ngc


subroutine atoms_get_ixcpsp(atoms, ixcpsp)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: ixcpsp

  ixcpsp => atoms%ixcpsp
END SUBROUTINE atoms_get_ixcpsp


subroutine atoms_get_amu(atoms, amu)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:), pointer :: amu

  amu => atoms%amu
END SUBROUTINE atoms_get_amu


!!$subroutine atoms_get_aocc(atoms, aocc)
!!$  use module_types
!!$  implicit none
!!$  type(atoms_data), intent(in) :: atoms
!!$  real(gp), dimension(:,:), pointer :: aocc
!!$
!!$  aocc => atoms%aocc
!!$END SUBROUTINE atoms_get_aocc


!> get radii_cf values
subroutine atoms_get_radii_cf(atoms, radii_cf)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: radii_cf

  radii_cf => atoms%radii_cf
END SUBROUTINE atoms_get_radii_cf


subroutine atoms_get_psppar(atoms, psppar)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:,:), pointer :: psppar

  psppar => atoms%psppar
END SUBROUTINE atoms_get_psppar

subroutine atoms_get_nlccpar(atoms, nlccpar)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: nlccpar

  nlccpar => atoms%nlccpar
END SUBROUTINE atoms_get_nlccpar

!subroutine atoms_get_ig_nlccpar(atoms, ig_nlccpar)
!  use module_types
!  implicit none
!  type(atoms_data), intent(in) :: atoms
!  real(gp), dimension(:,:), pointer :: ig_nlccpar
!
!  ig_nlccpar => atoms%ig_nlccpar
!END SUBROUTINE atoms_get_ig_nlccpar

!> For bindings only, use input_dicts module for Fortran usage.
subroutine astruct_merge_to_dict_binding(dict, astruct)
  use module_atoms, only: wrapper => astruct_merge_to_dict
  use module_types, only: atomic_structure
  use dictionaries, only: dictionary
  implicit none
  type(dictionary), pointer :: dict
  type(atomic_structure), intent(in) :: astruct

  call wrapper(dict, astruct,astruct%rxyz)
END SUBROUTINE astruct_merge_to_dict_binding


subroutine astruct_copy_geometry_data(astruct, geocode, format, units)
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  character(len = 1), intent(out) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len = 5), intent(out) :: format
  character(len = 20), intent(out) :: units

  write(geocode, "(A1)") astruct%geocode
  write(format,  "(A5)") astruct%inputfile_format
  write(units,  "(A20)") astruct%units
END SUBROUTINE astruct_copy_geometry_data


subroutine atoms_copy_psp_data(atoms, natsc, donlcc)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: natsc
  logical, intent(out) :: donlcc

  natsc = atoms%natsc
  donlcc = atoms%donlcc
END SUBROUTINE atoms_copy_psp_data


subroutine astruct_copy_name(astruct, ityp, name, ln)
  use module_types
  implicit none
  !Arguments
  type(atomic_structure), intent(in) :: astruct
  integer, intent(in) :: ityp
  character(len=1), dimension(20), intent(out) :: name
!  character(len=*), intent(out) :: name
  integer, intent(out) :: ln
  !Local variables 
  integer :: i,lname

  if (astruct%ntypes > 0) then
     lname = len(name)
     ln=min(len(trim(astruct%atomnames(ityp))),20)
     !print *,'lnt2',lnt
     do i = 1, ln, 1
     !name(i:i) = astruct%atomnames(ityp)(i:i)
     write(name(i),'(a1)') astruct%atomnames(ityp)(i:i)
     end do
     do i = ln + 1, lname, 1
        name(i) = ' '
     end do
  end if
END SUBROUTINE astruct_copy_name


subroutine astruct_copy_alat(astruct, alat)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(out) :: alat(3)

  alat(1) = astruct%cell_dim(1)
  alat(2) = astruct%cell_dim(2)
  alat(3) = astruct%cell_dim(3)
END SUBROUTINE astruct_copy_alat

!>Write the extra info necessary for the output file
subroutine write_extra_info(extra,natpol,ifrztyp)
  use module_atoms, only: frozen_itof
  use ao_inguess, only: charge_and_spol
  implicit none 
  integer, intent(in) :: natpol,ifrztyp
  character(len=50), intent(out) :: extra
  !local variables
  character(len=4) :: frzchain
  integer :: ispol,ichg

  call charge_and_spol(natpol,ichg,ispol)

  call frozen_itof(ifrztyp,frzchain)

  !takes into account the blocked atoms and the input polarisation
  if (ispol == 0 .and. ichg == 0 ) then
     write(extra,'(2x,a4)')frzchain
  else if (ispol /= 0 .and. ichg == 0) then
     write(extra,'(i7,2x,a4)')ispol,frzchain
  else if (ichg /= 0) then
     write(extra,'(2(i7),2x,a4)')ispol,ichg,frzchain
  else
     write(extra,'(2x,a4)') ''
  end if

END SUBROUTINE write_extra_info


!!$!> Module used for the input positions lines variables
!!$module position_files
!!$   implicit none
!!$   contains
!!$   subroutine directGetLine(line, ifile, eof)
!!$      !Arguments
!!$      integer, intent(in) :: ifile
!!$      character(len=150), intent(out) :: line
!!$      logical, intent(out) :: eof
!!$      !Local variables
!!$      integer :: i_stat
!!$
!!$      eof = .false.
!!$      read(ifile,'(a150)', iostat = i_stat) line
!!$      if (i_stat /= 0) eof = .true.
!!$   END SUBROUTINE directGetLine
!!$
!!$   subroutine archiveGetLine(line, ifile, eof)
!!$      !Arguments
!!$      integer, intent(in) :: ifile
!!$      character(len=150), intent(out) :: line
!!$      logical, intent(out) :: eof
!!$      !Local variables
!!$      integer :: i_stat
!!$      !The argument ifile is not used but it is used as argument routine
!!$      !eof = .false.
!!$      eof = (ifile /= ifile)
!!$      call extractNextLine(line, i_stat)
!!$      if (i_stat /= 0) eof = .true.
!!$   END SUBROUTINE archiveGetLine
!!$end module position_files
