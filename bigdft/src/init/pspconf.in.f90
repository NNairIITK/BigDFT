!> @file
!! Contains all data related to builtin pseudopotentials
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Routine containing the data for the builtin pseudopotentials
!! The file pspconf.f90 is generated from the file pspconf.in.f90.
subroutine psp_from_data(symbol, nzatom, nelpsp, npspcode, ixc, psppar, exists)
  use module_base
  use module_xc
  implicit none
  
  character(len = *), intent(in) :: symbol
  integer, intent(inout) :: ixc
  integer, intent(out) :: nzatom, nelpsp, npspcode
  real(gp), intent(out) :: psppar(0:4,0:6)
  logical, intent(out) :: exists

  character(len=500) :: name_ixc, name_xcpsp(3)

  exists      = .false.
  nzatom      = 0
  nelpsp      = 0
  npspcode    = 1
  psppar(:,:) = 0._gp

  if (ixc < 0) then
     call xc_get_name(name_ixc, ixc, XC_MIXED)
  else
     call xc_get_name(name_ixc, ixc, XC_ABINIT)
  end if
  call xc_get_name(name_xcpsp(1),  1, XC_ABINIT)
  call xc_get_name(name_xcpsp(2), 11, XC_ABINIT)

  if (trim(symbol) == "X") then
     return
!!PSP_TABLE!!
  end if
end subroutine psp_from_data
