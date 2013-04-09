!> @file
!! Check psuedopotentials
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to check pseudopotentials
program krachboum
  use module_base

  implicit none

  character(len = 1024) :: path
  integer :: istat, i, j
  logical :: exists, read_radii
  
  integer :: nzatom, nelpsp, npspcode, ixcpsp
  real(gp) :: psppar(0:4,0:6), radii_cf(3)

  integer :: nzatom_d, nelpsp_d, npspcode_d
  real(gp) :: psppar_d(0:4,0:6)

  call get_command_argument(1, value = path, status = istat)
  if (istat /= 0) stop "argument"

  call psp_from_file(path, nzatom, nelpsp, npspcode, ixcpsp, psppar, radii_cf, &
       & read_radii, exists)
  if (.not. exists) stop "psp file"

  i = index(path, "/", back = .true.)
  j = index(path, "-", back = .true.)

  write(*, "(1x,A6,1x)", advance = "NO") trim(path(i + 1:))

  call psp_from_data(path(i + 1:j - 1), nzatom_d, nelpsp_d, npspcode_d, &
       & ixcpsp, psppar_d, exists)
  if (.not. exists) stop
  if (exists .and. nelpsp /= nelpsp_d) then
     ! Try with semi-core
     call psp_from_data(path(i + 1:j - 1) // "_sc", nzatom_d, nelpsp_d, npspcode_d, &
          & ixcpsp, psppar_d, exists)
     if (exists .and. nelpsp /= nelpsp_d) then
        ! Try with semi-core +
        call psp_from_data(path(i + 1:j - 1) // "_sc+", nzatom_d, nelpsp_d, npspcode_d, &
             & ixcpsp, psppar_d, exists)
     end if
  end if

  if (nzatom /= nzatom_d) stop "nzatom"
  if (nelpsp /= nelpsp_d) stop "nelpsp"
  if (npspcode /= npspcode_d) stop "npspcode"
  if (sum(abs(psppar - psppar_d)) > 1.d-8) stop "psppar"

  write(*,"(A,F12.6,A)") "OK (checksum is:", sum(psppar) + nzatom + nelpsp + npspcode, ")"

end program krachboum
