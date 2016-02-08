!> @file
!! Check pseudopotentials
!! @author
!!    Copyright (C) 2011-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to check pseudopotentials
program psp_test

  use module_base
  use yaml_output
  use pseudopotentials
  implicit none

  character(len = 1024) :: path
  integer :: istat, i, j
  logical :: exists
  
  integer :: nzatom, nelpsp, npspcode, ixcpsp
  real(gp), dimension(0:4,0:6) :: psppar
  real(gp), dimension(3) :: radii_cf

  integer :: nzatom_d, nelpsp_d, npspcode_d
  real(gp), dimension(0:4,0:6) :: psppar_d
  !ALEX: If psp_from_file supports NLCC, it needs additional arguments:
  real(gp):: rcore, qcore
  logical:: donlcc, pawpatch
  type(io_stream) :: ios

  call f_lib_initialize()

  call get_command_argument(1, value = path, status = istat)
  !if (istat /= 0) stop "argument"
  if (istat /= 0) call f_err_throw('No argument',err_name='GENERIC_ERROR')

  call f_iostream_from_file(ios,path)
  call psp_from_stream(ios, nzatom, nelpsp, npspcode, ixcpsp, psppar, &
         donlcc, rcore, qcore, radii_cf,pawpatch)
  call f_iostream_release(ios)
  !if (.not. exists) stop "psp file"

  i = index(path, "/", back = .true.)
  j = index(path, "-", back = .true.)

  !write(*, "(1x,A6,1x)", advance = "NO") trim(path(i + 1:))

  call psp_from_data(path(i + 1:j - 1), nzatom_d, nelpsp_d, npspcode_d, &
       & ixcpsp, psppar_d, exists)
  !if (.not. exists) stop
  if (.not. exists) call f_err_throw('The psp data'+path(i+1:)+'does not exist!',err_name='GENERIC_ERROR')
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

  !if (nzatom /= nzatom_d) stop "nzatom"
  !if (nelpsp /= nelpsp_d) stop "nelpsp"
  !if (npspcode /= npspcode_d) stop "npspcode"
  !if (sum(abs(psppar - psppar_d)) > 1.d-8) stop "psppar"
  if (nzatom /= nzatom_d) &
     & call f_err_throw('Inconsistency for '+path(i+1:)+'nzatom /= nzatom_d)',err_name='GENERIC_ERROR')
  if (nelpsp /= nelpsp_d) call f_err_throw('Inconsistency for '+path(i+1:)+'nelpsp /= nelpsp_d)',err_name='GENERIC_ERROR')
  if (npspcode /= npspcode_d) call f_err_throw('Inconsistency for '+path(i+1:)+'npspcode /= npspcode_d)',err_name='GENERIC_ERROR')
  if (sum(abs(psppar - psppar_d)) > 1.d-8) &
     & call f_err_throw('Inconsistency (psppar /= psppar_d)',err_name='GENERIC_ERROR')

  !write(*,"(A,F12.6,A)") "OK (checksum is:", sum(psppar) + nzatom + nelpsp + npspcode, ")"
  call yaml_mapping_open(trim(path(i + 1:)),flow=.True.)
    call yaml_map('check',.True.)
    call yaml_map('checksum',sum(psppar) + nzatom + nelpsp + npspcode)
  call yaml_mapping_close()

  call f_lib_finalize_noreport()

end program psp_test
