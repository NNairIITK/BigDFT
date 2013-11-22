!> @file
!! Routine to flush a unit file 
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Routine to flush a unit file 
subroutine bigdft_utils_flush(unit)
  use yaml_output
  use dictionaries
  implicit none
  integer, intent(in) :: unit
  !local variables
  integer :: ierr
  logical :: unit_is_connected

  unit_is_connected=.false.
  inquire(unit=unit,opened=unit_is_connected,iostat=ierr)
  if (f_err_raise(ierr /=0,'error in unit inquiring, ierr='//trim(yaml_toa(ierr)),&
        err_name='YAML_INVALID')) return
  if (unit_is_connected) flush(unit=unit)
END SUBROUTINE bigdft_utils_flush
