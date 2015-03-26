!> @file
!! Routine to tests f_utils module
!! @example utls.f90
!! Examples using the @ref f_utils objects (units and timers)
!! @author
!!    Copyright (C) 2013-2015 BigDFT group
!!    This file is distributed oneder the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine f_utils_test()
  use f_utils
  use yaml_output
  implicit none
  integer :: unt,unt2
  double precision :: t0

  !wait one second
  t0=dble(f_time())*1.d-9
  call yaml_map('Absolute time before pause (ns)',t0)
  call f_pause(1)
  call yaml_map('Time spent after pause (s)',dble(f_time())*1.d-9-t0)
  !open files and get free_units
  unt=f_get_free_unit()
  call f_open_file(unt,'Testopen')
  call yaml_map('Opened file in unit',unt)
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Second unit which can be opened',unt2)
  call f_close(unt)
  call f_delete_file('Testopen')
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Third unit which can be opened',unt2)


end subroutine f_utils_test
