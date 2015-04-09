!> @file
!! Routine to tests f_utils module
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed oneder the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine f_utils_test()
  use f_utils
  use yaml_output
  implicit none
  integer :: unt,unt2,u
  double precision :: t0

  !wait one second
  t0=dble(f_time())*1.d-9
  call yaml_map('Absolute time before pause (ns)',t0)
  call f_pause(1)
  call yaml_map('Time spent after pause (s)',dble(f_time())*1.d-9-t0)
  !open files and get free_units
  unt=f_get_free_unit()
  call yaml_map('First unit which can be opened',unt)
  call f_open_file(unt,'Testopen')
  call yaml_map('Opened file in unit',unt)
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Second unit which can be opened',unt2)
  call f_open_file(unt2,'Testopen2')
  call yaml_map('Opened file in unit',unt2)
  call f_close(unt)
  call f_delete_file('Testopen')
  call f_delete_file('Testopen2')
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Third unit which can be opened',unt2)

  call yaml_mapping_open('Bastian test')
  u=f_get_free_unit()                                                   
  call yaml_map('First file',u)
  open(unit=u,file='test') 
  !call f_open_file(u,'test')
  call yaml_map('First file opened',u)
  unt=f_get_free_unit()                                                 
  call yaml_map('Second file',unt)
  open(unit=unt,file='test2')
  call yaml_map('Second file opened',unt)
  unt2=f_get_free_unit()                                      
  call yaml_map('Third file',unt2)
  open(unit=unt2,file='test3')     
  call yaml_map('Third file opened',unt2)
  call f_delete_file('test')
  call f_delete_file('test2')
  call f_delete_file('test3')
  call yaml_mapping_close()
  call yaml_map('If this value is 7 then all files have been correctly closed',f_get_free_unit())

end subroutine f_utils_test
