!> @file
!!  High level routines which needs more medium-level modules of the f_lib
!!  They should be external, in the sense that no interface should be needed to call them
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> routine which dump an error according to the arguments.
subroutine f_dump_error(newerror_code,err_msg)
  use dictionaries, only: f_get_error_dict
  use yaml_output, only: yaml_dict_dump,yaml_map
  implicit none
  integer, intent(in) :: newerror_code
  character(len=*), intent(in) :: err_msg

  !dump the error message
  call yaml_dict_dump(f_get_error_dict(newerror_code))
  if (len_trim(err_msg)/=0) call yaml_map('Additional Info',err_msg)

end subroutine f_dump_error

!>routine which initializes f_lib global pointers, to be called in general
subroutine f_lib_initialize()
  implicit none
  
end subroutine f_lib_initialize

!>routine which finalize f_lib 
subroutine f_lib_finalize()
  use dictionaries
  implicit none

  call f_err_finalize()
end subroutine f_lib_finalize

