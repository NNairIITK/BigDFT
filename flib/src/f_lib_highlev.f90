!> @file
!!  High level routines which needs more medium-level modules of the f_lib.
!!  They should be external, in the sense that no interface should be needed to call them.
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>print error information about last error
subroutine f_dump_last_error()
  use dictionaries, only: f_get_error_dict,f_get_last_error,max_field_length
  use yaml_output, only: yaml_dict_dump,yaml_map
  implicit none
  !local variables
  integer :: ierr
  character(len=max_field_length) :: add_msg

  ierr=f_get_last_error(add_msg)

  if (ierr /=0) then
     call yaml_dict_dump(f_get_error_dict(ierr))
     if (trim(add_msg)/= 'UNKNOWN') call yaml_map('Additional Info',add_msg)
  end if
end subroutine f_dump_last_error

!> dump the list of possible errors as they are defined at present
subroutine f_dump_possible_errors(extra_msg)
  use yaml_output
  use dictionaries, only: f_get_error_definitions
  implicit none
  character(len=*), intent(in) :: extra_msg
  
  call yaml_newline()
  call yaml_comment('Error list',hfill='~')
  call yaml_open_map('List of errors defined so far')
!  call yaml_dict_dump(f_get_error_definitions(),verbatim=.true.)
  call yaml_dict_dump(f_get_error_definitions())
  call yaml_close_map()
  call yaml_comment('End of error list',hfill='~')
  if (len_trim(extra_msg) > 0) then
     call yaml_map('Additional Info',trim(extra_msg))
  else
     call yaml_map('Dump ended',.true.)
  end if
end subroutine f_dump_possible_errors

!!$!> routine which dump an error according to the arguments.
!!$subroutine f_dump_error(newerror_code,err_msg)
!!$  use dictionaries, only: f_get_error_dict
!!$  use yaml_output, only: yaml_dict_dump,yaml_map
!!$  implicit none
!!$  integer, intent(in) :: newerror_code
!!$  character(len=*), intent(in) :: err_msg
!!$
!!$  !dump the error message
!!$  call yaml_dict_dump(f_get_error_dict(newerror_code))
!!$  if (len_trim(err_msg)/=0) call yaml_map('Additional Info',err_msg)
!!$
!!$end subroutine f_dump_error

!>routine which initializes f_lib global pointers, to be called in general
subroutine f_lib_initialize()
  implicit none
  
end subroutine f_lib_initialize

!>routine which finalize f_lib 
subroutine f_lib_finalize()
  use dictionaries
  use dynamic_memory
  implicit none
  call f_malloc_finalize()
  call f_err_finalize()
end subroutine f_lib_finalize

