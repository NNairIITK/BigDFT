!> @file
!!  High level routines which needs more medium-level modules of the f_lib
!!  They should be external, in the sense that no interface should be needed to call them
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!print error information about last error
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

subroutine f_dump_all_errors()
  use dictionaries, only: f_get_error_dict,f_get_past_error,&
       max_field_length,f_get_no_of_errors
  use yaml_output, only: yaml_dict_dump,yaml_map
  implicit none
  !local variables
  integer :: ierr
  character(len=max_field_length) :: add_msg
  
  do ierr=0,f_get_no_of_errors()-1
     call yaml_dict_dump(f_get_error_dict(f_get_past_error(ierr,add_msg)))
     if (trim(add_msg)/= 'UNKNOWN') call yaml_map('Additional Info',add_msg)
  end do
end subroutine f_dump_all_errors

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

subroutine initialize_flib_errors()
  use dictionaries, only: dictionaries_errors
  use yaml_output, only: yaml_output_errors
  use yaml_parse, only: yaml_parse_errors
  use dynamic_memory, only: dynamic_memory_errors
  implicit none

  call dictionaries_errors()
  call yaml_output_errors()
  call yaml_parse_errors()
  call dynamic_memory_errors()
  
end subroutine initialize_flib_errors

!>routine which initializes f_lib global pointers, to be called in general
subroutine f_lib_initialize()
  use dictionaries, only: f_err_initialize
  use dynamic_memory, only: f_malloc_initialize
  implicit none
  
  !general initialization, for lowest level f_lib calling
  call f_err_initialize()
  call initialize_flib_errors()
  !initializtion of memory allocation profiling
  call f_malloc_initialize()

end subroutine f_lib_initialize

!>routine which finalize f_lib 
subroutine f_lib_finalize()
  use dictionaries, only: f_err_finalize
  use dynamic_memory, only: f_malloc_finalize
  use yaml_output, only: yaml_close_all_streams
  implicit none
  call f_malloc_finalize()

  !general finalization, the f_lib should come back to uninitialized status
  call yaml_close_all_streams()
  call f_err_finalize()

end subroutine f_lib_finalize

