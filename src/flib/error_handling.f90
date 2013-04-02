!> @file
!!  Module defining handling of errors
!!  In the spirit of this module each error message is handled separately
!!  The user of this module is able to define the stopping routine, the new error and the error codes
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

module error_handling
  use dictionaries, dict_msg_len=>max_field_length
  implicit none

  private

  !some parameters
  integer, parameter :: ERR_SUCCESS=0
  character(len=*), parameter :: errid='Id'
  character(len=*), parameter :: errmsg='Message'
  character(len=*), parameter :: errunspec='UNSPECIFIED'
  character(len=*), parameter :: errundef='UNKNOWN'
  

  logical :: exceptions=.false. !<determine whether exception handling is active
  logical :: abort_has_been_defined=.false.
  integer :: last_err_code=ERR_SUCCESS

  type(dictionary), pointer :: dict_errors=>null() !< the global dictionaries of possible errors, nullified if not initialized
  type(dictionary), pointer :: dict_present_error=>null() !< local pointer of present error, nullified if success

  external :: customized_abort

contains

  !> initialize module
  subroutine err_initialize()
    implicit none
    
    call dict_init(dict_errors)
  end subroutine err_initialize

  subroutine err_finalize()
    implicit none
    
    call dict_free(dict_errors)
  end subroutine err_finalize

  !>define a error event, pointing the dictionary to it
  !! in case the error is absent, return the error code associated to it
  !! if exceptions is set to false, stops the code with the suitable stopping routine
  subroutine err_event(code,err_id,err_msg,err_code)
    implicit none
    integer, intent(in), optional :: code !< error code event, to be present in the dictionary
    character(len=*), intent(in), optional :: err_id !< ID of the error, in order of definition
    character(len=*), intent(in), optional :: err_msg !< message to be printed when the error is produced or retrieved
    integer, intent(out), optional :: err_code !<error code event to be retrieved if a new error is defined
    !local variables

    nullify(dict_present_error)
    if (present(code)) then
       if (code < dict_len(dict_errors)) then
          !the error has to be retrieved
          dict_present_error=>dict_errors//code
          last_err_code=code
       end if
    end if

    !define the error otherwise
    if (.not. associated(dict_present_error)) then

       call dict_init(dict_present_error)
       if (present(err_msg)) then
          call set(dict_present_error//errmsg,err_msg)
       else
          call set(dict_present_error//errmsg,errunspec)
       end if
       if (present(err_id)) then
          call set(dict_present_error//errid,err_id)
       else
          if (present(err_msg)) then
             call set(dict_present_error//errid,errunspec)
          else
             call set(dict_present_error//errid,errundef)
          end if
       end if
       if (associated(dict_errors)) call add(dict_errors,dict_present_error)
       last_err_code=dict_len(dict_errors)
    end if
    
    if (present(err_code)) then
       err_code=last_err_code
    end if

    call err_exception()

    !dry run in case of not associated dictionary
    if (.not. associated(dict_errors)) then
       call yaml_dict_free(dict_present_error)
    end if
  end subroutine err_event
 
  !> routine which makes the system crash if there is a problem and cpontinues in case of a exception raised
  subroutine err_exception()
    implicit none

    !stopping routine 
    if (.not. exceptions) then
       call yaml_dict_dump(dict_present_error)
       call err_abort()
    end if
    
  end subroutine err_exception

  !subroutine which defines the way the system stops
  subroutine err_abort()
    implicit none
    
    if (abort_has_been_defined) then
       call customized_abort()
    else
       stop
    end if
  
  end subroutine err_abort
  
  !> Defines the error routine which have to be used
  subroutine err_abort_function(customized_abort)
    implicit none
    interface
       subroutine customized_abort()
       end subroutine customized_abort
    end interface
  end subroutine err_abort_function

  function err_last()
    implicit none
    integer err_last
    err_last=last_err_code
  end function err_last

  function err_last_msg()
    implicit none
    character(len=dict_msg_len) :: err_last_msg
    
    if (associated(dict_errors)) then
       err_last_msg=dict_errors//last_err_code//errmsg
    else if (last_err_code /= ERR_SUCCESS) then
       err_last_msg=errundef
    else
       err_last_msg=repeat(' ',len(err_last_msg))
    end if
  end function err_last_msg

end module error_handling
