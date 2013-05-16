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
  use exception_callbacks
  use dictionaries, dict_msg_len=>max_field_length
  implicit none

  private

  !some parameters
  character(len=*), parameter :: errid='Id'
  character(len=*), parameter :: errmsg='Message'
  character(len=*), parameter :: erract='Action'
  character(len=*), parameter :: errclbk='Callback Procedure Address'
  character(len=*), parameter :: errclbkadd='Callback Procedure Data Address'

  character(len=*), parameter :: errunspec='UNSPECIFIED'
  character(len=*), parameter :: errundef='UNKNOWN'

  integer :: ERR_GENERIC,ERR_SUCCESS,ERR_NOT_DEFINED

  type(dictionary), pointer :: dict_errors=>null() !< the global dictionaries of possible errors, nullified if not initialized
  type(dictionary), pointer :: dict_present_error=>null() !< local pointer of present error, nullified if success

  public :: f_err_initialize,f_err_finalize
  public :: f_err_define,f_err_check,f_err_raise,f_err_clean

  !public variables of the callback module
  public :: f_err_set_callback,f_err_unset_callback
  public :: f_err_severe,f_err_severe_override,f_err_severe_restore
  public :: f_loc

contains

  !> initialize module
  subroutine f_err_initialize()
    implicit none
    !local variables
    call f_err_unset_callback()
    call dict_init(dict_present_error)
    call dict_init(dict_errors)
    !initialize the dictionary with the generic case
    call f_err_define('SUCCESS','Operation has succeeded',ERR_SUCCESS,err_action='No action')
    call f_err_define('GENERIC_ERROR',errunspec,ERR_GENERIC,err_action=errundef)
    call f_err_define('ERR_NOT_DEFINED','The error id or name is invalid',ERR_NOT_DEFINED,&
         err_action='Control if the err id exists')
  end subroutine f_err_initialize

  subroutine f_err_finalize()
    implicit none
    call f_err_unset_callback()
    call dict_free(dict_errors)
    call dict_free(dict_present_error)
  end subroutine f_err_finalize

  !> define a new error specification and returns the corresponding error code
  !! optionally, a error-specific callback function can be defined
  subroutine f_err_define(err_name,err_msg,err_id,err_action,callback,callback_data)
    implicit none
    character(len=*), intent(in) :: err_name,err_msg
    integer, intent(out) :: err_id
    integer(kind=8), intent(in), optional :: callback_data
    character(len=*), intent(in), optional :: err_action
    external :: callback
    optional :: callback
    !local variable
    type(dictionary), pointer :: dict_error

    !callback stuff to be done
    err_id=ERR_GENERIC
    if (associated(dict_errors)) then
       err_id=dict_len(dict_errors)
    end if

    call dict_init(dict_error)
    call set(dict_error//err_name//errid,err_id)
    call set(dict_error//err_name//errmsg,err_msg)
    if (present(err_action)) then
       if (len_trim(err_action) /=0)&
            call set(dict_error//err_name//erract,trim(err_action))
    end if
    if (present(callback)) then
       call set(dict_error//err_name//errclbk,f_loc(callback))
       if (present(callback_data)) &
            call set(dict_error//err_name//errclbkadd,callback_data)
    end if
    !callback stuff to be done
    if (associated(dict_errors)) call add(dict_errors,dict_error)

  end subroutine f_err_define

  !> this function returns true if a generic error has been raised  !! in case of specified errors, it returns true if an error of this kind has been raised
  function f_err_check(err_id,err_name)
    implicit none
    integer, intent(in), optional :: err_id !< the code of the error to be checked for
    character(len=*), intent(in), optional :: err_name !name of the error to search
    logical :: f_err_check
    include 'get_err-inc.f90'

    !check if a particular error has been found
    if (get_error==-1) then
       f_err_check = dict_len(dict_present_error) /= 0 
    else
       !othewise check is some error is present
       f_err_check =  get_error/=0
    end if

  end function f_err_check

  !> this routine should be generalized to allow the possiblity of addin customized message at the 
  !! raise of the error. Also customized callback should be allowed
  function f_err_raise(condition,err_msg,err_id,err_name,callback,callback_data)
    use yaml_output, only: yaml_dict_dump,yaml_map
    implicit none
    logical, intent(in), optional :: condition !< the condition which raise the error
    integer, intent(in), optional :: err_id !< the code of the error to be raised.
                                              !! it should already have been defined by f_err_define
    character(len=*), intent(in), optional :: err_name,err_msg !search for the error and add a message to it 
    integer(kind=8), intent(in), optional :: callback_data
    external :: callback
    optional :: callback
    logical :: f_err_raise
    !local variables
    integer :: new_errcode
    integer(kind=8) :: clbk_add,clbk_data_add
    character(len=dict_msg_len), dimension(1) :: keys
    type(dictionary), pointer :: dict_tmp

    if (present(condition)) then
       f_err_raise=condition
    else
       f_err_raise=.true.
    end if

    !once the error has been identified add it to the present errors and call callback function if needed
    if (f_err_raise) then
       !search the error ID (the generic error is always the first)
       new_errcode=ERR_GENERIC
       if (present(err_name)) then
          new_errcode= max(dict_errors .index. err_name,ERR_GENERIC)
!!$          call yaml_map('Error raised, name entered',err_name)
!!$          call yaml_map('Errorcode found',new_errcode)
!!$          call yaml_map('Index function returned',dict_errors .index. err_name)
!!$          call yaml_dict_dump(dict_errors)
!!$          call yaml_map('Dump ended',.true.)
       else if (present(err_id)) then
          new_errcode=ERR_GENERIC
          if (err_id < dict_len(dict_errors)) new_errcode=err_id
       end if

       call add(dict_present_error,new_errcode)
       !       call err_exception(new_errcode)
       !call yaml_dict_dump(dict_errors) !a routine to plot all created errors should be created
       !print *,'debug',new_errcode
       dict_tmp=>dict_errors//new_errcode
       call yaml_dict_dump(dict_tmp)
       if (present(err_msg)) call yaml_map('Additional Info',err_msg)
       !identify callback function 
       clbk_add=callback_add
       clbk_data_add=callback_data_add
       if (present(callback_data)) clbk_data_add=callback_data
       if (present(callback)) then
          clbk_add=f_loc(callback)
       else
          !find the callback in the error definition
          !that is how dict_keys function it should be called
          if (dict_size(dict_tmp) <= size(keys)) keys=dict_keys(dict_tmp)
          dict_tmp=>dict_tmp//trim(keys(1))
          if (has_key(dict_tmp,errclbk)) clbk_add=dict_tmp//errclbk
          if (has_key(dict_tmp,errclbkadd)) clbk_data_add=dict_tmp//errclbkadd
       end if
       call err_abort(clbk_add,clbk_data_add)
    end if
  end function f_err_raise

  !> clean the dictionary of present errors
  subroutine f_err_clean()
    implicit none
    call dict_free(dict_present_error)
    call dict_init(dict_present_error)
  end subroutine f_err_clean

end module error_handling
