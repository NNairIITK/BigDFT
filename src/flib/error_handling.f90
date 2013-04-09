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
  character(len=*), parameter :: errid='Id'
  character(len=*), parameter :: errmsg='Message'
  character(len=*), parameter :: erract='Action'
  character(len=*), parameter :: errclbk='Callback Procedure Address'
  character(len=*), parameter :: errclbkadd='Callback Procedure Data Address'

  character(len=*), parameter :: errunspec='UNSPECIFIED'
  character(len=*), parameter :: errundef='UNKNOWN'

  integer :: ERR_GENERIC,ERR_SUCCESS

  type(dictionary), pointer :: dict_errors=>null() !< the global dictionaries of possible errors, nullified if not initialized
  type(dictionary), pointer :: dict_present_error=>null() !< local pointer of present error, nullified if success

  integer(kind=8) :: callback_add=0
  integer(kind=8) :: callback_data_add=0

  interface f_err_set_callback
     module procedure err_set_callback_simple,err_set_callback_advanced
  end interface

  public :: f_err_initialize,f_err_finalize,f_err_set_callback,f_err_unset_callback,f_err_define,f_err_check,f_err_raise

contains

  function f_loc(routine)
    implicit none
    external :: routine
    integer(kind=8) :: f_loc

    call getlongaddress(routine,f_loc)

  end function f_loc

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
    !local variables
    integer :: ierr,nerr,jerr
    character(len=dict_msg_len) :: name

    nerr=dict_len(dict_present_error)
    f_err_check=nerr/=0
    if (present(err_name)) then
       f_err_check=.false.
       do ierr=0,nerr-1
          jerr=dict_present_error//ierr
          name=dict_key(dict_errors//jerr)
          if (trim(name)==trim(err_name)) then
             f_err_check=.true. 
             exit
          end if
       end do
    else if (present(err_id)) then
       f_err_check=.false.
       do ierr=0,nerr-1
          jerr=dict_present_error//ierr
          if (jerr==err_id) then
             f_err_check=.true. 
             exit
          end if
       end do
    end if
       
  end function f_err_check

  function f_err_raise(condition,err_id,err_name,err_msg,err_action,callback,callback_data)
    implicit none
    logical, intent(in) :: condition !< the condition which raise the error
    integer, intent(in), optional :: err_id !< the code of the error to be raised.
                                              !! it should already have been defined by f_err_define
    character(len=*), intent(in), optional :: err_name,err_msg,err_action !search for the error or define it in case it is absent 
    integer(kind=8), intent(in), optional :: callback_data
    external :: callback
    optional :: callback
    logical :: f_err_raise
    !local variables
    integer :: new_errcode
    character(len=dict_msg_len) :: error_name,error_message,error_action

    f_err_raise=condition

    error_name=repeat(' ',len(error_name))
    error_name=errunspec
    error_message=repeat(' ',len(error_message))
    error_action=repeat(' ',len(error_action))

    !if a name is specified search for it in the error dictionary
    if (present(err_name)) then
       new_errcode= dict_errors .index. err_name
       !print *,'index',new_errcode,'name',err_name
       if (new_errcode == 0) then !the error has not been found
          error_name=err_name
          if (present(err_msg)) error_message=err_msg
          if (present(err_action)) error_action=err_action
          if (present(callback)) then
             if (present(callback_data)) then
                call f_err_define(trim(error_name),trim(error_message),new_errcode,&
                     err_action=error_action,callback=callback,callback_data=callback_data)
             else
                call f_err_define(trim(error_name),trim(error_message),new_errcode,&
                     err_action=error_action,callback=callback)
             end if
          else
             call f_err_define(trim(error_name),trim(error_message),new_errcode,&
                  err_action=error_action)
          end if
       end if
    else if (present(err_id)) then
       !print *,'here',err_id,dict_len(dict_errors)
       if (err_id < dict_len(dict_errors)) then
          new_errcode=err_id
       else
          new_errcode=ERR_GENERIC
       end if
    else
       new_errcode=ERR_GENERIC
    end if

    !define error if there is information to do that
    if (f_err_raise) then
       call add(dict_present_error,new_errcode)
       call err_exception(new_errcode)
    end if
  end function f_err_raise

  !> clean the dictionary of present errors
  subroutine f_err_clean()
    implicit none
    call dict_free(dict_present_error)
    call dict_init(dict_present_error)
  end subroutine f_err_clean

  !> routine which makes the system crash if there is a problem and cpontinues in case of a exception raised
  subroutine err_exception(err_id)
    use yaml_output, only: yaml_dict_dump
    implicit none
    integer, intent(in) :: err_id
    !local variables
    integer(kind=8) :: clbk_add,clbk_data_add
    type(dictionary), pointer :: dict_tmp

    !call yaml_dict_dump(dict_present_error)
    !print *,'err_id',err_id
    !call yaml_dict_dump(dict_errors)
    !print *,'test'
    dict_tmp=>dict_errors//err_id
    call yaml_dict_dump(dict_tmp)
    if (has_key(dict_tmp,errclbk)) then
       clbk_add=dict_tmp//errclbk
       if (has_key(dict_tmp,errclbkadd)) then
          clbk_data_add=dict_tmp//errclbkadd
          call err_abort(clbk_add,clbk_data_add)
       else
          call err_abort(clbk_add,int(0,kind=8))
       end if
    else
       !global case
       call err_abort(callback_add,callback_data_add)
    end if

  end subroutine err_exception

  !subroutine which defines the way the system stops
  subroutine err_abort(callback,callback_data)
    !use metadata_interfaces
    implicit none
    integer(kind=8), intent(in) :: callback,callback_data

    if (callback_data /=0 .and. callback /=0) then
       call call_external_c_fromadd(callback) !for the moment data are ignored
    else if (callback /=0) then
       call call_external_c_fromadd(callback)
    else
       stop
    end if

 end subroutine err_abort

  !> Defines the error routine which have to be used
  subroutine err_set_callback_simple(callback)
    implicit none
    external :: callback

    callback_add=f_loc(callback)
    callback_data_add=0

  end subroutine err_set_callback_simple

  subroutine err_set_callback_advanced(callback,callback_data_address)
    implicit none
    integer(kind=8), intent(in) :: callback_data_address
    external :: callback

    callback_add=f_loc(callback)
    callback_data_add=callback_data_address

  end subroutine err_set_callback_advanced

  subroutine f_err_unset_callback()
    implicit none

    callback_add=0
    callback_data_add=0
  end subroutine f_err_unset_callback

end module error_handling
