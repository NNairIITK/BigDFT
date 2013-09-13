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

!!$module error_handling
!!$  use exception_callbacks
!!$  use dictionaries, dict_msg_len=>max_field_length
!!$  implicit none
!!$
!!$  private
!!$
!!$  !some parameters
!!$  character(len=*), parameter :: errid='Id'
!!$  character(len=*), parameter :: errmsg='Message'
!!$  character(len=*), parameter :: erract='Action'
!!$  character(len=*), parameter :: errclbk='Callback Procedure Address'
!!$  character(len=*), parameter :: errclbkadd='Callback Procedure Data Address'
!!$
!!$  character(len=*), parameter :: errunspec='UNSPECIFIED'
!!$  character(len=*), parameter :: errundef='UNKNOWN'
!!$
!!$  integer :: ERR_GENERIC,ERR_SUCCESS,ERR_NOT_DEFINED
!!$
!!$  type(dictionary), pointer :: dict_errors=>null() !< the global dictionaries of possible errors, nullified if not initialized
!!$  type(dictionary), pointer :: dict_present_error=>null() !< local pointer of present error, nullified if success
!!$
!!$  public :: f_err_initialize,f_err_finalize
!!$  public :: f_err_define,f_err_check,f_err_raise,f_err_clean,f_get_error_dict
!!$
!!$  !public variables of the callback module
!!$  public :: f_err_set_callback,f_err_unset_callback
!!$  public :: f_err_severe,f_err_severe_override,f_err_severe_restore
!!$  public :: f_loc
!!$
!!$contains

  !> initialize module, in case it has not be done so
  !! the double initialization should not be possible when library is
  !! correctly set up
  subroutine f_err_initialize()
    implicit none
    !local variables
    call f_err_unset_callback()
    if (associated(dict_present_error)) then
       call f_err_clean()
    else
       call dict_init(dict_present_error)
    end if
    if (.not. associated(dict_errors)) then
       call dict_init(dict_errors)
       !call dictionaries_errors()
    end if
  end subroutine f_err_initialize

  subroutine f_err_finalize()
    implicit none
    call f_err_unset_callback()
    call f_err_severe_restore()
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

    !assure initialization of the library in case of misuse
    !when f_lib library becomes beta, this should become a severe error
    if (.not. associated(dict_errors)) then !call f_err_initialize()
       write(0,*)'error_handling library not initialized, f_lib_initialized should be called'
       call f_err_severe()
    end if

    err_id=ERR_GENERIC
!    if (associated(dict_errors)) then
       err_id=dict_len(dict_errors)
!    end if

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

    !if (associated(dict_errors)) 
    call add(dict_errors,dict_error)

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
    use yaml_strings, only: yaml_toa
    !use yaml_output, only: yaml_dict_dump,yaml_map
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
!    integer :: new_errcode
    integer(kind=8) :: clbk_data_add
!    integer(kind=8) :: clbk_add
!    character(len=max_field_length), dimension(1) :: keys
!    type(dictionary), pointer :: dict_tmp
    character(len=max_field_length) :: message

    if (present(condition)) then
       f_err_raise=condition
    else
       f_err_raise=.true.
    end if

    !once the error has been identified add it to the present errors and call callback function if needed
    if (f_err_raise) then

       !trow the error with the annoying stuff of optional variables
       if (present(err_msg)) then
          message(1:len(message))=err_msg
       else
          message(1:len(message))='UNKNOWN'
       end if

       clbk_data_add=callback_data_add
       if (present(callback_data)) clbk_data_add=callback_data

       if (present(callback)) then
          if (present(err_id)) then
             call f_err_throw(message,err_id,callback=callback,callback_data=clbk_data_add)
          else if (present(err_name)) then
             call f_err_throw(message,err_name=err_name,callback=callback,callback_data=clbk_data_add) 
          else
             call f_err_throw(message,callback=callback,callback_data=clbk_data_add)
          end if
       else
          if (present(err_id)) then
             call f_err_throw(message,err_id,callback_data=clbk_data_add)
          else if (present(err_name)) then
             call f_err_throw(message,err_name=err_name,callback_data=clbk_data_add)
          else
             call f_err_throw(message,callback_data=clbk_data_add)
          end if
       end if
    end if
  end function f_err_raise

  !raise the error indicated
  subroutine f_err_throw(err_msg,err_id,err_name,callback,callback_data)
    use yaml_strings, only: yaml_toa
    implicit none
    integer, intent(in), optional :: err_id !< the code of the error to be raised.
                                      !! it should already have been defined by f_err_define
    character(len=*), intent(in), optional :: err_name,err_msg !<search for the error and add a message to it 
    integer(kind=8), intent(in), optional :: callback_data
    external :: callback
    optional :: callback
    !local variables
    integer :: new_errcode
    integer(kind=8) :: clbk_add,clbk_data_add
    character(len=max_field_length), dimension(1) :: keys
    type(dictionary), pointer :: dict_tmp
    character(len=max_field_length) :: message
    
    !to prevent infinite loop due to not association of the error handling
    if (.not. associated(dict_present_error)) then
       write(0,*)'error_handling library not initialized'
       call f_err_severe()
    end if

    !search the error ID (the generic error is always the first)
    new_errcode=ERR_GENERIC
    if (present(err_name)) then
       new_errcode= max(dict_errors .index. err_name,ERR_GENERIC)
       !add a potentially verbose error list in case the specified error name has not been found
       if ((dict_errors .index. err_name) < ERR_GENERIC) then
          call f_dump_possible_errors('Error raised, name entered= '//trim(err_name)//&
               '. Errorcode found='//trim(yaml_toa(new_errcode))//&
               '; index function returned'//trim(yaml_toa(dict_errors .index. err_name)))
       end if
       else if (present(err_id)) then
       new_errcode=ERR_GENERIC
       if (err_id < dict_len(dict_errors)) then
          new_errcode=err_id
       else
          !add a potentially verbose error list in case the specified error name has not been found
          call f_dump_possible_errors('Error raised, id='//trim(yaml_toa(err_id)))
       end if
    end if

    if (present(err_msg)) then
       message(1:len(message))=err_msg
    else
       message(1:len(message))='UNKNOWN'
    end if
    call add(dict_present_error,&
         dict_new((/ errid .is. yaml_toa(new_errcode),'Additional Info' .is. message/)))

    !identify callback function 
    clbk_add=callback_add
    clbk_data_add=callback_data_add
    if (present(callback_data)) clbk_data_add=callback_data
    if (present(callback)) then
       clbk_add=f_loc(callback)
    else
       !find the callback in the error definition
       !these data can be inserted in a function
       dict_tmp=>f_get_error_dict(new_errcode)
       !that is how dict_keys function should be called      
       if (dict_size(dict_tmp) <= size(keys)) keys=dict_keys(dict_tmp)
       dict_tmp=>dict_tmp//trim(keys(1))

       if (has_key(dict_tmp,errclbk)) clbk_add=dict_tmp//errclbk
       if (has_key(dict_tmp,errclbkadd)) clbk_data_add=dict_tmp//errclbkadd
    end if

    call err_abort(clbk_add,clbk_data_add)
    
  end subroutine f_err_throw

  function f_get_error_dict(icode)
    implicit none
    integer, intent(in) :: icode
    type(dictionary), pointer :: f_get_error_dict
    !local variables

    f_get_error_dict=>dict_errors//icode
  end function f_get_error_dict

  !> identify id of lastr error occured
  function f_get_last_error(add_msg)
    implicit none
    character(len=*), optional :: add_msg
    integer :: f_get_last_error
    !local variables
    integer :: ierr
    ierr=dict_len(dict_present_error)-1
    if (ierr >= 0) then
       f_get_last_error=dict_present_error//ierr//errid
       if (present(add_msg)) add_msg=dict_present_error//ierr//'Additional Info'
    else
       f_get_last_error=0
       if (present(add_msg)) add_msg=repeat(' ',len(add_msg))
    end if
  end function f_get_last_error

  !> clean the dictionary of present errors
   subroutine f_err_clean()
    implicit none
    call dict_free(dict_present_error)
    call dict_init(dict_present_error)
  end subroutine f_err_clean

  !> activate the exception handling for all errors
  subroutine f_err_open_try()
    implicit none
    call f_err_set_callback(f_err_ignore)
  end subroutine f_err_open_try
  
  !> close the try environment
  subroutine f_err_close_try()
    implicit none
    call f_err_clean() !no errors anymore
    call f_err_unset_callback()
  end subroutine f_err_close_try
  
  function f_get_error_definitions()
    implicit none
    type(dictionary), pointer :: f_get_error_definitions
    f_get_error_definitions => dict_errors
  end function f_get_error_definitions

!!$end module error_handling
