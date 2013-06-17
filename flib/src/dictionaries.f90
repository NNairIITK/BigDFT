!> @file
!!  Module defining a dictionary
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Define a dictionary and its basic usage rules
module dictionaries
  use exception_callbacks
  use dictionaries_base
!  use error_handling
  implicit none

  private

  !>public to be used in list_new() constructor.
  type, public :: list_container
     character(len=max_field_length) :: val=' '
     type(dictionary), pointer :: dict => null()
  end type list_container
  !>public to be used in dict_new() constructor.
  type, public :: dictionary_container
     character(len=max_field_length) :: key=' '
     character(len=max_field_length) :: value=' '
    type(dictionary), pointer :: child => null()
  end type dictionary_container


  !> Error codes
  integer, public :: DICT_KEY_ABSENT
  integer, public :: DICT_VALUE_ABSENT
  integer, public :: DICT_ITEM_NOT_VALID
  integer, public :: DICT_CONVERSION_ERROR
  integer, public :: DICT_INVALID_LIST
  integer, public :: DICT_INVALID

  interface operator(.index.)
     module procedure find_index
  end interface
  interface operator(.item.)
     module procedure item_char,item_dict
  end interface
  interface operator(.is.)
     module procedure dict_cont_new_with_value, dict_cont_new_with_dict
  end interface

  interface assignment(=)
     module procedure get_value,get_integer,get_real,get_double,get_long,get_bool!,get_dict, xlf does not like
  end interface
  interface pop
     module procedure pop_dict,pop_item,pop_last
  end interface

  interface set
     module procedure put_child,put_value,put_list,put_integer,put_real,put_double,put_long,put_bool
  end interface

  interface add
     module procedure add_char,add_dict,add_integer,add_real,add_double,add_long
  end interface

  interface dict_new
     module procedure dict_new,dict_new_single
  end interface

  integer(kind=8), external :: f_loc

  !> Public routines
  public :: operator(//),operator(.index.),assignment(=)
  public :: set,dict_init,dict_free,pop,append,prepend,add
  !Handle exceptions
  public :: find_key,dict_len,dict_size,dict_key,dict_next,has_key,dict_keys
  public :: dict_new,list_new
  !> Public elements of dictionary_base
  public :: operator(.is.),operator(.item.)
  public :: dictionary,max_field_length

  !header of error handling part
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

  public :: f_err_initialize,f_err_finalize,f_get_last_error,f_get_error_definitions
  public :: f_err_define,f_err_check,f_err_raise,f_err_clean,f_get_error_dict,f_err_throw

  !public variables of the callback module
  public :: f_err_set_callback,f_err_unset_callback,f_err_open_try,f_err_close_try
  public :: f_err_severe,f_err_severe_override,f_err_severe_restore,f_err_ignore
  public :: f_loc


contains

  !>define the errors of the dictionary module
  subroutine dictionary_errors()
    implicit none
    call f_err_define('DICT_KEY_ABSENT',&
         'The dictionary has no key',DICT_KEY_ABSENT,&
         err_action='Internal error, contact developers')
    call f_err_define('DICT_ITEM_NOT_VALID',&
         'The item of this list is not correct',DICT_ITEM_NOT_VALID,&
         err_action='Internal error, contact developers')
    call f_err_define('DICT_VALUE_ABSENT',&
         'The value for this key is absent',DICT_VALUE_ABSENT)
    call f_err_define('DICT_INVALID',&
         'Dictionary is not associated',DICT_INVALID)
    call f_err_define('DICT_INVALID_LIST',&
         'Current node is not a list',DICT_INVALID_LIST)
    call f_err_define('DICT_CONVERSION_ERROR',&
         'Conversion error of the dictionary value',DICT_CONVERSION_ERROR,&
         err_action='Check the nature of the conversion')

  end subroutine dictionary_errors

!!$  !> seems unused routine
!!$  subroutine reset_next(next,dict)
!!$    implicit none
!!$    type(dictionary), target :: next
!!$    type(dictionary) :: dict
!!$    !local variables
!!$    type(dictionary), pointer :: dict_all
!!$    
!!$    !do something only if needed
!!$    if (.not. associated(dict%next,target=next)) then
!!$       dict_all=>dict%next
!!$       dict%next=>next
!!$       deallocate(dict_all)
!!$    end if
!!$
!!$  end subroutine reset_next

  subroutine pop_dict(dict,key)
    implicit none
    type(dictionary), intent(inout), pointer :: dict 
    character(len=*), intent(in) :: key
    
    !check if we are at the first level
    call pop_dict_(dict%child,key)
    !if it is the last the dictionary should be empty
    if (.not. associated(dict%parent) .and. .not. associated(dict%child)) then
       call dict_free(dict)
    end if

  contains
    !> Eliminate a key from a dictionary if it exists
    recursive subroutine pop_dict_(dict,key)
      implicit none
      type(dictionary), intent(inout), pointer :: dict 
      character(len=*), intent(in) :: key
      !local variables
      type(dictionary), pointer :: dict_first !<in case of first occurrence

      if (associated(dict)) then
         !follow the chain, stop at the first occurence
         if (trim(dict%data%key) == trim(key)) then
            !          print *,'here',trim(key),associated(dict%next)
            if (associated(dict%parent)) then
               dict%parent%data%nelems=dict%parent%data%nelems-1
            else
               dict%data%nelems=dict%data%nelems-1
            end if
            if (associated(dict%next)) then
               call dict_free(dict%child)
               dict_first => dict
               !this is valid if we are not at the first element
               if (associated(dict%previous)) then
                  call define_brother(dict%previous,dict%next) 
                  dict%previous%next => dict%next
               else
                  nullify(dict%next%previous)
                  !the next should now become me
                  dict => dict%next
               end if
               deallocate(dict_first)
            else
               call dict_free(dict)
            end if
         else if (associated(dict%next)) then
            call pop_dict_(dict%next,key)
         else
            if (f_err_raise(err_msg='Key is '//trim(key),&
                 err_id=DICT_KEY_ABSENT)) return
         end if
      else
         if (f_err_raise(err_msg='Key is '//trim(key),&
              err_id=DICT_KEY_ABSENT)) return
      end if

    end subroutine pop_dict_
  end subroutine pop_dict

  !> add to a list
  subroutine add_char(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: val
    include 'dict_add-inc.f90'
  end subroutine add_char
  subroutine add_dict(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    type(dictionary), pointer :: val
    include 'dict_add-inc.f90'
  end subroutine add_dict
  subroutine add_integer(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    integer, intent(in) :: val
    include 'dict_add-inc.f90'
  end subroutine add_integer
  subroutine add_real(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    real, intent(in) :: val
    include 'dict_add-inc.f90'
  end subroutine add_real
  subroutine add_double(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    double precision, intent(in) :: val
    include 'dict_add-inc.f90'
  end subroutine add_double
  subroutine add_long(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    integer(kind=8), intent(in) :: val
    include 'dict_add-inc.f90'
  end subroutine add_long

  !defines a dictionary from a array of storage data
  function dict_new(dicts)
!    use yaml_output
!    type(storage), dimension(:), intent(in) :: st_arr
    type(dictionary_container), dimension(:), intent(in) :: dicts
    type(dictionary), pointer :: dict_new
    !local variables
    integer :: i_st,n_st
    type(dictionary), pointer :: dict_tmp

    !initialize dictionary
    call dict_init(dict_tmp)
    n_st=size(dicts)
    do i_st=1,n_st
       if (associated(dicts(i_st)%child)) then
          call set(dict_tmp//dicts(i_st)%key, dicts(i_st)%child)
       else
          call set(dict_tmp//dicts(i_st)%key, dicts(i_st)%value)
       end if
    end do
    dict_new => dict_tmp
  end function dict_new

  !defines a dictionary from a array of storage data
  function dict_new_single(cont)
    type(dictionary_container), intent(in), optional :: cont
    type(dictionary), pointer :: dict_new_single
    !local variables
    type(dictionary), pointer :: dict_tmp

    call dict_init(dict_tmp)
    if (present(cont)) then
       if (associated(cont%child)) then
          call set(dict_tmp//cont%key, cont%child)
       else
          call set(dict_tmp//cont%key, cont%value)
       end if
    end if
    dict_new_single => dict_tmp

  end function dict_new_single

  !defines a new dictionary from a key and a value
  function dict_cont_new_with_value(key, val)
    character(len = *), intent(in) :: key, val
    type(dictionary_container) :: dict_cont_new_with_value

    dict_cont_new_with_value%key(1:max_field_length) = key
    dict_cont_new_with_value%value(1:max_field_length) = val

  end function dict_cont_new_with_value

  function dict_cont_new_with_dict(key, val)
    character(len = *), intent(in) :: key
    type(dictionary), pointer, intent(in) :: val
    type(dictionary_container) :: dict_cont_new_with_dict

    dict_cont_new_with_dict%key(1:max_field_length) = key
    dict_cont_new_with_dict%child => val

  end function dict_cont_new_with_dict

  function dict_next(dict)
    implicit none
    type(dictionary), pointer, intent(in) :: dict
    type(dictionary), pointer :: dict_next
    
    if (associated(dict)) then
       if (associated(dict%parent)) then
          dict_next=>dict%next
       else
          dict_next=>dict%child
       end if
    else
       nullify(dict_next)
    end if
  end function dict_next

  !> Returns the position of the name in the dictionary
  !! returns 0 if the dictionary is nullified or the name is absent
  function find_index(dict,name)
    implicit none
    type(dictionary), pointer, intent(in) :: dict
    character(len=*), intent(in) :: name
    integer :: find_index
    !local variables
    integer :: ind
    type(dictionary), pointer :: dict_tmp

    find_index =0
    ind=-1
    if (associated(dict)) then
       dict_tmp=>dict_next(dict)
       loop_find: do while(associated(dict_tmp))
          ind=ind+1
          if (name_is(dict_tmp,name)) then
             find_index=ind
             exit loop_find
          end if
          dict_tmp=>dict_next(dict_tmp)
       end do loop_find
    end if

  end function find_index

  subroutine pop_last(dict)
    implicit none
    type(dictionary), intent(inout), pointer :: dict 
    !local variables
    integer :: nitems

    nitems=dict_len(dict)

    if (f_err_raise(nitems <= 0,'Pop not allowed for this node',&
         err_id=DICT_ITEM_NOT_VALID)) then
       return
    else
       call pop_item(dict,nitems-1)
    end if

  end subroutine pop_last
  
  subroutine pop_item(dict,item)
    implicit none
    type(dictionary), intent(inout), pointer :: dict 
    integer, intent(in) :: item

    !check if we are at the first level
 !TEST   if (associated(dict%parent)) then
       call pop_item_(dict%child,item)
       !if it is the last the dictionary should be empty
       if (.not. associated(dict%parent) .and. .not. associated(dict%child)) then
          call dict_free(dict)
       end if

!TEST    else
!TEST       call pop_item_(dict,item)
!TEST    end if
  contains
    !> Eliminate a key from a dictionary if it exists
    recursive subroutine pop_item_(dict,item)
      use yaml_strings, only: yaml_toa
      implicit none
      type(dictionary), intent(inout), pointer :: dict 
      integer, intent(in) :: item
      !local variables
      type(dictionary), pointer :: dict_first !<in case of first occurrence

      if (associated(dict)) then
!         print *,dict%data%item,trim(dict%data%key)
         !follow the chain, stop at the first occurence
         if (dict%data%item == item) then
            if (associated(dict%parent)) then
               dict%parent%data%nitems=dict%parent%data%nitems-1
            end if
            if (associated(dict%next)) then
               call dict_free(dict%child)
               dict_first => dict
               !this is valid if we are not at the first element
               if (associated(dict%previous)) then
                  call define_brother(dict%previous,dict%next) 
                  dict%previous%next => dict%next
               else
                  !the next should now become me
                  dict => dict%next
               end if
               deallocate(dict_first)
            else
               call dict_free(dict)
            end if
         else if (associated(dict%next)) then
            call pop_item_(dict%next,item)
         else
            if (f_err_raise(err_msg='Item No. '//trim(yaml_toa(item)),&
                 err_id=DICT_ITEM_NOT_VALID)) return
         end if
      else
         if (f_err_raise(err_msg='Item No. '//trim(yaml_toa(item)),&
              err_id=DICT_ITEM_NOT_VALID)) return
      end if

    end subroutine pop_item_
  end subroutine pop_item
 
  !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, search for it in the next chain 
  !! Key Must be already present 
  recursive function find_key(dict,key) result (dict_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: dict_ptr
    if (.not. associated(dict)) then
       nullify(dict_ptr)
       return
    end if
    !TEST 
    if (.not. associated(dict%parent)) then
       dict_ptr => find_key(dict%child,key)
       !print *,'parent'
       return
    end if

    !print *,'here ',trim(key),', key ',trim(dict%data%key)
    !follow the chain, stop at the first occurence
    if (trim(dict%data%key) == trim(key)) then
       dict_ptr => dict
    else if (associated(dict%next)) then
       dict_ptr => find_key(dict%next,key)
    else 
       nullify(dict_ptr)
    end if

  end function find_key

  function dict_keys(dict)
    implicit none
    type(dictionary), intent(in) :: dict !<the dictionary must be associated
    character(len=max_field_length), dimension(dict%data%nelems) :: dict_keys
    !local variables
    integer :: ikey
    type(dictionary), pointer :: dict_tmp

    !if (associated(dict)) then
       ikey=0
       dict_tmp=>dict%child
       do while(associated(dict_tmp))
          ikey=ikey+1
          dict_keys(ikey)=dict_key(dict_tmp)
          dict_tmp=>dict_tmp%next
       end do
    !end if

  end function dict_keys

  !> Search in the dictionary if some of the child has the given
  !! If the key does not exists, search for it in the next chain 
  !! Key Must be already present 
  !! 
  function has_key(dict,key)
    implicit none
    type(dictionary), intent(in), pointer :: dict 
    character(len=*), intent(in) :: key
    logical :: has_key

    if (.not. associated(dict) .or. trim(key) == "") then
       has_key=.false.
       return
    end if

    has_key=has_key_(dict%child,key)
  
  contains

    recursive function has_key_(dict,key) result(has)
      implicit none
      type(dictionary), intent(in), pointer :: dict 
      character(len=*), intent(in) :: key
      logical :: has
      if (.not. associated(dict)) then
         has=.false.
         return
      end if

      !print *,'here ',trim(key),', key ',trim(dict%data%key)
      !follow the chain, stop at the first occurence
      if (trim(dict%data%key) == trim(key)) then
         has=.true.
      else if (associated(dict%next)) then
         has=has_key_(dict%next,key)
      else 
         has=.false.
      end if

    end function has_key_
  end function has_key

  !> assign a child to the dictionary
  recursive subroutine put_child(dict,subd)
    implicit none
    type(dictionary), pointer :: dict
    type(dictionary), pointer :: subd

    !if the dictionary starts with a master tree, eliminate it and put the child
    if (.not. associated(subd%parent)) then
       call put_child(dict,subd%child)
       nullify(subd%child)
       call dict_free(subd)
       return
    end if

    if (f_err_raise(no_key(dict),err_id=DICT_KEY_ABSENT)) return
    !call check_key(dict)

    call set_field(repeat(' ',max_field_length),dict%data%value)
    if ( .not. associated(dict%child,target=subd) .and. &
         associated(dict%child)) then
       call dict_free(dict%child)
    end if
    dict%child=>subd
    !inherit the number of elements or items from subd's parent
    !which is guaranteed to be associated
    dict%data%nelems=subd%parent%data%nelems
    dict%data%nitems=subd%parent%data%nitems
    call define_parent(dict,dict%child)

  end subroutine put_child

  !> append another dictionary
  recursive subroutine append(dict,brother)
    implicit none
    type(dictionary), pointer :: dict
    type(dictionary), pointer :: brother

    if (.not. associated(dict)) then
       !this should be verifyed by passing a dictionary which is not in the beginning
       if (associated(brother%parent)) then
          call dict_init(dict)
          call set(dict,brother)
       else
          dict=>brother
       end if
    else if (.not. associated(dict%parent)) then
       call append(dict%child,brother)
    else if (.not. associated(brother%parent)) then
       call append(dict,brother%child)
       nullify(brother%child)
       call dict_free(brother)
    else if (associated(dict%next)) then
       call append(dict%next,brother)
    else
       dict%next=>brother
       call define_brother(dict,dict%next)
       dict%parent%data%nelems=dict%parent%data%nelems+brother%parent%data%nelems
       !print *,'appending',dict%parent%data%nelems,dict%data%nelems,brother%data%nelems,brother%parent%data%nelems
    end if
  end subroutine append

  !> append another dictionary
  recursive subroutine prepend(dict,brother)
    implicit none
    type(dictionary), pointer :: dict
    type(dictionary), pointer :: brother
    !local variables
    type(dictionary), pointer :: dict_tmp

    if (.not. associated(brother)) return

    if (.not. associated(dict)) then
       if (associated(brother%parent)) then
          call dict_init(dict)
          call set(dict,brother)
       else
          dict=>brother
       end if
    else if (.not. associated(dict%parent)) then
       call prepend(dict%child,brother)
    else if (.not. associated(brother%parent)) then
       call define_parent(dict%parent,brother%child)
       call prepend(dict,brother%child)
       nullify(brother%child)
       call dict_free(brother)
    else if (associated(dict%previous)) then
       call prepend(dict%previous,brother)
    else
       dict_tmp=>brother
       call append(brother,dict)
       dict=>dict_tmp
    end if
  end subroutine prepend

  !> assign the value to the dictionary
  subroutine put_value(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: val

    if (f_err_raise(no_key(dict),err_id=DICT_KEY_ABSENT)) return
    !call check_key(dict)

    if (associated(dict%child)) call dict_free(dict%child)

    call set_field(val,dict%data%value)

  end subroutine put_value
  
  !> assign the value to the dictionary (to be rewritten)
  subroutine put_list(dict,list)!,nitems)
    implicit none
    type(dictionary), pointer :: dict
!    integer, intent(in) :: nitems
    character(len=*), dimension(:), intent(in) :: list
    !local variables
    integer :: item,nitems

    nitems=size(list)
    do item=1,nitems
       call set(dict//(item-1),list(item))
    end do

  end subroutine put_list
  
  !> creates a dictionary which has only one entry as a list
  function item_char(val) result(elem)
    implicit none
    character(len=*), intent(in) :: val
    type(list_container) :: elem

    elem%val(1:max_field_length)=val

  end function item_char

  function item_dict(val) result(elem)
    implicit none
    type(dictionary), pointer, intent(in) :: val
    type(list_container) :: elem
    
    elem%dict=>val
  end function item_dict

  !creates a list from a table of dictionaries
  function list_new(dicts)
    implicit none
    type(list_container), dimension(:) :: dicts
    type(dictionary), pointer :: list_new
    !local variables
    integer :: i_st,n_st
    type(dictionary), pointer :: dict_tmp

    !initialize dictionary
    call dict_init(dict_tmp)
    
    n_st=size(dicts)
    do i_st=1,n_st
       if (associated(dicts(i_st)%dict)) then
          call add(dict_tmp,dicts(i_st)%dict)
       else if (len_trim(dicts(i_st)%val) > 0) then
          call add(dict_tmp,dicts(i_st)%val)
       end if
    end do

    list_new => dict_tmp

  end function list_new

  !> get the value from the dictionary
  subroutine get_value(val,dict)
    implicit none
    character(len=*), intent(out) :: val
    type(dictionary), intent(in) :: dict
    val(1:len(val))=' '
    if (f_err_raise(no_key(dict),err_id=DICT_KEY_ABSENT)) return
    if (f_err_raise(no_value(dict),err_id=DICT_VALUE_ABSENT)) return

    call get_field(dict%data%value,val)

  end subroutine get_value

  !> get the value from the dictionary
  !! This routine only works if the dictionary is associated
  !! the problem is solved if any of the routines have the dict variable as a pointer
  subroutine get_dict(dictval,dict)
    implicit none
    type(dictionary), pointer, intent(out) :: dictval
    type(dictionary), pointer, intent(in) :: dict

    dictval=>dict

  end subroutine get_dict

  !set and get routines for different types
  subroutine get_bool(lval,dict)
    logical, intent(out) :: lval
    type(dictionary), intent(in) :: dict
    !local variables
    integer :: ierror
    character(len=max_field_length) :: val

    !take value
    val=dict
    ierror = 0
    if (trim(val) == "Yes") then
       lval = .true.
    else if (trim(val) == "No") then
       lval = .false.
    else
       !look at conversion
       read(val,*,iostat=ierror)lval
    end if

    if (f_err_raise(ierror/=0,'Value '//val,err_id=DICT_CONVERSION_ERROR)) return    
  end subroutine get_bool

  !set and get routines for different types (this routine can be called from error_check also)
  recursive subroutine get_integer(ival,dict)
    integer, intent(out) :: ival
    type(dictionary), intent(in) :: dict
    !local variables
    integer :: ierror
    character(len=max_field_length) :: val

    !take value
    val=dict
    !look at conversion
    read(val,*,iostat=ierror)ival
    !is the value existing?
    if (ierror/=0) then
       if (f_err_check(err_id=DICT_VALUE_ABSENT))then
          ival=0
          return
       end if
    end if
    if (f_err_raise(ierror/=0,'Value '//val,err_id=DICT_CONVERSION_ERROR)) return    
  end subroutine get_integer

  !set and get routines for different types
  subroutine get_long(ival,dict)
    integer(kind=8), intent(out) :: ival
    type(dictionary), intent(in) :: dict
    !local variables
    integer :: ierror
    character(len=max_field_length) :: val

    !take value
    val=dict
    !look at conversion
    read(val,*,iostat=ierror)ival

    if (f_err_raise(ierror/=0,'Value '//val,err_id=DICT_CONVERSION_ERROR)) return
    
  end subroutine get_long

  !> Read a real or real/real, real:real 
  !! Here the fraction is indicated by the ':' or '/'
  !! The problem is that / is a separator for Fortran
  subroutine read_fraction_string(string,var,ierror)
    implicit none
    !Arguments
    character(len=*), intent(in) :: string
    double precision, intent(out) :: var
    integer, intent(out) :: ierror
    !Local variables
    character(len=max_field_length) :: tmp
    integer :: num,den,pfr,psp

    !First look at the first blank after trim
    tmp(1:max_field_length)=trim(string)
    psp = scan(tmp,' ')

    !see whether there is a fraction in the string
    if(psp==0) psp=len(tmp)
    pfr = scan(tmp(1:psp),':')
    if (pfr == 0) pfr = scan(tmp(1:psp),'/')
    !It is not a fraction
    if (pfr == 0) then
       read(tmp(1:psp),*,iostat=ierror) var
    else 
       read(tmp(1:pfr-1),*,iostat=ierror) num
       read(tmp(pfr+1:psp),*,iostat=ierror) den
       if (ierror == 0) var=real(num)/real(den)
    end if
    !Value by defaut
    if (ierror /= 0) var = huge(1.d0) 
  END SUBROUTINE read_fraction_string

  !set and get routines for different types
  subroutine get_real(rval,dict)
    real(kind=4), intent(out) :: rval
    type(dictionary), intent(in) :: dict
    !local variables
    integer :: ierror
    double precision :: dval
    character(len=max_field_length) :: val

    !take value
    val=dict
    !look at conversion
    call read_fraction_string(val, dval, ierror)
    rval = dval

    if (f_err_raise(ierror/=0,'Value '//val,err_id=DICT_CONVERSION_ERROR)) return
    
  end subroutine get_real

  !set and get routines for different types
  subroutine get_double(dval,dict)
    real(kind=8), intent(out) :: dval
    type(dictionary), intent(in) :: dict
    !local variables
    integer :: ierror
    character(len=max_field_length) :: val

    !take value
    val=dict
    !look at conversion
    call read_fraction_string(val, dval, ierror)

    if (f_err_raise(ierror/=0,'Value '//val,err_id=DICT_CONVERSION_ERROR)) return
     
  end subroutine get_double

  !> assign the value to the dictionary
  subroutine put_bool(dict,lval,fmt)
    use yaml_strings, only:yaml_toa
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(in) :: lval
    character(len=*), optional, intent(in) :: fmt

    if (present(fmt)) then
       call put_value(dict,trim(adjustl(yaml_toa(lval,fmt=fmt))))
    else
       call put_value(dict,trim(adjustl(yaml_toa(lval))))
    end if
  end subroutine put_bool

  !> assign the value to the dictionary
  subroutine put_integer(dict,ival,fmt)
    use yaml_strings, only:yaml_toa
    implicit none
    type(dictionary), pointer :: dict
    integer, intent(in) :: ival
    character(len=*), optional, intent(in) :: fmt

    if (present(fmt)) then
       call put_value(dict,trim(adjustl(yaml_toa(ival,fmt=fmt))))
    else
       call put_value(dict,trim(adjustl(yaml_toa(ival))))
    end if

  end subroutine put_integer

  !> assign the value to the dictionary
  subroutine put_double(dict,dval,fmt)
    use yaml_strings, only:yaml_toa
    implicit none
    type(dictionary), pointer :: dict
    real(kind=8), intent(in) :: dval
    character(len=*), optional, intent(in) :: fmt

    if (present(fmt)) then
       call put_value(dict,adjustl(trim(yaml_toa(dval,fmt=fmt))))
    else
       call put_value(dict,adjustl(trim(yaml_toa(dval))))
    end if

  end subroutine put_double

  !> assign the value to the dictionary
  subroutine put_real(dict,rval,fmt)
    use yaml_strings, only:yaml_toa
    implicit none
    type(dictionary), pointer :: dict
    real(kind=4), intent(in) :: rval
    character(len=*), optional, intent(in) :: fmt

    if (present(fmt)) then
       call put_value(dict,adjustl(trim(yaml_toa(rval,fmt=fmt))))
    else
       call put_value(dict,adjustl(trim(yaml_toa(rval))))
    end if

  end subroutine put_real

  !> assign the value to the dictionary
  subroutine put_long(dict,ilval,fmt)
    use yaml_strings, only:yaml_toa
    implicit none
    type(dictionary), pointer :: dict
    integer(kind=8), intent(in) :: ilval
    character(len=*), optional, intent(in) :: fmt

    if (present(fmt)) then
       call put_value(dict,adjustl(trim(yaml_toa(ilval,fmt=fmt))))
    else
       call put_value(dict,adjustl(trim(yaml_toa(ilval))))
    end if

  end subroutine put_long

  !include the module of error handling
  include 'error_handling.f90'

end module dictionaries

