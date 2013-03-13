!> @file
!!  Module defining a dictionary
!! @author Luigi Genovese
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Define a dictionary and its basic usage rules
module dictionaries
  implicit none

  private

  integer, parameter, public :: max_field_length = 256

  !> error codes
  integer, parameter :: DICT_SUCCESS=0
  integer, parameter :: DICT_KEY_ABSENT=1
  integer, parameter :: DICT_VALUE_ABSENT=2
  integer, parameter :: DICT_ITEM_NOT_VALID=3
  integer, parameter :: DICT_CONVERSION_ERROR=-1

  logical :: exceptions=.false.
  integer :: last_error = DICT_SUCCESS

  type, public :: storage
     integer :: item !< Id of the item associated to the list
     integer :: nitems !< No. of items in the list
     integer :: nelems !< No. of items in the dictionary
     character(len=max_field_length) :: key
     character(len=max_field_length) :: value
  end type storage

  !> structure of the dictionary element
  type, public :: dictionary
     type(storage) :: data
     type(dictionary), pointer :: parent,next,child,previous
  end type dictionary

  !> operators in the dictionary
!  interface operator(/)
! TEST     module procedure get_dict_ptr,get_item_ptr
!     module procedure get_child_ptr,get_list_ptr
!  end interface
  interface operator(//)
     module procedure get_child_ptr,get_list_ptr
  end interface
  interface assignment(=)
     module procedure get_value,get_integer,get_real,get_double,get_long,get_dict
  end interface
  interface pop
     module procedure pop_dict,pop_item,pop_last
  end interface

  interface set
     module procedure put_child,put_value,put_list,put_integer,put_real,put_double,put_long
  end interface

  public :: operator(//), assignment(=)
  public :: set,dict_init,dict_free,try,close_try,pop,append,prepend,find_key,length
  

contains

  !error handling routines
  subroutine try()
    implicit none

    exceptions=.true.
  end subroutine try

  subroutine close_try()
    implicit none

    exceptions=.false.
  end subroutine close_try

  pure function try_error() result(ierr)
    implicit none
    integer :: ierr
    ierr = last_error
  end function try_error

  function no_key(dict)
    implicit none
    type(dictionary), intent(in) :: dict
    logical :: no_key
    !TEST
    no_key=(len(trim(dict%data%key)) == 0 .and. dict%data%item == -1) .and. associated(dict%parent)
  end function no_key

  function no_value(dict)
    implicit none
    type(dictionary), intent(in) :: dict
    logical :: no_value

    no_value=len(trim(dict%data%value)) == 0 .and. .not. associated(dict%child)
  end function no_value

  subroutine check_key(dict)
    implicit none
    type(dictionary), intent(in) :: dict

    if (no_key(dict)) then
       if (exceptions) then
          last_error=DICT_KEY_ABSENT
       else
          write(*,*)'ERROR: key absent in dictionary'
          stop
       end if
    end if

  end subroutine check_key

  subroutine check_value(dict)
    implicit none
    type(dictionary), intent(in) :: dict

    if (no_value(dict)) then
       if (exceptions) then
          last_error=DICT_VALUE_ABSENT
       else
          write(*,*)'ERROR: value absent in dictionary'
          stop
       end if
    end if

  end subroutine check_value

  subroutine check_conversion(ierror)
    implicit none
    integer, intent(in) :: ierror
    if (ierror /= 0) then
       if (exceptions) then
          last_error=DICT_CONVERSION_ERROR
       else
          write(*,*)'ERROR: conversion error'
          stop
       end if
    end if

  end subroutine check_conversion

  subroutine set_item(dict,item)
    implicit none
    type(dictionary) :: dict
    integer, intent(in) :: item

    dict%data%item=item
    if (associated(dict%parent)) then
       dict%parent%data%nitems=dict%parent%data%nitems+1
       if (item+1 > dict%parent%data%nitems) then
          if (exceptions) then
             last_error=DICT_ITEM_NOT_VALID
          else
             write(*,*)'ERROR: item not valid',item,dict%parent%data%nitems
             stop
          end if
       end if
    end if

  end subroutine set_item

  recursive subroutine define_parent(dict,child)
    implicit none
    type(dictionary), target :: dict
    type(dictionary) :: child

    child%parent=>dict
    if (associated(child%next)) call define_parent(dict,child%next)
  end subroutine define_parent

  subroutine define_brother(brother,dict)
    implicit none
    type(dictionary), target :: brother
    type(dictionary) :: dict

    dict%previous=>brother
  end subroutine define_brother

  subroutine reset_next(next,dict)
    implicit none
    type(dictionary), target :: next
    type(dictionary) :: dict
    !local variables
    type(dictionary), pointer :: dict_all
    
    !do something only if needed
    if (.not. associated(dict%next,target=next)) then
       dict_all=>dict%next
       dict%next=>next
       deallocate(dict_all)
    end if

  end subroutine reset_next


  subroutine pop_dict(dict,key)
    implicit none
    type(dictionary), intent(inout), pointer :: dict 
    character(len=*), intent(in) :: key
    
    !check if we are at the first level
!TEST    if (associated(dict%parent)) then
       call pop_dict_(dict%child,key)
       !if it is the last the dictionary should be empty
       if (.not. associated(dict%parent) .and. .not. associated(dict%child)) then
          call dict_free(dict)
       end if

!TEST    else
!TEST       call pop_dict_(dict,key)
!TEST    end if
  contains
    !> Eliminate a key from a dictionary if it exists
    recursive subroutine pop_dict_(dict,key)
      implicit none
      type(dictionary), intent(inout), pointer :: dict 
      character(len=*), intent(in) :: key
      !local variables
      type(dictionary), pointer :: dict_first !<in case of first occurrence

      if (associated(dict)) then
         !follow the chain, stop at  first occurence
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
            if (exceptions) then
               last_error=DICT_KEY_ABSENT
            else
               write(*,*)'ERROR: key absent in dictionary'
               stop
            end if
         end if
      else
         if (exceptions) then
            last_error=DICT_KEY_ABSENT
         else
            write(*,*)'ERROR: key absent in dictionary'
            stop
         end if
      end if

    end subroutine pop_dict_
  end subroutine pop_dict

  function length(dict)
    implicit none
    type(dictionary), intent(in), pointer :: dict
    integer :: length
    
    if (associated(dict)) then
       length=dict%data%nitems
    else
       length=-1
    end if
  end function length

  subroutine pop_last(dict)
    implicit none
    type(dictionary), intent(inout), pointer :: dict 
    !local variables
    integer :: nitems

    nitems=dict%data%nitems
    if (nitems > 0) then
       call pop_item(dict,nitems-1)
    else
       if (exceptions) then
          last_error=DICT_ITEM_NOT_VALID
       else
          write(*,*)'ERROR: list empty, pop not possible'
          stop
       end if
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
      implicit none
      type(dictionary), intent(inout), pointer :: dict 
      integer, intent(in) :: item
      !local variables
      type(dictionary), pointer :: dict_first !<in case of first occurrence

      if (associated(dict)) then
!         print *,dict%data%item,trim(dict%data%key)
         !follow the chain, stop at  first occurence
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
            if (exceptions) then
               last_error=DICT_KEY_ABSENT
            else
               write(*,*)'ERROR: item absent in dictionary'
               stop
            end if
         end if
      else
         if (exceptions) then
            last_error=DICT_KEY_ABSENT
         else
            write(*,*)'ERROR: item absent in dictionary'
            stop
         end if
      end if

    end subroutine pop_item_
  end subroutine pop_item


  function get_ptr(dict,key) result (ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: ptr

    !if we are not at the topmost level check for the child
    if (associated(dict%parent)) then
       ptr=>get_child_ptr(dict,key)
    else
       ptr=>get_dict_ptr(dict,key)
    end if
  end function get_ptr
  
 !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, create it in the next chain 
  !! Key Must be already present 
  recursive function find_key(dict,key) result (dict_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: dict_ptr
    if (.not. associated(dict)) then
       nullify(dict_ptr)
    end if
!TEST 
if (.not. associated(dict%parent)) then
   dict_ptr =>  find_key(dict%child,key)
   return
end if

!    print *,'here',trim(key)
    !follow the chain, stop at  first occurence
    if (trim(dict%data%key) == trim(key)) then
       dict_ptr => dict
    else if (associated(dict%next)) then
       dict_ptr => find_key(dict%next,key)
    else 
       nullify(dict_ptr)
    end if

  end function find_key

  !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, create it in the next chain 
  !! Key Must be already present 
  recursive function get_dict_ptr(dict,key) result (dict_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: dict_ptr

!    print *,'here',trim(key)
    !follow the chain, stop at  first occurence
    if (trim(dict%data%key) == trim(key)) then
       dict_ptr => dict
    else if (associated(dict%next)) then
       dict_ptr => get_dict_ptr(dict%next,key)
    else if (no_key(dict)) then !this is useful for the first assignation
       call set_elem(dict,key)
       !call set_field(key,dict%data%key)
       dict_ptr => dict
    else
       call dict_init(dict%next)
       !call set_field(key,dict%next%data%key)
       call define_brother(dict,dict%next) !chain the list in both directions
       if (associated(dict%parent)) call define_parent(dict%parent,dict%next)
       call set_elem(dict%next,key)
       dict_ptr => dict%next
    end if

  end function get_dict_ptr

  !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, create it in the child chain
  function get_child_ptr(dict,key) result (subd_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: subd_ptr

    call check_key(dict)
    
    if (associated(dict%child)) then
       subd_ptr => get_dict_ptr(dict%child,key)
    else
       call dict_init(dict%child)
       !call set_field(key,dict%child%data%key)
       call define_parent(dict,dict%child)
       call set_elem(dict%child,key)
       subd_ptr => dict%child
    end if

  end function get_child_ptr

  !> Retrieve the pointer to the item of the list.
  !! If the list does not exists, create it in the child chain.
  !! If the list is too short, create it in the next chain
  recursive function get_item_ptr(dict,item) result (item_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    integer, intent(in) :: item
    type(dictionary), pointer :: item_ptr

    !follow the chain, stop at  first occurence
    if (dict%data%item == item) then
       item_ptr => dict
    else if (associated(dict%next)) then
       item_ptr => get_item_ptr(dict%next,item)
    else if (no_key(dict)) then
       call set_item(dict,item)
       item_ptr => dict
    else
       call dict_init(dict%next)
       call define_brother(dict,dict%next) !chain the list in both directions
       if (associated(dict%parent)) call define_parent(dict%parent,dict%next)
       call set_item(dict%next,item)
       item_ptr => dict%next
    end if

  end function get_item_ptr

  !> Retrieve the pointer to the item of the list.
  !! If the list does not exists, create it in the child chain.
  !! If the list is too short, create it in the next chain
  function get_list_ptr(dict,item) result (subd_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    integer, intent(in) :: item
    type(dictionary), pointer :: subd_ptr

    call check_key(dict)
    
    if (associated(dict%child)) then
       subd_ptr => get_item_ptr(dict%child,item)
    else
       call dict_init(dict%child)
       call define_parent(dict,dict%child)
       call set_item(dict%child,item)
       subd_ptr => dict%child
    end if

  end function get_list_ptr
!
  !> assign a child to the  dictionary
  recursive subroutine put_child(dict,subd)
    implicit none
    type(dictionary), pointer :: dict
    type(dictionary), pointer :: subd

    !TEST
!if the dictionary starts with a master tree, eliminate it and put the child
    if (.not. associated(subd%parent)) then
       call put_child(dict,subd%child)
       nullify(subd%child)
       call dict_free(subd)
       return
    end if

    call check_key(dict)

    call set_field(repeat(' ',max_field_length),dict%data%value)
    if ( .not. associated(dict%child,target=subd) .and. &
         associated(dict%child)) then
       call dict_free(dict%child)
    end if
    dict%child=>subd
    call define_parent(dict,dict%child)
    dict%data%nelems=dict%data%nelems+1
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
       dict%data%nelems=dict%data%nelems+brother%data%nelems
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

  !> assign the value to the  dictionary
  subroutine put_value(dict,val)
    implicit none
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: val

    call check_key(dict)

    if (associated(dict%child)) call dict_free(dict%child)

    call set_field(val,dict%data%value)

  end subroutine put_value

  !> assign the value to the  dictionary (to be rewritten)
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

  !> get the value from the  dictionary
  subroutine get_value(val,dict)
    implicit none
    character(len=*), intent(out) :: val
    type(dictionary), intent(in) :: dict

    call check_key(dict)
    call check_value(dict)

    call get_field(dict%data%value,val)

  end subroutine get_value

  !> get the value from the  dictionary
  !! This routine only works if the dictionary is associated
  !! the problem is solved if any of the routines have the dict variable as a pointer
  subroutine get_dict(dictval,dict)
    implicit none
    type(dictionary), pointer, intent(out) :: dictval
    type(dictionary), target, intent(in) :: dict

    call check_key(dict)
    
    !if (associated(dict%child)) then
    !   dictval=>dict%child
!    if (associated(dict)) then
       dictval=>dict
!    else
!       nullify(dictval)
!    end if

  end subroutine get_dict


  pure subroutine dictionary_nullify(dict)
    implicit none
    type(dictionary), intent(inout) :: dict

    dict%data%key=repeat(' ',max_field_length)
    dict%data%value=repeat(' ',max_field_length)
    dict%data%item=-1
    dict%data%nitems=0
    dict%data%nelems=0
    nullify(dict%child,dict%next,dict%parent,dict%previous)
  end subroutine dictionary_nullify

  pure subroutine dict_free(dict)
    type(dictionary), pointer :: dict

    if (associated(dict)) then
       call dict_free_(dict)
       deallocate(dict)
    end if

  contains

    pure recursive subroutine dict_free_(dict)
      implicit none
      type(dictionary), pointer :: dict

      !first destroy the children
      if (associated(dict%child)) then
         call dict_free_(dict%child)
         deallocate(dict%child)
         nullify(dict%child)
      end if
      !then destroy younger brothers
      if (associated(dict%next)) then
         call dict_free_(dict%next)
         deallocate(dict%next)
         nullify(dict%next)
      end if
      call dictionary_nullify(dict)

    end subroutine dict_free_

  end subroutine dict_free

  pure subroutine dict_init(dict)
    implicit none
    type(dictionary), pointer :: dict

    allocate(dict)
    call dictionary_nullify(dict)

  end subroutine dict_init
  
  pure subroutine set_elem(dict,key)
    implicit none
    type(dictionary), pointer :: dict !!TO BE VERIFIED
    character(len=*), intent(in) :: key

    call set_field(trim(key),dict%data%key)
    if (associated(dict%parent)) then
       dict%parent%data%nelems=dict%parent%data%nelems+1
    else
       dict%data%nelems=dict%data%nelems+1
    end if

  end subroutine set_elem

  pure subroutine set_field(input,output)
    implicit none
    character(len=*), intent(in) :: input !intent eliminated
    character(len=max_field_length), intent(out) :: output !intent eliminated
    !local variables
    integer :: ipos,i

    ipos=min(len(trim(input)),max_field_length)
    do i=1,ipos
       output(i:i)=input(i:i)
    end do
    do i=ipos+1,max_field_length
       output(i:i)=' ' 
    end do

  end subroutine set_field

  pure subroutine get_field(input,output)
    implicit none
    character(len=max_field_length), intent(in) :: input !intent eliminated
    character(len=*), intent(out) :: output !intent eliminated
    !local variables
    integer :: ipos,i

    ipos=min(len(output),max_field_length)
    do i=1,ipos
       output(i:i)=input(i:i)
    end do
    do i=ipos+1,len(output)
       output(i:i)=' ' 
    end do

  end subroutine get_field

  recursive function dict_list_size(dict) result(ntot)
    implicit none
    type(dictionary), intent(in) :: dict
    integer :: ntot
    !item
    integer :: npos=0
    
    ntot=-1
    !print *,field_to_integer(dict%key),trim(dict%key),npos
    if (associated(dict%parent)) npos=0 !beginning of the list
    if (field_to_integer(dict%data%key) == npos) then
       npos=npos+1
       ntot=npos
    else
       npos=0
       return
    end if

    if (associated(dict%next)) then
       ntot=dict_list_size(dict%next)
    end if
  end function dict_list_size
 
  function field_to_integer(input)
    implicit none
    character(len=max_field_length), intent(in) :: input
    integer :: field_to_integer
    !local variables
    integer :: iprobe,ierror

    !look at conversion
    read(input,*,iostat=ierror)iprobe
    !print *,trim(input),'test',ierror
    if (ierror /=0) then
       field_to_integer=-1
    else
       field_to_integer=iprobe
    end if

  end function field_to_integer

  !set and get routines for different types
  subroutine get_integer(ival,dict)
    integer, intent(out) :: ival
    type(dictionary), intent(in) :: dict
    !local variables
    integer :: ierror
    character(len=max_field_length) :: val

    !take value
    val=dict
    !look at conversion
    read(val,*,iostat=ierror)ival

    call check_conversion(ierror)
    
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

    call check_conversion(ierror)
    
  end subroutine get_long


  !set and get routines for different types
  subroutine get_real(rval,dict)
    real(kind=4), intent(out) :: rval
    type(dictionary), intent(in) :: dict
    !local variables
    integer :: ierror
    character(len=max_field_length) :: val

    !take value
    val=dict
    !look at conversion
    read(val,*,iostat=ierror)rval

    call check_conversion(ierror)
    
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
    read(val,*,iostat=ierror)dval

    call check_conversion(ierror)
    
  end subroutine get_double


  !> assign the value to the  dictionary
  subroutine put_integer(dict,ival,fmt)
    use yaml_strings, only:yaml_toa
    implicit none
    type(dictionary), pointer :: dict
    integer, intent(in) :: ival
    character(len=*), optional, intent(in) :: fmt

    if (present(fmt)) then
       call put_value(dict,adjustl(trim(yaml_toa(ival,fmt=fmt))))
    else
       call put_value(dict,adjustl(trim(yaml_toa(ival))))
    end if

  end subroutine put_integer

  !> assign the value to the  dictionary
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

  !> assign the value to the  dictionary
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

  !> assign the value to the  dictionary
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


end module dictionaries
