!> @file
!!  Module defining a dictionary
!! @ingroup flib
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module which defines a dictionary and the pure functions for its basic usage rules (no dependency)
module dictionaries_base

  implicit none

  integer, parameter, public :: max_field_length = 256
  character(len=max_field_length), parameter, public :: TYPE_DICT='__dict__'
  character(len=max_field_length), parameter, public :: TYPE_LIST='__list__'

  type, public :: storage
     integer :: item   !< Id of the item associated to the list
     integer :: nitems !< No. of items in the list
     integer :: nelems !< No. of items in the dictionary
     character(len=max_field_length) :: key
     character(len=max_field_length) :: value
  end type storage

  !> structure of the dictionary element (this internal structure is private)
  type, public :: dictionary
     type(storage) :: data
     type(dictionary), pointer :: parent,next,child,previous
  end type dictionary

  !> operator to access and create a key in the dictionary
  interface operator(//)
     module procedure get_child_ptr,get_list_ptr
  end interface
  interface operator(.is.)
     module procedure storage_data
  end interface

contains
  
  !> Test if keys are present
  pure function no_key(dict)
    implicit none
    type(dictionary), intent(in) :: dict
    logical :: no_key
    !TEST
    no_key=(len_trim(dict%data%key) == 0 .and. dict%data%item == -1) .and. &
         associated(dict%parent)
  end function no_key

  pure function no_value(dict)
    implicit none
    type(dictionary), intent(in) :: dict
    logical :: no_value

    no_value=len(trim(dict%data%value)) == 0 .and. .not. associated(dict%child)
  end function no_value

  pure function storage_null() result(st)
    type(storage) :: st
    st%key=repeat(' ',max_field_length)
    st%value=repeat(' ',max_field_length)
    st%item=-1
    st%nitems=0
    st%nelems=0
  end function storage_null

  pure subroutine dictionary_nullify(dict)
    implicit none
    type(dictionary), intent(inout) :: dict
    dict%data=storage_null()
    nullify(dict%child,dict%next,dict%parent,dict%previous)
  end subroutine dictionary_nullify

  pure subroutine dict_init(dict)
    implicit none
    type(dictionary), pointer :: dict

    allocate(dict)
    call dictionary_nullify(dict)

  end subroutine dict_init

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

  !> return the length of the list
  pure function dict_len(dict)
    implicit none
    type(dictionary), intent(in), pointer :: dict
    integer :: dict_len
    
    if (associated(dict)) then
       dict_len=dict%data%nitems
    else
       dict_len=-1
    end if
  end function dict_len

  !> return the size of the dictionary
  pure function dict_size(dict)
    implicit none
    type(dictionary), intent(in), pointer :: dict
    integer :: dict_size
    
    if (associated(dict)) then
       dict_size=dict%data%nelems
    else
       dict_size=-1
    end if
  end function dict_size

  !> this function returns the key if present otherwise the value of the element if in a list
  pure function name_is(dict,name)
    implicit none
    type(dictionary), pointer, intent(in) :: dict
    character(len=*), intent(in) :: name
    logical :: name_is

!!$    print *,'name is',trim(name),no_key(dict%child),no_value(dict%child),&
!!$         'value',trim(dict%data%value),'key',trim(dict%data%key),&
!!$         'item',dict%data%item,&
!!$         'child',associated(dict%child),'valueC',trim(dict%child%data%value),&
!!$         'keyC',trim(dict%child%data%key)

    name_is=.false.
    if (trim(name) == trim(dict%data%key)) then
       name_is=.true.
    else if (associated(dict%child)) then
       !most likely dict is a item list
       name_is=(trim(name) == trim(dict%child%data%key))
       !print *,'here',name_is,trim(dict%child%data%value),'ag',trim(name)
    else
       name_is=(trim(name) == trim(dict%data%value))
    end if
  end function name_is

  !> fill output with input and the rest with blanks
  !! this routine is only useful for its interface
  pure subroutine set_field(input,output)
    implicit none
    character(len=*), intent(in) :: input 
    character(len=max_field_length), intent(out) :: output 
    !local variables
    integer :: ipos,i

    !one could also write
    !output(1:len(output))=input

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
    character(len=max_field_length), intent(in) :: input
    character(len=*), intent(out) :: output
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

  !> returns the value of the key of the dictionary
  pure function dict_key(dict)
    type(dictionary), pointer, intent(in) :: dict
    character(len=max_field_length) :: dict_key
    
    if (associated(dict)) then 
       !call check_key(dict)
       dict_key=dict%data%key
    else
       dict_key=repeat(' ',len(dict_key))
    end if
  end function dict_key

  !> returns the value of the key of the dictionary
  pure function dict_item(dict)
    type(dictionary), pointer, intent(in) :: dict
    integer :: dict_item

    if (associated(dict)) then 
       dict_item=dict%data%item
    else
       dict_item=-1
    end if
  end function dict_item

  
  !>value of the dictionary, if present, otherwise empty
  !if the value is a dictionary, it returns __dict__ in the character
  pure function dict_value(dict)
    type(dictionary), pointer, intent(in) :: dict
    character(len=max_field_length) :: dict_value

    if (associated(dict)) then 
       !call check_key(dict)
       if (associated(dict%child)) then
          if (dict%data%nitems > 0) then
             dict_value=TYPE_LIST
          else if (dict%data%nelems > 0) then
             dict_value=TYPE_DICT
          end if
       else
          dict_value=dict%data%value
       end if
    else
       dict_value=repeat(' ',len(dict_value))
    end if

  end function dict_value

  !non-pure subroutines, due to pointer assignments

  !> define the same parent(dict) for any of the elements of the linked chain (child)
  recursive subroutine define_parent(dict,child)
    implicit none
    type(dictionary), target :: dict
    type(dictionary) :: child

    child%parent=>dict
    if (associated(child%next)) call define_parent(dict,child%next)
  end subroutine define_parent

  !> set brother as the previous element of dict
  subroutine define_brother(brother,dict)
    implicit none
    type(dictionary), target :: brother
    type(dictionary) :: dict

    dict%previous=>brother
  end subroutine define_brother

  !> this routine creates a key for the dictionary in case it is absent
  !! the it adds one to the number of elements of the parent dictionary
  subroutine set_elem(dict,key)
    implicit none
    type(dictionary), pointer :: dict !!TO BE VERIFIED
    character(len=*), intent(in) :: key

    !print *,'set_elem in ',trim(key),dict%data%nelems,dict%parent%data%nelems
    call set_field(trim(key),dict%data%key)
    if (associated(dict%parent)) then
       dict%parent%data%nelems=dict%parent%data%nelems+1
    else
       dict%data%nelems=dict%data%nelems+1
    end if
    !print *,'set_elem out ',trim(key),dict%data%nelems,dict%parent%data%nelems

  end subroutine set_elem

  !> associates an extra item to the dictionary and takes care of the 
  !! increment of the number of items
  !! this routine adds the check that the number of items is preserved
  !! for the moment this check is removed
  subroutine set_item(dict,item)
    implicit none
    type(dictionary) :: dict
    integer, intent(in) :: item

    dict%data%item=item
    if (associated(dict%parent)) then
       dict%parent%data%nitems=dict%parent%data%nitems+1
       if (item+1 > dict%parent%data%nitems) then
          !if (exceptions) then
          !   last_error=DICT_ITEM_NOT_VALID
          !else
             write(*,*)'ERROR: item not valid',item,dict%parent%data%nitems
             stop
          !end if
       end if
    end if

  end subroutine set_item

  !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, create it in the child chain
  function get_child_ptr(dict,key) result (subd_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: subd_ptr

    !!commented out, the key is checked only when retrieving
    !call check_key(dict)
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

    !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, create it in the next chain 
  !! Key Must be already present 
  recursive function get_dict_ptr(dict,key) result (dict_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: dict_ptr

!    print *,'here',trim(key)
    !follow the chain, stop at the first occurence
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

    !commented out for the moment
    !call check_key(dict)
    
    if (associated(dict%child)) then
       subd_ptr => get_item_ptr(dict%child,item)
    else
       call dict_init(dict%child)
       call define_parent(dict,dict%child)
       call set_item(dict%child,item)
       subd_ptr => dict%child
    end if

  end function get_list_ptr

  !> defines a storage structure with a key-value couple
  elemental pure function storage_data(key,val)
    character(len=*), intent(in) :: key,val
    type(storage) :: storage_data

    storage_data=storage_null()

    call set_field(key,storage_data%key)
    call set_field(val,storage_data%value)

  end function storage_data

  !test to see the g95 behaviour
  pure function stored_key(st) result(key)
    implicit none
    type(storage), intent(in) :: st
    character(len=max_field_length) :: key
    call get_field(st%key,key)
  end function stored_key

    !test to see the g95 behaviour
  pure function stored_value(st) result(val)
    implicit none
    type(storage), intent(in) :: st
    character(len=max_field_length) :: val
    call get_field(st%value,val)
  end function stored_value


end module dictionaries_base
