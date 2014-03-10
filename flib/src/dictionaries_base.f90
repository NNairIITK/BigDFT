!> @file
!!  Module defining a dictionary
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
  character(len=max_field_length), parameter :: TYPE_DICT='__dict__'
  character(len=max_field_length), parameter :: TYPE_LIST='__list__'

  character(len = max_field_length), parameter, private :: NOT_A_VALUE = "__not_a_value__"

  !global variables associated to the number of dictionaries allocated
  integer, private :: ndicts=0         !< number of dictionaries allocated simultaneously
  integer, private :: ndicts_max=0     !< maximum number of dictionaries allocated
  integer, private :: nfolders=0       !< number of libraries
  integer, private :: nfolders_max=0   !< maximum number of libraries allocated
  
  integer, private :: nfolder_size=10000 !< size of the folder of pre-allocated dictionaries. Should be about 7 MB.
  
  
  type, public :: storage
     sequence
     integer :: item   !< Id of the item associated to the list
     integer :: nitems !< No. of items in the list
     integer :: nelems !< No. of items in the dictionary
     character(len=max_field_length) :: key
     character(len=max_field_length) :: value
  end type storage

  !> structure of the dictionary element (this internal structure is in principle private)
  type, public :: dictionary
     type(storage) :: data
     type(dictionary), pointer :: parent => null()
     type(dictionary), pointer :: next => null()
     type(dictionary), pointer :: child => null()
     type(dictionary), pointer :: previous => null()
  end type dictionary

  !> the database book is a workspace of pre-allocated dictionaries, which is used to manage dictionary creation
  type, private :: database_book
     !> this is the registry of the dictionary. It keeps track of the address of the associated dictionaries
     !it is null when the corresponding file item is free
     integer(kind=8), dimension(:), pointer :: registry => null()
     !>this is the place occupied by a dictionary, which is pre-allocated by
     !! the allocation in database_book. The routine dict_init associates the dictionary to a database item
     type(dictionary), dimension(:), pointer :: files => null()
     !> this is useful to create new books of the library in case the book is full
     type(database_book), pointer :: next => null()
     type(database_book), pointer :: previous => null()
  end type database_book


  !>database of objects, associated contiguously and used for storage
  type(database_book), pointer, private :: library => null()

  !> operator to access and create a key in the dictionary
  interface operator(//)
     module procedure get_child_ptr,get_list_ptr
  end interface

  private :: allocate_library, allocate_file, deallocate_file, destroy_library

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

    no_value=trim(dict%data%value) == NOT_A_VALUE .and. .not. associated(dict%child)
  end function no_value

  pure function storage_null() result(st)
    type(storage) :: st
    st%key=repeat(' ',max_field_length)
    st%value(1:max_field_length)=NOT_A_VALUE
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

  pure subroutine allocate_library(library,previous)
    implicit none
    type(database_book), pointer :: previous 
    type(database_book), intent(out) :: library 
    !local variables
    integer :: ifolder

    allocate(library%registry(nfolder_size))
    allocate(library%files(nfolder_size))
    nullify(library%next)
    library%previous=>previous
    !all free places in the library
    do ifolder=1,nfolder_size
       library%registry(ifolder)=int(0,kind=8)
    end do
  end subroutine allocate_library

  !> assign a place in the library to the newcoming dictionary
  recursive function allocate_file(previous) result(dict)
    implicit none
    type(dictionary), pointer :: dict
    type(database_book), pointer, optional :: previous !> previous library to be linked with
    !local variables
    integer :: ifolder
    integer(kind=8), external :: f_loc
    type(database_book), pointer :: lib_tmp
        
    !first create the place for a folder if there is none
    if (.not. associated(library%registry)) then
       if (present(previous)) then
          lib_tmp=>previous
       else
          nullify(lib_tmp)
       end if
       call allocate_library(library,lib_tmp)
       nfolders=nfolders+1
       nfolders_max=max(nfolders_max,nfolders)
    end if
    !then search in the present library the first free place
    ifolder=0
    find_free_place: do
       ifolder=ifolder+1
       if (ifolder > nfolder_size) exit find_free_place
       if (library%registry(ifolder) == int(0,kind=8)) exit find_free_place
    end do find_free_place
    !if no free place has been found move to the next library
    if (ifolder > nfolder_size) then 
       !then the library which is used is the following one
       lib_tmp=>library
       if (.not. associated(library%next)) allocate(library%next) !the pointer
       library=>library%next
       dict=>allocate_file(previous=lib_tmp)
    else 
       !if the free place has been found the dictionary can be identified
       !the registry is updated
       library%registry(ifolder)=f_loc(library%files(ifolder))
       call dictionary_nullify(library%files(ifolder))
       dict=>library%files(ifolder)
    end if
  end function allocate_file

  recursive subroutine deallocate_file(file)
    implicit none
    integer(kind=8), intent(in) :: file !< the target is already inside

    !first, find last library
    if (associated(library%next)) then
       library=>library%next
       call deallocate_file(file)
    else
       call deallocate_file_(file)
       !then search the file to destroy
    end if

  contains
    !free the space from the library
    recursive subroutine deallocate_file_(file)
      implicit none
      integer(kind=8), intent(in) :: file !< the target is already inside
      !local variables
      integer :: ifolder
      !find the file in the registry
      !then search in the present library the first free place
      ifolder=0
      find_place: do
         ifolder=ifolder+1
         if (ifolder > nfolder_size) exit find_place
         if (library%registry(ifolder) == file) exit find_place
      end do find_place
      !if no place has been found go back to the previous library
      if (ifolder > nfolder_size) then 
         if (associated(library%previous)) then
            library=>library%previous
            call deallocate_file_(file)
         else
          !this might also mean that a variable has not been found
          stop 'dictionary has not been allocated by pool'
         end if
      else 
         !if the place has been found the dictionary can be identified
         !the registry is freed and the dictionary nullified
         library%registry(ifolder)=int(0,kind=8)
         call dictionary_nullify(library%files(ifolder))
      end if

    end subroutine deallocate_file_
  end subroutine deallocate_file

  !terminate the library and free all memory space
  recursive subroutine destroy_library()
    implicit none
    !first, find last library
    if (associated(library%next)) then
       library=>library%next
       call destroy_library()
    else
       !once we are at the last tree free memory spaces
       deallocate(library%registry)
       nullify(library%registry)
       deallocate(library%files)
       nullify(library%files)
    end if
    !then come back of one step
    if (associated(library%previous)) then
       library => library%previous
       !and free the memory space of the library
       deallocate(library%next)
       nullify(library%next)
       nfolders=nfolders-1
       call destroy_library()
    end if
  end subroutine destroy_library

  subroutine dict_init(dict)
    implicit none
    type(dictionary), pointer :: dict
    if (nfolder_size == 0) then
       !this is normal allocation, let the system do the allocation
       allocate(dict)
       call dictionary_nullify(dict)
    else
       !this strategy retrieves the allocation from the library
       if (.not. associated(library)) allocate(library) ! the pointer
       dict=>allocate_file()
    end if
    !global variables counting the number of dictionaries simultaneously used
    ndicts=ndicts+1
    ndicts_max=max(ndicts_max,ndicts)
  end subroutine dict_init

  !> destroy only one level
  subroutine dict_destroy(dict)
    implicit none
    type(dictionary), pointer :: dict
    !local variables
    integer(kind=8), external :: f_loc
    !strategy which uses systems' allocation
    if (nfolder_size==0) then
       deallocate(dict)
    else
       !free a space in the library and let the allocation lives
       call deallocate_file(f_loc(dict))
    end if

    ndicts=ndicts-1
    nullify(dict)
    !if there are no dictionaries anymore the library can be destroyed
    if (ndicts==0 .and. associated(library)) then
       call destroy_library()
       !otherwise free the last tree 
       deallocate(library)
       nullify(library)
    end if
  end subroutine dict_destroy

  subroutine dict_free(dict)
    type(dictionary), pointer :: dict

    if (associated(dict)) then
       call dict_free_(dict)
       call dict_destroy(dict)
       nullify(dict)
    end if

  contains

    recursive subroutine dict_free_(dict)
      implicit none
      type(dictionary), pointer :: dict

      !first destroy the children
      if (associated(dict%child)) then
         call dict_free_(dict%child)
         call dict_destroy(dict%child)
      end if
      !then destroy younger brothers
      if (associated(dict%next)) then
         call dict_free_(dict%next)
         call dict_destroy(dict%next)
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
          else
             dict_value= NOT_A_VALUE !illegal condition
          end if
       else if (trim(dict%data%value) == NOT_A_VALUE) then
          dict_value=repeat(' ',len(dict_value))
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
  pure subroutine set_elem(dict,key)
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
  pure subroutine set_item(dict,item)
    implicit none
    type(dictionary), pointer :: dict !< recently modified, to be checked
    integer, intent(in) :: item

    dict%data%item=item
    if (associated(dict%parent)) then
       if (len_trim(dict%parent%data%value) > 0) dict%parent%data%value=repeat(' ',max_field_length)
       dict%parent%data%nitems=dict%parent%data%nitems+1
    end if

  end subroutine set_item

  !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, create it in the child chain
  function get_child_ptr(dict,key) result(subd_ptr)
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
  subroutine item_ptr_find(dict,item,item_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    integer, intent(in) :: item
    type(dictionary), pointer :: item_ptr

    item_ptr=>dict
    find_item: do 
       if (item_ptr%data%item == item) exit find_item
       if (associated(item_ptr%next)) then
          item_ptr=>item_ptr%next
          cycle find_item
       end if
       if (no_key(item_ptr)) then
          call set_item(item_ptr,item)
       else
          call dict_init(item_ptr%next)
          call define_brother(item_ptr,item_ptr%next) !chain the list in both directions
          if (associated(item_ptr%parent)) &
               call define_parent(item_ptr%parent,item_ptr%next)
          call set_item(item_ptr%next,item)
          item_ptr=>item_ptr%next          
       end if
       exit find_item
    end do find_item
  end subroutine item_ptr_find


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
    
    !print *,'nelems',dict%data%nelems,dict%data%nitems
    !free previously existing dictionaries
    call clean_subdict(dict)
    
    if (associated(dict%child)) then         
       !subd_ptr => get_item_ptr(dict%child,item)
       call item_ptr_find(dict%child,item,subd_ptr)
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

  !> clean the child and put to zero the number of elements in the case of a dictionary
  !pure 
  subroutine clean_subdict(dict)
    implicit none
    type(dictionary), pointer :: dict
    
    if (associated(dict)) then
       if (dict%data%nelems > 0) then
          call dict_free(dict%child)
          dict%data%nelems=0
       end if
    end if
    
  end subroutine clean_subdict

  !!> export the number of dictionaries
  !! this routine is useful to understand the total usage of the 
  !! dictionaries in the f_lib module
  subroutine dict_get_num(ndict,ndict_max,nlibs,nlibs_max)
    implicit none
    integer, intent(out) :: ndict     !< actual number of dictionaries active
    integer, intent(out) :: ndict_max !< maximum number of dictionaries active during the session
    integer, intent(out) :: nlibs     !< actual number of libraries active at present
    integer, intent(out) :: nlibs_max !< maximum number of libraries active at present (should be nlibs)

    ndict=ndicts
    ndict_max=ndicts_max
    nlibs=nfolders
    nlibs_max=nfolders_max
  end subroutine dict_get_num

end module dictionaries_base

!> Routines for bindings only (external of module)
subroutine dict_new(dict)
  use dictionaries_base
  implicit none
  type(dictionary), pointer :: dict

  call dict_init(dict)
end subroutine dict_new

subroutine dict_free(dict)
  use dictionaries_base, mod_dict_free => dict_free
  implicit none
  type(dictionary), pointer :: dict

  call mod_dict_free(dict)
end subroutine dict_free
