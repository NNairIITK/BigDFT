!> @file
!! Define operations on enumerators
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_enums
  implicit none

  private

  integer, parameter, private :: NULL_INT=-1024
  character(len=*), parameter, private :: null_name='nullified enumerator'

  !> enumerator type, useful to define different modes
  type, public :: f_enumerator
     character(len=64) :: name=null_name
     integer :: id=null_int
     type(f_enumerator), pointer :: family => null()
  end type f_enumerator

  type(f_enumerator), parameter, private :: &
       f_enum_null=f_enumerator(null_name,NULL_INT,null())

  interface operator(==)
     module procedure enum_is_enum,enum_is_char,enum_is_int
  end interface operator(==)

  interface operator(/=)
     module procedure enum_is_not_enum,enum_is_not_char,enum_is_not_int
  end interface operator(/=)

  interface operator(.hasattr.)
     module procedure enum_has_attribute,enum_has_char,enum_has_int
  end interface operator(.hasattr.)

  interface int
     module procedure int_enum
  end interface int

  interface str
     module procedure char_enum
  end interface str

  public :: f_enum_attr,operator(.hasattr.),nullify_f_enum
  public :: int,str,f_enumerator_null,operator(==),operator(/=)

contains

  pure subroutine nullify_f_enum(en)
    implicit none
    type(f_enumerator), intent(out) :: en
    en=f_enumerator_null()
  end subroutine nullify_f_enum

  pure function f_enumerator_null() result(en)
    use yaml_strings, only: f_strcpy
    implicit none
    type(f_enumerator) :: en
    call f_strcpy(src=null_name,dest=en%name)
    en%id=null_int
    nullify(en%family)
  end function f_enumerator_null

  !>associate the enumerator to a family
  subroutine f_enum_attr(dest,attr)
    implicit none
    type(f_enumerator), intent(inout), target :: attr !< to avoid the parameter attribute
    type(f_enumerator), intent(inout) :: dest
    !local variables
    type(f_enumerator), pointer :: iter
    !print *,'ASSOCIATING',trim(char(dest)),trim(char(attr))
    !search the first iterator that has no family
    if (.not. associated(dest%family)) then
       !print *,'here'
       dest%family=>attr
    else
       !print *,'there'
       !if the enumerator already exists do not associate
       iter => dest%family
       do while(associated(iter%family) .and. (iter /= attr))
          iter => iter%family
       end do
       if (iter /= attr) iter%family=>attr
    end if
  end subroutine f_enum_attr

  function enum_has_attribute(en,family) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: family
    type(f_enumerator), intent(in) :: en
    logical :: ok
    !local variables
    type(f_enumerator), pointer :: iter

    ok=.false.
    iter=>en%family
    do while(associated(iter) .and. .not. ok)
       ok= iter == family
       iter => iter%family
    end do
  end function enum_has_attribute

  function enum_has_char(en,family) result(ok)
    implicit none
    character(len=*), intent(in) :: family
    type(f_enumerator), intent(in) :: en
    logical :: ok
    !local variables
    type(f_enumerator), pointer :: iter

    ok=.false.
    iter=>en%family
    do while(associated(iter) .and. .not. ok)
       ok= iter == family
       iter => iter%family
    end do
  end function enum_has_char

  function enum_has_int(en,family) result(ok)
    implicit none
    integer, intent(in) :: family
    type(f_enumerator), intent(in) :: en
    logical :: ok
    !local variables
    type(f_enumerator), pointer :: iter

    ok=.false.
    iter=>en%family
    do while(associated(iter) .and. .not. ok)
       ok= iter == family
       iter => iter%family
    end do
  end function enum_has_int

  elemental pure function enum_is_enum(en,en1) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    type(f_enumerator), intent(in) :: en1
    logical :: ok
    ok = (en == en1%id) .and. (en == en1%name)
  end function enum_is_enum

  elemental pure function enum_is_int(en,int) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    integer, intent(in) :: int
    logical :: ok
    ok = en%id == int
  end function enum_is_int

  elemental pure function enum_is_char(en,char) result(ok)
    use yaml_strings, only: operator(.eqv.)
    implicit none
    type(f_enumerator), intent(in) :: en
    character(len=*), intent(in) :: char
    logical :: ok
    ok = trim(en%name) .eqv. trim(char)
  end function enum_is_char

  elemental pure function enum_is_not_enum(en,en1) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    type(f_enumerator), intent(in) :: en1
    logical :: ok
    ok = .not. (en == en1)
  end function enum_is_not_enum

  elemental pure function enum_is_not_int(en,int) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    integer, intent(in) :: int
    logical :: ok
    !ok = .not. (en == int)
    ok = .not. enum_is_int(en,int)
  end function enum_is_not_int

  elemental pure function enum_is_not_char(en,char) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    character(len=*), intent(in) :: char
    logical :: ok
    ok = .not. (en == char)
  end function enum_is_not_char


  !>integer of f_enumerator type.
  elemental pure function int_enum(en)
    type(f_enumerator), intent(in) :: en
    integer :: int_enum
    int_enum=en%id
  end function int_enum

  !>char of f_enumerator type.
  elemental pure function char_enum(en)
    type(f_enumerator), intent(in) :: en
    character(len=len(en%name)) :: char_enum
    char_enum=en%name
  end function char_enum


end module f_enums
