!> @file
!!  template to be used to add different items to a list
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!! included in dictionaries.f90
!header of routine
!!$  subroutine add(dict,val)
!!$    implicit none
!!$    type(dictionary), pointer :: dict
!!$    <generic>, intent(in) :: val
  !local variables
  integer :: length,isize

  ! Consistency checks.
  isize=dict_size(dict)
  length=dict_len(dict)
  if (f_err_raise(isize > 0,'Add not allowed for this node',&
       err_id=DICT_INVALID_LIST)) return
  if (f_err_raise(length == -1,'Add not allowed for this node',&
       err_id=DICT_INVALID)) return

  if (present(last_item_ptr)) then
     ! Use sibling as last item to add to.
     if (associated(last_item_ptr)) then
        call init_next(last_item_ptr)
        call set_item(last_item_ptr, length)
     else
        last_item_ptr => dict // length
     end if
     call set(last_item_ptr, val)
  else
     ! Add new item at the end.
     call set(dict//length,val)
  end if

