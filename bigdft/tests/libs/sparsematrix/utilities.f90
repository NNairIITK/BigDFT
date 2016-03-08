!> @file
!! Test of the sparsematrix library
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


module utilities
  private

  public :: get_ccs_data_from_file

  contains
    subroutine get_ccs_data_from_file(filename, nfvctr, nvctr, row_ind, col_ptr)
      use module_base
      use sparsematrix_init, only: read_ccs_format
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nfvctr, nvctr
      integer,dimension(:),pointer,intent(out) :: row_ind
      integer,dimension(:),pointer,intent(out) :: col_ptr

      ! Local variables
      real(kind=8),dimension(:),pointer :: val

      call read_ccs_format(filename, nfvctr, nvctr, col_ptr, row_ind, val)

      call f_free_ptr(val)

    end subroutine get_ccs_data_from_file
end module utilities
