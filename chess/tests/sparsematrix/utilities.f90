!> @file
!!   File containing helper routines for the test drivers
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


module utilities
  private

  public :: get_ccs_data_from_file
  public :: operation_using_dense_lapack
  public :: check_deviation_from_unity_dense

  contains
    subroutine get_ccs_data_from_file(filename, nfvctr, nvctr, row_ind, col_ptr)
      use sparsematrix_init, only: read_ccs_format
      use dynamic_memory, only: f_free_ptr
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
