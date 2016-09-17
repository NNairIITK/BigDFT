!> @file
!!   File handling the sparse matrix errors
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


module sparsematrix_errorhandling
  implicit none

  private

  !> Errorcodes
  integer,save,public :: SPARSEMATRIX_ALLOCATION_ERROR
  integer,save,public :: SPARSEMATRIX_MANIPULATION_ERROR
  integer,save,public :: SPARSEMATRIX_RUNTIME_ERROR
  integer,save,public :: SPARSEMATRIX_INITIALIZATION_ERROR

  !> Public routine
  public :: sparsematrix_init_errors

  contains

    !> Define the sparsematrix errors
    subroutine sparsematrix_init_errors()
      use dictionaries
      implicit none

      call f_err_define('SPARSEMATRIX_ALLOCATION_ERROR',&
           'a problem occured during the allocation of a sparse matrix',&
           SPARSEMATRIX_ALLOCATION_ERROR,&
           err_action='Check the calling arguments of the allocation routine')

      call f_err_define('SPARSEMATRIX_MANIPULATION_ERROR',&
           'a problem occured during the manipulation ((un)compression,sparsity pattern transformation) of a sparse matrix',&
           SPARSEMATRIX_MANIPULATION_ERROR,&
           err_action='Check the calling arguments of the manipulation routine and the array sizes')

      call f_err_define('SPARSEMATRIX_RUNTIME_ERROR',&
           'a general problem related to sparse matrices occured during runtime',&
           SPARSEMATRIX_MANIPULATION_ERROR,&
           err_action='Check the dedicated error message')

      call f_err_define('SPARSEMATRIX_INITIALIZATION_ERROR',&
           'a problem related to the initialization of a sparse matrix occured',&
           SPARSEMATRIX_MANIPULATION_ERROR,&
           err_action='Check the calling arguments and the dedicated error message')
  
      end subroutine sparsematrix_init_errors

end module sparsematrix_errorhandling
