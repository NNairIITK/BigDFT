!> @file
!!   Basic file which allows to use the sparse matrices
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


!> Module defining the basic structures and routines related to the sparse matrix library
module sparsematrix_base
  !use dynamic_memory
  !use dictionaries
  !use yaml_output
  !use yaml_strings
  !use wrapper_MPI
  !use wrapper_linalg
  !use f_utils
  !use numerics
  !use time_profiling

  ! Very basic sparsematrix modules
  use sparsematrix_errorhandling
  use sparsematrix_timing
  use sparsematrix_types
 use sparsematrix_memory

  implicit none

  ! This module is public, such that all other modules using this one inherit all modules used in here

  ! Old and new version of the sparse matrix matrix multiplication
  integer,parameter :: MATMUL_NEW = 101
  integer,parameter :: MATMUL_OLD = 102
  integer,parameter :: matmul_version = MATMUL_NEW 


end module sparsematrix_base
