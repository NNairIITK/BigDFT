!> @file
!!  Basic file defining the structures to deal with the sparse matrices
!! @author
!!    Copyright (C) 2014-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining the basic structures and routines related to the sparse matrix library
module sparsematrix_base
  ! f_lib modules
  use dynamic_memory
  use dictionaries
  use yaml_output
  use yaml_strings
  use module_defs, only: UNINITIALIZED
  use wrapper_MPI
  use wrapper_linalg
  use f_utils
  use numerics
  use time_profiling

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
