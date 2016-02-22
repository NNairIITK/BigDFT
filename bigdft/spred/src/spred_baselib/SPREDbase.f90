!> @file
!!    Modulefile for handling basic definitions for SPRED
!!
!! @author
!!    B. Schaefer, L. Genovese
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module SPREDbase
  use f_precisions
!!  use wrapper_linalg
!!  use wrapper_MPI
  use numerics
  use dictionaries, dict_set => set !error_handling
  use dynamic_memory
!!  use time_profiling
  use f_utils
  use f_enums, f_char=>str, f_int => toi
  use f_refcnts
  use f_trees
  use yaml_strings
  implicit none
  ! General precision, density and the potential types, to be moved in a low-levle module
  integer, parameter, public :: gp=f_double  !< general-type precision
  integer, parameter, public :: dp=f_double  !< density-type precision
end module SPREDbase
