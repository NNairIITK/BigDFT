!> @file
!!  Module base of BigDFT package
!! @author Luigi Genovese
!!    Copyright (C) 2008-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Modules which contains the low level definitions, as well as some profiling procedures
module module_base 
  use wrapper_linalg
  use wrapper_MPI
  use module_defs
  use dictionaries, dict_set => set !error_handling
  use dynamic_memory
  use time_profiling
  use f_utils
  use yaml_strings
  implicit none  
end module module_base
