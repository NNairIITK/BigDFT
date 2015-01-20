!> @file
!!  Define main module for BigDFT API
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG, DC)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Module which contains all information about types (data structures) and interfaces
module BigDFT_API
  use module_interfaces
  use module_base
  use module_types
  use module_xc
  use module_atoms
  use module_input_dicts
  use psp_projectors
  use ao_inguess
  use communications_base
  use communications_init
  use gaussians
  use public_keys
  use public_enums
  use bigdft_run
  implicit none
end module BigDFT_API
