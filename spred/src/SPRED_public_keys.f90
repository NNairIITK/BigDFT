!> @file
!! Contains public parameters that have to be used here and there in the code
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module SPRED_public_keys
  implicit none

  public ! guess why?

  !> spred 
  character(len = *), parameter :: SPRED_VARIABLES = "spred"
  character(len = *), parameter :: FP_METHOD_KEY = "fpmethod"
  character(len = *), parameter :: FP_NATX_SPHERE = "natx_sphere"
  character(len = *), parameter :: FP_ANGMOM = "angmom"
end module SPRED_public_keys

!>module identifying constants that have to be used as enumerators
!! they can be used to define f_enumerator types or directly as integers
module SPRED_public_enums
  use f_enums
  implicit none
  
  private !as we use at least one module, private become compulsory

  !> spred
  type(f_enumerator), parameter, public :: OMF_FP_METHOD      =f_enumerator('OMF_FP_METHOD',-1000,null())
  type(f_enumerator), parameter, public :: OMP_FP_METHOD      =f_enumerator('OMP_FP_METHOD',-999,null())
end module SPRED_public_enums





