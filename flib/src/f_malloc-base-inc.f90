!> @file
!! Include fortran file for f_malloc routines
!! initialize the structure
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  logical, intent(in), optional :: profile
  character(len=*), intent(in), optional :: id,routine_id
  !local variables
  integer :: lgt
  call nullify_malloc_information(m)
  include 'f_malloc-inc.f90'
