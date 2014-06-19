!> @file
!! Include fortran file for deallocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!variable declaration, to be included in dynamic_memory.f90
  !local variables
  integer :: ierror
!!$  character(len=info_length) :: address
  logical :: use_global
  integer(kind=8) :: ilsize,jlsize,iadd
  character(len=namelen) :: array_id,routine_id
  type(dictionary), pointer :: dict_add

  if (f_err_raise(ictrl == 0,&
       'ERROR (f_free): the routine f_malloc_initialize has not been called',&
       ERR_MALLOC_INTERNAL)) return

  !here we should add a control of the OMP behaviour of allocation
  !in particular for what concerns the OMP nesting procedure

