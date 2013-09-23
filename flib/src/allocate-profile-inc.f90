!> @file
!! Include fortran file for allocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  integer :: ierror
  integer(kind=8) :: iadd

  if (f_err_raise(ictrl == 0,&
       'ERROR (f_malloc): the routine f_malloc_initialize has not been called',&
       ERR_MALLOC_INTERNAL)) return
