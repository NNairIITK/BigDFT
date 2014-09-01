!> @file
!! Include fortran file for broadcast operations
!! declarations for arrays
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in), optional :: root  !< @copydoc doc::root
  integer, intent(in), optional :: comm  !< @copydoc doc::comm  
  !local variables
  integer :: n,iroot,mpi_comm,ierr
  external :: MPI_BCAST

  n=size(buffer)
  !no bcast if size is zero
  if (n==0) return

