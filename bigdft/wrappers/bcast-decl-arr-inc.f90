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
  logical, intent(in), optional :: check !< performs the check of the arguments
  !local variables
  logical chk
  integer :: n,iroot,mpi_comm,ierr
  external :: MPI_BCAST

  chk=.false.
  n=size(buffer)

