!> @file
!! Include fortran file for broadcast operations
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  iroot=0
  if (present(root)) iroot=root
  mpi_comm=MPI_COMM_WORLD
  if (present(comm)) mpi_comm=comm

  call MPI_BCAST(buffer,n,mpitype(buffer),iroot,mpi_comm,ierr)
  if (ierr /=0) then
     call f_err_throw('An error in calling to MPI_BCAST occured',&
          err_id=ERR_MPI_WRAPPERS)
     return
  end if
