!> @file
!! Include fortran file for maxdiff operations with scalars
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (present(comm)) then
     mpi_comm=comm
  else
     mpi_comm=MPI_COMM_WORLD 
  end if
  if (present(root)) then
     iroot=root
  else
     iroot=0
  end if
  nproc=mpisize(mpi_comm)

  if (nproc == 1 .or. ndims == 0) return

  !check that the positions are identical for all the processes
  array_glob=f_malloc((/ndims,nproc/),id='array_glob')

  call mpigather(sendbuf=array,sendcount=ndims,recvbuf=array_glob,&
       root=iroot,comm=mpi_comm)

  if ( mpirank(mpi_comm) == iroot) then
     do jproc=2,nproc
        do i=1,ndims
           maxdiff=max(maxdiff,&
                abs(array_glob(i,jproc)-array_glob(i,1)))
        end do
     end do
  end if

  call f_free(array_glob)
