!!****f* PSEUDO/MPIfake
!! FUNCTION
!!    Fake functions for MPI in the case of serial version
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine  MPI_INIT(ierr)
  implicit none
  integer, intent(out) :: ierr
  ierr=0
END SUBROUTINE MPI_INIT
!!***
        
subroutine MPI_INITIALIZED(init,ierr)
  implicit none
  integer, intent(out) :: init,ierr
  init=1
  ierr=0
END SUBROUTINE  MPI_INITIALIZED

subroutine  MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: iproc,ierr
  iproc=0
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_RANK

subroutine  MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: nproc,ierr
  nproc=1
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_SIZE

subroutine  MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: MPI_GROUP,ierr
  MPI_GROUP=1
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_GROUP

subroutine  MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP,MPI_COMM,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: MPI_GROUP,MPI_COMM,ierr
  MPI_GROUP=1
  MPI_COMM=1
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_CREATE

subroutine  MPI_GROUP_INCL(GROUP,N,NRANKS,NEWGROUP,ierr)
  implicit none
  integer, intent(in) :: GROUP,N
  integer, intent(in) :: NRANKS(N)
  integer, intent(out) :: NEWGROUP,ierr
  NEWGROUP=1
  ierr=GROUP*0
END SUBROUTINE MPI_GROUP_INCL

!here we have routines which do not transform the argument for nproc==1
!these routines can be safely called also in the serial version
subroutine  MPI_FINALIZE(ierr)
  implicit none
  integer, intent(out) :: ierr
  ierr=0
END SUBROUTINE MPI_FINALIZE

subroutine MPI_BCAST()
  implicit none
END SUBROUTINE MPI_BCAST

subroutine  MPI_BARRIER(MPI_COMM_WORLD,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: ierr
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_BARRIER

subroutine MPI_REDUCE()
  implicit none
END SUBROUTINE MPI_REDUCE

subroutine  MPI_ALLREDUCE()
  implicit none
END SUBROUTINE MPI_ALLREDUCE


! These routines in serial version should not be called.
! A stop is added

subroutine  MPI_ALLGatherV()
  implicit none
  stop 'MPIFAKE: ALLGATHERV'
END SUBROUTINE  MPI_ALLGatherV

subroutine  MPI_ALLGATHER()
  implicit none
  stop 'MPIFAKE: ALLGATHER'
END SUBROUTINE  MPI_ALLGATHER

subroutine  MPI_GatherV()
  implicit none
  stop 'MPIFAKE: GATHERV'
END SUBROUTINE  MPI_GatherV

subroutine  MPI_Gather()
  implicit none
  stop 'MPIFAKE: GATHER'
END SUBROUTINE  MPI_Gather


subroutine  MPI_ALLTOALL()
  implicit none
  stop 'MPIFAKE: ALLTOALL'
END SUBROUTINE  MPI_ALLTOALL

subroutine  MPI_ALLTOALLV()
  implicit none
  stop 'MPIFAKE: ALLTOALLV'
END SUBROUTINE  MPI_ALLTOALLV

subroutine  MPI_REDUCE_SCATTER()
  implicit none
  stop 'MPIFAKE: REDUCE_SCATTER'
END SUBROUTINE  MPI_REDUCE_SCATTER

subroutine  MPI_ABORT()
  implicit none
  stop 'MPIFAKE: MPI_ABORT'
END SUBROUTINE  MPI_ABORT

subroutine  MPI_IRECV()
  implicit none
  stop 'MPIFAKE: IRECV'
END SUBROUTINE  MPI_IRECV

subroutine  MPI_RECV()
  implicit none
  stop 'MPIFAKE: RECV'
END SUBROUTINE  MPI_RECV

subroutine  MPI_ISEND()
  implicit none
  stop 'MPIFAKE: ISEND'
END SUBROUTINE  MPI_ISEND

subroutine  MPI_SEND()
  implicit none
  stop 'MPIFAKE: SEND'
END SUBROUTINE  MPI_SEND

subroutine  MPI_WAITALL()
  implicit none
  stop 'MPIFAKE: WAITALL'
END SUBROUTINE  MPI_WAITALL

subroutine MPI_GET_PROCESSOR_NAME()
  implicit none
  stop 'MPIFAKE: MPI_GET_PROCESSOR_NAME'
end subroutine  MPI_GET_PROCESSOR_NAME

subroutine  mpi_error_string()
  implicit none
  stop 'MPIFAKE: mpi_error_string'
END SUBROUTINE  MPI_ERROR_STRING

subroutine  MPI_SCATTERV()
  implicit none
  stop 'MPIFAKE: SCATTERV'
END SUBROUTINE  MPI_SCATTERV

subroutine mpi_attr_get ()
  implicit none
  stop 'MPIFAKE: mpi_attr_get'
END SUBROUTINE  MPI_ATTR_GET

subroutine mpi_type_size ()
  implicit none
  stop 'MPIFAKE: mpi_type_size'
END SUBROUTINE  MPI_TYPE_SIZE

subroutine mpi_comm_free ()
  implicit none
  stop 'MPIFAKE: mpi_comm_free'
END SUBROUTINE  MPI_COMM_FREE
