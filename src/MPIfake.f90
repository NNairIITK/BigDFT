!!****f* BigDFT/MPIfake
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

        subroutine  MPI_comm_group()
          implicit none
          stop 'MPIFAKE: comm_group'
        END SUBROUTINE  MPI_comm_group

        subroutine  MPI_comm_create()
          implicit none
          stop 'MPIFAKE: comm_create'
        END SUBROUTINE  MPI_comm_create

        subroutine  MPI_group_incl()
          implicit none
          stop 'MPIFAKE: comm_incl'
        END SUBROUTINE  MPI_group_incl

        subroutine  MPI_recv()
          implicit none
          stop 'MPIFAKE: recv'
        END SUBROUTINE  MPI_recv

        subroutine  MPI_send()
          implicit none
          stop 'MPIFAKE: send'
        END SUBROUTINE  MPI_send


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

        subroutine  MPI_SCATTERV()
        implicit none
        stop 'MPIFAKE: SCATTERV'
        END SUBROUTINE  MPI_SCATTERV

        subroutine  MPI_ABORT()
        implicit none
        stop 'MPIFAKE: MPI_ABORT'
        END SUBROUTINE  MPI_ABORT

        subroutine  MPI_IRECV()
        implicit none
        stop 'MPIFAKE: IRECV'
        END SUBROUTINE  MPI_IRECV
        
        subroutine  MPI_ISEND()
        implicit none
        stop 'MPIFAKE: ISEND'
        END SUBROUTINE  MPI_ISEND
        
        subroutine  MPI_WAITALL()
        implicit none
        stop 'MPIFAKE: WAITALL'
        END SUBROUTINE  MPI_WAITALL

        subroutine MPI_INITIALIZED(init,ierr)
          implicit none
          integer, intent(out) :: init,ierr
          init=1
          ierr=0
        END SUBROUTINE  MPI_INITIALIZED

        subroutine MPI_GET_PROCESSOR_NAME()
        implicit none
        stop 'MPIFAKE: MPI_GET_PROCESSOR_NAME'
        end subroutine  MPI_GET_PROCESSOR_NAME

        subroutine  mpi_error_string()
          implicit none
          stop 'MPIFAKE: mpi_error_string'
        END SUBROUTINE  MPI_ERROR_STRING

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
