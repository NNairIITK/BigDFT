!!****f* BigDFT/MPIfake
!! FUNCTION
!!    Fake functions for MPI in the case onf serial version
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2008 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
        subroutine  MPI_INIT(ierr)
        implicit none
        integer :: ierr
        end subroutine MPI_INIT
!!***
        
        subroutine  MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        implicit none
        integer :: MPI_COMM_WORLD,iproc,ierr
        iproc=0
        ierr=0
        end subroutine MPI_COMM_RANK

        subroutine  MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        implicit none
        integer :: MPI_COMM_WORLD,nproc,ierr
        nproc=1
        ierr=0
        end subroutine MPI_COMM_SIZE

        subroutine  MPI_FINALIZE(ierr)
        implicit none
        integer :: ierr
        end subroutine MPI_FINALIZE

        subroutine MPI_BCAST()
        implicit none
        end subroutine MPI_BCAST

        subroutine  MPI_BARRIER(MPI_COMM_WORLD,ierr)
        implicit none
        integer MPI_COMM_WORLD,ierr
        end subroutine MPI_BARRIER

        subroutine MPI_REDUCE()
        implicit none
        end subroutine MPI_REDUCE

! These routines in serial version should not be called.
! A stop is added

        subroutine  MPI_ALLREDUCE()
        implicit none
        stop 'MPIFAKE: ALLREDUCE'
        end subroutine MPI_ALLREDUCE

        subroutine  MPI_ALLGatherV()
        implicit none
        stop 'MPIFAKE: ALLGATHERV'
        end subroutine  MPI_ALLGatherV

        subroutine  MPI_ALLGather()
        implicit none
        stop 'MPIFAKE: ALLGATHER'
        end subroutine  MPI_ALLGather

        subroutine  MPI_GatherV()
        implicit none
        stop 'MPIFAKE: GATHERV'
        end subroutine  MPI_GatherV

        subroutine  MPI_Gather()
        implicit none
        stop 'MPIFAKE: GATHER'
        end subroutine  MPI_Gather

        subroutine  MPI_ALLTOALL()
        implicit none
        stop 'MPIFAKE: ALLTOALL'
        end subroutine  MPI_ALLTOALL

        subroutine  MPI_ALLTOALLV()
        implicit none
        stop 'MPIFAKE: ALLTOALLV'
        end subroutine  MPI_ALLTOALLV

        subroutine  MPI_REDUCE_SCATTER()
        implicit none
        stop 'MPIFAKE: REDUCE_SCATTER'
        end subroutine  MPI_REDUCE_SCATTER

        subroutine  MPI_ABORT()
        implicit none
        stop 'MPIFAKE: MPI_ABORT'
        end subroutine  MPI_ABORT

        subroutine  MPI_IRECV()
        implicit none
        stop 'MPIFAKE: IRECV'
        end subroutine  MPI_IRECV
        
        subroutine  MPI_ISEND()
        implicit none
        stop 'MPIFAKE: ISEND'
        end subroutine  MPI_ISEND
        
        subroutine  MPI_WAITALL()
        implicit none
        stop 'MPIFAKE: WAITALL'
        end subroutine  MPI_WAITALL
