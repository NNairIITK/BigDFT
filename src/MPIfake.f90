        subroutine  MPI_INIT(ierr)
        stop 'INIT'
        return
        end

        subroutine  MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        stop 'RANK'
        return
        end

        subroutine  MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        stop 'SIZE'
        return
        end

        subroutine  MPI_FINALIZE(ierr)
        stop 'FINAL'
        return
        end

        subroutine  MPI_ALLREDUCE(wrkallred,allred,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        stop 'ALLREDUCE'
        return
        end

        subroutine  MPI_BCAST(psi,nvctr,MPI_DOUBLE_PRECISION,jproc,MPI_COMM_WORLD,ierr)
        stop 'BCAST'
        return
        end

        subroutine  MPI_BARRIER(MPI_COMM_WORLD,ierr)
        stop 'BARRIER'
        return
        end

        subroutine MPI_REDUCE (t1,t2,norbe,MPI_DOUBLE_PRECISION,MPI_SUM,j,MPI_COMM_WORLD,ierr)
        stop 'REDUCE'
        return
        end

        subroutine  MPI_ALLGatherV()
        stop 'ALLGATHERV'
        return
        end

        subroutine  MPI_ALLGather()
        stop 'ALLGATHER'
        return
        end

        subroutine  MPI_ALLTOALL()
        stop 'ALLTOALL'
        return
        end

