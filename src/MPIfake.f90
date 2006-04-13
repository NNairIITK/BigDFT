        subroutine  MPI_INIT(ierr)
        return
        end

        subroutine  MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        return
        end

        subroutine  MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        return
        end

        subroutine  MPI_FINALIZE(ierr)
        return
        end

        subroutine  MPI_ALLREDUCE(wrkallred,allred,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        return
        end

        subroutine  MPI_BCAST(psi,nvctr,MPI_DOUBLE_PRECISION,jproc,MPI_COMM_WORLD,ierr)
        return
        end

        subroutine  MPI_BARRIER(MPI_COMM_WORLD,ierr)
        return
        end

        subroutine MPI_REDUCE (t1,t2,norbe,MPI_DOUBLE_PRECISION,MPI_SUM,j,MPI_COMM_WORLD,ierr)
        return
        end

        subroutine MPI_ALLGATHER(zf,md,MPI_double_precision,rhopot,md2,&
                 MPI_double_precision2,MPI_COMM_WORLD,ierr)
        return
        end

        subroutine MPI_ALLGATHERV(psi,sendcount,MPI_DOUBLE_PRECISION,  &
                          tpsi,recvcounts,displs, MPI_DOUBLE_PRECISION2,MPI_COMM_WORLD,ierr)
        return
        end
        
        subroutine MPI_ALLTOALL(zmpi2,n1, MPI_double_precision, &
                          zmpi1,n2, MPI_double_precision2,MPI_COMM_WORLD,ierr)
        return
        end
