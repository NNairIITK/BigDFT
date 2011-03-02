! simple MPI fake routines to compile a serial version of some MPI code without using an MPI compiler

subroutine  MPI_INIT(ierr)
integer::ierr
ierr=0
return
end subroutine

subroutine MPI_FINALIZE(ierr)
integer::ierr
ierr=0
return
end subroutine


subroutine  MPI_COMM_RANK(isMPI_COMM_WORLD,iproc,ierr)
integer::iproc,ierr
ierr=0
iproc=0
return
end subroutine


subroutine  MPI_COMM_SIZE(isMPI_COMM_WORLD,nproc,ierr)
integer::nproc,ierr
ierr=0
nproc=1
return
end subroutine

subroutine  MPI_BARRIER(isMPI_COMM_WORLD,ierr) 
integer::ierr
ierr=0
return
end subroutine


subroutine MPI_BCAST(sendbuf,n,isMPI_DOUBLE_PRECISION,iroot,&
                     isMPI_COMM_WORLD,ierr) 
integer::ierr
real(8)::arr(n)
ierr=0
return
end subroutine


subroutine MPI_ALLREDUCE(sendbuf,recbuf,n,isMPI_DOUBLE_PRECISION,&
                         isMPI_SUM,isMPI_COMM_WORLD,ierr)
REAL(8),dimension(n):: sendbuf,recbuf
integer::ierr
recbuf=sendbuf
ierr=0
return
end subroutine


subroutine MPI_ALLTOALL(sendbuf,nsend,isMPI_DOUBLE_PRECISION,&
                        recbuf,nrec,isalsoMPI_DOUBLE_PRECISION,&
                        isMPI_COMM_WORLD,ierr)
REAL(8),dimension(nsend):: sendbuf
REAL(8),dimension(nrec):: recbuf
integer::ierr
ierr=0
return
end subroutine
