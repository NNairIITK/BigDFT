program mpi_check
  use BigDFT_API
  implicit none
  logical :: failed
  integer :: ntot,ierr,i,iproc,nproc,nspin,nother
  integer, dimension(:,:), allocatable  :: isizes
  real(kind=8), dimension(:), allocatable :: copybuffer,buffer

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ntot=10000
  nspin=2
  nother=10

  allocate(buffer(ntot),copybuffer(ntot),isizes(nspin,nother))

  buffer=1.0d0
  copybuffer=0.d0
  isizes=1
  !last component indicates the size
  isizes(nspin,nother)=ntot

  !reference (supposed)
  call MPI_ALLREDUCE(buffer,copybuffer,ntot,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  call mpiallred(buffer(1),isizes(nspin,nother),MPI_SUM,MPI_COMM_WORLD,ierr)
  failed=.false.
  do i=1,ntot
     if (abs(copybuffer(i)-buffer(i))>1.d-12 .or. abs(copybuffer(i)-real(nproc,kind=8))> 1.d-12) then
        write(*,'(a,i7,2(1pe25.17),i6,i6)')'Failed',i,buffer(i),copybuffer(i),iproc,nproc
        failed=.true.
     end if
  end do
  
  deallocate(buffer,copybuffer,isizes)

  !wait all processes before finalisation
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  print *,'Ok, iproc',iproc   
  call MPI_FINALIZE(ierr)
end program mpi_check
