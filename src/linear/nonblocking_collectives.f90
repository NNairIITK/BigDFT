subroutine post_mpi_ialltoallv_dble(iproc, nproc, nsendbuf, sendbuf, sendcounts, senddspls, &
           nrecvbuf, recvbuf, recvcounts, recvdspls, comm, requests, communication_complete, messages_posted)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsendbuf, nrecvbuf, comm
  integer,dimension(0:nproc-1),intent(in) :: sendcounts, senddspls, recvcounts, recvdspls
  real(kind=8),dimension(nsendbuf),intent(in) :: sendbuf
  real(kind=8),dimension(nrecvbuf),intent(out) :: recvbuf
  integer,dimension(0:nproc-1,2),intent(out) :: requests
  logical,intent(inout) :: communication_complete
  logical,intent(out) :: messages_posted

  ! Local variables
  integer :: jproc, ierr, ist

  if(.not.communication_complete) stop 'ERROR: there is already a p2p communication going on...'

  do jproc=0,nproc-1
      ist=senddspls(jproc)+1
      call mpi_isend(sendbuf(ist), sendcounts(jproc), mpi_double_precision, &
           jproc, jproc, comm, requests(jproc,1), ierr)
  end do

  do jproc=0,nproc-1
      ist=recvdspls(jproc)+1
      call mpi_irecv(recvbuf(ist), recvcounts(jproc), mpi_double_precision, &
           jproc, iproc, comm, requests(jproc,2), ierr)
  end do

  messages_posted=.true.
  communication_complete=.false.

end subroutine post_mpi_ialltoallv_dble


subroutine wait_mpi_ialltoallv_dble(nproc, messages_posted, communication_complete, requests)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: nproc
  logical,intent(in) :: messages_posted
  logical,intent(inout) :: communication_complete
  integer,dimension(0:nproc-1,2),intent(inout) :: requests

  ! Local variables
  integer :: ierr

  if (.not.communication_complete) then

      if (.not.messages_posted) stop 'ERROR: trying to test for messages which have never been posted!'

      ! Tests the sends.
      call mpi_waitall(nproc, requests(0,1), mpi_statuses_ignore, ierr)
 
      ! Tests the sends.
      call mpi_waitall(nproc, requests(0,2), mpi_statuses_ignore, ierr)

  end if

  communication_complete=.true.

end subroutine wait_mpi_ialltoallv_dble
