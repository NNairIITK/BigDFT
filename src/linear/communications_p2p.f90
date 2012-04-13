subroutine post_p2p_communication(iproc, nproc, nsendbuf, sendbuf, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, nsendbuf
  real(8),dimension(nsendbuf),intent(in):: sendbuf
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  integer:: jproc, joverlap, nsends, nreceives, mpisource, istsource, ncount, mpidest, istdest, tag, ierr
  
  
  !!! Post the messages
  !!if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for p2p communication'
  
  
  ! First only post receives
  nreceives=0
  do jproc=0,nproc-1
      do joverlap=1,comm%noverlaps(jproc)
          mpisource=comm%comarr(1,joverlap,jproc)
          istsource=comm%comarr(2,joverlap,jproc)
          ncount=comm%comarr(3,joverlap,jproc)
          mpidest=comm%comarr(4,joverlap,jproc)
          istdest=comm%comarr(5,joverlap,jproc)
          tag=comm%comarr(6,joverlap,jproc)
          if(iproc==mpidest) then
              nreceives=nreceives+1
              call mpi_irecv(comm%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
                   comm%requests(nreceives,2), ierr)
          end if
      end do
  end do
  
  ! Now the sends
  nsends=0
  do jproc=0,nproc-1
      do joverlap=1,comm%noverlaps(jproc)
          mpisource=comm%comarr(1,joverlap,jproc)
          istsource=comm%comarr(2,joverlap,jproc)
          ncount=comm%comarr(3,joverlap,jproc)
          mpidest=comm%comarr(4,joverlap,jproc)
          istdest=comm%comarr(5,joverlap,jproc)
          tag=comm%comarr(6,joverlap,jproc)
          if(iproc==mpisource) then
              nsends=nsends+1
              call mpi_isend(sendbuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
                   comm%requests(nsends,1), ierr)
          end if
      end do
  end do
  
  !!if(iproc==0) write(*,'(a)') 'done.'
  
  comm%nsend=nsends
  comm%nrecv=nreceives
  
  if(nreceives/=comm%noverlaps(iproc)) then
      write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comm%noverlaps(iproc)',&
           nreceives, comm%noverlaps(iproc)
    stop
  end if


end subroutine post_p2p_communication


subroutine wait_p2p_communication(iproc, nproc, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  integer:: ierr, ind, i, nsend, nrecv
  
  
  ! Wait for the sends to complete.
  nsend=0
  if(comm%nsend>0) then
      wait_sends: do
         call mpi_waitany(comm%nsend-nsend, comm%requests(1,1), ind, mpi_status_ignore, ierr)
         nsend=nsend+1
         do i=ind,comm%nsend-nsend
            comm%requests(i,1)=comm%requests(i+1,1)
         end do
         if(nsend==comm%nsend) exit wait_sends
      end do wait_sends
  end if
 
 
  ! Wait for the receives to complete.
  nrecv=0
  if(comm%nrecv>0) then
      wait_recvs: do
         call mpi_waitany(comm%nrecv-nrecv, comm%requests(1,2), ind, mpi_status_ignore, ierr)
         nrecv=nrecv+1
         do i=ind,comm%nrecv-nrecv
            comm%requests(i,2)=comm%requests(i+1,2)
         end do
         if(nrecv==comm%nrecv) exit wait_recvs
      end do wait_recvs
  end if

end subroutine wait_p2p_communication
