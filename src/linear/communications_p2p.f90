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
  
  
  ! Post the messages
  if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for communicating the sendbufential... '
  
  
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
  
  if(iproc==0) write(*,'(a)') 'done.'
  
  comm%nsend=nsends
  comm%nrecv=nreceives
  
  if(nreceives/=comm%noverlaps(iproc)) then
      write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comm%noverlaps(iproc)',&
           nreceives, comm%noverlaps(iproc)
    stop
  end if


end subroutine post_p2p_communication
