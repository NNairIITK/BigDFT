subroutine post_p2p_communication(iproc, nproc, nsendbuf, sendbuf, nrecvbuf, recvbuf, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
  real(8),dimension(nsendbuf),intent(in):: sendbuf
  real(8),dimension(nrecvbuf),intent(out):: recvbuf
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  integer:: jproc, joverlap, nsends, nreceives, mpisource, istsource, ncount, mpidest, istdest, tag, ierr

  if(.not.comm%communication_complete) stop 'ERROR: there is already a p2p communication going on...'
  
  nreceives=0
  nsends=0
  ! First only post receives
  do jproc=0,nproc-1
      do joverlap=1,comm%noverlaps(jproc)
          mpisource=comm%comarr(1,joverlap,jproc)
          istsource=comm%comarr(2,joverlap,jproc)
          ncount=comm%comarr(3,joverlap,jproc)
          mpidest=comm%comarr(4,joverlap,jproc)
          istdest=comm%comarr(5,joverlap,jproc)
          tag=comm%comarr(6,joverlap,jproc)
          if(ncount>0) then
              if(nproc>1) then
                  if(iproc==mpidest) then
                      if(mpidest/=mpisource) then
                          nreceives=nreceives+1
                          call mpi_irecv(recvbuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
                               comm%requests(nreceives,2), ierr)
                      else
                          call dcopy(ncount, sendbuf(istsource), 1, recvbuf(istdest), 1)
                      end if
                  end if
              else
                  nsends=nsends+1
                  nreceives=nreceives+1
                  call dcopy(ncount, sendbuf(istsource), 1, recvbuf(istdest), 1)
              end if
          end if
      end do
  end do
  
  ! Now the sends
  do jproc=0,nproc-1
      do joverlap=1,comm%noverlaps(jproc)
          mpisource=comm%comarr(1,joverlap,jproc)
          istsource=comm%comarr(2,joverlap,jproc)
          ncount=comm%comarr(3,joverlap,jproc)
          mpidest=comm%comarr(4,joverlap,jproc)
          istdest=comm%comarr(5,joverlap,jproc)
          tag=comm%comarr(6,joverlap,jproc)
          if(ncount>0) then
              if(nproc>1) then
                  if(iproc==mpisource) then
                      if(mpisource/=mpidest) then
                          nsends=nsends+1
                          call mpi_isend(sendbuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
                           comm%requests(nsends,1), ierr)
                      end if
                  end if
              end if
          end if
      end do
  end do
  
  !!if(iproc==0) write(*,'(a)') 'done.'
  
  comm%nsend=nsends
  comm%nrecv=nreceives
  
  !!if(nreceives/=comm%noverlaps(iproc)) then
  !!    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comm%noverlaps(iproc)',&
  !!         nreceives, comm%noverlaps(iproc)
  !!  stop
  !!end if
  
  ! Flag indicating whether the communication is complete or not
  if(nproc>1) then
      comm%communication_complete=.false.
      comm%messages_posted=.true.
  else
      comm%communication_complete=.true.
      comm%messages_posted=.false.
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
  
  
  if(.not.comm%communication_complete) then

      if(.not.comm%messages_posted) stop 'ERROR: trying to wait for messages which have never been posted!'

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

  end if

  ! Flag indicating that the communication is complete
  comm%communication_complete=.true.

end subroutine wait_p2p_communication



module p2p_tags_data
  use module_base
  implicit none
  logical,save:: initialized
  integer,dimension(:),allocatable,save:: tags
  integer,save:: tag_max
end module p2p_tags_data


subroutine init_p2p_tags(nproc)
  use module_base
  use p2p_tags_data
  implicit none

  ! Calling arguments
  integer,intent(in):: nproc
  character(len=*),parameter:: subname='init_p2p_tags'

  ! Local variables
  integer:: jproc, istat, ierr
  logical:: success

  if(initialized) stop 'trying to initialize the counter for the p2p tags which is already running!'

  allocate(tags(0:nproc-1),stat=istat)
  call memocc(istat,tags,'tags',subname)
  do jproc=0,nproc-1
      tags(jproc)=0
  end do
  initialized=.true.

  ! Determine the largest possible tag
  if (nproc > 1 ) then
     call mpi_attr_get(mpi_comm_world, mpi_tag_ub, tag_max, success, ierr)
     if(.not.success) stop 'could not extract largest possible tag...'
  else
     tag_max=1
  end if
end subroutine init_p2p_tags


function p2p_tag(jproc)
  use module_base
  use p2p_tags_data
  implicit none

  ! Calling arguments
  integer,intent(in):: jproc
  integer:: p2p_tag

  if(.not.initialized) stop 'counter for tag was not properly initialized!'

  tags(jproc)=mod(tags(jproc)+1,tag_max)
  p2p_tag=tags(jproc)

end function p2p_tag


subroutine finalize_p2p_tags()
  use module_base
  use p2p_tags_data
  implicit none

  ! Local variables
  integer:: istat, iall
  character(len=*),parameter:: subname='finalize_p2p_tags'

  if(.not.initialized) stop 'trying to finalize the counter for the p2p tags which was not initialized!'

  iall=-product(shape(tags))*kind(tags)
  deallocate(tags,stat=istat)
  call memocc(istat,iall,'tags',subname)

  initialized=.false.

end subroutine finalize_p2p_tags
