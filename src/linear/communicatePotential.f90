subroutine initialize_communication_potential(iproc, nproc, nscatterarr, orbs, lzd, comgp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(p2pComms),intent(out):: comgp
  
  ! Local variables
  integer:: is1, ie1, is2, ie2, is3, ie3, ilr, ii, iorb, iiorb, jproc, kproc, istat, iall, istsource
  integer:: ioverlap, is3j, ie3j, is3k, ie3k, mpidest, istdest, ioffset, is3min, ie3max, tag, p2p_tag, ncount
  integer,dimension(:,:),allocatable:: iStartEnd
  character(len=*),parameter:: subname='setCommunicationPotential'
  
  call nullify_p2pComms(comgp)
  
  ! Determine the bounds of the potential that we need for
  ! the orbitals on this process.
  allocate(iStartEnd(6,0:nproc-1), stat=istat)
  call memocc(istat, iStartEnd, 'iStartEnd', subname)
  is1=0
  ie1=0
  is2=0
  ie2=0
  is3=0
  ie3=0
  iiorb=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc,0)
          
          iiorb=iiorb+1 
          ilr=orbs%inwhichlocreg(iiorb)
      
          ii=lzd%Llr(ilr)%nsi1
          if(ii < is1 .or. iorb==1) then
              is1=ii
          end if
          ii=lzd%Llr(ilr)%nsi1+lzd%Llr(ilr)%d%n1i
          if(ii > ie1 .or. iorb==1) then
              ie1=ii
          end if
      
          ii=lzd%Llr(ilr)%nsi2
          if(ii < is2 .or. iorb==1) then
              is2=ii
          end if
          ii=lzd%Llr(ilr)%nsi2+lzd%Llr(ilr)%d%n2i
          if(ii > ie2 .or. iorb==1) then
              ie2=ii
          end if
      
          ii=lzd%Llr(ilr)%nsi3
          if(ii < is3 .or. iorb==1) then
              is3=ii
          end if
          ii=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i
          if(ii > ie3 .or. iorb==1) then
              ie3=ii
          end if
      
      end do
      iStartEnd(1,jproc)=is1
      iStartEnd(2,jproc)=ie1
      iStartEnd(3,jproc)=is2
      iStartEnd(4,jproc)=ie2
      iStartEnd(5,jproc)=is3
      iStartEnd(6,jproc)=ie3
  end do
  
  ! Determine how many slices each process receives.
  allocate(comgp%noverlaps(0:nproc-1), stat=istat)
  call memocc(istat, comgp%noverlaps, 'comgp%noverlaps', subname)
  do jproc=0,nproc-1
      is3j=istartEnd(5,jproc)
      ie3j=istartEnd(6,jproc)
      mpidest=jproc
      ioverlap=0
      do kproc=0,nproc-1
          is3k=nscatterarr(kproc,3)+1
          ie3k=is3k+nscatterarr(kproc,2)-1
          if(is3j<=ie3k .and. ie3j>=is3k) then
              ioverlap=ioverlap+1
              !if(iproc==0) write(*,'(2(a,i0),a)') 'process ',jproc,' gets potential from process ',kproc,'.' 
          !TAKE INTO ACCOUNT THE PERIODICITY HERE
          else if(ie3j > lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F') then
              ie3j = istartEnd(6,jproc) - lzd%Glr%d%n3i
              if(ie3j>=is3k) then
                 ioverlap=ioverlap+1
              end if
              if(is3j <= ie3k)then
                 ioverlap=ioverlap+1
              end if
          end if
      end do
      comgp%noverlaps(jproc)=ioverlap
      !!if(iproc==0) write(*,'(2(a,i0),a)') 'Process ',jproc,' gets ',ioverlap,' potential slices.'
  end do
  
  ! Determine the parameters for the communications.
  allocate(comgp%overlaps(comgp%noverlaps(iproc)), stat=istat)
  call memocc(istat, comgp%overlaps, 'comgp%overlaps', subname)
  allocate(comgp%comarr(6,maxval(comgp%noverlaps),0:nproc-1))
  call memocc(istat, comgp%comarr, 'comgp%comarr', subname)
  allocate(comgp%ise3(2,0:nproc-1), stat=istat)
  call memocc(istat, comgp%ise3, 'comgp%ise3', subname)
  !allocate(comgp%requests(2,comgp%noverlaps(iproc)), stat=istat)
  allocate(comgp%requests(nproc,2), stat=istat) !nproc is in general too much
  call memocc(istat, comgp%requests, 'comgp%requests', subname)
call to_zero(2*nproc,comgp%requests(1,1))
comgp%nsend = 0 ; comgp%nrecv = 0
  comgp%nrecvBuf = 0
  is3min=0
  ie3max=0
  do jproc=0,nproc-1
      is3j=istartEnd(5,jproc)
      ie3j=istartEnd(6,jproc)
      mpidest=jproc
      ioverlap=0
      istdest=1
      do kproc=0,nproc-1
          is3k=nscatterarr(kproc,3)+1
          ie3k=is3k+nscatterarr(kproc,2)-1
  !SHOULD TAKE INTO ACCOUNT THE PERIODICITY HERE
  !Need to split the region
          if(is3j<=ie3k .and. ie3j>=is3k) then
              is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
              ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
              ioffset=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
              ioverlap=ioverlap+1
              !tag=tag+1
              tag=p2p_tag(jproc)
              if(is3<is3min .or. ioverlap==1) then
                  is3min=is3
              end if
              if(ie3>ie3max .or. ioverlap==1) then
                  ie3max=ie3
              end if
              !!call setCommunicationPotential(kproc, is3, ie3, ioffset, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
              !!     istdest, tag, comgp%comarr(1,ioverlap,jproc))
              istsource=ioffset*lzd%glr%d%n1i*lzd%glr%d%n2i+1
              ncount=(ie3-is3+1)*lzd%glr%d%n1i*lzd%glr%d%n2i
              call setCommsParameters(kproc, jproc, istsource, istdest, ncount, tag, comgp%comarr(1,ioverlap,jproc))
              istdest = istdest + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
              if(iproc==jproc) then
                  comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
              end if
          else if(ie3j > lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F')then
               ie3j = istartEnd(6,jproc) - lzd%Glr%d%n3i
               if(ie3j>=is3k) then
                   is3=max(0,is3k) ! starting index in z dimension for data to be sent
                   ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
                   ioffset=is3-0 ! starting index (in z direction) of data to be sent (actually it is the index -1)
                   ioverlap=ioverlap+1
                   !tag=tag+1
                   tag=p2p_tag(jproc)
                   if(is3<is3min .or. ioverlap==1) then
                       is3min=is3
                   end if
                   if(ie3>ie3max .or. ioverlap==1) then
                       ie3max=ie3
                   end if
                   !!call setCommunicationPotential(kproc, is3, ie3, ioffset, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
                   !!     istdest, tag, comgp%comarr(1,ioverlap,jproc))
                   istsource=ioffset*lzd%glr%d%n1i*lzd%glr%d%n2i+1
                   ncount=(ie3-is3+1)*lzd%glr%d%n1i*lzd%glr%d%n2i
                   call setCommsParameters(kproc, jproc, istsource, istdest, ncount, tag, comgp%comarr(1,ioverlap,jproc))
                   istdest = istdest + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
                   if(iproc==jproc) then
                       comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
                   end if
               end if
               if(is3j <= ie3k)then
                   is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
                   ie3=min(lzd%Glr%d%n3i,ie3k) ! ending index in z dimension for data to be sent
                   ioffset=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
                   ioverlap=ioverlap+1
                   !tag=tag+1
                   tag=p2p_tag(jproc)
                   if(is3<is3min .or. ioverlap==1) then
                       is3min=is3
                   end if
                   if(ie3>ie3max .or. ioverlap==1) then
                       ie3max=ie3
                   end if
                   !!call setCommunicationPotential(kproc, is3, ie3, ioffset, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
                   !!     istdest, tag, comgp%comarr(1,ioverlap,jproc))
                   istsource=ioffset*lzd%glr%d%n1i*lzd%glr%d%n2i+1
                   ncount=(ie3-is3+1)*lzd%glr%d%n1i*lzd%glr%d%n2i
                   call setCommsParameters(kproc, jproc, istsource, istdest, ncount, tag, comgp%comarr(1,ioverlap,jproc))
                   istdest = istdest + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
                   if(iproc==jproc) then
                       comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
                   end if
               end if
          end if
      end do
      comgp%ise3(1,jproc)=is3min
      comgp%ise3(2,jproc)=ie3max
      if(ioverlap/=comgp%noverlaps(jproc)) stop 'ioverlap/=comgp%noverlaps(jproc)'
  end do
  
  
  iall=-product(shape(iStartEnd))*kind(iStartEnd)
  deallocate(iStartEnd, stat=istat)
  call memocc(istat, iall, 'iStartEnd', subname)

  ! To indicate that to communication has been started
  comgp%communication_complete=.true.

end subroutine initialize_communication_potential



!!!!subroutine setCommunicationPotential(mpisource, is3, ie3, ioffset, n1i, n2i, mpidest, istdest, tag, comarr)
!!!!use module_base
!!!!implicit none
!!!!! Calling arguments
!!!!integer,intent(in):: mpisource, is3, ie3, ioffset, n1i, n2i, mpidest, istdest, tag
!!!!integer,dimension(6),intent(out):: comarr
!!!!
!!!!! Local variables
!!!!integer:: istsource, ncount
!!!!
!!!!! From which MPI process shall the slice be sent
!!!!comarr(1)=mpisource
!!!!
!!!!! Starting index on the sending process
!!!!istsource=ioffset*n1i*n2i+1
!!!!comarr(2)=istsource
!!!!
!!!!! Amount of data to be sent
!!!!ncount=(ie3-is3+1)*n1i*n2i
!!!!comarr(3)=ncount
!!!!
!!!!! To which MPI process shall the slice be sent
!!!!comarr(4)=mpidest
!!!!
!!!!! Starting index on the receiving index
!!!!comarr(5)=istdest
!!!!
!!!!! Tag for the communication
!!!!comarr(6)=tag
!!!!
!!!!! comarr(7): this entry is used as request for the mpi_isend.
!!!!
!!!!! comarr(8): this entry is used as request for the mpi_irecv.
!!!!
!!!!
!!!!end subroutine setCommunicationPotential

!!!subroutine gatherPotential(iproc, nproc, comgp)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!!type(p2pCommsGatherPot),intent(inout):: comgp
!!!type(p2pComms),intent(inout):: comgp
!!!
!!!! Local variables
!!!integer:: kproc, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc, ind, i, nsend, nrecv, ncomplete
!!!integer,dimension(mpi_status_size):: stat
!!!logical:: sendComplete, receiveComplete, received
!!!
!!!
!!!
!!!! Wait for the sends to complete.
!!!!!if (nproc > 1) then
!!!   nsend=0
!!!   if(comgp%nsend>0) then
!!!       waitLoopSend: do
!!!          !!call mpi_waitsome(comgp%nsend, comgp%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!!          !!nsend=nsend+ncomplete
!!!          !!if(nsend==comgp%nsend) exit waitLoopSend
!!!          call mpi_waitany(comgp%nsend-nsend, comgp%requests(1,1), ind, mpi_status_ignore, ierr)
!!!          nsend=nsend+1
!!!          do i=ind,comgp%nsend-nsend
!!!             comgp%requests(i,1)=comgp%requests(i+1,1)
!!!          end do
!!!          if(nsend==comgp%nsend) exit waitLoopSend
!!!       end do waitLoopSend
!!!   end if
!!!
!!!
!!!   nrecv=0
!!!   if(comgp%nrecv>0) then
!!!       waitLoopRecv: do
!!!          call mpi_waitany(comgp%nrecv-nrecv, comgp%requests(1,2), ind, mpi_status_ignore, ierr)
!!!          !call mpi_testany(comgp%nrecv-nrecv, comgp%requests(1,2), ind, received, mpi_status_ignore, ierr)
!!!          !ind=1
!!!          ncomplete=1
!!!          received=.true.
!!!          if(received) then
!!!             nrecv=nrecv+ncomplete
!!!             !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comgp%nrecv-nrecv
!!!             !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comgp%requests(1:comgp%nrecv,2)
!!!             do i=ind,comgp%nrecv-nrecv
!!!                comgp%requests(i,2)=comgp%requests(i+1,2)
!!!             end do
!!!             !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comgp%requests(1:comgp%nrecv,2)
!!!             if(nrecv==comgp%nrecv) exit waitLoopRecv
!!!          end if
!!!       end do waitLoopRecv
!!!   end if
!!!!!end if
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!!!if(iproc==0) write(*,'(1x,a)',advance='no') 'Gathering the potential... '
!!!!!! Check whether the communications have completed.
!!!!!nfast=0
!!!!!nsameproc=0
!!!!!testLoop: do 
!!!!!    do jproc=0,nproc-1
!!!!!        do kproc=1,comgp%noverlaps(jproc)
!!!!!           if(comgp%communComplete(kproc,jproc)) cycle
!!!!!            call mpi_test(comgp%comarr(7,kproc,jproc), sendComplete, stat, ierr)      
!!!!!            call mpi_test(comgp%comarr(8,kproc,jproc), receiveComplete, stat, ierr)   
!!!!!            ! Attention: mpi_test is a local function.
!!!!!            if(sendComplete .and. receiveComplete) comgp%communComplete(kproc,jproc)=.true.
!!!!!            !!if(comgp%communComplete(kproc,jproc)) then
!!!!!            !!    !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', korb
!!!!!            !!    mpisource=comgp%comarr(1,kproc,jproc)
!!!!!            !!    mpidest=comgp%comarr(4,kproc,jproc)
!!!!!            !!    if(mpisource/=mpidest) then
!!!!!            !!        nfast=nfast+1
!!!!!            !!    else
!!!!!            !!        nsameproc=nsameproc+1
!!!!!            !!    end if
!!!!!            !!end if
!!!!!        end do
!!!!!    end do
!!!!!    ! If we made it until here, either all all the communication is
!!!!!    ! complete or we better wait for each single orbital.
!!!!!    exit testLoop
!!!!!end do testLoop
!!!!!
!!!!!
!!!!!! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!!!!call mpiallred(comgp%communComplete(1,0), nproc*maxval(comgp%noverlaps), mpi_land, mpi_comm_world, ierr)
!!!!!
!!!!!! Wait for the communications that have not completed yet
!!!!!nslow=0
!!!!!do jproc=0,nproc-1
!!!!!    do kproc=1,comgp%noverlaps(jproc)
!!!!!        if(comgp%communComplete(kproc,jproc)) then
!!!!!            mpisource=comgp%comarr(1,kproc,jproc)
!!!!!            mpidest=comgp%comarr(4,kproc,jproc)
!!!!!            if(mpisource==mpidest) then
!!!!!                nsameproc=nsameproc+1
!!!!!            else
!!!!!                nfast=nfast+1
!!!!!            end if
!!!!!            cycle
!!!!!        end if
!!!!!        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
!!!!!        nslow=nslow+1
!!!!!        call mpi_wait(comgp%comarr(7,kproc,jproc), stat, ierr)  
!!!!!        call mpi_wait(comgp%comarr(8,kproc,jproc), stat, ierr) 
!!!!!        comgp%communComplete(kproc,jproc)=.true.
!!!!!    end do
!!!!!end do
!!!!!
!!!!!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!if (verbose > 3) then
!!!!!   if(iproc==0) write(*,'(a,f5.1,a)') 'done. Communication overlap ratio:',100.d0*dble(nfast)/(dble(nfast+nslow)),'%'
!!!!!else
!!!!!   if(iproc==0) write(*,'(a,f5.1,a)') 'done.'
!!!!!end if
!!!!!!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!!!!                       nfast, ' could be overlapped with computation.'
!!!!!!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!!
!!!
!!!end subroutine gatherPotential



!!subroutine cancelCommunicationPotential(iproc, nproc, comgp)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!!type(p2pCommsGatherPot),intent(inout):: comgp
!!type(p2pComms),intent(inout):: comgp
!!
!!! Local variables
!!integer:: jproc, kproc, ierr
!!integer,dimension(mpi_status_size):: stat
!!logical:: sendComplete, receiveComplete
!!
!!! Cancel all communications. 
!!! It gives errors, therefore simply wait for the communications to complete.
!!do jproc=0,nproc-1
!!    do kproc=1,comgp%noverlaps(jproc)
!!        !call mpi_test(comgp%comarr(7,kproc,jproc), sendComplete, stat, ierr)
!!        !call mpi_test(comgp%comarr(8,kproc,jproc), receiveComplete, stat, ierr)
!!        !if(sendComplete .and. receiveComplete) cycle
!!        !call mpi_cancel(comgp%comarr(7,kproc,jproc), ierr)
!!        !call mpi_cancel(comgp%comarr(8,kproc,jproc), ierr)
!!        call mpi_wait(comgp%comarr(7,kproc,jproc), stat, ierr)
!!        call mpi_wait(comgp%comarr(8,kproc,jproc), stat, ierr)
!!    end do
!!end do
!!
!!end subroutine cancelCommunicationPotential
