!!!!!!subroutine postCommunicationsPotential(iproc, nproc, ndimpot, pot, comgp)
!!!!!!use module_base
!!!!!!use module_types
!!!!!!implicit none
!!!!!!
!!!!!!! Calling arguments
!!!!!!integer,intent(in):: iproc, nproc, ndimpot
!!!!!!real(8),dimension(ndimpot),intent(in):: pot
!!!!!!!type(p2pCommsGatherPot),intent(inout):: comgp
!!!!!!type(p2pComms),intent(inout):: comgp
!!!!!!
!!!!!!! Local variables
!!!!!!integer:: jproc, kproc, nsends, nreceives, istat, mpisource, istsource, ncount, mpidest, istdest, tag, ierr
!!!!!!
!!!!!!
!!!!!!! Post the messages
!!!!!!if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for communicating the potential... '
!!!!!!
!!!!!!
!!!!!!! First only post receives
!!!!!!nreceives=0
!!!!!!do jproc=0,nproc-1
!!!!!!    do kproc=1,comgp%noverlaps(jproc)
!!!!!!        mpisource=comgp%comarr(1,kproc,jproc)
!!!!!!        istsource=comgp%comarr(2,kproc,jproc)
!!!!!!        ncount=comgp%comarr(3,kproc,jproc)
!!!!!!        mpidest=comgp%comarr(4,kproc,jproc)
!!!!!!        istdest=comgp%comarr(5,kproc,jproc)
!!!!!!        tag=comgp%comarr(6,kproc,jproc)
!!!!!!        if(iproc==mpidest) then
!!!!!!            nreceives=nreceives+1
!!!!!!            call mpi_irecv(comgp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
!!!!!!                 comgp%requests(nreceives,2), ierr)
!!!!!!        end if
!!!!!!    end do
!!!!!!end do
!!!!!!
!!!!!!! Now the sends
!!!!!!nsends=0
!!!!!!do jproc=0,nproc-1
!!!!!!    do kproc=1,comgp%noverlaps(jproc)
!!!!!!        mpisource=comgp%comarr(1,kproc,jproc)
!!!!!!        istsource=comgp%comarr(2,kproc,jproc)
!!!!!!        ncount=comgp%comarr(3,kproc,jproc)
!!!!!!        mpidest=comgp%comarr(4,kproc,jproc)
!!!!!!        istdest=comgp%comarr(5,kproc,jproc)
!!!!!!        tag=comgp%comarr(6,kproc,jproc)
!!!!!!        if(iproc==mpisource) then
!!!!!!            nsends=nsends+1
!!!!!!            call mpi_isend(pot(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
!!!!!!                 comgp%requests(nsends,1), ierr)
!!!!!!        end if
!!!!!!    end do
!!!!!!end do
!!!!!!
!!!!!!
!!!!!!!!!nreceives=0
!!!!!!!!!nsends=0
!!!!!!!!!comgp%communComplete=.false.
!!!!!!!!!destLoop: do jproc=0,nproc-1
!!!!!!!!!    sourceLoop: do kproc=1,comgp%noverlaps(jproc)
!!!!!!!!!        mpisource=comgp%comarr(1,kproc,jproc)
!!!!!!!!!        istsource=comgp%comarr(2,kproc,jproc)
!!!!!!!!!        ncount=comgp%comarr(3,kproc,jproc)
!!!!!!!!!        mpidest=comgp%comarr(4,kproc,jproc)
!!!!!!!!!        istdest=comgp%comarr(5,kproc,jproc)
!!!!!!!!!        tag=comgp%comarr(6,kproc,jproc)
!!!!!!!!!        if(mpisource/=mpidest) then
!!!!!!!!!            if(iproc==mpisource) then
!!!!!!!!!                !!write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', &
!!!!!!!!!                !!    istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
!!!!!!!!!                !!call mpi_isend(pot(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
!!!!!!!!!                !!     comgp%comarr(7,kproc,jproc), ierr)
!!!!!!!!!                nsends=nsends+1
!!!!!!!!!                call mpi_isend(pot(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
!!!!!!!!!                     comgp%requests(nsends,1), ierr)
!!!!!!!!!                comgp%comarr(8,kproc,jproc)=mpi_request_null !is this correct?
!!!!!!!!!            else if(iproc==mpidest) then
!!!!!!!!!               !!write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, &
!!!!!!!!!               !!    ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!!!!!!!                !!call mpi_irecv(comgp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
!!!!!!!!!                !!     comgp%comarr(8,kproc,jproc), ierr)
!!!!!!!!!                nreceives=nreceives+1
!!!!!!!!!                call mpi_irecv(comgp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
!!!!!!!!!                     comgp%requests(nreceives,2), ierr)
!!!!!!!!!                comgp%comarr(7,kproc,jproc)=mpi_request_null !is this correct?
!!!!!!!!!            else
!!!!!!!!!                comgp%comarr(7,kproc,jproc)=mpi_request_null
!!!!!!!!!                comgp%comarr(8,kproc,jproc)=mpi_request_null
!!!!!!!!!            end if
!!!!!!!!!        else
!!!!!!!!!            if(iproc==mpisource) then
!!!!!!!!!               !!write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, &
!!!!!!!!!               !!    ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!!!!!!!                nsends=nsends+1
!!!!!!!!!                call mpi_isend(pot(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
!!!!!!!!!                     comgp%requests(nsends,1), ierr)
!!!!!!!!!                nreceives=nreceives+1
!!!!!!!!!                    call mpi_irecv(comgp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
!!!!!!!!!                         comgp%requests(nreceives,2), ierr)
!!!!!!!!!            end if
!!!!!!!!!        end if
!!!!!!!!!    end do sourceLoop
!!!!!!!!!end do destLoop
!!!!!!!!!if(iproc==0) write(*,'(a)') 'done.'
!!!!!!
!!!!!!comgp%nsend=nsends
!!!!!!comgp%nrecv=nreceives
!!!!!!
!!!!!!if(nreceives/=comgp%noverlaps(iproc)) then
!!!!!!    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comgp%noverlaps(iproc)',&
!!!!!!         nreceives, comgp%noverlaps(iproc)
!!!!!!    stop
!!!!!!end if
!!!!!!
!!!!!!
!!!!!!end subroutine postCommunicationsPotential


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
