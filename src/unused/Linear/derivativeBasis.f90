!!!subroutine postCommsRepartition(iproc, nproc, orbs, comrp, nsendBuf, sendBuf, nrecvBuf, recvBuf)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
!!!type(orbitals_data),intent(in):: orbs
!!!!type(p2pCommsRepartition),intent(inout):: comrp
!!!type(p2pComms),intent(inout):: comrp
!!!real(8),dimension(nsendBuf),intent(in):: sendBuf
!!!real(8),dimension(nrecvBuf),intent(out):: recvBuf
!!!
!!!! Local variables
!!!integer:: nsends, nreceives, jproc, jorb, mpisource, mpidest, istsource, istdest, ncount, tag, ierr, isend, irecv
!!!
!!!
!!!
!!!
!!!! First only post receives
!!!irecv=0
!!!do jproc=0,nproc-1
!!!    !do jorb=1,comrp%noverlaps(jproc)
!!!    do jorb=1,4*orbs%norb_par(jproc,0)
!!!        mpisource=comrp%comarr(1,jorb,jproc)
!!!        istsource=comrp%comarr(2,jorb,jproc)
!!!        ncount=comrp%comarr(3,jorb,jproc)
!!!        mpidest=comrp%comarr(4,jorb,jproc)
!!!        istdest=comrp%comarr(5,jorb,jproc)
!!!        tag=comrp%comarr(6,jorb,jproc)
!!!        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
!!!        !irecv=irecv+1
!!!        if(iproc==mpidest .and. nproc >1) then
!!!            irecv=irecv+1
!!!            !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!            !call mpi_irecv(comrp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!!!            !     mpi_comm_world, comrp%comarr(8,iorb,jproc), ierr)
!!!            tag=jorb
!!!            !t1=mpi_wtime()
!!!            !write(*,'(2(a,i0))') 'iproc=',iproc,' receives data at ',istdest
!!!            !write(*,'(a,2(i0,a))') 'RECV: process ',iproc,' receives data from process ',mpisource,'.'
!!!            call mpi_irecv(recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, &
!!!                 mpi_comm_world, comrp%requests(irecv,2), ierr)
!!!            !t2=mpi_wtime()
!!!            !timecommun=timecommun+t2-t1
!!!        !else
!!!        !     comrp%requests(irecv,2)=mpi_request_null
!!!        end if
!!!    end do
!!!end do
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call timing(iproc,'p2pOrtho_post ','ON')
!!!write(*,*) 'old posts: irecv',irecv
!!!
!!!
!!!! Now the sends
!!!isend=0
!!!do jproc=0,nproc-1
!!!    !do jorb=1,comrp%noverlaps(jproc)
!!!    do jorb=1,4*orbs%norb_par(jproc,0)
!!!        mpisource=comrp%comarr(1,jorb,jproc)
!!!        istsource=comrp%comarr(2,jorb,jproc)
!!!        ncount=comrp%comarr(3,jorb,jproc)
!!!        mpidest=comrp%comarr(4,jorb,jproc)
!!!        istdest=comrp%comarr(5,jorb,jproc)
!!!        tag=comrp%comarr(6,jorb,jproc)
!!!        ! The orbitals are on different processes, so we need a point to point communication.
!!!        !isend=isend+1
!!!        if(iproc==mpisource .and. nproc >1) then
!!!            isend=isend+1
!!!            !t2=mpi_wtime()
!!!            !timeextract=timeextract+t2-t1
!!!            tag=jorb
!!!            !t1=mpi_wtime()
!!!            call mpi_irsend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, &
!!!                 mpi_comm_world, comrp%requests(isend,1), ierr)
!!!            !call mpi_rsend(comrp%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!!!            !     mpi_comm_world, ierr)
!!!            !t2=mpi_wtime()
!!!            !timecommun=timecommun+t2-t1
!!!         else if (nproc==1) then
!!!            call vcopy(ncount,sendBuf(istsource),1,recvBuf(istdest),1)
!!!        end if
!!!    end do
!!!end do
!!!
!!!
!!!
!!!
!!!
!!!!!!!nsends=0
!!!!!!!nreceives=0
!!!!!!!comrp%communComplete=.false.
!!!!!!!do jproc=0,nproc-1
!!!!!!!    do jorb=1,4*orbs%norb_par(jproc,0)
!!!!!!!        mpisource=comrp%comarr(1,jorb,jproc)
!!!!!!!        istsource=comrp%comarr(2,jorb,jproc)
!!!!!!!        ncount=comrp%comarr(3,jorb,jproc)
!!!!!!!        mpidest=comrp%comarr(4,jorb,jproc)
!!!!!!!        istdest=comrp%comarr(5,jorb,jproc)
!!!!!!!        tag=comrp%comarr(6,jorb,jproc)
!!!!!!!        !write(*,'(10(a,i0))') 'post on iproc=',iproc,': process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,'; tag=',tag,'  jproc=',jproc,' jorb=',jorb,' ncnt=',comrp%comarr(3,jorb,jproc)
!!!!!!!        if(mpisource/=mpidest) then
!!!!!!!            ! The orbitals are on different processes, so we need a point to point communication.
!!!!!!!            if(iproc==mpisource) then
!!!!!!!                !write(*,'(9(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag,'  jproc=',jproc,' jorb=',jorb,' ncnt=',comrp%comarr(3,jorb,jproc)
!!!!!!!                call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, &
!!!!!!!                     comrp%comarr(7,jorb,jproc), ierr)
!!!!!!!                !!do ierr=istsource,istsource+ncount-1
!!!!!!!                !!    write(tag,*) ierr, sendBuf(ierr)
!!!!!!!                !!end do
!!!!!!!                !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
!!!!!!!                comrp%comarr(8,jorb,jproc)=mpi_request_null !is this correct?
!!!!!!!                nsends=nsends+1
!!!!!!!            else if(iproc==mpidest) then
!!!!!!!                !write(*,'(9(a,i0))') 'process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag,'  jproc=',jproc,' jorb=',jorb,' ncnt=',comrp%comarr(3,jorb,jproc)
!!!!!!!                call mpi_irecv(recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world, &
!!!!!!!                     comrp%comarr(8,jorb,jproc), ierr)
!!!!!!!                comrp%comarr(7,jorb,jproc)=mpi_request_null !is this correct?
!!!!!!!                nreceives=nreceives+1
!!!!!!!            else
!!!!!!!                comrp%comarr(7,jorb,jproc)=mpi_request_null
!!!!!!!                comrp%comarr(8,jorb,jproc)=mpi_request_null
!!!!!!!            end if
!!!!!!!        else
!!!!!!!            ! The orbitals are on the same process, so simply copy them.
!!!!!!!            if(iproc==mpisource) then
!!!!!!!                call dcopy(ncount, sendBuf(istsource), 1, recvBuf(istdest), 1)
!!!!!!!                !!do ierr=istsource,istsource+ncount-1
!!!!!!!                !!    write(tag,*) ierr, sendBuf(ierr)
!!!!!!!                !!end do
!!!!!!!                !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
!!!!!!!                comrp%comarr(7,jorb,jproc)=mpi_request_null
!!!!!!!                comrp%comarr(8,jorb,jproc)=mpi_request_null
!!!!!!!                nsends=nsends+1
!!!!!!!                nreceives=nreceives+1
!!!!!!!                comrp%communComplete(jorb,iproc)=.true.
!!!!!!!            else
!!!!!!!                comrp%comarr(7,jorb,jproc)=mpi_request_null
!!!!!!!                comrp%comarr(8,jorb,jproc)=mpi_request_null
!!!!!!!                comrp%communComplete(jorb,iproc)=.true.
!!!!!!!            end if
!!!!!!!
!!!!!!!        end if
!!!!!!!    end do
!!!!!!!end do
!!!
!!!
!!!end subroutine postCommsRepartition




!!!subroutine gatherDerivativeOrbitals(iproc, nproc, orbs, comrp)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(orbitals_data),intent(in):: orbs
!!!!type(p2pCommsRepartition),intent(inout):: comrp
!!!type(p2pComms),intent(inout):: comrp
!!!
!!!! Local variables
!!!integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc, nsend, nrecv, ncomplete, ind, i
!!!integer,dimension(mpi_status_size):: stat
!!!logical:: sendComplete, receiveComplete, received
!!!
!!!
!!!
!!!
!!!
!!!!!allocate(indexarray(comrp%nrecv), stat=istat)
!!!!!call memocc(istat, indexarray, 'indexarray', subname)
!!!!!do i=1,comrp%nrecv
!!!!!   indexarray(i)=i
!!!!!end do
!!!
!!!!call timing(iproc,'p2pOrtho_wait ','ON')
!!!
!!!! Wait for the sends to complete.
!!!if (nproc > 1) then
!!!   !t1=mpi_wtime()
!!!   nsend=0
!!!   if(comrp%nsend>0) then
!!!       waitLoopSend: do
!!!          !!call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!!          !!nsend=nsend+ncomplete
!!!          !!if(nsend==comon%nsend) exit waitLoopSend
!!!          call mpi_waitany(comrp%nsend-nsend, comrp%requests(1,1), ind, mpi_status_ignore, ierr)
!!!          nsend=nsend+1
!!!          do i=ind,comrp%nsend-nsend
!!!             comrp%requests(i,1)=comrp%requests(i+1,1)
!!!          end do
!!!          if(nsend==comrp%nsend) exit waitLoopSend
!!!       end do waitLoopSend
!!!   end if
!!!   !t2=mpi_wtime()
!!!   !timecommunp2p=timecommunp2p+t2-t1
!!!
!!!
!!!   nrecv=0
!!!   if(comrp%nrecv>0) then
!!!       waitLoopRecv: do
!!!          !t1=mpi_wtime()
!!!          call mpi_waitany(comrp%nrecv-nrecv, comrp%requests(1,2), ind, mpi_status_ignore, ierr)
!!!          !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
!!!          !ind=1
!!!          !t2=mpi_wtime()
!!!          !timecommunp2p=timecommunp2p+t2-t1
!!!          ncomplete=1
!!!          received=.true.
!!!          if(received) then
!!!             nrecv=nrecv+ncomplete
!!!             !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
!!!             !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
!!!             do i=ind,comrp%nrecv-nrecv
!!!                comrp%requests(i,2)=comrp%requests(i+1,2)
!!!                !indexarray(i)=indexarray(i+1)
!!!             end do
!!!             !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
!!!             if(nrecv==comrp%nrecv) exit waitLoopRecv
!!!          end if
!!!       end do waitLoopRecv
!!!   end if
!!!end if
!!!!write(*,'(a,i0,a)') 'iproc=',iproc,' is here'
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!
!!!!!!!
!!!!!!!
!!!!!!!! Check whether the communications have completed.
!!!!!!!nfast=0
!!!!!!!nsameproc=0
!!!!!!!testLoop: do
!!!!!!!    do jproc=0,nproc-1
!!!!!!!        do jorb=1,4*orbs%norb_par(jproc,0)
!!!!!!!            if(comrp%communComplete(jorb,jproc)) cycle
!!!!!!!            call mpi_test(comrp%comarr(7,jorb,jproc), sendComplete, stat, ierr)
!!!!!!!            call mpi_test(comrp%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
!!!!!!!            ! Attention: mpi_test is a local function.
!!!!!!!            if(sendComplete .and. receiveComplete) comrp%communComplete(jorb,jproc)=.true.
!!!!!!!            !!if(comrp%communComplete(jorb,jproc)) then
!!!!!!!            !!    !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
!!!!!!!            !!    mpisource=comrp%comarr(1,jorb,jproc)
!!!!!!!            !!    mpidest=comrp%comarr(4,jorb,jproc)
!!!!!!!            !!    if(mpisource/=mpidest) then
!!!!!!!            !!        nfast=nfast+1
!!!!!!!            !!    else
!!!!!!!            !!        nsameproc=nsameproc+1
!!!!!!!            !!    end if
!!!!!!!            !!end if
!!!!!!!        end do
!!!!!!!    end do
!!!!!!!    ! If we made it until here, either all all the communication is
!!!!!!!    ! complete or we better wait for each single orbital.
!!!!!!!    exit testLoop
!!!!!!!end do testLoop
!!!!!!!
!!!!!!!! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!!!!!!call mpiallred(comrp%communComplete(1,0), nproc*4*maxval(orbs%norb_par(:,0)), mpi_land, mpi_comm_world, ierr)
!!!!!!!
!!!!!!!! Wait for the communications that have not completed yet
!!!!!!!nslow=0
!!!!!!!do jproc=0,nproc-1
!!!!!!!    do jorb=1,4*orbs%norb_par(jproc,0)
!!!!!!!        if(comrp%communComplete(jorb,jproc)) then
!!!!!!!            mpisource=comrp%comarr(1,jorb,jproc)
!!!!!!!            mpidest=comrp%comarr(4,jorb,jproc)
!!!!!!!            if(mpisource==mpidest) then
!!!!!!!                nsameproc=nsameproc+1
!!!!!!!            else
!!!!!!!                nfast=nfast+1
!!!!!!!            end if
!!!!!!!            cycle
!!!!!!!        end if
!!!!!!!        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
!!!!!!!        nslow=nslow+1
!!!!!!!        call mpi_wait(comrp%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!!!!!        call mpi_wait(comrp%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!!!!!        comrp%communComplete(jorb,jproc)=.true.
!!!!!!!    end do
!!!!!!!end do
!!!!!!!
!!!!!!!!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!!!!!                       nfast, ' could be overlapped with computation.'
!!!!!!!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!!
!!!
!!!
!!!
!!!
!!!end subroutine gatherDerivativeOrbitals
