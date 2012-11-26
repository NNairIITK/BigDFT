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


!!subroutine get_divergence(ndim, hgrid, lzd, lorbs, phider, phidiv)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => get_divergence
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: ndim
!!real(kind=8),intent(in) :: hgrid
!!type(local_zone_descriptors),intent(in) :: lzd
!!type(orbitals_data),intent(in) :: lorbs
!!real(kind=8),dimension(3*ndim),intent(in) :: phider !< Derivative Basis functions
!!real(kind=8),dimension(ndim),intent(inout) :: phidiv  !< Divergence of basis functions
!!
!!! Local variables
!!integer :: ist1_c, ist1_f, istdiv, nf, istat, iall, iorb, jproc
!!integer :: jlr, offset, ilr, iiorb, idir, iorbsmall
!!real(kind=8),dimension(0:3),parameter :: scal=1.d0
!!real(kind=8),dimension(:),allocatable :: w_f1, w_f2, w_f3
!!real(kind=8),dimension(:),allocatable :: phiX, phiY,phiZ
!!real(kind=8),dimension(:,:,:),allocatable :: w_c, phix_c, phiy_c, phiz_c
!!real(kind=8),dimension(:,:,:,:),allocatable :: w_f, phix_f, phiy_f, phiz_f
!!character(len=*),parameter :: subname='get_divergence'
!!
!!  if(mod(3*lorbs%isorb,3)+1 .ne. 1) STOP 'Derivative is not starting with X.'
!!  if(mod(3*lorbs%norbp+3*lorbs%isorb-1,3)+1 .ne. 3) STOP 'Derivative is not ending with Z.'
!!
!!  ist1_c=1
!!  istdiv=1
!!  do iorb=1,3*lorbs%norbp
!!
!!      iiorb=3*lorbs%isorb+iorb
!!      idir=mod(iiorb-1,3) + 1 ! get direction: x=1, y=2 or z=3 
!!      iorbsmall=ceiling(dble(iiorb)/3.d0)
!!      ilr=lorbs%inWhichLocreg(iorbsmall)
!!
!!      call allocateWorkarrays()
!!
!!      ist1_f=ist1_c+lzd%llr(ilr)%wfd%nvctr_c
!!
!!      ! Uncompress the wavefunction.
!!      call uncompress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, & 
!!           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!           scal, phider(ist1_c), phider(ist1_f), w_c, w_f, w_f1, w_f2, w_f3)
!!
!!      ist1_c = ist1_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f 
!!
!!      if(idir==1) then
!!          allocate(phiX(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f),stat=istat)
!!          call memocc(istat,phiX,'phiX',subname)
!!          allocate(phiY(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f),stat=istat)
!!          call memocc(istat,phiY,'phiY',subname)
!!          allocate(phiZ(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f),stat=istat)
!!          call memocc(istat,phiZ,'phiZ',subname)
!!
!!          call createDerivativeX(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
!!               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
!!               hgrid, lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibyz_f, &
!!               w_c, w_f, w_f1, phix_c, phix_f)
!!          ! Compress the x wavefunction.
!!          call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
!!               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
!!               lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
!!               lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!               lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!               lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!               scal, phix_c, phix_f, phiX(1), phiX(1+lzd%llr(ilr)%wfd%nvctr_c))
!!      else if (idir==2) then
!!         call createDerivativeY(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
!!               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
!!               hgrid, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxz_f, &
!!               w_c, w_f, w_f2, phiy_c, phiy_f)
!!         ! Compress the y wavefunction.
!!         call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!              lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
!!              lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
!!              lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
!!              lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!              lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!              lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!              scal, phiy_c, phiy_f, phiY(1), phiY(1+lzd%llr(ilr)%wfd%nvctr_c))
!!      else
!!         call createDerivativeZ(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
!!               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
!!               hgrid, lzd%llr(ilr)%bounds%kb%ibxy_c, lzd%llr(ilr)%bounds%kb%ibxy_f, &
!!               w_c, w_f, w_f3, phiz_c, phiz_f)
!!         ! Compress the z wavefunction.
!!         call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!              lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
!!              lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
!!              lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
!!              lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!              lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!              lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!              scal, phiz_c, phiz_f, phiZ(1),phiZ(1+lzd%llr(ilr)%wfd%nvctr_c ))
!!
!!         !Accumulate the divergence 
!!         call daxpy(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f, 1.0_dp, phiY, 1, phiX,1)
!!         call daxpy(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f, 1.0_dp, phiZ, 1, phiX,1)
!!         call dcopy(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f,phiX,1,phidiv(istdiv),1)
!!         istdiv = istdiv + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!
!!         iall = -product(shape(phiX))*kind(phiX)
!!         deallocate(phiX,stat=istat)
!!         call memocc(istat,iall,'phiX',subname)
!!         iall = -product(shape(phiY))*kind(phiY)
!!         deallocate(phiY,stat=istat)
!!         call memocc(istat,iall,'phiY',subname)
!!         iall = -product(shape(phiZ))*kind(phiZ)
!!         deallocate(phiZ,stat=istat)
!!         call memocc(istat,iall,'phiZ',subname)
!!      end if
!!
!!      call deallocateWorkarrays()                                
!!
!!  end do
!!
!!contains
!!  subroutine allocateWorkarrays()
!!
!!    ! THIS IS COPIED FROM allocate_work_arrays. Works only for free boundary.
!!    nf=(lzd%llr(ilr)%d%nfu1-lzd%llr(ilr)%d%nfl1+1)*(lzd%llr(ilr)%d%nfu2-lzd%llr(ilr)%d%nfl2+1)* &
!!       (lzd%llr(ilr)%d%nfu3-lzd%llr(ilr)%d%nfl3+1)
!!
!!    ! Allocate work arrays
!!    allocate(w_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3+ndebug), stat=istat)
!!    call memocc(istat, w_c, 'w_c', subname)
!!    !!w_c=0.d0
!!    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), w_c(0,0,0))
!!
!!    allocate(w_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
!!                 lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3+ndebug), stat=istat)
!!    call memocc(istat, w_f, 'w_f', subname)
!!    !!w_f=0.d0
!!    call to_zero(7*nf, w_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))
!!
!!
!!    allocate(w_f1(nf+ndebug), stat=istat)
!!    call memocc(istat, w_f1, 'w_f1', subname)
!!    !!w_f1=0.d0
!!    call to_zero(nf, w_f1(1))
!!
!!    allocate(w_f2(nf+ndebug), stat=istat)
!!    call memocc(istat, w_f2, 'w_f2', subname)
!!    !!w_f2=0.d0
!!    call to_zero(nf, w_f2(1))
!!
!!    allocate(w_f3(nf+ndebug), stat=istat)
!!    call memocc(istat, w_f3, 'w_f3', subname)
!!    !!w_f3=0.d0
!!    call to_zero(nf, w_f3(1))
!!
!!
!!    allocate(phix_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
!!                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
!!    call memocc(istat, phix_f, 'phix_f', subname)
!!    !!phix_f=0.d0
!!    call to_zero(7*nf, phix_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))
!!
!!    allocate(phix_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
!!    call memocc(istat, phix_c, 'phix_c', subname)
!!    !!phix_c=0.d0
!!    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phix_c(0,0,0))
!!
!!    allocate(phiy_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
!!                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
!!    call memocc(istat, phiy_f, 'phiy_f', subname)
!!    !!phiy_f=0.d0
!!    call to_zero(7*nf, phiy_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))
!!
!!    allocate(phiy_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
!!    call memocc(istat, phiy_c, 'phiy_c', subname)
!!    !!phiy_c=0.d0
!!    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiy_c(0,0,0))
!!
!!    allocate(phiz_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
!!                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
!!    call memocc(istat, phiz_f, 'phiz_f', subname)
!!    !!phiz_f=0.d0
!!    call to_zero(7*nf, phiz_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))
!!
!!    allocate(phiz_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
!!    call memocc(istat, phiz_c, 'phiz_c', subname)
!!    !!phiz_c=0.d0
!!    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiz_c(0,0,0))
!!
!!  end subroutine allocateWorkarrays
!!
!!
!!  subroutine deallocateWorkarrays
!!
!!    iall=-product(shape(w_c))*kind(w_c)
!!    deallocate(w_c, stat=istat)
!!    call memocc(istat, iall, 'w_c', subname)
!!
!!    iall=-product(shape(w_f))*kind(w_f)
!!    deallocate(w_f, stat=istat)
!!    call memocc(istat, iall, 'w_f', subname)
!!
!!    iall=-product(shape(w_f1))*kind(w_f1)
!!    deallocate(w_f1, stat=istat)
!!    call memocc(istat, iall, 'w_f1', subname)
!!
!!    iall=-product(shape(w_f2))*kind(w_f2)
!!    deallocate(w_f2, stat=istat)
!!    call memocc(istat, iall, 'w_f2', subname)
!!
!!    iall=-product(shape(w_f3))*kind(w_f3)
!!    deallocate(w_f3, stat=istat)
!!    call memocc(istat, iall, 'w_f3', subname)
!!
!!    iall=-product(shape(phix_f))*kind(phix_f)
!!    deallocate(phix_f, stat=istat)
!!    call memocc(istat, iall, 'phix_f', subname)
!!
!!    iall=-product(shape(phix_c))*kind(phix_c)
!!    deallocate(phix_c, stat=istat)
!!    call memocc(istat, iall, 'phix_c', subname)
!!
!!    iall=-product(shape(phiy_f))*kind(phiy_f)
!!    deallocate(phiy_f, stat=istat)
!!    call memocc(istat, iall, 'phiy_f', subname)
!!
!!    iall=-product(shape(phiy_c))*kind(phiy_c)
!!    deallocate(phiy_c, stat=istat)
!!    call memocc(istat, iall, 'phiy_c', subname)
!!
!!    iall=-product(shape(phiz_f))*kind(phiz_f)
!!    deallocate(phiz_f, stat=istat)
!!    call memocc(istat, iall, 'phiz_f', subname)
!!
!!    iall=-product(shape(phiz_c))*kind(phiz_c)
!!    deallocate(phiz_c, stat=istat)
!!    call memocc(istat, iall, 'phiz_c', subname)
!!
!!  end subroutine deallocateWorkarrays
!!
!!end subroutine get_divergence
