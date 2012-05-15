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

