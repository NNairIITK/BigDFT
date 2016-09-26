! Count for each orbital and each process the number of overlapping orbitals.
!!subroutine countOverlapsSphere(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(local_zone_descriptors),intent(in):: lzd
!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!  type(overlapParameters),intent(out):: op
!!  type(p2pComms),intent(out):: comon
!!
!!  ! Local variables
!!  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, iiorb
!!  integer ::  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
!!  logical:: ovrlpx, ovrlpy, ovrlpz
!!  real(8):: dx, dy, dz, rr
!!
!!  iiorb=0
!!  do jproc=0,nproc-1
!!     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!     ilrold=-1
!!     do iorb=1,orbs%norb_par(jproc,0)
!!        ioverlaporb=0 ! counts the overlaps for the given orbital.
!!        iiorb=iiorb+1 ! counts the total orbitals
!!        ilr=onWhichAtom(iiorb)
!!        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!        do jorb=1,orbs%norb
!!           jlr=onWhichAtom(jorb)
!!           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
!!           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
!!           !if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!           dx=(lzd%llr(ilr)%locregCenter(1)-lzd%llr(jlr)%locregCenter(1))**2
!!           dy=(lzd%llr(ilr)%locregCenter(2)-lzd%llr(jlr)%locregCenter(2))**2
!!           dz=(lzd%llr(ilr)%locregCenter(3)-lzd%llr(jlr)%locregCenter(3))**2
!!           rr=(lzd%llr(ilr)%locrad+lzd%llr(jlr)%locrad)**2
!!           if(dx+dy+dz<=rr) then
!!              ioverlaporb=ioverlaporb+1
!!              if(ilr/=ilrold) then
!!                 ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
!!                 ! would get the same orbitals again. Therefore the counter is not increased
!!                 ! in that case.
!!                 ioverlapMPI=ioverlapMPI+1
!!              end if
!!           end if
!!        end do
!!        op%noverlaps(iiorb)=ioverlaporb
!!        !!if(iproc==0) write(*,'(a,2i8)') 'iiorb, op%noverlaps(iiorb)', iiorb, op%noverlaps(iiorb)
!!        ilrold=ilr
!!     end do
!!     comon%noverlaps(jproc)=ioverlapMpi
!!     !if(iproc==0) write(*,'(a,2i8)') 'jproc, comon%noverlaps(jproc)', jproc, comon%noverlaps(jproc)
!!  end do
!!
!!end subroutine countOverlapsSphere

!!subroutine determineOverlapDescriptors(iproc, nproc, orbs, lzd, Glr, onWhichAtom, op)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(locreg_descriptors),intent(in):: Glr
!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!  type(overlapParameters),intent(inout):: op
!!
!!  ! Local variables
!!  integer:: iorb, jorb, jjorb, ilr, jlr, iiorb
!!
!!
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb_par(iproc)+iorb
!!     ilr=onWhichAtom(iiorb)
!!     !if(iproc==0) write(*,'(a,2i10)') 'iorb, op%noverlaps(iorb)', iorb, op%noverlaps(iorb)
!!     do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        jlr=onWhichAtom(jjorb)
!!        !write(*,*) 'calling get_overlap_region_periodic'
!!        !call get_overlap_region_periodic(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iorb))
!!        call get_overlap_region_periodic2(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iorb))
!!        !write(*,'(a,13i8)') 'iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, nvctr_c, nvctr_f, ncount, n1, n2, n3', iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, &
!!        !    op%olr(jorb,iorb)%wfd%nvctr_c, op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%d%n1, op%olr(jorb,iorb)%d%n2, op%olr(jorb,iorb)%d%n3
!!     end do
!!  end do
!!
!!end subroutine determineOverlapDescriptors




!!subroutine determineOverlapDescriptorsSphere(iproc, nproc, orbs, lzd, Glr, onWhichAtom, hx, hy, hz, op)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(locreg_descriptors),intent(in):: Glr
!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!  real(8):: hx, hy, hz
!!  type(overlapParameters),intent(inout):: op
!!
!!  ! Local variables
!!  integer:: iorb, jorb, jjorb, ilr, jlr, iiorb
!!
!!
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb_par(iproc)+iorb
!!     ilr=onWhichAtom(iiorb)
!!     !if(iproc==0) write(*,'(a,2i10)') 'iorb, op%noverlaps(iorb)', iorb, op%noverlaps(iorb)
!!     do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        jlr=onWhichAtom(jjorb)
!!        call determine_overlapdescriptors_from_descriptors(lzd%llr(ilr), lzd%llr(jlr), lzd%glr, op%olr(jorb,iorb))
!!        !write(*,'(a,13i8)') 'iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, nvctr_c, nvctr_f, ncount, n1, n2, n3', iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, &
!!        !    op%olr(jorb,iorb)%wfd%nvctr_c, op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%d%n1, op%olr(jorb,iorb)%d%n2, op%olr(jorb,iorb)%d%n3
!!     end do
!!  end do
!!
!!end subroutine determineOverlapDescriptorsSphere

!!subroutine gatherOrbitalsOverlapWithComput(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp, expanded)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => gatherOrbitalsOverlapWithComput
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(input_variables),intent(in):: input
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(in):: op
!!  type(p2pComms),intent(inout):: comon
!!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!  logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded
!!
!!  ! Local variables
!!  integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc, jjorb, orbsource, orbdest
!!  integer,dimension(mpi_status_size):: stat
!!  logical:: sendComplete, receiveComplete
!!
!!
!!  ! Check whether the communications have completed.
!!  nfast=0
!!  nsameproc=0
!!  testLoop: do
!!     do jproc=0,nproc-1
!!        do jorb=1,comon%noverlaps(jproc)
!!           mpisource=comon%comarr(1,jorb,jproc)
!!           mpidest=comon%comarr(4,jorb,jproc)
!!           orbsource=comon%comarr(9,jorb,jproc)
!!           orbdest=comon%comarr(10,jorb,jproc)
!!           if(mpisource==mpidest) then
!!              if(iproc==mpidest) then
!!                 call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
!!                      orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!                 expanded(orbsource,orbdest-orbs%isorb)=.true.
!!              end if
!!           end if
!!           if(comon%communComplete(jorb,jproc)) cycle
!!           ! Attention: mpi_test is a local function.
!!           call mpi_test(comon%comarr(7,jorb,jproc), sendComplete, stat, ierr)
!!           call mpi_test(comon%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
!!           if(sendComplete .and. receiveComplete) comon%communComplete(jorb,jproc)=.true.
!!           if(comon%communComplete(jorb,jproc)) then
!!              !if(iproc==jproc) write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
!!              mpisource=comon%comarr(1,jorb,jproc)
!!              mpidest=comon%comarr(4,jorb,jproc)
!!              orbsource=comon%comarr(9,jorb,jproc)
!!              orbdest=comon%comarr(10,jorb,jproc)
!!              if(iproc==mpidest) then
!!                 call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
!!                      orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!                 expanded(orbsource,orbdest-orbs%isorb)=.true.
!!              end if
!!              if(mpisource/=mpidest) then
!!                 !nfast=nfast+1
!!              else
!!                 !nsameproc=nsameproc+1
!!              end if
!!           end if
!!        end do
!!     end do
!!     ! If we made it until here, either all all the communication is
!!     ! complete or we better wait for each single orbital.
!!     exit testLoop
!!  end do testLoop
!!
!!  ! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!  call mpiallred(comon%communComplete(1,0), nproc*maxval(comon%noverlaps), mpi_land, mpi_comm_world, ierr)
!!
!!  ! Wait for the communications that have not completed yet
!!  nslow=0
!!  do jproc=0,nproc-1
!!     do jorb=1,comon%noverlaps(jproc)
!!        !!if(comon%communComplete(jorb,jproc)) cycle
!!        if(comon%communComplete(jorb,jproc)) then
!!           mpisource=comon%comarr(1,jorb,jproc)
!!           mpidest=comon%comarr(4,jorb,jproc)
!!           if(mpisource==mpidest) then
!!              nsameproc=nsameproc+1
!!           else
!!              nfast=nfast+1
!!           end if
!!           cycle
!!        end if
!!        !write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
!!        nslow=nslow+1
!!        call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!        call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!        comon%communComplete(jorb,jproc)=.true.
!!        mpidest=comon%comarr(4,jorb,jproc)
!!        orbsource=comon%comarr(9,jorb,jproc)
!!        orbdest=comon%comarr(10,jorb,jproc)
!!        if(iproc==mpidest) then
!!           !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!           expanded(orbsource,orbdest-orbs%isorb)=.false.
!!        end if
!!
!!
!!
!!
!!
!!
!!!!!write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
!!        !!nslow=nslow+1
!!        !!call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!        !!call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!        !!comon%communComplete(jorb,jproc)=.true.
!!        !!mpisource=comon%comarr(1,jorb,jproc)
!!        !!mpidest=comon%comarr(4,jorb,jproc)
!!        !!orbsource=comon%comarr(9,jorb,jproc)
!!        !!orbdest=comon%comarr(10,jorb,jproc)
!!        !!if(iproc==mpidest) then
!!        !!    !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!        !!    expanded(orbsource,orbdest-orbs%isorb)=.false.
!!        !!end if
!!!!!write(*,'(3(a,i0))') 'process ', iproc, ' has finally received orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
!!     end do
!!  end do
!!
!!  !call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!  !call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!  !call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!  !call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!  if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!       nfast, ' could be overlapped with computation.'
!!  if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!
!!
!!end subroutine gatherOrbitalsOverlapWithComput







!!subroutine gatherOrbitalsOverlapWithComput2(iproc, nproc, orbs, input, lzd, op, comon,&
!!     nsendbuf, sendbuf, nrecvbuf, recvbuf, lphiovrlp, expanded, ovrlp)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => gatherOrbitalsOverlapWithComput2
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
!!  type(orbitals_data),intent(in):: orbs
!!  type(input_variables),intent(in):: input
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(in):: op
!!  type(p2pComms),intent(inout):: comon
!!  real(8),dimension(nsendbuf),intent(in):: sendbuf
!!  real(8),dimension(nrecvbuf),intent(in):: recvbuf
!!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!  logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded
!!  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
!!
!!  ! Local variables
!!  integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc, jjorb, orbsource, orbdest, ncount, ist, jst
!!  integer,dimension(mpi_status_size):: stat
!!  logical:: sendComplete, receiveComplete
!!  real(8):: ddot
!!
!!  if(iproc==0) write(*,*) 'new subroutine'
!!
!!  ! Check whether the communications have completed. Only check the receives.
!!  nfast=0
!!  nsameproc=0
!!  testLoop: do
!!     do jproc=0,nproc-1
!!        do jorb=1,comon%noverlaps(jproc)
!!           mpisource=comon%comarr(1,jorb,jproc)
!!           mpidest=comon%comarr(4,jorb,jproc)
!!           orbsource=comon%comarr(9,jorb,jproc)
!!           orbdest=comon%comarr(10,jorb,jproc)
!!           if(iproc==mpidest) then
!!              if(mpisource==mpidest) then
!!                 call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
!!                      orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!                 expanded(orbsource,orbdest-orbs%isorb)=.true.
!!                 ! Calculate matrix element here.
!!                 call getStartingIndicesGlobal(orbsource, orbdest, op, orbs, ist, jst, ncount)
!!                 ovrlp(orbdest,orbsource)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
!!              end if
!!              if(comon%communComplete(jorb,jproc)) cycle
!!              ! Attention: mpi_test is a local function.
!!              ! Test whether the receive has completed.
!!              !call mpi_test(comon%comarr(7,jorb,jproc), sendComplete, stat, ierr)
!!              call mpi_test(comon%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
!!              !if(sendComplete .and. receiveComplete) comon%communComplete(jorb,jproc)=.true.
!!              if(receiveComplete) comon%communComplete(jorb,jproc)=.true.
!!              if(comon%communComplete(jorb,jproc)) then
!!                 !if(iproc==jproc) write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
!!                 mpisource=comon%comarr(1,jorb,jproc)
!!                 mpidest=comon%comarr(4,jorb,jproc)
!!                 orbsource=comon%comarr(9,jorb,jproc)
!!                 orbdest=comon%comarr(10,jorb,jproc)
!!                 if(iproc==mpidest) then
!!                    call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
!!                         orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!                    expanded(orbsource,orbdest-orbs%isorb)=.true.
!!                    ! Calculate matrix element here.
!!                    call getStartingIndicesGlobal(orbsource, orbdest, op, orbs, ist, jst, ncount)
!!                    ovrlp(orbdest,orbsource)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
!!                 end if
!!                 if(mpisource/=mpidest) then
!!                    !nfast=nfast+1
!!                 else
!!                    !nsameproc=nsameproc+1
!!                 end if
!!              end if
!!           end if
!!        end do
!!     end do
!!     ! If we made it until here, either all all the communication is
!!     ! complete or we better wait for each single orbital.
!!     exit testLoop
!!  end do testLoop
!!
!!  !!call mpi_barrier(mpi_comm_world, ierr)
!!  !!if(iproc==0) write(*,*) 'after test loop'
!!  !! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!  !call mpiallred(comon%communComplete(1,0), nproc*maxval(comon%noverlaps), mpi_land, mpi_comm_world, ierr)
!!
!!  ! Wait for the communications that have not completed yet
!!  nslow=0
!!  do jproc=0,nproc-1
!!     do jorb=1,comon%noverlaps(jproc)
!!        mpisource=comon%comarr(1,jorb,jproc)
!!        mpidest=comon%comarr(4,jorb,jproc)
!!        if(iproc==mpidest) then
!!           ! Only check the receive
!!           !!if(comon%communComplete(jorb,jproc)) cycle
!!           if(comon%communComplete(jorb,jproc)) then
!!              if(mpisource==mpidest) then
!!                 nsameproc=nsameproc+1
!!              else
!!                 nfast=nfast+1
!!              end if
!!              cycle
!!           end if
!!           !write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
!!           nslow=nslow+1
!!           !call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!           call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!           comon%communComplete(jorb,jproc)=.true.
!!           mpidest=comon%comarr(4,jorb,jproc)
!!           orbsource=comon%comarr(9,jorb,jproc)
!!           orbdest=comon%comarr(10,jorb,jproc)
!!           if(iproc==mpidest) then
!!              !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!              expanded(orbsource,orbdest-orbs%isorb)=.false.
!!              call getStartingIndicesGlobal(orbsource, orbdest, op, orbs, ist, jst, ncount)
!!              ovrlp(orbdest,orbsource)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
!!           end if
!!
!!        else if(iproc==mpisource) then
!!           ! Check the send
!!           call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)
!!        end if
!!
!!
!!
!!
!!
!!
!!!!!write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
!!        !!nslow=nslow+1
!!        !!call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!        !!call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!        !!comon%communComplete(jorb,jproc)=.true.
!!        !!mpisource=comon%comarr(1,jorb,jproc)
!!        !!mpidest=comon%comarr(4,jorb,jproc)
!!        !!orbsource=comon%comarr(9,jorb,jproc)
!!        !!orbdest=comon%comarr(10,jorb,jproc)
!!        !!if(iproc==mpidest) then
!!        !!    !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
!!        !!    expanded(orbsource,orbdest-orbs%isorb)=.false.
!!        !!end if
!!!!!write(*,'(3(a,i0))') 'process ', iproc, ' has finally received orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
!!     end do
!!  end do
!!
!!  !call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!  !call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!  !call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!  !call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!  if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!       nfast, ' could be overlapped with computation.'
!!  if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!
!!
!!end subroutine gatherOrbitalsOverlapWithComput2

!!subroutine expandRemainingOrbitals(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, expanded, lphiovrlp)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(input_variables),intent(in):: input
!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(in):: op
!!  type(p2pComms),intent(in):: comon
!!  logical,dimension(orbs%norb,orbs%norbp),intent(in):: expanded
!!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!
!!  ! Local variables
!!  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest
!!
!!
!!  !lphiovrlp=0.d0
!!
!!  ind=1
!!  ilrold=-1
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb+iorb
!!     ilr=onWhichAtom(iiorb)
!!     if(ilr==ilrold) cycle
!!     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!     do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        ! Starting index of orbital jjorb
!!        jst=op%indexInRecvBuf(iorb,jjorb)
!!        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
!!        ldim=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
!!        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
!!        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
!!        if(.not. expanded(jjorb,iorb)) then
!!           do i=0,ldim-1
!!              indDest=ind+op%indexExpand(jst+i)-1
!!              lphiovrlp(indDest)=comon%recvBuf(jst+i)
!!           end do
!!        end if
!!        ind=ind+gdim
!!     end do
!!     ilrold=ilr
!!  end do
!!
!!end subroutine expandRemainingOrbitals



!!subroutine expandOneOrbital(iproc, nproc, orbsource, orbdest, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc, orbsource, orbdest
!!  type(orbitals_data),intent(in):: orbs
!!  type(input_variables),intent(in):: input
!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(in):: op
!!  type(p2pComms),intent(in):: comon
!!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!
!!  ! Local variables
!!  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest
!!
!!
!!  !lphiovrlp=0.d0
!!
!!  ind=1
!!  ilrold=-1
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb+iorb
!!     ilr=onWhichAtom(iiorb)
!!     if(ilr==ilrold) cycle
!!     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!     do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        ! Starting index of orbital jjorb
!!        jst=op%indexInRecvBuf(iorb,jjorb)
!!        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
!!        ldim=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
!!        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
!!        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
!!        if(jjorb==orbsource .and. iorb==orbdest) then
!!           do i=0,ldim-1
!!              indDest=ind+op%indexExpand(jst+i)-1
!!              lphiovrlp(indDest)=comon%recvBuf(jst+i)
!!           end do
!!        end if
!!        ind=ind+gdim
!!     end do
!!     ilrold=ilr
!!  end do
!!
!!end subroutine expandOneOrbital



!!subroutine expandOneOrbital2(iproc, nproc, orbsource, orbdest, orbs, input, onWhichAtom, lzd, op, nrecvbuf, recvbuf, lphiovrlp)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc, orbsource, orbdest, nrecvbuf
!!  type(orbitals_data),intent(in):: orbs
!!  type(input_variables),intent(in):: input
!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(in):: op
!!  real(8),dimension(nrecvbuf),intent(in):: recvbuf
!!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!
!!  ! Local variables
!!  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest
!!
!!
!!  ind=1
!!  ilrold=-1
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb+iorb
!!     ilr=onWhichAtom(iiorb)
!!     if(ilr==ilrold) cycle
!!     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!     do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        ! Starting index of orbital jjorb
!!        jst=op%indexInRecvBuf(iorb,jjorb)
!!        ldim=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
!!        if(jjorb==orbsource .and. iorb==orbdest) then
!!           do i=0,ldim-1
!!              indDest=ind+op%indexExpand(jst+i)-1
!!              lphiovrlp(indDest)=recvbuf(jst+i)
!!              !lphiovrlp(indDest)=recvbuf(1+i)
!!           end do
!!        end if
!!        ind=ind+gdim
!!     end do
!!     ilrold=ilr
!!  end do
!!
!!end subroutine expandOneOrbital2

!!subroutine checkUnity(iproc, norb, ovrlp, maxError)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, norb
!!  real(8),dimension(norb,norb),intent(in):: ovrlp
!!  real(8),intent(out):: maxError
!!
!!  ! Local variables
!!  integer:: iorb, jorb
!!  real(8):: error
!!
!!  maxError=0.d0
!!  do iorb=1,norb
!!     do jorb=1,norb
!!        if(iorb==jorb) then
!!           error=abs(ovrlp(jorb,iorb)-1.d0)
!!        else
!!           error=abs(ovrlp(jorb,iorb))
!!        end if
!!        if(error>maxError) then
!!           maxError=error
!!        end if
!!        !write(20000+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!     end do
!!  end do
!!
!!end subroutine checkUnity

!!subroutine indicesForExpansion(iproc, nproc, nspin, orbs, onWhichAtom, lzd, op, comon)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => indicesForExpansion
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, nspin
!!type(orbitals_data),intent(in):: orbs
!!integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(inout):: op
!!type(p2pComms),intent(in):: comon
!!
!!! Local variables
!!integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, istat
!!character(len=*),parameter:: subname='indicesForExpansion'
!!
!!
!!ind=1
!!ilrold=-1
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=onWhichAtom(iiorb)
!!    if(ilr==ilrold) cycle
!!    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        ! Starting index of orbital jjorb
!!        jst=op%indexInRecvBuf(iorb,jjorb)
!!        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
!!        ldim=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
!!        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
!!        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
!!        call index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, nspin, &
!!             lzd%llr(ilr), op%olr(jorb,iorb), op%indexExpand(jst))
!!
!!        call countExpansionSegments(ldim, op%indexExpand(jst), op%expseg(jorb,iorb)%nseg)
!!        allocate(op%expseg(jorb,iorb)%segborders(2,op%expseg(jorb,iorb)%nseg), stat=istat)
!!        call memocc(istat, op%expseg(jorb,iorb)%segborders, 'op%expseg(jorb,iorb)%segborders', subname)
!!        call determineExpansionSegments(ldim, op%indexExpand(jst), op%expseg(jorb,iorb)%nseg, op%expseg(jorb,iorb)%segborders)
!!
!!        ind=ind+gdim
!!    end do
!!    ilrold=ilr
!!end do
!!
!!end subroutine indicesForExpansion



!!subroutine countExpansionSegments(ldim, indexExpand, nseg)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: ldim
!!integer,dimension(ldim),intent(in):: indexExpand
!!integer,intent(out):: nseg
!!integer,dimension(:,:),pointer:: segborders
!!
!!! Local variables
!!integer:: i
!!
!!! Count the numbers of segments
!!nseg=0
!!do i=2,ldim
!!    if(indexExpand(i)==indexExpand(i-1)+1) then
!!        !consecutive, same segment
!!    else
!!        !segment ended
!!        nseg=nseg+1
!!    end if
!!end do
!!nseg=nseg+1 !last segment
!!
!!end subroutine countExpansionSegments


!!subroutine determineExpansionSegments(ldim, indexExpand, nseg, segborders)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: ldim, nseg
!!integer,dimension(ldim),intent(in):: indexExpand
!!integer,dimension(2,nseg):: segborders
!!
!!! Local variables
!!integer:: iseg, i
!!character(len=*),parameter:: subname='determineExpansionSegments'
!!
!!iseg=1
!!segborders(1,iseg)=indexExpand(1)
!!do i=2,ldim
!!    if(indexExpand(i)==indexExpand(i-1)+1) then
!!        !consecutive, same segment
!!    else
!!        !segment ended
!!        segborders(2,iseg)=indexExpand(i-1)
!!        iseg=iseg+1
!!        segborders(1,iseg)=indexExpand(i)
!!    end if
!!end do
!!segborders(2,iseg)=indexExpand(ldim) !last segments
!!
!!
!!
!!end subroutine determineExpansionSegments




!!subroutine indicesForExtraction(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, comon)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => indicesForExtraction
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, sizePhi
!!type(orbitals_data),intent(in):: orbs
!!integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(inout):: op
!!type(p2pComms),intent(out):: comon
!!
!!! Local variables
!!integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
!!integer:: i, istat
!!character(len=*),parameter:: subname='indicesForExtraction'
!!
!!indovrlp=1
!!op%indexInSendBuf=0
!!
!!ilrold=-1
!!do iorb=1,orbs%norb
!!    ilr=onWhichAtom(iorb)
!!    iiproc=orbs%onWhichMPI(iorb)
!!    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
!!    do jorb=1,op%noverlaps(iorb)
!!        jjorb=op%overlaps(jorb,iorb)
!!        jjlr=onWhichAtom(jjorb)
!!        jjproc=orbs%onWhichMPI(jjorb)
!!        if(iproc==jjproc) then
!!            ! Get the correct descriptors
!!            korb=jjorb-orbs%isorb
!!            !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
!!            do i=1,op%noverlaps(jjorb)
!!                !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
!!                if(op%overlaps(i,jjorb)==iorb) then
!!                    lorb=i
!!                    exit
!!                end if
!!            end do
!!            !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
!!            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
!!            ldim=op%wfd_overlap(lorb,korb)%nvctr_c+7*op%wfd_overlap(lorb,korb)%nvctr_f
!!            ind=1
!!            do kkorb=orbs%isorb+1,jjorb-1
!!                klr=onWhichAtom(kkorb)
!!                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
!!            end do
!!            !write(*,'(5(a,i0))') 'process ',iproc,' adds ',op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f,' elements at position ',indovrlp,' from orbital ',jjorb,' for orbital ', iorb
!!            !call psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), phi(ind), comon%sendBuf(indovrlp))
!!            call index_of_psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), op%indexExtract(indovrlp))
!!
!!            call countExpansionSegments(ldim, op%indexExtract(indovrlp), op%extseg(lorb,korb)%nseg)
!!            allocate(op%extseg(lorb,korb)%segborders(2,op%extseg(lorb,korb)%nseg), stat=istat)
!!            call memocc(istat, op%extseg(lorb,korb)%segborders, 'op%extseg(lorb,korb)%segborders', subname)
!!            call determineExpansionSegments(ldim, op%indexExtract(indovrlp), op%extseg(lorb,korb)%nseg, &
!!                 op%extseg(lorb,korb)%segborders)
!!
!!            op%indexInSendBuf(jjorb-orbs%isorb,iorb)=indovrlp
!!            indovrlp=indovrlp+op%wfd_overlap(lorb,korb)%nvctr_c+7*op%wfd_overlap(lorb,korb)%nvctr_f
!!        end if
!!    end do
!!    ilrold=ilr
!!    iiprocold=iiproc
!!end do
!!
!!if(indovrlp/=comon%nsendBuf+1) then
!!    write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
!!    stop
!!end if
!!
!!
!!
!!end subroutine indicesForExtraction

!!!subroutine initCommsOrthoVariable(iproc, nproc, lzd, orbs, orbsig, onWhichAtomAll, input, op, comon, tag)
!!!use module_base
!!!use module_types
!!!use module_interfaces, exceptThisOne => initCommsOrthoVariable
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(local_zone_descriptors),intent(in):: lzd
!!!type(orbitals_data),intent(in):: orbs, orbsig
!!!integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
!!!type(input_variables),intent(in):: input
!!!type(overlapParameters),intent(out):: op
!!!type(p2pComms),intent(out):: comon
!!!integer,intent(inout):: tag
!!!
!!!! Local variables
!!!integer:: iorb, jorb, iiorb, jproc, ioverlaporb, ioverlapMPI, ilr, jlr
!!!integer:: ilrold, is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3
!!!integer::  je3, istat, i1, i2
!!!logical:: ovrlpx, ovrlpy, ovrlpz
!!!character(len=*),parameter:: subname='initCommsOrthoVariable'
!!!
!!!allocate(comon%noverlaps(0:nproc-1), stat=istat)
!!!call memocc(istat, comon%noverlaps, 'comon%noverlaps',subname)
!!!allocate(op%noverlaps(orbs%norb), stat=istat)
!!!call memocc(istat, op%noverlaps, 'op%noverlaps',subname)
!!!
!!!! Count how many overlaping regions each orbital / process has.
!!!call countOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
!!!
!!!op%noverlapsmax=maxval(op%noverlaps)
!!!allocate(op%overlaps(op%noverlapsmax,orbs%norb), stat=istat)
!!!call memocc(istat, op%overlaps, 'op%overlaps', subname)
!!!comon%noverlapsmax=maxval(comon%noverlaps)
!!!!!allocate(comon%overlaps(comon%noverlapsmax,0:nproc-1), stat=istat)
!!!!!call memocc(istat, comon%overlaps, 'comon%overlaps', subname)
!!!allocate(op%indexInRecvBuf(orbs%norbp,orbsig%norb), stat=istat)
!!!call memocc(istat, op%indexInRecvBuf, 'op%indexInRecvBuf', subname)
!!!allocate(op%indexInSendBuf(orbsig%norbp,orbsig%norb), stat=istat)
!!!call memocc(istat, op%indexInSendBuf, 'op%indexInSendBuf', subname)
!!!
!!!! Determine the overlapping orbitals.
!!!call determineOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
!!!
!!!!!allocate(op%olr(op%noverlapsmax,orbs%norb), stat=istat)
!!!!!!do i2=1,orbs%norbp
!!!!!do i2=1,orbs%norb
!!!!!    do i1=1,op%noverlapsmax
!!!!!        call nullify_locreg_descriptors(op%olr(i1,i2))
!!!!!    end do
!!!!!end do
!!!
!!!! Set the orbital descriptors for the overlap regions.
!!!!call determineOverlapDescriptorsVariable(iproc, nproc, orbs, orbsig, lzd, lzd%Glr, onWhichAtomAll, op)
!!! call determineOverlapDescriptorsVariable2(orbs, orbsig, lzd, op)
!!!
!!!allocate(comon%comarr(10,comon%noverlapsmax,0:nproc-1), stat=istat)
!!!call memocc(istat, comon%comarr, 'comon%comarr', subname)
!!!allocate(comon%communComplete(comon%noverlapsmax,0:nproc-1), stat=istat)
!!!call memocc(istat, comon%communComplete, 'comun%communComplete', subname)
!!!call setCommsOrthoVariable(iproc, nproc, orbs, orbsig, lzd, op, comon, tag)
!!!
!!!
!!!!DONT NEED THIS ANYMORE
!!!!!! Initialize the index arrays for the transformations from overlap region
!!!!!! to ordinary localization region.
!!!!!allocate(op%indexExpand(comon%nrecvBuf), stat=istat)
!!!!!call memocc(istat, op%indexExpand, 'op%indexExpand',subname)
!!!!!call indicesForExpansionVariable(iproc, nproc, orbs, input, lzd, op, comon)
!!!!!
!!!!!! Initialize the index arrays for the transformations from the ordinary localization region
!!!!!! to the overlap region.
!!!!!allocate(op%indexExtract(comon%nsendBuf), stat=istat)
!!!!!call memocc(istat, op%indexExtract, 'op%indexExtract',subname)
!!!!!call indicesForExtractionVariable(iproc, nproc, orbs, orbsig, orbs%npsidim_orbs, lzd, op, comon)
!!!
!!!end subroutine initCommsOrthoVariable



!!!! Count for each orbital and each process the number of overlapping orbitals.
!!!subroutine countOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(orbitals_data),intent(in):: orbs, orbsig
!!!type(local_zone_descriptors),intent(in):: lzd
!!!type(overlapParameters),intent(out):: op
!!!type(p2pComms),intent(out):: comon
!!!
!!!! Local variables
!!!integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, iiorb
!!!logical:: isoverlap
!!!!integer:: is1, ie1, is2, ie2, is3, ie3
!!!!integer:: js1, je1, js2, je2, js3, je3
!!!!logical:: ovrlpx, ovrlpy, ovrlpz
!!!
!!!iiorb=0
!!!do jproc=0,nproc-1
!!!    ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!!    ilrold=-1
!!!    do iorb=1,orbs%norb_par(jproc,0)
!!!        ioverlaporb=0 ! counts the overlaps for the given orbital.
!!!        iiorb=iiorb+1 ! counts the total orbitals
!!!        ilr=orbs%inWhichLocreg(iiorb)
!!!!        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!!        do jorb=1,orbsig%norb
!!!            jlr=orbsig%inWhichLocreg(jorb)
!!!!            call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!!!            ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!!            ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!!            ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!!!            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!!             call check_overlap_cubic_periodic(lzd%Glr,lzd%Llr(ilr),lzd%Llr(jlr),isoverlap)
!!!             if(isoverlap) then
!!!                ioverlaporb=ioverlaporb+1
!!!                if(ilr/=ilrold) then
!!!                    ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
!!!                    ! would get the same orbitals again. Therefore the counter is not increased
!!!                    ! in that case.
!!!                    ioverlapMPI=ioverlapMPI+1
!!!                end if
!!!            end if
!!!        end do 
!!!        op%noverlaps(iiorb)=ioverlaporb
!!!        !!if(iproc==0) write(*,'(a,2i8)') 'iiorb, op%noverlaps(iiorb)', iiorb, op%noverlaps(iiorb)
!!!        ilrold=ilr
!!!    end do
!!!    comon%noverlaps(jproc)=ioverlapMpi
!!!    !if(iproc==0) write(*,'(a,2i8)') 'jproc, comon%noverlaps(jproc)', jproc, comon%noverlaps(jproc)
!!!end do
!!!
!!!end subroutine countOverlapsVariable




!!subroutine determineOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(orbitals_data),intent(in):: orbs, orbsig
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(out):: op
!!type(p2pComms),intent(out):: comon
!!
!!! Local variables
!!integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, iiorb
!!logical:: isoverlap
!!!integer :: is1, ie1, is2, ie2, is3, ie3
!!!integer:: js1, je1, js2, je2, js3, je3
!!!logical:: ovrlpx, ovrlpy, ovrlpz
!!
!!
!!  ! Initialize to some value which will never be used.
!!  op%overlaps=-1
!!  !!comon%overlaps=-1
!!
!!  iiorb=0
!!  do jproc=0,nproc-1
!!      ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!      ilrold=-1
!!      do iorb=1,orbs%norb_par(jproc,0)
!!          ioverlaporb=0 ! counts the overlaps for the given orbital.
!!          iiorb=iiorb+1 ! counts the total orbitals
!!          ilr=orbs%inWhichLocreg(iiorb)
!!!          call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!          do jorb=1,orbsig%norb
!!              jlr=orbsig%inWhichLocreg(jorb)
!!               call check_overlap_cubic_periodic(lzd%Glr,lzd%Llr(ilr),lzd%Llr(jlr),isoverlap)
!!!              call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!!              ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!              ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!              ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!!              if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!               if(isoverlap) then
!!                  ioverlaporb=ioverlaporb+1
!!                  op%overlaps(ioverlaporb,iiorb)=jorb
!!                  if(ilr/=ilrold) then
!!                      ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
!!                      ! would get the same orbitals again. Therefore the counter is not increased
!!                      ! in that case.
!!                      ioverlapMPI=ioverlapMPI+1
!!                      !!comon%overlaps(ioverlapMPI,jproc)=jorb
!!                  end if
!!              end if
!!          end do 
!!          !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
!!          ilrold=ilr
!!      end do
!!      !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
!!  end do
!!
!!end subroutine determineOverlapsVariable





!!subroutine determineOverlapDescriptorsVariable(iproc, nproc, orbs, orbsig, lzd, Glr, onWhichAtom, op)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(orbitals_data),intent(in):: orbs, orbsig
!!type(local_zone_descriptors),intent(in):: lzd
!!type(locreg_descriptors),intent(in):: Glr
!!integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!type(overlapParameters),intent(inout):: op
!!
!!! Local variables
!!integer:: iorb, jorb, jjorb, ilr, jlr, iiorb
!!
!!do iiorb=1,orbs%norb
!!    ilr=orbs%inWhichLocreg(iiorb)
!!    do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        jlr=orbsig%inWhichLocreg(jjorb)
!!        call get_overlap_region_periodic2(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iiorb))
!!    end do
!!end do
!!
!!end subroutine determineOverlapDescriptorsVariable

!!subroutine determineOverlapDescriptorsVariable2(orbs, orbsig, lzd, op)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!type(orbitals_data),intent(in):: orbs, orbsig
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(inout):: op
!!
!!! Local variables
!!integer:: iorb, jorb, jjorb, ilr, jlr, iiorb, istat
!!logical :: isoverlap,isoverlapf
!!character(len=*),parameter:: subname='determineOverlapDescriptorsVariable2'
!!
!!allocate(op%wfd_overlap(orbs%norb,orbs%norbp), stat=istat)
!!
!!do iiorb=1,orbs%norb
!!    ilr=orbs%inWhichLocreg(iiorb)
!!    do jorb=1,op%noverlaps(iiorb)
!!        call nullify_wavefunctions_descriptors(op%wfd_overlap(iiorb,jorb))
!!        jjorb=op%overlaps(jorb,iiorb)
!!        jlr=orbsig%inWhichLocreg(jjorb)
!!        !Get number of coarse segments
!!        call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_c,&
!!                 lzd%llr(ilr)%wfd%keyglob(1,1), lzd%llr(jlr)%wfd%keyglob(1,1), &
!!                 isoverlap, op%wfd_overlap(iiorb,jorb)%nseg_c)
!!        !Get number of fine segments
!!        call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_f,&
!!                 lzd%llr(ilr)%wfd%keyglob(1,1), lzd%llr(jlr)%wfd%keyglob(1,1+lzd%llr(jlr)%wfd%nseg_c), &
!!                 isoverlapf, op%wfd_overlap(iiorb,jorb)%nseg_f)
!!        ! Allocate the keys
!!        call allocate_wfd(op%wfd_overlap(iiorb,jorb),subname)
!!        ! Get the coarse part
!!        call get_overlap_from_descriptors_periodic(lzd%Llr(ilr)%wfd%nseg_c, lzd%Llr(jlr)%wfd%nseg_c,&
!!           lzd%Llr(ilr)%wfd%keyglob(1,1), lzd%Llr(jlr)%wfd%keyglob(1,1), &                                                                                                      
!!           isoverlap,op%wfd_overlap(iiorb,jorb)%nseg_c, op%wfd_overlap(iiorb,jorb)%nvctr_c, &
!!           op%wfd_overlap(iiorb,jorb)%keyglob(1,1), op%wfd_overlap(iiorb,jorb)%keyvglob(1))
!!        !Get fine part
!!        if(op%wfd_overlap(iiorb,jorb)%nseg_f > 0)then
!!        call get_overlap_from_descriptors_periodic(lzd%Llr(ilr)%wfd%nseg_c, lzd%Llr(jlr)%wfd%nseg_f,&
!!           lzd%Llr(ilr)%wfd%keyglob(1,1), lzd%Llr(jlr)%wfd%keyglob(1,1+lzd%Llr(jlr)%wfd%nseg_c), &                                                                                                      
!!           isoverlapf,op%wfd_overlap(iiorb,jorb)%nseg_f, op%wfd_overlap(iiorb,jorb)%nvctr_f, &
!!           op%wfd_overlap(iiorb,jorb)%keyglob(1,1+op%wfd_overlap(iiorb,jorb)%nseg_c), &
!!           op%wfd_overlap(iiorb,jorb)%keyvglob(1+op%wfd_overlap(iiorb,jorb)%nseg_c))
!!        else
!!           op%wfd_overlap(iiorb,jorb)%nvctr_f = 0
!!        end if
!!    end do
!!end do
!!
!!end subroutine determineOverlapDescriptorsVariable2

!!!subroutine setCommsOrthoVariable(iproc, nproc, orbs, orbsig, lzd, op, comon, tag)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(orbitals_data),intent(in):: orbs, orbsig
!!!type(local_zone_descriptors),intent(in):: lzd
!!!type(overlapParameters),intent(inout):: op
!!!type(p2pComms),intent(out):: comon
!!!integer,intent(inout):: tag
!!!
!!!! Local variables
!!!integer:: jproc, iorb, jorb, iiorb, jjorb, mpisource, mpidest, istsource, istdest, ncount, istat, iall, ijorb
!!!integer:: ilr, ilrold, jprocold, ildim, ierr
!!!integer,dimension(:),allocatable:: istsourceArr, istdestArr
!!!character(len=*),parameter:: subname='setCommsOrtho'
!!!logical,dimension(:),allocatable:: receivedOrbital
!!!
!!!allocate(istsourceArr(0:nproc-1), stat=istat)
!!!call memocc(istat, istsourceArr, 'istsourceArr',subname)
!!!allocate(istdestArr(0:nproc-1), stat=istat)
!!!call memocc(istat, istdestArr, 'istdestArr',subname)
!!!allocate(receivedOrbital(orbsig%norb), stat=istat)
!!!call memocc(istat, receivedOrbital, 'receivedOrbital', subname)
!!!
!!!istsourceArr=1
!!!istdestArr=1
!!!
!!!comon%nsendBuf=0
!!!comon%nrecvBuf=0
!!!
!!!
!!!op%indexInRecvBuf=0
!!!op%ndim_lphiovrlp=0
!!!
!!!iiorb=0
!!!jprocold=-1
!!!do jproc=0,nproc-1
!!!    receivedOrbital=.false.
!!!    ijorb=0
!!!    ilrold=-1
!!!    do iorb=1,orbs%norb_par(jproc,0)
!!!       iiorb=iiorb+1 
!!!       ilr=orbs%inWhichLocreg(iiorb)
!!!       ! Check whether process jproc has already received orbital jjorb.
!!!       !if(iproc==0) write(*,'(a,5i8)') 'jproc, iorb, iiorb, ilr, ilrold', jproc, iorb, iiorb, ilr, ilrold
!!!       if(ilr==ilrold) cycle
!!!       ildim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!!       do jorb=1,op%noverlaps(iiorb)
!!!           jjorb=op%overlaps(jorb,iiorb)
!!!           !write(*,'(a,7i8)') 'iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold', iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold
!!!           ijorb=ijorb+1
!!!           mpisource=orbsig%onWhichMPI(jjorb)
!!!           mpidest=jproc
!!!           istsource=istsourceArr(mpisource)
!!!           istdest=istdestArr(mpidest)
!!!           if(iproc==jproc) then
!!!               ncount=op%wfd_overlap(jorb,iiorb)%nvctr_c+7*op%wfd_overlap(jorb,iiorb)%nvctr_f
!!!           end if
!!!           call mpi_bcast(ncount, 1, mpi_integer, jproc, mpi_comm_world, ierr)
!!!           tag=tag+1
!!!           receivedOrbital(jjorb)=.true.
!!!           call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comon%comarr(1,ijorb,jproc))
!!!           !if(iproc==0) write(*,'(6(a,i0))') 'process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,', tag=',tag
!!!           if(iproc==mpisource) then
!!!               !write(*,'(5(a,i0))') 'adding ',ncount,' elements from orbital ',jjorb,' for orbital ',iiorb,' to nsendBuf, iproc=',iproc,', jproc=',jproc
!!!               comon%nsendBuf=comon%nsendBuf+ncount
!!!print *,'2: comon%nsendbuf',comon%nsendBuf
!!!           end if
!!!           if(iproc==mpidest) then
!!!               !write(*,'(3(a,i0))') 'process ',iproc,' will get orbital ',jjorb,' at position ',istdest
!!!               op%indexInRecvBuf(iorb,jjorb)=istdest
!!!               comon%nrecvBuf=comon%nrecvBuf+ncount
!!!               op%ndim_lphiovrlp=op%ndim_lphiovrlp+ildim
!!!           end if
!!!           istsourceArr(mpisource) = istsourceArr(mpisource) + ncount
!!!           istdestArr(mpidest) = istdestArr(mpidest) + ncount
!!!       end do
!!!       ilrold=ilr
!!!    end do
!!!    jprocold=jproc
!!!end do
!!!   
!!!iall = -product(shape(istsourceArr))*kind(istsourceArr)
!!!deallocate(istsourceArr, stat=istat)
!!!call memocc(istat, iall, 'istsourceArr',subname)
!!!iall = -product(shape(istdestArr))*kind(istdestArr)
!!!deallocate(istdestArr, stat=istat)
!!!call memocc(istat, iall, 'istdestArr',subname)
!!!iall = -product(shape(receivedOrbital))*kind(receivedOrbital)
!!!deallocate(receivedOrbital, stat=istat)
!!!call memocc(istat, iall, 'receivedOrbital',subname)
!!!
!!!
!!!end subroutine setCommsOrthoVariable

!!subroutine indicesForExpansionVariable(iproc, nproc, orbs, input, lzd, op, comon)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => indicesForExpansionVariable
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(orbitals_data),intent(in):: orbs
!!type(input_variables),intent(in):: input
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(inout):: op
!!type(p2pComms),intent(in):: comon
!!
!!! Local variables
!!integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, ierr
!!
!!
!!
!!ind=1
!!ilrold=-1
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=orbs%inWhichLocreg(iiorb)
!!    if(ilr==ilrold) cycle
!!    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        ! Starting index of orbital jjorb
!!        jst=op%indexInRecvBuf(iorb,jjorb)
!!        ldim=op%wfd_overlap(jorb,iiorb)%nvctr_c+7*op%wfd_overlap(jorb,iiorb)%nvctr_f
!!        call index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, &
!!             input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), op%indexExpand(jst:jst+ldim-1))
!!        ind=ind+gdim
!!    end do
!!    ilrold=ilr
!!end do
!!if(jst+ldim/=comon%nrecvBuf+1) then
!!    write(*,*) 'ERROR on process ',iproc,': jst+ldim/=comon%nrecvBuf+1',jst+ldim,comon%nrecvBuf+1
!!    stop
!!end if
!!
!!end subroutine indicesForExpansionVariable



!!subroutine indicesForExtractionVariable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, comon)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, sizePhi
!!type(orbitals_data),intent(in):: orbs, orbsig
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(inout):: op
!!type(p2pComms),intent(out):: comon
!!
!!! Local variables
!!integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
!!integer:: i, jj, j, jjjlr, iilr, jlr, jjorb2, iiorb, ijorb
!!
!!indovrlp=1
!!op%indexInSendBuf=0
!!
!!ilrold=-1
!!ijorb=0
!!do iorb=1,orbs%norb
!!    ilr=orbs%inWhichLocreg(iorb)
!!    iiproc=orbs%onWhichMPI(iorb)
!!    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
!!    do jorb=1,op%noverlaps(iorb)
!!        jjorb=op%overlaps(jorb,iorb)
!!        jjlr=orbsig%inWhichLocreg(jjorb)
!!        jjproc=orbsig%onWhichMPI(jjorb)
!!        if(iproc==jjproc) then
!!            ijorb=ijorb+1
!!            ! Get the correct descriptors
!!            ! Get an orbs-orbital in the same locreg as the orbsig-orbital jjorb.
!!            do j=1,orbs%norb
!!                jjjlr=orbs%inWhichLocreg(j)
!!                if(jjjlr==jjlr) then
!!                    korb=j
!!                    exit
!!                end if
!!            end do
!!            ! Get an orbsig-orbital in the same locreg as the orbs-orbital iorb.
!!            lorb=0
!!            do i=1,op%noverlaps(korb)
!!                iiorb=op%overlaps(i,korb)
!!                iilr=orbsig%inWhichLocreg(iiorb)
!!                if(iilr==ilr) then
!!                    lorb=i
!!                    exit
!!                end if
!!            end do
!!            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
!!            ldim=op%wfd_overlap(lorb,korb)%nvctr_c+7*op%wfd_overlap(lorb,korb)%nvctr_f
!!            ind=1
!!            do kkorb=orbs%isorb+1,jjorb-1
!!                klr=orbsig%inWhichLocreg(kkorb)
!!                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
!!            end do
!!            call index_of_psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), op%indexExtract(indovrlp))
!!            op%indexInSendBuf(jjorb-orbsig%isorb,iorb)=indovrlp
!!            indovrlp=indovrlp+op%wfd_overlap(lorb,korb)%nvctr_c+7*op%wfd_overlap(lorb,korb)%nvctr_f
!!        end if
!!    end do
!!    ilrold=ilr
!!    iiprocold=iiproc
!!end do
!!
!!if(indovrlp/=comon%nsendBuf+1) then
!!    write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
!!    stop
!!end if
!!
!!
!!
!!end subroutine indicesForExtractionVariable



!!!!subroutine extractOrbital2Variable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, phi, comon)
!!!!use module_base
!!!!use module_types
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc, sizePhi
!!!!type(orbitals_data),intent(in):: orbs, orbsig
!!!!type(local_zone_descriptors),intent(in):: lzd
!!!!type(overlapParameters),intent(inout):: op
!!!!real(8),dimension(sizePhi),intent(in):: phi
!!!!type(p2pComms),intent(out):: comon
!!!!
!!!!! Local variables
!!!!integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
!!!!integer:: i, indSource, j, jjjlr, iiorb, iilr, iseg, istart, iend, ncount, kseg, kstart, kend, kold, start, jst, ifine, isend
!!!!integer:: igrid
!!!!
!!!!indovrlp=1
!!!!op%indexInSendBuf=0
!!!!
!!!!ilrold=-1
!!!!iiprocold=-1
!!!!do iorb=1,orbs%norb
!!!!    ilr=orbs%inWhichLocreg(iorb)
!!!!    iiproc=orbs%onWhichMPI(iorb)
!!!!    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
!!!!    do jorb=1,op%noverlaps(iorb)
!!!!        jjorb=op%overlaps(jorb,iorb)
!!!!        jjlr=orbsig%inWhichLocreg(jjorb)
!!!!        jjproc=orbsig%onWhichMPI(jjorb)
!!!!        if(iproc==jjproc) then
!!!!            ! Get the correct descriptors
!!!!            ! Get an orbs-orbital in the same locreg as the orbsig-orbital jjorb.
!!!!            do j=1,orbs%norb
!!!!                jjjlr=orbs%inWhichLocreg(j)
!!!!                if(jjjlr==jjlr) then
!!!!                    korb=j
!!!!                    exit
!!!!                end if
!!!!            end do
!!!!            ! Get an orbsig-orbital in the same locreg as the orbs-orbital iorb.
!!!!            lorb=0
!!!!            do i=1,op%noverlaps(korb)
!!!!                iiorb=op%overlaps(i,korb)
!!!!                iilr=orbsig%inWhichLocreg(iiorb)
!!!!                if(iilr==ilr) then
!!!!                    lorb=i
!!!!                    exit
!!!!                end if
!!!!            end do
!!!!            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
!!!!            ldim=op%wfd_overlap(lorb,korb)%nvctr_c+7*op%wfd_overlap(lorb,korb)%nvctr_f
!!!!            ind=1
!!!!            do kkorb=orbsig%isorb+1,jjorb-1
!!!!                klr=orbsig%inWhichLocreg(kkorb)
!!!!                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
!!!!            end do
!!!!!!            do i=0,ldim-1
!!!!!!                indSource=ind+op%indexExtract(indovrlp+i)-1
!!!!!!                comon%sendBuf(indovrlp+i)=phi(indSource)
!!!!!!            end do
!!!!              !! THIS IS THE NEWEST VERSION (NOT OPTIMIZED, could probably store the keyv once and for all)
!!!!              jst=0
!!!!              kold = 1
!!!!              do iseg=1,op%wfd_overlap(lorb,korb)%nseg_c
!!!!                  istart=op%wfd_overlap(lorb,korb)%keyglob(1,iseg)
!!!!                  iend=op%wfd_overlap(lorb,korb)%keyglob(2,iseg)
!!!!                  ncount=iend-istart+1  !the overlap segment is always contained inside the Llrs segments, so ncount is ok
!!!!                  inner_loop: do kseg=kold,lzd%llr(jjlr)%wfd%nseg_c
!!!!                     kstart = lzd%llr(jjlr)%wfd%keyglob(1,kseg)
!!!!                     kend   = lzd%llr(jjlr)%wfd%keyglob(2,kseg)
!!!!                     if(kstart <= iend .and. kend >= istart) then  
!!!!                        kold = kseg
!!!!                        start = lzd%llr(jjlr)%wfd%keyvglob(kseg) + max(0,istart-kstart)
!!!!                        call dcopy(ncount, phi(ind+start-1), 1, comon%sendBuf(indovrlp+jst), 1)
!!!!                        jst=jst+ncount
!!!!                        exit inner_loop
!!!!                     end if
!!!!                  end do inner_loop
!!!!              end do
!!!!              if(jst .ne. op%wfd_overlap(lorb,korb)%nvctr_c) then
!!!!                 print *,'extractOrbital2V:jst = ',jst,'not equal to ldim = ',op%wfd_overlap(lorb,korb)%nvctr_c
!!!!                 stop
!!!!              end if
!!!!              !DO FINE GRID 
!!!!              jst=0
!!!!              kold = 1
!!!!              do iseg=1,op%wfd_overlap(lorb,korb)%nseg_f
!!!!                  istart=op%wfd_overlap(lorb,korb)%keyglob(1,iseg+op%wfd_overlap(lorb,korb)%nseg_c)
!!!!                  iend=op%wfd_overlap(lorb,korb)%keyglob(2,iseg+op%wfd_overlap(lorb,korb)%nseg_c)
!!!!                  ncount=7*(iend-istart+1)  !the overlap segment is always contained inside the Llrs segments, so ncount is ok
!!!!                  inner_loop2: do kseg=kold,lzd%llr(jjlr)%wfd%nseg_f
!!!!                     kstart = lzd%llr(jjlr)%wfd%keyglob(1,kseg+lzd%llr(jjlr)%wfd%nseg_c)
!!!!                     kend   = lzd%llr(jjlr)%wfd%keyglob(2,kseg+lzd%llr(jjlr)%wfd%nseg_c)
!!!!                     if(kstart <= iend .and. kend >= istart) then  
!!!!                        kold = kseg
!!!!                        start = lzd%llr(jjlr)%wfd%nvctr_c+(lzd%llr(jjlr)%wfd%keyvglob(kseg+lzd%llr(jjlr)%wfd%nseg_c) +&
!!!!                                max(0,istart-kstart)-1)*7
!!!!                        call dcopy(ncount,phi(ind+start),1,&
!!!!                                comon%sendBuf(indovrlp+jst+op%wfd_overlap(lorb,korb)%nvctr_c),1)
!!!!                        jst=jst+ncount
!!!!                        exit inner_loop2
!!!!                     end if
!!!!                  end do inner_loop2
!!!!              end do
!!!!              if(jst .ne. 7*op%wfd_overlap(lorb,korb)%nvctr_f) then
!!!!                 print *,'extractOrbital2V:jst = ',jst,'not equal to ldim = ',7*op%wfd_overlap(lorb,korb)%nvctr_f
!!!!                 stop
!!!!              end if
!!!!            op%indexInSendBuf(jjorb-orbsig%isorb,iorb)=indovrlp
!!!!            indovrlp=indovrlp+op%wfd_overlap(lorb,korb)%nvctr_c+7*op%wfd_overlap(lorb,korb)%nvctr_f
!!!!        end if
!!!!    end do
!!!!    ilrold=ilr
!!!!    iiprocold=iiproc
!!!!end do
!!!!
!!!!if(indovrlp/=comon%nsendBuf+1) then
!!!!    write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
!!!!    stop
!!!!end if
!!!!
!!!!end subroutine extractOrbital2Variable



!!subroutine expandOrbital2Variable(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(inc, nproc
!!type(orbitals_data),intent(in):: orbs
!!type(input_variables),intent(in):: input
!!type(local_zone_descriptors),intent(in):: lzd
!!type(overlapParameters),intent(in):: op
!!type(p2pComms),intent(in):: comon
!!real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!
!!! Local variables
!!integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest
!!
!!
!!lphiovrlp=0.d0
!!
!!ind=1
!!ilrold=-1
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=orbs%inWhichLocreg(iiorb)
!!    if(ilr==ilrold) cycle
!!    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        ! Starting index of orbital jjorb
!!        jst=op%indexInRecvBuf(iorb,jjorb)
!!        ldim=op%wfd_overlap(jorb,iiorb)%nvctr_c+7*op%wfd_overlap(jorb,iiorb)%nvctr_f
!!        do i=0,ldim-1
!!            indDest=ind+op%indexExpand(jst+i)-1
!!            lphiovrlp(indDest)=comon%recvBuf(jst+i)
!!        end do
!!        ind=ind+gdim
!!    end do
!!    ilrold=ilr
!!end do
!!
!!end subroutine expandOrbital2Variable





!!subroutine getOrbitals(iproc, nproc, comon)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(p2pCommsOrthonormality),intent(inout):: comon
!!
!!! Local variables
!!integer:: jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, nsends, nreceives, ierr
!!integer:: win, sizeOfDouble
!!real(8),dimension(:),allocatable:: ttarr
!!
!!
!!call mpi_type_size(mpi_double_precision, sizeOfDouble, ierr)
!!write(*,*) 'iproc, comon%nsendBuf, sizeOfDouble', iproc, comon%nsendBuf, sizeOfDouble
!!
!!nsends=0
!!nreceives=0
!!comon%communComplete=.false.
!!do jproc=0,nproc-1
!!    !write(*,'(3(a,i0))') 'iproc=',iproc,', jproc=',jproc,', comon%noverlaps(jproc)=', comon%noverlaps(jproc)
!!    do iorb=1,comon%noverlaps(jproc)
!!        mpisource=comon%comarr(1,iorb,jproc)
!!        istsource=comon%comarr(2,iorb,jproc)
!!        ncount=comon%comarr(3,iorb,jproc)
!!        mpidest=comon%comarr(4,iorb,jproc)
!!        istdest=comon%comarr(5,iorb,jproc)
!!        tag=comon%comarr(6,iorb,jproc)
!!        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
!!        if(mpisource/=mpidest) then
!!            ! The orbitals are on different processes, so we need a point to point communication.
!!            if(iproc==mpisource) then
!!            call mpi_win_create(comon%sendBuf(istsource), ncount*sizeOfDouble, sizeOfDouble, mpi_info_null, mpi_comm_world, win, ierr)
!!            call mpi_win_fence(0, win, ierr)
!!            !call mpi_get(comon%recvBuf(1), ncount, mpi_double_precision, mpisource, istsource, ncount, mpi_double_precision, win, ierr)
!!            !if(iproc==mpidest) then
!!            !    call mpi_get(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, istsource, ncount, mpi_double_precision, win, ierr)
!!            !    !call mpi_get(comon%recvBuf(1), 0, mpi_double_precision, mpisource, istsource, 0, mpi_double_precision, win, ierr)
!!            !else
!!            !    call mpi_get(ttarr(1), ncount, mpi_double_precision, mpisource, istsource, ncount, mpi_double_precision, win, ierr)
!!            !    !call mpi_get(comon%recvBuf(1), 0, mpi_double_precision, mpisource, istsource, 0, mpi_double_precision, win, ierr)
!!            !end if
!!            call mpi_win_fence(0, win, ierr)
!!            call mpi_win_free(win, ierr)
!!            if(iproc==mpisource) then
!!                !write(*,'(6(a,i0))') 'overlap: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
!!                !!call mpi_isend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!!                !!     mpi_comm_world, comon%comarr(7,iorb,jproc), ierr)
!!                !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
!!                !!comon%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
!!                !!nsends=nsends+1
!!            else if(iproc==mpidest) then
!!                !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!                !!call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!!                !!     mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
!!                !!comon%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
!!                !!nreceives=nreceives+1
!!            else
!!                !comon%comarr(7,iorb,jproc)=mpi_request_null
!!                !comon%comarr(8,iorb,jproc)=mpi_request_null
!!            end if
!!        else
!!            ! The orbitals are on the same process, so simply copy them.
!!            if(iproc==mpisource) then
!!                call dcopy(ncount, comon%sendBuf(istsource), 1, comon%recvBuf(istdest), 1)
!!                !write(*,'(6(a,i0))') 'overlap: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
!!                comon%comarr(7,iorb,jproc)=mpi_request_null
!!                comon%comarr(8,iorb,jproc)=mpi_request_null
!!                nsends=nsends+1
!!                nreceives=nreceives+1
!!                comon%communComplete(iorb,iproc)=.true.
!!            else
!!                comon%comarr(7,iorb,jproc)=mpi_request_null
!!                comon%comarr(8,iorb,jproc)=mpi_request_null
!!            end if
!!
!!        end if
!!    end do
!!end do
!!
!!
!!
!!end subroutine getOrbitals

!!subroutine collectAndCalculateOverlap(iproc, nproc, comon, mad, op, orbs, input, lzd, &
!!           nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, lphiovrlp,timecommunp2p, &
!!           timecommuncoll, timeoverlap, timeexpand, timecompress)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => collectAndCalculateOverlap
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
!!type(p2pComms),intent(inout):: comon
!!type(matrixDescriptors),intent(in):: mad
!!type(overlapParameters),intent(in):: op
!!type(orbitals_data),intent(in):: orbs
!!type(input_variables),intent(in):: input
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),dimension(nsendbuf),intent(in):: sendbuf
!!real(8),dimension(nrecvbuf),intent(inout):: recvbuf
!!real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
!!real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!real(8),intent(inout):: timecommunp2p, timecommuncoll, timeoverlap, timeexpand, timecompress
!!
!!! Local variables
!!integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
!!real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
!!integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
!!real(8):: ddot, t1, t2
!!character(len=*),parameter:: subname='collectAndCalculateOverlap'
!!logical,dimension(:),allocatable:: done
!!logical:: received
!!
!!allocate(indcomplete(comon%nrecv), stat=istat)
!!call memocc(istat, indcomplete, 'indcomplete', subname)
!!allocate(done(comon%nrecv), stat=istat)
!!call memocc(istat, done, 'done', subname)
!!done=.false.
!!
!!allocate(indexarray(comon%nrecv), stat=istat)
!!call memocc(istat, indexarray, 'indexarray', subname)
!!do i=1,comon%nrecv
!!    indexarray(i)=i
!!end do
!!
!!
!!! Now the sends
!!if (nproc >1) then
!!   t1=mpi_wtime()
!!   nsend=0
!!   waitLoopSend: do
!!      !!call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!      !!nsend=nsend+ncomplete
!!      !!if(nsend==comon%nsend) exit waitLoopSend
!!      call mpi_waitany(comon%nsend-nsend, comon%requests(1,1), ind, mpi_status_ignore, ierr)
!!      nsend=nsend+1
!!      do i=ind,comon%nsend-nsend
!!         comon%requests(i,1)=comon%requests(i+1,1)
!!      end do
!!      if(nsend==comon%nsend) exit waitLoopSend
!!   end do waitLoopSend
!!   t2=mpi_wtime()
!!   timecommunp2p=timecommunp2p+t2-t1
!!end if
!!
!!ovrlp=0.d0
!!lphiovrlp=0.d0
!!nrecv=0
!!waitLoopRecv: do
!!   if (nproc > 1) then
!!      t1=mpi_wtime()
!!      call mpi_waitany(comon%nrecv-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
!!      !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
!!      !ind=1
!!      t2=mpi_wtime()
!!      timecommunp2p=timecommunp2p+t2-t1
!!   end if
!!   ncomplete=1
!!   received=.true.
!!   if(received) then
!!      do i=1,ncomplete
!!         ! Calculate overlap matrix
!!         !jorb=indcomplete(i)
!!         jorb=indexarray(ind)
!!         !!write(*,'(2(a,i0))') 'process ',iproc,' has received message ',jorb
!!         if(.not.done(jorb)) then
!!            orbsource=comon%comarr(9,jorb,iproc)
!!            do iorb=1,orbs%norbp
!!               orbdest=orbs%isorb+iorb
!!               call getStartingIndicesGlobal(orbdest, orbsource, op, orbs, ist, jst, ncount)
!!               !jst=comon%comarr(5,jorb,iproc)
!!               !ncount=comon%comarr(3,jorb,iproc)
!!               !write(*,'(a,i4,2i5,3x,2i8,4x,2i8)') 'iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc) ; jst, comon%comarr(5,jorb,iproc)', iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc), jst, comon%comarr(5,jorb,iproc)
!!               t1=mpi_wtime()
!!               ovrlp(orbsource,orbdest)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
!!               t2=mpi_wtime()
!!               timeoverlap=timeoverlap+t2-t1
!!               t1=mpi_wtime()
!!               call expandoneorbital2(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
!!                    orbs%inwhichlocreg, lzd, op, nrecvbuf, recvbuf, lphiovrlp)
!!               t2=mpi_wtime()
!!               timeexpand=timeexpand+t2-t1
!!            end do
!!            done(jorb)=.true.
!!         end if
!!      end do
!!      nrecv=nrecv+ncomplete
!!      !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
!!      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
!!      do i=ind,comon%nrecv-nrecv
!!         comon%requests(i,2)=comon%requests(i+1,2)
!!         indexarray(i)=indexarray(i+1)
!!      end do
!!      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
!!      if(nrecv==comon%nrecv) exit waitLoopRecv
!!   end if
!!end do waitLoopRecv
!!
!!!write(*,'(a,i0,a)') 'iproc=',iproc,' is here'
!!
!!iall=-product(shape(done))*kind(done)
!!deallocate(done, stat=istat)
!!call memocc(istat, iall, 'done', subname)
!!
!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!deallocate(indcomplete, stat=istat)
!!call memocc(istat, iall, 'indcomplete', subname)
!!
!!iall=-product(shape(indexarray))*kind(indexarray)
!!deallocate(indexarray, stat=istat)
!!call memocc(istat, iall, 'indexarray', subname)
!!
!!!!allocate(indcomplete(comon%nsend), stat=istat)
!!!!call memocc(istat, indcomplete, 'indcomplete', subname)
!!
!!!!do iorb=1,orbs%norb
!!!!  do jorb=1,orbs%norb
!!!!    write(400+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!!  end do
!!!!end do
!!
!!!!! Now the sends
!!!!nsend=0
!!!!waitLoopSend: do
!!!!    call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!!!    nsend=nsend+ncomplete
!!!!    if(nsend==comon%nsend) exit waitLoopSend
!!!!end do waitLoopSend
!!
!!!write(*,*) 'here: iproc', iproc
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call mpi_finalize(ierr)
!!!stop
!!
!!!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!!!deallocate(indcomplete, stat=istat)
!!!!call memocc(istat, iall, 'indcomplete', subname)
!!
!!!do iorb=1,orbs%norb
!!!  do jorb=1,orbs%norb
!!!    write(90+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!  end do
!!!end do
!!
!!
!!! Communicate the matrix
!!allocate(ovrlpCompressed(mad%nvctr), stat=istat)
!!call memocc(istat, ovrlpCompressed, 'ovrlpCompressed', subname)
!!allocate(sendcounts(0:nproc-1), stat=istat)
!!call memocc(istat, sendcounts, 'sendcounts', subname)
!!allocate(displs(0:nproc-1), stat=istat)
!!call memocc(istat, displs, 'displs', subname)
!!
!!t1=mpi_wtime()
!!call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed, sendcounts, displs)
!!t2=mpi_wtime()
!!timecompress=timecompress+t2-t1     
!!
!!allocate(ovrlpCompressed2(mad%nvctr), stat=istat)
!!call memocc(istat, ovrlpCompressed2, 'ovrlpCompressed2', subname)
!!
!!t1=mpi_wtime()
!!if (nproc >1) then
!!   call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
!!        sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
!!else
!!   call vcopy(sendcounts(iproc),ovrlpCompressed(displs(iproc)+1),1,ovrlpCompressed2(1+displs(iproc)),1)
!!end if
!!
!!t2=mpi_wtime()
!!timecommuncoll=timecommuncoll+t2-t1     
!!
!!t1=mpi_wtime()
!!call uncompress_matrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)
!!t2=mpi_wtime()
!!timecompress=timecompress+t2-t1     
!!
!!!do iorb=1,orbs%norb
!!!  do jorb=1,orbs%norb
!!!    write(100+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!  end do
!!!end do
!!
!!iall=-product(shape(ovrlpCompressed))*kind(ovrlpCompressed)
!!deallocate(ovrlpCompressed, stat=istat)
!!call memocc(istat, iall, 'ovrlpCompressed', subname)
!!iall=-product(shape(ovrlpCompressed2))*kind(ovrlpCompressed2)
!!deallocate(ovrlpCompressed2, stat=istat)
!!call memocc(istat, iall, 'ovrlpCompressed2', subname)
!!iall=-product(shape(sendcounts))*kind(sendcounts)
!!deallocate(sendcounts, stat=istat)
!!call memocc(istat, iall, 'sendcounts', subname)
!!iall=-product(shape(displs))*kind(displs)
!!deallocate(displs, stat=istat)
!!call memocc(istat, iall, 'displs', subname)
!!
!!
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call mpi_finalize(ierr)
!!
!!end subroutine collectAndCalculateOverlap

!!!subroutine collectnew2(iproc, nproc, comon, mad, op, orbs, input, lzd, &
!!!   nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, timecommunp2p, timecommuncoll, timecompress)
!!!use module_base
!!!use module_types
!!!use module_interfaces, exceptThisOne => collectnew
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
!!!type(p2pCommsOrthonormality),intent(inout):: comon
!!!type(matrixDescriptors),intent(in):: mad
!!!type(overlapParameters),intent(in):: op
!!!type(orbitals_data),intent(in):: orbs
!!!type(input_variables),intent(in):: input
!!!type(local_zone_descriptors),intent(in):: lzd
!!!real(8),dimension(nsendbuf),intent(in):: sendbuf
!!!real(8),dimension(nrecvbuf),intent(inout):: recvbuf
!!!real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
!!!real(8),intent(inout):: timecommunp2p, timecommuncoll, timecompress
!!!
!!!! Local variables
!!!integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
!!!real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
!!!integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
!!!real(8):: ddot, t1, t2
!!!character(len=*),parameter:: subname='collectAndCalculateOverlap'
!!!logical,dimension(:),allocatable:: done
!!!logical:: received
!!!
!!!
!!!allocate(indexarray(comon%nrecvtemp), stat=istat)
!!!call memocc(istat, indexarray, 'indexarray', subname)
!!!do i=1,comon%nrecvtemp
!!!   indexarray(i)=i
!!!end do
!!!
!!!if (nproc > 1) then
!!!   ! Wait for the sends to complete.
!!!   t1=mpi_wtime()
!!!   nsend=0
!!!   if(comon%nsendtemp/=0) then
!!!      waitLoopSend: do
!!!         !!call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!!         !!nsend=nsend+ncomplete
!!!         !!if(nsend==comon%nsend) exit waitLoopSend
!!!         call mpi_waitany(comon%nsendtemp-nsend, comon%requests(1,1), ind, mpi_status_ignore, ierr)
!!!         nsend=nsend+1
!!!         do i=ind,comon%nsendtemp-nsend
!!!            comon%requests(i,1)=comon%requests(i+1,1)
!!!         end do
!!!         if(nsend==comon%nsendtemp) exit waitLoopSend
!!!      end do waitLoopSend
!!!   end if
!!!   t2=mpi_wtime()
!!!   timecommunp2p=timecommunp2p+t2-t1
!!!
!!!
!!!
!!!   ovrlp=0.d0
!!!   nrecv=0
!!!   if(comon%nrecvtemp/=0) then
!!!      waitLoopRecv: do
!!!         t1=mpi_wtime()
!!!         call mpi_waitany(comon%nrecvtemp-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
!!!         !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
!!!         !ind=1
!!!         t2=mpi_wtime()
!!!         timecommunp2p=timecommunp2p+t2-t1
!!!         ncomplete=1
!!!         received=.true.
!!!         if(received) then
!!!            nrecv=nrecv+ncomplete
!!!            !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
!!!            !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
!!!            do i=ind,comon%nrecvtemp-nrecv
!!!               comon%requests(i,2)=comon%requests(i+1,2)
!!!               indexarray(i)=indexarray(i+1)
!!!            end do
!!!            !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
!!!            if(nrecv==comon%nrecvtemp) exit waitLoopRecv
!!!         end if
!!!      end do waitLoopRecv
!!!   end if
!!!end if
!!!
!!!iall=-product(shape(indexarray))*kind(indexarray)
!!!deallocate(indexarray, stat=istat)
!!!call memocc(istat, iall, 'indexarray', subname)
!!!
!!!
!!!end subroutine collectnew2








!!subroutine collectAndCalculateOverlap2(iproc, nproc, comon, mad, op, orbs, input, lzd, &
!!   nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, timecommunp2p, timecommuncoll, timeoverlap, timecompress)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => collectAndCalculateOverlap2
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
!!type(p2pComms),intent(inout):: comon
!!type(matrixDescriptors),intent(in):: mad
!!type(overlapParameters),intent(in):: op
!!type(orbitals_data),intent(in):: orbs
!!type(input_variables),intent(in):: input
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),dimension(nsendbuf),intent(in):: sendbuf
!!real(8),dimension(nrecvbuf),intent(inout):: recvbuf
!!real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
!!real(8),intent(inout):: timecommunp2p, timecommuncoll, timeoverlap, timecompress
!!
!!! Local variables
!!integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
!!real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
!!integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
!!real(8):: ddot, t1, t2
!!character(len=*),parameter:: subname='collectAndCalculateOverlap'
!!logical,dimension(:),allocatable:: done
!!logical:: received
!!
!!allocate(indcomplete(comon%nrecv), stat=istat)
!!call memocc(istat, indcomplete, 'indcomplete', subname)
!!allocate(done(comon%nrecv), stat=istat)
!!call memocc(istat, done, 'done', subname)
!!done=.false.
!!
!!allocate(indexarray(comon%nrecv), stat=istat)
!!call memocc(istat, indexarray, 'indexarray', subname)
!!do i=1,comon%nrecv
!!   indexarray(i)=i
!!end do
!!
!!
!!
!!ovrlp=0.d0
!!nrecv=0
!!waitLoopRecv: do
!!   if (nproc > 1) then
!!      t1=mpi_wtime()
!!      call mpi_waitany(comon%nrecv-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
!!      !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
!!      !ind=1
!!      t2=mpi_wtime()
!!      timecommunp2p=timecommunp2p+t2-t1
!!   end if
!!   ncomplete=1
!!   received=.true.
!!   if(received) then
!!      do i=1,ncomplete
!!         ! Calculate overlap matrix
!!         !jorb=indcomplete(i)
!!         jorb=indexarray(ind)
!!         !!write(*,'(2(a,i0))') 'process ',iproc,' has received message ',jorb
!!         if(.not.done(jorb)) then
!!            orbsource=comon%comarr(9,jorb,iproc)
!!            do iorb=1,orbs%norbp
!!               orbdest=orbs%isorb+iorb
!!               call getStartingIndicesGlobal(orbdest, orbsource, op, orbs, ist, jst, ncount)
!!               !jst=comon%comarr(5,jorb,iproc)
!!               !ncount=comon%comarr(3,jorb,iproc)
!!               !write(*,'(a,i4,2i5,3x,2i8,4x,2i8)') 'iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc) ; jst, comon%comarr(5,jorb,iproc)', iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc), jst, comon%comarr(5,jorb,iproc)
!!               t1=mpi_wtime()
!!               ovrlp(orbsource,orbdest)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
!!               t2=mpi_wtime()
!!               timeoverlap=timeoverlap+t2-t1
!!               !!t1=mpi_wtime()
!!               !!call expandoneorbital2(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
!!               !!     orbs%inwhichlocreg, lzd, op, nrecvbuf, recvbuf, lphiovrlp)
!!               !!t2=mpi_wtime()
!!               !!timeexpand=timeexpand+t2-t1
!!            end do
!!            done(jorb)=.true.
!!         end if
!!      end do
!!      nrecv=nrecv+ncomplete
!!      !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
!!      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
!!      do i=ind,comon%nrecv-nrecv
!!         comon%requests(i,2)=comon%requests(i+1,2)
!!         indexarray(i)=indexarray(i+1)
!!      end do
!!      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
!!      if(nrecv==comon%nrecv) exit waitLoopRecv
!!   end if
!!end do waitLoopRecv
!!
!!!write(*,'(a,i0,a)') 'iproc=',iproc,' is here'
!!
!!iall=-product(shape(done))*kind(done)
!!deallocate(done, stat=istat)
!!call memocc(istat, iall, 'done', subname)
!!
!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!deallocate(indcomplete, stat=istat)
!!call memocc(istat, iall, 'indcomplete', subname)
!!
!!iall=-product(shape(indexarray))*kind(indexarray)
!!deallocate(indexarray, stat=istat)
!!call memocc(istat, iall, 'indexarray', subname)
!!
!!!!allocate(indcomplete(comon%nsend), stat=istat)
!!!!call memocc(istat, indcomplete, 'indcomplete', subname)
!!
!!!!do iorb=1,orbs%norb
!!!!  do jorb=1,orbs%norb
!!!!    write(400+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!!  end do
!!!!end do
!!
!!!!! Now the sends
!!!!nsend=0
!!!!waitLoopSend: do
!!!!    call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!!!    nsend=nsend+ncomplete
!!!!    if(nsend==comon%nsend) exit waitLoopSend
!!!!end do waitLoopSend
!!
!!!write(*,*) 'here: iproc', iproc
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call mpi_finalize(ierr)
!!!stop
!!
!!!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!!!deallocate(indcomplete, stat=istat)
!!!!call memocc(istat, iall, 'indcomplete', subname)
!!
!!!do iorb=1,orbs%norb
!!!  do jorb=1,orbs%norb
!!!    write(90+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!  end do
!!!end do
!!
!!
!!! Communicate the matrix
!!allocate(ovrlpCompressed(mad%nvctr), stat=istat)
!!call memocc(istat, ovrlpCompressed, 'ovrlpCompressed', subname)
!!allocate(sendcounts(0:nproc-1), stat=istat)
!!call memocc(istat, sendcounts, 'sendcounts', subname)
!!allocate(displs(0:nproc-1), stat=istat)
!!call memocc(istat, displs, 'displs', subname)
!!
!!t1=mpi_wtime()
!!call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed, sendcounts, displs)
!!t2=mpi_wtime()
!!timecompress=timecompress+t2-t1     
!!
!!allocate(ovrlpCompressed2(mad%nvctr), stat=istat)
!!call memocc(istat, ovrlpCompressed2, 'ovrlpCompressed2', subname)
!!
!!t1=mpi_wtime()
!!if (nproc > 1) then
!!   call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
!!        sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
!!else
!!   call vcopy(sendcounts(iproc),ovrlpCompressed(displs(iproc)+1),1,ovrlpCompressed2(1+displs(iproc)),1)
!!end if
!!t2=mpi_wtime()
!!timecommuncoll=timecommuncoll+t2-t1     
!!
!!t1=mpi_wtime()
!!call uncompress_matrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)
!!t2=mpi_wtime()
!!timecompress=timecompress+t2-t1     
!!
!!!do iorb=1,orbs%norb
!!!  do jorb=1,orbs%norb
!!!    write(100+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!  end do
!!!end do
!!
!!iall=-product(shape(ovrlpCompressed))*kind(ovrlpCompressed)
!!deallocate(ovrlpCompressed, stat=istat)
!!call memocc(istat, iall, 'ovrlpCompressed', subname)
!!iall=-product(shape(ovrlpCompressed2))*kind(ovrlpCompressed2)
!!deallocate(ovrlpCompressed2, stat=istat)
!!call memocc(istat, iall, 'ovrlpCompressed2', subname)
!!iall=-product(shape(sendcounts))*kind(sendcounts)
!!deallocate(sendcounts, stat=istat)
!!call memocc(istat, iall, 'sendcounts', subname)
!!iall=-product(shape(displs))*kind(displs)
!!deallocate(displs, stat=istat)
!!call memocc(istat, iall, 'displs', subname)
!!
!!
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call mpi_finalize(ierr)
!!
!!end subroutine collectAndCalculateOverlap2

!subroutine extractOneOrbital(iproc, nproc, orbsource, orbdest, orbs, lzd, op, phi, nsendBuf, sendBuf)
!use module_base
!use module_types
!implicit none
!
!! Calling arguments
!integer,intent(in):: iproc, nproc, orbsource, orbdest, nsendbuf
!type(orbitals_data),intent(in):: orbs
!type(local_zone_descriptors),intent(in):: lzd
!type(overlapParameters),intent(inout):: op
!real(8),dimension(lzd%lpsidimtot),intent(in):: phi
!real(8),dimension(nsendBuf),intent(out):: sendBuf
!
!
!! Local variables
!integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
!integer:: i, indSource
!
!indovrlp=1
!op%indexInSendBuf=0
!
!ilrold=-1
!iiprocold=-1
!do iorb=1,orbs%norb
!    ilr=orbs%inWhichLocreg(iorb)
!    iiproc=orbs%onWhichMPI(iorb)
!    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
!    do jorb=1,op%noverlaps(iorb)
!        jjorb=op%overlaps(jorb,iorb)
!        jjlr=orbs%inWhichLocreg(jjorb)
!        jjproc=orbs%onWhichMPI(jjorb)
!        if(iproc==jjproc) then
!            ! Get the correct descriptors
!            korb=jjorb-orbs%isorb
!            !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
!            do i=1,op%noverlaps(jjorb)
!                !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
!                if(op%overlaps(i,jjorb)==iorb) then
!                    lorb=i
!                    exit
!                end if
!            end do
!            !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
!            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
!            ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
!            ind=1
!            do kkorb=orbs%isorb+1,jjorb-1
!                klr=orbs%inWhichLocreg(kkorb)
!                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
!            end do
!            !write(*,'(5(a,i0))') 'process ',iproc,' adds ',op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f,' elements at position ',indovrlp,' from orbital ',jjorb,' for orbital ', iorb
!            !call psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), phi(ind), comon%sendBuf(indovrlp))
!            if(jjorb==orbsource .and. iorb==orbdest) then
!                do i=0,ldim-1
!                    indSource=ind+op%indexExtract(indovrlp+i)-1
!                    sendBuf(indovrlp+i)=phi(indSource)
!                end do
!            end if
!            op%indexInSendBuf(jjorb-orbs%isorb,iorb)=indovrlp
!            indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
!        end if
!    end do
!    ilrold=ilr
!    iiprocold=iiproc
!end do
!
!if(indovrlp/=nsendBuf+1) then
!    write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=nsendBuf+1', indovrlp, nsendBuf+1
!    stop
!end if
!
!end subroutine extractOneOrbital




!subroutine overlapMatrixCubic(iproc, nproc, gorbs, orbs, comms, lzd, input, methTransformOverlap, blocksize_dsyev, &
!           blocksize_pdgemm, mad, lphi)
!use module_base
!use module_types
!use module_interfaces
!implicit none
!
!! Calling arguments
!integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_dsyev, blocksize_pdgemm
!type(orbitals_data),intent(in):: orbs, gorbs
!type(comms_cubic),intent(in):: comms
!type(local_zone_descriptors),intent(in):: lzd
!type(input_variables),intent(in):: input
!type(matrixDescriptors),intent(in):: mad
!real(8),dimension(orbs%npsidim_orbs), intent(inout):: lphi
!
!! Local variables
!integer:: ind1, ind2, iorb, ilr, ldim, gdim, istat, nvctrp, ierr, iall
!real(8),dimension(:),allocatable:: phi, temparr
!real(8),dimension(:),pointer:: phiWork
!real(8),dimension(:,:),allocatable:: ovrlp
!character(len=*),parameter:: subname='overlapMatrixCubic'
!
!allocate(phi(gorbs%npsidim), stat=istat)
!call memocc(istat, phi, 'phi',subname)
!allocate(phiWork(gorbs%npsidim), stat=istat)
!call memocc(istat, phiWork, 'phiWork',subname)
!allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, ovrlp, 'ovrlp', subname)
!
!
!  ind1=1
!  ind2=1
!  phi=0.d0
!  do iorb=1,orbs%norbp
!      ilr = orbs%inWhichLocregp(iorb)
!      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
!      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
!      call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norb, orbs%nspinor, input%nspin, lzd%Glr,&
!           lzd%Llr(ilr), lphi(ind2), phi(ind1))
!      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
!      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
!  end do
!  call transpose_v(iproc, nproc, orbs, lzd%glr%wfd%Glr%wfd, comms, phi, work=phiWork)
!
!
!  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!  call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi, nvctrp, phi, nvctrp, 0.d0, ovrlp, orbs%norb)
!
!  call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
!
!  !!do ind1=1,orbs%norb
!  !!    do ind2=1,orbs%norb
!  !!        write(600+iproc,*) ind1, ind2, ovrlp(ind2,ind1)
!  !!    end do
!  !!end do
!  call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, methTransformOverlap, blocksize_dsyev, &
!        blocksize_pdgemm, orbs%norb, mad, ovrlp)
!
!  allocate(temparr(gorbs%npsidim), stat=istat)
!  call memocc(istat, temparr, 'temparr', subname)
!  call dgemm('n', 'n', nvctrp, orbs%norb, orbs%norb, 1.d0, phi, &
!       nvctrp, ovrlp, orbs%norb, 0.d0, temparr, nvctrp)
!
!  call untranspose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, temparr, work=phiWork)
!
!  ind1=1
!  ind2=1
!  do iorb=1,orbs%norbp
!      ilr = orbs%inWhichLocregp(iorb)
!      ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!      gdim=lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
!      call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzd%glr, temparr(ind1), lphi(ind2))
!      ind1=ind1+lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
!      ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!  end do
!
!
!
!  iall=-product(shape(phi))*kind(phi)
!  deallocate(phi, stat=istat)
!  call memocc(istat, iall, 'phi', subname)
!
!  iall=-product(shape(phiWork))*kind(phiWork)
!  deallocate(phiWork, stat=istat)
!  call memocc(istat, iall, 'phiWork', subname)
!
!  iall=-product(shape(temparr))*kind(temparr)
!  deallocate(temparr, stat=istat)
!  call memocc(istat, iall, 'temparr', subname)
!
!  iall=-product(shape(ovrlp))*kind(ovrlp)
!  deallocate(ovrlp, stat=istat)
!  call memocc(istat, iall, 'ovrlp', subname)
!
!
!
!end subroutine overlapMatrixCubic





!!subroutine applyOrthoconstraintNonorthogonalCubic(iproc, nproc, methTransformOverlap,&
!!     blocksize_pdgemm, orbs, gorbs, comms, lzd, input, &
!!     op, ovrlp, mad, lphi, lhphi, trH)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => applyOrthoconstraintNonorthogonalCubic
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm
!!type(orbitals_data),intent(in):: orbs, gorbs
!!type(comms_cubic),intent(in):: comms
!!type(local_zone_descriptors),intent(in):: lzd
!!type(input_variables),intent(in):: input
!!type(overlapParameters),intent(in):: op
!!real(8),dimension(orbs%norb,orbs%norb),intent(inout):: ovrlp
!!type(matrixDescriptors),intent(in):: mad
!!real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lphi, lhphi
!!real(8),intent(out):: trH
!!
!!! Local variables
!!integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, info, i, istat, iall, ierr, iseg, nvctrp, ldim, gdim, ind1, ind2
!!real(8):: tt, t1, t2, time_dsymm, time_daxpy
!!real(8),dimension(:,:),allocatable:: ovrlp2
!!real(8),dimension(:,:),allocatable:: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
!!character(len=*),parameter:: subname='applyOrthoconstraintNonorthogonalCubic'
!!logical:: val, segment
!!real(8),dimension(:,:),allocatable:: lagmat
!!real(8),dimension(:),allocatable:: phi, hphi
!!real(8),dimension(:),pointer:: phiWork
!!
!!allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, lagmat, 'lagmat', subname)
!!
!!allocate(phi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)), stat=istat)
!!call memocc(istat, phi, 'phi', subname)
!!allocate(hphi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)), stat=istat)
!!call memocc(istat, hphi, 'hphi', subname)
!!allocate(phiWork(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)), stat=istat)
!!call memocc(istat, phiWork, 'phiWork', subname)
!!
!!! Calculate the overlap matrix
!!  ind1=1
!!  ind2=1
!!  phi=0.d0
!!  do iorb=1,orbs%norbp
!!      !ilr = orbs%inWhichLocregp(iorb)
!!      iiorb=orbs%isorb+iorb
!!      ilr = orbs%inWhichLocreg(iiorb)
!!      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
!!      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
!!      call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norb, orbs%nspinor, input%nspin, lzd%Glr,&
!!           lzd%Llr(ilr), lphi(ind2), phi(ind1))
!!      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
!!      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
!!  end do
!!  call transpose_v(iproc, nproc, orbs, lzd%glr%wfd%Glr%wfd, comms, phi, work=phiWork)
!!
!!
!!  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!  call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi, nvctrp, phi, nvctrp, 0.d0, ovrlp, orbs%norb)
!!
!!  call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
!!
!!
!!  ! Now calculate the Lagrange multiplier matrix
!!  ind1=1
!!  ind2=1
!!  hphi=0.d0
!!  do iorb=1,orbs%norbp
!!      !ilr = orbs%inWhichLocregp(iorb)
!!      iiorb=orbs%isorb+iorb
!!      ilr = orbs%inWhichLocreg(iiorb)
!!      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
!!      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
!!      call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norb, orbs%nspinor, input%nspin, lzd%Glr,&
!!           lzd%Llr(ilr), lhphi(ind2), hphi(ind1))
!!      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
!!      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
!!  end do
!!  call transpose_v(iproc, nproc, orbs, lzd%glr%wfd%Glr%wfd, comms, hphi, work=phiWork)
!!
!!  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!  call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi, nvctrp, hphi, nvctrp, 0.d0, lagmat, orbs%norb)
!!
!!  call mpiallred(lagmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
!!  trH=0.d0
!!  do iorb=1,orbs%norb
!!      trH=trH+lagmat(iorb,iorb)
!!  end do
!!
!!
!!
!!
!!
!!allocate(ovrlp_minus_one_lagmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
!!allocate(ovrlp_minus_one_lagmat_trans(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
!!allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, ovrlp2, 'ovrlp2', subname)
!!
!!call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)
!!
!!! Invert the overlap matrix
!!call mpi_barrier(mpi_comm_world, ierr)
!!call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, orbs%norb, mad, orbs, ovrlp2)
!!
!!
!!! Multiply the Lagrange multiplier matrix with S^-1/2.
!!! First fill the upper triangle.
!!do iorb=1,orbs%norb
!!    do jorb=1,iorb-1
!!        ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
!!    end do
!!end do
!!call cpu_time(t1)
!!if(blocksize_pdgemm<0) then
!!    !call dsymm('l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
!!    !     0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
!!    ovrlp_minus_one_lagmat=0.d0
!!    call dgemm_compressed2(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
!!         mad%keygmatmul, ovrlp2, lagmat, ovrlp_minus_one_lagmat)
!!    ! Transpose lagmat
!!    do iorb=1,orbs%norb
!!        do jorb=iorb+1,orbs%norb
!!            tt=lagmat(jorb,iorb)
!!            lagmat(jorb,iorb)=lagmat(iorb,jorb)
!!            lagmat(iorb,jorb)=tt
!!        end do
!!    end do
!!    !call dsymm('l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
!!    !     0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
!!    ovrlp_minus_one_lagmat_trans=0.d0
!!    call dgemm_compressed2(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
!!         mad%keygmatmul, ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
!!else
!!    call dsymm_parallel(iproc, nproc, blocksize_pdgemm, mpi_comm_world, 'l', 'l', orbs%norb, orbs%norb, 1.d0, &
!!         ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, 0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
!!    ! Transpose lagmat
!!    do iorb=1,orbs%norb
!!        do jorb=iorb+1,orbs%norb
!!            tt=lagmat(jorb,iorb)
!!            lagmat(jorb,iorb)=lagmat(iorb,jorb)
!!            lagmat(iorb,jorb)=tt
!!        end do
!!    end do
!!    call dsymm_parallel(iproc, nproc, blocksize_pdgemm, mpi_comm_world, 'l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), &
!!         orbs%norb, lagmat(1,1), orbs%norb, &
!!         0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
!!end if
!!call cpu_time(t2)
!!time_dsymm=t2-t1
!!
!!
!!call dgemm('n', 'n', nvctrp, orbs%norb, orbs%norb, -.5_wp, phi, nvctrp, lagmat, orbs%norb, 1.0_wp, hphi, nvctrp)
!!call dgemm('n', 't', nvctrp, orbs%norb, orbs%norb, -.5_wp, phi, nvctrp, lagmat, orbs%norb, 1.0_wp, hphi, nvctrp)
!!
!!  call untranspose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, hphi, work=phiWork)
!!
!!  ind1=1
!!  ind2=1
!!  do iorb=1,orbs%norbp
!!      !ilr = orbs%inWhichLocregp(iorb)
!!      iiorb=orbs%isorb+iorb
!!      ilr = orbs%inWhichLocreg(iiorb)
!!      ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!      gdim=lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
!!      call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzd%glr, hphi(ind1), lhphi(ind2))
!!      ind1=ind1+lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
!!      ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!  end do
!!
!!
!!call cpu_time(t1)
!!!!ist=1
!!!!jst=1
!!!!ilrold=-1
!!!!do iorb=1,orbs%norbp
!!!!    iiorb=orbs%isorb+iorb
!!!!    ilr=onWhichAtom(iiorb)
!!!!    if(ilr==ilrold) then
!!!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!!!        jst=jst-op%noverlaps(iiorb)*ncount
!!!!    end if
!!!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!!!    !call dscal(ncount, 1.5d0, lhphi(ist), 1)
!!!!    do jorb=1,op%noverlaps(iiorb)
!!!!        jjorb=op%overlaps(jorb,iiorb)
!!!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
!!!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
!!!!        jst=jst+ncount
!!!!    end do
!!!!    ist=ist+ncount
!!!!    ilrold=ilr
!!!!end do
!!call cpu_time(t2)
!!time_daxpy=t2-t1
!!
!!if (verbose > 2) then
!!   call mpiallred(time_dsymm, 1, mpi_sum, mpi_comm_world, ierr)
!!   call mpiallred(time_daxpy, 1, mpi_sum, mpi_comm_world, ierr)
!!   if(iproc==0) write(*,'(a,es15.6)') 'time for dsymm',time_dsymm/dble(nproc)
!!   if(iproc==0) write(*,'(a,es15.6)') 'time for daxpy',time_daxpy/dble(nproc)
!!end if
!!
!!
!!iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
!!deallocate(ovrlp_minus_one_lagmat, stat=istat)
!!call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
!!iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
!!deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
!!call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
!!iall=-product(shape(ovrlp2))*kind(ovrlp2)
!!deallocate(ovrlp2, stat=istat)
!!call memocc(istat, iall, 'ovrlp2', subname)
!!iall=-product(shape(lagmat))*kind(lagmat)
!!deallocate(lagmat, stat=istat)
!!call memocc(istat, iall, 'lagmat', subname)
!!iall=-product(shape(phi))*kind(phi)
!!deallocate(phi, stat=istat)
!!call memocc(istat, iall, 'phi', subname)
!!iall=-product(shape(hphi))*kind(hphi)
!!deallocate(hphi, stat=istat)
!!call memocc(istat, iall, 'hphi', subname)
!!iall=-product(shape(phiWork))*kind(phiWork)
!!deallocate(phiWork, stat=istat)
!!call memocc(istat, iall, 'phiWork', subname)
!!
!!
!!end subroutine applyOrthoconstraintNonorthogonalCubic






!!subroutine localloewdin(iproc, nproc, lzd, orbs, op, ovrlp, lphiovrlp, lphi)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => overlapPowerMinusOneHalf
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(orbitals_data),intent(in):: orbs
!!  type(overlapParameters),intent(in):: op
!!  real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
!!  real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
!!  real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
!!  
!!  ! Local variables
!!  integer:: lwork, istat, iall, iorb, jorb, info, iiorb, korb, ilr, jlr, norbinlocreg, iorbinlocreg, jorbinlocreg
!!  integer:: iiorbinlocreg, jjorbinlocreg, klr, ist, jst, korbinlocreg, iilr, ilrold, ncount, jjorb, jjlr
!!  character(len=*),parameter:: subname='localloewdin'
!!  real(8),dimension(:),allocatable:: eval, work
!!  real(8),dimension(:,:),allocatable:: ovrlp2, ovrlp3, ovrlplocal
!!  real(8),dimension(:,:,:),allocatable:: tempArr
!!  real(8):: tt, dnrm2, ddot
!!
!!
!!  lphi=0.d0
!!
!!  do ilr=1,lzd%nlr
!!
!!      ! Determine number of orbitals in this locreg
!!      norbinlocreg=0
!!      do jorb=1,orbs%norb
!!          jlr=orbs%inwhichlocreg(jorb)
!!          if(jlr==ilr) norbinlocreg=norbinlocreg+1
!!      end do
!!
!!      ! Exctract the corresponding part of the overlap matrix
!!      allocate(ovrlplocal(norbinlocreg,norbinlocreg), stat=istat)
!!      call memocc(istat, ovrlplocal, 'ovrlplocal', subname)
!!      jorbinlocreg=0
!!      do jorb=1,orbs%norb
!!          jlr=orbs%inwhichlocreg(jorb)
!!          if(jlr==ilr) then
!!              jorbinlocreg=jorbinlocreg+1
!!              korbinlocreg=0
!!              do korb=1,orbs%norb
!!                  klr=orbs%inwhichlocreg(korb)
!!                  if(klr==ilr) then
!!                      korbinlocreg=korbinlocreg+1
!!                      ovrlplocal(korbinlocreg,jorbinlocreg)=ovrlp(korb,jorb)
!!                  end if
!!              end do
!!          end if
!!      end do
!!  
!!      ! Exact calculation of ovrlp**(-1/2)
!!      allocate(eval(norbinlocreg), stat=istat)
!!      call memocc(istat, eval, 'eval', subname)
!!      allocate(tempArr(norbinlocreg,norbinlocreg,2), stat=istat)
!!      call memocc(istat, tempArr, 'tempArr', subname)
!!
!!      lwork=1000*norbinlocreg
!!      allocate(work(lwork), stat=istat)
!!      call memocc(istat, work, 'work', subname)
!!      call dsyev('v', 'l', norbinlocreg, ovrlplocal(1,1), norbinlocreg, eval, work, lwork, info)
!!      if(info/=0) then
!!          write(*,'(a,i0)') 'ERROR in dsyev, info=', info
!!          stop
!!      end if
!!      iall=-product(shape(work))*kind(work)
!!      deallocate(work, stat=istat)
!!      call memocc(istat, iall, 'work', subname)
!!      
!!      ! Calculate S^{-1/2}. 
!!      ! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
!!      ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
!!      do iorb=1,norbinlocreg
!!          do jorb=1,norbinlocreg
!!              !tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
!!              tempArr(jorb,iorb,1)=ovrlplocal(jorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
!!          end do
!!      end do
!!      
!!      ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
!!      ! This will give S^{-1/2}.
!!      call dgemm('n', 't', norbinlocreg, norbinlocreg, norbinlocreg, 1.d0, ovrlplocal(1,1), &
!!           norbinlocreg, tempArr(1,1,1), norbinlocreg, 0.d0, tempArr(1,1,2), norbinlocreg)
!!      call dcopy(norbinlocreg**2, tempArr(1,1,2), 1, ovrlplocal(1,1), 1)
!!      
!!      
!!      iall=-product(shape(eval))*kind(eval)
!!      deallocate(eval, stat=istat)
!!      call memocc(istat, iall, 'eval', subname)
!!      iall=-product(shape(tempArr))*kind(tempArr)
!!      deallocate(tempArr, stat=istat)
!!      call memocc(istat, iall, 'tempArr', subname)
!!
!!      ! Orthogonalize the orbitals in this locreg
!!      ist=1
!!      jst=1
!!      ilrold=-1
!!      do iorb=1,orbs%norbp
!!          iiorb=orbs%isorb+iorb
!!          iilr=orbs%inwhichlocreg(iiorb)
!!          if(iilr==ilr) then
!!              if(iilr==ilrold) then
!!                  ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!                  jst=jst-op%noverlaps(iiorb)*ncount
!!              end if
!!              ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
!!              do jorb=1,op%noverlaps(iiorb)
!!                  jjorb=op%overlaps(jorb,iiorb)
!!                  jjlr=orbs%inwhichlocreg(jjorb)
!!                  if(jjlr==iilr) then
!!                      ! First determine number of iiorb and jjorb in this locreg
!!                      iiorbinlocreg=0
!!                      do korb=1,orbs%norb
!!                          klr=orbs%inwhichlocreg(korb)
!!                          if(klr==iilr) iiorbinlocreg=iiorbinlocreg+1
!!                          if(korb==iiorb) exit
!!                      end do
!!                      jjorbinlocreg=0
!!                      do korb=1,orbs%norb
!!                          klr=orbs%inwhichlocreg(korb)
!!                          if(klr==iilr) jjorbinlocreg=jjorbinlocreg+1
!!                          if(korb==jjorb) exit
!!                      end do
!!                      call daxpy(ncount, ovrlplocal(jjorbinlocreg,iiorbinlocreg), lphiovrlp(jst), 1, lphi(ist), 1)
!!                  end if
!!                  jst=jst+ncount
!!              end do
!!      
!!              ! Normalize
!!              tt=dnrm2(ncount, lphi(ist), 1)
!!              call dscal(ncount, 1/tt, lphi(ist), 1)
!!      
!!              ist=ist+ncount
!!              ilrold=iilr
!!          else
!!              ! only increase indices
!!              if(iilr==ilrold) then
!!                  ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!                  jst=jst-op%noverlaps(iiorb)*ncount
!!              end if
!!              ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
!!              do jorb=1,op%noverlaps(iiorb)
!!                  jst=jst+ncount
!!              end do
!!              ist=ist+ncount
!!              ilrold=iilr
!!          end if
!!      end do
!!
!!
!!      iall=-product(shape(ovrlplocal))*kind(ovrlplocal)
!!      deallocate(ovrlplocal, stat=istat)
!!      call memocc(istat, iall, 'ovrlplocal', subname)
!!
!!
!!
!!  end do
!!
!!
!!endsubroutine localloewdin

!!subroutine applyOrthoconstraintlocal(iproc, nproc, lzd, orbs, op, lagmat, lphiovrlp, lhphi)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => applyOrthoconstraintlocal
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(local_zone_descriptors),intent(in):: lzd
!!type(orbitals_data),intent(in):: orbs
!!type(overlapParameters),intent(in):: op
!!real(8),dimension(orbs%norb,orbs%norb),intent(inout):: lagmat
!!real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
!!real(8),dimension(orbs%npsidim_orbs),intent(out):: lhphi
!!
!!! Local variables
!!integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, info, i, istat, iall, ierr, iilr
!!real(8):: tt, t1, t2, time_dsymm, time_daxpy
!!character(len=*),parameter:: subname='applyOrthoconstraintlocal'
!!
!!
!!call cpu_time(t1)
!!do ilr=1,lzd%nlr
!!
!!    ist=1
!!    jst=1
!!    ilrold=-1
!!    do iorb=1,orbs%norbp
!!        iiorb=orbs%isorb+iorb
!!        iilr=orbs%inwhichlocreg(iiorb)
!!        if(iilr==ilr) then
!!            if(iilr==ilrold) then
!!                ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!                jst=jst-op%noverlaps(iiorb)*ncount
!!            end if
!!            ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
!!            !call dscal(ncount, 1.5d0, lhphi(ist), 1)
!!            do jorb=1,op%noverlaps(iiorb)
!!                jjorb=op%overlaps(jorb,iiorb)
!!                call daxpy(ncount, -lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
!!                jst=jst+ncount
!!            end do
!!            ist=ist+ncount
!!            ilrold=iilr
!!        else
!!            ! only increase index
!!            if(iilr==ilrold) then
!!                ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!                jst=jst-op%noverlaps(iiorb)*ncount
!!            end if
!!            ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
!!            do jorb=1,op%noverlaps(iiorb)
!!                jst=jst+ncount
!!            end do
!!            ist=ist+ncount
!!            ilrold=iilr
!!        end if
!!    end do
!!end do
!!
!!
!!call cpu_time(t2)
!!time_daxpy=t2-t1
!!
!!call mpiallred(time_daxpy, 1, mpi_sum, mpi_comm_world, ierr)
!!if(iproc==0) write(*,'(a,es15.6)') 'time for daxpy',time_daxpy/dble(nproc)
!!
!!
!!end subroutine applyOrthoconstraintlocal

!!subroutine calculate_overlap_matrix(iproc, orbs, collComms, psi_i, psi_j, ovrlp)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc
!!type(orbitals_data),intent(in):: orbs
!!type(collectiveComms),intent(in):: collComms
!!real(8),dimension(orbs%npsidim_orbs),intent(in):: psi_i, psi_j
!!real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
!!
!!! Local variables
!!integer:: iorb, jorb, ist, jst, ncnt_iorb, ncnt_jorb, iloc, jloc, i, ii, jj, ierr, iist, jjst, istat
!!logical:: istop, jstop
!!real(8):: tt
!!
!!iist=1
!!do iorb=1,orbs%norb
!!    jjst=1
!!    do jorb=1,orbs%norb
!!       ncnt_iorb=collComms%nvctr_par(iorb,iproc) 
!!       ncnt_jorb=collComms%nvctr_par(jorb,iproc) 
!!
!!       ii=0
!!       jj=0
!!       istop=.false.
!!       jstop=.false.
!!       tt=0.d0
!!       ist=iist
!!       jst=jjst
!!       !do i=1,max(ncnt_iorb,ncnt_jorb)
!!       !write(*,'(a,2i8,i10,2i12,2l)') 'iproc, jorb, jst, ncnt_jorb, orbs%npsidim, istop, jstop', iproc, jorb, jst, ncnt_jorb, orbs%npsidim, istop, jstop
!!       !!if(iproc==0 .and. iorb==1) then
!!       !!    do istat=iist,iist+ncnt_iorb
!!       !!        write(110,*) istat, collComms%indexarray(istat)
!!       !!    end do
!!       !!end if
!!       !!if(iproc==0 .and. jorb==5) then
!!       !!    do istat=jjst,jjst+ncnt_jorb
!!       !!        write(510,*) istat, collComms%indexarray(istat)
!!       !!    end do
!!       !!end if
!!
!!       do
!!           !if(jstop .and. jst>orbs%npsidim) write(*,'(a,5i12,2l)') 'iproc, jj, ncnt_jorb, jst, orbs%npsidim, istop, jstop', iproc, jj, ncnt_jorb, jst, orbs%npsidim, istop, jstop
!!           if(.not.istop) iloc=collComms%indexarray(ist)  
!!           if(.not.jstop) jloc=collComms%indexarray(jst)  
!!           !!write(1000*(iproc+1)+700+jorb,'(a,2i5,2i12,3x,2i12)') 'iorb, jorb, ist, jst, iloc, jloc', &
!!           !!    iorb, jorb, ist, jst, iloc, jloc
!!           !!if(iorb==1 .and. jorb==5) then
!!           !!    write(980+iproc,'(2i9,4x,2i12,2l4,2i9)') ist, jst, collComms%indexarray(ist), collComms%indexarray(jst), istop, jstop, ncnt_iorb, ncnt_jorb
!!           !!end if
!!           if(iloc==jloc) then
!!               !if(iorb==1 .and. jorb==5) then
!!               !    write(980+iproc,'(a,2i9,4x,2i12,2l4,2i9)') 'HERE1: ',ist, jst, collComms%indexarray(ist), collComms%indexarray(jst), istop, jstop, ncnt_iorb, ncnt_jorb
!!               !end if
!!               !if(istop) write(*,*) 'STRANGE: istop is true...'
!!               !if(jstop) write(*,*) 'STRANGE: jstop is true...'
!!               !if(iorb==1 .and. jorb==5) then
!!               !    write(iproc*1000+999,*) iloc
!!               !end if
!!               tt=tt+psi_i(ist)*psi_j(jst)
!!               ist=ist+1
!!               jst=jst+1
!!               ii=ii+1
!!               jj=jj+1
!!               !if(iorb==1 .and. jorb==5) then
!!               !    write(980+iproc,'(a,2i9,4x,2i12,2l4,2i9)') 'HERE2: ',ist, jst, collComms%indexarray(ist), collComms%indexarray(jst), istop, jstop, ncnt_iorb, ncnt_jorb
!!               !end if
!!           else if((iloc<jloc .or. jstop) .and. .not.istop) then
!!               ist=ist+1
!!               ii=ii+1
!!           else if((jloc<iloc .or. istop) .and. .not.jstop) then
!!               jst=jst+1
!!               jj=jj+1
!!           end if
!!           if(ii==ncnt_iorb) istop=.true.
!!           if(jj==ncnt_jorb) jstop=.true.
!!           if(istop .and. jstop) then
!!               !!if(iorb==1 .and. jorb==5) write(980+iproc,'(a)') 'exit since both stops are true...'
!!               exit
!!           end if
!!       end do
!!
!!       ovrlp(jorb,iorb)=tt
!!       jjst=jjst+ncnt_jorb
!!
!!    end do
!!    iist=iist+ncnt_iorb
!!end do
!!
!!
!!call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
!!
!!
!!end subroutine calculate_overlap_matrix

!!!!subroutine postCommsOverlapNew(iproc, nproc, orbs, op, lzd, phi, comon, timecommun, timeextract)
!!!!use module_base
!!!!use module_types
!!!!use module_interfaces, exceptThisOne => postCommsOverlapNew
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc
!!!!type(orbitals_data),intent(in):: orbs
!!!!type(overlapParameters),intent(in):: op
!!!!type(local_zone_descriptors),intent(in):: lzd
!!!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: phi
!!!!type(p2pComms),intent(inout):: comon
!!!!real(8),intent(inout):: timecommun, timeextract
!!!!
!!!!! Local variables
!!!!integer:: jproc, jorb, mpisource, istsource, ncount, mpidest, istdest, tag, nsends, ierr, irecv, isend, iall, orbsource, orbdest
!!!!integer:: i, indsource, ind, klr, kkorb, istat, iorb
!!!!integer,dimension(:),allocatable:: istextracted
!!!!logical:: done
!!!!real(8):: t1, t2
!!!!
!!!!
!!!!
!!!!! First only post receives
!!!!irecv=0
!!!!do jproc=0,nproc-1
!!!!    do jorb=1,comon%noverlaps(jproc)
!!!!        mpisource=comon%comarr(1,jorb,jproc)
!!!!        istsource=comon%comarr(2,jorb,jproc)
!!!!        ncount=comon%comarr(3,jorb,jproc)
!!!!        mpidest=comon%comarr(4,jorb,jproc)
!!!!        istdest=comon%comarr(5,jorb,jproc)
!!!!        tag=comon%comarr(6,jorb,jproc)
!!!!        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
!!!!        !irecv=irecv+1
!!!!        if(iproc==mpidest .and. nproc >1) then
!!!!            irecv=irecv+1
!!!!            !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!!            !call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!!!!            !     mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
!!!!            tag=jorb
!!!!            t1=mpi_wtime()
!!!!            !write(*,'(2(a,i0))') 'iproc=',iproc,' receives data at ',istdest
!!!!            call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!!!!                 mpi_comm_world, comon%requests(irecv,2), ierr)
!!!!            t2=mpi_wtime()
!!!!            timecommun=timecommun+t2-t1
!!!!        !else
!!!!        !     comon%requests(irecv,2)=mpi_request_null
!!!!        end if
!!!!    end do
!!!!end do
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call timing(iproc,'p2pOrtho_post ','ON')
!!!!
!!!!! Now the sends
!!!!isend=0
!!!!do jproc=0,nproc-1
!!!!    do jorb=1,comon%noverlaps(jproc)
!!!!        mpisource=comon%comarr(1,jorb,jproc)
!!!!        istsource=comon%comarr(2,jorb,jproc)
!!!!        ncount=comon%comarr(3,jorb,jproc)
!!!!        mpidest=comon%comarr(4,jorb,jproc)
!!!!        istdest=comon%comarr(5,jorb,jproc)
!!!!        tag=comon%comarr(6,jorb,jproc)
!!!!        ! The orbitals are on different processes, so we need a point to point communication.
!!!!        !isend=isend+1
!!!!        if(iproc==mpisource .and. nproc >1) then
!!!!            isend=isend+1
!!!!            t2=mpi_wtime()
!!!!            timeextract=timeextract+t2-t1
!!!!            tag=jorb
!!!!            t1=mpi_wtime()
!!!!            call mpi_irsend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!!!!                 mpi_comm_world, comon%requests(isend,1), ierr)
!!!!            !call mpi_rsend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!!!!            !     mpi_comm_world, ierr)
!!!!            t2=mpi_wtime()
!!!!            timecommun=timecommun+t2-t1
!!!!         else if (nproc==1) then
!!!!            call vcopy(ncount,comon%sendBuf(istsource),1,comon%recvBuf(istdest),1)
!!!!        end if
!!!!    end do
!!!!end do
!!!!!!comon%nsend=isend
!!!!!!if(isend/=comon%nsend) then
!!!!!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR in process ',iproc,': irecv/=comon%nrecv',irecv,comon%nsend
!!!!!!    stop
!!!!!!end if
!!!!
!!!!call timing(iproc,'p2pOrtho_post ','OF')
!!!!
!!!!end subroutine postCommsOverlapNew

!!!!subroutine collectnew(iproc, nproc, comon, mad, op, orbs, lzd, &
!!!!   nsendbuf, sendbuf, nrecvbuf, recvbuf, timecommunp2p, timecommuncoll, timecompress)
!!!!use module_base
!!!!use module_types
!!!!use module_interfaces, exceptThisOne => collectnew
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
!!!!type(p2pComms),intent(inout):: comon
!!!!type(matrixDescriptors),intent(in):: mad
!!!!type(overlapParameters),intent(in):: op
!!!!type(orbitals_data),intent(in):: orbs
!!!!type(local_zone_descriptors),intent(in):: lzd
!!!!real(8),dimension(nsendbuf),intent(in):: sendbuf
!!!!real(8),dimension(nrecvbuf),intent(inout):: recvbuf
!!!!real(8),intent(inout):: timecommunp2p, timecommuncoll, timecompress
!!!!
!!!!! Local variables
!!!!integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
!!!!real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
!!!!integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
!!!!real(8):: ddot, t1, t2
!!!!character(len=*),parameter:: subname='collectnew'
!!!!logical,dimension(:),allocatable:: done
!!!!logical:: received
!!!!
!!!!allocate(indcomplete(comon%nrecv), stat=istat)
!!!!call memocc(istat, indcomplete, 'indcomplete', subname)
!!!!allocate(done(comon%nrecv), stat=istat)
!!!!call memocc(istat, done, 'done', subname)
!!!!done=.false.
!!!!
!!!!allocate(indexarray(comon%nrecv), stat=istat)
!!!!call memocc(istat, indexarray, 'indexarray', subname)
!!!!do i=1,comon%nrecv
!!!!   indexarray(i)=i
!!!!end do
!!!!
!!!!call timing(iproc,'p2pOrtho_wait ','ON')
!!!!
!!!!! Wait for the sends to complete.
!!!!if (nproc > 1) then
!!!!   t1=mpi_wtime()
!!!!   nsend=0
!!!!   if(comon%nsend>0) then
!!!!       waitLoopSend: do
!!!!          !!call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!!!          !!nsend=nsend+ncomplete
!!!!          !!if(nsend==comon%nsend) exit waitLoopSend
!!!!          call mpi_waitany(comon%nsend-nsend, comon%requests(1,1), ind, mpi_status_ignore, ierr)
!!!!          nsend=nsend+1
!!!!          do i=ind,comon%nsend-nsend
!!!!             comon%requests(i,1)=comon%requests(i+1,1)
!!!!          end do
!!!!          if(nsend==comon%nsend) exit waitLoopSend
!!!!       end do waitLoopSend
!!!!   end if
!!!!   t2=mpi_wtime()
!!!!   timecommunp2p=timecommunp2p+t2-t1
!!!!
!!!!
!!!!   nrecv=0
!!!!   if(comon%nrecv>0) then
!!!!       waitLoopRecv: do
!!!!          t1=mpi_wtime()
!!!!          call mpi_waitany(comon%nrecv-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
!!!!          !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
!!!!          !ind=1
!!!!          t2=mpi_wtime()
!!!!          timecommunp2p=timecommunp2p+t2-t1
!!!!          ncomplete=1
!!!!          received=.true.
!!!!          if(received) then
!!!!             nrecv=nrecv+ncomplete
!!!!             !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
!!!!             !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
!!!!             do i=ind,comon%nrecv-nrecv
!!!!                comon%requests(i,2)=comon%requests(i+1,2)
!!!!                indexarray(i)=indexarray(i+1)
!!!!             end do
!!!!             !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
!!!!             if(nrecv==comon%nrecv) exit waitLoopRecv
!!!!          end if
!!!!       end do waitLoopRecv
!!!!   end if
!!!!end if
!!!!!write(*,'(a,i0,a)') 'iproc=',iproc,' is here'
!!!!
!!!!iall=-product(shape(done))*kind(done)
!!!!deallocate(done, stat=istat)
!!!!call memocc(istat, iall, 'done', subname)
!!!!
!!!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!!!deallocate(indcomplete, stat=istat)
!!!!call memocc(istat, iall, 'indcomplete', subname)
!!!!
!!!!iall=-product(shape(indexarray))*kind(indexarray)
!!!!deallocate(indexarray, stat=istat)
!!!!call memocc(istat, iall, 'indexarray', subname)
!!!!
!!!!
!!!!call timing(iproc,'p2pOrtho_wait ','OF')
!!!!
!!!!
!!!!!!allocate(indcomplete(comon%nsend), stat=istat)
!!!!!!call memocc(istat, indcomplete, 'indcomplete', subname)
!!!!
!!!!!!do iorb=1,orbs%norb
!!!!!!  do jorb=1,orbs%norb
!!!!!!    write(400+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!!!!  end do
!!!!!!end do
!!!!
!!!!!!! Now the sends
!!!!!!nsend=0
!!!!!!waitLoopSend: do
!!!!!!    call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!!!!!    nsend=nsend+ncomplete
!!!!!!    if(nsend==comon%nsend) exit waitLoopSend
!!!!!!end do waitLoopSend
!!!!
!!!!!write(*,*) 'here: iproc', iproc
!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!call mpi_finalize(ierr)
!!!!!stop
!!!!
!!!!!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!!!!!deallocate(indcomplete, stat=istat)
!!!!!!call memocc(istat, iall, 'indcomplete', subname)
!!!!
!!!!!do iorb=1,orbs%norb
!!!!!  do jorb=1,orbs%norb
!!!!!    write(90+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!!!  end do
!!!!!end do
!!!!
!!!!
!!!!!!! Communicate the matrix
!!!!!!allocate(ovrlpCompressed(mad%nvctr), stat=istat)
!!!!!!call memocc(istat, ovrlpCompressed, 'ovrlpCompressed', subname)
!!!!!!allocate(sendcounts(0:nproc-1), stat=istat)
!!!!!!call memocc(istat, sendcounts, 'sendcounts', subname)
!!!!!!allocate(displs(0:nproc-1), stat=istat)
!!!!!!call memocc(istat, displs, 'displs', subname)
!!!!!!
!!!!!!t1=mpi_wtime()
!!!!!!call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed, sendcounts, displs)
!!!!!!t2=mpi_wtime()
!!!!!!timecompress=timecompress+t2-t1     
!!!!!!
!!!!!!allocate(ovrlpCompressed2(mad%nvctr), stat=istat)
!!!!!!call memocc(istat, ovrlpCompressed2, 'ovrlpCompressed2', subname)
!!!!!!
!!!!!!t1=mpi_wtime()
!!!!!!call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
!!!!!!     sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
!!!!!!t2=mpi_wtime()
!!!!!!timecommuncoll=timecommuncoll+t2-t1     
!!!!!!
!!!!!!t1=mpi_wtime()
!!!!!!call uncompress_matrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)
!!!!!!t2=mpi_wtime()
!!!!!!timecompress=timecompress+t2-t1     
!!!!!!
!!!!!!!do iorb=1,orbs%norb
!!!!!!!  do jorb=1,orbs%norb
!!!!!!!    write(100+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!!!!!!  end do
!!!!!!!end do
!!!!!!
!!!!!!iall=-product(shape(ovrlpCompressed))*kind(ovrlpCompressed)
!!!!!!deallocate(ovrlpCompressed, stat=istat)
!!!!!!call memocc(istat, iall, 'ovrlpCompressed', subname)
!!!!!!iall=-product(shape(ovrlpCompressed2))*kind(ovrlpCompressed2)
!!!!!!deallocate(ovrlpCompressed2, stat=istat)
!!!!!!call memocc(istat, iall, 'ovrlpCompressed2', subname)
!!!!!!iall=-product(shape(sendcounts))*kind(sendcounts)
!!!!!!deallocate(sendcounts, stat=istat)
!!!!!!call memocc(istat, iall, 'sendcounts', subname)
!!!!!!iall=-product(shape(displs))*kind(displs)
!!!!!!deallocate(displs, stat=istat)
!!!!!!call memocc(istat, iall, 'displs', subname)
!!!!
!!!!
!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!call mpi_finalize(ierr)
!!!!
!!!!end subroutine collectnew

!!!subroutine gatherOrbitals2(iproc, nproc, comon)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  integer,intent(in):: iproc, nproc
!!!  type(p2pComms),intent(inout):: comon
!!!
!!!  ! Local variables
!!!  integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
!!!  integer,dimension(mpi_status_size):: stat
!!!  logical:: sendComplete, receiveComplete
!!!
!!!
!!!!!! Check whether the communications have completed.
!!!  !!nfast=0
!!!  !!nsameproc=0
!!!  !!testLoop: do
!!!  !!    do jproc=0,nproc-1
!!!  !!        do jorb=1,comon%noverlaps(jproc)
!!!  !!            if(comon%communComplete(jorb,jproc)) cycle
!!!  !!            call mpi_test(comon%comarr(7,jorb,jproc), sendComplete, stat, ierr)
!!!  !!            call mpi_test(comon%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
!!!  !!            ! Attention: mpi_test is a local function.
!!!  !!            if(sendComplete .and. receiveComplete) comon%communComplete(jorb,jproc)=.true.
!!!  !!        end do
!!!  !!    end do
!!!  !!    ! If we made it until here, either all all the communication is
!!!  !!    ! complete or we better wait for each single orbital.
!!!  !!    exit testLoop
!!!  !!end do testLoop
!!!  !!
!!!!!! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!!  !!call mpiallred(comon%communComplete(1,0), nproc*maxval(comon%noverlaps), mpi_land, mpi_comm_world, ierr)
!!!
!!!
!!!  ! Wait for the communications that have not completed yet
!!!  nslow=0
!!!  do jproc=0,nproc-1
!!!     do jorb=1,comon%noverlaps(jproc)
!!!        if(comon%communComplete(jorb,jproc)) then
!!!           mpisource=comon%comarr(1,jorb,jproc)
!!!           mpidest=comon%comarr(4,jorb,jproc)
!!!           if(mpisource==mpidest) then
!!!              !nsameproc=nsameproc+1
!!!           else
!!!              !nfast=nfast+1
!!!           end if
!!!           cycle
!!!        end if
!!!        !write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
!!!        nslow=nslow+1
!!!        !if(iproc==mpidest) call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!        call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!        !if(iproc==mpisource) call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!        call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!        comon%communComplete(jorb,jproc)=.true.
!!!     end do
!!!  end do
!!!
!!!  !call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!!  !call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!!  !call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!!  !call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!!  !nfast=nint(dble(nfast)/dble(nproc))
!!!  !nfast=nint(dble(nslow)/dble(nproc))
!!!  !nfast=nint(dble(nsameproc)/dble(nproc))
!!!  !if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!  !                       nfast, ' could be overlapped with computation.'
!!!  !if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!!
!!!
!!!end subroutine gatherOrbitals2

!!!!subroutine postCommsOverlap(iproc, nproc, comon)
!!!!  use module_base
!!!!  use module_types
!!!!  use module_interfaces, exceptThisOne => postCommsOverlap
!!!!  implicit none
!!!!
!!!!  ! Calling arguments
!!!!  integer,intent(in):: iproc, nproc
!!!!  type(p2pComms),intent(inout):: comon
!!!!
!!!!  ! Local variables
!!!!  integer:: jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, nsends, nreceives, ierr
!!!!
!!!!!!! THIS IS THE ORIGINAL ##################################
!!!!  !nsends=0
!!!!  !nreceives=0
!!!!  !comon%communComplete=.false.
!!!!  !do jproc=0,nproc-1
!!!!  !    !write(*,'(3(a,i0))') 'iproc=',iproc,', jproc=',jproc,', comon%noverlaps(jproc)=', comon%noverlaps(jproc)
!!!!  !    do iorb=1,comon%noverlaps(jproc)
!!!!  !        mpisource=comon%comarr(1,iorb,jproc)
!!!!  !        istsource=comon%comarr(2,iorb,jproc)
!!!!  !        ncount=comon%comarr(3,iorb,jproc)
!!!!  !        mpidest=comon%comarr(4,iorb,jproc)
!!!!  !        istdest=comon%comarr(5,iorb,jproc)
!!!!  !        tag=comon%comarr(6,iorb,jproc)
!!!!  !        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
!!!!  !        if(mpisource/=mpidest) then
!!!!  !            ! The orbitals are on different processes, so we need a point to point communication.
!!!!  !            if(iproc==mpisource) then
!!!!  !                !write(*,'(6(a,i0))') 'overlap: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
!!!!  !                call mpi_isend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!!!!  !                     mpi_comm_world, comon%comarr(7,iorb,jproc), ierr)
!!!!  !                !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
!!!!  !                comon%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
!!!!  !                nsends=nsends+1
!!!!  !            else if(iproc==mpidest) then
!!!!  !                !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!!  !                call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!!!!  !                     mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
!!!!  !                comon%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
!!!!  !                nreceives=nreceives+1
!!!!  !            else
!!!!  !                comon%comarr(7,iorb,jproc)=mpi_request_null
!!!!  !                comon%comarr(8,iorb,jproc)=mpi_request_null
!!!!  !            end if
!!!!  !        else
!!!!  !            ! The orbitals are on the same process, so simply copy them.
!!!!  !            if(iproc==mpisource) then
!!!!  !                call dcopy(ncount, comon%sendBuf(istsource), 1, comon%recvBuf(istdest), 1)
!!!!  !                !write(*,'(6(a,i0))') 'overlap: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
!!!!  !                comon%comarr(7,iorb,jproc)=mpi_request_null
!!!!  !                comon%comarr(8,iorb,jproc)=mpi_request_null
!!!!  !                nsends=nsends+1
!!!!  !                nreceives=nreceives+1
!!!!  !                comon%communComplete(iorb,mpisource)=.true.
!!!!  !            else
!!!!  !                comon%comarr(7,iorb,jproc)=mpi_request_null
!!!!  !                comon%comarr(8,iorb,jproc)=mpi_request_null
!!!!  !                comon%communComplete(iorb,mpisource)=.true.
!!!!  !            end if
!!!!  !
!!!!  !        end if
!!!!  !    end do
!!!!  !end do
!!!!
!!!!
!!!!
!!!!  !! NEW ####################################
!!!!  ! First only post receives
!!!!  nsends=0
!!!!  nreceives=0
!!!!  comon%communComplete=.false.
!!!!  do jproc=0,nproc-1
!!!!     !write(*,'(3(a,i0))') 'iproc=',iproc,', jproc=',jproc,', comon%noverlaps(jproc)=', comon%noverlaps(jproc)
!!!!     do iorb=1,comon%noverlaps(jproc)
!!!!        mpisource=comon%comarr(1,iorb,jproc)
!!!!        istsource=comon%comarr(2,iorb,jproc)
!!!!        ncount=comon%comarr(3,iorb,jproc)
!!!!        mpidest=comon%comarr(4,iorb,jproc)
!!!!        istdest=comon%comarr(5,iorb,jproc)
!!!!        tag=comon%comarr(6,iorb,jproc)
!!!!        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
!!!!        if(mpisource/=mpidest) then
!!!!           ! The orbitals are on different processes, so we need a point to point communication.
!!!!           if(iproc==mpisource) then
!!!!           else if(iproc==mpidest) then
!!!!              !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount,&
!!!!              !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!!              call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!!!!                   mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
!!!!              comon%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
!!!!              nreceives=nreceives+1
!!!!           end if
!!!!        end if
!!!!     end do
!!!!  end do
!!!!
!!!!
!!!!  ! Now the rest.
!!!!  do jproc=0,nproc-1
!!!!     !write(*,'(3(a,i0))') 'iproc=',iproc,', jproc=',jproc,', comon%noverlaps(jproc)=', comon%noverlaps(jproc)
!!!!     do iorb=1,comon%noverlaps(jproc)
!!!!        mpisource=comon%comarr(1,iorb,jproc)
!!!!        istsource=comon%comarr(2,iorb,jproc)
!!!!        ncount=comon%comarr(3,iorb,jproc)
!!!!        mpidest=comon%comarr(4,iorb,jproc)
!!!!        istdest=comon%comarr(5,iorb,jproc)
!!!!        tag=comon%comarr(6,iorb,jproc)
!!!!        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
!!!!        if(mpisource/=mpidest) then
!!!!           ! The orbitals are on different processes, so we need a point to point communication.
!!!!           if(iproc==mpisource) then
!!!!              !write(*,'(6(a,i0))') 'overlap: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
!!!!              call mpi_isend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!!!!                   mpi_comm_world, comon%comarr(7,iorb,jproc), ierr)
!!!!              !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
!!!!              comon%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
!!!!              nsends=nsends+1
!!!!           else if(iproc==mpidest) then
!!!!           else
!!!!              comon%comarr(7,iorb,jproc)=mpi_request_null
!!!!              comon%comarr(8,iorb,jproc)=mpi_request_null
!!!!           end if
!!!!        else
!!!!           ! The orbitals are on the same process, so simply copy them.
!!!!           if(iproc==mpisource) then
!!!!              call dcopy(ncount, comon%sendBuf(istsource), 1, comon%recvBuf(istdest), 1)
!!!!              !write(*,'(6(a,i0))') 'overlap: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
!!!!              comon%comarr(7,iorb,jproc)=mpi_request_null
!!!!              comon%comarr(8,iorb,jproc)=mpi_request_null
!!!!              nsends=nsends+1
!!!!              nreceives=nreceives+1
!!!!              comon%communComplete(iorb,mpisource)=.true.
!!!!           else
!!!!              comon%comarr(7,iorb,jproc)=mpi_request_null
!!!!              comon%comarr(8,iorb,jproc)=mpi_request_null
!!!!              comon%communComplete(iorb,mpisource)=.true.
!!!!           end if
!!!!
!!!!        end if
!!!!     end do
!!!!  end do
!!!!
!!!!
!!!!
!!!!end subroutine postCommsOverlap




!!subroutine determine_overlapParameters_fast(iproc, nproc, orbs, lzd, op)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(out):: op
!!
!!  ! Local variables
!!  integer:: istat, iall, n0, iorb, ii, iiorb, ilr, jj, j0, j1, i3, i2, i1, i0, ii1, ii2, ii3, iseg, i, inumber, ierr
!!  integer,dimension(:),allocatable:: numbers, overlaps_tmp
!!  integer,dimension(:,:,:),allocatable:: orbitalnumbers
!!  character(len=*),parameter:: subname='determine_overlapParameters_fast'
!!
!!
!!
!!  allocate(numbers(orbs%norb), stat=istat)
!!  call memocc(istat, numbers, 'numbers', subname)
!!
!!
!!  ! Determine the reference number n0
!!  n0=(1+orbs%norb)*orbs%norb/2+1
!!  if(n0>2**30) then
!!      stop 'ERROR: this will give an integer overflow!'
!!  end if
!!
!!  ! Determine all numbers
!!  ii=0
!!  do iorb=1,orbs%norb
!!      ii=ii+iorb
!!      numbers(iorb)=n0+ii
!!      !if(iproc==0) write(*,*) 'iorb, numbers(iorb)', iorb, numbers(iorb)
!!  end do
!!
!!
!!  allocate(orbitalnumbers(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
!!  call memocc(istat, orbitalnumbers, 'orbitalnumbers', subname)
!!  call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1), orbitalnumbers(0,0,0))
!!  
!!  ! Assign the orbital numbers to the all grid points where they extend. Fine part not needed
!!  ! since the fine region is contained in the coarse region.
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      ilr=orbs%inwhichlocreg(iiorb)
!!      do iseg=1,lzd%llr(ilr)%wfd%nseg_c
!!          jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
!!          j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
!!          j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
!!          ii=j0-1
!!          i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
!!          ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
!!          i2=ii/(lzd%llr(ilr)%d%n1+1)
!!          i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
!!          i1=i0+j1-j0
!!          do i=i0,i1
!!              ii1=i+lzd%llr(ilr)%ns1
!!              ii2=i2+lzd%llr(ilr)%ns2
!!              ii3=i3+lzd%llr(ilr)%ns3
!!              orbitalnumbers(ii1,ii2,ii3)=numbers(iiorb)
!!          end do
!!      end do
!!  end do
!!
!!
!!  ! Communicate the grid among all processes
!!  call mpiallred(orbitalnumbers(0,0,0), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1), &
!!       mpi_sum, mpi_comm_world, ierr)
!!
!!
!!  allocate(overlaps_tmp(orbs%norb), stat=istat)
!!  call memocc(istat, overlaps_tmp, 'overlaps_tmp', subname)
!! 
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      ilr=orbs%inwhichlocreg(iiorb)
!!      do iseg=1,lzd%llr(ilr)%wfd%nseg_c
!!          jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
!!          j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
!!          j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
!!          ii=j0-1
!!          i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
!!          ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
!!          i2=ii/(lzd%llr(ilr)%d%n1+1)
!!          i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
!!          i1=i0+j1-j0
!!          do i=i0,i1
!!              ii1=i+lzd%llr(ilr)%ns1
!!              ii2=i2+lzd%llr(ilr)%ns2
!!              ii3=i3+lzd%llr(ilr)%ns3
!!              inumber=orbitalnumbers(ii1,ii2,ii3)
!!          end do
!!      end do
!!  end do
!!
!!
!!  iall=-product(shape(numbers))*kind(numbers)
!!  deallocate(numbers, stat=istat)
!!  call memocc(istat, iall, 'numbers', subname)
!!
!!  iall=-product(shape(orbitalnumbers))*kind(orbitalnumbers)
!!  deallocate(orbitalnumbers, stat=istat)
!!  call memocc(istat, iall, 'orbitalnumbers', subname)
!!
!!  iall=-product(shape(overlaps_tmp))*kind(overlaps_tmp)
!!  deallocate(overlaps_tmp, stat=istat)
!!  call memocc(istat, iall, 'overlaps_tmp', subname)
!!
!!end subroutine determine_overlapParameters_fast




!!!subroutine expandOrbital2(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  integer,intent(in):: iproc, nproc
!!!  type(orbitals_data),intent(in):: orbs
!!!  type(input_variables),intent(in):: input
!!!  type(local_zone_descriptors),intent(in):: lzd
!!!  type(overlapParameters),intent(in):: op
!!!  type(p2pComms),intent(in):: comon
!!!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
!!!
!!!  ! Local variables
!!!  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest, m, ii
!!!  integer:: start, iseg, istart, iend, ncount, kseg, kstart, kend, jst2, kold, ifine, isend
!!!  integer:: igrid 
!!!!!real(8):: t1, t2, time1, time2
!!!
!!!!!t1=mpi_wtime()
!!!  lphiovrlp=0.d0
!!!!!t2=mpi_wtime()
!!!!!time1=t2-t1
!!!
!!!!!t1=mpi_wtime()
!!!  ind=1
!!!  ilrold=-1
!!!  do iorb=1,orbs%norbp
!!!     iiorb=orbs%isorb+iorb
!!!     ilr=orbs%inwhichlocreg(iiorb)
!!!     if(ilr==ilrold) cycle
!!!     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!!     do jorb=1,op%noverlaps(iiorb)
!!!        jjorb=op%overlaps(jorb,iiorb)
!!!        ! Starting index of orbital jjorb
!!!        jst=op%indexInRecvBuf(iorb,jjorb)
!!!        ldim=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
!!!        !! THIS IS THE OLD VERSION
!!!        !!do i=0,ldim-1
!!!        !!    indDest=ind+op%indexExpand(jst+i)-1
!!!        !!    lphiovrlp(indDest)=comon%recvBuf(jst+i)
!!!        !!end do
!!!        !! THIS IS NEW
!!!        !!m=mod(ldim,7)
!!!        !!if(m/=0) then
!!!        !!   do i=0,m-1
!!!        !!        ii=jst+i
!!!        !!        indDest=ind+op%indexExpand(ii)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii)
!!!        !!   end do
!!!        !!end if
!!!        !!do i=m,ldim-1,7
!!!        !!        ii=jst+i
!!!        !!        indDest=ind+op%indexExpand(ii+0)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii+0)
!!!        !!        indDest=ind+op%indexExpand(ii+1)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii+1)
!!!        !!        indDest=ind+op%indexExpand(ii+2)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii+2)
!!!        !!        indDest=ind+op%indexExpand(ii+3)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii+3)
!!!        !!        indDest=ind+op%indexExpand(ii+4)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii+4)
!!!        !!        indDest=ind+op%indexExpand(ii+5)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii+5)
!!!        !!        indDest=ind+op%indexExpand(ii+6)-1
!!!        !!        lphiovrlp(indDest)=comon%recvBuf(ii+6)
!!!        !!end do
!!!        !! THIS IS THE NEWEST (NOT OPTIMIZED, could probably store the keyv once and for all)
!!!        jst2=0
!!!        kold = 1
!!!        do iseg=1,op%wfd_overlap(jorb,iorb)%nseg_c
!!!            istart=op%wfd_overlap(jorb,iorb)%keyglob(1,iseg)
!!!            iend=op%wfd_overlap(jorb,iorb)%keyglob(2,iseg)
!!!            ncount=iend-istart+1  !the overlap segment is always contained inside the Llrs segments, so ncount is ok
!!!            inner_loop: do kseg=kold,lzd%llr(ilr)%wfd%nseg_c
!!!               kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg)
!!!               kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg)
!!!               if(kstart <= iend .and. kend >= istart) then
!!!                  kold = kseg
!!!                  start = lzd%llr(ilr)%wfd%keyvglob(kseg) + max(0,istart-kstart)
!!!                  call dcopy(ncount, comon%recvBuf(jst+jst2), 1, lphiovrlp(start+ind-1), 1)
!!!                  jst2=jst2+ncount
!!!                  exit inner_loop
!!!               end if
!!!            end do inner_loop
!!!        end do
!!!        if(jst2 .ne. op%wfd_overlap(jorb,iorb)%nvctr_c) then
!!!          print *, 'jst2 = ',jst2,'not equal to ldim = ',op%wfd_overlap(jorb,iorb)%nvctr_c
!!!          stop
!!!        end if
!!!        !Do Fine grid
!!!        jst2=0
!!!        kold = 1
!!!        do iseg=1,op%wfd_overlap(jorb,iorb)%nseg_f
!!!            istart=op%wfd_overlap(jorb,iorb)%keyglob(1,iseg+op%wfd_overlap(jorb,iorb)%nseg_c)
!!!            iend=op%wfd_overlap(jorb,iorb)%keyglob(2,iseg+op%wfd_overlap(jorb,iorb)%nseg_c)
!!!            ncount=7*(iend-istart+1)  !the overlap segment is always contained inside the Llrs segments, so ncount is ok
!!!            inner_loop2: do kseg=kold,lzd%llr(ilr)%wfd%nseg_f
!!!               kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg+lzd%llr(ilr)%wfd%nseg_c)
!!!               kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg+lzd%llr(ilr)%wfd%nseg_c)
!!!               if(kstart <= iend .and. kend >= istart) then
!!!                  kold = kseg
!!!                  start = lzd%llr(ilr)%wfd%nvctr_c+(lzd%llr(ilr)%wfd%keyvglob(kseg+lzd%llr(ilr)%wfd%nseg_c) +&
!!!                          max(0,istart-kstart)-1)*7
!!!                  call dcopy(ncount,comon%recvBuf(jst+jst2+op%wfd_overlap(jorb,iorb)%nvctr_c),1,&
!!!                          lphiovrlp(ind+start),1)
!!!                  jst2=jst2+ncount
!!!                  exit inner_loop2
!!!               end if
!!!            end do inner_loop2
!!!        end do
!!!        if(jst2 .ne. 7*op%wfd_overlap(jorb,iorb)%nvctr_f) then
!!!          print *, 'jst2 = ',jst2,'not equal to ldim = ',7*op%wfd_overlap(jorb,iorb)%nvctr_f
!!!          stop
!!!        end if
!!!        ind=ind+gdim
!!!     end do
!!!     ilrold=ilr
!!!  end do
!!!!!t2=mpi_wtime()
!!!!!time2=t2-t1
!!!!!if(iproc==0) write(*,*) 'time1, time2', time1, time2
!!!
!!!end subroutine expandOrbital2



!!!subroutine determineOverlapsSphere(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  integer,intent(in):: iproc, nproc
!!!  type(orbitals_data),intent(in):: orbs
!!!  type(local_zone_descriptors),intent(in):: lzd
!!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!!  type(overlapParameters),intent(out):: op
!!!  type(p2pComms),intent(out):: comon
!!!
!!!  ! Local variables
!!!  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
!!!  integer:: js1, je1, js2, je2, js3, je3, iiorb
!!!  logical:: ovrlpx, ovrlpy, ovrlpz
!!!  real(8):: dx, dy, dz, rr
!!!
!!!  ! Initialize to some value which will never be used.
!!!  op%overlaps=-1
!!!  !!comon%overlaps=-1
!!!
!!!  iiorb=0
!!!  do jproc=0,nproc-1
!!!     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!!     ilrold=-1
!!!     do iorb=1,orbs%norb_par(jproc,0)
!!!        ioverlaporb=0 ! counts the overlaps for the given orbital.
!!!        iiorb=iiorb+1 ! counts the total orbitals
!!!        ilr=onWhichAtom(iiorb)
!!!        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!!        do jorb=1,orbs%norb
!!!           jlr=onWhichAtom(jorb)
!!!           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!!           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!!           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
!!!           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
!!!           dx=(lzd%llr(ilr)%locregCenter(1)-lzd%llr(jlr)%locregCenter(1))**2
!!!           dy=(lzd%llr(ilr)%locregCenter(2)-lzd%llr(jlr)%locregCenter(2))**2
!!!           dz=(lzd%llr(ilr)%locregCenter(3)-lzd%llr(jlr)%locregCenter(3))**2
!!!           rr=(lzd%llr(ilr)%locrad+lzd%llr(jlr)%locrad)**2
!!!           if(dx+dy+dz<=rr) then
!!!              !if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!!              ioverlaporb=ioverlaporb+1
!!!              op%overlaps(ioverlaporb,iiorb)=jorb
!!!              if(ilr/=ilrold) then
!!!                 ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
!!!                 ! would get the same orbitals again. Therefore the counter is not increased
!!!                 ! in that case.
!!!                 ioverlapMPI=ioverlapMPI+1
!!!                 !!comon%overlaps(ioverlapMPI,jproc)=jorb
!!!              end if
!!!           end if
!!!        end do
!!!        !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
!!!        ilrold=ilr
!!!     end do
!!!     !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
!!!  end do
!!!
!!!end subroutine determineOverlapsSphere



!!subroutine getStartingIndicesGlobal(iiorbx, jjorbx, op, orbs, ist, jst, ncount)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => getStartingIndicesGlobal
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in):: iiorbx, jjorbx
!!  type(overlapParameters),intent(in):: op
!!  type(orbitals_data),intent(in):: orbs
!!  integer,intent(out):: ist, jst, ncount
!!  
!!  ! Local variables
!!  integer:: iiorb, jjorb, iorb, jorb
!!
!!  ! This only works localy on a process, i.e. we can only get the informations related
!!  ! to the currect process (but this is enough)
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      do jorb=1,op%noverlaps(iiorb)
!!          jjorb=op%overlaps(jorb,iiorb)
!!          if(iiorb==iiorbx .and. jjorb==jjorbx) then
!!              call getStartingIndices(iorb, jorb, op, orbs, ist, jst)
!!              ncount=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
!!          end if
!!      end do
!!  end do
!!
!!end subroutine getStartingIndicesGlobal


!!subroutine determineOverlaps(iproc, nproc, orbs, lzd, op, comon)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(out):: op
!!  type(p2pComms),intent(out):: comon
!!
!!  ! Local variables
!!  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, iiorb
!!!  integer ::  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
!!!  logical:: ovrlpx, ovrlpy, ovrlpz
!!  logical :: isoverlap
!!
!!  ! Initialize to some value which will never be used.
!!  op%overlaps=-1
!!  !!comon%overlaps=-1
!!
!!  iiorb=0
!!  do jproc=0,nproc-1
!!     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!     ilrold=-1
!!     do iorb=1,orbs%norb_par(jproc,0)
!!        ioverlaporb=0 ! counts the overlaps for the given orbital.
!!        iiorb=iiorb+1 ! counts the total orbitals
!!        ilr=orbs%inwhichlocreg(iiorb)
!!!        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!        do jorb=1,orbs%norb
!!           jlr=orbs%inwhichlocreg(jorb)
!!!           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!           call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzd%llr(jlr),isoverlap)
!!!           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
!!           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
!!!           if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!           if(isoverlap) then
!!              ioverlaporb=ioverlaporb+1
!!              op%overlaps(ioverlaporb,iiorb)=jorb
!!              if(ilr/=ilrold) then
!!                 ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
!!                 ! would get the same orbitals again. Therefore the counter is not increased
!!                 ! in that case.
!!                 ioverlapMPI=ioverlapMPI+1
!!                 !!comon%overlaps(ioverlapMPI,jproc)=jorb
!!              end if
!!           end if
!!        end do
!!        !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
!!        ilrold=ilr
!!     end do
!!     !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
!!  end do
!!
!!end subroutine determineOverlaps


!!! Count for each orbital and each process the number of overlapping orbitals.
!!subroutine countOverlaps(iproc, nproc, orbs, lzd, op, comon)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(overlapParameters),intent(out):: op
!!  type(p2pComms),intent(out):: comon
!!
!!  ! Local variables
!!  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, iiorb
!!   logical :: isoverlap
!!
!!  iiorb=0
!!  do jproc=0,nproc-1
!!     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!     ilrold=-1
!!     do iorb=1,orbs%norb_par(jproc,0)
!!        ioverlaporb=0 ! counts the overlaps for the given orbital.
!!        iiorb=iiorb+1 ! counts the total orbitals
!!        ilr=orbs%inwhichlocreg(iiorb)
!! !       call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!        do jorb=1,orbs%norb
!!           jlr=orbs%inwhichlocreg(jorb)
!!!           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!           call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzd%llr(jlr),isoverlap) 
!!!           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
!!           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
!!!           if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!            if(isoverlap) then
!!              ioverlaporb=ioverlaporb+1
!!              if(ilr/=ilrold) then
!!                 ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
!!                 ! would get the same orbitals again. Therefore the counter is not increased
!!                 ! in that case.
!!                 ioverlapMPI=ioverlapMPI+1
!!              end if
!!           end if
!!        end do
!!        op%noverlaps(iiorb)=ioverlaporb
!!        !!if(iproc==0) write(*,'(a,2i8)') 'iiorb, op%noverlaps(iiorb)', iiorb, op%noverlaps(iiorb)
!!        ilrold=ilr
!!     end do
!!     comon%noverlaps(jproc)=ioverlapMPI
!!     !if(iproc==0) write(*,'(a,2i8)') 'jproc, comon%noverlaps(jproc)', jproc, comon%noverlaps(jproc)
!!  end do
!!
!!end subroutine countOverlaps


subroutine set_comms_ortho(iproc, nproc, orbs, lzd, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(overlapParameters),intent(inout) :: op
  type(p2pComms),intent(out) :: comon

  ! Local variables
  integer :: jproc, iorb, jorb, iiorb, jjorb, mpisource, mpidest, istsource, istdest, ncount, istat, iall, ijorb
  integer :: ilr, ilrold, jprocold, ildim, ierr, isend, irecv, p2p_tag, tag
  integer,dimension(:),allocatable :: istsourceArr, istdestArr
  character(len=*),parameter :: subname='set_comms_ortho'
  logical,dimension(:),allocatable :: receivedOrbital

  allocate(istsourceArr(0:nproc-1), stat=istat)
  call memocc(istat, istsourceArr, 'istsourceArr',subname)
  allocate(istdestArr(0:nproc-1), stat=istat)
  call memocc(istat, istdestArr, 'istdestArr',subname)
  allocate(receivedOrbital(orbs%norb), stat=istat)
  call memocc(istat, receivedOrbital, 'receivedOrbital', subname)

  istsourceArr=1
  istdestArr=1

  comon%nsendBuf=0
  comon%nrecvBuf=0


  op%indexInRecvBuf=0
  op%ndim_lphiovrlp=0

  iiorb=0
  jprocold=-1
  do jproc=0,nproc-1
     receivedOrbital=.false.
     ijorb=0
     ilrold=-1
     do iorb=1,orbs%norb_par(jproc,0)
        iiorb=iiorb+1 
        ilr=orbs%inwhichlocreg(iiorb)
        ! Check whether process jproc has already received orbital jjorb.
        !if(iproc==0) write(*,'(a,5i8)') 'jproc, iorb, iiorb, ilr, ilrold', jproc, iorb, iiorb, ilr, ilrold
        if(ilr==ilrold) cycle
        ildim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
        do jorb=1,op%noverlaps(iiorb)
           jjorb=op%overlaps(jorb,iiorb)
           !write(*,'(a,7i8)') 'iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold', iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold
           ijorb=ijorb+1
           mpisource=orbs%onWhichMPI(jjorb)
           mpidest=jproc
           istsource=istsourceArr(mpisource)
           istdest=istdestArr(mpidest)
           if(iproc==jproc) then
              ncount=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
              !write(*,'(a,4i9)') 'iproc, iorb, jorb, ncount', iproc, iorb, jorb, ncount
           end if
           call mpi_bcast(ncount, 1, mpi_integer, jproc, mpi_comm_world, ierr)
           !tag=tag+1
           tag=p2p_tag(jproc)
           receivedOrbital(jjorb)=.true.
           call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comon%comarr(1,ijorb,jproc))
           !!comon%comarr(9,ijorb,jproc)=jjorb
           !!comon%comarr(10,ijorb,jproc)=iiorb
           !if(iproc==0) write(*,'(6(a,i0))') 'process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,', tag=',tag
           if(iproc==mpisource) then
              !write(*,'(4(a,i0))') 'adding ',ncount,' elements for orbital ',iiorb,' to nsendBuf, iproc=',iproc,', jproc=',jproc
              comon%nsendBuf=comon%nsendBuf+ncount
           end if
           if(iproc==mpidest) then
              !write(*,'(3(a,i0))') 'process ',iproc,' will get orbital ',jjorb,' at position ',istdest
              op%indexInRecvBuf(iorb,jjorb)=istdest
              comon%nrecvBuf=comon%nrecvBuf+ncount
              op%ndim_lphiovrlp=op%ndim_lphiovrlp+ildim
           end if
           istsourceArr(mpisource) = istsourceArr(mpisource) + ncount
           istdestArr(mpidest) = istdestArr(mpidest) + ncount
        end do
        ilrold=ilr
     end do
     jprocold=jproc
  end do

  iall = -product(shape(istsourceArr))*kind(istsourceArr)
  deallocate(istsourceArr, stat=istat)
  call memocc(istat, iall, 'istsourceArr',subname)
  iall = -product(shape(istdestArr))*kind(istdestArr)
  deallocate(istdestArr, stat=istat)
  call memocc(istat, iall, 'istdestArr',subname)
  iall = -product(shape(receivedOrbital))*kind(receivedOrbital)
  deallocate(receivedOrbital, stat=istat)
  call memocc(istat, iall, 'receivedOrbital',subname)


  ! Determine comon%nrecv and comon%nsend - maybe can be integrated in the above loop.
  ! First receives
  irecv=0
  do jproc=0,nproc-1
     do jorb=1,comon%noverlaps(jproc)
        mpidest=comon%comarr(4,jorb,jproc)
        if(iproc==mpidest) then
           irecv=irecv+1
        end if
     end do
  end do
  comon%nrecv=irecv

  ! Now the sends
  isend=0
  do jproc=0,nproc-1
     do jorb=1,comon%noverlaps(jproc)
        mpisource=comon%comarr(1,jorb,jproc)
        if(iproc==mpisource) then
           isend=isend+1
        end if
     end do
  end do
  comon%nsend=isend

  ! Allocate the requests for the point to point communication.
  allocate(comon%requests(max(comon%nsend,comon%nrecv),2), stat=istat)
  call memocc(istat, comon%requests, 'comon%requests', subname)

  ! To indicate that no communication is going on.
  comon%communication_complete=.true.
  comon%messages_posted=.false.


end subroutine set_comms_ortho



subroutine getMatrixElements2(iproc, nproc, lzd, orbs, op_lb, comon_lb, lphi, lhphi, mad, matrixElements)
use module_base
use module_types
use module_interfaces, exceptThisOne => getMatrixElements2
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc
type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(in) :: orbs
type(overlapParameters),intent(inout) :: op_lb
type(p2pComms),intent(inout) :: comon_lb
real(kind=8),dimension(orbs%npsidim_orbs),intent(in) :: lphi, lhphi
type(matrixDescriptors),intent(in) :: mad
real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: matrixElements

! Local variables
character(len=*),parameter :: subname='getMatrixElements2'
!! real(kind=8),dimension(:),allocatable :: lphiovrlp
!! real(kind=8) :: tt1, tt2, tt3
!! type(input_variables) :: input


  ! Put lphi in the sendbuffer,i.e. lphi will be sent to other processes' receive buffer.
  call extractOrbital3(iproc,nproc,orbs,orbs,orbs%npsidim_orbs,lzd,lzd,&
       op_lb,op_lb,lphi,comon_lb%nsendBuf,comon_lb%sendBuf)
  !!call postCommsOverlapNew(iproc,nproc,orbs,op_lb,lzd,lphi,comon_lb,tt1,tt2)
  call post_p2p_communication(iproc, nproc, comon_lb%nsendbuf, comon_lb%sendbuf, &
       comon_lb%nrecvbuf, comon_lb%recvbuf, comon_lb)
  !!call collectnew(iproc,nproc,comon_lb,mad,op_lb,orbs,lzd,comon_lb%nsendbuf,&
  !!     comon_lb%sendbuf,comon_lb%nrecvbuf,comon_lb%recvbuf,tt1,tt2,tt3)
  call wait_p2p_communication(iproc, nproc, comon_lb)
  ! Put lhphi to the sendbuffer,so we can the calculate <lphi|lhphi>
  call extractOrbital3(iproc,nproc,orbs,orbs,orbs%npsidim_orbs,lzd,lzd,&
       op_lb,op_lb,lhphi,comon_lb%nsendBuf,comon_lb%sendBuf)
  call calculateOverlapMatrix3(iproc,nproc,orbs,op_lb,comon_lb%nsendBuf,&
                               comon_lb%sendBuf,comon_lb%nrecvBuf,comon_lb%recvBuf,mad,matrixElements)

end subroutine getMatrixElements2



subroutine getOverlapMatrix2(iproc, nproc, lzd, orbs, comon_lb, op_lb, lphi, mad, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => getOverlapMatrix2
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  type(p2pComms),intent(inout) :: comon_lb
  type(overlapParameters),intent(inout) :: op_lb
  real(kind=8),dimension(orbs%npsidim_orbs),intent(inout) :: lphi
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: ovrlp

  ! Local variables
  character(len=*),parameter :: subname='getOverlapMatrix2'
!!  real(kind=8) :: tt1, tt2, tt3
  call allocateCommuncationBuffersOrtho(comon_lb, subname)
  call extractOrbital3(iproc,nproc,orbs,orbs,orbs%npsidim_orbs,&
       lzd,lzd,op_lb,op_lb,lphi,comon_lb%nsendBuf,comon_lb%sendBuf)
  call post_p2p_communication(iproc, nproc, comon_lb%nsendbuf, comon_lb%sendbuf, &
       comon_lb%nrecvbuf, comon_lb%recvbuf, comon_lb)
  call wait_p2p_communication(iproc, nproc, comon_lb)
  call calculateOverlapMatrix3(iproc, nproc, orbs, op_lb, comon_lb%nsendBuf, &
       comon_lb%sendBuf, comon_lb%nrecvBuf, comon_lb%recvBuf, mad, ovrlp)
  call deallocateCommuncationBuffersOrtho(comon_lb, subname)


end subroutine getOverlapMatrix2



subroutine extractOrbital3(iproc, nproc, orbs, orbsig, sizePhi, lzd, lzdig, op, opig, phi, nsendBuf, sendBuf)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, sizePhi
  type(orbitals_data),intent(in) :: orbs, orbsig
  type(local_zone_descriptors),intent(in) :: lzd, lzdig
  type(overlapParameters),intent(inout) :: op, opig
  real(kind=8),dimension(sizePhi),intent(in) :: phi
  integer,intent(in) :: nsendBuf
  real(kind=8),dimension(nsendBuf),intent(out) :: sendBuf


  ! Local variables
  integer :: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
  integer :: i, jst, istart, iend, ncount, iseg, kstart, kend, start, kold, kseg, knseg
!! integer :: indSource

  call timing(iproc,'extract_orbs  ','ON')

  indovrlp=1
  op%indexInSendBuf=0
  sendbuf = 0.0_dp
  ilrold=-1
  iiprocold=-1
  do iorb=1,orbs%norb
     ilr=orbs%inwhichlocreg(iorb)
     iiproc=orbs%onWhichMPI(iorb)
     if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
     do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=orbsig%inwhichlocreg(jjorb)
        jjproc=orbsig%onWhichMPI(jjorb)
        if(iproc==jjproc) then
           ! Get the correct descriptors
           korb=jjorb-orbsig%isorb
           !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
           do i=1,opig%noverlaps(jjorb)
              !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
              if(opig%overlaps(i,jjorb)==iorb) then
                 lorb=i
                 exit
              end if
           end do
           !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
           gdim=lzdig%llr(jjlr)%wfd%nvctr_c+7*lzdig%llr(jjlr)%wfd%nvctr_f
           ldim=opig%wfd_overlap(lorb,korb)%nvctr_c+7*opig%wfd_overlap(lorb,korb)%nvctr_f
           ind=1
           do kkorb=orbsig%isorb+1,jjorb-1
              klr=orbs%inwhichlocreg(kkorb)
              ind = ind + lzdig%llr(klr)%wfd%nvctr_c + 7*lzdig%llr(klr)%wfd%nvctr_f
           end do
           !! THIS IS THE OLD VERSION
           !!do i=0,ldim-1
           !!    indSource=ind+op%indexExtract(indovrlp+i)-1
           !!    sendBuf(indovrlp+i)=phi(indSource)
           !!end do

            !! THIS IS THE NEWEST VERSION (NOT OPTIMIZED, could probably store the keyv once and for all)
            jst=0
            kold = 1
            knseg = 0
            do iseg=1,opig%wfd_overlap(lorb,korb)%nseg_c
                istart=opig%wfd_overlap(lorb,korb)%keyglob(1,iseg)
                iend=opig%wfd_overlap(lorb,korb)%keyglob(2,iseg)
                ncount=iend-istart+1  
                inner_loop: do kseg=kold,lzdig%llr(jjlr)%wfd%nseg_c
                   kstart = lzdig%llr(jjlr)%wfd%keyglob(1,kseg)
                   kend   = lzdig%llr(jjlr)%wfd%keyglob(2,kseg)
                   if(kstart <= iend .and. kend >= istart) then 
                      knseg = knseg + 1 
                      kold = kseg + 1
                      start = lzdig%llr(jjlr)%wfd%keyvglob(kseg) + max(0,istart-kstart)
                      call dcopy(ncount, phi(ind+start-1), 1, sendBuf(indovrlp+jst), 1)
                      jst=jst+ncount
                      exit inner_loop
                   end if
                end do inner_loop
            end do
            if(jst .ne. opig%wfd_overlap(lorb,korb)%nvctr_c) then
               print *,'ilr,jjlr',ilr,jjlr
               print *,'knseg: ',knseg,'onseg:',opig%wfd_overlap(lorb,korb)%nseg_c
               print *, 'extractOrbital3: jst = ',jst,'not equal to ldim = ',opig%wfd_overlap(lorb,korb)%nvctr_c
               stop
            end if
            !MUST DO THE FINE GRID NOW
            jst=0
            kold = 1
            knseg = 0
            do iseg=1,opig%wfd_overlap(lorb,korb)%nseg_f
                istart=opig%wfd_overlap(lorb,korb)%keyglob(1,iseg+opig%wfd_overlap(lorb,korb)%nseg_c)
                iend=opig%wfd_overlap(lorb,korb)%keyglob(2,iseg+opig%wfd_overlap(lorb,korb)%nseg_c)
                ncount=7*(iend-istart+1) 
                inner_loop2: do kseg=kold,lzdig%llr(jjlr)%wfd%nseg_f
                   kstart = lzdig%llr(jjlr)%wfd%keyglob(1,kseg+lzdig%llr(jjlr)%wfd%nseg_c)
                   kend   = lzdig%llr(jjlr)%wfd%keyglob(2,kseg+lzdig%llr(jjlr)%wfd%nseg_c)
                   if(kstart <= iend .and. kend >= istart) then 
                      knseg = knseg + 1 
                      kold = kseg + 1
                      start = lzdig%llr(jjlr)%wfd%nvctr_c+(lzdig%llr(jjlr)%wfd%keyvglob(kseg+lzdig%llr(jjlr)%wfd%nseg_c) +&
                              max(0,istart-kstart)-1)*7
                      call dcopy(ncount,phi(ind+start),1,&
                              sendBuf(indovrlp+jst+opig%wfd_overlap(lorb,korb)%nvctr_c),1)
                      jst=jst+ncount
                      exit inner_loop2
                   end if
                end do inner_loop2
            end do
            if(jst .ne. 7*opig%wfd_overlap(lorb,korb)%nvctr_f) then
               print *,'ilr,jjlr',ilr,jjlr
               print *,'knseg: ',knseg,'onseg:',opig%wfd_overlap(lorb,korb)%nseg_f
               print *, 'extractOrbital3: jst = ',jst,'not equal to ldim = ',7*opig%wfd_overlap(lorb,korb)%nvctr_f
               stop
            end if
            
            !!!! THIS IS THE NEW VERSION
            !!m=mod(ldim,7)
            !!if(m/=0) then
            !!    do i=0,m-1
            !!        ii=indovrlp+i
            !!        indSource=ind+op%indexExtract(ii)-1
            !!        sendBuf(ii)=phi(indSource)
            !!    end do
            !!end if
            !!do i=m,ldim-1,7
            !!    ii=indovrlp+i
            !!    indSource=ind+op%indexExtract(ii+0)-1
            !!    sendBuf(ii+0)=phi(indSource)
            !!    indSource=ind+op%indexExtract(ii+1)-1
            !!    sendBuf(ii+1)=phi(indSource)
            !!    indSource=ind+op%indexExtract(ii+2)-1
            !!    sendBuf(ii+2)=phi(indSource)
            !!    indSource=ind+op%indexExtract(ii+3)-1
            !!    sendBuf(ii+3)=phi(indSource)
            !!    indSource=ind+op%indexExtract(ii+4)-1
            !!    sendBuf(ii+4)=phi(indSource)
            !!    indSource=ind+op%indexExtract(ii+5)-1
            !!    sendBuf(ii+5)=phi(indSource)
            !!    indSource=ind+op%indexExtract(ii+6)-1
            !!    sendBuf(ii+6)=phi(indSource)
            !!end do
           if(orbs%norb==orbsig%norb) then
               op%indexInSendBuf(jjorb-orbsig%isorb,iorb)=indovrlp
           else if(orbsig%norb>orbs%norb) then
               opig%indexInSendBuf(jjorb-orbsig%isorb,iorb)=indovrlp
           else
               stop 'orbsig%norb<orbs%norb not yet implemented'
           end if
           indovrlp=indovrlp+opig%wfd_overlap(lorb,korb)%nvctr_c+7*opig%wfd_overlap(lorb,korb)%nvctr_f
        end if
     end do
     ilrold=ilr
     iiprocold=iiproc

  end do

  if(indovrlp/=nsendBuf+1) then
     write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=nsendBuf+1', indovrlp, nsendBuf+1
     stop
  end if

  call timing(iproc,'extract_orbs  ','OF')

end subroutine extractOrbital3



subroutine calculateOverlapMatrix3(iproc, nproc, orbs, op, nsendBuf, sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsendBuf, nrecvBuf
  type(orbitals_data),intent(in) :: orbs
  type(overlapParameters),intent(in) :: op
  real(kind=8),dimension(nsendBuf),intent(in) :: sendBuf
  real(kind=8),dimension(nrecvBuf),intent(in) :: recvBuf
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: ovrlp

  ! Local variables
  integer :: iorb, jorb, iiorb, jjorb, ist, jst, ncount, ierr, istat, iall
  real(kind=8) :: ddot ! , tt, ttmax
  real(kind=8),dimension(:),allocatable :: ovrlpCompressed_send, ovrlpCompressed_receive
  character(len=*),parameter :: subname='calculateOverlapMatrix3'
  integer,dimension(:),allocatable :: sendcounts, displs


  call timing(iproc,'lovrlp_comp   ','ON')

  !!ovrlp=0.d0
  call to_zero(orbs%norb**2, ovrlp(1,1))

  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        call getStartingIndices(iorb, jorb, op, orbs, ist, jst)
        ncount=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
        ovrlp(jjorb,iiorb)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
        !!write(300+iproc,*) iiorb,jjorb,ovrlp(jjorb,iiorb)
     end do
  end do

  call timing(iproc,'lovrlp_comp   ','OF')


  allocate(ovrlpCompressed_send(mad%nvctr), stat=istat)
  call memocc(istat, ovrlpCompressed_send, 'ovrlpCompressed_send', subname)

  allocate(sendcounts(0:nproc-1), stat=istat)
  call memocc(istat, sendcounts, 'sendcounts', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  call timing(iproc,'lovrlp_compr  ','ON')
  call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed_send, sendcounts, displs)
  call timing(iproc,'lovrlp_compr  ','OF')
  allocate(ovrlpCompressed_receive(mad%nvctr), stat=istat)
  call memocc(istat, ovrlpCompressed_receive, 'ovrlpCompressed_receive', subname)
  call timing(iproc,'lovrlp_comm   ','ON')
  !!do istat=1,mad%nvctr
  !!    write(400+iproc,*) istat, ovrlpCompressed_send(istat)
  !!end do
  !!write(*,'(a,i5,2x,100i5)') 'iproc, sendcounts', iproc, sendcounts
  !!write(*,'(a,i5,2x,100i5)') 'iproc, displs', iproc, displs
  if (nproc >1) then
     call mpi_allgatherv(ovrlpCompressed_send(displs(iproc)+1), sendcounts(iproc),&
          mpi_double_precision, ovrlpCompressed_receive(1), &
          sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
  else
     call vcopy(sendcounts(iproc),ovrlpCompressed_send(displs(iproc)+1),1,ovrlpCompressed_receive(1+displs(iproc)),1)
  end if
  !!do istat=1,mad%nvctr
  !!    write(500+iproc,*) istat, ovrlpCompressed_receive(istat)
  !!end do
  call timing(iproc,'lovrlp_comm   ','OF')

  call timing(iproc,'lovrlp_uncompr','ON')
  call uncompress_matrix(orbs%norb, mad, ovrlpCompressed_receive, ovrlp)
  call timing(iproc,'lovrlp_uncompr','OF')

  iall=-product(shape(ovrlpCompressed_send))*kind(ovrlpCompressed_send)
  deallocate(ovrlpCompressed_send, stat=istat)
  call memocc(istat, iall, 'ovrlpCompressed_send', subname)
  iall=-product(shape(ovrlpCompressed_receive))*kind(ovrlpCompressed_receive)
  deallocate(ovrlpCompressed_receive, stat=istat)
  call memocc(istat, iall, 'ovrlpCompressed_receive', subname)
  iall=-product(shape(sendcounts))*kind(sendcounts)
  deallocate(sendcounts, stat=istat)
  call memocc(istat, iall, 'sendcounts', subname)
  iall=-product(shape(displs))*kind(displs)
  deallocate(displs, stat=istat)
  call memocc(istat, iall, 'displs', subname)

  !ttmax=0.d0
  !do iorb=1,orbs%norb
  !    do jorb=iorb,orbs%norb
  !        tt=abs(ovrlp(jorb,iorb)-ovrlp(iorb,jorb))
  !        if(tt>ttmax) ttmax=tt
  !    end do
  !end do
  !if(iproc==0) write(*,*) 'in calculateOverlapMatrix3: max dev from symmetry:', ttmax



end subroutine calculateOverlapMatrix3



subroutine calculateOverlapMatrix3Partial(iproc, nproc, orbs, op, nsendBuf, &
     sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsendBuf, nrecvBuf
  type(orbitals_data),intent(in) :: orbs
  type(overlapParameters),intent(in) :: op
  real(kind=8),dimension(nsendBuf),intent(in) :: sendBuf
  real(kind=8),dimension(nrecvBuf),intent(in) :: recvBuf
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: ovrlp

  ! Local variables
  integer :: iorb, jorb, iiorb, jjorb, ist, jst, ncount
  real(kind=8) :: ddot
  character(len=*),parameter :: subname='calculateOverlapMatrix3'


  !!ovrlp=0.d0
  call to_zero(orbs%norb**2, ovrlp(1,1))

  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        call getStartingIndices(iorb, jorb, op, orbs, ist, jst)
        ncount=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
        ovrlp(jjorb,iiorb)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
     end do
  end do

end subroutine calculateOverlapMatrix3Partial



subroutine globalLoewdin(iproc, nproc, orbs, lzd, op, comon, ovrlp, lphiovrlp, lphi)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(overlapParameters),intent(in) :: op
  type(p2pComms),intent(in) :: comon
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: ovrlp
  real(kind=8),dimension(op%ndim_lphiovrlp),intent(in) :: lphiovrlp
  real(kind=8),dimension(orbs%npsidim_orbs),intent(out) :: lphi

  ! Local variables
  integer :: iorb, iiorb, ilr, ist, ncount
  real(kind=8) :: tt, dnrm2


  call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, comon%recvbuf, ovrlp, .true., lphi)

  ! Normalize
  ist=1
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr=orbs%inwhichlocreg(iiorb)
     ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f

     ! Normalize
     tt=dnrm2(ncount, lphi(ist), 1)
     call dscal(ncount, 1/tt, lphi(ist), 1)

     ist=ist+ncount
  end do


end subroutine globalLoewdin



subroutine getStartingIndices(iorb, jorb, op, orbs, ist, jst)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iorb, jorb
  type(overlapParameters),intent(in) :: op
  type(orbitals_data),intent(in) :: orbs
  integer,intent(out) :: ist, jst
  
  ! Local variables
  integer :: jjlr, jjproc, jj, iilr, iiproc, korb, kkorb, ii, iiorb, jjorb


  iiorb=orbs%isorb+iorb
  jjorb=op%overlaps(jorb,iiorb)
  ! Starting index of orbital iorb, already transformed to overlap region with jjorb.
  ! We have to find the first orbital on the same MPI and in the same locreg as jjorb.
  jjlr=orbs%inWhichLocreg(jjorb)
  jjproc=orbs%onWhichMPI(jjorb)
  jj=orbs%isorb_par(jjproc)+1
  do
      if(orbs%inWhichLocreg(jj)==jjlr) exit
      jj=jj+1
  end do
  ist=op%indexInSendBuf(iorb,jj)
  ! Starting index of orbital jjorb.
  ! We have to find the first orbital on the same MPI and in the same locreg as iiorb.
  iilr=orbs%inWhichLocreg(iiorb)
  iiproc=orbs%onWhichMPI(iiorb)
  do korb=1,orbs%norbp
      kkorb=orbs%isorb_par(iiproc)+korb
      if(orbs%inWhichLocreg(kkorb)==iilr) then
          ii=korb
          exit
      end if
  end do
  jst=op%indexInRecvBuf(ii,jjorb)


end subroutine getStartingIndices


!> Returns the starting and ending indices (on the coarse grid) of a given localization region.
subroutine getIndices(lr, is1, ie1, is2, ie2, is3, ie3)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(in) :: lr
  integer,intent(out) :: is1, ie1, is2, ie2, is3, ie3

  is1=lr%ns1!+1 !PB: should not have a +1
  ie1=lr%ns1+lr%d%n1
  is2=lr%ns2!+1
  ie2=lr%ns2+lr%d%n2
  is3=lr%ns3!+1
  ie3=lr%ns3+lr%d%n3

end subroutine getIndices





subroutine applyOrthoconstraintNonorthogonal2(iproc, nproc, methTransformOverlap, blocksize_pdgemm, &
           correction_orthoconstraint, orbs, lagmat, ovrlp, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => applyOrthoconstraintNonorthogonal2
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, methTransformOverlap, blocksize_pdgemm, correction_orthoconstraint
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: ovrlp
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: lagmat
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans

  ! Local variables
  integer :: iorb, jorb, istat, iall
  real(kind=8),dimension(:,:),allocatable :: ovrlp2, lagmat_trans
  character(len=*),parameter :: subname='applyOrthoconstraintNonorthogonal2'

  call timing(iproc,'lagmat_orthoco','ON')

  allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp2, 'ovrlp2', subname)

  correctionIf: if(correction_orthoconstraint==0) then
  
    !call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)

    allocate(lagmat_trans(orbs%norb,orbs%norb), stat=istat)
    call memocc(istat, lagmat_trans, 'lagmat_trans', subname)

    call dcopy(orbs%norb**2, lagmat(1,1), 1, lagmat_trans(1,1), 1)
  
    ! Invert the overlap matrix
    call timing(iproc,'lagmat_orthoco','OF')
    call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, blocksize_pdgemm, orbs%norb, ovrlp, ovrlp2)
    call timing(iproc,'lagmat_orthoco','ON')
  
  
    ! Multiply the Lagrange multiplier matrix with S^-1
    ! First fill the upper triangle.
    do iorb=1,orbs%norb
       do jorb=1,iorb-1
          ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
       end do
    end do
    !if(blocksize_pdgemm<0) then
       !!ovrlp_minus_one_lagmat=0.d0
       call to_zero(orbs%norb**2, ovrlp_minus_one_lagmat(1,1))
       call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
            0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
  
       ! Transpose lagmat
       do iorb=1,orbs%norb
           do jorb=iorb+1,orbs%norb
               lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
               lagmat_trans(iorb,jorb)=lagmat(jorb,iorb)
           end do
       end do
       call to_zero(orbs%norb**2, ovrlp_minus_one_lagmat_trans(1,1))
       call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat_trans(1,1), orbs%norb, &
            0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
    ! THIS SHOULD BE DGEMM NOT DSYMM_PARALLEL
    !else
    !  call dsymm_parallel(iproc, nproc, blocksize_pdgemm, bigdft_mpi%mpi_comm, 'l', 'l', orbs%norb, orbs%norb, 1.d0, &
    !       ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, 0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
    !  ! Transpose lagmat
    !  do iorb=1,orbs%norb
    !      do jorb=iorb+1,orbs%norb
    !          lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
    !          lagmat_trans(iorb,jorb)=lagmat(jorb,iorb)
    !      end do
    !  end do
    !  call dsymm_parallel(iproc, nproc, blocksize_pdgemm, bigdft_mpi%mpi_comm, 'l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), &
    !       orbs%norb, lagmat_trans(1,1), orbs%norb, &
    !       0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
    !end if

    iall=-product(shape(lagmat_trans))*kind(lagmat_trans)
    deallocate(lagmat_trans, stat=istat)
    call memocc(istat, iall, 'lagmat_trans', subname)
  
  !!else if(correction_orthoconstraint==1) then correctionIf
  !!    do iorb=1,orbs%norb
  !!        do jorb=1,orbs%norb
  !!            ovrlp_minus_one_lagmat(jorb,iorb)=lagmat(jorb,iorb)
  !!            ovrlp_minus_one_lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
  !!        end do
  !!    end do
  end if correctionIf
  
  call timing(iproc,'lagmat_orthoco','OF')
  
  
  iall=-product(shape(ovrlp2))*kind(ovrlp2)
  deallocate(ovrlp2, stat=istat)
  call memocc(istat, iall, 'ovrlp2', subname)


end subroutine applyOrthoconstraintNonorthogonal2

subroutine overlapPowerPlusMinusOneHalf(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, &
           blocksize_pdgemm, norb, ovrlp, inv_ovrlp_half, plusminus)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb
  type(sparse_matrix),intent(inout) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half
  logical, intent(in) :: plusminus ! if true S^1/2, if false S^-1/2

  ! Local variables
  integer :: lwork, istat, iall, iorb, jorb, info, iiorb, jjorb, ii, ii_inv!, iseg
  integer :: matrixindex_in_compressed
  character(len=*),parameter :: subname='overlapPowerMinusOneHalf'
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:,:),allocatable :: tempArr
  !*real(8),dimension(:,:), allocatable :: vr,vl ! for non-symmetric LAPACK
  !*real(8),dimension(:),allocatable:: eval1 ! for non-symmetric LAPACK
  !*real(dp) :: temp
  !*real(dp), allocatable, dimension(:) :: temp_vec

  call timing(iproc,'lovrlp^-1/2   ','ON')

  if(methTransformOrder==0) then

      ! Exact calculation of ovrlp**(-1/2)
      allocate(eval(norb), stat=istat)
      call memocc(istat, eval, 'eval', subname)
      allocate(tempArr(norb,norb,2), stat=istat)
      call memocc(istat, tempArr, 'tempArr', subname)

      allocate(ovrlp%matrix(norb,norb), stat=istat)
      call memocc(istat, ovrlp%matrix, 'ovrlp%matrix', subname)

      call uncompress_matrix(iproc,ovrlp)
  
      allocate(inv_ovrlp_half%matrix(norb,norb), stat=istat)
      call memocc(istat, inv_ovrlp_half%matrix, 'inv_ovrlp_half%matrix', subname)

      call vcopy(norb*norb,ovrlp%matrix(1,1),1,inv_ovrlp_half%matrix(1,1),1)

      iall=-product(shape(ovrlp%matrix))*kind(ovrlp%matrix)
      deallocate(ovrlp%matrix, stat=istat)
      call memocc(istat, iall, 'ovrlp%matrix', subname)
    
      if(blocksize_dsyev>0) then
          call dsyev_parallel(iproc, nproc, min(blocksize_dsyev,norb), comm, 'v', 'l', norb, inv_ovrlp_half%matrix(1,1), &
               norb, eval(1), info)
          if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
              !stop
          end if
      else
          
          allocate(work(100), stat=istat)
          call dsyev('v', 'l', norb, inv_ovrlp_half%matrix(1,1), norb, eval, work, -1, info)
          lwork = int(work(1))
          deallocate(work, stat=istat)
          allocate(work(lwork), stat=istat)
          call memocc(istat, work, 'work', subname)
          call dsyev('v', 'l', norb, inv_ovrlp_half%matrix(1,1), norb, eval, work, lwork, info)

          !*!lwork=1000*norb
          !*allocate(work(1), stat=istat)
          !*call DGEEV( 'v','v', norb, inv_ovrlp_half%matrix(1,1), norb, eval, eval1, VL, norb, VR,&
          !*     norb, WORK, -1, info )
          !*lwork = nint(work(1))
          !*deallocate(work, stat=istat)
          !*allocate(work(lwork), stat=istat)
          !*call memocc(istat, work, 'work', subname)
          !*! lr408 - see if LAPACK is stil to blame for convergence issues
          !*allocate(vl(1:norb,1:norb))
          !*allocate(vr(1:norb,1:norb))
          !*allocate(eval1(1:norb))
          !*call DGEEV( 'v','v', norb, inv_ovrlp_half%matrix(1,1), norb, eval, eval1, VL, norb, VR,&
          !*     norb, WORK, LWORK, info )
          !*inv_ovrlp_half%matrix=vl
          !*write(14,*) eval1
          !*deallocate(eval1)
          !*deallocate(vr)
          !*deallocate(vl)
          !*allocate(temp_vec(1:norb))
          !*do iorb=1,norb
          !*   do jorb=iorb+1,norb
          !*      if (eval(jorb) < eval(iorb)) then
          !*         temp = eval(iorb)
          !*         temp_vec = inv_ovrlp_half%matrix(:,iorb)
          !*         eval(iorb) = eval(jorb)
          !*         eval(jorb) = temp
          !*         inv_ovrlp_half%matrix(:,iorb) = inv_ovrlp_half%matrix(:,jorb)
          !*         inv_ovrlp_half%matrix(:,jorb) = temp_vec
          !*      end if
          !*   end do
          !*end do
          !*deallocate(temp_vec)

          !  lr408 - see if LAPACK is stil to blame for convergence issues
          if(info/=0) then
              write(*,'(a,i0,2x,i0)') 'ERROR in dsyev (overlapPowerPlusMinusOneHalf), info, norb=', info, norb
              stop
          end if
          iall=-product(shape(work))*kind(work)
          deallocate(work, stat=istat)
          call memocc(istat, iall, 'work', subname)
      end if

      ! Calculate S^{-1/2}. 
      ! First calculate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
      ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
      do iorb=1,norb
          do jorb=1,norb
              if (plusminus) then
                  tempArr(jorb,iorb,1)=inv_ovrlp_half%matrix(jorb,iorb)*sqrt(abs(eval(iorb)))
              else
                  tempArr(jorb,iorb,1)=inv_ovrlp_half%matrix(jorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
              end if
          end do
      end do
      
      ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
      ! This will give S^{-1/2}.
      if(blocksize_pdgemm<0) then
          call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half%matrix(1,1), &
               norb, tempArr(1,1,1), norb, 0.d0, tempArr(1,1,2), norb)
      else
          call dgemm_parallel(iproc, nproc, blocksize_pdgemm, comm, 'n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half%matrix(1,1), &
               norb, tempArr(1,1,1), norb, 0.d0, tempArr(1,1,2), norb)
      end if
      call dcopy(norb**2, tempArr(1,1,2), 1, inv_ovrlp_half%matrix(1,1), 1)

      call compress_matrix(iproc,inv_ovrlp_half)

      iall=-product(shape(eval))*kind(eval)
      deallocate(eval, stat=istat)
      call memocc(istat, iall, 'eval', subname)
      iall=-product(shape(tempArr))*kind(tempArr)
      deallocate(tempArr, stat=istat)
      call memocc(istat, iall, 'tempArr', subname)

      iall=-product(shape(inv_ovrlp_half%matrix))*kind(inv_ovrlp_half%matrix)
      deallocate(inv_ovrlp_half%matrix, stat=istat)
      call memocc(istat, iall, 'inv_ovrlp_half%matrix', subname)

  else if(methTransformOrder==1) then
      ! inv_ovrlp_half can be less sparse than ovrlp, so pad with zeros first
      call to_zero(inv_ovrlp_half%nvctr,inv_ovrlp_half%matrix_compr(1))
      ! Taylor expansion up to first order.
      if (plusminus) then
           do ii=1,ovrlp%nvctr
              iiorb = ovrlp%orb_from_index(1,ii)
              jjorb = ovrlp%orb_from_index(2,ii)

              ii_inv = matrixindex_in_compressed(inv_ovrlp_half,iiorb,jjorb) ! double check this order
              if(iiorb==jjorb) then
                  inv_ovrlp_half%matrix_compr(ii_inv)=0.5d0+0.5d0*ovrlp%matrix_compr(ii)
              else
                  inv_ovrlp_half%matrix_compr(ii_inv)=0.5d0*ovrlp%matrix_compr(ii)
              end if
          end do
      else
           do ii=1,ovrlp%nvctr
              iiorb = ovrlp%orb_from_index(1,ii)
              jjorb = ovrlp%orb_from_index(2,ii)

              ii_inv = matrixindex_in_compressed(inv_ovrlp_half,iiorb,jjorb) ! double check this order
              if(iiorb==jjorb) then
                  inv_ovrlp_half%matrix_compr(ii_inv)=1.5d0-.5d0*ovrlp%matrix_compr(ii)
              else
                  inv_ovrlp_half%matrix_compr(ii_inv)=-.5d0*ovrlp%matrix_compr(ii)
              end if
          end do
      end if
  else
      
      stop 'deprecated'

  end if

  call timing(iproc,'lovrlp^-1/2   ','OF')

end subroutine overlapPowerPlusMinusOneHalf




! Should be used if sparsemat is not available... to be cleaned
subroutine overlapPowerPlusMinusOneHalf_old(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, &
           blocksize_pdgemm, norb, ovrlp, inv_ovrlp_half, plusminus, orbs)
  use module_base
  use module_types
  use module_interfaces, fake_name => overlapPowerPlusMinusOneHalf_old
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, norb, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp
  real(kind=8),dimension(norb,norb),intent(inout) :: inv_ovrlp_half
  type(orbitals_data), optional, intent(in) :: orbs
  logical, intent(in) :: plusminus ! if true, calculates S^1/2, if false S^-1/2
  
  ! Local variables
  integer :: lwork, istat, iall, iorb, jorb, info, iiorb, ierr
  character(len=*),parameter :: subname='overlapPowerMinusOneHalf'
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_halfp
  real(kind=8),dimension(:,:,:),allocatable :: tempArr

  call timing(iproc,'lovrlp^-1/2old','ON')

  ! this first option is called from coeff_weight_analysis and for constrained DFT, needs optimizing!
  if(methTransformOrder==0) then
  
      call vcopy(norb*norb,ovrlp(1,1),1,inv_ovrlp_half(1,1),1)

      ! Exact calculation of ovrlp**(-1/2)
      allocate(eval(norb), stat=istat)
      call memocc(istat, eval, 'eval', subname)
      allocate(tempArr(norb,norb,2), stat=istat)
      call memocc(istat, tempArr, 'tempArr', subname)
      
      
      if(blocksize_dsyev>0) then
          call dsyev_parallel(iproc, nproc, min(blocksize_dsyev,norb), comm, 'v', 'l', norb, inv_ovrlp_half(1,1), &
               norb, eval(1), info)
          if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
              !stop
          end if
      else
          
          allocate(work(10000), stat=istat)
          call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, -1, info)

          lwork = nint(work(1))

          deallocate(work, stat=istat)
          allocate(work(lwork), stat=istat)
          call memocc(istat, work, 'work', subname)
          call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, lwork, info)

          !*allocate(work(1), stat=istat)
          !*call DGEEV( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
          !*     norb, WORK, -1, info )
          !*lwork = work(1)
          !*deallocate(work, stat=istat)
          !*allocate(work(lwork), stat=istat)
          !*call memocc(istat, work, 'work', subname)
          !*! lr408 - see if LAPACK is stil to blame for convergence issues
          !*allocate(vl(1:norb,1:norb))
          !*allocate(vr(1:norb,1:norb))
          !*allocate(eval1(1:norb))
          !*call DGEEV( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
          !*     norb, WORK, LWORK, info )
          !*inv_ovrlp_half=vl
          !*write(14,*) eval1
          !*deallocate(eval1)
          !*deallocate(vr)
          !*deallocate(vl)
          !*allocate(temp_vec(1:norb))
          !*do iorb=1,norb
          !*   do jorb=iorb+1,norb
          !*      if (eval(jorb) < eval(iorb)) then
          !*         temp = eval(iorb)
          !*         temp_vec = inv_ovrlp_half(:,iorb)
          !*         eval(iorb) = eval(jorb)
          !*         eval(jorb) = temp
          !*         inv_ovrlp_half(:,iorb) = inv_ovrlp_half(:,jorb)
          !*         inv_ovrlp_half(:,jorb) = temp_vec
          !*      end if
          !*   end do
          !*end do
          !*deallocate(temp_vec)

          !  lr408 - see if LAPACK is stil to blame for convergence issues
          if(info/=0) then
              write(*,'(a,i0,2x,i0)') 'ERROR in dsyev (overlapPowerPlusMinusOneHalf_old), info, norb=', info, norb
              stop
          end if
          iall=-product(shape(work))*kind(work)
          deallocate(work, stat=istat)
          call memocc(istat, iall, 'work', subname)
      end if

      ! Calculate S^{-1/2}. 
      ! First calculate ovrlp*diag(1/sqrt(eval)) (ovrlp is the diagonalized overlap
      ! matrix and diag(1/sqrt(eval)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
      !$omp parallel do default(private) shared(tempArr,inv_ovrlp_half,norb,plusminus,eval)
      do iorb=1,norb
          do jorb=1,norb
              if (plusminus) then
                 tempArr(jorb,iorb,1)=inv_ovrlp_half(jorb,iorb)*sqrt(abs(eval(iorb)))
              else
                 tempArr(jorb,iorb,1)=inv_ovrlp_half(jorb,iorb)/sqrt(abs(eval(iorb)))
              end if
          end do
      end do
      !$omp end parallel do
      
      ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
      ! This will give S^{-1/2}.
      if(blocksize_pdgemm<0) then
          call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half(1,1), &
               norb, tempArr(1,1,1), norb, 0.d0, tempArr(1,1,2), norb)
      else
          call dgemm_parallel(iproc, nproc, blocksize_pdgemm, comm, 'n', 't', norb, norb, norb, 1.d0, &
               inv_ovrlp_half(1,1), norb, tempArr(1,1,1), norb, 0.d0, tempArr(1,1,2), norb)
      end if
      call dcopy(norb**2, tempArr(1,1,2), 1, inv_ovrlp_half(1,1), 1)

      iall=-product(shape(eval))*kind(eval)
      deallocate(eval, stat=istat)
      call memocc(istat, iall, 'eval', subname)
      iall=-product(shape(tempArr))*kind(tempArr)
      deallocate(tempArr, stat=istat)
      call memocc(istat, iall, 'tempArr', subname)

  else if(methTransformOrder==1) then
     if (present(orbs)) then ! parallel version

        !if (norb/=orbs%norb) add warning later

        ! Taylor expansion up to first order.
        allocate(inv_ovrlp_halfp(orbs%norb,orbs%norbp), stat=istat)
        call memocc(istat, inv_ovrlp_halfp, 'inv_ovrlp_halfp', subname)

        ! No matrix compression available
        if (plusminus) then
           !!$omp parallel do default(private) shared(inv_ovrlp_half,ovrlp,orbs)
           do iorb=1,orbs%norbp
              iiorb=orbs%isorb+iorb
              do jorb=1,orbs%norb
                 if(iiorb==jorb) then
                    inv_ovrlp_halfp(jorb,iorb)=0.5d0+0.5d0*ovrlp(jorb,iiorb)
                 else
                    inv_ovrlp_halfp(jorb,iorb)=0.5d0*ovrlp(jorb,iiorb)
                 end if
              end do
           end do
           !!$omp end parallel do
        else
            !!$omp parallel do default(private) shared(inv_ovrlp_half,ovrlp,orbs)
           do iorb=1,orbs%norbp
              iiorb=orbs%isorb+iorb
              do jorb=1,orbs%norb
                 if(iiorb==jorb) then
                    inv_ovrlp_halfp(jorb,iorb)=1.5d0-.5d0*ovrlp(jorb,iiorb)
                 else
                    inv_ovrlp_halfp(jorb,iorb)=-.5d0*ovrlp(jorb,iiorb)
                 end if
              end do
           end do
           !!$omp end parallel do
        end if


        call timing(iproc,'lovrlp^-1/2old','OF')
        call timing(iproc,'lovrlp^-1/2com','ON')
        ! gather together
        if(nproc > 1) then
           call mpi_allgatherv(inv_ovrlp_halfp, orbs%norb*orbs%norbp, mpi_double_precision, inv_ovrlp_half, &
                orbs%norb*orbs%norb_par(:,0), orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
        else
           call dcopy(orbs%norbp*orbs%norb,inv_ovrlp_halfp(1,1),1,inv_ovrlp_half(1,1),1)
        end if
        call timing(iproc,'lovrlp^-1/2com','OF')
        call timing(iproc,'lovrlp^-1/2old','ON')

        iall=-product(shape(inv_ovrlp_halfp))*kind(inv_ovrlp_halfp)
        deallocate(inv_ovrlp_halfp, stat=istat)
        call memocc(istat, iall, 'inv_ovrlp_halfp', subname)

     else ! no orbs present, use serial version
        ! No matrix compression available
        if (plusminus) then
           !$omp parallel do default(private) shared(inv_ovrlp_half,ovrlp,norb)
           do iorb=1,norb
              do jorb=iorb,norb
                 if(iorb==jorb) then
                    inv_ovrlp_half(jorb,iorb)=0.5d0+0.5d0*ovrlp(jorb,iorb)
                 else
                    inv_ovrlp_half(jorb,iorb)=0.5d0*ovrlp(jorb,iorb)
                    inv_ovrlp_half(iorb,jorb)=0.5d0*ovrlp(jorb,iorb)
                 end if
              end do
           end do
           !$omp end parallel do
        else
           !$omp parallel do default(private) shared(inv_ovrlp_half,ovrlp,norb)
           do iorb=1,norb
              do jorb=iorb,norb
                 if(iorb==jorb) then
                    inv_ovrlp_half(jorb,iorb)=1.5d0-.5d0*ovrlp(jorb,iorb)
                 else
                    inv_ovrlp_half(jorb,iorb)=-.5d0*ovrlp(jorb,iorb)
                    inv_ovrlp_half(iorb,jorb)=-.5d0*ovrlp(jorb,iorb)
                 end if
              end do
           end do
           !$omp end parallel do
        end if
     end if
  else

     ! Taylor expansion up to 4th order
      stop 'must fix this'

  end if

  call timing(iproc,'lovrlp^-1/2old','OF')

end subroutine overlapPowerPlusMinusOneHalf_old


!!subroutine first_order_taylor_sparse(power,ovrlp,inv_ovrlp)
!!  use module_base
!!  use module_types
!!  use sparsematrix_base, only: sparse_matrix
!!  use sparsematrix_init, only: matrixindex_in_compressed
!!  implicit none
!!  integer, intent(in) :: power
!!  type(sparse_matrix),intent(in) :: ovrlp
!!  type(sparse_matrix),intent(out) :: inv_ovrlp
!!
!!  integer :: ii,iii,iorb,jorb,ii_inv,ierr,iii_inv
!!
!!  if (power/=1 .and. power/=2 .and. power/=-2) stop 'Error in first_order_taylor_sparse'
!!
!!  if (inv_ovrlp%parallel_compression==0.or.bigdft_mpi%nproc==1) then
!!     ! inv_ovrlp can be less sparse than ovrlp, so pad with zeros first
!!     call to_zero(inv_ovrlp%nvctr,inv_ovrlp%matrix_compr(1))
!!     if (power==1) then
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii=1,ovrlp%nvctr
!!           iorb = ovrlp%orb_from_index(1,ii)
!!           jorb = ovrlp%orb_from_index(2,ii)
!!           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
!!           if(iorb==jorb) then
!!              inv_ovrlp%matrix_compr(ii_inv) = 2.0d0 - ovrlp%matrix_compr(ii)
!!           else
!!              inv_ovrlp%matrix_compr(ii_inv) = -ovrlp%matrix_compr(ii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     else if (power==2) then
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii=1,ovrlp%nvctr
!!           iorb = ovrlp%orb_from_index(1,ii)
!!           jorb = ovrlp%orb_from_index(2,ii)
!!           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
!!           if(iorb==jorb) then
!!              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0 + 0.5d0*ovrlp%matrix_compr(ii)
!!           else
!!              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0*ovrlp%matrix_compr(ii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     else
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii=1,ovrlp%nvctr
!!           iorb = ovrlp%orb_from_index(1,ii)
!!           jorb = ovrlp%orb_from_index(2,ii)
!!           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
!!           if(iorb==jorb) then
!!              inv_ovrlp%matrix_compr(ii_inv) = 1.5d0 - 0.5d0*ovrlp%matrix_compr(ii)
!!           else
!!              inv_ovrlp%matrix_compr(ii_inv) = -0.5d0*ovrlp%matrix_compr(ii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     end if
!!  else if (inv_ovrlp%parallel_compression==1) then
!!     call to_zero(inv_ovrlp%nvctr,inv_ovrlp%matrix_compr(1))
!!     if (power==1) then
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii=1,ovrlp%nvctrp
!!           iii=ii+ovrlp%isvctr
!!           iorb = ovrlp%orb_from_index(1,iii)
!!           jorb = ovrlp%orb_from_index(2,iii)
!!           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
!!           if(iorb==jorb) then
!!              inv_ovrlp%matrix_compr(ii_inv) = 2.0d0 - ovrlp%matrix_compr(iii)
!!           else
!!              inv_ovrlp%matrix_compr(ii_inv) = -ovrlp%matrix_compr(iii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     else if (power==2) then
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii=1,ovrlp%nvctrp
!!           iii=ii+ovrlp%isvctr
!!           iorb = ovrlp%orb_from_index(1,iii)
!!           jorb = ovrlp%orb_from_index(2,iii)
!!           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
!!           if(iorb==jorb) then
!!              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0 + 0.5d0*ovrlp%matrix_compr(iii)
!!           else
!!              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0*ovrlp%matrix_compr(iii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     else
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii=1,ovrlp%nvctrp
!!           iii=ii+ovrlp%isvctr
!!           iorb = ovrlp%orb_from_index(1,iii)
!!           jorb = ovrlp%orb_from_index(2,iii)
!!           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
!!           if(iorb==jorb) then
!!              inv_ovrlp%matrix_compr(ii_inv) = 1.5d0 - 0.5d0*ovrlp%matrix_compr(iii)
!!           else
!!              inv_ovrlp%matrix_compr(ii_inv) = -0.5d0*ovrlp%matrix_compr(iii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     end if
!!  else
!!     inv_ovrlp%matrix_comprp=f_malloc_ptr((inv_ovrlp%nvctrp),id='inv_ovrlp%matrix_comprp')
!!     if (power==1) then
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii_inv=1,inv_ovrlp%nvctrp
!!           iii_inv=ii_inv+inv_ovrlp%isvctr
!!           iorb = inv_ovrlp%orb_from_index(1,iii_inv)
!!           jorb = inv_ovrlp%orb_from_index(2,iii_inv)
!!           ii = matrixindex_in_compressed(ovrlp,iorb,jorb)
!!           if (ii==0) then
!!              inv_ovrlp%matrix_comprp(ii_inv) = 0.0d0
!!           else if(iorb==jorb) then
!!              inv_ovrlp%matrix_comprp(ii_inv) = 2.0d0 - ovrlp%matrix_compr(ii)
!!           else
!!              inv_ovrlp%matrix_comprp(ii_inv) = -ovrlp%matrix_compr(ii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     else if (power==2) then
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii_inv=1,inv_ovrlp%nvctrp
!!           iii_inv=ii_inv+inv_ovrlp%isvctr
!!           iorb = inv_ovrlp%orb_from_index(1,iii_inv)
!!           jorb = inv_ovrlp%orb_from_index(2,iii_inv)
!!           ii = matrixindex_in_compressed(ovrlp,iorb,jorb)
!!           if (ii==0) then
!!              inv_ovrlp%matrix_comprp(ii_inv) = 0.0d0
!!           else if(iorb==jorb) then
!!              inv_ovrlp%matrix_comprp(ii_inv) = 0.5d0 + 0.5d0*ovrlp%matrix_compr(ii)
!!           else
!!              inv_ovrlp%matrix_comprp(ii_inv) = 0.5d0*ovrlp%matrix_compr(ii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     else
!!        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
!!        do ii_inv=1,inv_ovrlp%nvctrp
!!           iii_inv=ii_inv+inv_ovrlp%isvctr
!!           iorb = inv_ovrlp%orb_from_index(1,iii_inv)
!!           jorb = inv_ovrlp%orb_from_index(2,iii_inv)
!!           ii = matrixindex_in_compressed(ovrlp,iorb,jorb)
!!           if (ii==0) then
!!              inv_ovrlp%matrix_comprp(ii_inv) = 0.0d0
!!           else if(iorb==jorb) then
!!              inv_ovrlp%matrix_comprp(ii_inv) = 1.5d0 - 0.5d0*ovrlp%matrix_compr(ii)
!!           else
!!              inv_ovrlp%matrix_comprp(ii_inv) = -0.5d0*ovrlp%matrix_compr(ii)
!!           end if
!!        end do
!!        !$omp end parallel do
!!     end if
!!  end if
!!
!!end subroutine first_order_taylor_sparse






!!subroutine check_accur_overlap_minus_one_sparse(iproc, nproc, smat, norb, norbp, isorb, nseq, nout, &
!!           ivectorindex, amat_seq, bmatp, power, &
!!           max_error, mean_error, dmat_seq, cmatp)
!!  use module_base
!!  use sparsematrix_base, only: sparse_matrix
!!  use sparsematrix, only: sparsemm
!!  implicit none
!!  integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseq, nout, power
!!  type(sparse_matrix) :: smat
!!  integer,dimension(nseq),intent(in) :: ivectorindex
!!  real(kind=8),dimension(nseq),intent(in) :: amat_seq
!!  real(kind=8),dimension(norb,norbp),intent(in) :: bmatp
!!  real(kind=8),intent(out) :: max_error, mean_error
!!  real(kind=8),dimension(nseq),intent(in),optional :: dmat_seq
!!  real(kind=8),dimension(norb,norbp),intent(in),optional :: cmatp
!!
!!  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
!!  real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
!!  integer :: ierr, i,j
!!
!!  call f_routine(id='check_accur_overlap_minus_one_sparse')
!!
!!  tmpp=f_malloc0((/norb,norbp/),id='tmpp')
!!  if (power==1) then
!!     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
!!     !!     norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
!!     call sparsemm(smat, amat_seq, bmatp, tmpp)
!!     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, smat, max_error, mean_error)
!!  else if (power==2) then
!!      if (.not.present(cmatp)) stop 'cmatp not present'
!!     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
!!     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
!!     !write(*,*) 'iproc, sum(amat_seq), sum(bmatp)', iproc, sum(amat_seq), sum(bmatp)
!!     call sparsemm(smat, amat_seq, bmatp, tmpp)
!!     call max_matrix_diff_parallel(iproc, norb, norbp, isorb, tmpp, cmatp, smat, max_error, mean_error)
!!     !max_error=0.5d0*max_error
!!     !mean_error=0.5d0*mean_error
!!  else if (power==-2) then
!!     if (.not.present(dmat_seq)) stop 'dmat_seq not present'
!!     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
!!     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
!!     call sparsemm(smat, amat_seq, bmatp, tmpp)
!!     tmp2p=f_malloc0((/norb,norbp/),id='tmp2p')
!!     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
!!     !!     norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
!!     call sparsemm(smat, dmat_seq, tmpp, tmp2p)
!!     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, smat, max_error, mean_error)
!!     !max_error=0.5d0*max_error
!!     !mean_error=0.5d0*mean_error
!!     call f_free(tmp2p)
!!  else
!!     stop 'Error in check_accur_overlap_minus_one_sparse'
!!  end if
!!  call f_free(tmpp)
!!
!!  call f_release_routine()
!!
!!end subroutine check_accur_overlap_minus_one_sparse









!!subroutine deviation_from_unity(iproc, norb, ovrlp, deviation)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, norb
!!  real(8),dimension(norb,norb),intent(in):: ovrlp
!!  real(8),intent(out):: deviation
!!
!!  ! Local variables
!!  integer:: iorb, jorb
!!  real(8):: error
!!
!!  call timing(iproc,'dev_from_unity','ON') 
!!  deviation=0.d0
!!  do iorb=1,norb
!!     do jorb=1,norb
!!        if(iorb==jorb) then
!!           error=abs(ovrlp(jorb,iorb)-1.d0)
!!        else
!!           error=abs(ovrlp(jorb,iorb))
!!        end if
!!        deviation=max(error,deviation)
!!     end do
!!  end do
!!  call timing(iproc,'dev_from_unity','OF') 
!!
!!end subroutine deviation_from_unity





!!subroutine deviation_from_symmetry(iproc, norb, ovrlp, deviation)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, norb
!!  real(8),dimension(norb,norb),intent(in):: ovrlp
!!  real(8),intent(out):: deviation
!!
!!  ! Local variables
!!  integer:: iorb, jorb
!!  real(8):: error
!!
!!  call timing(iproc,'dev_from_unity','ON') 
!!  deviation=0.d0
!!  do iorb=1,norb
!!     do jorb=iorb+1,norb
!!        error=abs(ovrlp(jorb,iorb)-ovrlp(iorb,jorb))
!!        deviation=max(error,deviation)
!!     end do
!!  end do
!!  call timing(iproc,'dev_from_unity','OF') 
!!
!!end subroutine deviation_from_symmetry



!!subroutine diagonalize_subset(iproc, nproc, orbs, ovrlp, ovrlp_mat, ham, ham_mat)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  use sparsematrix_base, only: sparse_matrix, matrices!, matrices_null, &
!!                               !allocate_matrices, deallocate_matrices
!!  use sparsematrix_init, only: matrixindex_in_compressed
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in) :: iproc, nproc
!!  type(orbitals_data),intent(inout) :: orbs
!!  type(sparse_matrix),intent(in) :: ovrlp
!!  type(matrices),intent(in) :: ovrlp_mat
!!  type(sparse_matrix),intent(in) :: ham
!!  type(matrices),intent(in) :: ham_mat
!!
!!  ! Local variables
!!  integer :: iend, i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
!!  integer :: iiorb, ierr, ii, iseg, ind, lwork, idiag, ind_ham, ind_ovrlp, info
!!  real(kind=8) :: error
!!  real(kind=8),dimension(:,:),pointer :: ovrlp_tmp, ham_tmp
!!  real(kind=8),dimension(:),allocatable :: eval, work
!!  logical,dimension(:),allocatable :: in_neighborhood
!!  !type(matrices) :: inv_ovrlp_half_
!!
!!  call f_routine('diagonalize_subset')
!!
!!  in_neighborhood = f_malloc(orbs%norb,id='in_neighborhood')
!!
!!  call f_zero(orbs%norb, orbs%eval(1))
!!
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb+iorb
!!     !ilr=orbs%inwhichlocreg(iiorb)
!!     ! We are at the start of a new atom
!!     ! Count all orbitals that are in the neighborhood
!!
!!     iseg=ham%istsegline(iiorb)
!!     iend=iiorb*orbs%norb
!!     n=0
!!     in_neighborhood(:)=.false.
!!     do 
!!        do i=ham%keyg(1,iseg),ham%keyg(2,iseg)
!!           ii=i-(iiorb-1)*orbs%norb
!!           in_neighborhood(ii)=.true.
!!           n=n+1
!!           if (ii==iiorb) then
!!               !this is the diagonal element
!!               idiag=n
!!           end if
!!        end do
!!        iseg=iseg+1
!!        if (iseg>ham%nseg) exit
!!        if (ham%keyg(1,iseg)>iend) exit
!!     end do
!!
!!     ham_tmp = f_malloc0_ptr((/n,n/),id='ovrlp_tmp')
!!     ovrlp_tmp = f_malloc0_ptr((/n,n/),id='ovrlp_tmp')
!!
!!     jjorb=0
!!     do jorb=1,orbs%norb
!!        if (.not.in_neighborhood(jorb)) cycle
!!        jjorb=jjorb+1
!!        kkorb=0
!!        do korb=1,orbs%norb
!!           if (.not.in_neighborhood(korb)) cycle
!!           kkorb=kkorb+1
!!           ind_ham = matrixindex_in_compressed(ham,korb,jorb)
!!           ind_ovrlp = matrixindex_in_compressed(ovrlp,korb,jorb)
!!           if (ind_ham>0) then
!!               ham_tmp(kkorb,jjorb)=ham_mat%matrix_compr(ind_ham)
!!           else
!!               ham_tmp(kkorb,jjorb)=0.d0
!!           end if
!!           if (ind_ovrlp>0) then
!!              ovrlp_tmp(kkorb,jjorb)=ovrlp_mat%matrix_compr(ind_ovrlp)
!!           else
!!              ovrlp_tmp(kkorb,jjorb)=0.d0
!!           end if
!!           !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
!!        end do
!!     end do
!!
!!     lwork=100*n
!!     work = f_malloc(lwork,id='work')
!!     eval = f_malloc(n,id='eval')
!!     call dsygv(1, 'n', 'l', n, ham_tmp, n, ovrlp_tmp, n, eval, work, lwork, info)
!!     orbs%eval(iiorb)=eval(idiag)
!!     call f_free(work)
!!     call f_free(eval)
!!     call f_free_ptr(ham_tmp)
!!     call f_free_ptr(ovrlp_tmp)
!!
!!
!!
!!
!! end do
!!
!! call f_free(in_neighborhood)
!!
!! if (nproc>1) then
!!     call mpiallred(orbs%eval(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm)
!! end if
!!
!!  call f_release_routine()
!!
!!end subroutine diagonalize_subset
