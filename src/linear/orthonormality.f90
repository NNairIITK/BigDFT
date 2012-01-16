subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, &
     blocksize_pdgemm, orbs, op, comon, lzd, onWhichAtomAll, convCritOrtho, input, mad, lphi, ovrlp, method)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc,nproc,methTransformOverlap,nItOrtho,blocksize_dsyev,blocksize_pdgemm
  !type(linearParameters),intent(inout):: lin
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(inout):: op
  type(p2pCommsOrthonormality),intent(inout):: comon
  type(local_zone_descriptors),intent(in):: lzd
  integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
  real(8),intent(in):: convCritOrtho
  type(input_variables),intent(in):: input
  !real(8),dimension(lin%lorbs%npsidim),intent(inout):: lphi
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(inout) :: lphi
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
character(len=3),intent(in):: method

  ! Local variables
  integer:: it, istat, iall, iorb, jorb, ierr, ind1, ind2, maxvaloverlap
  real(8),dimension(:),allocatable:: lphiovrlp, sendbufcopy
  character(len=*),parameter:: subname='orthonormalize'
  logical:: converged
  real(8):: maxError, t1, t2, timeCommun, timeComput, timeCalcOvrlp, t3, t4, timeExpand, timeLoewdin, timeTransform, timeExtract
  real(8):: timecommunp2p, timecommuncoll, timeoverlap, timecompress



  timeComput=0.d0
  timeCommun=0.d0
  timeLoewdin=0.d0
  timeTransform=0.d0
  timeextract=0.d0
  timecommunp2p=0.d0
  timecommuncoll=0.d0
  timeoverlap=0.d0
  timeexpand=0.d0
  timecompress=0.d0
  converged=.false.
  !write(*,*) 'iproc, comon%nsend', iproc, comon%nsend
  !do it=1,nItOrtho+orbs%norb
  do it=1,nItOrtho
     t1=mpi_wtime()

     ! Allocate the send and receive buffers for the communication.
     call allocateSendBufferOrtho(comon, subname)
     call allocateRecvBufferOrtho(comon, subname)

     ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
     call extractOrbital3(iproc, nproc, orbs, orbs%npsidim_orbs, onWhichAtomAll, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)

     t2=mpi_wtime()
     timeExtract=timeExtract+t2-t1
     timeComput=timeComput+t2-t1

     t1=mpi_wtime()
     ! Post the send messages.
     call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, timecommunp2p, timeextract)
     allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
     call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
     !call collectAndCalculateOverlap(iproc, nproc, comon, mad, op, orbs, input, lzd, comon%nsendbuf, &
     !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ovrlp, lphiovrlp, timecommunp2p, timecommuncoll, timeoverlap, timeexpand, timecompress)
     call collectnew(iproc, nproc, comon, mad, op, orbs, input, lzd, comon%nsendbuf, &
          comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ovrlp, timecommunp2p, timecommuncoll, timecompress)

!!! THIS IS NEW #####################
     !!maxvaloverlap=maxval(comon%noverlaps)
     !!comon%nstepoverlap=50000
     !!comon%isoverlap=1
     !!do 
     !!    call postCommsOverlapNew2(iproc, nproc, orbs, op, lzd, lphi, comon, timecommunp2p, timeextract)
     !!    !write(*,'(a,3i9)') 'iproc, comon%nrecv, comon%nsend', iproc, comon%nrecv, comon%nsend
     !!    call collectnew2(iproc, nproc, comon, mad, op, orbs, input, lzd, comon%nsendbuf, &
     !!         comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ovrlp, timecommunp2p, timecommuncoll, timecompress)
     !!    comon%isoverlap=comon%isoverlap+comon%nstepoverlap
     !!    if(comon%isoverlap>=maxvaloverlap) exit
     !!end do
!!! #################################

     !!do ind1=1,orbs%norb
     !!    do ind2=1,orbs%norb
     !!        write(500+iproc,*) ind1, ind2, ovrlp(ind2,ind1)
     !!    end do
     !!end do
     t2=mpi_wtime()
     timeCommun=timeCommun+t2-t1
     t1=mpi_wtime()
     !call deallocateSendBufferOrtho(comon, subname)
      !call expandOrbital2(iproc, nproc, orbs, input, onWhichAtomAll, lzd, op, comon, lphiovrlp)
     t2=mpi_wtime()
     timeexpand=timeexpand+t2-t1
     !call checkUnity(iproc, orbs%norb, ovrlp, maxError)
     !if(iproc==0) write(*,'(3x,a,es12.4)') 'maximal deviation from unity:', maxError
     t3=mpi_wtime()

     !call overlapMatrixCubic(iproc, nproc, gorbs, orbs, comms, lzd, input, lphi, ovrlp)
     call calculateOverlapMatrix3(iproc, nproc, orbs, op, orbs%inWhichLocreg, comon%nsendBuf, &
          comon%sendBuf, comon%nrecvBuf, comon%recvBuf, mad, ovrlp)
     t4=mpi_wtime()
     call checkUnity(iproc, orbs%norb, ovrlp, maxError)
      if(iproc==0) write(*,*) 'deviation from unity:', maxError
      !!if(maxError<1.d-2) then
      !!    if(iproc==0) write(*,*) 'error is small enough, no need for orthogonalization, only normalization...'
      !!    ist=1
      !!    do iorb=1,orbs%norbp
      !!        iiorb=orbs%isorb+iorb
      !!        ilr=orbs%inwhichlocreg(iiorb)
      !!        ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      !!        ! Normalize
      !!        tt=dnrm2(ncount, lphi(ist), 1)
      !!        call dscal(ncount, 1/tt, lphi(ist), 1)
      !!        ist=ist+ncount
      !!    end do
      !!    call deallocateRecvBufferOrtho(comon, subname)
      !!    call deallocateSendBufferOrtho(comon, subname)
      !!    iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
      !!    deallocate(lphiovrlp, stat=istat)
      !!    call memocc(istat, iall, 'lphiovrlp', subname)
      !!    exit
      !!end if
     timeoverlap=timeoverlap+t4-t3

     t3=mpi_wtime()
      !!if(method=='old' .or. maxError>2.d-20) then
          call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, methTransformOverlap, blocksize_dsyev, &
                blocksize_pdgemm, orbs%norb, mad, ovrlp)
      !!end if
     t4=mpi_wtime()
     timeTransform=timeTransform+t4-t3
     t3=mpi_wtime()


     t4=mpi_wtime()
     timeExpand=timeExpand+t4-t3
     t3=mpi_wtime()
      !!if(method=='old' .or. maxError>2.d-20) then
          call globalLoewdin(iproc, nproc, orbs, orbs, onWhichAtomAll, lzd, op, comon, ovrlp, lphiovrlp, lphi)
      !!end if
      call deallocateRecvBufferOrtho(comon, subname)
      call deallocateSendBufferOrtho(comon, subname)

      !!if(method=='new' .and. maxError<=2.d-20) then
      !!    call localloewdin(iproc, nproc, lzd, orbs, op, ovrlp, lphiovrlp, lphi)
      !!end if

      !!call mpi_barrier(mpi_comm_world, ierr)
      !!call gramschmidt_forlinear(iproc, nproc, it, orbs, orbs, onWhichAtomAll, lzd, op, ovrlp, lphiovrlp, lphi)
      !!call mpi_barrier(mpi_comm_world, ierr)
     iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
     deallocate(lphiovrlp, stat=istat)
     call memocc(istat, iall, 'lphiovrlp', subname)
     t4=mpi_wtime()
     timeLoewdin=timeLoewdin+t4-t3
     t2=mpi_wtime()
     timeComput=timeComput+t2-t1

      !if(it==nItOrtho+orbs%norb) then
     if(it==nItOrtho) then
        converged=.false.
        exit
     end if
     !if(maxError<convCritOrtho) then
     !    converged=.true.
     !    exit
     !else if(it==nItOrtho) then
     !    exit
     !end if
  end do

  if(converged) then
     if(iproc==0) write(*,'(3x,a,i0,a)') 'done in ', it, ' iterations.'
  else 
     if(iproc==0) write(*,'(3x,a,i0,a)') 'WARNING: orthonormalization not converged within ', nItOrtho, ' iterations.'
  end if



  !write(*,*)'iproc, timecommunp2p', iproc, timecommunp2p
  !write(*,*)'iproc, timecommuncoll', iproc, timecommuncoll

  if (verbose > 2) then
     timeComput=timeLoewdin+timeTransform+timeextract+timeoverlap+timeexpand+timecompress
     timeCommun=timecommunp2p+timecommuncoll
     call mpiallred(timeComput, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeCommun, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeLoewdin, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeTransform, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeextract, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timecommunp2p, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timecommuncoll, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeoverlap, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeexpand, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timecompress, 1, mpi_sum, mpi_comm_world, ierr)
     timeComput=timeComput/dble(nproc)
     timeCommun=timeCommun/dble(nproc)
     timeLoewdin=timeLoewdin/dble(nproc)
     timeTransform=timeTransform/dble(nproc)
     timeExtract=timeExtract/dble(nproc)
     timecommunp2p=timecommunp2p/dble(nproc)
     timecommuncoll=timecommuncoll/dble(nproc)
     timeoverlap=timeoverlap/dble(nproc)
     timeexpand=timeexpand/dble(nproc)
     timecompress=timecompress/dble(nproc)
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for computation:', timeComput, '=',&
          100.d0*timeComput/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for communication:', timeCommun, '=',&
          100.d0*timeCommun/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for calculating overlap:', timeoverlap, '=',&
          100.d0*timeoverlap/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for expansion:', timeexpand, '=',&
          100.d0*timeexpand/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for Loewdin:', timeLoewdin, '=',&
          100.d0*timeLoewdin/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for transform:', timeTransform, '=',&
          100.d0*timeTransform/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for extract:', timeExtract, '=',&
          100.d0*timeExtract/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for p2p communication:', timecommunp2p, '=',&
          100.d0*timecommunp2p/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for collective communication:', timecommuncoll, '=',&
          100.d0*timecommuncoll/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for expand:', timeExpand, '=',&
          100.d0*timeExpand/(timeComput+timeCommun), '%'
  end if


end subroutine orthonormalizeLocalized




!subroutine orthonormalizeLocalized2(iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, &
!     blocksize_pdgemm, orbs, op, comon, lzd, gorbs, comms, onWhichAtomAll, convCritOrtho, input, mad, lphi, ovrlp)
!  use module_base
!  use module_types
!  use module_interfaces, exceptThisOne => orthonormalizeLocalized2
!  implicit none
!
!  ! Calling arguments
!  integer,intent(in):: iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, blocksize_pdgemm
!  !type(linearParameters),intent(inout):: lin
!  type(orbitals_data),intent(in):: orbs, gorbs
!  type(overlapParameters),intent(inout):: op
!  type(p2pCommsOrthonormality),intent(inout):: comon
!  type(local_zone_descriptors),intent(in):: lzd
!  type(communications_arrays),intent(in):: comms
!  integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
!  real(8),intent(in):: convCritOrtho
!  type(input_variables),intent(in):: input
!  !real(8),dimension(lin%lorbs%npsidim),intent(inout):: lphi
!  type(matrixDescriptors),intent(in):: mad
!  real(8),dimension(orbs%npsidim),intent(inout):: lphi
!  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
!
!  ! Local variables
!  integer:: it, istat, iall, iorb, jorb, ierr
!  real(8),dimension(:),allocatable:: lphiovrlp, sendbufcopy
!  character(len=*),parameter:: subname='orthonormalize'
!  logical:: converged
!  real(8):: maxError, t1, t2, timeCommun, timeComput, timeCalcOvrlp, t3, t4, timeExpand, timeLoewdin, timeTransform, timeExtract
!  real(8):: timecommunp2p, timecommuncoll, timeoverlap, timecompress
!
!
!
!  timeComput=0.d0
!  timeCommun=0.d0
!  timeLoewdin=0.d0
!  timeTransform=0.d0
!  timeextract=0.d0
!  timecommunp2p=0.d0
!  timecommuncoll=0.d0
!  timeoverlap=0.d0
!  timeexpand=0.d0
!  timecompress=0.d0
!  converged=.false.
!  call overlapMatrixCubic(iproc, nproc, gorbs, orbs, comms, lzd, input, methTransformOverlap, blocksize_dsyev, &
!       blocksize_pdgemm, mad, lphi)
!  !do it=1,nItOrtho
!  !    t1=mpi_wtime()
!
!  !    ! Allocate the send and receive buffers for the communication.
!  !    call allocateSendBufferOrtho(comon, subname)
!  !    call allocateRecvBufferOrtho(comon, subname)
!  !    
!  !    ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
!  !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, onWhichAtomAll, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
!
!  !    t2=mpi_wtime()
!  !    timeExtract=timeExtract+t2-t1
!  !    timeComput=timeComput+t2-t1
!
!  !    t1=mpi_wtime()
!  !    ! Post the send messages.
!  !    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, timecommunp2p, timeextract)
!  !    allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
!  !    call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
!  !    call expandOrbital2(iproc, nproc, orbs, input, onWhichAtomAll, lzd, op, comon, lphiovrlp)
!  !    !call collectAndCalculateOverlap(iproc, nproc, comon, mad, op, orbs, input, lzd, comon%nsendbuf, &
!  !    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ovrlp, lphiovrlp, timecommunp2p, timecommuncoll, timeoverlap, timeexpand, timecompress)
!  !    t2=mpi_wtime()
!  !    timeCommun=timeCommun+t2-t1
!  !    t1=mpi_wtime()
!  !    call deallocateSendBufferOrtho(comon, subname)
!  !    t2=mpi_wtime()
!  !    !call checkUnity(iproc, orbs%norb, ovrlp, maxError)
!  !    !if(iproc==0) write(*,'(3x,a,es12.4)') 'maximal deviation from unity:', maxError
!  !    t3=mpi_wtime()
!
!  !    !!call overlapMatrixCubic(iproc, nproc, gorbs, orbs, comms, lzd, input, lphi, ovrlp)
!
!  !    call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, methTransformOverlap, blocksize_dsyev, &
!  !          blocksize_pdgemm, orbs%norb, mad, ovrlp)
!  !    t4=mpi_wtime()
!  !    timeTransform=timeTransform+t4-t3
!  !    t3=mpi_wtime()
!  !    call deallocateRecvBufferOrtho(comon, subname)
!  !    t4=mpi_wtime()
!  !    timeExpand=timeExpand+t4-t3
!  !    t3=mpi_wtime()
!  !    call globalLoewdin(iproc, nproc, orbs, orbs, onWhichAtomAll, lzd, op, ovrlp, lphiovrlp, lphi)
!  !    iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
!  !    deallocate(lphiovrlp, stat=istat)
!  !    call memocc(istat, iall, 'lphiovrlp', subname)
!  !    t4=mpi_wtime()
!  !    timeLoewdin=timeLoewdin+t4-t3
!  !    t2=mpi_wtime()
!  !    timeComput=timeComput+t2-t1
!
!  !    if(it==nItOrtho) then
!  !        converged=.false.
!  !        exit
!  !    end if
!  !    !if(maxError<convCritOrtho) then
!  !    !    converged=.true.
!  !    !    exit
!  !    !else if(it==nItOrtho) then
!  !    !    exit
!  !    !end if
!  !end do
!
!  !if(converged) then
!  !    if(iproc==0) write(*,'(3x,a,i0,a)') 'done in ', it, ' iterations.'
!  !else 
!  !    if(iproc==0) write(*,'(3x,a,i0,a)') 'WARNING: orthonormalization not converged within ', nItOrtho, ' iterations.'
!  !end if
!
!
!
!  !write(*,*)'iproc, timecommunp2p', iproc, timecommunp2p
!  !write(*,*)'iproc, timecommuncoll', iproc, timecommuncoll
!
!  if (verbose > 2) then
!     timeComput=timeLoewdin+timeTransform+timeextract+timeoverlap+timeexpand+timecompress
!     timeCommun=timecommunp2p+timecommuncoll
!     call mpiallred(timeComput, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timeCommun, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timeLoewdin, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timeTransform, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timeextract, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timecommunp2p, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timecommuncoll, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timeoverlap, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timeexpand, 1, mpi_sum, mpi_comm_world, ierr)
!     call mpiallred(timecompress, 1, mpi_sum, mpi_comm_world, ierr)
!     timeComput=timeComput/dble(nproc)
!     timeCommun=timeCommun/dble(nproc)
!     timeLoewdin=timeLoewdin/dble(nproc)
!     timeTransform=timeTransform/dble(nproc)
!     timeExtract=timeExtract/dble(nproc)
!     timecommunp2p=timecommunp2p/dble(nproc)
!     timecommuncoll=timecommuncoll/dble(nproc)
!     timeoverlap=timeoverlap/dble(nproc)
!     timeexpand=timeexpand/dble(nproc)
!     timecompress=timecompress/dble(nproc)
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for computation:', timeComput, '=',&
!          100.d0*timeComput/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for communication:', timeCommun, '=',&
!          100.d0*timeCommun/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for calculating overlap:', timeoverlap, '=',&
!          100.d0*timeoverlap/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for expansion:', timeexpand, '=',&
!          100.d0*timeexpand/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for Loewdin:', timeLoewdin, '=',&
!          100.d0*timeLoewdin/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for transform:', timeTransform, '=',&
!          100.d0*timeTransform/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for extract:', timeExtract, '=',&
!          100.d0*timeExtract/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for p2p communication:', timecommunp2p, '=',&
!          100.d0*timecommunp2p/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for collective communication:', timecommuncoll, '=',&
!          100.d0*timecommuncoll/(timeComput+timeCommun), '%'
!     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for expand:', timeExpand, '=',&
!          100.d0*timeExpand/(timeComput+timeCommun), '%'
!
!  end if
!
!end subroutine orthonormalizeLocalized2







subroutine orthoconstraintNonorthogonal(iproc, nproc, lin, input, ovrlp, lphi, lhphi, mad, trH, W, eval)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(linearParameters),intent(inout):: lin
  type(input_variables),intent(in):: input
  real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: ovrlp
  !real(8),dimension(lin%lorbs%npsidim),intent(in):: lphi
  real(8),dimension(max(lin%orbs%npsidim_comp,lin%orbs%npsidim_orbs)),intent(in):: lphi
  !real(8),dimension(lin%lorbs%npsidim),intent(inout):: lhphi
  real(8),dimension(max(lin%orbs%npsidim_comp,lin%orbs%npsidim_orbs)),intent(inout):: lhphi
  type(matrixDescriptors),intent(in):: mad
  real(8),intent(out):: trH
  real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(out),optional:: W
  real(8),dimension(lin%orbs%norb),intent(out),optional:: eval

  ! Local variables
  integer:: it, istat, iall, iorb, ierr, i, maxvaloverlap, lwork, info, jorb, k
  real(8),dimension(:),allocatable:: lphiovrlp, sendBuf, ovrlpCompressed, ovrlpCompressed2, work
  real(8),dimension(:,:),allocatable:: lagmat
  character(len=*),parameter:: subname='orthoconstraintLocalized'
  real(8):: t1, t2, timeExtract, timeExpand, timeApply, timeCalcMatrix, timeCommun, timeComput, tt
  real(8):: timecommunp2p, timecommuncoll, timecompress
  logical,dimension(:,:),allocatable:: expanded
  integer,dimension(:),allocatable:: sendcounts, displs
  logical:: present_W, present_eval



  allocate(lagmat(lin%orbs%norb,lin%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat',subname)
  allocate(lphiovrlp(lin%op%ndim_lphiovrlp), stat=istat)
  call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
  allocate(expanded(lin%orbs%norb,lin%orbs%norbp), stat=istat)
  call memocc(istat, expanded, 'expanded',subname)
  allocate(sendBuf(lin%comon%nsendBuf), stat=istat)
  call memocc(istat, sendBuf, 'sendBuf',subname)
  lphiovrlp=0.d0


  timeComput=0.d0
  timeCommun=0.d0
  timeCalcMatrix=0.d0
  timeExpand=0.d0
  timeApply=0.d0
  timeExtract=0.d0

  !!allocate(lhphiovrlp(lin%op%ndim_lphiovrlp), stat=istat)
  !!call memocc(istat, lhphiovrlp, 'lhphiovrlp',subname)

  call allocateCommuncationBuffersOrtho(lin%comon, subname)

  ! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
  !call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
  t1=mpi_wtime()
  !call extractOrbital2(iproc, nproc, lin%orbs, lin%orbs%npsidim, lin%orbs%inWhichLocreg, lin%lzd, lin%op, lphi, lin%comon)
  call extractOrbital3(iproc, nproc, lin%orbs, lin%orbs%npsidim_orbs, lin%orbs%inWhichLocreg, lin%lzd, lin%op, &
       lphi, lin%comon%nsendBuf, lin%comon%sendBuf)
  t2=mpi_wtime()
  timeExtract=t2-t1
  !call postCommsOverlap(iproc, nproc, lin%comon)
  call postCommsOverlapNew(iproc, nproc, lin%orbs, lin%op, lin%lzd, lphi, lin%comon, timecommunp2p, timeextract)
  !call extractOrbital3(iproc, nproc, lin%orbs, lin%orbs%npsidim, lin%orbs%inWhichLocreg, lin%lzd, lin%op, &
  !                     lhphi, lin%comon%nsendBuf, sendBuf)
  !call gatherOrbitals2(iproc, nproc, lin%comon)
  !call gatherOrbitalsOverlapWithComput(iproc, nproc, lin%orbs, input, lin%lzd, lin%op, lin%comon, lphiovrlp, expanded)
  !!call gatherOrbitalsOverlapWithComput2(iproc, nproc, lin%orbs, input, lin%lzd, lin%op, lin%comon, lin%comon%nsendbuf, sendbuf, lin%comon%nrecvbuf, lin%comon%recvbuf, lphiovrlp, expanded, lagmat)
  call collectnew(iproc, nproc, lin%comon, lin%mad, lin%op, lin%orbs, input, lin%lzd, lin%comon%nsendbuf, &
       lin%comon%sendbuf, lin%comon%nrecvbuf, lin%comon%recvbuf, lagmat, timecommunp2p, timecommuncoll, timecompress)
!!!! THIS IS NEW #####################
  !!maxvaloverlap=maxval(lin%comon%noverlaps)
  !!lin%comon%nstepoverlap=50000
  !!lin%comon%isoverlap=1
  !!do 
  !!    call postCommsOverlapNew2(iproc, nproc, lin%orbs, lin%op, lin%lzd, lphi, lin%comon, timecommunp2p, timeextract)
  !!    !write(*,'(a,3i9)') 'iproc, lin%comon%nrecv, lin%comon%nsend', iproc, lin%comon%nrecv, lin%comon%nsend
  !!    call collectnew2(iproc, nproc, lin%comon, lin%mad, lin%op, lin%orbs, input, lin%lzd, lin%comon%nsendbuf, &
  !!         lin%comon%sendbuf, lin%comon%nrecvbuf, lin%comon%recvbuf, lagmat, timecommunp2p, timecommuncoll, timecompress)
  !!    lin%comon%isoverlap=lin%comon%isoverlap+lin%comon%nstepoverlap
  !!    if(lin%comon%isoverlap>=maxvaloverlap) exit
  !!end do
  !! #################################
  ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
  !call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lhphi, lin%comon)
  !!do iall=1,lin%comon%nrecvBuf
  !!    !write(3000+iproc,*) iall, lin%comon%sendBuf(iall)
  !!    !write(3000+iproc,*) iall, sendBuf(iall)
  !!    write(4000+iproc,*) iall, lin%comon%recvBuf(iall)
  !!end do
  !call calculateOverlapMatrix2(iproc, nproc, lin%orbs, lin%op, lin%comon, lin%orbs%inWhichLocreg, lagmat)
  t1=mpi_wtime()
  !call extractOrbital3(iproc, nproc, lin%orbs, lin%orbs%npsidim, lin%orbs%inWhichLocreg, lin%lzd, lin%op, &
  !                     lhphi, lin%comon%nsendBuf, sendBuf)
  call extractOrbital3(iproc, nproc, lin%orbs, lin%orbs%npsidim_orbs, lin%orbs%inWhichLocreg, lin%lzd, lin%op, &
       lhphi, lin%comon%nsendBuf, lin%comon%sendBuf)
  t2=mpi_wtime()
  timeextract=timeextract+t2-t1
  t1=mpi_wtime()
  call calculateOverlapMatrix3(iproc, nproc, lin%orbs, lin%op, lin%orbs%inWhichLocreg, lin%comon%nsendBuf, &
       lin%comon%sendBuf, lin%comon%nrecvBuf, lin%comon%recvBuf, mad, lagmat)
  !!do iorb=1,lin%orbs%norb
  !!    do jorb=1,lin%orbs%norb
  !!        if(iproc==0) write(332,'(2i9,es20.10)') iorb,jorb,lagmat(jorb,iorb)
  !!    end do 
  !!end do 
  t2=mpi_wtime()
  timecalcmatrix=timecalcmatrix+t2-t1

  present_W=present(W)
  present_eval=present(eval)
  if(present_W .and. .not.present_eval .or. present_eval .and. .not.present_W) then
      write(*,*) 'ERROR W and eval must be present at the same time'
      stop
  end if
  if(present_W .and. present_eval) then
      !! Just copy lagmat... used for debugging
      if(iproc==0) write(*,*) 'copy lagmat to W'
      W=lagmat
  end if


  !call cpu_time(t2)
  !timeCalcMatrix=t2-t1
  trH=0.d0
  do iorb=1,lin%orbs%norb
     trH=trH+lagmat(iorb,iorb)
  end do
  ! Expand the receive buffer, i.e. lphi
  t1=mpi_wtime()
  !call expandOrbital2(iproc, nproc, lin%orbs, input, lin%orbs%inWhichLocreg, lin%lzd, lin%op, lin%comon, lphiovrlp)
  !call expandRemainingOrbitals(iproc, nproc, lin%orbs, input, lin%orbs%inWhichLocreg, lin%lzd, lin%op, &
  !     lin%comon, expanded, lphiovrlp)
  t2=mpi_wtime()
  timeExpand=t2-t1
  !!do iall=1,lin%op%ndim_lphiovrlp
  !!    write(2000+iproc,*) iall, lphiovrlp(iall)
  !!end do

  !! I think this is not needed??
!!!! Now we also have to send lhphi
!!!!call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lhphi, lin%comon)
!!!call extractOrbital2(iproc, nproc, lin%orbs, lin%orbs%npsidim, lin%orbs%inWhichLocreg, lin%lzd, lin%op, lhphi, lin%comon)
!!!call postCommsOverlap(iproc, nproc, lin%comon)
!!!call gatherOrbitals2(iproc, nproc, lin%comon)
!!!! Expand the receive buffer, i.e. lhphi
!!!call expandOrbital2(iproc, nproc, lin%orbs, input, lin%orbs%inWhichLocreg, lin%lzd, lin%op, lin%comon, lhphiovrlp)
  
  !!!!!! TEST #####################################
  !!!!if(iproc==0) write(*,*) 'WARNING: testing: deleting all off-diagonal elements of lagmat!!!'
  !!!!do iorb=1,lin%orbs%norb
  !!!!    do jorb=1,lin%orbs%norb
  !!!!        if(iorb/=jorb) lagmat(jorb,iorb)=0.d0
  !!!!    end do
  !!!!end do
  !!!!!! END TEST ##############################

  t1=mpi_wtime()

  if(iproc==0) write(*,*) 'old version of orthoconstraint'
  call applyOrthoconstraintNonorthogonal2(iproc, nproc, lin%methTransformOverlap, lin%blocksize_pdgemm, lin%orbs, lin%orbs, &
       lin%orbs%inWhichLocreg, lin%lzd, lin%op, lin%comon, lagmat, ovrlp, lphiovrlp, mad, lhphi)
  !!if(iproc==0) write(*,*) 'new version of orthoconstraint'
  !!call applyOrthoconstraintlocal(iproc, nproc, lin%lzd, lin%orbs, lin%op, lagmat, lphiovrlp, lhphi)

  t2=mpi_wtime()
  timeApply=t2-t1

  call deallocateCommuncationBuffersOrtho(lin%comon, subname)

  present_W=present(W)
  present_eval=present(eval)
  !present_eval=.false.
  if(present_W .and. .not.present_eval .or. present_eval .and. .not.present_W) then
     write(*,*) 'ERROR W and eval must be present at the same time'
     stop
  end if
  if(present_W .and. present_eval) then
      !! Just copy lagmat... used for debugging
      !!W=lagmat
      !!! Diagonalize lagmat -- EXPERIMENTAL

      !!! Copy lagmat to W and symmetrize it.
      !!do iorb=1,lin%orbs%norb
      !!    do jorb=1,lin%orbs%norb
      !!        W(jorb,iorb)=.5d0*lagmat(jorb,iorb)+.5d0*lagmat(iorb,jorb)
      !!    end do
      !!end do
      !!allocate(work(1), stat=istat)
      !!call memocc(istat, work, 'work', subname)
      !!call dsyev('v', 'l', lin%orbs%norb, W(1,1), lin%orbs%norb, eval, work, -1, info)
      !!if(info/=0) stop 'ERROR in dsyev'
      !!lwork=work(1)
      !!iall=-product(shape(work))*kind(work)
      !!deallocate(work, stat=istat)
      !!call memocc(istat, iall, 'work', subname)
      !!allocate(work(lwork), stat=istat)
      !!call memocc(istat, work, 'work', subname)
      !!call dsyev('v', 'l', lin%orbs%norb, W(1,1), lin%orbs%norb, eval, work, lwork, info)
      !!if(info/=0) stop 'ERROR in dsyev'
      !!iall=-product(shape(work))*kind(work)
      !!deallocate(work, stat=istat)
      !!call memocc(istat, iall, 'work', subname)

      !!! Transpose W
      !!do iorb=1,lin%orbs%norb
      !!    do jorb=iorb,lin%orbs%norb
      !!        tt=W(jorb,iorb)
      !!        W(jorb,iorb)=W(iorb,jorb)
      !!        W(iorb,jorb)=tt
      !!    end do
      !!end do

      !!! Check
      !!!do iorb=1,lin%orbs%norb
      !!!    do jorb=1,lin%orbs%norb
      !!!        tt=0.d0
      !!!        do k=1,lin%orbs%norb
      !!!            tt=tt+eval(k)*W(k,iorb)*W(k,jorb)
      !!!        end do
      !!!        if(iproc==0) then
      !!!            write(*,'(a,2i7,2es16.8)') 'iproc, jorb, tt, .5d0*lagmat(jorb,iorb)+.5d0*lagmat(iorb,jorb)', iproc, jorb, tt, .5d0*lagmat(jorb,iorb)+.5d0*lagmat(iorb,jorb)
      !!!        end if
      !!!    end do
      !!!end do
  end if

  !!

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)
  iall=-product(shape(expanded))*kind(expanded)
  deallocate(expanded, stat=istat)
  call memocc(istat, iall, 'expanded', subname)
  iall=-product(shape(sendBuf))*kind(sendBuf)
  deallocate(sendBuf, stat=istat)
  call memocc(istat, iall, 'sendBuf', subname)
  !!iall=-product(shape(lhphiovrlp))*kind(lhphiovrlp)
  !!deallocate(lhphiovrlp, stat=istat)
  !!call memocc(istat, iall, 'lhphiovrlp', subname)

  if (verbose > 2) then
     timeComput=timeExtract+timeExpand+timeApply+timecalcmatrix
     timeCommun=timecommunp2p+timecommuncoll
     call mpiallred(timeComput, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeCommun, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeCalcMatrix, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeExpand, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeApply, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(timeExtract, 1, mpi_sum, mpi_comm_world, ierr)
     timeComput=timeComput/dble(nproc)
     timeCommun=timeCommun/dble(nproc)
     timeCalcMatrix=timeCalcMatrix/dble(nproc)
     timeExpand=timeExpand/dble(nproc)
     timeApply=timeApply/dble(nproc)
     timeExtract=timeExtract/dble(nproc)
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for computation:', timeComput, '=', &
          100.d0*timeComput/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for communication:', timeCommun, '=', &
          100.d0*timeCommun/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for calculating overlap:', timeCalcMatrix, &
          '=', 100.d0*timeCalcMatrix/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for expansion:', timeExpand, '=', &
          100.d0*timeExpand/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for applying orthoconstraint:', timeApply, &
          '=', 100.d0*timeApply/(timeComput+timeCommun), '%'
     if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for extract:', timeExtract, '=', &
          100.d0*timeExtract/(timeComput+timeCommun), '%'
  end if

end subroutine orthoconstraintNonorthogonal





!!subroutine getOverlapMatrix(iproc, nproc, lin, input, lphi, mad, ovrlp)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => getOverlapMatrix
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(linearParameters),intent(inout):: lin
!!type(input_variables),intent(in):: input
!!!real(8),dimension(lin%lorbs%npsidim),intent(inout):: lphi
!!real(8),dimension(lin%orbs%npsidim),intent(inout):: lphi
!!type(matrixDescriptors),intent(in):: mad
!!real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(out):: ovrlp
!!
!!! Local variables
!!integer:: it, istat, iall, iorb, jorb, ierr
!!real(8),dimension(:),allocatable:: lphiovrlp
!!character(len=*),parameter:: subname='orthonormalize'
!!logical:: converged
!!
!!  call allocateCommuncationBuffersOrtho(lin%comon, subname)
!!
!!  !call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
!!  call extractOrbital2(iproc, nproc, lin%orbs, lin%orbs%npsidim, lin%orbs%inWhichLocreg, lin%lzd, lin%op, lphi, lin%comon)
!!  call postCommsOverlap(iproc, nproc, lin%comon)
!!  call gatherOrbitals2(iproc, nproc, lin%comon)
!!  !!call calculateOverlapMatrix2(iproc, nproc, lin%orbs, lin%op, lin%comon, lin%orbs%inWhichLocreg, mad, ovrlp)
!!  call calculateOverlapMatrix3(iproc, nproc, lin%orbs, lin%op, lin%orbs%inWhichLocreg, lin%comon%nsendBuf, &
!!                               lin%comon%sendBuf, lin%comon%nrecvBuf, lin%comon%recvBuf, mad, ovrlp)
!!  call deallocateCommuncationBuffersOrtho(lin%comon, subname)
!!
!!
!!end subroutine getOverlapMatrix


subroutine getOverlapMatrix2(iproc, nproc, lzd, orbs, comon_lb, op_lb, lphi, mad, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => getOverlapMatrix2
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  type(p2pCommsOrthonormality),intent(inout):: comon_lb
  type(overlapParameters),intent(inout):: op_lb
  real(8),dimension(orbs%npsidim_orbs),intent(inout):: lphi
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

  ! Local variables
  integer:: it, istat, iall, iorb, jorb, ierr
  real(8),dimension(:),allocatable:: lphiovrlp
  character(len=*),parameter:: subname='orthonormalize'
  logical:: converged
  real(8):: tt1, tt2, tt3
  type(input_variables):: input

  !!allocate(lphiovrlp(lin%op_lb%ndim_lphiovrlp), stat=istat)
  !!call memocc(istat, lphiovrlp, 'lphiovrlp',subname)

  call allocateCommuncationBuffersOrtho(comon_lb, subname)
  !call extractOrbital2(iproc,nproc,orbs,orbs%npsidim,orbs%inWhichLocreg,lzd,op_lb,lphi,comon_lb)
  call extractOrbital3(iproc,nproc,orbs,orbs%npsidim_orbs,orbs%inWhichLocreg,&
       lzd,op_lb,lphi,comon_lb%nsendBuf,comon_lb%sendBuf)
  !call postCommsOverlap(iproc,nproc,comon_lb)
  call postCommsOverlapNew(iproc,nproc,orbs,op_lb,lzd,lphi,comon_lb,tt1,tt2)
  !call gatherOrbitals2(iproc,nproc,comon_lb)
  call collectnew(iproc,nproc,comon_lb,mad,op_lb,orbs,input,lzd,comon_lb%nsendbuf,&
       comon_lb%sendbuf,comon_lb%nrecvbuf,comon_lb%recvbuf,ovrlp,tt1,tt2,tt3)
  !!call calculateOverlapMatrix2(iproc,nproc,orbs,op_lb,comon_lb,orbs%inWhichLocreg,mad,ovrlp)
  call calculateOverlapMatrix3(iproc, nproc, orbs, op_lb, orbs%inWhichLocreg, comon_lb%nsendBuf, &
       comon_lb%sendBuf, comon_lb%nrecvBuf, comon_lb%recvBuf, mad, ovrlp)
  call deallocateCommuncationBuffersOrtho(comon_lb, subname)




end subroutine getOverlapMatrix2




subroutine initCommsOrtho(iproc, nproc, lzd, orbs, onWhichAtomAll, input, locregShape, op, comon, tag)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initCommsOrtho
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
  type(input_variables),intent(in):: input
  character(len=1),intent(in):: locregShape
  type(overlapParameters),intent(out):: op
  type(p2pCommsOrthonormality),intent(out):: comon
  integer,intent(inout):: tag

  ! Local variables
  integer:: iorb, jorb, iiorb, jproc, ioverlaporb, ioverlapMPI, ilr, jlr
  integer:: ilrold, is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3
  integer::  je3, istat, i1, i2, irecv, isend, mpidest, mpisource
  logical:: ovrlpx, ovrlpy, ovrlpz
  character(len=*),parameter:: subname='initCommsOrtho'

  call nullify_overlapParameters(op)

  ! Allocate the arrays that count the number of overlaps per process (comon%noverlaps)
  ! and per orbital (op%noverlaps)
  allocate(comon%noverlaps(0:nproc-1), stat=istat)
  call memocc(istat, comon%noverlaps, 'comon%noverlaps',subname)
  allocate(op%noverlaps(orbs%norb), stat=istat)
  call memocc(istat, op%noverlaps, 'op%noverlaps',subname)

  ! Allocate the arrays holding the starting indices of the data to communicate in the
  ! send and receive buffers.
  allocate(op%indexInRecvBuf(orbs%norbp,orbs%norb), stat=istat)
  call memocc(istat, op%indexInRecvBuf, 'op%indexInRecvBuf', subname)
  allocate(op%indexInSendBuf(orbs%norbp,orbs%norb), stat=istat)
  call memocc(istat, op%indexInSendBuf, 'op%indexInSendBuf', subname)


  ! Count how many overlaping regions each orbital / process has.
  if(locregShape=='c') then
     call countOverlaps(iproc, nproc, orbs, lzd, onWhichAtomAll, op, comon)
     allocate(op%overlaps(maxval(op%noverlaps),orbs%norb), stat=istat)
     call memocc(istat, op%overlaps, 'op%overlaps', subname)
     allocate(comon%overlaps(maxval(comon%noverlaps),0:nproc-1), stat=istat)
     call memocc(istat, comon%overlaps, 'comon%overlaps', subname)
  else if(locregShape=='s') then
     call determine_overlap_from_descriptors(iproc, nproc, orbs, lzd, op, comon)
  end if

  ! Determine the overlapping orbitals.
  if(locregShape=='c') then
     call determineOverlaps(iproc, nproc, orbs, lzd, onWhichAtomAll, op, comon)
  else if(locregShape=='s') then
     ! already done in determine_overlap_from_descriptors
  end if

  ! Allocate the types describing the overlap localization regions.
  op%noverlapsmaxp=maxval(op%noverlaps(orbs%isorb+1:orbs%isorb+orbs%norbp))
  allocate(op%olr(op%noverlapsmaxp,orbs%norbp), stat=istat)
  do i2=1,orbs%norbp
     do i1=1,op%noverlapsmaxp
        call nullify_locreg_descriptors(op%olr(i1,i2))
     end do
  end do

  ! Set the orbital descriptors for the overlap regions.
  if(locregShape=='c') then
     call determineOverlapDescriptors(iproc, nproc, orbs, lzd, lzd%Glr, onWhichAtomAll, op)
  else if(locregShape=='s') then
     call determineOverlapDescriptorsSphere(iproc, nproc, orbs, lzd, lzd%Glr, onWhichAtomAll, input%hx, input%hy, input%hz, op)
  end if


  ! Initialize the communications.
  comon%noverlapsmax=maxval(comon%noverlaps)
  allocate(comon%comarr(10,comon%noverlapsmax,0:nproc-1), stat=istat)
  call memocc(istat, comon%comarr, 'comon%comarr', subname)
  allocate(comon%communComplete(comon%noverlapsmax,0:nproc-1), stat=istat)
  call memocc(istat, comon%communComplete, 'comun%communComplete', subname)
  call setCommsOrtho(iproc, nproc, orbs, onWhichAtomAll, lzd, op, comon, tag)


  ! Initialize the index arrays for the transformations from overlap region
  ! to ordinary localization region.
  allocate(op%indexExpand(comon%nrecvBuf), stat=istat)
  call memocc(istat, op%indexExpand, 'op%indexExpand',subname)

  allocate(op%expseg(op%noverlapsmaxp,orbs%norbp), stat=istat)
  do i1=1,orbs%norbp
      do i2=1,op%noverlapsmaxp
          call nullify_expansionSegments(op%expseg(i2,i1))
      end do
  end do
  call indicesForExpansion(iproc, nproc, orbs, input, onWhichAtomAll, lzd, op, comon)

  ! Initialize the index arrays for the transformations from the ordinary localization region
  ! to the overlap region.
  allocate(op%indexExtract(comon%nsendBuf), stat=istat)
  call memocc(istat, op%indexExtract, 'op%indexExtract',subname)
  allocate(op%extseg(op%noverlapsmaxp,orbs%norbp), stat=istat)
  do i1=1,orbs%norbp
      do i2=1,op%noverlapsmaxp
          call nullify_expansionSegments(op%extseg(i2,i1))
      end do
  end do
  call indicesForExtraction(iproc, nproc, orbs, orbs%npsidim_orbs, onWhichAtomAll, lzd, op, comon)

end subroutine initCommsOrtho


subroutine allocateCommuncationBuffersOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pCommsOrthonormality),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat

  allocate(comon%recvBuf(comon%nrecvBuf), stat=istat)
  call memocc(istat, comon%recvBuf, 'comon%recvBuf', subname)
  allocate(comon%sendBuf(comon%nsendBuf), stat=istat)
  call memocc(istat, comon%sendBuf, 'comon%sendBuf', subname)

end subroutine allocateCommuncationBuffersOrtho





subroutine deallocateCommuncationBuffersOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pCommsOrthonormality),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat, iall

  iall = -product(shape(comon%recvBuf))*kind(comon%recvBuf)
  deallocate(comon%recvBuf, stat=istat)
  call memocc(istat, iall, 'comon%recvBuf',subname)
  iall = -product(shape(comon%sendBuf))*kind(comon%sendBuf)
  deallocate(comon%sendBuf, stat=istat)
  call memocc(istat, iall, 'comon%sendBuf',subname)

end subroutine deallocateCommuncationBuffersOrtho



subroutine allocateSendBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pCommsOrthonormality),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat

  allocate(comon%sendBuf(comon%nsendBuf), stat=istat)
  call memocc(istat, comon%sendBuf, 'comon%sendBuf', subname)

end subroutine allocateSendBufferOrtho


subroutine deallocateSendBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pCommsOrthonormality),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat, iall

  iall = -product(shape(comon%sendBuf))*kind(comon%sendBuf)
  deallocate(comon%sendBuf, stat=istat)
  call memocc(istat, iall, 'comon%sendBuf',subname)

end subroutine deallocateSendBufferOrtho



subroutine allocateRecvBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pCommsOrthonormality),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat

  allocate(comon%recvBuf(comon%nrecvBuf), stat=istat)
  call memocc(istat, comon%recvBuf, 'comon%recvBuf', subname)

end subroutine allocateRecvBufferOrtho



subroutine deallocateRecvBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pCommsOrthonormality),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat, iall

  iall = -product(shape(comon%recvBuf))*kind(comon%recvBuf)
  deallocate(comon%recvBuf, stat=istat)
  call memocc(istat, iall, 'comon%recvBuf',subname)

end subroutine deallocateRecvBufferOrtho



!> Returns the starting and ending indices (on the coarse grid) of a given localization region.
subroutine getIndices(lr, is1, ie1, is2, ie2, is3, ie3)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(in):: lr
  integer,intent(out):: is1, ie1, is2, ie2, is3, ie3

  is1=lr%ns1+1
  ie1=lr%ns1+lr%d%n1
  is2=lr%ns2+1
  ie2=lr%ns2+lr%d%n2
  is3=lr%ns3+1
  ie3=lr%ns3+lr%d%n3

end subroutine getIndices




! Count for each orbital and each process the number of overlapping orbitals.
subroutine countOverlaps(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(overlapParameters),intent(out):: op
  type(p2pCommsOrthonormality),intent(out):: comon

  ! Local variables
  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
  integer:: js1, je1, js2, je2, js3, je3, iiorb
  logical:: ovrlpx, ovrlpy, ovrlpz

  iiorb=0
  do jproc=0,nproc-1
     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
     ilrold=-1
     do iorb=1,orbs%norb_par(jproc,0)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=onWhichAtom(iiorb)
        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
           jlr=onWhichAtom(jorb)
           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
           if(ovrlpx .and. ovrlpy .and. ovrlpz) then
              ioverlaporb=ioverlaporb+1
              if(ilr/=ilrold) then
                 ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
                 ! would get the same orbitals again. Therefore the counter is not increased
                 ! in that case.
                 ioverlapMPI=ioverlapMPI+1
              end if
           end if
        end do
        op%noverlaps(iiorb)=ioverlaporb
        !!if(iproc==0) write(*,'(a,2i8)') 'iiorb, op%noverlaps(iiorb)', iiorb, op%noverlaps(iiorb)
        ilrold=ilr
     end do
     comon%noverlaps(jproc)=ioverlapMPI
     !if(iproc==0) write(*,'(a,2i8)') 'jproc, comon%noverlaps(jproc)', jproc, comon%noverlaps(jproc)
  end do

end subroutine countOverlaps




! Count for each orbital and each process the number of overlapping orbitals.
subroutine countOverlapsSphere(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(overlapParameters),intent(out):: op
  type(p2pCommsOrthonormality),intent(out):: comon

  ! Local variables
  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
  integer:: js1, je1, js2, je2, js3, je3, iiorb
  logical:: ovrlpx, ovrlpy, ovrlpz
  real(8):: dx, dy, dz, rr

  iiorb=0
  do jproc=0,nproc-1
     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
     ilrold=-1
     do iorb=1,orbs%norb_par(jproc,0)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=onWhichAtom(iiorb)
        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
           jlr=onWhichAtom(jorb)
           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
           !if(ovrlpx .and. ovrlpy .and. ovrlpz) then
           dx=(lzd%llr(ilr)%locregCenter(1)-lzd%llr(jlr)%locregCenter(1))**2
           dy=(lzd%llr(ilr)%locregCenter(2)-lzd%llr(jlr)%locregCenter(2))**2
           dz=(lzd%llr(ilr)%locregCenter(3)-lzd%llr(jlr)%locregCenter(3))**2
           rr=(lzd%llr(ilr)%locrad+lzd%llr(jlr)%locrad)**2
           if(dx+dy+dz<=rr) then
              ioverlaporb=ioverlaporb+1
              if(ilr/=ilrold) then
                 ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
                 ! would get the same orbitals again. Therefore the counter is not increased
                 ! in that case.
                 ioverlapMPI=ioverlapMPI+1
              end if
           end if
        end do
        op%noverlaps(iiorb)=ioverlaporb
        !!if(iproc==0) write(*,'(a,2i8)') 'iiorb, op%noverlaps(iiorb)', iiorb, op%noverlaps(iiorb)
        ilrold=ilr
     end do
     comon%noverlaps(jproc)=ioverlapMpi
     !if(iproc==0) write(*,'(a,2i8)') 'jproc, comon%noverlaps(jproc)', jproc, comon%noverlaps(jproc)
  end do

end subroutine countOverlapsSphere





subroutine determineOverlaps(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(overlapParameters),intent(out):: op
  type(p2pCommsOrthonormality),intent(out):: comon

  ! Local variables
  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
  integer:: js1, je1, js2, je2, js3, je3, iiorb
  logical:: ovrlpx, ovrlpy, ovrlpz

  ! Initialize to some value which will never be used.
  op%overlaps=-1
  comon%overlaps=-1

  iiorb=0
  do jproc=0,nproc-1
     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
     ilrold=-1
     do iorb=1,orbs%norb_par(jproc,0)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=onWhichAtom(iiorb)
        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
           jlr=onWhichAtom(jorb)
           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
           if(ovrlpx .and. ovrlpy .and. ovrlpz) then
              ioverlaporb=ioverlaporb+1
              op%overlaps(ioverlaporb,iiorb)=jorb
              if(ilr/=ilrold) then
                 ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
                 ! would get the same orbitals again. Therefore the counter is not increased
                 ! in that case.
                 ioverlapMPI=ioverlapMPI+1
                 comon%overlaps(ioverlapMPI,jproc)=jorb
              end if
           end if
        end do
        !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
        ilrold=ilr
     end do
     !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
  end do

end subroutine determineOverlaps


subroutine determineOverlapsSphere(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(overlapParameters),intent(out):: op
  type(p2pCommsOrthonormality),intent(out):: comon

  ! Local variables
  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
  integer:: js1, je1, js2, je2, js3, je3, iiorb
  logical:: ovrlpx, ovrlpy, ovrlpz
  real(8):: dx, dy, dz, rr

  ! Initialize to some value which will never be used.
  op%overlaps=-1
  comon%overlaps=-1

  iiorb=0
  do jproc=0,nproc-1
     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
     ilrold=-1
     do iorb=1,orbs%norb_par(jproc,0)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=onWhichAtom(iiorb)
        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
           jlr=onWhichAtom(jorb)
           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
           dx=(lzd%llr(ilr)%locregCenter(1)-lzd%llr(jlr)%locregCenter(1))**2
           dy=(lzd%llr(ilr)%locregCenter(2)-lzd%llr(jlr)%locregCenter(2))**2
           dz=(lzd%llr(ilr)%locregCenter(3)-lzd%llr(jlr)%locregCenter(3))**2
           rr=(lzd%llr(ilr)%locrad+lzd%llr(jlr)%locrad)**2
           if(dx+dy+dz<=rr) then
              !if(ovrlpx .and. ovrlpy .and. ovrlpz) then
              ioverlaporb=ioverlaporb+1
              op%overlaps(ioverlaporb,iiorb)=jorb
              if(ilr/=ilrold) then
                 ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
                 ! would get the same orbitals again. Therefore the counter is not increased
                 ! in that case.
                 ioverlapMPI=ioverlapMPI+1
                 comon%overlaps(ioverlapMPI,jproc)=jorb
              end if
           end if
        end do
        !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
        ilrold=ilr
     end do
     !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
  end do

end subroutine determineOverlapsSphere




subroutine determineOverlapDescriptors(iproc, nproc, orbs, lzd, Glr, onWhichAtom, op)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(locreg_descriptors),intent(in):: Glr
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(overlapParameters),intent(inout):: op

  ! Local variables
  integer:: iorb, jorb, jjorb, ilr, jlr, iiorb


  do iorb=1,orbs%norbp
     iiorb=orbs%isorb_par(iproc)+iorb
     ilr=onWhichAtom(iiorb)
     !if(iproc==0) write(*,'(a,2i10)') 'iorb, op%noverlaps(iorb)', iorb, op%noverlaps(iorb)
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        jlr=onWhichAtom(jjorb)
        !write(*,*) 'calling get_overlap_region_periodic'
        !call get_overlap_region_periodic(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iorb))
        call get_overlap_region_periodic2(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iorb))
        !write(*,'(a,13i8)') 'iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, nvctr_c, nvctr_f, ncount, n1, n2, n3', iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, &
        !    op%olr(jorb,iorb)%wfd%nvctr_c, op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%d%n1, op%olr(jorb,iorb)%d%n2, op%olr(jorb,iorb)%d%n3
     end do
  end do

end subroutine determineOverlapDescriptors




subroutine determineOverlapDescriptorsSphere(iproc, nproc, orbs, lzd, Glr, onWhichAtom, hx, hy, hz, op)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(locreg_descriptors),intent(in):: Glr
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  real(8):: hx, hy, hz
  type(overlapParameters),intent(inout):: op

  ! Local variables
  integer:: iorb, jorb, jjorb, ilr, jlr, iiorb


  do iorb=1,orbs%norbp
     iiorb=orbs%isorb_par(iproc)+iorb
     ilr=onWhichAtom(iiorb)
     !if(iproc==0) write(*,'(a,2i10)') 'iorb, op%noverlaps(iorb)', iorb, op%noverlaps(iorb)
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        jlr=onWhichAtom(jjorb)
        call determine_overlapdescriptors_from_descriptors(lzd%llr(ilr), lzd%llr(jlr), lzd%glr, op%olr(jorb,iorb))
        !write(*,'(a,13i8)') 'iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, nvctr_c, nvctr_f, ncount, n1, n2, n3', iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, &
        !    op%olr(jorb,iorb)%wfd%nvctr_c, op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%d%n1, op%olr(jorb,iorb)%d%n2, op%olr(jorb,iorb)%d%n3
     end do
  end do

end subroutine determineOverlapDescriptorsSphere





subroutine setCommsOrtho(iproc, nproc, orbs, onWhichAtom, lzd, op, comon, tag)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(inout):: op
  type(p2pCommsOrthonormality),intent(out):: comon
  integer,intent(inout):: tag

  ! Local variables
  integer:: jproc, iorb, jorb, iiorb, jjorb, mpisource, mpidest, istsource, istdest, ncount, istat, iall, ijorb
  integer:: ilr, ilrold, jprocold, ildim, ierr, isend, irecv
  integer,dimension(:),allocatable:: istsourceArr, istdestArr
  character(len=*),parameter:: subname='setCommsOrtho'
  logical,dimension(:),allocatable:: receivedOrbital

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
        ilr=onWhichAtom(iiorb)
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
              ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
              !write(*,'(a,4i9)') 'iproc, iorb, jorb, ncount', iproc, iorb, jorb, ncount
           end if
           call mpi_bcast(ncount, 1, mpi_integer, jproc, mpi_comm_world, ierr)
           tag=tag+1
           receivedOrbital(jjorb)=.true.
           call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comon%comarr(1,ijorb,jproc))
           comon%comarr(9,ijorb,jproc)=jjorb
           comon%comarr(10,ijorb,jproc)=iiorb
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

end subroutine setCommsOrtho



subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: mpisource, mpidest, istsource, istdest, ncount, tag
  integer,dimension(8),intent(out):: comarr


  ! From which MPI process shall the orbital be sent.
  comarr(1)=mpisource

  ! Starting index on the sending process.
  comarr(2)=istsource

  ! Amount of datat to be sent
  comarr(3)=ncount

  ! To which MPI process shall the orbital be sent.
  comarr(4)=mpidest

  ! Starting index on the receiving process.
  comarr(5)=istdest

  ! Tag for the communication
  comarr(6)=tag

  ! comarr(7): this entry is used as request for the mpi_isend.

  ! comarr(8): this entry is used as request for the mpi_irecv.



end subroutine setCommsParameters




subroutine postCommsOverlap(iproc, nproc, comon)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => postCommsOverlap
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(p2pCommsOrthonormality),intent(inout):: comon

  ! Local variables
  integer:: jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, nsends, nreceives, ierr

!!! THIS IS THE ORIGINAL ##################################
  !nsends=0
  !nreceives=0
  !comon%communComplete=.false.
  !do jproc=0,nproc-1
  !    !write(*,'(3(a,i0))') 'iproc=',iproc,', jproc=',jproc,', comon%noverlaps(jproc)=', comon%noverlaps(jproc)
  !    do iorb=1,comon%noverlaps(jproc)
  !        mpisource=comon%comarr(1,iorb,jproc)
  !        istsource=comon%comarr(2,iorb,jproc)
  !        ncount=comon%comarr(3,iorb,jproc)
  !        mpidest=comon%comarr(4,iorb,jproc)
  !        istdest=comon%comarr(5,iorb,jproc)
  !        tag=comon%comarr(6,iorb,jproc)
  !        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
  !        if(mpisource/=mpidest) then
  !            ! The orbitals are on different processes, so we need a point to point communication.
  !            if(iproc==mpisource) then
  !                !write(*,'(6(a,i0))') 'overlap: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
  !                call mpi_isend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
  !                     mpi_comm_world, comon%comarr(7,iorb,jproc), ierr)
  !                !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
  !                comon%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
  !                nsends=nsends+1
  !            else if(iproc==mpidest) then
  !                !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
  !                call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
  !                     mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
  !                comon%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
  !                nreceives=nreceives+1
  !            else
  !                comon%comarr(7,iorb,jproc)=mpi_request_null
  !                comon%comarr(8,iorb,jproc)=mpi_request_null
  !            end if
  !        else
  !            ! The orbitals are on the same process, so simply copy them.
  !            if(iproc==mpisource) then
  !                call dcopy(ncount, comon%sendBuf(istsource), 1, comon%recvBuf(istdest), 1)
  !                !write(*,'(6(a,i0))') 'overlap: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
  !                comon%comarr(7,iorb,jproc)=mpi_request_null
  !                comon%comarr(8,iorb,jproc)=mpi_request_null
  !                nsends=nsends+1
  !                nreceives=nreceives+1
  !                comon%communComplete(iorb,mpisource)=.true.
  !            else
  !                comon%comarr(7,iorb,jproc)=mpi_request_null
  !                comon%comarr(8,iorb,jproc)=mpi_request_null
  !                comon%communComplete(iorb,mpisource)=.true.
  !            end if
  !
  !        end if
  !    end do
  !end do



  !! NEW ####################################
  ! First only post receives
  nsends=0
  nreceives=0
  comon%communComplete=.false.
  do jproc=0,nproc-1
     !write(*,'(3(a,i0))') 'iproc=',iproc,', jproc=',jproc,', comon%noverlaps(jproc)=', comon%noverlaps(jproc)
     do iorb=1,comon%noverlaps(jproc)
        mpisource=comon%comarr(1,iorb,jproc)
        istsource=comon%comarr(2,iorb,jproc)
        ncount=comon%comarr(3,iorb,jproc)
        mpidest=comon%comarr(4,iorb,jproc)
        istdest=comon%comarr(5,iorb,jproc)
        tag=comon%comarr(6,iorb,jproc)
        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
        if(mpisource/=mpidest) then
           ! The orbitals are on different processes, so we need a point to point communication.
           if(iproc==mpisource) then
           else if(iproc==mpidest) then
              !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount,&
              !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
              call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
                   mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
              comon%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
              nreceives=nreceives+1
           end if
        end if
     end do
  end do


  ! Now the rest.
  do jproc=0,nproc-1
     !write(*,'(3(a,i0))') 'iproc=',iproc,', jproc=',jproc,', comon%noverlaps(jproc)=', comon%noverlaps(jproc)
     do iorb=1,comon%noverlaps(jproc)
        mpisource=comon%comarr(1,iorb,jproc)
        istsource=comon%comarr(2,iorb,jproc)
        ncount=comon%comarr(3,iorb,jproc)
        mpidest=comon%comarr(4,iorb,jproc)
        istdest=comon%comarr(5,iorb,jproc)
        tag=comon%comarr(6,iorb,jproc)
        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
        if(mpisource/=mpidest) then
           ! The orbitals are on different processes, so we need a point to point communication.
           if(iproc==mpisource) then
              !write(*,'(6(a,i0))') 'overlap: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
              call mpi_isend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
                   mpi_comm_world, comon%comarr(7,iorb,jproc), ierr)
              !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
              comon%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
              nsends=nsends+1
           else if(iproc==mpidest) then
           else
              comon%comarr(7,iorb,jproc)=mpi_request_null
              comon%comarr(8,iorb,jproc)=mpi_request_null
           end if
        else
           ! The orbitals are on the same process, so simply copy them.
           if(iproc==mpisource) then
              call dcopy(ncount, comon%sendBuf(istsource), 1, comon%recvBuf(istdest), 1)
              !write(*,'(6(a,i0))') 'overlap: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
              comon%comarr(7,iorb,jproc)=mpi_request_null
              comon%comarr(8,iorb,jproc)=mpi_request_null
              nsends=nsends+1
              nreceives=nreceives+1
              comon%communComplete(iorb,mpisource)=.true.
           else
              comon%comarr(7,iorb,jproc)=mpi_request_null
              comon%comarr(8,iorb,jproc)=mpi_request_null
              comon%communComplete(iorb,mpisource)=.true.
           end if

        end if
     end do
  end do



end subroutine postCommsOverlap




subroutine extractOrbital3(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, phi, nsendBuf, sendBuf)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, sizePhi
  type(orbitals_data),intent(in):: orbs
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(inout):: op
  real(8),dimension(sizePhi),intent(in):: phi
  integer,intent(in):: nsendBuf
  real(8),dimension(nsendBuf),intent(out):: sendBuf


  ! Local variables
  integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
  integer:: i, indSource, m, ii, jst, istart, iend, ncount, iseg

  indovrlp=1
  op%indexInSendBuf=0

  ilrold=-1
  iiprocold=-1
  do iorb=1,orbs%norb
     ilr=onWhichAtom(iorb)
     iiproc=orbs%onWhichMPI(iorb)
     if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
     do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=onWhichAtom(jjorb)
        jjproc=orbs%onWhichMPI(jjorb)
        if(iproc==jjproc) then
           ! Get the correct descriptors
           korb=jjorb-orbs%isorb
           !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
           do i=1,op%noverlaps(jjorb)
              !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
              if(op%overlaps(i,jjorb)==iorb) then
                 lorb=i
                 exit
              end if
           end do
           !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
           gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
           ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
           ind=1
           do kkorb=orbs%isorb+1,jjorb-1
              klr=onWhichAtom(kkorb)
              ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
           end do
           !! THIS IS THE OLD VERSION
           !!do i=0,ldim-1
           !!    indSource=ind+op%indexExtract(indovrlp+i)-1
           !!    sendBuf(indovrlp+i)=phi(indSource)
           !!end do

            !! THIS IS THE NEWEST VERSION
            jst=0
            do iseg=1,op%extseg(lorb,korb)%nseg
                istart=op%extseg(lorb,korb)%segborders(1,iseg)
                iend=op%extseg(lorb,korb)%segborders(2,iseg)
                ncount=iend-istart+1
                !call dcopy(ncount,  comon%recvBuf(jst), 1, lphi(indout+istart-1), 1)
                call dcopy(ncount, phi(ind+istart-1), 1, sendBuf(indovrlp+jst), 1)
                jst=jst+ncount
           end do


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
           op%indexInSendBuf(jjorb-orbs%isorb,iorb)=indovrlp
           indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
        end if
     end do
     ilrold=ilr
     iiprocold=iiproc

  end do
!!do i=1,nsendBuf
!!     write(1000+iproc,*) i, sendbuf(i)
!!end do

  if(indovrlp/=nsendBuf+1) then
     write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=nsendBuf+1', indovrlp, nsendBuf+1
     stop
  end if

end subroutine extractOrbital3





!!subroutine gatherOrbitals(iproc, nproc, comon)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(p2pCommsOrthonormality),intent(inout):: comon
!!
!!! Local variables
!!integer:: kproc, mpisource, mpidest, nfast, nslow, nsameproc, ierr
!!integer,dimension(mpi_status_size):: stat
!!logical:: sendComplete, receiveComplete
!!
!!
!!! Check whether the communications have completed.
!!nfast=0
!!nsameproc=0
!!testLoop: do
!!    do kproc=1,comon%noverlaps(iproc)
!!        if(comon%communComplete(kproc,iproc)) cycle
!!        call mpi_test(comon%comarr(7,kproc,iproc), sendComplete, stat, ierr)
!!        call mpi_test(comon%comarr(8,kproc,iproc), receiveComplete, stat, ierr)
!!        if(sendComplete .and. receiveComplete) comon%communComplete(kproc,iproc)=.true.
!!        if(comon%communComplete(kproc,iproc)) then
!!            !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', kproc
!!            mpisource=comon%comarr(1,kproc,iproc)
!!            mpidest=comon%comarr(4,kproc,iproc)
!!            if(mpisource/=mpidest) then
!!                nfast=nfast+1
!!            else
!!                nsameproc=nsameproc+1
!!            end if
!!        end if
!!    end do
!!    ! If we made it until here, either all all the communication is
!!    ! complete or we better wait for each single orbital.
!!    exit testLoop
!!end do testLoop
!!
!!
!!! Wait for the communications that have not completed yet
!!nslow=0
!!do kproc=1,comon%noverlaps(iproc)
!!    if(comon%communComplete(kproc,iproc)) cycle
!!    !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', kproc
!!    nslow=nslow+1
!!    call mpi_wait(comon%comarr(7,kproc,iproc), stat, ierr)   !COMMENTED BY PB
!!    call mpi_wait(comon%comarr(8,kproc,iproc), stat, ierr)   !COMMENTED BY PB
!!    comon%communComplete(kproc,iproc)=.true.
!!end do
!!
!!!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!                       nfast, ' could be overlapped with computation.'
!!!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!
!!
!!
!!end subroutine gatherOrbitals





subroutine gatherOrbitals2(iproc, nproc, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(p2pCommsOrthonormality),intent(inout):: comon

  ! Local variables
  integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
  integer,dimension(mpi_status_size):: stat
  logical:: sendComplete, receiveComplete


!!! Check whether the communications have completed.
  !!nfast=0
  !!nsameproc=0
  !!testLoop: do
  !!    do jproc=0,nproc-1
  !!        do jorb=1,comon%noverlaps(jproc)
  !!            if(comon%communComplete(jorb,jproc)) cycle
  !!            call mpi_test(comon%comarr(7,jorb,jproc), sendComplete, stat, ierr)
  !!            call mpi_test(comon%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
  !!            ! Attention: mpi_test is a local function.
  !!            if(sendComplete .and. receiveComplete) comon%communComplete(jorb,jproc)=.true.
  !!        end do
  !!    end do
  !!    ! If we made it until here, either all all the communication is
  !!    ! complete or we better wait for each single orbital.
  !!    exit testLoop
  !!end do testLoop
  !!
!!! Since mpi_test is a local function, check whether the communication has completed on all processes.
  !!call mpiallred(comon%communComplete(1,0), nproc*maxval(comon%noverlaps), mpi_land, mpi_comm_world, ierr)


  ! Wait for the communications that have not completed yet
  nslow=0
  do jproc=0,nproc-1
     do jorb=1,comon%noverlaps(jproc)
        if(comon%communComplete(jorb,jproc)) then
           mpisource=comon%comarr(1,jorb,jproc)
           mpidest=comon%comarr(4,jorb,jproc)
           if(mpisource==mpidest) then
              !nsameproc=nsameproc+1
           else
              !nfast=nfast+1
           end if
           cycle
        end if
        !write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
        nslow=nslow+1
        !if(iproc==mpidest) call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        !if(iproc==mpisource) call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        comon%communComplete(jorb,jproc)=.true.
     end do
  end do

  !call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
  !nfast=nint(dble(nfast)/dble(nproc))
  !nfast=nint(dble(nslow)/dble(nproc))
  !nfast=nint(dble(nsameproc)/dble(nproc))
  !if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
  !                       nfast, ' could be overlapped with computation.'
  !if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


end subroutine gatherOrbitals2





subroutine gatherOrbitalsOverlapWithComput(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp, expanded)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => gatherOrbitalsOverlapWithComput
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(input_variables),intent(in):: input
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  type(p2pCommsOrthonormality),intent(inout):: comon
  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
  logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded

  ! Local variables
  integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc, jjorb, orbsource, orbdest
  integer,dimension(mpi_status_size):: stat
  logical:: sendComplete, receiveComplete


  ! Check whether the communications have completed.
  nfast=0
  nsameproc=0
  testLoop: do
     do jproc=0,nproc-1
        do jorb=1,comon%noverlaps(jproc)
           mpisource=comon%comarr(1,jorb,jproc)
           mpidest=comon%comarr(4,jorb,jproc)
           orbsource=comon%comarr(9,jorb,jproc)
           orbdest=comon%comarr(10,jorb,jproc)
           if(mpisource==mpidest) then
              if(iproc==mpidest) then
                 call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
                      orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
                 expanded(orbsource,orbdest-orbs%isorb)=.true.
              end if
           end if
           if(comon%communComplete(jorb,jproc)) cycle
           ! Attention: mpi_test is a local function.
           call mpi_test(comon%comarr(7,jorb,jproc), sendComplete, stat, ierr)
           call mpi_test(comon%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
           if(sendComplete .and. receiveComplete) comon%communComplete(jorb,jproc)=.true.
           if(comon%communComplete(jorb,jproc)) then
              !if(iproc==jproc) write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
              mpisource=comon%comarr(1,jorb,jproc)
              mpidest=comon%comarr(4,jorb,jproc)
              orbsource=comon%comarr(9,jorb,jproc)
              orbdest=comon%comarr(10,jorb,jproc)
              if(iproc==mpidest) then
                 call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
                      orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
                 expanded(orbsource,orbdest-orbs%isorb)=.true.
              end if
              if(mpisource/=mpidest) then
                 !nfast=nfast+1
              else
                 !nsameproc=nsameproc+1
              end if
           end if
        end do
     end do
     ! If we made it until here, either all all the communication is
     ! complete or we better wait for each single orbital.
     exit testLoop
  end do testLoop

  ! Since mpi_test is a local function, check whether the communication has completed on all processes.
  call mpiallred(comon%communComplete(1,0), nproc*maxval(comon%noverlaps), mpi_land, mpi_comm_world, ierr)

  ! Wait for the communications that have not completed yet
  nslow=0
  do jproc=0,nproc-1
     do jorb=1,comon%noverlaps(jproc)
        !!if(comon%communComplete(jorb,jproc)) cycle
        if(comon%communComplete(jorb,jproc)) then
           mpisource=comon%comarr(1,jorb,jproc)
           mpidest=comon%comarr(4,jorb,jproc)
           if(mpisource==mpidest) then
              nsameproc=nsameproc+1
           else
              nfast=nfast+1
           end if
           cycle
        end if
        !write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
        nslow=nslow+1
        call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        comon%communComplete(jorb,jproc)=.true.
        mpidest=comon%comarr(4,jorb,jproc)
        orbsource=comon%comarr(9,jorb,jproc)
        orbdest=comon%comarr(10,jorb,jproc)
        if(iproc==mpidest) then
           !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
           expanded(orbsource,orbdest-orbs%isorb)=.false.
        end if






!!!write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
        !!nslow=nslow+1
        !!call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        !!call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        !!comon%communComplete(jorb,jproc)=.true.
        !!mpisource=comon%comarr(1,jorb,jproc)
        !!mpidest=comon%comarr(4,jorb,jproc)
        !!orbsource=comon%comarr(9,jorb,jproc)
        !!orbdest=comon%comarr(10,jorb,jproc)
        !!if(iproc==mpidest) then
        !!    !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
        !!    expanded(orbsource,orbdest-orbs%isorb)=.false.
        !!end if
!!!write(*,'(3(a,i0))') 'process ', iproc, ' has finally received orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
     end do
  end do

  !call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
  if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
       nfast, ' could be overlapped with computation.'
  if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


end subroutine gatherOrbitalsOverlapWithComput







subroutine gatherOrbitalsOverlapWithComput2(iproc, nproc, orbs, input, lzd, op, comon,&
     nsendbuf, sendbuf, nrecvbuf, recvbuf, lphiovrlp, expanded, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => gatherOrbitalsOverlapWithComput2
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
  type(orbitals_data),intent(in):: orbs
  type(input_variables),intent(in):: input
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  type(p2pCommsOrthonormality),intent(inout):: comon
  real(8),dimension(nsendbuf),intent(in):: sendbuf
  real(8),dimension(nrecvbuf),intent(in):: recvbuf
  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
  logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

  ! Local variables
  integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc, jjorb, orbsource, orbdest, ncount, ist, jst
  integer,dimension(mpi_status_size):: stat
  logical:: sendComplete, receiveComplete
  real(8):: ddot

  if(iproc==0) write(*,*) 'new subroutine'

  ! Check whether the communications have completed. Only check the receives.
  nfast=0
  nsameproc=0
  testLoop: do
     do jproc=0,nproc-1
        do jorb=1,comon%noverlaps(jproc)
           mpisource=comon%comarr(1,jorb,jproc)
           mpidest=comon%comarr(4,jorb,jproc)
           orbsource=comon%comarr(9,jorb,jproc)
           orbdest=comon%comarr(10,jorb,jproc)
           if(iproc==mpidest) then
              if(mpisource==mpidest) then
                 call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
                      orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
                 expanded(orbsource,orbdest-orbs%isorb)=.true.
                 ! Calculate matrix element here.
                 call getStartingIndicesGlobal(orbsource, orbdest, op, orbs, ist, jst, ncount)
                 ovrlp(orbdest,orbsource)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
              end if
              if(comon%communComplete(jorb,jproc)) cycle
              ! Attention: mpi_test is a local function.
              ! Test whether the receive has completed.
              !call mpi_test(comon%comarr(7,jorb,jproc), sendComplete, stat, ierr)
              call mpi_test(comon%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
              !if(sendComplete .and. receiveComplete) comon%communComplete(jorb,jproc)=.true.
              if(receiveComplete) comon%communComplete(jorb,jproc)=.true.
              if(comon%communComplete(jorb,jproc)) then
                 !if(iproc==jproc) write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
                 mpisource=comon%comarr(1,jorb,jproc)
                 mpidest=comon%comarr(4,jorb,jproc)
                 orbsource=comon%comarr(9,jorb,jproc)
                 orbdest=comon%comarr(10,jorb,jproc)
                 if(iproc==mpidest) then
                    call expandOneOrbital(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
                         orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
                    expanded(orbsource,orbdest-orbs%isorb)=.true.
                    ! Calculate matrix element here.
                    call getStartingIndicesGlobal(orbsource, orbdest, op, orbs, ist, jst, ncount)
                    ovrlp(orbdest,orbsource)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
                 end if
                 if(mpisource/=mpidest) then
                    !nfast=nfast+1
                 else
                    !nsameproc=nsameproc+1
                 end if
              end if
           end if
        end do
     end do
     ! If we made it until here, either all all the communication is
     ! complete or we better wait for each single orbital.
     exit testLoop
  end do testLoop

  !!call mpi_barrier(mpi_comm_world, ierr)
  !!if(iproc==0) write(*,*) 'after test loop'
  !! Since mpi_test is a local function, check whether the communication has completed on all processes.
  !call mpiallred(comon%communComplete(1,0), nproc*maxval(comon%noverlaps), mpi_land, mpi_comm_world, ierr)

  ! Wait for the communications that have not completed yet
  nslow=0
  do jproc=0,nproc-1
     do jorb=1,comon%noverlaps(jproc)
        mpisource=comon%comarr(1,jorb,jproc)
        mpidest=comon%comarr(4,jorb,jproc)
        if(iproc==mpidest) then
           ! Only check the receive
           !!if(comon%communComplete(jorb,jproc)) cycle
           if(comon%communComplete(jorb,jproc)) then
              if(mpisource==mpidest) then
                 nsameproc=nsameproc+1
              else
                 nfast=nfast+1
              end if
              cycle
           end if
           !write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
           nslow=nslow+1
           !call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
           call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
           comon%communComplete(jorb,jproc)=.true.
           mpidest=comon%comarr(4,jorb,jproc)
           orbsource=comon%comarr(9,jorb,jproc)
           orbdest=comon%comarr(10,jorb,jproc)
           if(iproc==mpidest) then
              !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
              expanded(orbsource,orbdest-orbs%isorb)=.false.
              call getStartingIndicesGlobal(orbsource, orbdest, op, orbs, ist, jst, ncount)
              ovrlp(orbdest,orbsource)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
           end if

        else if(iproc==mpisource) then
           ! Check the send
           call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)
        end if






!!!write(*,'(3(a,i0))') 'process ', iproc, ' is waiting for orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
        !!nslow=nslow+1
        !!call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        !!call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        !!comon%communComplete(jorb,jproc)=.true.
        !!mpisource=comon%comarr(1,jorb,jproc)
        !!mpidest=comon%comarr(4,jorb,jproc)
        !!orbsource=comon%comarr(9,jorb,jproc)
        !!orbdest=comon%comarr(10,jorb,jproc)
        !!if(iproc==mpidest) then
        !!    !call expandOneOrbital(iproc, nproc, jjorb, orbs, input, orbs%inWhichLocreg, lzd, op, comon, lphiovrlp)
        !!    expanded(orbsource,orbdest-orbs%isorb)=.false.
        !!end if
!!!write(*,'(3(a,i0))') 'process ', iproc, ' has finally received orbital ',jorb,'; tag=',comon%comarr(6,jorb,jproc)
     end do
  end do

  !call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
  !call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
  if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
       nfast, ' could be overlapped with computation.'
  if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


end subroutine gatherOrbitalsOverlapWithComput2









subroutine calculateOverlapMatrix3(iproc, nproc, orbs, op, onWhichAtom, nsendBuf, sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(in):: op
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  real(8),dimension(nsendBuf),intent(in):: sendBuf
  real(8),dimension(nrecvBuf),intent(in):: recvBuf
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

  ! Local variables
  integer:: iorb, jorb, iiorb, jjorb, jjproc, ii, ist, jst, jjlr, ncount, ierr, jj, korb, kkorb, iilr, iiproc, istat, iall
  real(8):: ddot, tt, ttmax
  real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2
  character(len=*),parameter:: subname='calculateOverlapMatrix3'
  integer,dimension(:),allocatable:: sendcounts, displs


  ovrlp=0.d0

  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
!!! Starting index of orbital iorb, already transformed to overlap region with jjorb.
!!! We have to find the first orbital on the same MPI and in the same locreg as jjorb.
        !!jjlr=onWhichAtom(jjorb)
        !!jjproc=orbs%onWhichMPI(jjorb)
        !!jj=orbs%isorb_par(jjproc)+1
        !!do
        !!    if(onWhichAtom(jj)==jjlr) exit
        !!    jj=jj+1
        !!end do
        !!ist=op%indexInSendBuf(iorb,jj)
!!! Starting index of orbital jjorb.
!!! We have to find the first orbital on the same MPI and in the same locreg as iiorb.
        !!iilr=onWhichAtom(iiorb)
        !!iiproc=orbs%onWhichMPI(iiorb)
        !!do korb=1,orbs%norbp
        !!    kkorb=orbs%isorb_par(iiproc)+korb
        !!    if(onWhichAtom(kkorb)==iilr) then
        !!        ii=korb
        !!        exit
        !!    end if
        !!end do
        !!jst=op%indexInRecvBuf(ii,jjorb)
        call getStartingIndices(iorb, jorb, op, orbs, ist, jst)
        !write(*,'(5(a,i0))') 'process ',iproc,' calculates overlap of ',iiorb,' and ',jjorb,'. ist=',ist,' jst=',jst 
        !ncount=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !write(*,'(a,4i8)') 'iproc, iiorb, jjorb, ncount', iproc, iiorb, jjorb, ncount
        !ovrlp(iiorb,jjorb)=ddot(ncount, comon%sendBuf(ist), 1, comon%recvBuf(jst), 1)
        !ovrlp(iiorb,jjorb)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
        ovrlp(jjorb,iiorb)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
     end do
  end do


  allocate(ovrlpCompressed(mad%nvctr), stat=istat)
  call memocc(istat, ovrlpCompressed, 'ovrlpCompressed', subname)

  !call compressMatrix(orbs%norb, mad, ovrlp, ovrlpCompressed)
  allocate(sendcounts(0:nproc-1), stat=istat)
  call memocc(istat, sendcounts, 'sendcounts', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed, sendcounts, displs)
  !call mpiallred(ovrlpCompressed(1), mad%nvctr, mpi_sum, mpi_comm_world, ierr)
  allocate(ovrlpCompressed2(mad%nvctr), stat=istat)
  call memocc(istat, ovrlpCompressed2, 'ovrlpCompressed2', subname)
  if (nproc >1) then
     call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
          sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
  else
     call vcopy(sendcounts(iproc),ovrlpCompressed(displs(iproc)+1),1,ovrlpCompressed2(1+displs(iproc)),1)
  end if

  !call uncompressMatrix(orbs%norb, mad, ovrlpCompressed, ovrlp)
  call uncompressMatrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)

  iall=-product(shape(ovrlpCompressed))*kind(ovrlpCompressed)
  deallocate(ovrlpCompressed, stat=istat)
  call memocc(istat, iall, 'ovrlpCompressed', subname)
  iall=-product(shape(ovrlpCompressed2))*kind(ovrlpCompressed2)
  deallocate(ovrlpCompressed2, stat=istat)
  call memocc(istat, iall, 'ovrlpCompressed2', subname)
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


  !!call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)

end subroutine calculateOverlapMatrix3






subroutine calculateOverlapMatrix3Partial(iproc, nproc, orbs, op, onWhichAtom, nsendBuf, &
     sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(in):: op
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  real(8),dimension(nsendBuf),intent(in):: sendBuf
  real(8),dimension(nrecvBuf),intent(in):: recvBuf
  type(matrixDescriptors),intent(in):: mad
  !logical,dimension(0:nproc-1),intent(in):: skip
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

  ! Local variables
  integer:: iorb, jorb, iiorb, jjorb, jjproc, ii, ist, jst, jjlr, ncount, ierr, jj, korb, kkorb, iilr, iiproc, istat, iall
  real(8):: ddot, tt, ttmax
  real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2
  character(len=*),parameter:: subname='calculateOverlapMatrix3'
  integer,dimension(:),allocatable:: sendcounts, displs


  ovrlp=0.d0

  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
!!! Starting index of orbital iorb, already transformed to overlap region with jjorb.
!!! We have to find the first orbital on the same MPI and in the same locreg as jjorb.
        !!jjlr=onWhichAtom(jjorb)
        !!jjproc=orbs%onWhichMPI(jjorb)
        !!jj=orbs%isorb_par(jjproc)+1
        !!do
        !!    if(onWhichAtom(jj)==jjlr) exit
        !!    jj=jj+1
        !!end do
        !!ist=op%indexInSendBuf(iorb,jj)
!!! Starting index of orbital jjorb.
!!! We have to find the first orbital on the same MPI and in the same locreg as iiorb.
        !!iilr=onWhichAtom(iiorb)
        !!iiproc=orbs%onWhichMPI(iiorb)
        !!do korb=1,orbs%norbp
        !!    kkorb=orbs%isorb_par(iiproc)+korb
        !!    if(onWhichAtom(kkorb)==iilr) then
        !!        ii=korb
        !!        exit
        !!    end if
        !!end do
        !!jst=op%indexInRecvBuf(ii,jjorb)
        call getStartingIndices(iorb, jorb, op, orbs, ist, jst)
        !write(*,'(5(a,i0))') 'process ',iproc,' calculates overlap of ',iiorb,' and ',jjorb,'. ist=',ist,' jst=',jst 
        !ncount=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !write(*,'(a,4i8)') 'iproc, iiorb, jjorb, ncount', iproc, iiorb, jjorb, ncount
        !ovrlp(iiorb,jjorb)=ddot(ncount, comon%sendBuf(ist), 1, comon%recvBuf(jst), 1)
        !ovrlp(iiorb,jjorb)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
        ovrlp(jjorb,iiorb)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
     end do
  end do


  !!allocate(ovrlpCompressed(mad%nvctr), stat=istat)
  !!call memocc(istat, ovrlpCompressed, 'ovrlpCompressed', subname)
  !!
!!!call compressMatrix(orbs%norb, mad, ovrlp, ovrlpCompressed)
  !!allocate(sendcounts(0:nproc-1), stat=istat)
  !!call memocc(istat, sendcounts, 'sendcounts', subname)
  !!allocate(displs(0:nproc-1), stat=istat)
  !!call memocc(istat, displs, 'displs', subname)
  !!call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed, sendcounts, displs)
!!!call mpiallred(ovrlpCompressed(1), mad%nvctr, mpi_sum, mpi_comm_world, ierr)
  !!allocate(ovrlpCompressed2(mad%nvctr), stat=istat)
  !!call memocc(istat, ovrlpCompressed2, 'ovrlpCompressed2', subname)
  !!call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
  !!     sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
!!!call uncompressMatrix(orbs%norb, mad, ovrlpCompressed, ovrlp)
  !!call uncompressMatrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)
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

  !ttmax=0.d0
  !do iorb=1,orbs%norb
  !    do jorb=iorb,orbs%norb
  !        tt=abs(ovrlp(jorb,iorb)-ovrlp(iorb,jorb))
  !        if(tt>ttmax) ttmax=tt
  !    end do
  !end do
  !if(iproc==0) write(*,*) 'in calculateOverlapMatrix3: max dev from symmetry:', ttmax


  !!call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)

end subroutine calculateOverlapMatrix3Partial






subroutine expandOrbital2(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(input_variables),intent(in):: input
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  type(p2pCommsOrthonormality),intent(in):: comon
  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp

  ! Local variables
  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest, m, ii
!!real(8):: t1, t2, time1, time2

!!t1=mpi_wtime()
  lphiovrlp=0.d0
!!t2=mpi_wtime()
!!time1=t2-t1

!!t1=mpi_wtime()
  ind=1
  ilrold=-1
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr=onWhichAtom(iiorb)
     if(ilr==ilrold) cycle
     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !! THIS IS THE OLD VERSION
        !!do i=0,ldim-1
        !!    indDest=ind+op%indexExpand(jst+i)-1
        !!    lphiovrlp(indDest)=comon%recvBuf(jst+i)
        !!end do
        !! THIS IS NEW
        m=mod(ldim,7)
        if(m/=0) then
           do i=0,m-1
                ii=jst+i
                indDest=ind+op%indexExpand(ii)-1
                lphiovrlp(indDest)=comon%recvBuf(ii)
           end do
        end if
        do i=m,ldim-1,7
                ii=jst+i
                indDest=ind+op%indexExpand(ii+0)-1
                !if(iproc==0) write(333,*) indDest
                lphiovrlp(indDest)=comon%recvBuf(ii+0)
                indDest=ind+op%indexExpand(ii+1)-1
                !if(iproc==0) write(333,*) indDest
                lphiovrlp(indDest)=comon%recvBuf(ii+1)
                indDest=ind+op%indexExpand(ii+2)-1
                !if(iproc==0) write(333,*) indDest
                lphiovrlp(indDest)=comon%recvBuf(ii+2)
                indDest=ind+op%indexExpand(ii+3)-1
                !if(iproc==0) write(333,*) indDest
                lphiovrlp(indDest)=comon%recvBuf(ii+3)
                indDest=ind+op%indexExpand(ii+4)-1
                !if(iproc==0) write(333,*) indDest
                lphiovrlp(indDest)=comon%recvBuf(ii+4)
                indDest=ind+op%indexExpand(ii+5)-1
                !if(iproc==0) write(333,*) indDest
                lphiovrlp(indDest)=comon%recvBuf(ii+5)
                indDest=ind+op%indexExpand(ii+6)-1
                !if(iproc==0) write(333,*) indDest
                lphiovrlp(indDest)=comon%recvBuf(ii+6)
        end do
        ind=ind+gdim
     end do
     ilrold=ilr
  end do
!!t2=mpi_wtime()
!!time2=t2-t1
!!if(iproc==0) write(*,*) 'time1, time2', time1, time2

end subroutine expandOrbital2



subroutine expandRemainingOrbitals(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, expanded, lphiovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(input_variables),intent(in):: input
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  type(p2pCommsOrthonormality),intent(in):: comon
  logical,dimension(orbs%norb,orbs%norbp),intent(in):: expanded
  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp

  ! Local variables
  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest


  !lphiovrlp=0.d0

  ind=1
  ilrold=-1
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr=onWhichAtom(iiorb)
     if(ilr==ilrold) cycle
     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
        if(.not. expanded(jjorb,iorb)) then
           do i=0,ldim-1
              indDest=ind+op%indexExpand(jst+i)-1
              lphiovrlp(indDest)=comon%recvBuf(jst+i)
           end do
        end if
        ind=ind+gdim
     end do
     ilrold=ilr
  end do

end subroutine expandRemainingOrbitals



subroutine expandOneOrbital(iproc, nproc, orbsource, orbdest, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, orbsource, orbdest
  type(orbitals_data),intent(in):: orbs
  type(input_variables),intent(in):: input
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  type(p2pCommsOrthonormality),intent(in):: comon
  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp

  ! Local variables
  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest


  !lphiovrlp=0.d0

  ind=1
  ilrold=-1
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr=onWhichAtom(iiorb)
     if(ilr==ilrold) cycle
     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
        if(jjorb==orbsource .and. iorb==orbdest) then
           do i=0,ldim-1
              indDest=ind+op%indexExpand(jst+i)-1
              lphiovrlp(indDest)=comon%recvBuf(jst+i)
           end do
        end if
        ind=ind+gdim
     end do
     ilrold=ilr
  end do

end subroutine expandOneOrbital



subroutine expandOneOrbital2(iproc, nproc, orbsource, orbdest, orbs, input, onWhichAtom, lzd, op, nrecvbuf, recvbuf, lphiovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, orbsource, orbdest, nrecvbuf
  type(orbitals_data),intent(in):: orbs
  type(input_variables),intent(in):: input
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  real(8),dimension(nrecvbuf),intent(in):: recvbuf
  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp

  ! Local variables
  integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest


  ind=1
  ilrold=-1
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr=onWhichAtom(iiorb)
     if(ilr==ilrold) cycle
     gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
     do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        if(jjorb==orbsource .and. iorb==orbdest) then
           do i=0,ldim-1
              indDest=ind+op%indexExpand(jst+i)-1
              lphiovrlp(indDest)=recvbuf(jst+i)
              !lphiovrlp(indDest)=recvbuf(1+i)
           end do
        end if
        ind=ind+gdim
     end do
     ilrold=ilr
  end do

end subroutine expandOneOrbital2








subroutine checkUnity(iproc, norb, ovrlp, maxError)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: ovrlp
  real(8),intent(out):: maxError

  ! Local variables
  integer:: iorb, jorb
  real(8):: error

  maxError=0.d0
  do iorb=1,norb
     do jorb=1,norb
        if(iorb==jorb) then
           error=abs(ovrlp(jorb,iorb)-1.d0)
        else
           error=abs(ovrlp(jorb,iorb))
        end if
        if(error>maxError) then
           maxError=error
        end if
          write(20000+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
     end do
  end do

end subroutine checkUnity




subroutine globalLoewdin(iproc, nproc, orbs, lorbs, onWhichAtom, lzd, op, comon, ovrlp, lphiovrlp, lphi)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs, lorbs
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(in):: comon
  real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
  real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
  real(8),dimension(lorbs%npsidim_orbs),intent(out):: lphi

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount
  real(8):: tt, dnrm2


call build_new_linear_combinations(lzd, orbs, op, comon%nrecvbuf, comon%recvbuf, ovrlp, .true., lphi)

! Normalize
  ist=1
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr=onWhichAtom(iiorb)
     ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f

     ! Normalize
     tt=dnrm2(ncount, lphi(ist), 1)
     call dscal(ncount, 1/tt, lphi(ist), 1)

     ist=ist+ncount
  end do

!!lphi=0.d0
!!
!!ist=1
!!jst=1
!!ilrold=-1
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=onWhichAtom(iiorb)
!!    if(ilr==ilrold) then
!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!        jst=jst-op%noverlaps(iiorb)*ncount
!!    end if
!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        call daxpy(ncount, ovrlp(jjorb,iiorb), lphiovrlp(jst), 1, lphi(ist), 1)
!!        jst=jst+ncount
!!    end do
!!
!!    ! Normalize
!!    tt=dnrm2(ncount, lphi(ist), 1)
!!    call dscal(ncount, 1/tt, lphi(ist), 1)
!!
!!    ist=ist+ncount
!!    ilrold=ilr
!!
!!end do


end subroutine globalLoewdin



!!!!!!
!!!!!!subroutine gramschmidt_forlinear(iproc, nproc, iorbwork, orbs, lorbs, onWhichAtom, lzd, op, ovrlp, lphiovrlp, lphi)
!!!!!!use module_base
!!!!!!use module_types
!!!!!!implicit none
!!!!!!
!!!!!!! Calling arguments
!!!!!!integer,intent(in):: iproc, nproc, iorbwork
!!!!!!type(orbitals_data),intent(in):: orbs, lorbs
!!!!!!integer,dimension(orbs%norb),intent(in):: onWhichAtom
!!!!!!type(local_zone_descriptors),intent(in):: lzd
!!!!!!type(overlapParameters),intent(in):: op
!!!!!!real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
!!!!!!real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
!!!!!!real(8),dimension(lorbs%npsidim),intent(out):: lphi
!!!!!!
!!!!!!! Local variables
!!!!!!integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, jst2
!!!!!!real(8):: tt, dnrm2, ddot
!!!!!!
!!!!!!!lphi=0.d0
!!!!!!
!!!!!!ist=1
!!!!!!jst=1
!!!!!!jst2=1
!!!!!!ilrold=-1
!!!!!!do iorb=1,orbs%norbp
!!!!!!    iiorb=orbs%isorb+iorb
!!!!!!    ilr=onWhichAtom(iiorb)
!!!!!!    if(ilr==ilrold) then
!!!!!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!!!!!        jst=jst-op%noverlaps(iiorb)*ncount
!!!!!!        jst2=jst2-op%noverlaps(iiorb)*ncount
!!!!!!    end if
!!!!!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!!!!!    if(iiorb==iorbwork) then
!!!!!!        !write(*,*) 'before: sum(lphi)',sum(lphi)
!!!!!!        do jorb=1,op%noverlaps(iiorb)
!!!!!!            jjorb=op%overlaps(jorb,iiorb)
!!!!!!            if(jjorb<iiorb) then
!!!!!!                !write(*,'(3(a,i0),a,es15.8,a,es15.8)') 'iproc ',iproc,' works on orbital ', iiorb,'; subtracting orbital ',jjorb,'. ddot before:',&
!!!!!!                !    ddot(ncount, lphiovrlp(jst), 1, lphi(ist), 1), '    ddot before2', ovrlp(jjorb,iiorb)
!!!!!!                !write(*,'(a,2i11,2es15.6)') 'ist, jst, ddot(lphi), ddot(lphiovrlp)', ist, jst, ddot(ncount, lphi(ist), 1, lphi(ist), 1), &
!!!!!!                !    ddot(ncount, lphiovrlp(jst), 1, lphiovrlp(jst), 1)
!!!!!!                tt=ddot(ncount, lphiovrlp(jst), 1, lphi(ist), 1)/ddot(ncount, lphiovrlp(jst), 1, lphiovrlp(jst), 1)
!!!!!!                call daxpy(ncount, -tt, lphiovrlp(jst), 1, lphi(ist), 1)
!!!!!!                write(*,'(a,es15.8)') 'ddot :after', ddot(ncount, lphiovrlp(jst), 1, lphi(ist), 1)
!!!!!!            end if
!!!!!!            jst=jst+ncount
!!!!!!        end do
!!!!!!        do jorb=1,op%noverlaps(iiorb)
!!!!!!            jjorb=op%overlaps(jorb,iiorb)
!!!!!!            if(jjorb<iiorb) then
!!!!!!                write(*,'(a,es15.8)') 'overview: ddot :after', ddot(ncount, lphiovrlp(jst2), 1, lphi(ist), 1)
!!!!!!            end if
!!!!!!            jst2=jst2+ncount
!!!!!!        end do
!!!!!!        !write(*,*) 'after1: sum(lphi)',sum(lphi)
!!!!!!
!!!!!!        ! Normalize
!!!!!!        tt=dnrm2(ncount, lphi(ist), 1)
!!!!!!        !call dscal(ncount, 1/tt, lphi(ist), 1)
!!!!!!        !write(*,*) 'after2: sum(lphi)',sum(lphi)
!!!!!!    else
!!!!!!        ! only increase index
!!!!!!        do jorb=1,op%noverlaps(iiorb)
!!!!!!            jjorb=op%overlaps(jorb,iiorb)
!!!!!!            jst=jst+ncount
!!!!!!            jst2=jst2+ncount
!!!!!!        end do
!!!!!!    end if
!!!!!!
!!!!!!    ist=ist+ncount
!!!!!!    ilrold=ilr
!!!!!!
!!!!!!end do
!!!!!!
!!!!!!
!!!!!!end subroutine gramschmidt_forlinear








subroutine applyOrthoconstraintNonorthogonal2(iproc, nproc, methTransformOverlap, blocksize_pdgemm, orbs, lorbs, onWhichAtom, lzd, &
     op, comon, lagmat, ovrlp, lphiovrlp, mad, lhphi)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => applyOrthoconstraintNonorthogonal2
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm
  type(orbitals_data),intent(in):: orbs, lorbs
  integer,dimension(orbs%norb),intent(in):: onWhichAtom
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  type(p2pCommsOrthonormality),intent(in):: comon
  real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
  real(8),dimension(orbs%norb,orbs%norb),intent(inout):: lagmat
  real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(out):: lhphi

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, info, i, istat, iall, ierr, iseg
  real(8):: tt, t1, t2, time_dsymm, time_daxpy
  real(8),dimension(:,:),allocatable:: ovrlp2
  real(8),dimension(:,:),allocatable:: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
  character(len=*),parameter:: subname='applyOrthoconstraintNonorthogonal2'
  logical:: val, segment

  allocate(ovrlp_minus_one_lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
  allocate(ovrlp_minus_one_lagmat_trans(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
  allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp2, 'ovrlp2', subname)

  call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)

  ! Invert the overlap matrix
  call mpi_barrier(mpi_comm_world, ierr)
  call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, orbs%norb, mad, orbs, ovrlp2)

!!!! test: replace ovrlp2 by unity ##########################
!!!write(*,*) 'ATTENTION: test: replace ovrlp2 by unity'
!!!write(*,*) 'ATTENTION: test: delete off-diagonal elements of lagmat'
!!!do iorb=1,orbs%norb
!!!    do jorb=1,orbs%norb
!!!        if(iorb==jorb) then
!!!            ovrlp2(jorb,iorb)=1.d0
!!!        else
!!!            ovrlp2(jorb,iorb)=0.d0
!!!            lagmat(jorb,iorb)=0.d0
!!!        end if
!!!    end do
!!!end do
!!!! end test ###############################################


  ! Multiply the Lagrange multiplier matrix with S^-1/2.
  ! First fill the upper triangle.
  do iorb=1,orbs%norb
     do jorb=1,iorb-1
        ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
     end do
  end do
  if(blocksize_pdgemm<0) then
     !! ATTENTION: HERE IT IS ASSUMED THAT THE INVERSE OF THE OVERLAP MATRIX HAS THE SAME SPARITY
     !! AS THE OVERLAP MATRIX ITSELF. CHECK THIS!!

     !call dsymm('l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
     !     0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
     t1=mpi_wtime()
     ovrlp_minus_one_lagmat=0.d0
     !!call dgemm_compressed2(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
     !!     mad%keygmatmul, ovrlp2, lagmat, ovrlp_minus_one_lagmat)
     !!do iorb=1,orbs%norb
     !!    do jorb=1,orbs%norb
     !!        if(iproc==0) write(200,*) iorb, jorb, ovrlp_minus_one_lagmat(jorb,iorb)
     !!    end do
     !!end do
    !call dgemm_compressed_parallel(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
    !     mad%keygmatmul, orbs%norb_par, orbs%isorb_par, orbs%norbp, ovrlp2, lagmat, ovrlp_minus_one_lagmat)
    call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
         0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
     !!do iorb=1,orbs%norb
     !!    do jorb=1,orbs%norb
     !!        if(iproc==0) write(201,*) iorb, jorb, ovrlp_minus_one_lagmat(jorb,iorb)
     !!    end do
     !!end do

     if (verbose > 2) then
        t2=mpi_wtime()
        if(iproc==0) write(*,*) 'time for first dgemm_compressed_parallel', t2-t1
        t1=mpi_wtime()
     end if
    ! Transpose lagmat
    do iorb=1,orbs%norb
        do jorb=iorb+1,orbs%norb
            tt=lagmat(jorb,iorb)
            lagmat(jorb,iorb)=lagmat(iorb,jorb)
            lagmat(iorb,jorb)=tt
        end do
    end do
    if (verbose >2) then
       t2=mpi_wtime()
       if(iproc==0) write(*,*) 'time for transposing', t2-t1
    !call dsymm('l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
    !     0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
       t1=mpi_wtime()
    end if
    ovrlp_minus_one_lagmat_trans=0.d0
    !!call dgemm_compressed2(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
    !!     mad%keygmatmul, ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
    !call dgemm_compressed_parallel(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
    !     mad%keygmatmul, orbs%norb_par, orbs%isorb_par, orbs%norbp, ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
    call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
         0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
    if (verbose >2) then
       t2=mpi_wtime()
       if(iproc==0) write(*,*) 'time for first dgemm_compressed_parallel', t2-t1
    end if
else
    call dsymm_parallel(iproc, nproc, blocksize_pdgemm, mpi_comm_world, 'l', 'l', orbs%norb, orbs%norb, 1.d0, &
         ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, 0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
    ! Transpose lagmat
    do iorb=1,orbs%norb
        do jorb=iorb+1,orbs%norb
            tt=lagmat(jorb,iorb)
            lagmat(jorb,iorb)=lagmat(iorb,jorb)
            lagmat(iorb,jorb)=tt
        end do
    end do
    call dsymm_parallel(iproc, nproc, blocksize_pdgemm, mpi_comm_world, 'l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), &
         orbs%norb, lagmat(1,1), orbs%norb, &
         0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
end if
call cpu_time(t2)
time_dsymm=t2-t1




do iorb=1,orbs%norb
    do jorb=1,orbs%norb
        ovrlp2(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat(jorb,iorb)
    end do
end do
call build_new_linear_combinations(lzd, orbs, op, comon%nrecvbuf, comon%recvbuf, ovrlp2, .false., lhphi)

do iorb=1,orbs%norb
    do jorb=1,orbs%norb
        ovrlp2(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat_trans(jorb,iorb)
    end do
end do
call build_new_linear_combinations(lzd, orbs, op, comon%nrecvbuf, comon%recvbuf, ovrlp2, .false., lhphi)


!!call cpu_time(t1)
!!ist=1
!!jst=1
!!ilrold=-1
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=onWhichAtom(iiorb)
!!    if(ilr==ilrold) then
!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!        jst=jst-op%noverlaps(iiorb)*ncount
!!    end if
!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    !call dscal(ncount, 1.5d0, lhphi(ist), 1)
!!    do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
!!        jst=jst+ncount
!!    end do
!!    ist=ist+ncount
!!    ilrold=ilr
!!end do
!!call cpu_time(t2)
!!time_daxpy=t2-t1

if (verbose > 2) then
   call mpiallred(time_dsymm, 1, mpi_sum, mpi_comm_world, ierr)
   call mpiallred(time_daxpy, 1, mpi_sum, mpi_comm_world, ierr)
   if(iproc==0) write(*,'(a,es15.6)') 'time for dsymm',time_dsymm/dble(nproc)
   if(iproc==0) write(*,'(a,es15.6)') 'time for daxpy',time_daxpy/dble(nproc)
end if


iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
deallocate(ovrlp_minus_one_lagmat, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
iall=-product(shape(ovrlp2))*kind(ovrlp2)
deallocate(ovrlp2, stat=istat)
call memocc(istat, iall, 'ovrlp2', subname)


end subroutine applyOrthoconstraintNonorthogonal2


subroutine overlapPowerMinusOne(iproc, nproc, iorder, norb, mad, orbs, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOne
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, iorder, norb
  type(orbitals_data),intent(in):: orbs
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(norb,norb),intent(inout):: ovrlp
  
  ! Local variables
  integer:: lwork, istat, iall, iorb, jorb, info
  character(len=*),parameter:: subname='overlapPowerMinusOne'
  real(8),dimension(:,:),allocatable:: ovrlp2, ovrlp3


  if(iorder==0) then

      ! Exact inversion
      call dpotrf('l', norb, ovrlp(1,1), norb, info)
      if(info/=0) then
          write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
          stop
      end if
      call dpotri('l', norb, ovrlp(1,1), norb, info)
      if(info/=0) then
          write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
          stop
      end if
  
  else if(iorder==1) then
      ! Taylor expansion up to first order.
      do iorb=1,norb
          do jorb=1,norb
              if(iorb==jorb) then
                  ovrlp(jorb,iorb) = 2.d0 - ovrlp(jorb,iorb)
              else
                  ovrlp(jorb,iorb) = -ovrlp(jorb,iorb)
              end if
          end do
      end do
  else if(iorder==2) then
      ! Taylor expansion up to second order.
  
      ! Calculate ovrlp**2
      allocate(ovrlp2(norb,norb), stat=istat)
      call memocc(istat, ovrlp2, 'ovrlp2', subname)
      !!call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, &
      !!     mad%nsegmatmul, mad%keygmatmul, ovrlp, ovrlp, ovrlp2)
      call dgemm_compressed_parallel(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, &
           mad%nsegmatmul, mad%keygmatmul, orbs%norb_par, orbs%isorb_par, orbs%norbp, ovrlp, ovrlp, ovrlp2)
      
      ! Build ovrlp**(-1) with a Taylor expansion up to second order.  
      do iorb=1,norb
          do jorb=1,norb
              if(iorb==jorb) then
                  ovrlp(jorb,iorb) = 3.d0 - 3.d0*ovrlp(jorb,iorb) + ovrlp2(jorb,iorb)
              else
                  ovrlp(jorb,iorb) = - 3.d0*ovrlp(jorb,iorb) + ovrlp2(jorb,iorb)
              end if
          end do
      end do
      
      iall=-product(shape(ovrlp2))*kind(ovrlp2)
      deallocate(ovrlp2, stat=istat)
      call memocc(istat, iall, 'ovrlp2', subname)
  !!else if(iorder==3) then
  !!    ! Taylor expansion up to third iorder.
  !!
  !!    ! Calculate ovrlp**2
  !!    allocate(ovrlp2(norb,norb), stat=istat)
  !!    call memocc(istat, ovrlp2, 'ovrlp2', subname)
  !!    allocate(ovrlp3(norb,norb), stat=istat)
  !!    call memocc(istat, ovrlp3, 'ovrlp3', subname)
  !!    call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, &
  !!         mad%nsegmatmul, mad%keygmatmul, ovrlp, ovrlp, ovrlp2)
  !!    call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, &
  !!         mad%nsegmatmul, mad%keygmatmul, ovrlp, ovrlp2, ovrlp3)
  !!    
  !!    ! Build ovrlp**(-1/2) with a Taylor expansion up to third order.  
  !!    do iorb=1,norb
  !!        do jorb=1,norb
  !!            if(iorb==jorb) then
  !!                ovrlp(jorb,iorb) = 4.d0 - 6.d0*ovrlp(jorb,iorb) + 4.d0*ovrlp2(jorb,iorb) - ovrlp3(jorb,iorb)
  !!            else
  !!                ovrlp(jorb,iorb) = - 6.d0*ovrlp(jorb,iorb) + 4.d0*ovrlp2(jorb,iorb) - ovrlp3(jorb,iorb)
  !!            end if
  !!        end do
  !!    end do
  !!    
  !!    iall=-product(shape(ovrlp2))*kind(ovrlp2)
  !!    deallocate(ovrlp2, stat=istat)
  !!    call memocc(istat, iall, 'ovrlp2', subname)
  !!    iall=-product(shape(ovrlp3))*kind(ovrlp3)
  !!    deallocate(ovrlp3, stat=istat)
  !!    call memocc(istat, iall, 'ovrlp3', subname)
  else
      write(*,'(1x,a)') 'ERROR: iorder must be 0,1 or 2!'
      stop
end if


end subroutine overlapPowerMinusOne




subroutine overlapPowerMinusOneHalf(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb, mad, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOneHalf
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(norb,norb),intent(inout):: ovrlp
  
  ! Local variables
  integer:: lwork, istat, iall, iorb, jorb, info
  character(len=*),parameter:: subname='overlapPowerMinusOneHalf'
  real(8),dimension(:),allocatable:: eval, work
  real(8),dimension(:,:),allocatable:: ovrlp2, ovrlp3
  real(8),dimension(:,:,:),allocatable:: tempArr
  
  if(methTransformOrder==0) then

      ! Exact calculation of ovrlp**(-1/2)

      allocate(eval(norb), stat=istat)
      call memocc(istat, eval, 'eval', subname)
      allocate(tempArr(norb,norb,2), stat=istat)
      call memocc(istat, tempArr, 'tempArr', subname)
      
      
      if(blocksize_dsyev>0) then
          call dsyev_parallel(iproc, nproc, min(blocksize_dsyev,norb), comm, 'v', 'l', norb, ovrlp(1,1), norb, eval(1), info)
          if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
              !stop
          end if
      else
          lwork=1000*norb
          allocate(work(lwork), stat=istat)
          call memocc(istat, work, 'work', subname)
          call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, lwork, info)
          if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev, info=', info
              stop
          end if
          iall=-product(shape(work))*kind(work)
          deallocate(work, stat=istat)
          call memocc(istat, iall, 'work', subname)
      end if
      !!do iorb=1,norb
      !!    do jorb=1,norb
      !!        if(iproc==0) write(1402,'(2i6,es26.17)') iorb, jorb, ovrlp(iorb,jorb)
      !!    end do
      !!    if(iproc==0) write(1450,'(i6,es25.12)') iorb, eval(iorb)
      !!end do
      
      ! Calculate S^{-1/2}. 
      ! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
      ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
      do iorb=1,norb
          do jorb=1,norb
              !tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
              tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
          end do
      end do
      !do iorb=1,norb
      !    do jorb=1,norb
      !        if(iproc==0) write(1403,'(2i6,es26.17)') iorb, jorb, temparr(iorb,jorb,1)
      !    end do
      !end do
      
      ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
      ! This will give S^{-1/2}.
      if(blocksize_pdgemm<0) then
          call dgemm('n', 't', norb, norb, norb, 1.d0, ovrlp(1,1), &
               norb, tempArr(1,1,1), norb, 0.d0, tempArr(1,1,2), norb)
      else
          call dgemm_parallel(iproc, nproc, blocksize_pdgemm, comm, 'n', 't', norb, norb, norb, 1.d0, ovrlp(1,1), &
               norb, tempArr(1,1,1), norb, 0.d0, tempArr(1,1,2), norb)
      end if
      call dcopy(norb**2, tempArr(1,1,2), 1, ovrlp(1,1), 1)
      !do iorb=1,norb
      !    do jorb=1,norb
      !        if(iproc==0) write(1405,'(2i6,es26.17)') iorb, jorb, ovrlp(iorb,jorb)
      !    end do
      !end do
      
      
      iall=-product(shape(eval))*kind(eval)
      deallocate(eval, stat=istat)
      call memocc(istat, iall, 'eval', subname)
      iall=-product(shape(tempArr))*kind(tempArr)
      deallocate(tempArr, stat=istat)
      call memocc(istat, iall, 'tempArr', subname)

  else if(methTransformOrder==1) then

      ! Taylor expansion up to first order.
      do iorb=1,norb
          do jorb=1,norb
              if(iorb==jorb) then
                  ovrlp(jorb,iorb)=1.5d0-.5d0*ovrlp(jorb,iorb)
              else
                  ovrlp(jorb,iorb)=-.5d0*ovrlp(jorb,iorb)
              end if
          end do
      end do

  else if(methTransformOrder==2) then

      ! Taylor expansion up to second order.
  
      ! Calculate ovrlp**2
      allocate(ovrlp2(norb,norb), stat=istat)
      call memocc(istat, ovrlp2, 'ovrlp2', subname)
      call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, &
           mad%nsegmatmul, mad%keygmatmul, ovrlp, ovrlp, ovrlp2)
      
      ! Build ovrlp**(-1/2) with a Taylor expansion up to second order.  
      do iorb=1,norb
          do jorb=1,norb
              if(iorb==jorb) then
                  ovrlp(jorb,iorb) = 1.125d0 + .25d0*ovrlp(jorb,iorb) - .375d0*ovrlp2(jorb,iorb)
              else
                  ovrlp(jorb,iorb) = .25d0*ovrlp(jorb,iorb) - .375d0*ovrlp2(jorb,iorb)
              end if
          end do
      end do
      
      iall=-product(shape(ovrlp2))*kind(ovrlp2)
      deallocate(ovrlp2, stat=istat)
      call memocc(istat, iall, 'ovrlp2', subname)

  !!else if(methTransformOrder==3) then

  !!    ! Taylor expansion up to third order.
  !!
  !!    ! Calculate ovrlp**2
  !!    allocate(ovrlp2(norb,norb), stat=istat)
  !!    call memocc(istat, ovrlp2, 'ovrlp2', subname)
  !!    allocate(ovrlp3(norb,norb), stat=istat)
  !!    call memocc(istat, ovrlp3, 'ovrlp3', subname)
  !!    call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, &
  !!         mad%nsegmatmul, mad%keygmatmul, ovrlp, ovrlp, ovrlp2)
  !!    call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, &
  !!         mad%nsegmatmul, mad%keygmatmul, ovrlp, ovrlp2, ovrlp3)
  !!    
  !!    ! Build ovrlp**(-1/2) with a Taylor expansion up to third order.  
  !!    do iorb=1,norb
  !!        do jorb=1,norb
  !!            if(iorb==jorb) then
  !!                ovrlp(jorb,iorb) = 1.4375d0 - 0.6875d0*ovrlp(jorb,iorb) + 0.5625d0*ovrlp2(jorb,iorb) - 0.3125d0*ovrlp3(jorb,iorb)
  !!            else
  !!                ovrlp(jorb,iorb) = - 0.6875d0*ovrlp(jorb,iorb) + 0.5625d0*ovrlp2(jorb,iorb) - 0.3125d0*ovrlp3(jorb,iorb)
  !!            end if
  !!        end do
  !!    end do
  !!    
  !!    iall=-product(shape(ovrlp2))*kind(ovrlp2)
  !!    deallocate(ovrlp2, stat=istat)
  !!    call memocc(istat, iall, 'ovrlp2', subname)
  !!    iall=-product(shape(ovrlp3))*kind(ovrlp3)
  !!    deallocate(ovrlp3, stat=istat)
  !!    call memocc(istat, iall, 'ovrlp3', subname)

  else

      write(*,'(1x,a)') 'ERROR: methTransformOrder must be 0,1 or 2!'
      stop

end if


endsubroutine overlapPowerMinusOneHalf





subroutine indicesForExpansion(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon)
use module_base
use module_types
use module_interfaces, exceptThisOne => indicesForExpansion
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(in):: comon

! Local variables
integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, istat
character(len=*),parameter:: subname='indicesForExpansion'


ind=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) cycle
    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
        call index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, &
             lzd%llr(ilr), op%olr(jorb,iorb), op%indexExpand(jst))

        call countExpansionSegments(ldim, op%indexExpand(jst), op%expseg(jorb,iorb)%nseg)
        allocate(op%expseg(jorb,iorb)%segborders(2,op%expseg(jorb,iorb)%nseg), stat=istat)
        call memocc(istat, op%expseg(jorb,iorb)%segborders, 'op%expseg(jorb,iorb)%segborders', subname)
        call determineExpansionSegments(ldim, op%indexExpand(jst), op%expseg(jorb,iorb)%nseg, op%expseg(jorb,iorb)%segborders)

        ind=ind+gdim
    end do
    ilrold=ilr
end do

end subroutine indicesForExpansion



subroutine countExpansionSegments(ldim, indexExpand, nseg)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: ldim
integer,dimension(ldim),intent(in):: indexExpand
integer,intent(out):: nseg
integer,dimension(:,:),pointer:: segborders

! Local variables
integer:: i

! Count the numbers of segments
nseg=0
do i=2,ldim
    if(indexExpand(i)==indexExpand(i-1)+1) then
        !consecutive, same segment
    else
        !segment ended
        nseg=nseg+1
    end if
end do
nseg=nseg+1 !last segment

end subroutine countExpansionSegments


subroutine determineExpansionSegments(ldim, indexExpand, nseg, segborders)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: ldim, nseg
integer,dimension(ldim),intent(in):: indexExpand
integer,dimension(2,nseg):: segborders

! Local variables
integer:: iseg, i
character(len=*),parameter:: subname='determineExpansionSegments'

iseg=1
segborders(1,iseg)=indexExpand(1)
do i=2,ldim
    if(indexExpand(i)==indexExpand(i-1)+1) then
        !consecutive, same segment
    else
        !segment ended
        segborders(2,iseg)=indexExpand(i-1)
        iseg=iseg+1
        segborders(1,iseg)=indexExpand(i)
    end if
end do
segborders(2,iseg)=indexExpand(ldim) !last segments



end subroutine determineExpansionSegments




subroutine indicesForExtraction(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, comon)
use module_base
use module_types
use module_interfaces, exceptThisOne => indicesForExtraction
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, sizePhi
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
integer:: i, istat
character(len=*),parameter:: subname='indicesForExtraction'

indovrlp=1
op%indexInSendBuf=0

ilrold=-1
do iorb=1,orbs%norb
    ilr=onWhichAtom(iorb)
    iiproc=orbs%onWhichMPI(iorb)
    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
    do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=onWhichAtom(jjorb)
        jjproc=orbs%onWhichMPI(jjorb)
        if(iproc==jjproc) then
            ! Get the correct descriptors
            korb=jjorb-orbs%isorb
            !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
            do i=1,op%noverlaps(jjorb)
                !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
                if(op%overlaps(i,jjorb)==iorb) then
                    lorb=i
                    exit
                end if
            end do
            !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
            ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
            ind=1
            do kkorb=orbs%isorb+1,jjorb-1
                klr=onWhichAtom(kkorb)
                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
            end do
            !write(*,'(5(a,i0))') 'process ',iproc,' adds ',op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f,' elements at position ',indovrlp,' from orbital ',jjorb,' for orbital ', iorb
            !call psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), phi(ind), comon%sendBuf(indovrlp))
            call index_of_psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), op%indexExtract(indovrlp))

            call countExpansionSegments(ldim, op%indexExtract(indovrlp), op%extseg(lorb,korb)%nseg)
            allocate(op%extseg(lorb,korb)%segborders(2,op%extseg(lorb,korb)%nseg), stat=istat)
            call memocc(istat, op%extseg(lorb,korb)%segborders, 'op%extseg(lorb,korb)%segborders', subname)
            call determineExpansionSegments(ldim, op%indexExtract(indovrlp), op%extseg(lorb,korb)%nseg, &
                 op%extseg(lorb,korb)%segborders)

            op%indexInSendBuf(jjorb-orbs%isorb,iorb)=indovrlp
            indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(indovrlp/=comon%nsendBuf+1) then
    write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
    stop
end if



end subroutine indicesForExtraction






subroutine getMatrixElements2(iproc, nproc, lzd, orbs, op_lb, comon_lb, lphi, lhphi, mad, matrixElements)
use module_base
use module_types
use module_interfaces, exceptThisOne => getMatrixElements2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(inout):: op_lb
type(p2pCommsOrthonormality),intent(inout):: comon_lb
real(8),dimension(orbs%npsidim_orbs),intent(in):: lphi, lhphi
type(matrixDescriptors),intent(in):: mad
real(8),dimension(orbs%norb,orbs%norb),intent(out):: matrixElements

! Local variables
integer:: it, istat, iall, iorb
real(8),dimension(:),allocatable:: lphiovrlp
character(len=*),parameter:: subname='getMatrixElements2'
real(8),dimension(:,:),allocatable:: ttmat
real(8):: tt1, tt2, tt3
type(input_variables):: input



  !!call allocateCommuncationBuffersOrtho(comon_lb, subname)
  ! Put lphi in the sendbuffer,i.e. lphi will be sent to other processes' receive buffer.
  !call extractOrbital2(iproc,nproc,orbs,orbs%npsidim,orbs%inWhichLocreg,lzd,op_lb,lphi,comon_lb)
  call extractOrbital3(iproc,nproc,orbs,orbs%npsidim_orbs,orbs%inWhichLocreg,lzd,op_lb,lphi,comon_lb%nsendBuf,comon_lb%sendBuf)
  !call postCommsOverlap(iproc,nproc,comon_lb)
  call postCommsOverlapNew(iproc,nproc,orbs,op_lb,lzd,lphi,comon_lb,tt1,tt2)
  !call gatherOrbitals2(iproc,nproc,comon_lb)
  allocate(ttmat(orbs%norb,orbs%norb))
  call collectnew(iproc,nproc,comon_lb,mad,op_lb,orbs,input,lzd,comon_lb%nsendbuf,&
       comon_lb%sendbuf,comon_lb%nrecvbuf,comon_lb%recvbuf,ttmat,tt1,tt2,tt3)
  deallocate(ttmat)
  ! Put lhphi to the sendbuffer,so we can the calculate <lphi|lhphi>
  !call extractOrbital2(iproc,nproc,orbs,orbs%npsidim,orbs%inWhichLocreg,lzd,op_lb,lhphi,comon_lb)
  call extractOrbital3(iproc,nproc,orbs,orbs%npsidim_orbs,orbs%inWhichLocreg,lzd,op_lb,lhphi,comon_lb%nsendBuf,comon_lb%sendBuf)
  !!call calculateOverlapMatrix2(iproc,nproc,orbs,op_lb,comon_lb,orbs%inWhichLocreg,mad,matrixElements)
  call calculateOverlapMatrix3(iproc,nproc,orbs,op_lb,orbs%inWhichLocreg,comon_lb%nsendBuf,&
                               comon_lb%sendBuf,comon_lb%nrecvBuf,comon_lb%recvBuf,mad,matrixElements)
  !!if(iproc==0) then
  !!    do iall=1,orbs%norb
  !!        do istat=1,orbs%norb
  !!            write(27000+iproc,*) iall, istat, matrixElements(istat,iall)
  !!        end do
  !!    end do
  !!    write(27000+iproc,*) '=============================='
  !!end if
  !!call deallocateCommuncationBuffersOrtho(comon_lb, subname)

  !!if(iproc==0) then
  !!    do iall=1,orbs%norb
  !!        do istat=1,orbs%norb
  !!            write(26000+iproc,*) iall, istat, matrixElements(istat,iall)
  !!        end do
  !!    end do
  !!    write(26000+iproc,*) '=============================='
  !!end if

end subroutine getMatrixElements2




subroutine initCommsOrthoVariable(iproc, nproc, lzd, orbs, orbsig, onWhichAtomAll, input, op, comon, tag)
use module_base
use module_types
use module_interfaces, exceptThisOne => initCommsOrthoVariable
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs, orbsig
integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
type(input_variables),intent(in):: input
type(overlapParameters),intent(out):: op
type(p2pCommsOrthonormality),intent(out):: comon
integer,intent(inout):: tag

! Local variables
integer:: iorb, jorb, iiorb, jproc, ioverlaporb, ioverlapMPI, ilr, jlr
integer:: ilrold, is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3
integer::  je3, istat, i1, i2
logical:: ovrlpx, ovrlpy, ovrlpz
character(len=*),parameter:: subname='initCommsOrthoVariable'

allocate(comon%noverlaps(0:nproc-1), stat=istat)
call memocc(istat, comon%noverlaps, 'comon%noverlaps',subname)
allocate(op%noverlaps(orbs%norb), stat=istat)
call memocc(istat, op%noverlaps, 'op%noverlaps',subname)

! Count how many overlaping regions each orbital / process has.
call countOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)

op%noverlapsmax=maxval(op%noverlaps)
allocate(op%overlaps(op%noverlapsmax,orbs%norb), stat=istat)
call memocc(istat, op%overlaps, 'op%overlaps', subname)
comon%noverlapsmax=maxval(comon%noverlaps)
allocate(comon%overlaps(comon%noverlapsmax,0:nproc-1), stat=istat)
call memocc(istat, comon%overlaps, 'comon%overlaps', subname)
allocate(op%indexInRecvBuf(orbs%norbp,orbsig%norb), stat=istat)
call memocc(istat, op%indexInRecvBuf, 'op%indexInRecvBuf', subname)
allocate(op%indexInSendBuf(orbsig%norbp,orbsig%norb), stat=istat)
call memocc(istat, op%indexInSendBuf, 'op%indexInSendBuf', subname)

! Determine the overlapping orbitals.
call determineOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)

allocate(op%olr(op%noverlapsmax,orbs%norb), stat=istat)
!do i2=1,orbs%norbp
do i2=1,orbs%norb
    do i1=1,op%noverlapsmax
        call nullify_locreg_descriptors(op%olr(i1,i2))
    end do
end do

! Set the orbital descriptors for the overlap regions.
call determineOverlapDescriptorsVariable(iproc, nproc, orbs, orbsig, lzd, lzd%Glr, onWhichAtomAll, op)


allocate(comon%comarr(10,comon%noverlapsmax,0:nproc-1), stat=istat)
call memocc(istat, comon%comarr, 'comon%comarr', subname)
allocate(comon%communComplete(comon%noverlapsmax,0:nproc-1), stat=istat)
call memocc(istat, comon%communComplete, 'comun%communComplete', subname)
call setCommsOrthoVariable(iproc, nproc, orbs, orbsig, lzd, op, comon, tag)



! Initialize the index arrays for the transformations from overlap region
! to ordinary localization region.
allocate(op%indexExpand(comon%nrecvBuf), stat=istat)
call memocc(istat, op%indexExpand, 'op%indexExpand',subname)
call indicesForExpansionVariable(iproc, nproc, orbs, input, lzd, op, comon)

! Initialize the index arrays for the transformations from the ordinary localization region
! to the overlap region.
allocate(op%indexExtract(comon%nsendBuf), stat=istat)
call memocc(istat, op%indexExtract, 'op%indexExtract',subname)
call indicesForExtractionVariable(iproc, nproc, orbs, orbsig, orbs%npsidim_orbs, lzd, op, comon)

end subroutine initCommsOrthoVariable



! Count for each orbital and each process the number of overlapping orbitals.
subroutine countOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, orbsig
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(out):: op
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
integer:: js1, je1, js2, je2, js3, je3, iiorb
logical:: ovrlpx, ovrlpy, ovrlpz

iiorb=0
do jproc=0,nproc-1
    ioverlapMPI=0 ! counts the overlaps for the given MPI process.
    ilrold=-1
    do iorb=1,orbs%norb_par(jproc,0)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=orbs%inWhichLocreg(iiorb)
        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbsig%norb
            jlr=orbsig%inWhichLocreg(jorb)
            call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
            ovrlpx = ( is1<=je1 .and. ie1>=js1 )
            ovrlpy = ( is2<=je2 .and. ie2>=js2 )
            ovrlpz = ( is3<=je3 .and. ie3>=js3 )
            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                ioverlaporb=ioverlaporb+1
                if(ilr/=ilrold) then
                    ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
                    ! would get the same orbitals again. Therefore the counter is not increased
                    ! in that case.
                    ioverlapMPI=ioverlapMPI+1
                end if
            end if
        end do 
        op%noverlaps(iiorb)=ioverlaporb
        !!if(iproc==0) write(*,'(a,2i8)') 'iiorb, op%noverlaps(iiorb)', iiorb, op%noverlaps(iiorb)
        ilrold=ilr
    end do
    comon%noverlaps(jproc)=ioverlapMpi
    !if(iproc==0) write(*,'(a,2i8)') 'jproc, comon%noverlaps(jproc)', jproc, comon%noverlaps(jproc)
end do

end subroutine countOverlapsVariable




subroutine determineOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, orbsig
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(out):: op
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
integer:: js1, je1, js2, je2, js3, je3, iiorb
logical:: ovrlpx, ovrlpy, ovrlpz

  ! Initialize to some value which will never be used.
  op%overlaps=-1
  comon%overlaps=-1

  iiorb=0
  do jproc=0,nproc-1
      ioverlapMPI=0 ! counts the overlaps for the given MPI process.
      ilrold=-1
      do iorb=1,orbs%norb_par(jproc,0)
          ioverlaporb=0 ! counts the overlaps for the given orbital.
          iiorb=iiorb+1 ! counts the total orbitals
          ilr=orbs%inWhichLocreg(iiorb)
          call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
          do jorb=1,orbsig%norb
              jlr=orbsig%inWhichLocreg(jorb)
              call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
              ovrlpx = ( is1<=je1 .and. ie1>=js1 )
              ovrlpy = ( is2<=je2 .and. ie2>=js2 )
              ovrlpz = ( is3<=je3 .and. ie3>=js3 )
              if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                  ioverlaporb=ioverlaporb+1
                  op%overlaps(ioverlaporb,iiorb)=jorb
                  if(ilr/=ilrold) then
                      ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
                      ! would get the same orbitals again. Therefore the counter is not increased
                      ! in that case.
                      ioverlapMPI=ioverlapMPI+1
                      comon%overlaps(ioverlapMPI,jproc)=jorb
                  end if
              end if
          end do 
          !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
          ilrold=ilr
      end do
      !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
  end do

end subroutine determineOverlapsVariable





subroutine determineOverlapDescriptorsVariable(iproc, nproc, orbs, orbsig, lzd, Glr, onWhichAtom, op)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, orbsig
type(local_zone_descriptors),intent(in):: lzd
type(locreg_descriptors),intent(in):: Glr
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(overlapParameters),intent(inout):: op

! Local variables
integer:: iorb, jorb, jjorb, ilr, jlr, iiorb

do iiorb=1,orbs%norb
    ilr=orbs%inWhichLocreg(iiorb)
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        jlr=orbsig%inWhichLocreg(jjorb)
        call get_overlap_region_periodic2(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iiorb))
    end do
end do

end subroutine determineOverlapDescriptorsVariable



subroutine setCommsOrthoVariable(iproc, nproc, orbs, orbsig, lzd, op, comon, tag)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, orbsig
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(out):: comon
integer,intent(inout):: tag

! Local variables
integer:: jproc, iorb, jorb, iiorb, jjorb, mpisource, mpidest, istsource, istdest, ncount, istat, iall, ijorb
integer:: ilr, ilrold, jprocold, ildim, ierr
integer,dimension(:),allocatable:: istsourceArr, istdestArr
character(len=*),parameter:: subname='setCommsOrtho'
logical,dimension(:),allocatable:: receivedOrbital

allocate(istsourceArr(0:nproc-1), stat=istat)
call memocc(istat, istsourceArr, 'istsourceArr',subname)
allocate(istdestArr(0:nproc-1), stat=istat)
call memocc(istat, istdestArr, 'istdestArr',subname)
allocate(receivedOrbital(orbsig%norb), stat=istat)
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
       ilr=orbs%inWhichLocreg(iiorb)
       ! Check whether process jproc has already received orbital jjorb.
       !if(iproc==0) write(*,'(a,5i8)') 'jproc, iorb, iiorb, ilr, ilrold', jproc, iorb, iiorb, ilr, ilrold
       if(ilr==ilrold) cycle
       ildim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
       do jorb=1,op%noverlaps(iiorb)
           jjorb=op%overlaps(jorb,iiorb)
           !write(*,'(a,7i8)') 'iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold', iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold
           ijorb=ijorb+1
           mpisource=orbsig%onWhichMPI(jjorb)
           mpidest=jproc
           istsource=istsourceArr(mpisource)
           istdest=istdestArr(mpidest)
           if(iproc==jproc) then
               ncount=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
           end if
           call mpi_bcast(ncount, 1, mpi_integer, jproc, mpi_comm_world, ierr)
           tag=tag+1
           receivedOrbital(jjorb)=.true.
           call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comon%comarr(1,ijorb,jproc))
           !if(iproc==0) write(*,'(6(a,i0))') 'process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,', tag=',tag
           if(iproc==mpisource) then
               !write(*,'(5(a,i0))') 'adding ',ncount,' elements from orbital ',jjorb,' for orbital ',iiorb,' to nsendBuf, iproc=',iproc,', jproc=',jproc
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


end subroutine setCommsOrthoVariable




subroutine indicesForExpansionVariable(iproc, nproc, orbs, input, lzd, op, comon)
use module_base
use module_types
use module_interfaces, exceptThisOne => indicesForExpansionVariable
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(in):: comon

! Local variables
integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, ierr



ind=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    if(ilr==ilrold) cycle
    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        call index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, &
             input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), op%indexExpand(jst:jst+ldim-1))
        ind=ind+gdim
    end do
    ilrold=ilr
end do
if(jst+ldim/=comon%nrecvBuf+1) then
    write(*,*) 'ERROR on process ',iproc,': jst+ldim/=comon%nrecvBuf+1',jst+ldim,comon%nrecvBuf+1
    stop
end if

end subroutine indicesForExpansionVariable



subroutine indicesForExtractionVariable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, sizePhi
type(orbitals_data),intent(in):: orbs, orbsig
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
integer:: i, jj, j, jjjlr, iilr, jlr, jjorb2, iiorb, ijorb

indovrlp=1
op%indexInSendBuf=0

ilrold=-1
ijorb=0
do iorb=1,orbs%norb
    ilr=orbs%inWhichLocreg(iorb)
    iiproc=orbs%onWhichMPI(iorb)
    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
    do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=orbsig%inWhichLocreg(jjorb)
        jjproc=orbsig%onWhichMPI(jjorb)
        if(iproc==jjproc) then
            ijorb=ijorb+1
            ! Get the correct descriptors
            ! Get an orbs-orbital in the same locreg as the orbsig-orbital jjorb.
            do j=1,orbs%norb
                jjjlr=orbs%inWhichLocreg(j)
                if(jjjlr==jjlr) then
                    korb=j
                    exit
                end if
            end do
            ! Get an orbsig-orbital in the same locreg as the orbs-orbital iorb.
            lorb=0
            do i=1,op%noverlaps(korb)
                iiorb=op%overlaps(i,korb)
                iilr=orbsig%inWhichLocreg(iiorb)
                if(iilr==ilr) then
                    lorb=i
                    exit
                end if
            end do
            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
            ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
            ind=1
            do kkorb=orbs%isorb+1,jjorb-1
                klr=orbsig%inWhichLocreg(kkorb)
                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
            end do
            call index_of_psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), op%indexExtract(indovrlp))
            op%indexInSendBuf(jjorb-orbsig%isorb,iorb)=indovrlp
            indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(indovrlp/=comon%nsendBuf+1) then
    write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
    stop
end if



end subroutine indicesForExtractionVariable



subroutine extractOrbital2Variable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, phi, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, sizePhi
type(orbitals_data),intent(in):: orbs, orbsig
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
real(8),dimension(sizePhi),intent(in):: phi
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
integer:: i, indSource, j, jjjlr, iiorb, iilr

indovrlp=1
op%indexInSendBuf=0

ilrold=-1
iiprocold=-1
do iorb=1,orbs%norb
    ilr=orbs%inWhichLocreg(iorb)
    iiproc=orbs%onWhichMPI(iorb)
    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
    do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=orbsig%inWhichLocreg(jjorb)
        jjproc=orbsig%onWhichMPI(jjorb)
        if(iproc==jjproc) then
            ! Get the correct descriptors
            ! Get an orbs-orbital in the same locreg as the orbsig-orbital jjorb.
            do j=1,orbs%norb
                jjjlr=orbs%inWhichLocreg(j)
                if(jjjlr==jjlr) then
                    korb=j
                    exit
                end if
            end do
            ! Get an orbsig-orbital in the same locreg as the orbs-orbital iorb.
            lorb=0
            do i=1,op%noverlaps(korb)
                iiorb=op%overlaps(i,korb)
                iilr=orbsig%inWhichLocreg(iiorb)
                if(iilr==ilr) then
                    lorb=i
                    exit
                end if
            end do
            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
            ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
            ind=1
            do kkorb=orbsig%isorb+1,jjorb-1
                klr=orbsig%inWhichLocreg(kkorb)
                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
            end do
            do i=0,ldim-1
                indSource=ind+op%indexExtract(indovrlp+i)-1
                comon%sendBuf(indovrlp+i)=phi(indSource)
            end do
            op%indexInSendBuf(jjorb-orbsig%isorb,iorb)=indovrlp
            indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(indovrlp/=comon%nsendBuf+1) then
    write(*,'(1x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
    stop
end if

end subroutine extractOrbital2Variable



subroutine expandOrbital2Variable(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(in):: comon
real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp

! Local variables
integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest


lphiovrlp=0.d0

ind=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    if(ilr==ilrold) cycle
    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        do i=0,ldim-1
            indDest=ind+op%indexExpand(jst+i)-1
            lphiovrlp(indDest)=comon%recvBuf(jst+i)
        end do
        ind=ind+gdim
    end do
    ilrold=ilr
end do

end subroutine expandOrbital2Variable





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


subroutine getStartingIndices(iorb, jorb, op, orbs, ist, jst)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iorb, jorb
  type(overlapParameters),intent(in):: op
  type(orbitals_data),intent(in):: orbs
  integer,intent(out):: ist, jst
  
  ! Local variables
  integer:: jjlr, jjproc, jj, iilr, iiproc, korb, kkorb, ii, iiorb, jjorb


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




subroutine getStartingIndicesGlobal(iiorbx, jjorbx, op, orbs, ist, jst, ncount)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => getStartingIndicesGlobal
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iiorbx, jjorbx
  type(overlapParameters),intent(in):: op
  type(orbitals_data),intent(in):: orbs
  integer,intent(out):: ist, jst, ncount
  
  ! Local variables
  integer:: iiorb, jjorb, iorb, jorb

  ! This only works localy on a process, i.e. we can only get the informations related
  ! to the currect process (but this is enough)
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jorb=1,op%noverlaps(iiorb)
          jjorb=op%overlaps(jorb,iiorb)
          if(iiorb==iiorbx .and. jjorb==jjorbx) then
              call getStartingIndices(iorb, jorb, op, orbs, ist, jst)
              ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
          end if
      end do
  end do

end subroutine getStartingIndicesGlobal






subroutine postCommsOverlapNew(iproc, nproc, orbs, op, lzd, phi, comon, timecommun, timeextract)
use module_base
use module_types
use module_interfaces, exceptThisOne => postCommsOverlapNew
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(in):: op
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: phi
type(p2pCommsOrthonormality),intent(inout):: comon
real(8),intent(inout):: timecommun, timeextract

! Local variables
integer:: jproc, jorb, mpisource, istsource, ncount, mpidest, istdest, tag, nsends, ierr, irecv, isend, iall, orbsource, orbdest
integer:: i, indsource, ind, klr, kkorb, istat, iorb
integer,dimension(:),allocatable:: istextracted
logical:: done
real(8):: t1, t2



! First only post receives
irecv=0
do jproc=0,nproc-1
    do jorb=1,comon%noverlaps(jproc)
        mpisource=comon%comarr(1,jorb,jproc)
        istsource=comon%comarr(2,jorb,jproc)
        ncount=comon%comarr(3,jorb,jproc)
        mpidest=comon%comarr(4,jorb,jproc)
        istdest=comon%comarr(5,jorb,jproc)
        tag=comon%comarr(6,jorb,jproc)
        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
        !irecv=irecv+1
        if(iproc==mpidest .and. nproc >1) then
            irecv=irecv+1
            !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
            !call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
            !     mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
            tag=jorb
            t1=mpi_wtime()
            !write(*,'(2(a,i0))') 'iproc=',iproc,' receives data at ',istdest
            call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
                 mpi_comm_world, comon%requests(irecv,2), ierr)
            t2=mpi_wtime()
            timecommun=timecommun+t2-t1
        !else
        !     comon%requests(irecv,2)=mpi_request_null
        end if
    end do
end do
call mpi_barrier(mpi_comm_world, ierr)





! Now the sends
isend=0
do jproc=0,nproc-1
    do jorb=1,comon%noverlaps(jproc)
        mpisource=comon%comarr(1,jorb,jproc)
        istsource=comon%comarr(2,jorb,jproc)
        ncount=comon%comarr(3,jorb,jproc)
        mpidest=comon%comarr(4,jorb,jproc)
        istdest=comon%comarr(5,jorb,jproc)
        tag=comon%comarr(6,jorb,jproc)
        ! The orbitals are on different processes, so we need a point to point communication.
        !isend=isend+1
        if(iproc==mpisource .and. nproc >1) then
            isend=isend+1
            t2=mpi_wtime()
            timeextract=timeextract+t2-t1
            tag=jorb
            t1=mpi_wtime()
            call mpi_irsend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
                 mpi_comm_world, comon%requests(isend,1), ierr)
            !call mpi_rsend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
            !     mpi_comm_world, ierr)
            t2=mpi_wtime()
            timecommun=timecommun+t2-t1
         else if (nproc==1) then
            call vcopy(ncount,comon%sendBuf(istsource),1,comon%recvBuf(istdest),1)
        end if
    end do
end do
!!comon%nsend=isend
!!if(isend/=comon%nsend) then
!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR in process ',iproc,': irecv/=comon%nrecv',irecv,comon%nsend
!!    stop
!!end if







end subroutine postCommsOverlapNew






!subroutine postCommsOverlapNew2(iproc, nproc, orbs, op, lzd, phi, comon, timecommun, timeextract)
!use module_base
!use module_types
!use module_interfaces, exceptThisOne => postCommsOverlapNew
!implicit none
!
!! Calling arguments
!integer,intent(in):: iproc, nproc
!type(orbitals_data),intent(in):: orbs
!type(overlapParameters),intent(in):: op
!type(local_zone_descriptors),intent(in):: lzd
!real(8),dimension(lzd%lpsidimtot),intent(in):: phi
!type(p2pCommsOrthonormality),intent(inout):: comon
!real(8),intent(inout):: timecommun, timeextract
!
!! Local variables
!integer:: jproc, jorb, mpisource, istsource, ncount, mpidest, istdest, tag, nsends, ierr, irecv, isend, iall, orbsource, orbdest
!integer:: i, indsource, ind, klr, kkorb, istat, iorb
!integer,dimension(:),allocatable:: istextracted
!logical:: done
!real(8):: t1, t2
!
!
!
!! First only post receives
!irecv=0
!do jproc=0,nproc-1
!    !do jorb=1,comon%noverlaps(jproc)
!    do jorb=max(comon%isoverlap,1),min(comon%isoverlap+comon%nstepoverlap-1,comon%noverlaps(jproc))
!        mpisource=comon%comarr(1,jorb,jproc)
!        istsource=comon%comarr(2,jorb,jproc)
!        ncount=comon%comarr(3,jorb,jproc)
!        mpidest=comon%comarr(4,jorb,jproc)
!        istdest=comon%comarr(5,jorb,jproc)
!        tag=comon%comarr(6,jorb,jproc)
!        !write(*,'(6(a,i0))') 'iproc=',iproc,', tag=',tag,', mpisource=',mpisource,', mpidest=',mpidest,' jproc=',jproc,', iorb=',iorb
!        !irecv=irecv+1
!        if(iproc==mpidest) then
!            irecv=irecv+1
!            !write(*,'(6(a,i0))') 'overlap: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!            !call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!            !     mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
!            tag=jorb
!            t1=mpi_wtime()
!            !write(*,'(2(a,i0))') 'iproc=',iproc,' receives data at ',istdest
!            call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
!                 mpi_comm_world, comon%requests(irecv,2), ierr)
!            t2=mpi_wtime()
!            timecommun=timecommun+t2-t1
!        !else
!        !     comon%requests(irecv,2)=mpi_request_null
!        end if
!    end do
!end do
!comon%nrecvtemp=irecv
!
!
!
!
!
!! Now the sends
!isend=0
!do jproc=0,nproc-1
!    !do jorb=1,comon%noverlaps(jproc)
!    do jorb=max(comon%isoverlap,1),min(comon%isoverlap+comon%nstepoverlap-1,comon%noverlaps(jproc))
!        mpisource=comon%comarr(1,jorb,jproc)
!        istsource=comon%comarr(2,jorb,jproc)
!        ncount=comon%comarr(3,jorb,jproc)
!        mpidest=comon%comarr(4,jorb,jproc)
!        istdest=comon%comarr(5,jorb,jproc)
!        tag=comon%comarr(6,jorb,jproc)
!        ! The orbitals are on different processes, so we need a point to point communication.
!        !isend=isend+1
!        if(iproc==mpisource) then
!            isend=isend+1
!            tag=jorb
!            t1=mpi_wtime()
!            !call mpi_irsend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!            !     mpi_comm_world, comon%requests(isend,1), ierr)
!            call mpi_isend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!                 mpi_comm_world, comon%requests(isend,1), ierr)
!            !call mpi_rsend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
!            !     mpi_comm_world, ierr)
!            t2=mpi_wtime()
!            timecommun=timecommun+t2-t1
!        end if
!    end do
!end do
!comon%nsendtemp=isend
!!!comon%nsend=isend
!!!if(isend/=comon%nsend) then
!!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR in process ',iproc,': irecv/=comon%nrecv',irecv,comon%nsend
!!!    stop
!!!end if
!
!
!
!
!
!
!
!end subroutine postCommsOverlapNew2
!

subroutine collectAndCalculateOverlap(iproc, nproc, comon, mad, op, orbs, input, lzd, &
           nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, lphiovrlp,timecommunp2p, &
           timecommuncoll, timeoverlap, timeexpand, timecompress)
use module_base
use module_types
use module_interfaces, exceptThisOne => collectAndCalculateOverlap
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
type(p2pCommsOrthonormality),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
type(overlapParameters),intent(in):: op
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(nsendbuf),intent(in):: sendbuf
real(8),dimension(nrecvbuf),intent(inout):: recvbuf
real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
real(8),intent(inout):: timecommunp2p, timecommuncoll, timeoverlap, timeexpand, timecompress

! Local variables
integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
real(8):: ddot, t1, t2
character(len=*),parameter:: subname='collectAndCalculateOverlap'
logical,dimension(:),allocatable:: done
logical:: received

allocate(indcomplete(comon%nrecv), stat=istat)
call memocc(istat, indcomplete, 'indcomplete', subname)
allocate(done(comon%nrecv), stat=istat)
call memocc(istat, done, 'done', subname)
done=.false.

allocate(indexarray(comon%nrecv), stat=istat)
call memocc(istat, indexarray, 'indexarray', subname)
do i=1,comon%nrecv
    indexarray(i)=i
end do


! Now the sends
if (nproc >1) then
   t1=mpi_wtime()
   nsend=0
   waitLoopSend: do
      !!call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
      !!nsend=nsend+ncomplete
      !!if(nsend==comon%nsend) exit waitLoopSend
      call mpi_waitany(comon%nsend-nsend, comon%requests(1,1), ind, mpi_status_ignore, ierr)
      nsend=nsend+1
      do i=ind,comon%nsend-nsend
         comon%requests(i,1)=comon%requests(i+1,1)
      end do
      if(nsend==comon%nsend) exit waitLoopSend
   end do waitLoopSend
   t2=mpi_wtime()
   timecommunp2p=timecommunp2p+t2-t1
end if

ovrlp=0.d0
lphiovrlp=0.d0
nrecv=0
waitLoopRecv: do
   if (nproc > 1) then
      t1=mpi_wtime()
      call mpi_waitany(comon%nrecv-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
      !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
      !ind=1
      t2=mpi_wtime()
      timecommunp2p=timecommunp2p+t2-t1
   end if
   ncomplete=1
   received=.true.
   if(received) then
      do i=1,ncomplete
         ! Calculate overlap matrix
         !jorb=indcomplete(i)
         jorb=indexarray(ind)
         !!write(*,'(2(a,i0))') 'process ',iproc,' has received message ',jorb
         if(.not.done(jorb)) then
            orbsource=comon%comarr(9,jorb,iproc)
            do iorb=1,orbs%norbp
               orbdest=orbs%isorb+iorb
               call getStartingIndicesGlobal(orbdest, orbsource, op, orbs, ist, jst, ncount)
               !jst=comon%comarr(5,jorb,iproc)
               !ncount=comon%comarr(3,jorb,iproc)
               !write(*,'(a,i4,2i5,3x,2i8,4x,2i8)') 'iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc) ; jst, comon%comarr(5,jorb,iproc)', iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc), jst, comon%comarr(5,jorb,iproc)
               t1=mpi_wtime()
               ovrlp(orbsource,orbdest)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
               t2=mpi_wtime()
               timeoverlap=timeoverlap+t2-t1
               t1=mpi_wtime()
               call expandoneorbital2(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
                    orbs%inwhichlocreg, lzd, op, nrecvbuf, recvbuf, lphiovrlp)
               t2=mpi_wtime()
               timeexpand=timeexpand+t2-t1
            end do
            done(jorb)=.true.
         end if
      end do
      nrecv=nrecv+ncomplete
      !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
      do i=ind,comon%nrecv-nrecv
         comon%requests(i,2)=comon%requests(i+1,2)
         indexarray(i)=indexarray(i+1)
      end do
      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
      if(nrecv==comon%nrecv) exit waitLoopRecv
   end if
end do waitLoopRecv

!write(*,'(a,i0,a)') 'iproc=',iproc,' is here'

iall=-product(shape(done))*kind(done)
deallocate(done, stat=istat)
call memocc(istat, iall, 'done', subname)

iall=-product(shape(indcomplete))*kind(indcomplete)
deallocate(indcomplete, stat=istat)
call memocc(istat, iall, 'indcomplete', subname)

iall=-product(shape(indexarray))*kind(indexarray)
deallocate(indexarray, stat=istat)
call memocc(istat, iall, 'indexarray', subname)

!!allocate(indcomplete(comon%nsend), stat=istat)
!!call memocc(istat, indcomplete, 'indcomplete', subname)

!!do iorb=1,orbs%norb
!!  do jorb=1,orbs%norb
!!    write(400+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!  end do
!!end do

!!! Now the sends
!!nsend=0
!!waitLoopSend: do
!!    call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!    nsend=nsend+ncomplete
!!    if(nsend==comon%nsend) exit waitLoopSend
!!end do waitLoopSend

!write(*,*) 'here: iproc', iproc
!call mpi_barrier(mpi_comm_world, ierr)
!call mpi_finalize(ierr)
!stop

!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!deallocate(indcomplete, stat=istat)
!!call memocc(istat, iall, 'indcomplete', subname)

!do iorb=1,orbs%norb
!  do jorb=1,orbs%norb
!    write(90+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!  end do
!end do


! Communicate the matrix
allocate(ovrlpCompressed(mad%nvctr), stat=istat)
call memocc(istat, ovrlpCompressed, 'ovrlpCompressed', subname)
allocate(sendcounts(0:nproc-1), stat=istat)
call memocc(istat, sendcounts, 'sendcounts', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)

t1=mpi_wtime()
call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed, sendcounts, displs)
t2=mpi_wtime()
timecompress=timecompress+t2-t1     

allocate(ovrlpCompressed2(mad%nvctr), stat=istat)
call memocc(istat, ovrlpCompressed2, 'ovrlpCompressed2', subname)

t1=mpi_wtime()
if (nproc >1) then
   call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
        sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
else
   call vcopy(sendcounts(iproc),ovrlpCompressed(displs(iproc)+1),1,ovrlpCompressed2(1+displs(iproc)),1)
end if

t2=mpi_wtime()
timecommuncoll=timecommuncoll+t2-t1     

t1=mpi_wtime()
call uncompressMatrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)
t2=mpi_wtime()
timecompress=timecompress+t2-t1     

!do iorb=1,orbs%norb
!  do jorb=1,orbs%norb
!    write(100+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!  end do
!end do

iall=-product(shape(ovrlpCompressed))*kind(ovrlpCompressed)
deallocate(ovrlpCompressed, stat=istat)
call memocc(istat, iall, 'ovrlpCompressed', subname)
iall=-product(shape(ovrlpCompressed2))*kind(ovrlpCompressed2)
deallocate(ovrlpCompressed2, stat=istat)
call memocc(istat, iall, 'ovrlpCompressed2', subname)
iall=-product(shape(sendcounts))*kind(sendcounts)
deallocate(sendcounts, stat=istat)
call memocc(istat, iall, 'sendcounts', subname)
iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)


!!call mpi_barrier(mpi_comm_world, ierr)
!!call mpi_finalize(ierr)

end subroutine collectAndCalculateOverlap





subroutine collectnew(iproc, nproc, comon, mad, op, orbs, input, lzd, &
   nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, timecommunp2p, timecommuncoll, timecompress)
use module_base
use module_types
use module_interfaces, exceptThisOne => collectnew
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
type(p2pCommsOrthonormality),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
type(overlapParameters),intent(in):: op
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(nsendbuf),intent(in):: sendbuf
real(8),dimension(nrecvbuf),intent(inout):: recvbuf
real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
real(8),intent(inout):: timecommunp2p, timecommuncoll, timecompress

! Local variables
integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
real(8):: ddot, t1, t2
character(len=*),parameter:: subname='collectAndCalculateOverlap'
logical,dimension(:),allocatable:: done
logical:: received

allocate(indcomplete(comon%nrecv), stat=istat)
call memocc(istat, indcomplete, 'indcomplete', subname)
allocate(done(comon%nrecv), stat=istat)
call memocc(istat, done, 'done', subname)
done=.false.

allocate(indexarray(comon%nrecv), stat=istat)
call memocc(istat, indexarray, 'indexarray', subname)
do i=1,comon%nrecv
   indexarray(i)=i
end do


! Wait for the sends to complete.
if (nproc > 1) then
   t1=mpi_wtime()
   nsend=0
   waitLoopSend: do
      !!call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
      !!nsend=nsend+ncomplete
      !!if(nsend==comon%nsend) exit waitLoopSend
      call mpi_waitany(comon%nsend-nsend, comon%requests(1,1), ind, mpi_status_ignore, ierr)
      nsend=nsend+1
      do i=ind,comon%nsend-nsend
         comon%requests(i,1)=comon%requests(i+1,1)
      end do
      if(nsend==comon%nsend) exit waitLoopSend
   end do waitLoopSend
   t2=mpi_wtime()
   timecommunp2p=timecommunp2p+t2-t1


   ovrlp=0.d0
   nrecv=0
   waitLoopRecv: do
      t1=mpi_wtime()
      call mpi_waitany(comon%nrecv-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
      !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
      !ind=1
      t2=mpi_wtime()
      timecommunp2p=timecommunp2p+t2-t1
      ncomplete=1
      received=.true.
      if(received) then
         nrecv=nrecv+ncomplete
         !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
         !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
         do i=ind,comon%nrecv-nrecv
            comon%requests(i,2)=comon%requests(i+1,2)
            indexarray(i)=indexarray(i+1)
         end do
         !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
         if(nrecv==comon%nrecv) exit waitLoopRecv
      end if
   end do waitLoopRecv
end if
!write(*,'(a,i0,a)') 'iproc=',iproc,' is here'

iall=-product(shape(done))*kind(done)
deallocate(done, stat=istat)
call memocc(istat, iall, 'done', subname)

iall=-product(shape(indcomplete))*kind(indcomplete)
deallocate(indcomplete, stat=istat)
call memocc(istat, iall, 'indcomplete', subname)

iall=-product(shape(indexarray))*kind(indexarray)
deallocate(indexarray, stat=istat)
call memocc(istat, iall, 'indexarray', subname)

!!allocate(indcomplete(comon%nsend), stat=istat)
!!call memocc(istat, indcomplete, 'indcomplete', subname)

!!do iorb=1,orbs%norb
!!  do jorb=1,orbs%norb
!!    write(400+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!  end do
!!end do

!!! Now the sends
!!nsend=0
!!waitLoopSend: do
!!    call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!    nsend=nsend+ncomplete
!!    if(nsend==comon%nsend) exit waitLoopSend
!!end do waitLoopSend

!write(*,*) 'here: iproc', iproc
!call mpi_barrier(mpi_comm_world, ierr)
!call mpi_finalize(ierr)
!stop

!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!deallocate(indcomplete, stat=istat)
!!call memocc(istat, iall, 'indcomplete', subname)

!do iorb=1,orbs%norb
!  do jorb=1,orbs%norb
!    write(90+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!  end do
!end do


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
!!call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
!!     sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
!!t2=mpi_wtime()
!!timecommuncoll=timecommuncoll+t2-t1     
!!
!!t1=mpi_wtime()
!!call uncompressMatrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)
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


!!call mpi_barrier(mpi_comm_world, ierr)
!!call mpi_finalize(ierr)

end subroutine collectnew





subroutine collectnew2(iproc, nproc, comon, mad, op, orbs, input, lzd, &
   nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, timecommunp2p, timecommuncoll, timecompress)
use module_base
use module_types
use module_interfaces, exceptThisOne => collectnew
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
type(p2pCommsOrthonormality),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
type(overlapParameters),intent(in):: op
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(nsendbuf),intent(in):: sendbuf
real(8),dimension(nrecvbuf),intent(inout):: recvbuf
real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
real(8),intent(inout):: timecommunp2p, timecommuncoll, timecompress

! Local variables
integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
real(8):: ddot, t1, t2
character(len=*),parameter:: subname='collectAndCalculateOverlap'
logical,dimension(:),allocatable:: done
logical:: received


allocate(indexarray(comon%nrecvtemp), stat=istat)
call memocc(istat, indexarray, 'indexarray', subname)
do i=1,comon%nrecvtemp
   indexarray(i)=i
end do

if (nproc > 1) then
   ! Wait for the sends to complete.
   t1=mpi_wtime()
   nsend=0
   if(comon%nsendtemp/=0) then
      waitLoopSend: do
         !!call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
         !!nsend=nsend+ncomplete
         !!if(nsend==comon%nsend) exit waitLoopSend
         call mpi_waitany(comon%nsendtemp-nsend, comon%requests(1,1), ind, mpi_status_ignore, ierr)
         nsend=nsend+1
         do i=ind,comon%nsendtemp-nsend
            comon%requests(i,1)=comon%requests(i+1,1)
         end do
         if(nsend==comon%nsendtemp) exit waitLoopSend
      end do waitLoopSend
   end if
   t2=mpi_wtime()
   timecommunp2p=timecommunp2p+t2-t1



   ovrlp=0.d0
   nrecv=0
   if(comon%nrecvtemp/=0) then
      waitLoopRecv: do
         t1=mpi_wtime()
         call mpi_waitany(comon%nrecvtemp-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
         !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
         !ind=1
         t2=mpi_wtime()
         timecommunp2p=timecommunp2p+t2-t1
         ncomplete=1
         received=.true.
         if(received) then
            nrecv=nrecv+ncomplete
            !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
            !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
            do i=ind,comon%nrecvtemp-nrecv
               comon%requests(i,2)=comon%requests(i+1,2)
               indexarray(i)=indexarray(i+1)
            end do
            !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
            if(nrecv==comon%nrecvtemp) exit waitLoopRecv
         end if
      end do waitLoopRecv
   end if
end if

iall=-product(shape(indexarray))*kind(indexarray)
deallocate(indexarray, stat=istat)
call memocc(istat, iall, 'indexarray', subname)


end subroutine collectnew2








subroutine collectAndCalculateOverlap2(iproc, nproc, comon, mad, op, orbs, input, lzd, &
   nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, timecommunp2p, timecommuncoll, timeoverlap, timecompress)
use module_base
use module_types
use module_interfaces, exceptThisOne => collectAndCalculateOverlap2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
type(p2pCommsOrthonormality),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
type(overlapParameters),intent(in):: op
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(nsendbuf),intent(in):: sendbuf
real(8),dimension(nrecvbuf),intent(inout):: recvbuf
real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
real(8),intent(inout):: timecommunp2p, timecommuncoll, timeoverlap, timecompress

! Local variables
integer:: iorb, orbsource, orbdest, nrecv, nsend, ist, jst, ncount, ierr, ncomplete, i, istat, iall, jorb, ind
real(8),dimension(:),allocatable:: ovrlpCompressed, ovrlpCompressed2, sendbuf2, temparr
integer,dimension(:),allocatable:: sendcounts, displs, indcomplete, indexarray
real(8):: ddot, t1, t2
character(len=*),parameter:: subname='collectAndCalculateOverlap'
logical,dimension(:),allocatable:: done
logical:: received

allocate(indcomplete(comon%nrecv), stat=istat)
call memocc(istat, indcomplete, 'indcomplete', subname)
allocate(done(comon%nrecv), stat=istat)
call memocc(istat, done, 'done', subname)
done=.false.

allocate(indexarray(comon%nrecv), stat=istat)
call memocc(istat, indexarray, 'indexarray', subname)
do i=1,comon%nrecv
   indexarray(i)=i
end do



ovrlp=0.d0
nrecv=0
waitLoopRecv: do
   if (nproc > 1) then
      t1=mpi_wtime()
      call mpi_waitany(comon%nrecv-nrecv, comon%requests(1,2), ind, mpi_status_ignore, ierr)
      !call mpi_testany(comon%nrecv-nrecv, comon%requests(1,2), ind, received, mpi_status_ignore, ierr)
      !ind=1
      t2=mpi_wtime()
      timecommunp2p=timecommunp2p+t2-t1
   end if
   ncomplete=1
   received=.true.
   if(received) then
      do i=1,ncomplete
         ! Calculate overlap matrix
         !jorb=indcomplete(i)
         jorb=indexarray(ind)
         !!write(*,'(2(a,i0))') 'process ',iproc,' has received message ',jorb
         if(.not.done(jorb)) then
            orbsource=comon%comarr(9,jorb,iproc)
            do iorb=1,orbs%norbp
               orbdest=orbs%isorb+iorb
               call getStartingIndicesGlobal(orbdest, orbsource, op, orbs, ist, jst, ncount)
               !jst=comon%comarr(5,jorb,iproc)
               !ncount=comon%comarr(3,jorb,iproc)
               !write(*,'(a,i4,2i5,3x,2i8,4x,2i8)') 'iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc) ; jst, comon%comarr(5,jorb,iproc)', iproc, ind, jorb, ncount, comon%comarr(3,jorb,iproc), jst, comon%comarr(5,jorb,iproc)
               t1=mpi_wtime()
               ovrlp(orbsource,orbdest)=ddot(ncount, sendBuf(ist), 1, recvBuf(jst), 1)
               t2=mpi_wtime()
               timeoverlap=timeoverlap+t2-t1
               !!t1=mpi_wtime()
               !!call expandoneorbital2(iproc, nproc, orbsource, orbdest-orbs%isorb, orbs, input, &
               !!     orbs%inwhichlocreg, lzd, op, nrecvbuf, recvbuf, lphiovrlp)
               !!t2=mpi_wtime()
               !!timeexpand=timeexpand+t2-t1
            end do
            done(jorb)=.true.
         end if
      end do
      nrecv=nrecv+ncomplete
      !write(*,'(5(a,i0))') 'iproc=',iproc,': communication ',ind,' corresponding to jorb=',jorb,') has completed; moving requests from ',ind,' to ',comon%nrecv-nrecv
      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests before: ',comon%requests(1:comon%nrecv,2)
      do i=ind,comon%nrecv-nrecv
         comon%requests(i,2)=comon%requests(i+1,2)
         indexarray(i)=indexarray(i+1)
      end do
      !write(*,'(a,i0,a,4x,40i7)') 'iproc=',iproc,': requests after: ',comon%requests(1:comon%nrecv,2)
      if(nrecv==comon%nrecv) exit waitLoopRecv
   end if
end do waitLoopRecv

!write(*,'(a,i0,a)') 'iproc=',iproc,' is here'

iall=-product(shape(done))*kind(done)
deallocate(done, stat=istat)
call memocc(istat, iall, 'done', subname)

iall=-product(shape(indcomplete))*kind(indcomplete)
deallocate(indcomplete, stat=istat)
call memocc(istat, iall, 'indcomplete', subname)

iall=-product(shape(indexarray))*kind(indexarray)
deallocate(indexarray, stat=istat)
call memocc(istat, iall, 'indexarray', subname)

!!allocate(indcomplete(comon%nsend), stat=istat)
!!call memocc(istat, indcomplete, 'indcomplete', subname)

!!do iorb=1,orbs%norb
!!  do jorb=1,orbs%norb
!!    write(400+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!  end do
!!end do

!!! Now the sends
!!nsend=0
!!waitLoopSend: do
!!    call mpi_waitsome(comon%nsend, comon%requests(1,1), ncomplete, indcomplete, mpi_statuses_ignore, ierr)
!!    nsend=nsend+ncomplete
!!    if(nsend==comon%nsend) exit waitLoopSend
!!end do waitLoopSend

!write(*,*) 'here: iproc', iproc
!call mpi_barrier(mpi_comm_world, ierr)
!call mpi_finalize(ierr)
!stop

!!iall=-product(shape(indcomplete))*kind(indcomplete)
!!deallocate(indcomplete, stat=istat)
!!call memocc(istat, iall, 'indcomplete', subname)

!do iorb=1,orbs%norb
!  do jorb=1,orbs%norb
!    write(90+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!  end do
!end do


! Communicate the matrix
allocate(ovrlpCompressed(mad%nvctr), stat=istat)
call memocc(istat, ovrlpCompressed, 'ovrlpCompressed', subname)
allocate(sendcounts(0:nproc-1), stat=istat)
call memocc(istat, sendcounts, 'sendcounts', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)

t1=mpi_wtime()
call compressMatrix2(iproc, nproc, orbs, mad, ovrlp, ovrlpCompressed, sendcounts, displs)
t2=mpi_wtime()
timecompress=timecompress+t2-t1     

allocate(ovrlpCompressed2(mad%nvctr), stat=istat)
call memocc(istat, ovrlpCompressed2, 'ovrlpCompressed2', subname)

t1=mpi_wtime()
if (nproc > 1) then
   call mpi_allgatherv(ovrlpCompressed(displs(iproc)+1), sendcounts(iproc), mpi_double_precision, ovrlpCompressed2(1), &
        sendcounts, displs, mpi_double_precision, mpi_comm_world, ierr)
else
   call vcopy(sendcounts(iproc),ovrlpCompressed(displs(iproc)+1),1,ovrlpCompressed2(1+displs(iproc)),1)
end if
t2=mpi_wtime()
timecommuncoll=timecommuncoll+t2-t1     

t1=mpi_wtime()
call uncompressMatrix(orbs%norb, mad, ovrlpCompressed2, ovrlp)
t2=mpi_wtime()
timecompress=timecompress+t2-t1     

!do iorb=1,orbs%norb
!  do jorb=1,orbs%norb
!    write(100+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!  end do
!end do

iall=-product(shape(ovrlpCompressed))*kind(ovrlpCompressed)
deallocate(ovrlpCompressed, stat=istat)
call memocc(istat, iall, 'ovrlpCompressed', subname)
iall=-product(shape(ovrlpCompressed2))*kind(ovrlpCompressed2)
deallocate(ovrlpCompressed2, stat=istat)
call memocc(istat, iall, 'ovrlpCompressed2', subname)
iall=-product(shape(sendcounts))*kind(sendcounts)
deallocate(sendcounts, stat=istat)
call memocc(istat, iall, 'sendcounts', subname)
iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)


!!call mpi_barrier(mpi_comm_world, ierr)
!!call mpi_finalize(ierr)

end subroutine collectAndCalculateOverlap2

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
!type(communications_arrays),intent(in):: comms
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
!  call transpose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, phi, work=phiWork)
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





subroutine applyOrthoconstraintNonorthogonalCubic(iproc, nproc, methTransformOverlap,&
     blocksize_pdgemm, orbs, gorbs, comms, lzd, input, &
     op, ovrlp, mad, lphi, lhphi, trH)
use module_base
use module_types
use module_interfaces, exceptThisOne => applyOrthoconstraintNonorthogonalCubic
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm
type(orbitals_data),intent(in):: orbs, gorbs
type(communications_arrays),intent(in):: comms
type(local_zone_descriptors),intent(in):: lzd
type(input_variables),intent(in):: input
type(overlapParameters),intent(in):: op
real(8),dimension(orbs%norb,orbs%norb),intent(inout):: ovrlp
type(matrixDescriptors),intent(in):: mad
real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lphi, lhphi
real(8),intent(out):: trH

! Local variables
integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, info, i, istat, iall, ierr, iseg, nvctrp, ldim, gdim, ind1, ind2
real(8):: tt, t1, t2, time_dsymm, time_daxpy
real(8),dimension(:,:),allocatable:: ovrlp2
real(8),dimension(:,:),allocatable:: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
character(len=*),parameter:: subname='applyOrthoconstraintNonorthogonalCubic'
logical:: val, segment
real(8),dimension(:,:),allocatable:: lagmat
real(8),dimension(:),allocatable:: phi, hphi
real(8),dimension(:),pointer:: phiWork

allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, lagmat, 'lagmat', subname)

allocate(phi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)), stat=istat)
call memocc(istat, phi, 'phi', subname)
allocate(hphi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)), stat=istat)
call memocc(istat, hphi, 'hphi', subname)
allocate(phiWork(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)), stat=istat)
call memocc(istat, phiWork, 'phiWork', subname)

! Calculate the overlap matrix
  ind1=1
  ind2=1
  phi=0.d0
  do iorb=1,orbs%norbp
      !ilr = orbs%inWhichLocregp(iorb)
      iiorb=orbs%isorb+iorb
      ilr = orbs%inWhichLocreg(iiorb)
      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norb, orbs%nspinor, input%nspin, lzd%Glr,&
           lzd%Llr(ilr), lphi(ind2), phi(ind1))
      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
  end do
  call transpose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, phi, work=phiWork)


  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi, nvctrp, phi, nvctrp, 0.d0, ovrlp, orbs%norb)

  call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)


  ! Now calculate the Lagrange multiplier matrix
  ind1=1
  ind2=1
  hphi=0.d0
  do iorb=1,orbs%norbp
      !ilr = orbs%inWhichLocregp(iorb)
      iiorb=orbs%isorb+iorb
      ilr = orbs%inWhichLocreg(iiorb)
      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norb, orbs%nspinor, input%nspin, lzd%Glr,&
           lzd%Llr(ilr), lhphi(ind2), hphi(ind1))
      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
  end do
  call transpose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, hphi, work=phiWork)

  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi, nvctrp, hphi, nvctrp, 0.d0, lagmat, orbs%norb)

  call mpiallred(lagmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
  trH=0.d0
  do iorb=1,orbs%norb
      trH=trH+lagmat(iorb,iorb)
  end do





allocate(ovrlp_minus_one_lagmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
allocate(ovrlp_minus_one_lagmat_trans(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp2, 'ovrlp2', subname)

call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)

! Invert the overlap matrix
call mpi_barrier(mpi_comm_world, ierr)
call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, orbs%norb, mad, orbs, ovrlp2)


! Multiply the Lagrange multiplier matrix with S^-1/2.
! First fill the upper triangle.
do iorb=1,orbs%norb
    do jorb=1,iorb-1
        ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
    end do
end do
call cpu_time(t1)
if(blocksize_pdgemm<0) then
    !call dsymm('l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
    !     0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
    ovrlp_minus_one_lagmat=0.d0
    call dgemm_compressed2(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
         mad%keygmatmul, ovrlp2, lagmat, ovrlp_minus_one_lagmat)
    ! Transpose lagmat
    do iorb=1,orbs%norb
        do jorb=iorb+1,orbs%norb
            tt=lagmat(jorb,iorb)
            lagmat(jorb,iorb)=lagmat(iorb,jorb)
            lagmat(iorb,jorb)=tt
        end do
    end do
    !call dsymm('l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
    !     0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
    ovrlp_minus_one_lagmat_trans=0.d0
    call dgemm_compressed2(iproc, nproc, orbs%norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, &
         mad%keygmatmul, ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
else
    call dsymm_parallel(iproc, nproc, blocksize_pdgemm, mpi_comm_world, 'l', 'l', orbs%norb, orbs%norb, 1.d0, &
         ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, 0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
    ! Transpose lagmat
    do iorb=1,orbs%norb
        do jorb=iorb+1,orbs%norb
            tt=lagmat(jorb,iorb)
            lagmat(jorb,iorb)=lagmat(iorb,jorb)
            lagmat(iorb,jorb)=tt
        end do
    end do
    call dsymm_parallel(iproc, nproc, blocksize_pdgemm, mpi_comm_world, 'l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), &
         orbs%norb, lagmat(1,1), orbs%norb, &
         0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
end if
call cpu_time(t2)
time_dsymm=t2-t1


call dgemm('n', 'n', nvctrp, orbs%norb, orbs%norb, -.5_wp, phi, nvctrp, lagmat, orbs%norb, 1.0_wp, hphi, nvctrp)
call dgemm('n', 't', nvctrp, orbs%norb, orbs%norb, -.5_wp, phi, nvctrp, lagmat, orbs%norb, 1.0_wp, hphi, nvctrp)

  call untranspose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, hphi, work=phiWork)

  ind1=1
  ind2=1
  do iorb=1,orbs%norbp
      !ilr = orbs%inWhichLocregp(iorb)
      iiorb=orbs%isorb+iorb
      ilr = orbs%inWhichLocreg(iiorb)
      ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      gdim=lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
      call psi_to_locreg2(iproc, nproc, ldim, gdim, lzd%llr(ilr), lzd%glr, hphi(ind1), lhphi(ind2))
      ind1=ind1+lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
      ind2=ind2+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
  end do


call cpu_time(t1)
!!ist=1
!!jst=1
!!ilrold=-1
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=onWhichAtom(iiorb)
!!    if(ilr==ilrold) then
!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!        jst=jst-op%noverlaps(iiorb)*ncount
!!    end if
!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    !call dscal(ncount, 1.5d0, lhphi(ist), 1)
!!    do jorb=1,op%noverlaps(iiorb)
!!        jjorb=op%overlaps(jorb,iiorb)
!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
!!        jst=jst+ncount
!!    end do
!!    ist=ist+ncount
!!    ilrold=ilr
!!end do
call cpu_time(t2)
time_daxpy=t2-t1

if (verbose > 2) then
   call mpiallred(time_dsymm, 1, mpi_sum, mpi_comm_world, ierr)
   call mpiallred(time_daxpy, 1, mpi_sum, mpi_comm_world, ierr)
   if(iproc==0) write(*,'(a,es15.6)') 'time for dsymm',time_dsymm/dble(nproc)
   if(iproc==0) write(*,'(a,es15.6)') 'time for daxpy',time_daxpy/dble(nproc)
end if


iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
deallocate(ovrlp_minus_one_lagmat, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
iall=-product(shape(ovrlp2))*kind(ovrlp2)
deallocate(ovrlp2, stat=istat)
call memocc(istat, iall, 'ovrlp2', subname)
iall=-product(shape(lagmat))*kind(lagmat)
deallocate(lagmat, stat=istat)
call memocc(istat, iall, 'lagmat', subname)
iall=-product(shape(phi))*kind(phi)
deallocate(phi, stat=istat)
call memocc(istat, iall, 'phi', subname)
iall=-product(shape(hphi))*kind(hphi)
deallocate(hphi, stat=istat)
call memocc(istat, iall, 'hphi', subname)
iall=-product(shape(phiWork))*kind(phiWork)
deallocate(phiWork, stat=istat)
call memocc(istat, iall, 'phiWork', subname)


end subroutine applyOrthoconstraintNonorthogonalCubic






subroutine localloewdin(iproc, nproc, lzd, orbs, op, ovrlp, lphiovrlp, lphi)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOneHalf
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(in):: op
  real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
  real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
  real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
  
  ! Local variables
  integer:: lwork, istat, iall, iorb, jorb, info, iiorb, korb, ilr, jlr, norbinlocreg, iorbinlocreg, jorbinlocreg
  integer:: iiorbinlocreg, jjorbinlocreg, klr, ist, jst, korbinlocreg, iilr, ilrold, ncount, jjorb, jjlr
  character(len=*),parameter:: subname='localloewdin'
  real(8),dimension(:),allocatable:: eval, work
  real(8),dimension(:,:),allocatable:: ovrlp2, ovrlp3, ovrlplocal
  real(8),dimension(:,:,:),allocatable:: tempArr
  real(8):: tt, dnrm2, ddot


  lphi=0.d0

  do ilr=1,lzd%nlr

      ! Determine number of orbitals in this locreg
      norbinlocreg=0
      do jorb=1,orbs%norb
          jlr=orbs%inwhichlocreg(jorb)
          if(jlr==ilr) norbinlocreg=norbinlocreg+1
      end do

      ! Exctract the corresponding part of the overlap matrix
      allocate(ovrlplocal(norbinlocreg,norbinlocreg), stat=istat)
      call memocc(istat, ovrlplocal, 'ovrlplocal', subname)
      jorbinlocreg=0
      do jorb=1,orbs%norb
          jlr=orbs%inwhichlocreg(jorb)
          if(jlr==ilr) then
              jorbinlocreg=jorbinlocreg+1
              korbinlocreg=0
              do korb=1,orbs%norb
                  klr=orbs%inwhichlocreg(korb)
                  if(klr==ilr) then
                      korbinlocreg=korbinlocreg+1
                      ovrlplocal(korbinlocreg,jorbinlocreg)=ovrlp(korb,jorb)
                  end if
              end do
          end if
      end do
  
      ! Exact calculation of ovrlp**(-1/2)
      allocate(eval(norbinlocreg), stat=istat)
      call memocc(istat, eval, 'eval', subname)
      allocate(tempArr(norbinlocreg,norbinlocreg,2), stat=istat)
      call memocc(istat, tempArr, 'tempArr', subname)

      lwork=1000*norbinlocreg
      allocate(work(lwork), stat=istat)
      call memocc(istat, work, 'work', subname)
      call dsyev('v', 'l', norbinlocreg, ovrlplocal(1,1), norbinlocreg, eval, work, lwork, info)
      if(info/=0) then
          write(*,'(a,i0)') 'ERROR in dsyev, info=', info
          stop
      end if
      iall=-product(shape(work))*kind(work)
      deallocate(work, stat=istat)
      call memocc(istat, iall, 'work', subname)
      
      ! Calculate S^{-1/2}. 
      ! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
      ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
      do iorb=1,norbinlocreg
          do jorb=1,norbinlocreg
              !tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
              tempArr(jorb,iorb,1)=ovrlplocal(jorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
          end do
      end do
      
      ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
      ! This will give S^{-1/2}.
      call dgemm('n', 't', norbinlocreg, norbinlocreg, norbinlocreg, 1.d0, ovrlplocal(1,1), &
           norbinlocreg, tempArr(1,1,1), norbinlocreg, 0.d0, tempArr(1,1,2), norbinlocreg)
      call dcopy(norbinlocreg**2, tempArr(1,1,2), 1, ovrlplocal(1,1), 1)
      
      
      iall=-product(shape(eval))*kind(eval)
      deallocate(eval, stat=istat)
      call memocc(istat, iall, 'eval', subname)
      iall=-product(shape(tempArr))*kind(tempArr)
      deallocate(tempArr, stat=istat)
      call memocc(istat, iall, 'tempArr', subname)

      ! Orthogonalize the orbitals in this locreg
      ist=1
      jst=1
      ilrold=-1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          iilr=orbs%inwhichlocreg(iiorb)
          if(iilr==ilr) then
              if(iilr==ilrold) then
                  ! Set back the index of lphiovrlp, since we again need the same orbitals.
                  jst=jst-op%noverlaps(iiorb)*ncount
              end if
              ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
              do jorb=1,op%noverlaps(iiorb)
                  jjorb=op%overlaps(jorb,iiorb)
                  jjlr=orbs%inwhichlocreg(jjorb)
                  if(jjlr==iilr) then
                      ! First determine number of iiorb and jjorb in this locreg
                      iiorbinlocreg=0
                      do korb=1,orbs%norb
                          klr=orbs%inwhichlocreg(korb)
                          if(klr==iilr) iiorbinlocreg=iiorbinlocreg+1
                          if(korb==iiorb) exit
                      end do
                      jjorbinlocreg=0
                      do korb=1,orbs%norb
                          klr=orbs%inwhichlocreg(korb)
                          if(klr==iilr) jjorbinlocreg=jjorbinlocreg+1
                          if(korb==jjorb) exit
                      end do
                      call daxpy(ncount, ovrlplocal(jjorbinlocreg,iiorbinlocreg), lphiovrlp(jst), 1, lphi(ist), 1)
                  end if
                  jst=jst+ncount
              end do
      
              ! Normalize
              tt=dnrm2(ncount, lphi(ist), 1)
              call dscal(ncount, 1/tt, lphi(ist), 1)
      
              ist=ist+ncount
              ilrold=iilr
          else
              ! only increase indices
              if(iilr==ilrold) then
                  ! Set back the index of lphiovrlp, since we again need the same orbitals.
                  jst=jst-op%noverlaps(iiorb)*ncount
              end if
              ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
              do jorb=1,op%noverlaps(iiorb)
                  jst=jst+ncount
              end do
              ist=ist+ncount
              ilrold=iilr
          end if
      end do


      iall=-product(shape(ovrlplocal))*kind(ovrlplocal)
      deallocate(ovrlplocal, stat=istat)
      call memocc(istat, iall, 'ovrlplocal', subname)



  end do


endsubroutine localloewdin














subroutine applyOrthoconstraintlocal(iproc, nproc, lzd, orbs, op, lagmat, lphiovrlp, lhphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => applyOrthoconstraintlocal
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(in):: op
real(8),dimension(orbs%norb,orbs%norb),intent(inout):: lagmat
real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
real(8),dimension(orbs%npsidim_orbs),intent(out):: lhphi

! Local variables
integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, info, i, istat, iall, ierr, iilr
real(8):: tt, t1, t2, time_dsymm, time_daxpy
character(len=*),parameter:: subname='applyOrthoconstraintlocal'


call cpu_time(t1)
do ilr=1,lzd%nlr

    ist=1
    jst=1
    ilrold=-1
    do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        iilr=orbs%inwhichlocreg(iiorb)
        if(iilr==ilr) then
            if(iilr==ilrold) then
                ! Set back the index of lphiovrlp, since we again need the same orbitals.
                jst=jst-op%noverlaps(iiorb)*ncount
            end if
            ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
            !call dscal(ncount, 1.5d0, lhphi(ist), 1)
            do jorb=1,op%noverlaps(iiorb)
                jjorb=op%overlaps(jorb,iiorb)
                call daxpy(ncount, -lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
                jst=jst+ncount
            end do
            ist=ist+ncount
            ilrold=iilr
        else
            ! only increase index
            if(iilr==ilrold) then
                ! Set back the index of lphiovrlp, since we again need the same orbitals.
                jst=jst-op%noverlaps(iiorb)*ncount
            end if
            ncount=lzd%llr(iilr)%wfd%nvctr_c+7*lzd%llr(iilr)%wfd%nvctr_f
            do jorb=1,op%noverlaps(iiorb)
                jst=jst+ncount
            end do
            ist=ist+ncount
            ilrold=iilr
        end if
    end do
end do


call cpu_time(t2)
time_daxpy=t2-t1

call mpiallred(time_daxpy, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(a,es15.6)') 'time for daxpy',time_daxpy/dble(nproc)


end subroutine applyOrthoconstraintlocal




subroutine calculate_overlap_matrix(iproc, orbs, collComms, psi_i, psi_j, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc
type(orbitals_data),intent(in):: orbs
type(collectiveComms),intent(in):: collComms
real(8),dimension(orbs%npsidim_orbs),intent(in):: psi_i, psi_j
real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

! Local variables
integer:: iorb, jorb, ist, jst, ncnt_iorb, ncnt_jorb, iloc, jloc, i, ii, jj, ierr, iist, jjst, istat
logical:: istop, jstop
real(8):: tt

iist=1
do iorb=1,orbs%norb
    jjst=1
    do jorb=1,orbs%norb
       ncnt_iorb=collComms%nvctr_par(iorb,iproc) 
       ncnt_jorb=collComms%nvctr_par(jorb,iproc) 

       ii=0
       jj=0
       istop=.false.
       jstop=.false.
       tt=0.d0
       ist=iist
       jst=jjst
       !do i=1,max(ncnt_iorb,ncnt_jorb)
       !write(*,'(a,2i8,i10,2i12,2l)') 'iproc, jorb, jst, ncnt_jorb, orbs%npsidim, istop, jstop', iproc, jorb, jst, ncnt_jorb, orbs%npsidim, istop, jstop
       if(iproc==0 .and. iorb==1) then
           do istat=iist,iist+ncnt_iorb
               write(110,*) istat, collComms%indexarray(istat)
           end do
       end if
       if(iproc==0 .and. jorb==5) then
           do istat=jjst,jjst+ncnt_jorb
               write(510,*) istat, collComms%indexarray(istat)
           end do
       end if

       do
           !if(jstop .and. jst>orbs%npsidim) write(*,'(a,5i12,2l)') 'iproc, jj, ncnt_jorb, jst, orbs%npsidim, istop, jstop', iproc, jj, ncnt_jorb, jst, orbs%npsidim, istop, jstop
           if(.not.istop) iloc=collComms%indexarray(ist)  
           if(.not.jstop) jloc=collComms%indexarray(jst)  
           write(1000*(iproc+1)+700+jorb,'(a,2i5,2i12,3x,2i12)') 'iorb, jorb, ist, jst, iloc, jloc', &
               iorb, jorb, ist, jst, iloc, jloc
           if(iorb==1 .and. jorb==5) then
               write(980+iproc,'(2i9,4x,2i12,2l4,2i9)') ist, jst, collComms%indexarray(ist), collComms%indexarray(jst), istop, jstop, ncnt_iorb, ncnt_jorb
           end if
           if(iloc==jloc) then
               !if(iorb==1 .and. jorb==5) then
               !    write(980+iproc,'(a,2i9,4x,2i12,2l4,2i9)') 'HERE1: ',ist, jst, collComms%indexarray(ist), collComms%indexarray(jst), istop, jstop, ncnt_iorb, ncnt_jorb
               !end if
               !if(istop) write(*,*) 'STRANGE: istop is true...'
               !if(jstop) write(*,*) 'STRANGE: jstop is true...'
               !if(iorb==1 .and. jorb==5) then
               !    write(iproc*1000+999,*) iloc
               !end if
               tt=tt+psi_i(ist)*psi_j(jst)
               ist=ist+1
               jst=jst+1
               ii=ii+1
               jj=jj+1
               !if(iorb==1 .and. jorb==5) then
               !    write(980+iproc,'(a,2i9,4x,2i12,2l4,2i9)') 'HERE2: ',ist, jst, collComms%indexarray(ist), collComms%indexarray(jst), istop, jstop, ncnt_iorb, ncnt_jorb
               !end if
           else if((iloc<jloc .or. jstop) .and. .not.istop) then
               ist=ist+1
               ii=ii+1
           else if((jloc<iloc .or. istop) .and. .not.jstop) then
               jst=jst+1
               jj=jj+1
           end if
           if(ii==ncnt_iorb) istop=.true.
           if(jj==ncnt_jorb) jstop=.true.
           if(istop .and. jstop) then
               if(iorb==1 .and. jorb==5) write(980+iproc,'(a)') 'exit since both stops are true...'
               exit
           end if
       end do

       ovrlp(jorb,iorb)=tt
       jjst=jjst+ncnt_jorb

    end do
    iist=iist+ncnt_iorb
end do


call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)


end subroutine calculate_overlap_matrix
