subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, nItOrtho, &
           orbs, op, comon, lzd, mad, collcom, orthpar, bpo, lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc,nproc,methTransformOverlap,nItOrtho
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(inout):: op
  type(p2pComms),intent(inout):: comon
  type(local_zone_descriptors),intent(in):: lzd
  type(matrixDescriptors),intent(in):: mad
  type(collective_comms),intent(in):: collcom
  type(orthon_data),intent(in):: orthpar
  type(basis_performance_options),intent(in):: bpo
  real(8),dimension(orbs%npsidim_orbs), intent(inout) :: lphi
  real(8),dimension(:),pointer,intent(out):: psit_c, psit_f
  logical,intent(out):: can_use_transposed

  ! Local variables
  integer:: it, istat, iall, ierr, iorb, jorb, ilr, ncount, ist, iiorb
  real(8),dimension(:),allocatable:: lphiovrlp, psittemp_c, psittemp_f
  character(len=*),parameter:: subname='orthonormalizeLocalized'
  real(8):: maxError
  real(8),dimension(:,:),allocatable:: ovrlp


  allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)

  if(nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'
  can_use_transposed=.false.

  do it=1,nItOrtho

      if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
          allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, psit_c, 'psit_c', subname)
          allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, psit_f, 'psit_f', subname)
          call transpose_localized(iproc, nproc, orbs, collcom, lphi, psit_c, psit_f, lzd)
          call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)
      else if (bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
          ! Allocate the send and receive buffers for the communication.
          call allocateSendBufferOrtho(comon, subname)
          call allocateRecvBufferOrtho(comon, subname)
          ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
          call extractOrbital3(iproc, nproc, orbs, orbs, orbs%npsidim_orbs, lzd, lzd, op, op, &
               lphi, comon%nsendBuf, comon%sendBuf)
          ! Post the send messages.
          call post_p2p_communication(iproc, nproc, comon%nsendbuf, comon%sendbuf, comon%nrecvbuf, comon%recvbuf, comon)
          allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
          call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
          call wait_p2p_communication(iproc, nproc, comon)
          call calculateOverlapMatrix3(iproc, nproc, orbs, op, comon%nsendBuf, &
               comon%sendBuf, comon%nrecvBuf, comon%recvBuf, mad, ovrlp)
          !call checkUnity(iproc, orbs%norb, ovrlp, maxError)
          !if(iproc==0) write(*,*) 'deviation from unity:', maxError
      end if

      call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, methTransformOverlap, orthpar%blocksize_pdsyev, &
          orthpar%blocksize_pdgemm, orbs%norb, mad, ovrlp)

      if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
          allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, psittemp_c, 'psittemp_c', subname)
          allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, psittemp_f, 'psittemp_f', subname)
          call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
          call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)
          call build_linear_combination_transposed(orbs%norb, ovrlp, collcom, psittemp_c, psittemp_f, .true., psit_c, psit_f)
          call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f)
          call untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, lphi, lzd)
          can_use_transposed=.true.
          iall=-product(shape(psittemp_c))*kind(psittemp_c)
          deallocate(psittemp_c, stat=istat)
          call memocc(istat, iall, 'psittemp_c', subname)
          iall=-product(shape(psittemp_f))*kind(psittemp_f)
          deallocate(psittemp_f, stat=istat)
          call memocc(istat, iall, 'psittemp_f', subname)
          !!iall=-product(shape(psit_c))*kind(psit_c)
          !!deallocate(psit_c, stat=istat)
          !!call memocc(istat, iall, 'psit_c', subname)
          !!iall=-product(shape(psit_f))*kind(psit_f)
          !!deallocate(psit_f, stat=istat)
          !!call memocc(istat, iall, 'psit_f', subname)

          !!! Normalize... could this be done in the tranposed layout?
          !!ist=1
          !!do iorb=1,orbs%norbp
          !!   iiorb=orbs%isorb+iorb
          !!   ilr=orbs%inwhichlocreg(iiorb)
          !!   ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f

          !!   ! Normalize
          !!   tt=dnrm2(ncount, lphi(ist), 1)
          !!   call dscal(ncount, 1/tt, lphi(ist), 1)

          !!   ist=ist+ncount
          !!end do
          !!call transpose_localized(iproc, nproc, orbs, collcom, lphi, psit_c, psit_f, lzd)

          !call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f)


      else if (bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
          call globalLoewdin(iproc, nproc, orbs, lzd, op, comon, ovrlp, lphiovrlp, lphi)

          call deallocateSendBufferOrtho(comon, subname)
          call deallocateRecvBufferOrtho(comon, subname)

         iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
         deallocate(lphiovrlp, stat=istat)
         call memocc(istat, iall, 'lphiovrlp', subname)
      end if
  end do


  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)


end subroutine orthonormalizeLocalized




subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, orbs, op, comon, mad, collcom, orthpar, bpo, &
           lphi, lhphi, lagmat, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_Data),intent(in):: orbs
  type(overlapParameters),intent(inout):: op
  type(p2pComms),intent(inout):: comon
  type(matrixDescriptors),intent(in):: mad
  type(collective_comms),intent(in):: collcom
  type(orthon_data),intent(in):: orthpar
  type(basis_performance_options),intent(in):: bpo
  real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lphi !inout due to tranposition...
  real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lhphi
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: lagmat
  real(8),dimension(:),pointer,intent(inout):: psit_c, psit_f
  logical,intent(inout):: can_use_transposed

  ! Local variables
  integer:: istat, iall, ierr, iorb, jorb
  real(8),dimension(:),allocatable:: lphiovrlp, hpsit_c, hpsit_f
  real(8),dimension(:,:),allocatable:: ovrlp2, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
  character(len=*),parameter:: subname='orthoconstraintNonorthogonal'
integer :: i, j
real(8) :: diff_frm_ortho, diff_frm_sym ! lr408
  allocate(ovrlp_minus_one_lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
  allocate(ovrlp_minus_one_lagmat_trans(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
  allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp2, 'ovrlp2', subname)



  if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
      !!write(*,*) 'can_use_transposed',can_use_transposed
      if(.not. can_use_transposed) then
          allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, psit_c, 'psit_c', subname)
          allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, psit_f, 'psit_f', subname)
          !!write(*,*) 'transposing...'
          call transpose_localized(iproc, nproc, orbs, collcom, lphi, psit_c, psit_f, lzd)
      end if
      allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
      allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
      call transpose_localized(iproc, nproc, orbs, collcom, lhphi, hpsit_c, hpsit_f, lzd)
      call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat)
      call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp2)
  else if (bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
      allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
      call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
      !!lphiovrlp=0.d0
      call to_zero(op%ndim_lphiovrlp, lphiovrlp(1))
      call allocateCommuncationBuffersOrtho(comon, subname)
      ! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
      call extractOrbital3(iproc, nproc, orbs, orbs, orbs%npsidim_orbs, lzd, lzd, op, op, &
           lphi, comon%nsendBuf, comon%sendBuf)
      call post_p2p_communication(iproc, nproc, comon%nsendbuf, comon%sendbuf, comon%nrecvbuf, comon%recvbuf, comon)
      call wait_p2p_communication(iproc, nproc, comon)
      call calculateOverlapMatrix3(iproc, nproc, orbs, op, comon%nsendBuf, &
           comon%sendBuf, comon%nrecvBuf, comon%recvBuf, mad, ovrlp2)
      call extractOrbital3(iproc, nproc, orbs, orbs, orbs%npsidim_orbs, lzd, lzd, op, op, &
           lhphi, comon%nsendBuf, comon%sendBuf)
      call calculateOverlapMatrix3(iproc, nproc, orbs, op, comon%nsendBuf, &
           comon%sendBuf, comon%nrecvBuf, comon%recvBuf, mad, lagmat)
 end if

  !if (iproc==0) then ! debugging moved to getlocbasis
  !  open(30,file='ovrlp.dat',status='replace')
  !  diff_frm_ortho = 0.0_dp
  !  diff_frm_sym = 0.0_dp
  !  do i=1,orbs%norb
  !    do j=1,orbs%norb
  !       if (i==j) then
  !         diff_frm_ortho = diff_frm_ortho + abs(1.0_dp - ovrlp2(i,j))
  !       else
  !         diff_frm_ortho = diff_frm_ortho + abs(ovrlp2(i,j))
  !       end if
  !       diff_frm_sym = diff_frm_sym + abs(ovrlp2(i,j) - ovrlp2(j,i))
  !       if (iproc == 0) write(30,*) i,j,ovrlp2(i,j)
  !    end do
  !    write(30,*) ''
  !  end do
  !  close(30)
  !  if (iproc == 0) write(51,*) 'diff from ortho',diff_frm_ortho / (orbs%norb **2),&
  !       diff_frm_sym / (orbs%norb **2)
  !end if



  call applyOrthoconstraintNonorthogonal2(iproc, nproc, orthpar%methTransformOverlap, orthpar%blocksize_pdgemm, 0, &
       orbs, lagmat, ovrlp2, mad, &
       ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
  !!call applyOrthoconstraintNonorthogonal2(iproc, nproc, orthpar%methTransformOverlap, orthpar%blocksize_pdgemm, 1, &
  !!     orbs, lagmat, ovrlp2, mad, &
  !!     ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)


  if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then

      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              ovrlp2(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat(jorb,iorb)
          end do
      end do
      call build_linear_combination_transposed(orbs%norb, ovrlp2, collcom, psit_c, psit_f, .false., hpsit_c, hpsit_f)
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              ovrlp2(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat_trans(jorb,iorb)
          end do
      end do
      call build_linear_combination_transposed(orbs%norb, ovrlp2, collcom, psit_c, psit_f, .false., hpsit_c, hpsit_f)

      !call untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, lphi, lzd)
      call untranspose_localized(iproc, nproc, orbs, collcom, hpsit_c, hpsit_f, lhphi, lzd)
      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)
      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)
      iall=-product(shape(psit_c))*kind(psit_c)
      deallocate(psit_c, stat=istat)
      call memocc(istat, iall, 'psit_c', subname)
      can_use_transposed=.false.
      iall=-product(shape(psit_f))*kind(psit_f)
      deallocate(psit_f, stat=istat)
      call memocc(istat, iall, 'psit_f', subname)
  else if (bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              ovrlp2(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat(jorb,iorb)
          end do
      end do
      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, comon%recvbuf, ovrlp2, .false., lhphi)
      
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              ovrlp2(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat_trans(jorb,iorb)
          end do
      end do
      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, comon%recvbuf, ovrlp2, .false., lhphi)
      iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
      deallocate(lphiovrlp, stat=istat)
      call memocc(istat, iall, 'lphiovrlp', subname)
      call deallocateCommuncationBuffersOrtho(comon, subname)
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

  !!if (verbose > 2) then
  !!   timeComput=timeExtract+timeExpand+timeApply+timecalcmatrix
  !!   timeCommun=timecommunp2p+timecommuncoll
  !!   call mpiallred(timeComput, 1, mpi_sum, mpi_comm_world, ierr)
  !!   call mpiallred(timeCommun, 1, mpi_sum, mpi_comm_world, ierr)
  !!   call mpiallred(timeCalcMatrix, 1, mpi_sum, mpi_comm_world, ierr)
  !!   call mpiallred(timeExpand, 1, mpi_sum, mpi_comm_world, ierr)
  !!   call mpiallred(timeApply, 1, mpi_sum, mpi_comm_world, ierr)
  !!   call mpiallred(timeExtract, 1, mpi_sum, mpi_comm_world, ierr)
  !!   timeComput=timeComput/dble(nproc)
  !!   timeCommun=timeCommun/dble(nproc)
  !!   timeCalcMatrix=timeCalcMatrix/dble(nproc)
  !!   timeExpand=timeExpand/dble(nproc)
  !!   timeApply=timeApply/dble(nproc)
  !!   timeExtract=timeExtract/dble(nproc)
  !!   if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for computation:', timeComput, '=', &
  !!        100.d0*timeComput/(timeComput+timeCommun), '%'
  !!   if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for communication:', timeCommun, '=', &
  !!        100.d0*timeCommun/(timeComput+timeCommun), '%'
  !!   if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for calculating overlap:', timeCalcMatrix, &
  !!        '=', 100.d0*timeCalcMatrix/(timeComput+timeCommun), '%'
  !!   if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for expansion:', timeExpand, '=', &
  !!        100.d0*timeExpand/(timeComput+timeCommun), '%'
  !!   if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for applying orthoconstraint:', timeApply, &
  !!        '=', 100.d0*timeApply/(timeComput+timeCommun), '%'
  !!   if(iproc==0) write(*,'(3x,a,es9.3,a,f5.1,a)') 'time for extract:', timeExtract, '=', &
  !!        100.d0*timeExtract/(timeComput+timeCommun), '%'
  !!end if

end subroutine orthoconstraintNonorthogonal




subroutine getOverlapMatrix2(iproc, nproc, lzd, orbs, comon_lb, op_lb, lphi, mad, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => getOverlapMatrix2
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  type(p2pComms),intent(inout):: comon_lb
  type(overlapParameters),intent(inout):: op_lb
  real(8),dimension(orbs%npsidim_orbs),intent(inout):: lphi
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

  ! Local variables
  character(len=*),parameter:: subname='getOverlapMatrix2'
  real(8):: tt1, tt2, tt3

  call allocateCommuncationBuffersOrtho(comon_lb, subname)
  call extractOrbital3(iproc,nproc,orbs,orbs,orbs%npsidim_orbs,&
       lzd,lzd,op_lb,op_lb,lphi,comon_lb%nsendBuf,comon_lb%sendBuf)
  !call postCommsOverlapNew(iproc,nproc,orbs,op_lb,lzd,lphi,comon_lb,tt1,tt2)
  call post_p2p_communication(iproc, nproc, comon_lb%nsendbuf, comon_lb%sendbuf, &
       comon_lb%nrecvbuf, comon_lb%recvbuf, comon_lb)
  !!call collectnew(iproc,nproc,comon_lb,mad,op_lb,orbs,lzd,comon_lb%nsendbuf,&
  !!     comon_lb%sendbuf,comon_lb%nrecvbuf,comon_lb%recvbuf,tt1,tt2,tt3)
  call wait_p2p_communication(iproc, nproc, comon_lb)
  call calculateOverlapMatrix3(iproc, nproc, orbs, op_lb, comon_lb%nsendBuf, &
       comon_lb%sendBuf, comon_lb%nrecvBuf, comon_lb%recvBuf, mad, ovrlp)
  call deallocateCommuncationBuffersOrtho(comon_lb, subname)




end subroutine getOverlapMatrix2




subroutine initCommsOrtho(iproc, nproc, nspin, hx, hy, hz, lzd, lzdig, orbs, locregShape, op, comon) 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initCommsOrtho
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, nspin
  real(8),intent(in):: hx, hy, hz
  type(local_zone_descriptors),intent(in):: lzd, lzdig
  type(orbitals_data),intent(in):: orbs
  character(len=1),intent(in):: locregShape
  type(overlapParameters),intent(out):: op
  type(p2pComms),intent(out):: comon

  ! Local variables
  integer:: iorb, jorb, iiorb
  integer::  istat, i1, i2, jjorb, nsub, ierr
  character(len=*),parameter:: subname='initCommsOrtho'


  call timing(iproc,'init_commOrtho','ON')

  call nullify_overlapParameters(op)
  call nullify_p2pComms(comon)

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
!  allocate(overlaps_nseg(orbs%norb,orbs%norbp), stat=istat)
!  call memocc(istat, overlaps_nseg, 'overlaps_nseg', subname)


  ! Count how many overlaping regions each orbital / process has.
  if(locregShape=='c') then
     call countOverlaps(iproc, nproc, orbs, lzd, op, comon)
     allocate(op%overlaps(maxval(op%noverlaps),orbs%norb), stat=istat)
     call memocc(istat, op%overlaps, 'op%overlaps', subname)
     call determineOverlaps(iproc, nproc, orbs, lzd, op, comon)
  else if(locregShape=='s') then
     call determine_overlap_from_descriptors(iproc, nproc, orbs, orbs, lzd, lzd, op, comon)
  end if

  ! OLRs NOT NEEDED ANYMORE
  ! Allocate the types describing the overlap localization regions.
  !!op%noverlapsmaxp=maxval(op%noverlaps(orbs%isorb+1:orbs%isorb+orbs%norbp))
  !!allocate(op%olr(op%noverlapsmaxp,orbs%norbp), stat=istat)
  !!do i2=1,orbs%norbp
  !!   do i1=1,op%noverlapsmaxp
  !!      call nullify_locreg_descriptors(op%olr(i1,i2))
  !!   end do
  !!end do

  ! Set the orbital descriptors for the overlap regions.
  !!if(locregShape=='c') then
  !!   call determineOverlapDescriptors(iproc, nproc, orbs, lzd, lzd%Glr, onWhichAtomAll, op)
  !!else if(locregShape=='s') then
  !!   call determineOverlapDescriptorsSphere(iproc, nproc, orbs, lzd, lzd%Glr, onWhichAtomAll, hx, hy, hz, op)
  !!end if

  ! Initialize the communications.
  !!comon%noverlapsmax=maxval(comon%noverlaps)
  allocate(comon%comarr(6,maxval(comon%noverlaps),0:nproc-1), stat=istat)
  call memocc(istat, comon%comarr, 'comon%comarr', subname)
  !!allocate(comon%communComplete(comon%noverlapsmax,0:nproc-1), stat=istat)
  !!call memocc(istat, comon%communComplete, 'comun%communComplete', subname)
  call set_comms_ortho(iproc, nproc, orbs, lzd, op, comon)

  !DON'T need this anymore
  ! Initialize the index arrays for the transformations from overlap region
  ! to ordinary localization region.
  !!allocate(op%indexExpand(comon%nrecvBuf), stat=istat)
  !!call memocc(istat, op%indexExpand, 'op%indexExpand',subname)

  !!allocate(op%expseg(op%noverlapsmaxp,orbs%norbp), stat=istat)
  !!!call memocc(istat, op%expseg,'op%expseg',subname)
  !!
  !!do i1=1,orbs%norbp
  !!    do i2=1,op%noverlapsmaxp
  !!        call nullify_expansionSegments(op%expseg(i2,i1))
  !!    end do
  !!end do
  !!call indicesForExpansion(iproc, nproc, nspin, orbs, onWhichAtomAll, lzd, op, comon)

  !!allocate(op%indexExtract(comon%nsendBuf), stat=istat)
  !!call memocc(istat, op%indexExtract, 'op%indexExtract',subname)
  !!allocate(op%extseg(op%noverlapsmaxp,orbs%norbp), stat=istat)
  !!do i1=1,orbs%norbp
  !!    do i2=1,op%noverlapsmaxp
  !!        call nullify_expansionSegments(op%extseg(i2,i1))
  !!    end do
  !!end do
  !!call indicesForExtraction(iproc, nproc, orbs, orbs%npsidim_orbs, onWhichAtomAll, lzd, op, comon)

  ! Determine the number of non subdiagonals that the overlap matrix / overlap matrix will have.
  op%nsubmax=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jorb=1,op%noverlaps(iiorb)
          jjorb=op%overlaps(jorb,iiorb)
          nsub=jjorb-iiorb
          op%nsubmax=max(op%nsubmax,nsub)
      end do
  end do
  call mpiallred(op%nsubmax, 1, mpi_max, mpi_comm_world, ierr)
  !if(iproc==0) write(*,*) 'op%nsubmax', op%nsubmax

  call timing(iproc,'init_commOrtho','OF')

end subroutine initCommsOrtho




!> Returns the starting and ending indices (on the coarse grid) of a given localization region.
subroutine getIndices(lr, is1, ie1, is2, ie2, is3, ie3)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(in):: lr
  integer,intent(out):: is1, ie1, is2, ie2, is3, ie3

  is1=lr%ns1!+1 !PB: should not have a +1
  ie1=lr%ns1+lr%d%n1
  is2=lr%ns2!+1
  ie2=lr%ns2+lr%d%n2
  is3=lr%ns3!+1
  ie3=lr%ns3+lr%d%n3

end subroutine getIndices




! Count for each orbital and each process the number of overlapping orbitals.
subroutine countOverlaps(iproc, nproc, orbs, lzd, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(out):: op
  type(p2pComms),intent(out):: comon

  ! Local variables
  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, iiorb
   logical :: isoverlap

  iiorb=0
  do jproc=0,nproc-1
     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
     ilrold=-1
     do iorb=1,orbs%norb_par(jproc,0)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=orbs%inwhichlocreg(iiorb)
 !       call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jorb)
!           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
           call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzd%llr(jlr),isoverlap) 
!           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
!           if(ovrlpx .and. ovrlpy .and. ovrlpz) then
            if(isoverlap) then
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

subroutine determineOverlaps(iproc, nproc, orbs, lzd, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(out):: op
  type(p2pComms),intent(out):: comon

  ! Local variables
  integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, iiorb
!  integer ::  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
!  logical:: ovrlpx, ovrlpy, ovrlpz
  logical :: isoverlap

  ! Initialize to some value which will never be used.
  op%overlaps=-1
  !!comon%overlaps=-1

  iiorb=0
  do jproc=0,nproc-1
     ioverlapMPI=0 ! counts the overlaps for the given MPI process.
     ilrold=-1
     do iorb=1,orbs%norb_par(jproc,0)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=orbs%inwhichlocreg(iiorb)
!        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jorb)
!           call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
           call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzd%llr(jlr),isoverlap)
!           ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!           ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!           ovrlpz = ( is3<=je3 .and. ie3>=js3 )
           !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
           !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
!           if(ovrlpx .and. ovrlpy .and. ovrlpz) then
           if(isoverlap) then
              ioverlaporb=ioverlaporb+1
              op%overlaps(ioverlaporb,iiorb)=jorb
              if(ilr/=ilrold) then
                 ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
                 ! would get the same orbitals again. Therefore the counter is not increased
                 ! in that case.
                 ioverlapMPI=ioverlapMPI+1
                 !!comon%overlaps(ioverlapMPI,jproc)=jorb
              end if
           end if
        end do
        !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
        ilrold=ilr
     end do
     !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
  end do

end subroutine determineOverlaps


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

subroutine set_comms_ortho(iproc, nproc, orbs, lzd, op, comon)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(inout):: op
  type(p2pComms),intent(out):: comon

  ! Local variables
  integer:: jproc, iorb, jorb, iiorb, jjorb, mpisource, mpidest, istsource, istdest, ncount, istat, iall, ijorb
  integer:: ilr, ilrold, jprocold, ildim, ierr, isend, irecv, p2p_tag, tag
  integer,dimension(:),allocatable:: istsourceArr, istdestArr
  character(len=*),parameter:: subname='set_comms_ortho'
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



subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: mpisource, mpidest, istsource, istdest, ncount, tag
  integer,dimension(6),intent(out):: comarr


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

subroutine extractOrbital3(iproc, nproc, orbs, orbsig, sizePhi, lzd, lzdig, op, opig, phi, nsendBuf, sendBuf)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, sizePhi
  type(orbitals_data),intent(in):: orbs, orbsig
  type(local_zone_descriptors),intent(in):: lzd, lzdig
  type(overlapParameters),intent(inout):: op, opig
  real(8),dimension(sizePhi),intent(in):: phi
  integer,intent(in):: nsendBuf
  real(8),dimension(nsendBuf),intent(out):: sendBuf


  ! Local variables
  integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
  integer:: i, indSource, jst, istart, iend, ncount, iseg, kstart, kend, start, kold, kseg, knseg

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
  integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(in):: op
  real(8),dimension(nsendBuf),intent(in):: sendBuf
  real(8),dimension(nrecvBuf),intent(in):: recvBuf
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

  ! Local variables
  integer:: iorb, jorb, iiorb, jjorb, ist, jst, ncount, ierr, istat, iall
  real(8):: ddot, tt, ttmax
  real(8),dimension(:),allocatable:: ovrlpCompressed_send, ovrlpCompressed_receive
  character(len=*),parameter:: subname='calculateOverlapMatrix3'
  integer,dimension(:),allocatable:: sendcounts, displs


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
  call uncompressMatrix(orbs%norb, mad, ovrlpCompressed_receive, ovrlp)
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
  integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
  type(orbitals_data),intent(in):: orbs
  type(overlapParameters),intent(in):: op
  real(8),dimension(nsendBuf),intent(in):: sendBuf
  real(8),dimension(nrecvBuf),intent(in):: recvBuf
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

  ! Local variables
  integer:: iorb, jorb, iiorb, jjorb, ist, jst, ncount
  real(8):: ddot
  character(len=*),parameter:: subname='calculateOverlapMatrix3'


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

subroutine globalLoewdin(iproc, nproc, orbs, lzd, op, comon, ovrlp, lphiovrlp, lphi)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(overlapParameters),intent(in):: op
  type(p2pComms),intent(in):: comon
  real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
  real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
  real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, ist, ncount
  real(8):: tt, dnrm2


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

subroutine applyOrthoconstraintNonorthogonal2(iproc, nproc, methTransformOverlap, blocksize_pdgemm, &
           correction_orthoconstraint, orbs, &
           lagmat, ovrlp, mad, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => applyOrthoconstraintNonorthogonal2
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm, correction_orthoconstraint
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
  real(8),dimension(orbs%norb,orbs%norb),intent(inout):: lagmat
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, info, i, istat, iall, ierr, iseg
  real(8):: tt, t1, t2, time_dsymm, time_daxpy
  real(8),dimension(:,:),allocatable:: ovrlp2
  character(len=*),parameter:: subname='applyOrthoconstraintNonorthogonal2'

  call timing(iproc,'lagmat_orthoco','ON')

  allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp2, 'ovrlp2', subname)

correctionIf: if(correction_orthoconstraint==0) then

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
  if(blocksize_pdgemm<0) then
     !! ATTENTION: HERE IT IS ASSUMED THAT THE INVERSE OF THE OVERLAP MATRIX HAS THE SAME SPARSITY
     !! AS THE OVERLAP MATRIX ITSELF. CHECK THIS!!

     !call dsymm('l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
     !     0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
     t1=mpi_wtime()
     !!ovrlp_minus_one_lagmat=0.d0
     call to_zero(orbs%norb**2, ovrlp_minus_one_lagmat(1,1))
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
    !!ovrlp_minus_one_lagmat_trans=0.d0
    call to_zero(orbs%norb**2, ovrlp_minus_one_lagmat_trans(1,1))
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

else if(correction_orthoconstraint==1) then correctionIf
    do iorb=1,orbs%norb
        do jorb=1,orbs%norb
            ovrlp_minus_one_lagmat(jorb,iorb)=lagmat(jorb,iorb)
            ovrlp_minus_one_lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
        end do
    end do
end if correctionIf

  call timing(iproc,'lagmat_orthoco','OF')




!!if (verbose > 2) then
!!   call mpiallred(time_dsymm, 1, mpi_sum, mpi_comm_world, ierr)
!!   call mpiallred(time_daxpy, 1, mpi_sum, mpi_comm_world, ierr)
!!   if(iproc==0) write(*,'(a,es15.6)') 'time for dsymm',time_dsymm/dble(nproc)
!!   if(iproc==0) write(*,'(a,es15.6)') 'time for daxpy',time_daxpy/dble(nproc)
!!end if


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

  call timing(iproc,'lovrlp^-1     ','ON')

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
  else
      write(*,'(1x,a)') 'ERROR: iorder must be 0,1 or 2!'
      stop
end if

  call timing(iproc,'lovrlp^-1     ','OF')

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
  real(8),dimension(:,:), allocatable :: vr,vl ! for non-symmetric LAPACK
  real(8),dimension(:),allocatable:: eval1 ! for non-symmetric LAPACK
  call timing(iproc,'lovrlp^-1/2   ','ON')
  
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
          !call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, lwork, info)
          ! lr408 - see if LAPACK is stil to blame for convergence issues
          allocate(vl(1:norb,1:norb))
          allocate(vr(1:norb,1:norb))
          allocate(eval1(1:norb))
          call DGEEV( 'v','v', norb, ovrlp(1,1), norb, eval, eval1, VL, norb, VR,&
               norb, WORK, LWORK, info )
          ovrlp=vl
          deallocate(eval1)
          deallocate(vr)
          deallocate(vl)
          !  lr408 - see if LAPACK is stil to blame for convergence issues
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
  else

      write(*,'(1x,a)') 'ERROR: methTransformOrder must be 0,1 or 2!'
      stop

end if

call timing(iproc,'lovrlp^-1/2   ','OF')

end subroutine overlapPowerMinusOneHalf

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
type(p2pComms),intent(inout):: comon_lb
real(8),dimension(orbs%npsidim_orbs),intent(in):: lphi, lhphi
type(matrixDescriptors),intent(in):: mad
real(8),dimension(orbs%norb,orbs%norb),intent(out):: matrixElements

! Local variables
integer:: it, istat, iall, iorb
real(8),dimension(:),allocatable:: lphiovrlp
character(len=*),parameter:: subname='getMatrixElements2'
real(8):: tt1, tt2, tt3
type(input_variables):: input


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
              ncount=op%wfd_overlap(jorb,iorb)%nvctr_c+7*op%wfd_overlap(jorb,iorb)%nvctr_f
          end if
      end do
  end do

end subroutine getStartingIndicesGlobal


