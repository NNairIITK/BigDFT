!> @file
!! Orthonormalization
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, nItOrtho, &
           orbs, op, comon, lzd, mad, collcom, orthpar, bpo, lphi, psit_c, psit_f, &
           can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,nItOrtho
  type(orbitals_data),intent(in) :: orbs
  type(overlapParameters),intent(inout) :: op
  type(p2pComms),intent(inout) :: comon
  type(local_zone_descriptors),intent(in) :: lzd
  type(matrixDescriptors),intent(in) :: mad
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  type(basis_performance_options),intent(in) :: bpo
  real(kind=8),dimension(orbs%npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(out):: can_use_transposed

  ! Local variables
  integer :: it, istat, iall
  !integer :: ilr, iorb, i, jlr, jorb, j
  real(kind=8),dimension(:),allocatable :: lphiovrlp, psittemp_c, psittemp_f
  character(len=*),parameter :: subname='orthonormalizeLocalized'
  !real(kind=8) :: maxError
  real(kind=8),dimension(:,:),allocatable :: ovrlp


  allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)
  if(nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'
  !can_use_transposed=.false.
  do it=1,nItOrtho

      if(.not.can_use_transposed) then
          if(associated(psit_c)) then
              iall=-product(shape(psit_c))*kind(psit_c)
              deallocate(psit_c, stat=istat)
              call memocc(istat, iall, 'psit_c', subname)
          end if
          if(associated(psit_f)) then
              iall=-product(shape(psit_f))*kind(psit_f)
              deallocate(psit_f, stat=istat)
              call memocc(istat, iall, 'psit_f', subname)
          end if
          allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, psit_c, 'psit_c', subname)
          allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, psit_f, 'psit_f', subname)

          call transpose_localized(iproc, nproc, orbs, collcom, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.

      end if
      call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)

      !DEBUG: alternative calculation of overlap - 
      !i=1
      !do iorb=1,orbs%norbp
      !   j=1
      !   do jorb=1,orbs%norbp
      !      ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
      !      jlr=orbs%inwhichlocreg(jorb+orbs%isorb)
      !      call wpdot_wrap(1,lzd%llr(ilr)%wfd%nvctr_c,lzd%llr(ilr)%wfd%nvctr_f,lzd%llr(ilr)%wfd%nseg_c,lzd%llr(ilr)%wfd%nseg_f,&
      !            lzd%llr(ilr)%wfd%keyvglob,lzd%llr(ilr)%wfd%keyglob,lphi(i),  &
      !            lzd%llr(jlr)%wfd%nvctr_c,lzd%llr(jlr)%wfd%nvctr_f,lzd%llr(jlr)%wfd%nseg_c,lzd%llr(jlr)%wfd%nseg_f,&
      !            lzd%llr(jlr)%wfd%keyvglob,lzd%llr(jlr)%wfd%keyglob,lphi(j),ovrlp(iorb,jorb))
      !      j=j+lzd%llr(jlr)%wfd%nvctr_c+7*lzd%llr(jlr)%wfd%nvctr_f
      !   end do
      !   i=i+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      !end do

      call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, methTransformOverlap, orthpar%blocksize_pdsyev, &
          orthpar%blocksize_pdgemm, orbs%norb, ovrlp)

      allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psittemp_c, 'psittemp_c', subname)
      allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psittemp_f, 'psittemp_f', subname)

      call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
      call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)
      call build_linear_combination_transposed(orbs%norb, ovrlp, collcom, psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
!      call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f)
      call untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, lphi, lzd)
      can_use_transposed=.true.

      ! alternative normalization - would need to switch back if keeping the transposed form for further use in eg calculating overlap
      i=1
      do iorb=1,orbs%norbp
         ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
         call normalizevector(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f,lphi(i))
         i=i+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      end do

      iall=-product(shape(psittemp_c))*kind(psittemp_c)
      deallocate(psittemp_c, stat=istat)
      call memocc(istat, iall, 'psittemp_c', subname)
      iall=-product(shape(psittemp_f))*kind(psittemp_f)
      deallocate(psittemp_f, stat=istat)
      call memocc(istat, iall, 'psittemp_f', subname)

  end do


  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)



end subroutine orthonormalizeLocalized




subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, orbs, op, comon, mad, collcom, orthpar, bpo, bs, &
           lphi, lhphi, lagmat, ovrlp, psit_c, psit_f, hpsit_c, hpsit_f, can_use_transposed, overlap_calculated)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_Data),intent(in) :: orbs
  type(overlapParameters),intent(inout) :: op
  type(p2pComms),intent(inout) :: comon
  type(matrixDescriptors),intent(in) :: mad
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  type(basis_performance_options),intent(in) :: bpo
  type(basis_specifications),intent(in):: bs
  real(kind=8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout) :: lphi,lhphi
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: lagmat, ovrlp
  real(8),dimension(:),pointer:: psit_c, psit_f, hpsit_c, hpsit_f
  logical,intent(inout):: can_use_transposed, overlap_calculated

  ! Local variables
  integer :: istat, iall, iorb, jorb
  real(kind=8),dimension(:),allocatable :: lphiovrlp
  real(kind=8),dimension(:,:),allocatable :: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans, lagmat_tmp
  character(len=*),parameter :: subname='orthoconstraintNonorthogonal'
  allocate(ovrlp_minus_one_lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
  allocate(ovrlp_minus_one_lagmat_trans(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
  allocate(lagmat_tmp(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, lagmat_tmp, 'lagmat_tmp', subname)

  if(.not. can_use_transposed) then
      allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psit_c, 'psit_c', subname)
      allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psit_f, 'psit_f', subname)
      call transpose_localized(iproc, nproc, orbs, collcom, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.
  end if
  ! It is assumed that this routine is called with the transposed gradient ready if it is associated...
  if(.not.associated(hpsit_c)) then
      allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
      allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
      call transpose_localized(iproc, nproc, orbs, collcom, lhphi, hpsit_c, hpsit_f, lzd)
  end if
  call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat)
  if(.not. overlap_calculated) then
      call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)
  end if
  overlap_calculated=.true.


  call dcopy(orbs%norb**2, lagmat(1,1), 1, lagmat_tmp(1,1), 1)
  call applyOrthoconstraintNonorthogonal2(iproc, nproc, orthpar%methTransformOverlap, orthpar%blocksize_pdgemm, &
       bs%correction_orthoconstraint, orbs, lagmat_tmp, ovrlp, mad, &
       ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)


  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          ovrlp(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat(jorb,iorb)
      end do
  end do
  call build_linear_combination_transposed(orbs%norb, ovrlp, collcom, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          ovrlp(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat_trans(jorb,iorb)
      end do
  end do
  call build_linear_combination_transposed(orbs%norb, ovrlp, collcom, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)

  call untranspose_localized(iproc, nproc, orbs, collcom, hpsit_c, hpsit_f, lhphi, lzd)


  iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
  deallocate(ovrlp_minus_one_lagmat, stat=istat)
  call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
  iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
  deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
  call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
  iall=-product(shape(lagmat_tmp))*kind(lagmat_tmp)
  deallocate(lagmat_tmp, stat=istat)
  call memocc(istat, iall, 'lagmat_tmp', subname)


end subroutine orthoconstraintNonorthogonal





subroutine initCommsOrtho(iproc, nproc, nspin, hx, hy, hz, lzd, lzdig, orbs, locregShape, bpo, op, comon) 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initCommsOrtho
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nspin
  real(kind=8),intent(in) :: hx, hy, hz
  type(local_zone_descriptors),intent(in) :: lzd, lzdig
  type(orbitals_data),intent(in) :: orbs
  character(len=1),intent(in) :: locregShape
  type(basis_performance_options),intent(in):: bpo
  type(overlapParameters),intent(out) :: op
  type(p2pComms),intent(out) :: comon

  ! Local variables
  integer :: iorb, jorb, iiorb
  integer ::  istat, jjorb, nsub, ierr
  character(len=*),parameter :: subname='initCommsOrtho'


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

  ! Count how many overlaping regions each orbital / process has.
  if(locregShape=='c') then
     stop "ERROR: locregShape=='c' is deprecated!"
  else if(locregShape=='s') then
     call determine_overlap_from_descriptors(iproc, nproc, orbs, orbs, lzd, lzd, op, comon)
  end if


  ! Initialize the communications.
  allocate(comon%comarr(6,maxval(comon%noverlaps),0:nproc-1), stat=istat)
  call memocc(istat, comon%comarr, 'comon%comarr', subname)

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

  call timing(iproc,'init_commOrtho','OF')

end subroutine initCommsOrtho






subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: mpisource, mpidest, istsource, istdest, ncount, tag
  integer,dimension(6),intent(out) :: comarr


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




subroutine applyOrthoconstraintNonorthogonal2(iproc, nproc, methTransformOverlap, blocksize_pdgemm, &
           correction_orthoconstraint, orbs, &
           lagmat, ovrlp, mad, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => applyOrthoconstraintNonorthogonal2
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, methTransformOverlap, blocksize_pdgemm, correction_orthoconstraint
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: ovrlp
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(inout) :: lagmat
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans

  ! Local variables
  integer :: iorb, jorb, istat, iall, ierr
  real(kind=8) :: tt, t1, t2, time_dsymm
  real(kind=8),dimension(:,:),allocatable :: ovrlp2
  character(len=*),parameter :: subname='applyOrthoconstraintNonorthogonal2'

  call timing(iproc,'lagmat_orthoco','ON')

  allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp2, 'ovrlp2', subname)

  correctionIf: if(correction_orthoconstraint==0) then
  
    call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)
  
    ! Invert the overlap matrix
    call timing(iproc,'lagmat_orthoco','OF')
    call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, blocksize_pdgemm, orbs%norb, mad, orbs, ovrlp2)
    call timing(iproc,'lagmat_orthoco','ON')
  
  
    ! Multiply the Lagrange multiplier matrix with S^-1/2.
    ! First fill the upper triangle.
    do iorb=1,orbs%norb
       do jorb=1,iorb-1
          ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
       end do
    end do
    if(blocksize_pdgemm<0) then
       t1=mpi_wtime()
       !!ovrlp_minus_one_lagmat=0.d0
       call to_zero(orbs%norb**2, ovrlp_minus_one_lagmat(1,1))
       call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
            0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
  
       ! Transpose lagmat
       do iorb=1,orbs%norb
           do jorb=iorb+1,orbs%norb
               tt=lagmat(jorb,iorb)
               lagmat(jorb,iorb)=lagmat(iorb,jorb)
               lagmat(iorb,jorb)=tt
           end do
       end do
       call to_zero(orbs%norb**2, ovrlp_minus_one_lagmat_trans(1,1))
       call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
            0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
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
  
  else if(correction_orthoconstraint==1) then correctionIf
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              ovrlp_minus_one_lagmat(jorb,iorb)=lagmat(jorb,iorb)
              ovrlp_minus_one_lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
          end do
      end do
  end if correctionIf
  
  call timing(iproc,'lagmat_orthoco','OF')
  
  
  iall=-product(shape(ovrlp2))*kind(ovrlp2)
  deallocate(ovrlp2, stat=istat)
  call memocc(istat, iall, 'ovrlp2', subname)


end subroutine applyOrthoconstraintNonorthogonal2


subroutine overlapPowerMinusOne(iproc, nproc, iorder, blocksize, norb, mad, orbs, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOne
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, norb
  type(orbitals_data),intent(in) :: orbs
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(norb,norb),intent(inout) :: ovrlp
  
  ! Local variables
  integer :: iorb, jorb, info
  character(len=*),parameter :: subname='overlapPowerMinusOne'

  call timing(iproc,'lovrlp^-1     ','ON')

  if(iorder==0) then

      ! Exact inversion
      if (blocksize<0) then
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
      else
          call dpotrf_parallel(iproc, nproc, blocksize, mpi_comm_world, 'l', norb, ovrlp(1,1), norb)
          call dpotri_parallel(iproc, nproc, blocksize, mpi_comm_world, 'l', norb, ovrlp(1,1), norb)
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

      stop 'overlapPowerMinusOne: iorder==2 is deprecated!'

  else

      write(*,'(1x,a)') 'ERROR: iorder must be 0,1 or 2!'

      stop
  end if

  call timing(iproc,'lovrlp^-1     ','OF')

end subroutine overlapPowerMinusOne




subroutine overlapPowerMinusOneHalf(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOneHalf
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb
  real(kind=8),dimension(norb,norb),intent(inout) :: ovrlp
  
  ! Local variables
  integer :: lwork, istat, iall, iorb, jorb, info
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
      
      
      if(blocksize_dsyev>0) then
          call dsyev_parallel(iproc, nproc, min(blocksize_dsyev,norb), comm, 'v', 'l', norb, ovrlp(1,1), norb, eval(1), info)
          if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
              !stop
          end if
      else
          
          !lwork=1000*norb
          allocate(work(1), stat=istat)
          call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, -1, info)
          lwork = work(1)
          deallocate(work, stat=istat)
          allocate(work(lwork), stat=istat)
          call memocc(istat, work, 'work', subname)
          call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, lwork, info)

          !*!lwork=1000*norb
          !*allocate(work(1), stat=istat)
          !*call DGEEV( 'v','v', norb, ovrlp(1,1), norb, eval, eval1, VL, norb, VR,&
          !*     norb, WORK, -1, info )
          !*lwork = work(1)
          !*deallocate(work, stat=istat)
          !*allocate(work(lwork), stat=istat)
          !*call memocc(istat, work, 'work', subname)
          !*! lr408 - see if LAPACK is stil to blame for convergence issues
          !*allocate(vl(1:norb,1:norb))
          !*allocate(vr(1:norb,1:norb))
          !*allocate(eval1(1:norb))
          !*call DGEEV( 'v','v', norb, ovrlp(1,1), norb, eval, eval1, VL, norb, VR,&
          !*     norb, WORK, LWORK, info )
          !*ovrlp=vl
          !*write(14,*) eval1
          !*deallocate(eval1)
          !*deallocate(vr)
          !*deallocate(vl)
          !*allocate(temp_vec(1:norb))
          !*do iorb=1,norb
          !*   do jorb=iorb+1,norb
          !*      if (eval(jorb) < eval(iorb)) then
          !*         temp = eval(iorb)
          !*         temp_vec = ovrlp(:,iorb)
          !*         eval(iorb) = eval(jorb)
          !*         eval(jorb) = temp
          !*         ovrlp(:,iorb) = ovrlp(:,jorb)
          !*         ovrlp(:,jorb) = temp_vec
          !*      end if
          !*   end do
          !*end do
          !*deallocate(temp_vec)

          !  lr408 - see if LAPACK is stil to blame for convergence issues
          if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev, info=', info
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
              tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
          end do
      end do
      
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

  else

      write(*,'(1x,a)') 'ERROR: methTransformOrder must be 0 or 1!'
      stop

  end if

  call timing(iproc,'lovrlp^-1/2   ','OF')

end subroutine overlapPowerMinusOneHalf


subroutine deviation_from_unity(iproc, norb, ovrlp, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: ovrlp
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, jorb
  real(8):: error

  deviation=0.d0
  do iorb=1,norb
     do jorb=1,norb
        if(iorb==jorb) then
           error=abs(ovrlp(jorb,iorb)-1.d0)
        else
           error=abs(ovrlp(jorb,iorb))
        end if
        if(error>deviation) then
           deviation=error
        end if
     end do
  end do

end subroutine deviation_from_unity

