!> @file
!! Orthonormalization
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
           orbs, lzd, sparsemat, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparseMatrix),intent(in) :: sparsemat
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(orbs%npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout):: can_use_transposed

  ! Local variables
  integer :: it, istat, iall
  !!integer :: ii, iiorb, jjorb, ilr, iorb, i, iseg, jlr, jorb, j
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm, ovrlp_compr 
  !!real(kind=8),dimension(:),allocatable :: ovrlp_compr2
  !real(kind=8),dimension(:,:),allocatable :: ovrlp
  character(len=*),parameter :: subname='orthonormalizeLocalized'
  !!real(kind=8) :: maxError, tt
  !real(kind=8) :: maxError
  !integer :: ilr,iorb,i,jlr,jorb,j


  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'
  !can_use_transposed=.false.
  do it=1,orthpar%nItOrtho

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
      allocate(ovrlp_compr(sparsemat%nvctr), stat=istat)
      call memocc(istat, ovrlp_compr, 'ovrlp_compr', subname)
      call calculate_overlap_transposed(iproc, nproc, orbs, sparsemat, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp_compr)
      !!allocate(ovrlp_compr2(sparsemat%nvctr))
      !!ovrlp_compr2=ovrlp_compr

      if (methTransformOverlap==-1) then
          !allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
          !call memocc(istat, ovrlp, 'ovrlp', subname)
          !call uncompressMatrix(orbs%norb, sparsemat, ovrlp_compr, ovrlp)
          call overlap_power_minus_one_half_per_atom(iproc, nproc, bigdft_mpi%mpi_comm, orbs, lzd, sparsemat, collcom, ovrlp_compr)
          !call compress_matrix_for_allreduce(orbs%norb, sparsemat, ovrlp, ovrlp_compr)
          !iall=-product(shape(ovrlp))*kind(ovrlp)
          !deallocate(ovrlp, stat=istat)
          !call memocc(istat, iall, 'ovrlp', subname)
      else
          call overlapPowerMinusOneHalf(iproc, nproc, bigdft_mpi%mpi_comm, methTransformOverlap, orthpar%blocksize_pdsyev, &
              orthpar%blocksize_pdgemm, orbs%norb, orbs%norbp, orbs%isorb, sparsemat, ovrlp_compr)
          !!call uncompressMatrix(orbs%norb, sparsemat, ovrlp_compr, ovrlp)
      end if

      !!if (iproc==0) then
      !!    do istat=1,sparsemat%nvctr
      !!        write(333,'(i9,2es20.10)') istat, ovrlp_compr(istat), ovrlp_compr2(istat)
      !!    end do
      !!end if
      !!ovrlp_compr2=abs(ovrlp_compr2-ovrlp_compr)
      !!if (iproc==0) write(*,'(a,es20.10,i9)') 'maxval(ovrlp_compr2), maxloc(ovrlp_compr2)',maxval(ovrlp_compr2), maxloc(ovrlp_compr2)


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

      allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psittemp_c, 'psittemp_c', subname)
      allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psittemp_f, 'psittemp_f', subname)

      call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
      call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)
      call build_linear_combination_transposed(orbs%norb, ovrlp_compr, collcom, sparsemat, &
           psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
      allocate(norm(orbs%norb), stat=istat)
      call memocc(istat, norm, 'norm', subname)
      call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)
      !call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)
      iall=-product(shape(norm))*kind(norm)
      deallocate(norm, stat=istat)
      call memocc(istat, iall, 'norm', subname)
      call untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, lphi, lzd)

      !!! TEST ##################################
      !!call calculate_overlap_transposed(iproc, nproc, orbs, sparsemat, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp_compr)
      !!!!do iorb=1,orbs%norb
      !!!!    do jorb=1,orbs%norb
      !!!!        write(10000+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
      !!!!    end do
      !!!!end do
      !!if (iproc==0)  then
      !!    ii=0
      !!    maxError=0.d0
      !!    do iseg=1,sparsemat%nseg
      !!        do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
      !!            ii=ii+1
      !!            iiorb = (jorb-1)/orbs%norb + 1
      !!            jjorb = jorb - (iiorb-1)*orbs%norb
      !!            if (jjorb==iiorb) then
      !!                tt=abs(1.d0-ovrlp_compr(ii))
      !!            else
      !!                tt=abs(ovrlp_compr(ii))
      !!            end if
      !!            if (tt>maxError) then
      !!                maxError=tt
      !!            end if
      !!            !write(999,*) iiorb, jjorb, ovrlp_compr(ii)
      !!        end do
      !!    end do
      !!    write(888,*) 'maxError', maxError
      !!end if
      !!! END TEST ##############################

      iall=-product(shape(ovrlp_compr))*kind(ovrlp_compr)
      deallocate(ovrlp_compr, stat=istat)
      call memocc(istat, iall, 'ovrlp_compr', subname)

      ! alternative normalization - would need to switch back if keeping the transposed form for further use in eg calculating overlap
      !i=1
      !do iorb=1,orbs%norbp
      !   ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
      !   call normalizevector(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f,lphi(i))
      !   i=i+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      !end do

      iall=-product(shape(psittemp_c))*kind(psittemp_c)
      deallocate(psittemp_c, stat=istat)
      call memocc(istat, iall, 'psittemp_c', subname)
      iall=-product(shape(psittemp_f))*kind(psittemp_f)
      deallocate(psittemp_f, stat=istat)
      call memocc(istat, iall, 'psittemp_f', subname)

  end do


end subroutine orthonormalizeLocalized


! can still tidy this up more when tmblarge is removed
! use sparsity of density kernel for all inverse quantities
subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, orbs, collcom, orthpar, correction_orthoconstraint, &
           linmat, lphi, lhphi, lagmat, psit_c, psit_f, hpsit_c, hpsit_f, can_use_transposed, overlap_calculated)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_Data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  integer,intent(in) :: correction_orthoconstraint
  real(kind=8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(in) :: lphi
  real(kind=8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout) :: lhphi
  type(SparseMatrix),intent(inout) :: lagmat
  real(kind=8),dimension(:),pointer :: psit_c, psit_f, hpsit_c, hpsit_f
  logical,intent(inout) :: can_use_transposed, overlap_calculated
  type(linear_matrices),intent(inout) :: linmat

  ! Local variables
  integer :: istat, iall, iorb, jorb, ii, ii_trans
  real(kind=8),dimension(:),pointer :: tmp_mat_compr
  type(SparseMatrix) :: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
  character(len=*),parameter :: subname='orthoconstraintNonorthogonal'


  ! ASSUME denskern sparsity pattern is symmetric
  ! create ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans with sparsity pattern of denskern here
  ! this is slight overkill for no orthoconstraint correction, think about going back to just matrices
  call nullify_sparsematrix(ovrlp_minus_one_lagmat)
  call sparse_copy_pattern(linmat%denskern, ovrlp_minus_one_lagmat, subname)
  call nullify_sparsematrix(ovrlp_minus_one_lagmat_trans)
  call sparse_copy_pattern(linmat%denskern, ovrlp_minus_one_lagmat_trans, subname)

  if (correction_orthoconstraint==0) then
      allocate(ovrlp_minus_one_lagmat%matrix_compr(ovrlp_minus_one_lagmat%nvctr), stat=istat)
      call memocc(istat, ovrlp_minus_one_lagmat%matrix_compr, 'ovrlp_minus_one_lagmat%matrix_compr', subname)
      allocate(ovrlp_minus_one_lagmat_trans%matrix_compr(ovrlp_minus_one_lagmat_trans%nvctr), stat=istat)
      call memocc(istat, ovrlp_minus_one_lagmat_trans%matrix_compr, 'ovrlp_minus_one_lagmat_trans%matrix_compr', subname)
  else
      ovrlp_minus_one_lagmat%matrix_compr => lagmat%matrix_compr
      ovrlp_minus_one_lagmat_trans%matrix_compr => lagmat%matrix_compr
  end if

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

  call calculate_overlap_transposed(iproc, nproc, orbs, lagmat, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat%matrix_compr)

  if (correction_orthoconstraint==0) then
      if(overlap_calculated) stop 'overlap_calculated should be wrong... To be modified later'

      call calculate_overlap_transposed(iproc, nproc, orbs, linmat%ovrlp, collcom, psit_c, psit_c, psit_f, psit_f, &
           linmat%ovrlp%matrix_compr)

      allocate(linmat%ovrlp%matrix(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, linmat%ovrlp%matrix, 'linmat%ovrlp%matrix', subname)

      call uncompressMatrix(orbs%norb, linmat%ovrlp, linmat%ovrlp%matrix_compr, linmat%ovrlp%matrix)

      allocate(lagmat%matrix(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, lagmat%matrix, 'lagmat%matrix', subname)

      call uncompressMatrix(orbs%norb, lagmat, lagmat%matrix_compr, lagmat%matrix)

      allocate(ovrlp_minus_one_lagmat%matrix(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, ovrlp_minus_one_lagmat%matrix, 'ovrlp_minus_one_lagmat%matrix', subname)
 
      allocate(ovrlp_minus_one_lagmat_trans%matrix(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, ovrlp_minus_one_lagmat_trans%matrix, 'ovrlp_minus_one_lagmat_trans%matrix', subname)

      call applyOrthoconstraintNonorthogonal2(iproc, nproc, orthpar%methTransformOverlap, orthpar%blocksize_pdgemm, &
           correction_orthoconstraint, orbs, lagmat%matrix, linmat%ovrlp%matrix, &
           ovrlp_minus_one_lagmat%matrix, ovrlp_minus_one_lagmat_trans%matrix)

      iall=-product(shape(linmat%ovrlp%matrix))*kind(linmat%ovrlp%matrix)
      deallocate(linmat%ovrlp%matrix, stat=istat)
      call memocc(istat, iall, 'linmat%ovrlp%matrix', subname)

      iall=-product(shape(lagmat%matrix))*kind(lagmat%matrix)
      deallocate(lagmat%matrix, stat=istat)
      call memocc(istat, iall, 'lagmat%matrix', subname)

      call compress_matrix_for_allreduce(orbs%norb, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat%matrix, &
           ovrlp_minus_one_lagmat%matrix_compr)

      iall=-product(shape(ovrlp_minus_one_lagmat%matrix))*kind(ovrlp_minus_one_lagmat%matrix)
      deallocate(ovrlp_minus_one_lagmat%matrix, stat=istat)
      call memocc(istat, iall, 'ovrlp_minus_one_lagmat%matrix', subname)
 
      call compress_matrix_for_allreduce(orbs%norb, ovrlp_minus_one_lagmat_trans, ovrlp_minus_one_lagmat_trans%matrix, &
           ovrlp_minus_one_lagmat_trans%matrix_compr)

      iall=-product(shape(ovrlp_minus_one_lagmat_trans%matrix))*kind(ovrlp_minus_one_lagmat_trans%matrix)
      deallocate(ovrlp_minus_one_lagmat_trans%matrix, stat=istat)
      call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans%matrix', subname)
  end if

  allocate(tmp_mat_compr(ovrlp_minus_one_lagmat%nvctr), stat=istat)
  call memocc(istat, tmp_mat_compr, 'tmp_mat_compr', subname)

  do jorb=1,orbs%norb
     do iorb=1,orbs%norb
          ii_trans = ovrlp_minus_one_lagmat_trans%matrixindex_in_compressed(jorb, iorb)
          ii = ovrlp_minus_one_lagmat%matrixindex_in_compressed(iorb, jorb)
          if (ii==0.or.ii_trans==0) cycle
          tmp_mat_compr(ii)=-0.5d0*ovrlp_minus_one_lagmat%matrix_compr(ii) &
               -0.5d0*ovrlp_minus_one_lagmat_trans%matrix_compr(ii_trans)
      end do
  end do

  call build_linear_combination_transposed(orbs%norb, tmp_mat_compr, collcom, &
       ovrlp_minus_one_lagmat, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)

  iall=-product(shape(tmp_mat_compr))*kind(tmp_mat_compr)
  deallocate(tmp_mat_compr, stat=istat)
  call memocc(istat, iall, 'tmp_mat_compr', subname)

  call untranspose_localized(iproc, nproc, orbs, collcom, hpsit_c, hpsit_f, lhphi, lzd)

  if (correction_orthoconstraint/=0) then
      nullify(ovrlp_minus_one_lagmat%matrix_compr)
      nullify(ovrlp_minus_one_lagmat_trans%matrix_compr)
  end if

  call deallocate_sparseMatrix(ovrlp_minus_one_lagmat, subname)
  call deallocate_sparseMatrix(ovrlp_minus_one_lagmat_trans, subname)

  overlap_calculated=.false.

end subroutine orthoconstraintNonorthogonal




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
  
    call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)

    allocate(lagmat_trans(orbs%norb,orbs%norb), stat=istat)
    call memocc(istat, lagmat_trans, 'lagmat_trans', subname)

    call dcopy(orbs%norb**2, lagmat(1,1), 1, lagmat_trans(1,1), 1)
  
    ! Invert the overlap matrix
    call timing(iproc,'lagmat_orthoco','OF')
    call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, blocksize_pdgemm, orbs%norb, ovrlp2)
    call timing(iproc,'lagmat_orthoco','ON')
  
  
    ! Multiply the Lagrange multiplier matrix with S^-1/2.
    ! First fill the upper triangle.
    do iorb=1,orbs%norb
       do jorb=1,iorb-1
          ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
       end do
    end do
    if(blocksize_pdgemm<0) then
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
    else
      call dsymm_parallel(iproc, nproc, blocksize_pdgemm, bigdft_mpi%mpi_comm, 'l', 'l', orbs%norb, orbs%norb, 1.d0, &
           ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, 0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
      ! Transpose lagmat
      do iorb=1,orbs%norb
          do jorb=iorb+1,orbs%norb
              lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
              lagmat_trans(iorb,jorb)=lagmat(jorb,iorb)
          end do
      end do
      call dsymm_parallel(iproc, nproc, blocksize_pdgemm, bigdft_mpi%mpi_comm, 'l', 'l', orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), &
           orbs%norb, lagmat_trans(1,1), orbs%norb, &
           0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)
    end if

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


subroutine overlapPowerMinusOne(iproc, nproc, iorder, blocksize, norb, ovrlp)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOne
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, norb
  real(kind=8),dimension(norb,norb),intent(inout) :: ovrlp
  
  ! Local variables
  integer :: iorb, jorb, info

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
          call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, ovrlp(1,1), norb)
          call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, ovrlp(1,1), norb)
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




subroutine overlapPowerMinusOneHalf(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, &
           blocksize_pdgemm, norb, norbp, isorb, sparsemat, ovrlp_compr)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOneHalf
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb, norbp, isorb
  type(sparseMatrix),intent(in) :: sparsemat
  real(kind=8),dimension(sparsemat%nvctr),intent(inout) :: ovrlp_compr

  
  ! Local variables
  integer :: lwork, istat, iall, iorb, jorb, info, iseg, iiorb, jjorb, ii
  character(len=*),parameter :: subname='overlapPowerMinusOneHalf'
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:,:),allocatable :: tempArr
  real(kind=8),dimension(:,:),allocatable :: ovrlp
  !!real(kind=8),dimension(:,:),allocatable :: temp1, temp2
  !!real(kind=8),dimension(5) :: cc
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

      allocate(ovrlp(norb,norb), stat=istat)
      call memocc(istat, ovrlp, 'ovrlp', subname)
      call uncompressMatrix(norb, sparsemat, ovrlp_compr, ovrlp)
      
      
      if(blocksize_dsyev>0) then
          call dsyev_parallel(iproc, nproc, min(blocksize_dsyev,norb), comm, 'v', 'l', norb, ovrlp(1,1), norb, eval(1), info)
          if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
              !stop
          end if
      else
          
          !lwork=1000*norb
          allocate(work(100), stat=istat)
          call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, -1, info)
          lwork = int(work(1))
          deallocate(work, stat=istat)
          allocate(work(lwork), stat=istat)
          call memocc(istat, work, 'work', subname)
          call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, lwork, info)

          !*!lwork=1000*norb
          !*allocate(work(1), stat=istat)
          !*call DGEEV( 'v','v', norb, ovrlp(1,1), norb, eval, eval1, VL, norb, VR,&
          !*     norb, WORK, -1, info )
          !*lwork = nint(work(1))
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

      call compress_matrix_for_allreduce(norb, sparsemat, ovrlp, ovrlp_compr)

      iall=-product(shape(eval))*kind(eval)
      deallocate(eval, stat=istat)
      call memocc(istat, iall, 'eval', subname)
      iall=-product(shape(tempArr))*kind(tempArr)
      deallocate(tempArr, stat=istat)
      call memocc(istat, iall, 'tempArr', subname)

      iall=-product(shape(ovrlp))*kind(ovrlp)
      deallocate(ovrlp, stat=istat)
      call memocc(istat, iall, 'ovrlp', subname)

  else if(methTransformOrder==1) then


      ! Taylor expansion up to first order.
      !!!$omp parallel do default(private) shared(norb,sparsemat,ovrlp)
      ii=0
      do iseg=1,sparsemat%nseg
          do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/norb + 1
              jjorb = jorb - (iiorb-1)*norb
              if(iiorb==jjorb) then
                  ovrlp_compr(ii)=1.5d0-.5d0*ovrlp_compr(ii)
              else
                  ovrlp_compr(ii)=-.5d0*ovrlp_compr(ii)
              end if
          end do
      end do
      !!!$omp end parallel do

      !!if (present(sparsemat)) then

      !!    ! Matrix compression is availabale
 
      !!    !$omp parallel do default(private) shared(norb,sparsemat,ovrlp)
      !!    do iseg=1,sparsemat%nseg
      !!        do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
      !!            iiorb = (jorb-1)/norb + 1
      !!            jjorb = jorb - (iiorb-1)*norb
      !!            if(iiorb==jjorb) then
      !!                ovrlp(jjorb,iiorb)=1.5d0-.5d0*ovrlp(jjorb,iiorb)
      !!            else
      !!                ovrlp(jjorb,iiorb)=-.5d0*ovrlp(jjorb,iiorb)
      !!            end if
      !!        end do
      !!    end do
      !!   !$omp end parallel do
 
      !!else

      !!    ! No matrix compression available
      !!    !$omp parallel do default(private) shared(ovrlp,norb)
      !!    do iorb=1,norb
      !!        do jorb=1,norb
      !!            if(iorb==jorb) then
      !!                ovrlp(jorb,iorb)=1.5d0-.5d0*ovrlp(jorb,iorb)
      !!            else
      !!                ovrlp(jorb,iorb)=-.5d0*ovrlp(jorb,iorb)
      !!            end if
      !!        end do
      !!    end do
      !!    !$omp end parallel do

      !!end if


  else
      
      stop 'deprecated'

      !!! Taylor expansion up to 4th order
      !!if (norbp>0) then
      !!    cc(1)=-1.d0/2.d0
      !!    cc(2)=3.d0/8.d0
      !!    cc(3)=-15.d0/48.d0
      !!    cc(4)=105.d0/384.d0
      !!    cc(5)=-945.d0/3840.d0


      !!    do iorb=1,norb
      !!        do jorb=1,norb
      !!            !!write(400+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
      !!            if (iorb==jorb) then
      !!                ovrlp(jorb,iorb)=ovrlp(jorb,iorb)-1.d0
      !!            else
      !!                ovrlp(jorb,iorb)=ovrlp(jorb,iorb)
      !!            end if
      !!            !!write(420+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
      !!        end do
      !!    end do
      !!    allocate(ovrlp_compr(sparsemat%nvctr), stat=istat)
      !!    call memocc(istat, ovrlp_compr, 'ovrlp_compr', subname)
      !!    call compress_matrix_for_allreduce(norb, sparsemat, ovrlp, ovrlp_compr)

      !!    call to_zero(norb*norb, ovrlp(1,1))
      !!    do iorb=isorb+1,isorb+norbp
      !!        ovrlp(iorb,iorb)=1.d0
      !!    end do

      !!    allocate(temp1(norb,norbp), stat=istat)
      !!    call memocc(istat, temp1, 'temp1', subname)
      !!    allocate(temp2(norb,norbp), stat=istat)
      !!    call memocc(istat, temp2, 'temp2', subname)
      !!    call to_zero(norb*norbp, temp2(1,1))

      !!    do iorb=1,norbp
      !!        iiorb=isorb+iorb
      !!        do jorb=1,norb
      !!            if (iiorb==jorb) then
      !!                temp1(jorb,iorb)=1.d0
      !!            else
      !!                temp1(jorb,iorb)=0.d0
      !!            end if
      !!        end do
      !!    end do


      !!    !!do iorb=1,norb
      !!    !!    do jorb=1,norb
      !!    !!        write(460+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
      !!    !!    end do
      !!    !!end do
      !!    do i=1,5
      !!        call sparsemm(ovrlp_compr, temp1, temp2, norb, norbp, sparsemat)
      !!        do iorb=1,norbp
      !!            !!do jorb=1,norb
      !!            !!    write(440+iproc,*) iorb,jorb,temp2(jorb,iorb)
      !!            !!end do
      !!        end do
      !!        call daxpy(norb*norbp, cc(i), temp2, 1, ovrlp(1,isorb+1), 1) 
      !!        call dcopy(norb*norbp, temp2, 1, temp1, 1)
      !!    end do
      !!    !!do iorb=1,norb
      !!    !!    do jorb=1,norb
      !!    !!        write(480+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
      !!    !!    end do
      !!    !!end do

      !!    iall=-product(shape(temp1))*kind(temp1)
      !!    deallocate(temp1, stat=istat)
      !!    call memocc(istat, iall, 'temp1', subname)
      !!    iall=-product(shape(temp2))*kind(temp2)
      !!    deallocate(temp2, stat=istat)
      !!    call memocc(istat, iall, 'temp2', subname)

      !!end if

      !!call mpiallred(ovrlp(1,1), norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)



      !!write(*,'(1x,a)') 'ERROR: methTransformOrder must be 0 or 1!'
      !!stop

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


subroutine overlap_power_minus_one_half_per_atom(iproc, nproc, comm, orbs, lzd, sparsemat, collcom, ovrlp_compr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => overlap_power_minus_one_half_per_atom
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, comm
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparseMatrix),intent(in) :: sparsemat
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(sparsemat%nvctr),intent(inout) :: ovrlp_compr

  ! Local variables
  integer :: iend, i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb, ilr
  integer :: iiorb, ierr, ii, iseg, ind
  !!integer :: llorb, j, jlr, iiorbold, istart, ia, iaold
  !!real(kind=8) :: tt
  real(kind=8),dimension(:,:),allocatable :: ovrlp_tmp
  !!real(kind=8),dimension(:,:),allocatable :: ovrlp_old
  real(kind=8),dimension(:),allocatable :: ovrlp_compr_old
  logical,dimension(:),allocatable :: in_neighborhood
  character(len=*),parameter :: subname='overlap_power_minus_one_half_per_atom'


      !!! after input guess: only orthonormalize the orbitals on the same atom
      !!istart=1
      !!iaold=orbs%onwhichatom(istart)
      !!do iorb=1,orbs%norb
      !!    ia=orbs%onwhichatom(iorb)
      !!    if (ia/=iaold) then
      !!        ! We are at the start of a new atom
      !!        ! End of previous atom
      !!        iend=iorb-1
      !!        
      !!        ! Calculate S^-1/2 for the small overlap matrix of this atom
      !!        n=iend-istart+1
      !!        allocate(ovrlp_tmp(n,n), stat=istat)
      !!        call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
      !!        do i=1,n
      !!            do j=1,n
      !!                ovrlp_tmp(j,i)=ovrlp(istart+j-1,istart+i-1)
      !!            end do
      !!        end do
      !!        call overlapPowerMinusOneHalf(iproc, nproc, comm, 0, -8, -8, n, 0, 0, ovrlp_tmp)

      !!        ! Fill it back to the original matrix, filling the remaining parts of the column with zeros.
      !!        call to_zero(n*orbs%norb, ovrlp(1,istart))
      !!        do i=1,n
      !!            do j=1,n
      !!                ovrlp(istart+j-1,istart+i-1)=ovrlp_tmp(j,i)
      !!            end do
      !!        end do

      !!        iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
      !!        deallocate(ovrlp_tmp, stat=istat)
      !!        call memocc(istat, iall, 'ovrlp_tmp', subname)
      !!        
      !!        istart=iorb
      !!        iaold=ia
      !!    end if
      !!end do

      allocate(in_neighborhood(orbs%norb), stat=istat)
      call memocc(istat, in_neighborhood, 'in_neighborhood', subname)
      !!allocate(ovrlp_old(orbs%norb,orbs%norb), stat=istat)
      !!call memocc(istat, ovrlp_old, 'ovrlp_old', subname)
      !!call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp_old(1,1), 1)
      !!call to_zero(orbs%norb**2, ovrlp(1,1))

      allocate(ovrlp_compr_old(sparsemat%nvctr), stat=istat)
      call memocc(istat, ovrlp_compr_old, 'ovrlp_compr_old', subname)
      call vcopy(sparsemat%nvctr, ovrlp_compr(1), 1, ovrlp_compr_old(1), 1)

      call to_zero(sparsemat%nvctr, ovrlp_compr(1))

      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          ! We are at the start of a new atom
          ! Count all orbitals that are in the neighborhood

          iseg=sparsemat%istsegline(iiorb)
          iend =iiorb*orbs%norb
          n=0
          in_neighborhood(:)=.false.
          do 
              do i=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
                  ii=i-(iiorb-1)*orbs%norb
                  in_neighborhood(ii)=.true.
                  n=n+1
              end do
              iseg=iseg+1
              if (iseg>sparsemat%nseg) exit
              if (sparsemat%keyg(1,iseg)>iend) exit
          end do

          !!n=0
          !!do jorb=1,orbs%norb
          !!    jlr=orbs%inwhichlocreg(jorb)
          !!    tt = (lzd%llr(ilr)%locregcenter(1)-lzd%llr(jlr)%locregcenter(1))**2 &
          !!         + (lzd%llr(ilr)%locregcenter(2)-lzd%llr(jlr)%locregcenter(2))**2 &
          !!         + (lzd%llr(ilr)%locregcenter(3)-lzd%llr(jlr)%locregcenter(3))**2
          !!    tt=sqrt(tt)
          !!    if (tt<10.d0) then
          !!        n=n+1
          !!        in_neighborhood(jorb)=.true.
          !!    else
          !!        in_neighborhood(jorb)=.false.
          !!    end if
          !!end do
          allocate(ovrlp_tmp(n,n), stat=istat)
          call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
          call to_zero(n*n, ovrlp_tmp(1,1))

          !!jjorb=0
          !!do jorb=1,orbs%norb
          !!    if (.not.in_neighborhood(jorb)) cycle
          !!    jjorb=jjorb+1
          !!    kkorb=0
          !!    do korb=1,orbs%norb
          !!        if (.not.in_neighborhood(korb)) cycle
          !!        kkorb=kkorb+1
          !!        ovrlp_tmp(kkorb,jjorb)=ovrlp_old(korb,jorb)
          !!        write(1100+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
          !!    end do
          !!end do

          jjorb=0
          do jorb=1,orbs%norb
              if (.not.in_neighborhood(jorb)) cycle
              jjorb=jjorb+1
              kkorb=0
              do korb=1,orbs%norb
                  if (.not.in_neighborhood(korb)) cycle
                  kkorb=kkorb+1
                  ind = compressed_index(korb, jorb, orbs%norb, sparsemat)
                  if (ind>0) then
                      ovrlp_tmp(kkorb,jjorb)=ovrlp_compr_old(ind)
                  else
                      ovrlp_tmp(kkorb,jjorb)=0.d0
                  end if
                  !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
              end do
          end do

          !!iiorbold=0
          !!llorb=0
          !!ii=0
          !!do iseg=1,sparsemat%nseg
          !!    do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
          !!        ii=ii+1
          !!        iiorb = (jorb-1)/orbs%norb + 1
          !!        jjorb = jorb - (iiorb-1)*orbs%norb
          !!        if (.not.in_neighborhood(iiorb)) cycle
          !!        if (.not.in_neighborhood(jjorb)) cycle
          !!        if (iiorb/=iiorbold) then
          !!            llorb=llorb+1
          !!            kkorb=0
          !!        end if
          !!        kkorb=kkorb+1
          !!        iiorbold=iiorb
          !!        !if (iproc==0) write(*,'(2(a,i0),a,es14.6)') 'fill entry ',kkorb,',',llorb,' with ', ovrlp_compr(ii)
          !!        ovrlp_tmp(kkorb,llorb)=ovrlp_compr(ii)
          !!        write(1000+iproc,'(4i8,2es20.10)') iiorb, jjorb, kkorb, llorb, ovrlp_tmp(kkorb,llorb), ovrlp_old(jjorb,iiorb)
          !!    end do
          !!end do
          
          ! Calculate S^-1/2 for the small overlap matrix
          call overlapPowerMinusOneHalf_old(iproc, nproc, comm, 0, -8, -8, n, 0, 0, ovrlp_tmp)

          ! Fill it back to the original matrix, filling the remaining parts of the column with zeros.
          !!jjorb=0
          !!do jorb=1,orbs%norb
          !!    if (.not.in_neighborhood(jorb)) cycle
          !!    jjorb=jjorb+1
          !!    kkorb=0
          !!    if (jorb==iiorb) then
          !!        do korb=1,orbs%norb
          !!            if (.not.in_neighborhood(korb)) cycle
          !!            kkorb=kkorb+1
          !!            ovrlp(korb,jorb)=ovrlp_tmp(kkorb,jjorb)
          !!        end do
          !!    end if
          !!end do


          jjorb=0
          do jorb=1,orbs%norb
              if (.not.in_neighborhood(jorb)) cycle
              jjorb=jjorb+1
              kkorb=0
              if (jorb==iiorb) then
                  do korb=1,orbs%norb
                      if (.not.in_neighborhood(korb)) cycle
                      kkorb=kkorb+1
                      !ind = compressed_index(korb, jorb, orbs%norb, sparsemat)
                      ind = compressed_index(jorb, korb, orbs%norb, sparsemat)
                      !ind = collcom%matrixindex_in_compressed(korb,jorb)
                      !ind = collcom%matrixindex_in_compressed(jorb,korb)
                      if (ind>0) then
                          ovrlp_compr(ind)=ovrlp_tmp(kkorb,jjorb)
                      end if
                      !write(1300+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
                  end do
              end if
          end do



          !!iiorbold=0
          !!llorb=0
          !!ii=0
          !!do iseg=1,sparsemat%nseg
          !!    do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
          !!        ii=ii+1
          !!        iiorb = (jorb-1)/orbs%norb + 1
          !!        jjorb = jorb - (iiorb-1)*orbs%norb
          !!        if (.not.in_neighborhood(iiorb)) cycle
          !!        if (.not.in_neighborhood(jjorb)) cycle
          !!        if (iiorb/=iiorbold) then
          !!            llorb=llorb+1
          !!            kkorb=0
          !!        end if
          !!        kkorb=kkorb+1
          !!        iiorbold=iiorb
          !!        ovrlp_tmp(kkorb,llorb)=ovrlp_compr(ii)
          !!    end do
          !!end do

          iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
          deallocate(ovrlp_tmp, stat=istat)
          call memocc(istat, iall, 'ovrlp_tmp', subname)
      end do

      iall=-product(shape(ovrlp_compr_old))*kind(ovrlp_compr_old)
      deallocate(ovrlp_compr_old, stat=istat)
      call memocc(istat, iall, 'ovrlp_compr_old', subname)

      !call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      call mpiallred(ovrlp_compr(1), sparsemat%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)

      iall=-product(shape(in_neighborhood))*kind(in_neighborhood)
      deallocate(in_neighborhood, stat=istat)
      call memocc(istat, iall, 'in_neighborhood', subname)
      !!iall=-product(shape(ovrlp_old))*kind(ovrlp_old)
      !!deallocate(ovrlp_old, stat=istat)
      !!call memocc(istat, iall, 'ovrlp_old', subname)

 contains

    function compressed_index(iiorb, jjorb, norb, sparsemat)
      use module_base
      use module_types
      implicit none

      ! Calling arguments
      integer,intent(in) :: iiorb, jjorb, norb
      type(sparseMatrix),intent(in) :: sparsemat
      integer :: compressed_index
      logical :: notfound

      ! Local variables
      integer :: ii, iseg, isegend

      ii=(iiorb-1)*norb+jjorb

      iseg=sparsemat%istsegline(iiorb)
      ! isegend is the last possible segments where the index might be
      if(iiorb<norb) then
          isegend=sparsemat%istsegline(iiorb+1)-1
      else
          isegend=sparsemat%nseg
      end if

      notfound=.false.
      do
      !write(*,'(a,6i9)') 'iiorb, jjorb, ii, iseg, sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)', iiorb, jjorb, ii, iseg, sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
          if (ii>=sparsemat%keyg(1,iseg) .and. ii<=sparsemat%keyg(2,iseg)) then
              ! The matrix element is in this segment
              exit
          end if
          iseg=iseg+1
          if (iseg>isegend) then
              notfound=.true.
              exit
          end if
      end do

      if (notfound) then
          compressed_index=-1
      else
          compressed_index = sparsemat%keyv(iseg) + ii - sparsemat%keyg(1,iseg)
      end if

    end function compressed_index


end subroutine overlap_power_minus_one_half_per_atom



! Should be used if sparsemat is not available... to be cleaned
subroutine overlapPowerMinusOneHalf_old(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, &
           blocksize_pdgemm, norb, norbp, isorb, ovrlp, sparsemat)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => overlapPowerMinusOneHalf_old
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb, norbp, isorb
  real(kind=8),dimension(norb,norb),intent(inout) :: ovrlp
  type(sparseMatrix),intent(in),optional :: sparsemat

  
  ! Local variables
  integer :: lwork, istat, iall, iorb, jorb, info, iseg, iiorb, jjorb
  !!integer :: i, ierr
  character(len=*),parameter :: subname='overlapPowerMinusOneHalf'
  real(kind=8),dimension(:),allocatable :: eval, work
  !!real(kind=8),dimension(:),allocatable :: ovrlp_compr
  real(kind=8),dimension(:,:,:),allocatable :: tempArr
  !!real(kind=8),dimension(:,:),allocatable :: temp1, temp2
  !!real(kind=8),dimension(5) :: cc
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
          allocate(work(10), stat=istat)
          call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, -1, info)

          lwork = nint(work(1))

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
              tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)/sqrt(abs(eval(iorb)))
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
      if (present(sparsemat)) then

          ! Matrix compression is availabale
 
          !$omp parallel do default(private) shared(norb,sparsemat,ovrlp)
          do iseg=1,sparsemat%nseg
              do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
                  iiorb = (jorb-1)/norb + 1
                  jjorb = jorb - (iiorb-1)*norb
                  if(iiorb==jjorb) then
                      ovrlp(jjorb,iiorb)=1.5d0-.5d0*ovrlp(jjorb,iiorb)
                  else
                      ovrlp(jjorb,iiorb)=-.5d0*ovrlp(jjorb,iiorb)
                  end if
              end do
          end do
         !$omp end parallel do
 
      else

          ! No matrix compression available
          !$omp parallel do default(private) shared(ovrlp,norb)
          do iorb=1,norb
              do jorb=1,norb
                  if(iorb==jorb) then
                      ovrlp(jorb,iorb)=1.5d0-.5d0*ovrlp(jorb,iorb)
                  else
                      ovrlp(jorb,iorb)=-.5d0*ovrlp(jorb,iorb)
                  end if
              end do
          end do
          !$omp end parallel do

      end if


  else

      stop 'must fix this'

!!      ! Taylor expansion up to 4th order
!!      if (norbp>0) then
!!          cc(1)=-1.d0/2.d0
!!          cc(2)=3.d0/8.d0
!!          cc(3)=-15.d0/48.d0
!!          cc(4)=105.d0/384.d0
!!          cc(5)=-945.d0/3840.d0
!!
!!
!!          do iorb=1,norb
!!              do jorb=1,norb
!!                  !!write(400+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
!!                  if (iorb==jorb) then
!!                      ovrlp(jorb,iorb)=ovrlp(jorb,iorb)-1.d0
!!                  else
!!                      ovrlp(jorb,iorb)=ovrlp(jorb,iorb)
!!                  end if
!!                  !!write(420+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
!!              end do
!!          end do
!!          allocate(ovrlp_compr(sparsemat%nvctr), stat=istat)
!!          call memocc(istat, ovrlp_compr, 'ovrlp_compr', subname)
!!          call compress_matrix_for_allreduce(norb, sparsemat, ovrlp, ovrlp_compr)
!!
!!          call to_zero(norb*norb, ovrlp(1,1))
!!          do iorb=isorb+1,isorb+norbp
!!              ovrlp(iorb,iorb)=1.d0
!!          end do
!!
!!          allocate(temp1(norb,norbp), stat=istat)
!!          call memocc(istat, temp1, 'temp1', subname)
!!          allocate(temp2(norb,norbp), stat=istat)
!!          call memocc(istat, temp2, 'temp2', subname)
!!          call to_zero(norb*norbp, temp2(1,1))
!!
!!          do iorb=1,norbp
!!              iiorb=isorb+iorb
!!              do jorb=1,norb
!!                  if (iiorb==jorb) then
!!                      temp1(jorb,iorb)=1.d0
!!                  else
!!                      temp1(jorb,iorb)=0.d0
!!                  end if
!!              end do
!!          end do
!!
!!
!!          !!do iorb=1,norb
!!          !!    do jorb=1,norb
!!          !!        write(460+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
!!          !!    end do
!!          !!end do
!!          do i=1,5
!!              call sparsemm(ovrlp_compr, temp1, temp2, norb, norbp, sparsemat)
!!              do iorb=1,norbp
!!                  !!do jorb=1,norb
!!                  !!    write(440+iproc,*) iorb,jorb,temp2(jorb,iorb)
!!                  !!end do
!!              end do
!!              call daxpy(norb*norbp, cc(i), temp2, 1, ovrlp(1,isorb+1), 1) 
!!              call dcopy(norb*norbp, temp2, 1, temp1, 1)
!!          end do
!!          !!do iorb=1,norb
!!          !!    do jorb=1,norb
!!          !!        write(480+iproc,*) iorb,jorb,ovrlp(jorb,iorb)
!!          !!    end do
!!          !!end do
!!
!!          iall=-product(shape(temp1))*kind(temp1)
!!          deallocate(temp1, stat=istat)
!!          call memocc(istat, iall, 'temp1', subname)
!!          iall=-product(shape(temp2))*kind(temp2)
!!          deallocate(temp2, stat=istat)
!!          call memocc(istat, iall, 'temp2', subname)
!!
!!      end if
!!
!!      call mpiallred(ovrlp(1,1), norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!
!!
!!
!!      !!write(*,'(1x,a)') 'ERROR: methTransformOrder must be 0 or 1!'
!!      !!stop

  end if

  call timing(iproc,'lovrlp^-1/2   ','OF')

end subroutine overlapPowerMinusOneHalf_old
