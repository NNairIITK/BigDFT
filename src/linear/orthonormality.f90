!> @file
!! Orthonormalization
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, lzd, ovrlp, inv_ovrlp, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparseMatrix),intent(inout) :: ovrlp
  type(sparseMatrix),intent(in) :: inv_ovrlp ! just using structure for now
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall, iorb, jorb
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  type(sparseMatrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='orthonormalizeLocalized'

  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'

  call nullify_sparsematrix(inv_ovrlp_half)
  call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, subname)
  allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)

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

          call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.

      end if
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)

      if (methTransformOverlap==-1) then
          call overlap_power_minus_one_half_parallel(iproc, nproc, bigdft_mpi%mpi_comm, orbs, ovrlp, inv_ovrlp_half)
          !DEBUG
          !allocate(inv_ovrlp_half%matrix(orbs%norb,orbs%norb))
          !call uncompressMatrix(inv_ovrlp_half)
          !if (iproc==0) then
          !do iorb=1,orbs%norb
          !do jorb=1,orbs%norb
          !write(10,*) iorb,jorb,inv_ovrlp_half%matrix(iorb,jorb)
          !end do
          !end do
          !end if
          !deallocate(inv_ovrlp_half%matrix)
          !call overlapPowerMinusOneHalf(iproc, nproc, bigdft_mpi%mpi_comm, 0, orthpar%blocksize_pdsyev, &
          !    orthpar%blocksize_pdgemm, orbs%norb, ovrlp, inv_ovrlp_half)

          !allocate(inv_ovrlp_half%matrix(orbs%norb,orbs%norb))
          !call uncompressMatrix(inv_ovrlp_half)
          !if (iproc==0) then
          !do iorb=1,orbs%norb
          !do jorb=1,orbs%norb
          !write(11,*) iorb,jorb,inv_ovrlp_half%matrix(iorb,jorb)
          !end do
          !end do
          !end if
          !deallocate(inv_ovrlp_half%matrix)
          !call mpi_finalize(istat)
          !stop
          !END DEBUG
      else
          call overlapPowerMinusOneHalf(iproc, nproc, bigdft_mpi%mpi_comm, methTransformOverlap, orthpar%blocksize_pdsyev, &
              orthpar%blocksize_pdgemm, orbs%norb, ovrlp, inv_ovrlp_half)
      end if

      allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psittemp_c, 'psittemp_c', subname)
      allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psittemp_f, 'psittemp_f', subname)

      call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
      call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)
      call build_linear_combination_transposed(collcom, inv_ovrlp_half, &
           psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
      allocate(norm(orbs%norb), stat=istat)
      call memocc(istat, norm, 'norm', subname)
      call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)

      iall=-product(shape(norm))*kind(norm)
      deallocate(norm, stat=istat)
      call memocc(istat, iall, 'norm', subname)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

      iall=-product(shape(psittemp_c))*kind(psittemp_c)
      deallocate(psittemp_c, stat=istat)
      call memocc(istat, iall, 'psittemp_c', subname)
      iall=-product(shape(psittemp_f))*kind(psittemp_f)
      deallocate(psittemp_f, stat=istat)
      call memocc(istat, iall, 'psittemp_f', subname)
  end do

  call deallocate_sparseMatrix(inv_ovrlp_half, subname)

end subroutine orthonormalizeLocalized


! can still tidy this up more when tmblarge is removed
! use sparsity of density kernel for all inverse quantities
subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, npsidim_orbs, npsidim_comp, orbs, collcom, orthpar, &
           correction_orthoconstraint, linmat, lphi, lhphi, lagmat, psit_c, psit_f, hpsit_c, hpsit_f, &
           can_use_transposed, overlap_calculated)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim_orbs, npsidim_comp
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_Data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  integer,intent(in) :: correction_orthoconstraint
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(in) :: lphi
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(inout) :: lhphi
  type(SparseMatrix),intent(inout) :: lagmat
  real(kind=8),dimension(:),pointer :: psit_c, psit_f, hpsit_c, hpsit_f
  logical,intent(inout) :: can_use_transposed, overlap_calculated
  type(linear_matrices),intent(inout) :: linmat ! change to ovrlp and inv_ovrlp, and use inv_ovrlp instead of denskern

  ! Local variables
  integer :: istat, iall, iorb, jorb, ii, ii_trans
  type(SparseMatrix) :: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans, tmp_mat
  character(len=*),parameter :: subname='orthoconstraintNonorthogonal'


  ! ASSUME denskern sparsity pattern is symmetric
  ! create ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans with sparsity pattern of denskern here
  ! this is slight overkill for no orthoconstraint correction, think about going back to just matrices
  ! also isn't going to work unless denskern sparsity = lagmat sparsity...
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

      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.
  end if

  ! It is assumed that this routine is called with the transposed gradient ready if it is associated...
  if(.not.associated(hpsit_c)) then
      allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
 
      allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
 
     call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lhphi, hpsit_c, hpsit_f, lzd)
  end if

  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat)

  if (correction_orthoconstraint==0) then !not correctly working
      if(overlap_calculated) stop 'overlap_calculated should be wrong... To be modified later'

      ! problem here as psit match tmblarge whereas linmat%ovrlp matches tmb
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, linmat%ovrlp)

      allocate(linmat%ovrlp%matrix(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, linmat%ovrlp%matrix, 'linmat%ovrlp%matrix', subname)

      call uncompressMatrix(iproc,linmat%ovrlp)

      allocate(lagmat%matrix(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, lagmat%matrix, 'lagmat%matrix', subname)

      call uncompressMatrix(iproc,lagmat)

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

      call compress_matrix_for_allreduce(iproc,ovrlp_minus_one_lagmat)

      iall=-product(shape(ovrlp_minus_one_lagmat%matrix))*kind(ovrlp_minus_one_lagmat%matrix)
      deallocate(ovrlp_minus_one_lagmat%matrix, stat=istat)
      call memocc(istat, iall, 'ovrlp_minus_one_lagmat%matrix', subname)
 
      call compress_matrix_for_allreduce(iproc,ovrlp_minus_one_lagmat_trans)

      iall=-product(shape(ovrlp_minus_one_lagmat_trans%matrix))*kind(ovrlp_minus_one_lagmat_trans%matrix)
      deallocate(ovrlp_minus_one_lagmat_trans%matrix, stat=istat)
      call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans%matrix', subname)
  end if

  call nullify_sparseMatrix(tmp_mat)
  call sparse_copy_pattern(ovrlp_minus_one_lagmat,tmp_mat,subname)

  allocate(tmp_mat%matrix_compr(tmp_mat%nvctr), stat=istat)
  call memocc(istat, tmp_mat%matrix_compr, 'tmp_mat%matrix_compr', subname)

  do jorb=1,orbs%norb
     do iorb=1,orbs%norb
          ii_trans = ovrlp_minus_one_lagmat_trans%matrixindex_in_compressed(jorb, iorb)
          ii = ovrlp_minus_one_lagmat%matrixindex_in_compressed(iorb, jorb)
          if (ii==0.or.ii_trans==0) cycle
          tmp_mat%matrix_compr(ii)=-0.5d0*ovrlp_minus_one_lagmat%matrix_compr(ii) &
               -0.5d0*ovrlp_minus_one_lagmat_trans%matrix_compr(ii_trans)
      end do
  end do

  call build_linear_combination_transposed(collcom, tmp_mat, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)
  call deallocate_sparseMatrix(tmp_mat, subname)

  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, hpsit_c, hpsit_f, lhphi, lzd)

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


subroutine overlapPowerMinusOne(iproc, nproc, iorder, blocksize, norb, ovrlp, inv_ovrlp)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, norb
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp
  real(kind=8),dimension(norb,norb),intent(out) :: inv_ovrlp
  
  ! Local variables
  integer :: iorb, jorb, info

  call timing(iproc,'lovrlp^-1     ','ON')

  ! don't want to destroy ovrlp, and inv_ovrlp will eventually use sparse matrices
  call vcopy(norb*norb,ovrlp(1,1),1,inv_ovrlp(1,1),1)

  if(iorder==0) then

      ! Exact inversion
      if (blocksize<0) then
          call dpotrf('l', norb, inv_ovrlp(1,1), norb, info)
          if(info/=0) then
              write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
              stop
          end if
          call dpotri('l', norb, inv_ovrlp(1,1), norb, info)
          if(info/=0) then
              write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
              stop
          end if
      else
          call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp(1,1), norb)
          call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp(1,1), norb)
      end if
  
  else if(iorder==1) then
       
      ! Taylor expansion up to first order.
      do iorb=1,norb
          do jorb=1,norb
              if(iorb==jorb) then
                  inv_ovrlp(jorb,iorb) = 2.d0 - ovrlp(jorb,iorb)
              else
                  inv_ovrlp(jorb,iorb) = -ovrlp(jorb,iorb)
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
           blocksize_pdgemm, norb, ovrlp, inv_ovrlp_half)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb
  type(sparseMatrix),intent(inout) :: ovrlp
  type(sparseMatrix),intent(inout) :: inv_ovrlp_half

  ! Local variables
  integer :: lwork, istat, iall, iorb, jorb, info, iseg, iiorb, jjorb, ii, ii_inv
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

      call uncompressMatrix(iproc,ovrlp)
  
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
              tempArr(jorb,iorb,1)=inv_ovrlp_half%matrix(jorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
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

      call compress_matrix_for_allreduce(iproc,inv_ovrlp_half)

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
      !ii=0
      !do iseg=1,ovrlp%nseg
      !    do jorb=ovrlp%keyg(1,iseg),ovrlp%keyg(2,iseg)
      !        ii=ii+1
      !        iiorb = (jorb-1)/norb + 1
      !        jjorb = jorb - (iiorb-1)*norb
           do ii=1,ovrlp%nvctr
              iiorb = ovrlp%orb_from_index(ii,1)
              jjorb = ovrlp%orb_from_index(ii,2)

              ii_inv = inv_ovrlp_half%matrixindex_in_compressed(iiorb,jjorb) ! double check this order
              if(iiorb==jjorb) then
                  inv_ovrlp_half%matrix_compr(ii_inv)=1.5d0-.5d0*ovrlp%matrix_compr(ii)
              else
                  inv_ovrlp_half%matrix_compr(ii_inv)=-.5d0*ovrlp%matrix_compr(ii)
              end if
          end do
      !end do
  else
      
      stop 'deprecated'

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


subroutine overlap_power_minus_one_half_parallel(iproc, nproc, comm, orbs, ovrlp, inv_ovrlp_half)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, comm
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: ovrlp
  type(sparseMatrix),intent(inout) :: inv_ovrlp_half

  ! Local variables
  integer :: iend, i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
  integer :: iiorb, ierr, ii, iseg, ind
  real(kind=8),dimension(:,:),allocatable :: ovrlp_tmp, ovrlp_tmp_inv_half
  logical,dimension(:),allocatable :: in_neighborhood
  character(len=*),parameter :: subname='overlap_power_minus_one_half_parallel'


  allocate(in_neighborhood(orbs%norb), stat=istat)
  call memocc(istat, in_neighborhood, 'in_neighborhood', subname)
  call to_zero(inv_ovrlp_half%nvctr, inv_ovrlp_half%matrix_compr(1))

  !DEBUG
  !if (iproc==0) then
  !   jjorb=0
  !   do jorb=1,orbs%norb
  !      do korb=1,orbs%norb
  !         ind = ovrlp%matrixindex_in_compressed(korb, jorb)
  !         if (ind>0) then
  !            print*, korb,jorb,ovrlp%matrix_compr(ind)
  !         else
  !            print*, korb,jorb,0.0d0
  !         end if
  !         !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
  !      end do
  !   end do
  !end if
  !call mpi_barrier(bigdft_mpi%mpi_comm,istat)

  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     !ilr=orbs%inwhichlocreg(iiorb)
     ! We are at the start of a new atom
     ! Count all orbitals that are in the neighborhood

     iseg=inv_ovrlp_half%istsegline(iiorb)
     iend=iiorb*orbs%norb
     n=0
     in_neighborhood(:)=.false.
     do 
        do i=inv_ovrlp_half%keyg(1,iseg),inv_ovrlp_half%keyg(2,iseg)
           ii=i-(iiorb-1)*orbs%norb
           in_neighborhood(ii)=.true.
           n=n+1
        end do
        iseg=iseg+1
        if (iseg>inv_ovrlp_half%nseg) exit
        if (inv_ovrlp_half%keyg(1,iseg)>iend) exit
     end do

     allocate(ovrlp_tmp(n,n), stat=istat)
     call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
     call to_zero(n*n, ovrlp_tmp(1,1))

     jjorb=0
     do jorb=1,orbs%norb
        if (.not.in_neighborhood(jorb)) cycle
        jjorb=jjorb+1
        kkorb=0
        do korb=1,orbs%norb
           if (.not.in_neighborhood(korb)) cycle
           kkorb=kkorb+1
           ind = ovrlp%matrixindex_in_compressed(korb, jorb)
           if (ind>0) then
              ovrlp_tmp(kkorb,jjorb)=ovrlp%matrix_compr(ind)
           else
              ovrlp_tmp(kkorb,jjorb)=0.d0
           end if
           !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
        end do
     end do
          
     allocate(ovrlp_tmp_inv_half(n,n))
     call memocc(istat, ovrlp_tmp_inv_half, 'ovrlp_tmp_inv_half', subname)

     !if (iiorb==orbs%norb) then
     !print*,''
     !print*,'ovrlp_tmp',n,iiorb
     !do jorb=1,n
     !print*,jorb,ovrlp_tmp(:,jorb)
     !end do
     !end if

     ! Calculate S^-1/2 for the small overlap matrix
     call overlapPowerMinusOneHalf_old(iproc, nproc, comm, 0, -8, -8, n, ovrlp_tmp, ovrlp_tmp_inv_half)

     !if (iiorb==orbs%norb) then
     !print*,''
     !print*,'inv_ovrlp_tmp',n,iiorb
     !do jorb=1,n
     !print*,jorb,ovrlp_tmp_inv_half(:,jorb)
     !end do
     !end if

     jjorb=0
     do jorb=1,orbs%norb
        if (.not.in_neighborhood(jorb)) cycle
        jjorb=jjorb+1
        kkorb=0
        if (jorb==iiorb) then
           do korb=1,orbs%norb
              if (.not.in_neighborhood(korb)) cycle
              kkorb=kkorb+1
              ind = inv_ovrlp_half%matrixindex_in_compressed(korb,jorb)
              if (ind>0) then
                 inv_ovrlp_half%matrix_compr(ind)=ovrlp_tmp_inv_half(kkorb,jjorb)
                 !if (iiorb==orbs%norb) print*,'problem here?!',iiorb,kkorb,jjorb,korb,jorb,ind,ovrlp_tmp_inv_half(kkorb,jjorb)
              end if
              !write(1300+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
           end do
           exit !no need to keep looping
        end if
     end do

     iall=-product(shape(ovrlp_tmp_inv_half))*kind(ovrlp_tmp_inv_half)
     deallocate(ovrlp_tmp_inv_half, stat=istat)
     call memocc(istat, iall, 'ovrlp_tmp_inv_half', subname)

     iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
     deallocate(ovrlp_tmp, stat=istat)
     call memocc(istat, iall, 'ovrlp_tmp', subname)
  end do

  call mpiallred(inv_ovrlp_half%matrix_compr(1), inv_ovrlp_half%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  iall=-product(shape(in_neighborhood))*kind(in_neighborhood)
  deallocate(in_neighborhood, stat=istat)
  call memocc(istat, iall, 'in_neighborhood', subname)

  !if (iproc==0) then
  !   jjorb=0
  !   do jorb=1,orbs%norb
  !      do korb=1,orbs%norb
  !         ind = inv_ovrlp_half%matrixindex_in_compressed(korb, jorb)
  !         if (ind>0) then
  !            print*, korb,jorb,inv_ovrlp_half%matrix_compr(ind)
  !         else
  !            print*, korb,jorb,0.0d0
  !         end if
  !         !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
  !      end do
  !   end do
  !end if
  !call mpi_finalize(ind)
  !stop

end subroutine overlap_power_minus_one_half_parallel


! Should be used if sparsemat is not available... to be cleaned
subroutine overlapPowerMinusOneHalf_old(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, &
           blocksize_pdgemm, norb, ovrlp, inv_ovrlp_half, orbs)
  use module_base
  use module_types
  use module_interfaces, fake_name => overlapPowerMinusOneHalf_old
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, norb, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp
  real(kind=8),dimension(norb,norb),intent(inout) :: inv_ovrlp_half
  type(orbitals_data), optional, intent(in) :: orbs
  
  ! Local variables
  integer :: lwork, istat, iall, iorb, jorb, info, iiorb, ierr
  character(len=*),parameter :: subname='overlapPowerMinusOneHalf'
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_halfp
  real(kind=8),dimension(:,:,:),allocatable :: tempArr

  call timing(iproc,'lovrlp^-1/2old','ON')


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
          
          allocate(work(10), stat=istat)
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
              write(*,'(a,i0)') 'ERROR in dsyev, info=', info
              stop
          end if
          iall=-product(shape(work))*kind(work)
          deallocate(work, stat=istat)
          call memocc(istat, iall, 'work', subname)
      end if

      ! Calculate S^{-1/2}. 
      ! First calculate ovrlp*diag(1/sqrt(eval)) (ovrlp is the diagonalized overlap
      ! matrix and diag(1/sqrt(eval)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
      do iorb=1,norb
          do jorb=1,norb
              tempArr(jorb,iorb,1)=inv_ovrlp_half(jorb,iorb)/sqrt(abs(eval(iorb)))
          end do
      end do
      
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
  else

     ! Taylor expansion up to 4th order
      stop 'must fix this'

  end if

  call timing(iproc,'lovrlp^-1/2old','OF')

end subroutine overlapPowerMinusOneHalf_old
