!> @file
!! Orthonormalization
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Orthonormalized the localized orbitals
subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed, foe_obj)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix
  use foe_base, only: foe_data
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparse_matrix),intent(inout) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed
  type(foe_data),intent(in) :: foe_obj

  ! Local variables
  integer :: it, istat, iall
  !integer :: irow, ii, iorb, jcol, jorb
  real(kind=8), dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparse_matrix) :: inv_ovrlp_half
  character(len=*), parameter :: subname='orthonormalizeLocalized'
  !real(kind=8), dimension(orbs%norb,orbs%norb) :: tempmat
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: error
  logical :: ovrlp_associated, inv_ovrlp_associated
  type(matrices) :: ovrlp_, inv_ovrlp_half_


  !inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')

  inv_ovrlp_half_ = matrices_null()
  call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_)



  if(.not.can_use_transposed) then
      if(associated(psit_c)) then
          call f_free_ptr(psit_c)
      end if
      if(associated(psit_f)) then
          call f_free_ptr(psit_f)
      end if
      psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
      psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')

      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.

  end if

  ovrlp_ = matrices_null()
  call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
  ! This can then be deleted if the transition to the new type has been completed.


  if (methTransformOverlap==-1) then
      call overlap_power_minus_one_half_parallel(iproc, nproc, 0, orbs, ovrlp, ovrlp_, inv_ovrlp_half, inv_ovrlp_half_)
      !!inv_ovrlp_half_%matrix_compr = inv_ovrlp_half%matrix_compr
  else
      !ovrlp%matrix_compr=ovrlp_%matrix_compr
      call overlapPowerGeneral(iproc, nproc, methTransformOverlap, -2, &
           orthpar%blocksize_pdgemm, &
           imode=1, ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
           ovrlp_mat=ovrlp_, inv_ovrlp_mat=inv_ovrlp_half_, &
           check_accur=.false.)!!, &
  end if

  call deallocate_matrices(ovrlp_)

  psittemp_c = f_malloc(sum(collcom%nrecvcounts_c),id='psittemp_c')
  psittemp_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psittemp_f')

  call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
  call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)

  call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_, &
       psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)

  norm = f_malloc(orbs%norb,id='norm')
  call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)

  call f_free(norm)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

  call f_free(psittemp_c)
  call f_free(psittemp_f)

  call f_free_ptr(inv_ovrlp_half%matrix_compr)

  call deallocate_matrices(inv_ovrlp_half_)

end subroutine orthonormalizeLocalized


! can still tidy this up more when tmblarge is removed
! use sparsity of density kernel for all inverse quantities
subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, npsidim_orbs, npsidim_comp, orbs, collcom, orthpar, &
           correction_orthoconstraint, linmat, lphi, lhphi, lagmat, lagmat_, psit_c, psit_f, hpsit_c, hpsit_f, &
           can_use_transposed, overlap_calculated, experimental_mode, norder_taylor)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  use yaml_output
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: matrices_null, allocate_matrices, deallocate_matrices, sparsematrix_malloc, &
                               sparsematrix_malloc_ptr, DENSE_FULL, DENSE_PARALLEL, SPARSE_FULL, SPARSEMM_SEQ, &
                               assignment(=)
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: uncompress_matrix, uncompress_matrix_distributed, compress_matrix_distributed, &
                          sequential_acces_matrix_fast, sparsemm
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim_orbs, npsidim_comp
  type(local_zone_descriptors),intent(in) :: lzd
  !type(orbitals_Data),intent(in) :: orbs
  type(orbitals_Data),intent(inout) :: orbs !temporary inout
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  integer,intent(in) :: correction_orthoconstraint
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(in) :: lphi
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(inout) :: lhphi
  type(sparse_matrix),intent(inout) :: lagmat
  type(matrices),intent(out) :: lagmat_
  real(kind=8),dimension(:),pointer :: psit_c, psit_f, hpsit_c, hpsit_f
  logical,intent(inout) :: can_use_transposed, overlap_calculated
  type(linear_matrices),intent(inout) :: linmat ! change to ovrlp and inv_ovrlp, and use inv_ovrlp instead of denskern
  logical,intent(in) :: experimental_mode
  integer,intent(in) :: norder_taylor

  ! Local variables
  integer :: istat, iall, iorb, jorb, ii, ii_trans, irow, jcol, info, lwork, jj
  real(kind=8) :: error
  real(kind=8),dimension(:),allocatable :: tmp_mat_compr, lagmat_tmp_compr, work
  character(len=*),parameter :: subname='orthoconstraintNonorthogonal'
  real(kind=8),dimension(:,:),allocatable :: tmp_mat, tmp_mat2, tmp_mat3
  integer,dimension(:),allocatable :: ipiv
  type(matrices) :: inv_ovrlp_
  real(8),dimension(:),allocatable :: inv_ovrlp_seq
  real(8),dimension(:,:),allocatable :: lagmatp, inv_lagmatp

  ! removed option for correction orthoconstrain for now
  !if (correction_orthoconstraint==0) stop 'correction_orthoconstraint not working'

  if(.not. can_use_transposed) then
      psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
      psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')

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

  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat, lagmat_)
  ! This can then be deleted if the transition to the new type has been completed.
  lagmat%matrix_compr=lagmat_%matrix_compr


  tmp_mat_compr = sparsematrix_malloc(lagmat,iaction=SPARSE_FULL,id='tmp_mat_compr')
call timing(iproc,'misc','ON')
  do ii=1,lagmat%nvctr
     iorb = lagmat%orb_from_index(1,ii)
     jorb = lagmat%orb_from_index(2,ii)
     ii_trans=matrixindex_in_compressed(lagmat,jorb, iorb)
     tmp_mat_compr(ii)=-0.5d0*lagmat_%matrix_compr(ii)-0.5d0*lagmat_%matrix_compr(ii_trans)
     ! SM: This is a hack, should use another variable
     !if (.false..or.correction_orthoconstraint==2) then
     if (.false..and.correction_orthoconstraint==2) then
         if (iproc==0 .and. ii==1) write(*,*) 'only normalization constraint'
         if (iorb/=jorb) then
             tmp_mat_compr(ii)=0.d0
         end if
     end if
     if (experimental_mode) then
         if (iorb==jorb) then
             if (iproc==0 .and. iorb==1) then
                 call yaml_warning('EXPERIMENTAL: modify eval')
                 call yaml_newline()
             end if
             orbs%eval(iorb)=lagmat_%matrix_compr(ii)
         end if
     end if
  end do

  ! NEW: reactivate correction for non-orthogonality ##########
  if (correction_orthoconstraint==0) then
      !@NEW
      if (iproc==0) call yaml_map('correction orthoconstraint',.true.)
      inv_ovrlp_ = matrices_null()
      call allocate_matrices(linmat%l, allocate_full=.false., &
           matname='inv_ovrlp_', mat=inv_ovrlp_)
      call overlapPowerGeneral(iproc, nproc, norder_taylor, 1, -1, &
           imode=1, ovrlp_smat=linmat%s, inv_ovrlp_smat=linmat%l, &
           ovrlp_mat=linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp_, &
           check_accur=.true., error=error)

      inv_ovrlp_seq = sparsematrix_malloc(linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_seq')
      lagmatp = sparsematrix_malloc(linmat%m, iaction=DENSE_PARALLEL, id='lagmatp')
      inv_lagmatp = sparsematrix_malloc(linmat%m, iaction=DENSE_PARALLEL, id='inv_lagmatp')
      call sequential_acces_matrix_fast(linmat%l, inv_ovrlp_%matrix_compr, inv_ovrlp_seq)
      call uncompress_matrix_distributed(iproc, linmat%m, tmp_mat_compr, lagmatp)
      call sparsemm(linmat%l, inv_ovrlp_seq, lagmatp, inv_lagmatp)
      call compress_matrix_distributed(iproc, linmat%m, inv_lagmatp, tmp_mat_compr)
      call f_free(inv_ovrlp_seq)
      call f_free(lagmatp)
      call f_free(inv_lagmatp)
      call deallocate_matrices(inv_ovrlp_)
      !@ENDNEW


   !!!   ! WARNING: it is mandatory that the overlap matrix has been calculated before
   !!!   !!call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, linmat%ovrlp)
   !!!   !if (iproc==0) write(*,*) 'correction orthoconstraint'
   !!!   if (iproc==0) call yaml_map('correction orthoconstraint',.true.)
   !!!   linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(linmat%s, DENSE_FULL, id='linmat%ovrlp_%matrix')
   !!!   call uncompress_matrix(iproc, linmat%s, inmat=linmat%ovrlp_%matrix_compr, outmat=linmat%ovrlp_%matrix)
   !!!   allocate(tmp_mat(orbs%norb,orbs%norb))
   !!!   allocate(tmp_mat2(orbs%norb,orbs%norb))
   !!!   call to_zero(orbs%norb**2, tmp_mat(1,1))
   !!!   do ii=1,lagmat%nvctr
   !!!      irow = lagmat%orb_from_index(1,ii)
   !!!      jcol = lagmat%orb_from_index(2,ii)
   !!!      tmp_mat(irow,jcol)=tmp_mat_compr(ii)
   !!!   end do
   !!!   allocate(ipiv(orbs%norb))
   !!!   lwork=10*orbs%norb
   !!!   allocate(work(lwork))
   !!!   allocate(tmp_mat3(orbs%norb,orbs%norb))
   !!!   tmp_mat3=linmat%ovrlp_%matrix

   !!!   call dgetrf(orbs%norb, orbs%norb, linmat%ovrlp_%matrix, orbs%norb, ipiv, info)
   !!!   call dgetri(orbs%norb, linmat%ovrlp_%matrix, orbs%norb, ipiv, work, lwork, info)

   !!!   !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, linmat%ovrlp%matrix, orbs%norb, tmp_mat3, orbs%norb, 0.d0, tmp_mat2, orbs%norb)
   !!!   !!if (iproc==0) then
   !!!   !!  do iorb=1,orbs%norb
   !!!   !!    do jorb=1,orbs%norb
   !!!   !!      write(*,'(a,2i8,es14.5)') 'iorb, jorb, tmp_mat2(iorb,jorb)', iorb, jorb, tmp_mat2(iorb,jorb)
   !!!   !!    end do
   !!!   !!  end do
   !!!   !!end if

   !!!   ! This is the original
   !!!   call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, linmat%ovrlp_%matrix, orbs%norb, &
   !!!        tmp_mat, orbs%norb, 0.d0, tmp_mat2, orbs%norb)

   !!!   !!! Test
   !!!   !!call dgemm('n', 't', orbs%norb, orbs%norb, orbs%norb, 1.d0, tmp_mat, orbs%norb, linmat%ovrlp%matrix, orbs%norb, 0.d0, tmp_mat2, orbs%norb)

   !!!   do jj=1,lagmat%nvctr
   !!!      irow = lagmat%orb_from_index(1,jj)
   !!!      jcol = lagmat%orb_from_index(2,jj)
   !!!      !!if (iproc==0) write(*,'(a,3i8,2es16.6)') 'jj, irow, jcol, tmp_mat_compr(jj), tmp_mat2(irow,jcol)', &
   !!!      !!                                          jj, irow, jcol, tmp_mat_compr(jj), tmp_mat2(irow,jcol)
   !!!      tmp_mat_compr(jj)=tmp_mat2(irow,jcol)
   !!!   end do
   !!!   call f_free_ptr(linmat%ovrlp_%matrix)
   !!!   deallocate(tmp_mat)
   !!!   deallocate(tmp_mat2)
   !!!   deallocate(ipiv)
   !!!   deallocate(work)
  end if
  !! ##########################################################

  lagmat_tmp_compr = sparsematrix_malloc(lagmat,iaction=SPARSE_FULL,id='lagmat_tmp_compr')

  call vcopy(lagmat%nvctr,lagmat_%matrix_compr(1),1,lagmat_tmp_compr(1),1) ! need to keep a copy
  call vcopy(lagmat%nvctr,tmp_mat_compr(1),1,lagmat_%matrix_compr(1),1)

  call f_free(tmp_mat_compr)

call timing(iproc,'misc','OF')

  !lagmat_ = matrices_null()
  !call allocate_matrices(lagmat, allocate_full=.false., matname='lagmat_', mat=lagmat_)
  !lagmat_%matrix_compr = lagmat%matrix_compr
  call build_linear_combination_transposed(collcom, lagmat, lagmat_, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)
  !call deallocate_matrices(lagmat_)






  !call build_linear_combination_transposed(collcom, tmp_mat, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)

  !!! TEST ORTHOGONALITY OF GRADIENT AND TMBs ##############################
  !!call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, linmat%ovrlp)
  !!allocate(linmat%ovrlp%matrix(orbs%norb,orbs%norb))
  !!call uncompress_matrix(iproc,linmat%ovrlp)
  !!if (iproc==0) then
  !!  do iorb=1,orbs%norb
  !!    do jorb=1,orbs%norb
  !!      write(*,'(a,2i8,es16.6)') 'iorb, jorb, linmat%ovrlp%matrix(jorb,iorb)', iorb, jorb, linmat%ovrlp%matrix(jorb,iorb)
  !!    end do
  !!  end do
  !!end if
  !!deallocate(linmat%ovrlp%matrix)
  !!! END TEST #############################################################

  call vcopy(lagmat%nvctr,lagmat_tmp_compr(1),1,lagmat_%matrix_compr(1),1)

  call f_free(lagmat_tmp_compr)

  !call deallocate_sparse_matrix(tmp_mat, subname)

  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, hpsit_c, hpsit_f, lhphi, lzd)

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


!S^-1 exact only works for symmetric matrices
!BOTH sparse matrices must be present together and inv_ovrlp should be nullified pointer, NOT inv_ovrlp_smat%matrix
!when sparse matrices present, check is performed to see whether %matrix is allocated so that its allocated status remains unchanged
!contents of %matrix not guaranteed to be correct though - inv_ovrlp_smat%can_use_dense set accordingly
!power: -2 -> S^-1/2, 2 -> S^1/2, 1 -> S^-1
subroutine overlapPowerGeneral(iproc, nproc, iorder, power, blocksize, imode, &
           ovrlp_smat, inv_ovrlp_smat, ovrlp_mat, inv_ovrlp_mat, check_accur, &
           error)
     !!foe_nseg, foe_kernel_nsegline, foe_istsegline, foe_keyg)
  use module_base
  use module_types
  use module_interfaces, except_this_one => overlapPowerGeneral
  use sparsematrix_base, only: sparse_matrix, &
                          sparsematrix_malloc_ptr, sparsematrix_malloc, sparsematrix_malloc0, sparsematrix_malloc0_ptr, &
                          assignment(=), &
                          SPARSE_FULL, DENSE_PARALLEL, DENSE_FULL, SPARSEMM_SEQ
  use sparsematrix, only: compress_matrix, uncompress_matrix, transform_sparse_matrix, &
                          compress_matrix_distributed, uncompress_matrix_distributed, &
                          sequential_acces_matrix_fast, sparsemm
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, power
  integer,intent(in) :: imode
  type(sparse_matrix),intent(inout) :: ovrlp_smat, inv_ovrlp_smat
  type(matrices),intent(inout) :: ovrlp_mat, inv_ovrlp_mat
  logical,intent(in) :: check_accur
  real(kind=8),intent(out),optional :: error
  
  ! Local variables
  integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i, its, maxits
  integer :: matrixindex_in_compressed, nmaxvalk
  real(kind=8), dimension(:,:), pointer :: ovrlpminonep, ovrlpminone, inv_ovrlpp, ovrlppowerp, ovrlppoweroldp
  real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp, ovrlp_local, inv_ovrlp_local
  real(kind=8) :: factor, newfactor
  logical :: ovrlp_allocated, inv_ovrlp_allocated

  ! new for sparse taylor
  integer :: nout, nseq, nmaxsegk, nmaxval
  integer,dimension(:),allocatable :: ivectorindex
  integer,dimension(:,:),pointer :: onedimindices
  integer,dimension(:,:,:),allocatable :: istindexarr
  real(kind=8),dimension(:),pointer :: ovrlpminone_sparse
  real(kind=8),dimension(:),allocatable :: ovrlp_compr_seq, ovrlpminone_sparse_seq, ovrlp_large_compr
  real(kind=8),dimension(:),allocatable :: invovrlp_compr_seq
  real(kind=8),dimension(:,:),allocatable :: ovrlpminoneoldp, invovrlpp, ovrlp_largep
  real(kind=8),dimension(:,:),allocatable :: Amat12p, Amat21p, Amat21
  real(kind=8),dimension(:,:),pointer :: Amat12, Amat11p, Amat22p
  real(kind=8),dimension(:),pointer :: Amat12_compr
  real(kind=8),dimension(:),allocatable :: Amat21_compr, Amat12_seq, Amat21_seq
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2


  !!write(*,*) 'iorder',iorder


  call f_routine(id='overlapPowerGeneral')
  call timing(iproc,'lovrlp^-1     ','ON')


  if (iproc==0) then
      call yaml_newline()
      call yaml_open_sequence('overlap manipulation routine')
      if (imode==SPARSE) then
          call yaml_map('mode','sparse')
      else if (imode==DENSE) then
          call yaml_map('mode','dense')
      end if
      call yaml_map('power',power)
      call yaml_map('order',iorder)
      call yaml_close_sequence()
  end if




  ! Perform a check of the arguments

  if (imode/=SPARSE .and. imode/=DENSE) stop 'wrong imode'

  if (imode==DENSE) then
      if (.not.associated(ovrlp_mat%matrix)) stop 'ovrlp_mat%matrix not associated'
      if (.not.associated(inv_ovrlp_mat%matrix)) stop 'inv_ovrlp_mat%matrix not associated'
  end if
  
  if (check_accur) then
      if (.not.present(error)) stop 'error not present'
  end if

  if (power/=-2 .and. power/=1 .and. power/=2) stop 'wrong value of power'

  if (nproc/=1 .and. nproc/=bigdft_mpi%nproc) stop 'wrong value of nproc'

  ! Decide whether this routine is called in parallel or in serial.
  ! If parallel, take the default values from orbs, otherwise adjust them.
  !if (nproc>1) then
      norbp=ovrlp_smat%nfvctrp
      isorb=ovrlp_smat%isfvctr
  !else
  !    norbp=norb
  !    isorb=0
  !end if

  sparse_dense: if (imode==DENSE) then
      if (iorder==0) then
          call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr,ovrlp_mat%matrix(1,1),1,inv_ovrlp_mat%matrix(1,1),1)
          if (power==1) then
             if (blocksize<0) then
                call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_mat%matrix)
             else
                stop 'check if working - upper half may not be filled'
                call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                     ovrlp_smat%nfvctr, inv_ovrlp_mat%matrix(1,1), ovrlp_smat%nfvctr)
                call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                     ovrlp_smat%nfvctr, inv_ovrlp_mat%matrix(1,1), ovrlp_smat%nfvctr)
             end if
          else if (power==2) then
             !if (nproc>1) then
                 call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                      blocksize,.true.,inv_ovrlp_mat%matrix,inv_ovrlp_smat)
             !else
             !    call overlap_plus_minus_one_half_exact(ovrlp_smat%nfvctr,blocksize,.true.,inv_ovrlp_mat%matrix)
             !end if
          else if (power==-2) then
              !if (nproc>1) then
                  call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                       blocksize,.false.,inv_ovrlp_mat%matrix,inv_ovrlp_smat)
              !else
              !    call overlap_plus_minus_one_half_exact(ovrlp_smat%nfvctr,blocksize,.false.,inv_ovrlp_mat%matrix)
              !end if
          end if
      else if (iorder<0) then
          ! sign approach as used in CP2K
          ! use 4 submatrices
          !if (nproc>1) then
              !Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat12p')
              !Amat21p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat21p')
              !Amat11p = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat11p')
          !else
              Amat12p = f_malloc((/ovrlp_smat%nfvctr,norbp/), id='Amat12p')
              Amat21p = f_malloc((/ovrlp_smat%nfvctr,norbp/), id='Amat21p')
              Amat11p = f_malloc_ptr((/ovrlp_smat%nfvctr,norbp/), id='Amat11p')
          !end if
          ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
          Amat22p=>Amat11p
          Amat12=>inv_ovrlp_mat%matrix
          !if (nproc>1) then
          !    Amat21=sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat21')
          !else
              Amat21=f_malloc0((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr/), id='Amat21')
          !end if

          call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr,ovrlp_mat%matrix(1,1),1,Amat12(1,1),1)
          do iorb=1,ovrlp_smat%nfvctr
              Amat21(iorb,iorb)=1.0d0
          end do

          ! calculate Xn+1=0.5*Xn*(3I-Xn**2)
          do its=1,abs(iorder)
              if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, -0.5d0, Amat12(1,1), &
                   ovrlp_smat%nfvctr, Amat21(1,isorb+1), ovrlp_smat%nfvctr, 0.0d0, Amat11p(1,1), ovrlp_smat%nfvctr)
              !call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, -0.5d0, Amat21(1,1), &
              !     ovrlp_smat%nfvctr, Amat12(1,isorb+1), ovrlp_smat%nfvctr, 0.0d0, Amat22p(1,1), ovrlp_smat%nfvctr)
              do iorb=1,norbp
                  Amat11p(iorb+isorb,iorb)=Amat11p(iorb+isorb,iorb)+1.5d0
              !    Amat22p(iorb+isorb,iorb)=Amat22p(iorb+isorb,iorb)+1.5d0
              end do
              if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.0d0, Amat12(1,1), &
                   ovrlp_smat%nfvctr, Amat22p(1,1), ovrlp_smat%nfvctr, 0.0d0, Amat12p(1,1), ovrlp_smat%nfvctr)
              if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.0d0, Amat21(1,1), &
                   ovrlp_smat%nfvctr, Amat11p(1,1), ovrlp_smat%nfvctr, 0.0d0, Amat21p(1,1), ovrlp_smat%nfvctr)
              !if(nproc > 1) then
                  call mpi_allgatherv(Amat12p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, Amat12, &
                       ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                       mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                  call mpi_allgatherv(Amat21p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, Amat21, &
                       ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                       mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
              !else
              !    call vcopy(ovrlp_smat%nfvctr**2,Amat12p(1,1),1,Amat12(1,1),1)
              !    call vcopy(ovrlp_smat%nfvctr**2,Amat21p(1,1),1,Amat21(1,1),1)
              !end if
          end do

          nullify(Amat22p)
          call f_free_ptr(Amat11p)

          if (power==1) then
              if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.0d0, Amat21(1,1), &
                   ovrlp_smat%nfvctr, Amat21p(1,1), ovrlp_smat%nfvctr, 0.0d0, Amat12p(1,1), ovrlp_smat%nfvctr)
              !if (nproc>1) then
                  call mpi_allgatherv(Amat12p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, inv_ovrlp_mat%matrix, &
                       ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                       mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
              !else
              !    call vcopy(ovrlp_smat%nfvctr**2, Amat12p(1,1), 1, inv_ovrlp_mat%matrix(1,1), 1)
              !end if
          !else if (power==2) then
          !   call vcopy(ovrlp_smat%nfvctr**2,Amat12(1,1),1,inv_ovrlp_mat%matrix(1,1),1)
          else if (power==-2) then
              call vcopy(ovrlp_smat%nfvctr**2,Amat21(1,1),1,inv_ovrlp_mat%matrix(1,1),1)
          end if

          call f_free(Amat12p)
          call f_free(Amat21p)
          nullify(Amat12)
          call f_free(Amat21)

      else
          if (iorder>1) then
              if (nproc>1) then
                  ovrlpminone => inv_ovrlp_mat%matrix
                  !ovrlpminonep = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlpminonep')
                  ovrlpminonep = f_malloc_ptr((/ovrlp_smat%nfvctr,norbp/), id='ovrlpminonep')
              else
                  ovrlpminone = f_malloc_ptr((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr/), id='ovrlpminone')
                  ovrlpminonep => ovrlpminone
              end if

              if (norbp>0) call matrix_minus_identity_dense(ovrlp_smat%nfvctr,isorb,norbp,ovrlp_mat%matrix(1,isorb+1),ovrlpminonep)

              !if (nproc>1) then
                  !ovrlppoweroldp = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlppoweroldp')
              !else
                  ovrlppoweroldp = f_malloc_ptr((/ovrlp_smat%nfvctr,norbp/), id='ovrlppoweroldp')
              !end if

              if (norbp>0) call vcopy(ovrlp_smat%nfvctr*norbp,ovrlpminonep(1,1),1,ovrlppoweroldp(1,1),1)

              if(nproc > 1) then
                  call mpi_allgatherv(ovrlpminonep, ovrlp_smat%nfvctr*norbp, mpi_double_precision, ovrlpminone, &
                       ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                       mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                  call f_free_ptr(ovrlpminonep)
              else
                  nullify(ovrlpminonep)
              end if

              !if (nproc>1) then
              !    ovrlppowerp = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlppowerp')
              !else
                  ovrlppowerp = f_malloc_ptr((/ovrlp_smat%nfvctr,norbp/), id='ovrlppowerp')
              !end if

              if (power==1) then
                  factor=-1.0d0
              else if (power==2) then
                  factor=0.5d0
              else if (power==-2) then
                  factor=-0.5d0
              end if
          end if

          if (nproc>1) then
              !inv_ovrlpp = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='inv_ovrlpp')
              inv_ovrlpp = f_malloc_ptr((/ovrlp_smat%nfvctr,norbp/), id='inv_ovrlpp')
          else
              inv_ovrlpp => inv_ovrlp_mat%matrix
          end if

          if (norbp>0) call first_order_taylor_dense(ovrlp_smat%nfvctr,isorb,norbp,power,ovrlp_mat%matrix(1,isorb+1),inv_ovrlpp)

          do i=2,iorder
              if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.d0, ovrlpminone(1,1), &
                   ovrlp_smat%nfvctr, ovrlppoweroldp(1,1), ovrlp_smat%nfvctr, 0.d0, ovrlppowerp(1,1), ovrlp_smat%nfvctr)
              factor=newfactor(power,i,factor)
              call daxpy(ovrlp_smat%nfvctr*norbp,factor,ovrlppowerp,1,inv_ovrlpp,1)
              if (i/=iorder.and.norbp>0) call vcopy(ovrlp_smat%nfvctr*norbp,ovrlppowerp(1,1),1,ovrlppoweroldp(1,1),1)
          end do

          if (iorder>1) then
              if(nproc > 1) then
                  nullify(ovrlpminone)
              else
                  call f_free_ptr(ovrlpminone)
              end if

              call f_free_ptr(ovrlppowerp)
              call f_free_ptr(ovrlppoweroldp)

          end if

          if(nproc > 1) then
              call mpi_allgatherv(inv_ovrlpp, ovrlp_smat%nfvctr*norbp, mpi_double_precision, inv_ovrlp_mat%matrix, &
                   ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                   mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
              call f_free_ptr(inv_ovrlpp)
          else
              nullify(inv_ovrlpp)
          end if
      end if

      if (check_accur) then
          call check_accur_overlap_minus_one(iproc,nproc,ovrlp_smat%nfvctr,norbp,isorb,power,&
               ovrlp_mat%matrix,inv_ovrlp_mat%matrix,error)
      end if
  else if (imode==SPARSE) then
      if (iorder==0) then
          !!if (iproc==0) call yaml_warning('The compressed matrix will not be filled! You should know what you do.')
          ovrlp_local = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='ovrlp_local')
          inv_ovrlp_local = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='inv_ovrlp_local')
          call uncompress_matrix(iproc, ovrlp_smat, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_local)
          !!write(*,*) ovrlp_smat%matrix_compr
          !!write(*,*) '==============='
          !!write(*,*) ovrlp_local
          call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr,ovrlp_local(1,1),1,inv_ovrlp_local(1,1),1)
          !!do iorb=1,orbs%ovrlp_smat%nfvctr
          !!    do jorb=1,orbs%ovrlp_smat%nfvctr
          !!        inv_ovrlp_local(jorb,iorb)=0.5d0*(ovrlp_local(jorb,iorb)+ovrlp_local(iorb,jorb))
          !!    end do
          !!end do
          if (power==1) then
             if (blocksize<0) then
                call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_local)
             else
                stop 'check if working - upper half may not be filled'
                call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                     ovrlp_smat%nfvctr, inv_ovrlp_local(1,1), ovrlp_smat%nfvctr)
                call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                     ovrlp_smat%nfvctr, inv_ovrlp_local(1,1), ovrlp_smat%nfvctr)
             end if
          else if (power==2) then
              !if (nproc>1) then
                  call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                       blocksize,.true.,inv_ovrlp_local,inv_ovrlp_smat)
              !else
              !    call overlap_plus_minus_one_half_exact(ovrlp_smat%nfvctr,blocksize,.true.,inv_ovrlp_local)
              !end if
          else if (power==-2) then
              !if (nproc>1) then
                  call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                       blocksize,.false.,inv_ovrlp_local,inv_ovrlp_smat)
              !else
              !    call overlap_plus_minus_one_half_exact(ovrlp_smat%nfvctr,blocksize,.false.,inv_ovrlp_local)
              !end if
          end if
          !!! These two lines can be deleted as soon as the tests are stabilized ##########
          !!inv_ovrlp_smat%matrix=inv_ovrlp_mat%matrix
          call compress_matrix(iproc, inv_ovrlp_smat, inmat=inv_ovrlp_local, outmat=inv_ovrlp_mat%matrix_compr)
          call f_free_ptr(ovrlp_local)
          call f_free_ptr(inv_ovrlp_local)
          ! #############################################################################
      else if (iorder<0) then ! could improve timing for checking, but for now just making sure it works
          ! use 4 submatrices
          if (nproc>0) then
              Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat12p')
              Amat21p = sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat21p')
              Amat11p = sparsematrix_malloc0_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat11p')
          else
              Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat12p')
              Amat21p = sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat21p')
              Amat11p = sparsematrix_malloc0_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat11p')
          end if
          ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
          Amat22p=>Amat11p
          Amat12_compr=>inv_ovrlp_mat%matrix_compr

          call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
               ovrlp_mat%matrix_compr, Amat12_compr, 'small_to_large')
          Amat12_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='Amat12_seq')
          call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat12_compr, Amat12_seq)
          call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, Amat12_compr, Amat12p)

          do iorb=1,norbp
              Amat21p(iorb+isorb,iorb)=1.0d0
          end do
          Amat21_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='Amat21_compr')
          call compress_matrix_distributed(iproc, inv_ovrlp_smat, Amat21p, Amat21_compr)
          Amat21_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='Amat21_seq')
          call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat21_compr, Amat21_seq)

          ! calculate Xn+1=0.5*Xn*(3I-Xn**2)
          do its=1,abs(iorder)
              call sparsemm(inv_ovrlp_smat, Amat12_seq, Amat21p, Amat11p)

              if (norbp>0) call vscal(ovrlp_smat%nfvctr*norbp,-0.5d0,Amat11p(1,1),1)
              !call vscal(ovrlp_smat%nfvctr*norbp,-0.5d0,Amat22p(1,1),1)
              do iorb=1,norbp
                  Amat11p(iorb+isorb,iorb)=Amat11p(iorb+isorb,iorb)+1.5d0
              !    Amat22p(iorb+isorb,iorb)=Amat22p(iorb+isorb,iorb)+1.5d0
              end do

              call sparsemm(inv_ovrlp_smat, Amat12_seq, Amat22p, Amat12p)
              call sparsemm(inv_ovrlp_smat, Amat21_seq, Amat11p, Amat21p)

              if (its/=abs(iorder).or.power/=2) then
                  call compress_matrix_distributed(iproc, inv_ovrlp_smat, Amat21p, Amat21_compr)
              end if
              if (its/=abs(iorder).or.power==1) then
                  call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat21_compr, Amat21_seq)
              end if
              if (its/=abs(iorder).or.power==2) then
                  call compress_matrix_distributed(iproc, inv_ovrlp_smat, Amat12p, Amat12_compr)
              end if
              if (its/=abs(iorder)) then
                  call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat12_compr, Amat12_seq)
              end if
          end do

          call f_free(Amat12_seq)
          nullify(Amat22p)
          call f_free_ptr(Amat11p)

          if (power==1) then
              call sparsemm(inv_ovrlp_smat, Amat21_seq, Amat21p, Amat12p)
              call compress_matrix_distributed(iproc, inv_ovrlp_smat, Amat12p, inv_ovrlp_mat%matrix_compr)
          !else if (power==2) then
          !    call vcopy(inv_ovrlp_smat%nvctr,Amat12_compr(1),1,inv_ovrlp_smat%matrix_compr(1),1)
          else if (power==-2) then
              call vcopy(inv_ovrlp_smat%nvctr,Amat21_compr(1),1,inv_ovrlp_mat%matrix_compr(1),1)
          end if

          nullify(Amat12_compr)
          call f_free(Amat21_compr)
          call f_free(Amat12p)
          call f_free(Amat21p)
          call f_free(Amat21_seq)

      else
          ovrlp_large_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='ovrlp_large_compr')
          call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
               ovrlp_mat%matrix_compr, ovrlp_large_compr, 'small_to_large')

          if (iorder>1) then
              ovrlpminone_sparse_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlpminone_sparse_seq')
              ovrlpminone_sparse = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=SPARSE_FULL, id='ovrlpminone_sparse')
              ovrlpminoneoldp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlpminoneoldp')

              call matrix_minus_identity_sparse(ovrlp_smat%nfvctr, inv_ovrlp_smat, ovrlp_large_compr, ovrlpminone_sparse)
              call sequential_acces_matrix_fast(inv_ovrlp_smat, ovrlpminone_sparse, ovrlpminone_sparse_seq)
              call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, ovrlpminone_sparse, ovrlpminoneoldp)

              call f_free_ptr(ovrlpminone_sparse)

              if (power==1) then
                  factor=-1.0d0
              else if (power==2) then
                  factor=0.5d0
              else if (power==-2) then
                  factor=-0.5d0
              end if
          end if

          ovrlpminonep = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlpminonep')
          invovrlpp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='invovrlpp')

          if (norbp>0) then
              call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, ovrlp_large_compr, ovrlpminonep)
              if (.not.check_accur) call f_free(ovrlp_large_compr)
              call first_order_taylor_dense(ovrlp_smat%nfvctr,isorb,norbp,power,ovrlpminonep,invovrlpp)
          end if

          do i=2,iorder
               call sparsemm(inv_ovrlp_smat, ovrlpminone_sparse_seq, ovrlpminoneoldp, ovrlpminonep)
              factor=newfactor(power,i,factor)
              call daxpy(ovrlp_smat%nfvctr*norbp,factor,ovrlpminonep,1,invovrlpp,1)
              if (i/=iorder.and.norbp>0) call vcopy(ovrlp_smat%nfvctr*norbp,ovrlpminonep(1,1),1,ovrlpminoneoldp(1,1),1)
          end do
          !!call to_zero(inv_ovrlp_smat%nvctr, inv_ovrlp_smat%matrix_compr(1))
            call compress_matrix_distributed(iproc, inv_ovrlp_smat, invovrlpp, inv_ovrlp_mat%matrix_compr)

          if (iorder>1) then
              call f_free(ovrlpminone_sparse_seq)
              call f_free(ovrlpminoneoldp)
              !!if (.not.check_accur) call f_free(istindexarr)
              !!if (.not.check_accur) call f_free(ivectorindex)
              !!if (.not.check_accur) call f_free_ptr(onedimindices)
          end if

          if (.not.check_accur) call f_free(invovrlpp)
          call f_free_ptr(ovrlpminonep)
      end if

      if (check_accur) then
          ! HERE STARTS LINEAR CHECK ##########################
          if (iorder<1) then
              invovrlpp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='invovrlpp')
              ovrlp_large_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='ovrlp_large_compr')
              call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
                   ovrlp_mat%matrix_compr, ovrlp_large_compr, 'small_to_large')
          end if
          invovrlp_compr_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlp_large_compr_seq')
          ovrlp_largep = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlp_largep')
          call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, ovrlp_large_compr, ovrlp_largep)
          call sequential_acces_matrix_fast(inv_ovrlp_smat, inv_ovrlp_mat%matrix_compr, invovrlp_compr_seq)

          if (power==1) then
              call check_accur_overlap_minus_one_sparse(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, norbp, isorb, &
                   inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                   inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                   invovrlp_compr_seq, ovrlp_largep, power, &
                   error)
          else if (power==2) then
              call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, inv_ovrlp_mat%matrix_compr, invovrlpp)
              call check_accur_overlap_minus_one_sparse(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, norbp, isorb, &
                   inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                   inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                   invovrlp_compr_seq, invovrlpp, power, &
                   error,cmatp=ovrlp_largep)
          else if (power==-2) then
              ovrlp_compr_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlp_compr_seq') 
              call sequential_acces_matrix_fast(inv_ovrlp_smat, ovrlp_large_compr, ovrlp_compr_seq)
              call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, inv_ovrlp_mat%matrix_compr, invovrlpp)
              call check_accur_overlap_minus_one_sparse(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, norbp, isorb, &
                    inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                    inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                    invovrlp_compr_seq, invovrlpp, power, &
                    error, &
                    ovrlp_compr_seq)
              call f_free(ovrlp_compr_seq)
          else
              stop 'wrong power'
          end if
          call f_free(invovrlp_compr_seq)
          call f_free(ovrlp_largep)
          call f_free(invovrlpp)
          call f_free(ovrlp_large_compr)
          !HERE ENDS LINEAR CHECK #############################
      end if
  end if sparse_dense


  call timing(iproc,'lovrlp^-1     ','OF')
  call f_release_routine()


end subroutine overlapPowerGeneral


pure function newfactor(power,order,factor)
  implicit none
  integer, intent(in) :: power, order
  real(kind=8), intent(in) :: factor
  real(kind=8) :: newfactor

  if (power==1) then
     newfactor=-factor
  else if (power==2) then
     newfactor=factor*(1.5d0-real(order,kind=8))/real(order,kind=8)
  else if (power==-2) then
     newfactor=factor*(0.5d0-real(order,kind=8))/real(order,kind=8)
  end if

end function newfactor


subroutine matrix_minus_identity_dense(norb,isorb,norbp,matinp,matoutp)
  implicit none
  integer,intent(in) :: norb, isorb, norbp
  real(kind=8),dimension(norb,norbp),intent(in) :: matinp
  real(kind=8),dimension(norb,norbp),intent(out) :: matoutp

  integer :: iorb, jorb, iiorb

  !$omp parallel do default(private) shared(matinp,matoutp,norb,isorb,norbp)
  do iorb=1,norbp
     iiorb=isorb+iorb
     do jorb=1,norb
        if(iiorb==jorb) then
           matoutp(jorb,iorb) = matinp(jorb,iorb) - 1.0d0
        else
           matoutp(jorb,iorb) = matinp(jorb,iorb)
        end if
     end do
  end do
  !$omp end parallel do

end subroutine matrix_minus_identity_dense

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

subroutine first_order_taylor_dense(norb,isorb,norbp,power,ovrlpp,inv_ovrlpp)
use module_base
  implicit none
  integer,intent(in) :: norb, isorb, norbp, power
  real(kind=8),dimension(norb,norbp),intent(in) :: ovrlpp
  real(kind=8),dimension(norb,norbp),intent(out) :: inv_ovrlpp

  integer :: iorb, jorb, iiorb

  if (power==1) then
     !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
     do iorb=1,norbp
        iiorb=isorb+iorb
        do jorb=1,norb
           if(iiorb==jorb) then
              inv_ovrlpp(jorb,iorb) = 2.d0 - ovrlpp(jorb,iorb)
           else
              inv_ovrlpp(jorb,iorb) = -ovrlpp(jorb,iorb)
           end if
        end do
     end do
     !$omp end parallel do
  else if (power==2) then
     !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
     do iorb=1,norbp
        iiorb=isorb+iorb
        do jorb=1,norb
           if(iiorb==jorb) then
              inv_ovrlpp(jorb,iorb) = 0.5d0 + 0.5d0*ovrlpp(jorb,iorb)
           else
              inv_ovrlpp(jorb,iorb) = 0.5d0*ovrlpp(jorb,iorb)
           end if
        end do
     end do
     !$omp end parallel do
  else if (power==-2) then
     !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
     do iorb=1,norbp
        iiorb=isorb+iorb
        do jorb=1,norb
           if(iiorb==jorb) then
              inv_ovrlpp(jorb,iorb) = 1.5d0 - 0.5d0*ovrlpp(jorb,iorb)
           else
              inv_ovrlpp(jorb,iorb) = -0.5d0*ovrlpp(jorb,iorb)
           end if
        end do
     end do
     !$omp end parallel do
  else
     stop 'Error in first_order_taylor_dense'
  end if

end subroutine first_order_taylor_dense


subroutine overlap_minus_one_exact_serial(norb,inv_ovrlp)
  use module_base
  implicit none
  integer,intent(in) :: norb
  real(kind=8),dimension(norb,norb),intent(inout) :: inv_ovrlp

  integer :: info, iorb, jorb

  call dpotrf('u', norb, inv_ovrlp(1,1), norb, info)
  if(info/=0) then
     write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
     stop
  end if
  call dpotri('u', norb, inv_ovrlp(1,1), norb, info)
  if(info/=0) then
     write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
     stop
  end if

  ! fill lower half
  !$omp parallel do default(private) shared(inv_ovrlp,norb)
  do iorb=1,norb
     do jorb=1,iorb-1
        inv_ovrlp(iorb,jorb)=inv_ovrlp(jorb,iorb)
     end do
  end do
  !$omp end parallel do 

end subroutine overlap_minus_one_exact_serial


subroutine overlap_plus_minus_one_half_exact(nproc,norb,blocksize,plusminus,inv_ovrlp_half,smat)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  implicit none
  integer,intent(in) :: nproc,norb,blocksize
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_half
  logical, intent(in) :: plusminus
  type(sparse_matrix),intent(in),optional :: smat

  integer :: info, iorb, jorb, ierr, iiorb, isorb, norbp, lwork, jjorb
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:),allocatable :: tempArr, orig_ovrlp
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_halfp
  real(kind=8),dimension(:,:), allocatable :: vr,vl ! for non-symmetric LAPACK
  real(kind=8),dimension(:),allocatable:: eval1 ! for non-symmetric LAPACK
  real(dp) :: temp, error
  real(dp), allocatable, dimension(:) :: temp_vec
  logical, parameter :: symmetric=.true.
  logical, parameter :: check_lapack=.true.

  if (nproc>1) then
      if (.not.present(smat)) then 
          call f_err_throw('overlap_plus_minus_one_half_exact: for nproc>1, smat must be present!', &
               err_name='BIGDFT_RUNTIME_ERROR')
      end if
  end if
           

  eval=f_malloc(norb,id='eval')
  if(blocksize>0) then
     call dsyev_parallel(bigdft_mpi%iproc, bigdft_mpi%nproc, min(blocksize,norb), bigdft_mpi%mpi_comm, 'v', 'l', norb, &
          inv_ovrlp_half(1,1), norb, eval(1), info)
     if(info/=0) then
        write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
     end if
  else
     if (symmetric) then
        if (check_lapack) then
           orig_ovrlp=f_malloc((/norb,norb/),id='orig_ovrlp')
           call vcopy(norb*norb,inv_ovrlp_half(1,1),1,orig_ovrlp(1,1),1)
        end if
        work=f_malloc(1000,id='work')
        call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, -1, info)
        lwork = int(work(1))
        call f_free(work)
        work=f_malloc(lwork,id='work')
           !!do iorb=1,norb
           !!   do jorb=1,norb
           !!      write(2000+bigdft_mpi%iproc,'(a,3i8,es16.7)') 'iproc, iorb, jorb, val', bigdft_mpi%iproc, iorb, jorb, inv_ovrlp_half(jorb,iorb)
           !!   end do
           !!end do
        call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, lwork, info)
        if (check_lapack) then
           tempArr=f_malloc((/norb,norb/), id='tempArr')
           do iorb=1,norb
              do jorb=1,norb
                 tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
              end do
           end do
           inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
           call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
           call f_free(tempArr)
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, error)
           if (bigdft_mpi%iproc==0.and.abs(error)>1.0d-8) then
              print*,'LAPACK error for dsyev in overlap_plus_minus_one_half_exact',error
              open(99,file='dsyev_input.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,orig_ovrlp(iorb,jorb)
                 end do
              end do
              close(99)
              open(99,file='dsyev_output.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,inv_ovrlp_halfp(iorb,jorb)
                 end do
              end do
              close(99)
              call mpi_finalize(bigdft_mpi%mpi_comm)
              stop
           end if
           call f_free_ptr(inv_ovrlp_halfp)
           call f_free(orig_ovrlp)
        end if
     else
        if (check_lapack) then
           orig_ovrlp=f_malloc((/norb,norb/),id='orig_ovrlp')
           call vcopy(norb*norb,inv_ovrlp_half(1,1),1,orig_ovrlp(1,1),1)
        end if
        work=f_malloc(1000,id='work')
        eval1=f_malloc(norb,id='eval1')
        vl=f_malloc((/norb,norb/),id='vl')
        vr=f_malloc((/norb,norb/),id='vr')
        call dgeev( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
             norb, WORK, -1, info )
        lwork = nint(work(1))
        call f_free(work)
        work=f_malloc(lwork,id='work')
        call DGEEV( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
             norb, WORK, LWORK, info )
        call vcopy(norb*norb,vl(1,1),1,inv_ovrlp_half(1,1),1)
        call f_free(eval1)
        call f_free(vr)
        call f_free(vl)
        temp_vec=f_malloc(norb,id='temp_vec')
        do iorb=1,norb
           do jorb=iorb+1,norb
              if (eval(jorb) < eval(iorb)) then
                 temp = eval(iorb)
                 temp_vec = inv_ovrlp_half(:,iorb)
                 eval(iorb) = eval(jorb)
                 eval(jorb) = temp
                 inv_ovrlp_half(:,iorb) = inv_ovrlp_half(:,jorb)
                 inv_ovrlp_half(:,jorb) = temp_vec
              end if
           end do
        end do
        call f_free(temp_vec)
        if (check_lapack) then
           tempArr=f_malloc((/norb,norb/), id='tempArr')
           do iorb=1,norb
              do jorb=1,norb
                 tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
              end do
           end do
           inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
           call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
           call f_free(tempArr)
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, error)
           if (bigdft_mpi%iproc==0.and.abs(error)>1.0d-8) then
              print*,'LAPACK error for dgeev in overlap_plus_minus_one_half_exact',error
              open(99,file='dgeev_input.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,orig_ovrlp(iorb,jorb)
                 end do
              end do
              close(99)
              open(99,file='dgeev_output.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,inv_ovrlp_halfp(iorb,jorb)
                 end do
              end do
              close(99)
              call mpi_finalize(bigdft_mpi%mpi_comm)
              stop
           end if
           call f_free_ptr(inv_ovrlp_halfp)
           call f_free(orig_ovrlp)
        end if
     end if

     if(info/=0) then
        write(*,'(a,i0,2x,i0)') 'ERROR in dsyev (overlap_plus_minus_one_half_exact), info, norb=', info, norb
        stop
     end if
     call f_free(work)
  end if

  ! Calculate S^{-1/2}. 
  ! First calculate ovrlp*diag(1/sqrt(eval)) (ovrlp is the diagonalized overlap
  ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
  !write(*,*) 'iproc, present(orbs), norb, orbs%norb', bigdft_mpi%iproc, present(orbs), norb, orbs%norb
  !if (present(orbs).and.bigdft_mpi%nproc>1.and.blocksize<0) then
     !if (norb/=orbs%norb) stop 'Error with orbs%norb in overlap_plus_minus_one_half_exact'
     if (nproc>1) then
         norbp=smat%nfvctrp
         isorb=smat%isfvctr
     else
         norbp=norb
         isorb=0
     end if
  !else
  !   norbp=norb
  !   isorb=0
  !end if
  tempArr=f_malloc((/norbp,norb/), id='tempArr')
  !$omp parallel do default(private) shared(tempArr,inv_ovrlp_half,eval,plusminus,norb,norbp,isorb)
  do iorb=1,norb
     do jorb=1,norbp
        jjorb=jorb+isorb
        if (plusminus) then
           tempArr(jorb,iorb)=inv_ovrlp_half(jjorb,iorb)*sqrt(abs(eval(iorb)))
        else
           tempArr(jorb,iorb)=inv_ovrlp_half(jjorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
        end if
     end do
  end do
  !$omp end parallel do


  call f_free(eval)
  inv_ovrlp_halfp=f_malloc_ptr((/norb,norbp/), id='inv_ovrlp_halfp')
  ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
  ! This will give S^{+/-1/2}.
  if(blocksize<0) then
     if (norbp>0) call dgemm('n', 't', norb, norbp, norb, 1.d0, inv_ovrlp_half, &
          norb, tempArr, norbp, 0.d0, inv_ovrlp_halfp, norb)
     !if (present(orbs).and.bigdft_mpi%nproc>1) then
     if (nproc>1) then
        call mpi_allgatherv(inv_ovrlp_halfp, norb*norbp, mpi_double_precision, inv_ovrlp_half, &
                   norb*smat%nfvctr_par(:), norb*smat%isfvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call vcopy(norb*norbp,inv_ovrlp_halfp(1,1),1,inv_ovrlp_half(1,1),1)
     end if
  else
     call dgemm_parallel(bigdft_mpi%iproc, bigdft_mpi%nproc, blocksize, bigdft_mpi%mpi_comm, 'n', 't', norb, norb, norb, &
          1.d0, inv_ovrlp_half, norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
     call vcopy(norb*norbp,inv_ovrlp_halfp(1,1),1,inv_ovrlp_half(1,1),1)
  end if

  call f_free_ptr(inv_ovrlp_halfp)
  call f_free(tempArr)

end subroutine overlap_plus_minus_one_half_exact



subroutine check_accur_overlap_minus_one_sparse(iproc, nproc, smat, norb, norbp, isorb, nseq, nout, &
           ivectorindex, onedimindices, amat_seq, bmatp, power, &
           error, &
           dmat_seq, cmatp)
  use module_base
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only: sparsemm
  implicit none
  integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseq, nout, power
  type(sparse_matrix) :: smat
  integer,dimension(nseq),intent(in) :: ivectorindex
  integer,dimension(4,nout) :: onedimindices
  real(kind=8),dimension(nseq),intent(in) :: amat_seq
  real(kind=8),dimension(norb,norbp),intent(in) :: bmatp
  real(kind=8),intent(out) :: error
  real(kind=8),dimension(nseq),intent(in),optional :: dmat_seq
  real(kind=8),dimension(norb,norbp),intent(in),optional :: cmatp

  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
  real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
  integer :: ierr, i,j

  tmpp=f_malloc0((/norb,norbp/),id='tmpp')
  if (power==1) then
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(smat, amat_seq, bmatp, tmpp)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, error)
  else if (power==2) then
      if (.not.present(cmatp)) stop 'cmatp not present'
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(smat, amat_seq, bmatp, tmpp)
     call max_matrix_diff_parallel(iproc, norb, norbp, tmpp, cmatp, error)
     error=0.5d0*error
  else if (power==-2) then
     if (.not.present(dmat_seq)) stop 'dmat_seq not present'
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(smat, amat_seq, bmatp, tmpp)
     tmp2p=f_malloc0((/norb,norbp/),id='tmp2p')
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
     !!     norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
     call sparsemm(smat, dmat_seq, tmpp, tmp2p)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, error)
     error=0.5d0*error
     call f_free(tmp2p)
  else
     stop 'Error in check_accur_overlap_minus_one_sparse'
  end if
  call f_free(tmpp)

end subroutine check_accur_overlap_minus_one_sparse




subroutine check_accur_overlap_minus_one(iproc,nproc,norb,norbp,isorb,power,ovrlp,inv_ovrlp,error)
  use module_base
  implicit none
  integer,intent(in) :: iproc, nproc, norb, norbp, isorb, power
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp, inv_ovrlp
  real(kind=8),intent(out) :: error

  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
  real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
  integer :: ierr, i,j

  tmpp=f_malloc((/norb,norbp/),id='tmpp')
  if (power==1) then
     if (norbp>0) then
        call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
             norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     end if
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, error)
  else if (power==2) then
     if (norbp>0) then
        call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
             norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     end if
     call max_matrix_diff_parallel(iproc, norb, norbp, tmpp, ovrlp(1,isorb+1), error)
     error=0.5d0*error
  else if (power==-2) then
     if (norbp>0) then
        call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
             norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     end if
     tmp2p=f_malloc((/norb,norbp/),id='tmp2p')
     if (norbp>0) then
        call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
             norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
     end if
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, error)
     error=0.5d0*error
     call f_free(tmp2p)
  else
     stop 'Error in check_accur_overlap_minus_one'
  end if
  call f_free(tmpp)

end subroutine check_accur_overlap_minus_one


subroutine max_matrix_diff(iproc, norb, mat1, mat2, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: mat1, mat2
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, jorb
  real(8):: error

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norb
     do jorb=1,norb
        error=abs(mat1(jorb,iorb)-mat2(jorb,iorb))
        deviation=max(error,deviation)
     end do
  end do
  call timing(iproc,'dev_from_unity','OF') 

end subroutine max_matrix_diff


subroutine max_matrix_diff_parallel(iproc, norb, norbp, mat1, mat2, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb, norbp
  real(8),dimension(norb,norbp),intent(in):: mat1, mat2
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, jorb, ierr
  real(8):: error

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norbp
     !$omp parallel default(private) shared(iorb, norb, mat1, mat2, deviation)
     !$omp do reduction(max:deviation)
     do jorb=1,norb
        error=abs(mat1(jorb,iorb)-mat2(jorb,iorb))
        deviation=max(error,deviation)
     end do
     !$omp end do
     !$omp end parallel
  end do

  if (bigdft_mpi%nproc > 1) then
      call mpiallred(deviation, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if

  call timing(iproc,'dev_from_unity','OF') 

end subroutine max_matrix_diff_parallel



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

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norb
     do jorb=1,norb
        if(iorb==jorb) then
           error=abs(ovrlp(jorb,iorb)-1.d0)
        else
           error=abs(ovrlp(jorb,iorb))
        end if
        deviation=max(error,deviation)
     end do
  end do
  call timing(iproc,'dev_from_unity','OF') 

end subroutine deviation_from_unity


subroutine deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, ovrlp, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb, norbp, isorb
  real(8),dimension(norb,norbp),intent(in):: ovrlp
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, iiorb, jorb, ierr
  real(8):: error


  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norbp
     iiorb=iorb+isorb
     !$omp parallel default(private) shared(norb, iiorb, ovrlp, iorb, deviation)
     !$omp do reduction(max:deviation)
     do jorb=1,norb
        if(iiorb==jorb) then
           error=abs(ovrlp(jorb,iorb)-1.d0)
        else
           error=abs(ovrlp(jorb,iorb))
        end if
        deviation=max(error,deviation)
     end do
     !$omp end do
     !$omp end parallel
  end do
  if (nproc>1) then
      call mpiallred(deviation, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if
  call timing(iproc,'dev_from_unity','OF') 

end subroutine deviation_from_unity_parallel


subroutine deviation_from_symmetry(iproc, norb, ovrlp, deviation)
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

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norb
     do jorb=iorb+1,norb
        error=abs(ovrlp(jorb,iorb)-ovrlp(iorb,jorb))
        deviation=max(error,deviation)
     end do
  end do
  call timing(iproc,'dev_from_unity','OF') 

end subroutine deviation_from_symmetry

subroutine overlap_power_minus_one_half_parallel(iproc, nproc, meth_overlap, orbs, ovrlp, ovrlp_mat, &
           inv_ovrlp_half, inv_ovrlp_half_)
  use module_base
  use module_types
  use module_interfaces
  use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, &
                               allocate_matrices, deallocate_matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, meth_overlap
  type(orbitals_data),intent(in) :: orbs
  type(sparse_matrix),intent(inout) :: ovrlp
  type(matrices),intent(inout) :: ovrlp_mat
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half
  type(matrices),intent(inout) :: inv_ovrlp_half_

  ! Local variables
  integer :: iend, i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
  integer :: iiorb, ierr, ii, iseg, ind
  real(kind=8) :: error
  real(kind=8),dimension(:,:),pointer :: ovrlp_tmp, ovrlp_tmp_inv_half
  logical,dimension(:),allocatable :: in_neighborhood
  character(len=*),parameter :: subname='overlap_power_minus_one_half_parallel'
  !type(matrices) :: inv_ovrlp_half_

  call timing(iproc,'lovrlp^-1/2par','ON')

  in_neighborhood = f_malloc(orbs%norb,id='in_neighborhood')

  !inv_ovrlp_half_ = matrices_null()
  !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_)
  call to_zero(inv_ovrlp_half%nvctr, inv_ovrlp_half_%matrix_compr(1))

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

     ovrlp_tmp = f_malloc0_ptr((/n,n/),id='ovrlp_tmp')

     jjorb=0
     do jorb=1,orbs%norb
        if (.not.in_neighborhood(jorb)) cycle
        jjorb=jjorb+1
        kkorb=0
        do korb=1,orbs%norb
           if (.not.in_neighborhood(korb)) cycle
           kkorb=kkorb+1
           ind = matrixindex_in_compressed(ovrlp,korb, jorb)
           if (ind>0) then
              ovrlp_tmp(kkorb,jjorb)=ovrlp_mat%matrix_compr(ind)
           else
              ovrlp_tmp(kkorb,jjorb)=0.d0
           end if
           !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
        end do
     end do
          
     ovrlp_tmp_inv_half = f_malloc_ptr((/n,n/),id='ovrlp_tmp_inv_half')
     call vcopy(n*n, ovrlp_tmp(1,1), 1, ovrlp_tmp_inv_half(1,1), 1)

     !if (iiorb==orbs%norb) then
     !print*,''
     !print*,'ovrlp_tmp',n,iiorb
     !do jorb=1,n
     !print*,jorb,ovrlp_tmp(:,jorb)
     !end do
     !end if


     ! Calculate S^-1/2 for the small overlap matrix
     !!call overlapPowerGeneral(iproc, nproc, meth_overlap, -2, -8, n, orbs, imode=2, check_accur=.true.,&
     !!     ovrlp=ovrlp_tmp, inv_ovrlp=ovrlp_tmp_inv_half, error=error)
     !!call overlapPowerGeneral(iproc, 1, meth_overlap, -2, -8, n, orbs, imode=2, &
     !!     ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
     !!     ovrlp_mat=ovrlp_mat, inv_ovrlp_mat=inv_ovrlp_half_, check_accur=.true., &
     !!     ovrlp=ovrlp_tmp, inv_ovrlp=ovrlp_tmp_inv_half, error=error)
     call overlap_plus_minus_one_half_exact(1, n, -8, .false., ovrlp_tmp_inv_half)


     !if (iiorb==orbs%norb) then
     !print*,''
     !print*,'inv_ovrlp_tmp',n,iiorb,error
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
              ind = matrixindex_in_compressed(inv_ovrlp_half,korb,jorb)
              if (ind>0) then
                 inv_ovrlp_half_%matrix_compr(ind)=ovrlp_tmp_inv_half(kkorb,jjorb)
                 !if (iiorb==orbs%norb) print*,'problem here?!',iiorb,kkorb,jjorb,korb,jorb,ind,ovrlp_tmp_inv_half(kkorb,jjorb)
              end if
              !write(1300+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
           end do
           exit !no need to keep looping
        end if
     end do


     call f_free_ptr(ovrlp_tmp_inv_half)
     call f_free_ptr(ovrlp_tmp)

  end do

  if (nproc>1)then
      call mpiallred(inv_ovrlp_half_%matrix_compr(1), inv_ovrlp_half%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  call f_free(in_neighborhood)

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

  !call deallocate_matrices(inv_ovrlp_half_)

  call timing(iproc,'lovrlp^-1/2par','OF')

end subroutine overlap_power_minus_one_half_parallel

!> Orthonormalize a subset of orbitals
subroutine orthonormalize_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalize_subset
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparse_matrix),intent(inout) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall, iorb, jorb, iat, jat, ii
  logical :: iout, jout
  integer,dimension(:),allocatable :: icount_norb, jcount_norb
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparse_matrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='orthonormalize_subset'
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: error
  type(matrices) :: ovrlp_, inv_ovrlp_half_


  !call nullify_sparse_matrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
  inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')

  inv_ovrlp_half_ = matrices_null()
  call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_)


  if(.not.can_use_transposed) then
      if(associated(psit_c)) then
          call f_free_ptr(psit_c)
      end if
      if(associated(psit_f)) then
          call f_free_ptr(psit_f)
      end if
      psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
      psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')

      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.

  end if

  ovrlp_ = matrices_null()
  call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
  ! This can then be deleted if the transition to the new type has been completed.

  ! For the "higher" TMBs: delete off-diagonal elements and
  ! set diagonal elements to 1
  icount_norb = f_malloc(at%astruct%nat,id='icount_norb')
  jcount_norb = f_malloc(at%astruct%nat,id='jcount_norb')
  icount_norb=0
  do iorb=1,orbs%norb
      iat=orbs%onwhichatom(iorb)
      icount_norb(iat)=icount_norb(iat)+1
      if (icount_norb(iat)<minorbs_type(at%astruct%iatype(iat)) .or. &
          icount_norb(iat)>maxorbs_type(at%astruct%iatype(iat))) then
          iout=.true.
      else
          iout=.false.
      end if
      jcount_norb=0
      do jorb=1,orbs%norb
          jat=orbs%onwhichatom(jorb)
          jcount_norb(jat)=jcount_norb(jat)+1
          if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat)) .or. &
              jcount_norb(jat)>maxorbs_type(at%astruct%iatype(jat))) then
              jout=.true.
          else
              jout=.false.
          end if
          ii=matrixindex_in_compressed(ovrlp,jorb,iorb)
          !!if (iproc==0) write(444,'(a,2i7,2x,2i7,3x,2l4,3x,3i6)') 'iorb, jorb, iat, jat, iout, jout, icount_norb(iat), minorbs_type(at%iatype(iat)), maxorbs_type(at%iatype(iat))', &
          !!                                                         iorb, jorb, iat, jat, iout, jout, icount_norb(iat), minorbs_type(at%iatype(iat)), maxorbs_type(at%iatype(iat))
          if (ii/=0 .and. (iout .or. jout)) then
              if (jorb==iorb) then
                  ovrlp_%matrix_compr(ii)=1.d0
              else
                  ovrlp_%matrix_compr(ii)=0.d0
              end if
          end if
      end do
  end do
  call f_free(icount_norb)
  call f_free(jcount_norb)


  if (methTransformOverlap==-1) then
      !ovrlp_ = matrices_null()
      !call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
      !ovrlp_%matrix_compr=ovrlp%matrix_compr
      call overlap_power_minus_one_half_parallel(iproc, nproc, methTransformOverlap, &
           orbs, ovrlp, ovrlp_, inv_ovrlp_half, inv_ovrlp_half_)
      !call deallocate_matrices(ovrlp_)
  else
      nullify(inv_ovrlp_null)
      ! do sparse.. check later
      !ovrlp%matrix_compr=ovrlp_%matrix_compr
      call overlapPowerGeneral(iproc, nproc, methTransformOverlap, -2, &
           orthpar%blocksize_pdsyev, &
           imode=1, check_accur=.true., &
           ovrlp_mat=ovrlp_, inv_ovrlp_mat=inv_ovrlp_half_, &
           error=error, ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half)
  end if

  ! For the "higher" TMBs: delete off-diagonal elements and
  ! set diagonal elements to 1
  icount_norb = f_malloc(at%astruct%nat,id='icount_norb')
  jcount_norb = f_malloc(at%astruct%nat,id='jcount_norb')
  icount_norb=0
  do iorb=1,orbs%norb
      iat=orbs%onwhichatom(iorb)
      icount_norb(iat)=icount_norb(iat)+1
      if (icount_norb(iat)<minorbs_type(at%astruct%iatype(iat)) .or. &
          icount_norb(iat)>maxorbs_type(at%astruct%iatype(iat))) then
          iout=.true.
      else
          iout=.false.
      end if
      jcount_norb=0
      do jorb=1,orbs%norb
          jat=orbs%onwhichatom(jorb)
          jcount_norb(jat)=jcount_norb(jat)+1
          if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat)) .or. &
              jcount_norb(jat)>maxorbs_type(at%astruct%iatype(jat))) then
              jout=.true.
          else
              jout=.false.
          end if
          ii=matrixindex_in_compressed(ovrlp,jorb,iorb)
          if (ii/=0 .and. (iout .or. jout)) then
              if (jorb==iorb) then
                  ovrlp_%matrix_compr(ii)=1.d0
              else
                  ovrlp_%matrix_compr(ii)=0.d0
              end if
          end if
      end do
  end do
  call f_free(icount_norb)
  call f_free(jcount_norb)

  psittemp_c = f_malloc(sum(collcom%nrecvcounts_c),id='psittemp_c')
  psittemp_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psittemp_f')

  call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
  call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)

  !inv_ovrlp_half_%matrix_compr = inv_ovrlp_half%matrix_compr
  call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_, &
       psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)



  call deallocate_matrices(ovrlp_)

  call deallocate_matrices(inv_ovrlp_half_)


  norm = f_malloc(orbs%norb,id='norm')
  call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)
  call f_free(norm)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

  call f_free(psittemp_c)
  call f_free(psittemp_f)
  call f_free_ptr(inv_ovrlp_half%matrix_compr)

end subroutine orthonormalize_subset



subroutine gramschmidt_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => gramschmidt_subset
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparse_matrix),intent(inout) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall, iorb, jorb, iat, jat, ii
  logical :: iout, jout
  integer,dimension(:),allocatable :: icount_norb, jcount_norb
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparse_matrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='gramschmidt_subset'
  type(matrices) :: ovrlp_


  !call nullify_sparse_matrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
  inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')


  if(.not.can_use_transposed) then
      if(associated(psit_c)) then
          call f_free_ptr(psit_c)
      end if
      if(associated(psit_f)) then
          call f_free_ptr(psit_f)
      end if
      psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
      psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')

      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.

  end if


  ovrlp_ = matrices_null()
  call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
  ! This can then be deleted if the transition to the new type has been completed.
  !ovrlp%matrix_compr=ovrlp_%matrix_compr

  ! For the "higher" TMBs: delete off-diagonal elements and
  ! set diagonal elements to 1
  icount_norb = f_malloc0(at%astruct%nat,id='icount_norb')
  jcount_norb = f_malloc(at%astruct%nat,id='jcount_norb')
  do iorb=1,orbs%norb
      iat=orbs%onwhichatom(iorb)
      icount_norb(iat)=icount_norb(iat)+1
      if (icount_norb(iat)<minorbs_type(at%astruct%iatype(iat)) .or. &
          icount_norb(iat)>maxorbs_type(at%astruct%iatype(iat))) then
          iout=.true.
      else
          iout=.false.
      end if
      call to_zero(at%astruct%nat,jcount_norb(1))
      do jorb=1,orbs%norb
          jat=orbs%onwhichatom(jorb)
          jcount_norb(jat)=jcount_norb(jat)+1
          !!if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat)) .or. &
          !!    jcount_norb(jat)>maxorbs_type(at%astruct%iatype(jat))) then
          if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat))) then
              jout=.true.
          else
              jout=.false.
          end if
          ii=matrixindex_in_compressed(ovrlp,jorb,iorb)
          if (ii/=0) then
              if (iout) then
                  ovrlp_%matrix_compr(ii)=0.d0
              else
                  if (jout) then
                      ovrlp_%matrix_compr(ii)=-ovrlp_%matrix_compr(ii)
                  else
                      ovrlp_%matrix_compr(ii)=0.d0
                  end if
              end if
          end if
          !!if (iout .or. jout) then
          !!    if (jorb==iorb) then
          !!        ovrlp%matrix_compr(ii)=1.d0
          !!    else
          !!        ovrlp%matrix_compr(ii)=0.d0
          !!    end if
          !!end if
      end do
  end do
  call f_free(icount_norb)
  call f_free(jcount_norb)


  !!if (methTransformOverlap==-1) then
  !!    call overlap_power_minus_one_half_parallel(iproc, nproc, methTransformOverlap, orbs, ovrlp, inv_ovrlp_half)
  !!else
  !!    call overlapPowerMinusOneHalf(iproc, nproc, bigdft_mpi%mpi_comm, methTransformOverlap, orthpar%blocksize_pdsyev, &
  !!        orthpar%blocksize_pdgemm, orbs%norb, ovrlp, inv_ovrlp_half)
  !!end if

  !!! For the "higher" TMBs: delete off-diagonal elements and
  !!! set diagonal elements to 1
  !!allocate(icount_norb(at%nat),stat=istat)
  !!call memocc(istat,icount_norb,'icount_norb',subname)
  !!allocate(jcount_norb(at%nat),stat=istat)
  !!call memocc(istat,jcount_norb,'jcount_norb',subname)
  !!do iorb=1,orbs%norb
  !!    iat=orbs%onwhichatom(iorb)
  !!    icount_norb(iat)=icount_norb(iat)+1
  !!    if (icount_norb(iat)<minorbs_type(at%iatype(iat)) .or. &
  !!        icount_norb(iat)>maxorbs_type(at%iatype(iat))) then
  !!        iout=.true.
  !!    else
  !!        iout=.false.
  !!    end if
  !!    do jorb=1,orbs%norb
  !!        jat=orbs%onwhichatom(jorb)
  !!        jcount_norb(jat)=jcount_norb(jat)+1
  !!        if (jcount_norb(jat)>maxorbs_type(at%iatype(jat))) then
  !!            jout=.true.
  !!        else
  !!            jout=.false.
  !!        end if
  !!        ii=ovrlp%matrixindex_in_compressed(jorb,iorb)
  !!        if (iout .or. jout) then
  !!            if (jorb==iorb) then
  !!                ovrlp%matrix_compr(ii)=1.d0
  !!            else
  !!                ovrlp%matrix_compr(ii)=0.d0
  !!            end if
  !!        end if
  !!    end do
  !!end do
  !!iall=-product(shape(icount_norb))*kind(icount_norb)
  !!deallocate(icount_norb, stat=istat)
  !!call memocc(istat, iall, 'icount_norb', subname)
  !!iall=-product(shape(jcount_norb))*kind(jcount_norb)
  !!deallocate(jcount_norb, stat=istat)
  !!call memocc(istat, iall, 'jcount_norb', subname)

  psittemp_c = f_malloc(sum(collcom%nrecvcounts_c),id='psittemp_c')
  psittemp_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psittemp_f')

  call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
  call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)
  !!call build_linear_combination_transposed(collcom, inv_ovrlp_half, &
  !!     psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)



  ovrlp_ = matrices_null()
  call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  ovrlp_%matrix_compr = inv_ovrlp_half%matrix_compr
  call build_linear_combination_transposed(collcom, ovrlp, ovrlp_, &
       psittemp_c, psittemp_f, .false., psit_c, psit_f, iproc)
  call deallocate_matrices(ovrlp_)


  call deallocate_matrices(ovrlp_)


  norm = f_malloc(orbs%norb,id='norm')
  !call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)

  call f_free(norm)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

  call f_free(psittemp_c)
  call f_free(psittemp_f)
  call f_free_ptr(inv_ovrlp_half%matrix_compr)

end subroutine gramschmidt_subset









   subroutine matrix_minus_identity_sparse(norb, smat, ovrlp_compr, ovrlpminone_compr)
     use module_base
     use module_types
     use sparsematrix_base, only: sparse_matrix
     implicit none

     ! Calling arguments
     integer,intent(in) :: norb
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nvctr),intent(in) :: ovrlp_compr
     real(kind=8),dimension(smat%nvctr),intent(out) :: ovrlpminone_compr

     ! Local variables
     integer :: ii, iseg, jorb, iiorb, jjorb

     !$omp parallel do default(private) &
     !$omp shared(smat, norb, ovrlpminone_compr, ovrlp_compr)
     do iseg=1,smat%nseg
         ii=smat%keyv(iseg)-1
         do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
             iiorb = (jorb-1)/norb + 1
             jjorb = jorb - (iiorb-1)*norb
             ii=ii+1
             if (iiorb==jjorb) then
                 ovrlpminone_compr(ii)=ovrlp_compr(ii)-1.d0
             else
                 ovrlpminone_compr(ii)=ovrlp_compr(ii)
             end if
         end do
     end do
     !$omp end parallel do

   end subroutine matrix_minus_identity_sparse



subroutine overlap_minus_one_half_serial(iproc, nproc, iorder, power, blocksize, &
           norb, ovrlp_matrix, inv_ovrlp_matrix, check_accur, &
           error)
  use module_base
  use module_types
  use module_interfaces, except_this_one => overlap_minus_one_half_serial
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, power, norb
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp_matrix
  real(kind=8),dimension(:,:),pointer,intent(out) :: inv_ovrlp_matrix
  logical,intent(in) :: check_accur
  real(kind=8),intent(out),optional :: error
  
  ! Local variables
  integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i, its, maxits
  integer :: matrixindex_in_compressed, nmaxvalk
  real(kind=8), dimension(:,:), pointer :: ovrlpminonep, ovrlpminone, inv_ovrlpp, ovrlppowerp, ovrlppoweroldp
  real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp, ovrlp_local, inv_ovrlp_local
  real(kind=8) :: factor, newfactor
  logical :: ovrlp_allocated, inv_ovrlp_allocated

  ! new for sparse taylor
  integer :: nout, nseq, nmaxsegk, nmaxval
  integer,dimension(:),allocatable :: ivectorindex
  integer,dimension(:,:),pointer :: onedimindices
  integer,dimension(:,:,:),allocatable :: istindexarr
  real(kind=8),dimension(:),pointer :: ovrlpminone_sparse
  real(kind=8),dimension(:),allocatable :: ovrlp_compr_seq, ovrlpminone_sparse_seq, ovrlp_large_compr
  real(kind=8),dimension(:),allocatable :: invovrlp_compr_seq
  real(kind=8),dimension(:,:),allocatable :: ovrlpminoneoldp, invovrlpp, ovrlp_largep
  real(kind=8),dimension(:,:),allocatable :: Amat12p, Amat21p, Amat21
  real(kind=8),dimension(:,:),pointer :: Amat12, Amat11p, Amat22p
  real(kind=8),dimension(:),pointer :: Amat12_compr
  real(kind=8),dimension(:),allocatable :: Amat21_compr, Amat12_seq, Amat21_seq
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2


  !!write(*,*) 'iorder',iorder


  call f_routine(id='overlapPowerGeneral')
  call timing(iproc,'lovrlp^-1     ','ON')

  if (nproc>1) then
      stop 'this routine only works in serial'
  end if

  
  if (check_accur) then
      if (.not.present(error)) stop 'error not present'
  end if

  if (power/=-2 .and. power/=1 .and. power/=2) stop 'wrong value of power'


      if (iorder==0) then
          call vcopy(norb*norb,ovrlp_matrix(1,1),1,inv_ovrlp_matrix(1,1),1)
          if (power==1) then
             call overlap_minus_one_exact_serial(norb,inv_ovrlp_matrix)
          else if (power==2) then
             call overlap_plus_minus_one_half_exact(1,norb,blocksize,.true.,inv_ovrlp_matrix)
          else if (power==-2) then
             call overlap_plus_minus_one_half_exact(1,norb,blocksize,.false.,inv_ovrlp_matrix)
          end if
      else if (iorder<0) then
          Amat12p = f_malloc((/norb,norb/), id='Amat12p')
          Amat21p = f_malloc((/norb,norb/), id='Amat21p')
          Amat11p = f_malloc_ptr((/norb,norb/), id='Amat11p')
          ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
          Amat22p=>Amat11p
          Amat12=>inv_ovrlp_matrix
          Amat21=f_malloc0((/norb,norb/), id='Amat21')

          call vcopy(norb*norb,ovrlp_matrix(1,1),1,Amat12(1,1),1)
          do iorb=1,norb
              Amat21(iorb,iorb)=1.0d0
          end do

          ! calculate Xn+1=0.5*Xn*(3I-Xn**2)
          do its=1,abs(iorder)
              call dgemm('n', 'n', norb, norb, norb, -0.5d0, Amat12(1,1), &
                   norb, Amat21(1,1), norb, 0.0d0, Amat11p(1,1), norb)
              do iorb=1,norb
                  Amat11p(iorb,iorb)=Amat11p(iorb,iorb)+1.5d0
              end do
              call dgemm('n', 'n', norb, norb, norb, 1.0d0, Amat12(1,1), &
                   norb, Amat22p(1,1), norb, 0.0d0, Amat12p(1,1), norb)
              call dgemm('n', 'n', norb, norb, norb, 1.0d0, Amat21(1,1), &
                   norb, Amat11p(1,1), norb, 0.0d0, Amat21p(1,1), norb)
              call vcopy(norb**2,Amat12p(1,1),1,Amat12(1,1),1)
              call vcopy(norb**2,Amat21p(1,1),1,Amat21(1,1),1)
          end do

          nullify(Amat22p)
          call f_free_ptr(Amat11p)

          if (power==1) then
              call dgemm('n', 'n', norb, norb, norb, 1.0d0, Amat21(1,1), &
                   norb, Amat21p(1,1), norb, 0.0d0, Amat12p(1,1), norb)
              call vcopy(norb**2, Amat12p(1,1), 1, inv_ovrlp_matrix(1,1), 1)
          else if (power==-2) then
              call vcopy(norb**2,Amat21(1,1),1,inv_ovrlp_matrix(1,1),1)
          end if

          call f_free(Amat12p)
          call f_free(Amat21p)
          nullify(Amat12)
          call f_free(Amat21)

      else
          if (iorder>1) then
              ovrlpminone = f_malloc_ptr((/norb,norb/), id='ovrlpminone')
              ovrlpminonep => ovrlpminone

              call matrix_minus_identity_dense(norb,0,norb,ovrlp_matrix(1,1),ovrlpminonep)

                  ovrlppoweroldp = f_malloc_ptr((/norb,norb/), id='ovrlppoweroldp')

              call vcopy(norb*norb,ovrlpminonep(1,1),1,ovrlppoweroldp(1,1),1)

              nullify(ovrlpminonep)

              ovrlppowerp = f_malloc_ptr((/norb,norb/), id='ovrlppowerp')

              if (power==1) then
                  factor=-1.0d0
              else if (power==2) then
                  factor=0.5d0
              else if (power==-2) then
                  factor=-0.5d0
              end if
          end if

          if (nproc>1) then
              inv_ovrlpp = f_malloc_ptr((/norb,norb/), id='inv_ovrlpp')
          else
              inv_ovrlpp => inv_ovrlp_matrix
          end if

          call first_order_taylor_dense(norb,0,norb,power,ovrlp_matrix(1,1),inv_ovrlpp)

          do i=2,iorder
              call dgemm('n', 'n', norb, norb, norb, 1.d0, ovrlpminone(1,1), &
                   norb, ovrlppoweroldp(1,1), norb, 0.d0, ovrlppowerp(1,1), norb)
              factor=newfactor(power,i,factor)
              call daxpy(norb*norb,factor,ovrlppowerp,1,inv_ovrlpp,1)
              if (i/=iorder) call vcopy(norb*norb,ovrlppowerp(1,1),1,ovrlppoweroldp(1,1),1)
          end do

          if (iorder>1) then
              if(nproc > 1) then
                  nullify(ovrlpminone)
              else
                  call f_free_ptr(ovrlpminone)
              end if

              call f_free_ptr(ovrlppowerp)
              call f_free_ptr(ovrlppoweroldp)

          end if

          nullify(inv_ovrlpp)
      end if

      if (check_accur) then
          call check_accur_overlap_minus_one(iproc,nproc,norb,norb,0,power,ovrlp_matrix,inv_ovrlp_matrix,error)
      end if


  call timing(iproc,'lovrlp^-1     ','OF')
  call f_release_routine()


end subroutine overlap_minus_one_half_serial
