!> @file
!! Orthonormalization
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Orthonormalized the localized orbitals
subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, max_inversion_error, npsidim_orbs, &
           orbs, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed, foe_obj)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices, &
                               assignment(=), sparsematrix_malloc_ptr, SPARSE_TASKGROUP
  use sparsematrix, only: compress_matrix, uncompress_matrix, gather_matrix_from_taskgroups_inplace
  use foe_base, only: foe_data
  use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed, &
                                   normalize_transposed
  use matrix_operations, only: overlapPowerGeneral
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,npsidim_orbs
  integer,intent(inout) :: methTransformOverlap
  real(kind=8),intent(in) :: max_inversion_error
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
  real(kind=8), dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  character(len=*), parameter :: subname='orthonormalizeLocalized'
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: mean_error, max_error
  logical :: ovrlp_associated, inv_ovrlp_associated
  type(matrices) :: ovrlp_
  type(matrices),dimension(1) :: inv_ovrlp_half_
  integer :: ii, i, ispin


  call f_routine(id='orthonormalizeLocalized')

  inv_ovrlp_half_(1) = matrices_null()
  call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_(1))
  !inv_ovrlp_half_(1)%matrix_compr = sparsematrix_malloc_ptr(inv_ovrlp_half, &
  !                                  iaction=SPARSE_TASKGROUP, id='inv_ovrlp_half_(1)%matrix_compr')



  if(.not.can_use_transposed) then
      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_FULL, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.
      !!do i=1,collcom%ndimind_c
      !!    write(750+iproc,'(a,2i8,es14.5)') 'i, mod(i-1,ndimind_c/2)+1, val', i, mod(i-1,collcom%ndimind_c/2)+1, psit_c(i)
      !!end do
  end if

  ovrlp_ = matrices_null()
  call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
  !!ii=0
  !!do ispin=1,ovrlp%nspin
  !!    do i=1,ovrlp%nvctr
  !!        ii=ii+1
  !!        write(930+iproc,*) 'ii, i, val', ii, i, ovrlp_%matrix_compr(ii)
  !!    end do
  !!end do


  if (methTransformOverlap==-1) then
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, ovrlp, ovrlp_)
      call overlap_power_minus_one_half_parallel(iproc, nproc, 0, orbs, ovrlp, ovrlp_, inv_ovrlp_half, inv_ovrlp_half_(1))
  else
      call overlapPowerGeneral(iproc, nproc, methTransformOverlap, 1, (/-2/), &
           orthpar%blocksize_pdgemm, &
           imode=1, ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
           ovrlp_mat=ovrlp_, inv_ovrlp_mat=inv_ovrlp_half_, &
           check_accur=.true., mean_error=mean_error, max_error=max_error)!!, &
      !if (iproc==0) call yaml_map('max error',max_error)
      !if (iproc==0) call yaml_map('mean error',mean_error)
      call check_taylor_order(mean_error, max_inversion_error, methTransformOverlap)
      !!ii=0
      !!do ispin=1,inv_ovrlp_half%nspin
      !!    do i=1,inv_ovrlp_half%nvctr
      !!        ii=ii+1
      !!        write(1930+iproc,*) 'ii, i, val', ii, i, inv_ovrlp_half_%matrix_compr(ii)
      !!    end do
      !!end do
  end if

  call deallocate_matrices(ovrlp_)

  psittemp_c = f_malloc(sum(collcom%nrecvcounts_c),id='psittemp_c')
  psittemp_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psittemp_f')

  call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
  call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)

  if (methTransformOverlap==-1) then
      ! this is only because overlap_power_minus_one_half_parallel still needs the entire array without taskgroups... to be improved
  call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_(1), &
       psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
  else
  call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_(1), &
       psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
  end if
  !if (iproc==0) then
  !    do i=1,size(inv_ovrlp_half_(1)%matrix_compr)
  !        write(*,*) 'i, val', i, inv_ovrlp_half_(1)%matrix_compr(i)
  !    end do
  !end if


  norm = f_malloc(orbs%norb,id='norm')
  call normalize_transposed(iproc, nproc, orbs, ovrlp%nspin, collcom, psit_c, psit_f, norm)

  call f_free(norm)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
       TRANSPOSE_FULL, psit_c, psit_f, lphi, lzd)

  call f_free(psittemp_c)
  call f_free(psittemp_f)

  !call f_free_ptr(inv_ovrlp_half%matrix_compr)

  call deallocate_matrices(inv_ovrlp_half_(1))

  call f_release_routine()

end subroutine orthonormalizeLocalized


!> Can still tidy this up more when tmblarge is removed
!! use sparsity of density kernel for all inverse quantities
subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, npsidim_orbs, npsidim_comp, orbs, collcom, orthpar, &
           correction_orthoconstraint, linmat, lphi, lhphi, lagmat, lagmat_, psit_c, psit_f, &
           hpsit_c, hpsit_f, &
           can_use_transposed, overlap_calculated, experimental_mode, calculate_inverse, norder_taylor, max_inversion_error, &
           npsidim_orbs_small, lzd_small, hpsi_noprecond, wt_philarge, wt_hphi, wt_hpsinoprecond)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  use yaml_output
  use communications_base, only: work_transpose, TRANSPOSE_POST, TRANSPOSE_FULL, TRANSPOSE_GATHER
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: matrices_null, allocate_matrices, deallocate_matrices, sparsematrix_malloc, &
                               sparsematrix_malloc_ptr, DENSE_FULL, DENSE_MATMUL, SPARSE_FULL, SPARSEMM_SEQ, &
                               assignment(=), SPARSE_TASKGROUP
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: uncompress_matrix, &
                          sequential_acces_matrix_fast2, transform_sparse_matrix, &
                          gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace, &
                          transform_sparse_matrix_local, uncompress_matrix_distributed2, &
                          matrix_matrix_mult_wrapper
  use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed
  use matrix_operations, only: overlapPowerGeneral
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim_orbs, npsidim_comp, npsidim_orbs_small
  type(local_zone_descriptors),intent(in) :: lzd, lzd_small
  !type(orbitals_Data),intent(in) :: orbs
  type(orbitals_Data),intent(inout) :: orbs !temporary inout
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  integer,intent(in) :: correction_orthoconstraint
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(in) :: lphi
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(inout) :: lhphi
  type(sparse_matrix),intent(inout) :: lagmat
  type(matrices),intent(out) :: lagmat_
  real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: hpsit_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: hpsit_f
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed, overlap_calculated
  type(linear_matrices),intent(inout) :: linmat ! change to ovrlp and inv_ovrlp, and use inv_ovrlp instead of denskern
  logical,intent(in) :: experimental_mode, calculate_inverse
  integer,intent(inout) :: norder_taylor
  real(kind=8),intent(in) :: max_inversion_error
  real(kind=8),dimension(npsidim_orbs_small),intent(out) :: hpsi_noprecond
  type(work_transpose),intent(inout) :: wt_philarge
  type(work_transpose),intent(out) :: wt_hphi, wt_hpsinoprecond

  ! Local variables
  real(kind=8) :: max_error, mean_error
  real(kind=8),dimension(:),allocatable :: tmp_mat_compr, hpsit_tmp_c, hpsit_tmp_f, hphi_nococontra
  integer,dimension(:),allocatable :: ipiv
  type(matrices),dimension(1) :: inv_ovrlp_
  real(8),dimension(:),allocatable :: inv_ovrlp_seq, lagmat_large, tmpmat, tmparr
  real(8),dimension(:,:),allocatable :: lagmatp, inv_lagmatp
  integer,dimension(2) :: irowcol
  real(kind=8),dimension(:),pointer :: matrix_local
  integer,parameter :: GLOBAL_MATRIX=101, SUBMATRIX=102
  integer,parameter :: data_strategy_main=SUBMATRIX!GLOBAL_MATRIX
  !type(work_transpose) :: wt_

  call f_routine(id='orthoconstraintNonorthogonal')



  inv_ovrlp_seq = sparsematrix_malloc(linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_seq')
  inv_lagmatp = sparsematrix_malloc(linmat%l, iaction=DENSE_MATMUL, id='inv_lagmatp')
  lagmatp = sparsematrix_malloc(linmat%l, iaction=DENSE_MATMUL, id='lagmatp')
  !!inv_ovrlp_(1) = matrices_null()
  !!inv_ovrlp_(1)%matrix_compr = sparsematrix_malloc_ptr(linmat%l,iaction=SPARSE_FULL,id='inv_ovrlp_(1)%matrix_compr')

  if (calculate_inverse) then
      ! Invert the overlap matrix
      if (iproc==0) call yaml_map('calculation of S^-1','direct calculation')
      !!tmparr = sparsematrix_malloc(linmat%s,iaction=SPARSE_FULL,id='tmparr')
      !!call vcopy(linmat%s%nvctr*linmat%s%nspin, linmat%ovrlp_%matrix_compr(1), 1, tmparr(1), 1)
      !!call extract_taskgroup_inplace(linmat%s, linmat%ovrlp_)
      call overlapPowerGeneral(iproc, nproc, norder_taylor, 1, (/1/), -1, &
           imode=1, ovrlp_smat=linmat%s, inv_ovrlp_smat=linmat%l, &
           ovrlp_mat=linmat%ovrlp_, inv_ovrlp_mat=linmat%ovrlppowers_(3), &
           check_accur=.true., max_error=max_error, mean_error=mean_error)
      !!call vcopy(linmat%s%nvctr*linmat%s%nspin, tmparr(1), 1, linmat%ovrlp_%matrix_compr(1), 1)
      !!call f_free(tmparr)
      call check_taylor_order(mean_error, max_inversion_error, norder_taylor)
  else
      if (iproc==0) call yaml_map('calculation of S^-1','from memory')
  end if

  ! Gather together the data (has been posted in getLocalizedBasis
  call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
       TRANSPOSE_GATHER, lphi, psit_c, psit_f, lzd, wt_philarge)
  can_use_transposed=.true.

  ! Calculate <phi_alpha|g_beta>
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat, lagmat_)
  !call gather_matrix_from_taskgroups_inplace(iproc, nproc, lagmat, lagmat_)

  lagmat_large = sparsematrix_malloc(linmat%l, iaction=SPARSE_TASKGROUP, id='lagmat_large')

  call symmetrize_matrix()


  ! Apply S^-1
  !!call sequential_acces_matrix_fast2(linmat%l, linmat%ovrlppowers_(3)%matrix_compr, inv_ovrlp_seq)
  ! Transform the matrix to the large sparsity pattern (necessary for the following uncompress_matrix_distributed)
  if (correction_orthoconstraint==0) then
      if (data_strategy_main==GLOBAL_MATRIX) then
          stop 'deprecated'
          call transform_sparse_matrix(linmat%m, linmat%l, lagmat_%matrix_compr, lagmat_large, 'small_to_large')
      end if
      if (iproc==0) call yaml_map('correction orthoconstraint',.true.)
      !!call uncompress_matrix_distributed2(iproc, linmat%l, DENSE_MATMUL, lagmat_large, lagmatp)
      !!call sparsemm(linmat%l, inv_ovrlp_seq, lagmatp, inv_lagmatp)
      !!write(*,*) 'iproc, sum(inv_lagmatp)', iproc, sum(inv_lagmatp)
      !!call compress_matrix_distributed(iproc, nproc, linmat%l, DENSE_MATMUL, &
      !!     inv_lagmatp, lagmat_large)
      call matrix_matrix_mult_wrapper(iproc, nproc, linmat%l, &
           linmat%ovrlppowers_(3)%matrix_compr, lagmat_large, lagmat_large)
  end if
  if (data_strategy_main==SUBMATRIX) then
      call transform_sparse_matrix_local(linmat%m, linmat%l, lagmat_%matrix_compr, lagmat_large, 'large_to_small')
  end if
  call f_free(lagmat_large)
  call f_free(inv_ovrlp_seq)
  call f_free(lagmatp)
  call f_free(inv_lagmatp)

  !!tmparr = sparsematrix_malloc(lagmat,iaction=SPARSE_FULL,id='tmparr')
  !!call vcopy(lagmat%nvctr, lagmat_%matrix_compr(1), 1, tmparr(1), 1)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, lagmat, lagmat_)
  call build_linear_combination_transposed(collcom, lagmat, lagmat_, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)
  !!call vcopy(lagmat%nvctr, tmparr(1), 1, lagmat_%matrix_compr(1), 1)
  !!call f_free(tmparr)


  ! Start the untranspose process (will be gathered together in
  ! calculate_energy_and_gradient_linear)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
       TRANSPOSE_POST, hpsit_c, hpsit_f, lhphi, lzd, wt_hphi)

  ! The symmetrized Lagrange multiplier matrix has now the wrong sign
  call dscal(lagmat%nvctrp_tg*lagmat%nspin, -1.d0, lagmat_%matrix_compr(1), 1)



  overlap_calculated=.false.
  

  ! @NEW apply S^-1 to the gradient
  hpsit_tmp_c = f_malloc(collcom%ndimind_c,id='psit_tmp_c')
  hpsit_tmp_f = f_malloc(7*collcom%ndimind_f,id='psit_tmp_f')
  hphi_nococontra = f_malloc(npsidim_orbs,id='hphi_nococontra')
  call build_linear_combination_transposed(collcom, linmat%l, linmat%ovrlppowers_(3), &
       hpsit_c, hpsit_f, .true., hpsit_tmp_c, hpsit_tmp_f, iproc)

  ! Start the untranspose process (will be gathered together in
  ! getLocalizedBasis)
  ! Pass hphi_nococontra even if it is not used
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
       TRANSPOSE_POST, hpsit_tmp_c, hpsit_tmp_f, hphi_nococontra, lzd, wt_hpsinoprecond)
  !!call large_to_small_locreg(iproc, npsidim_orbs_small, npsidim_orbs, lzd_small, lzd, &
  !!     orbs, hphi_nococontra, hpsi_noprecond)
  ! END @NEW

  !call deallocate_matrices(inv_ovrlp_(1))
  call f_free(hpsit_tmp_c)
  call f_free(hpsit_tmp_f)
  call f_free(hphi_nococontra)

  call f_release_routine()


  contains


    subroutine symmetrize_matrix()
      implicit none
      integer :: ishift, itg, iitg
      integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
      integer,parameter :: comm_strategy=GET
      integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX
      integer :: iorb, jorb, ii, ii_trans, irow, jcol, info, lwork, jj, ispin, iseg, i
      integer :: isegstart, isegend, ierr


      call f_routine(id='symmetrize_matrix')

      ! Just to check the consistency
      if (data_strategy_main/=data_strategy) then
          stop 'data_strategy_main/=data_strategy'
      end if

      if (lagmat%nvctrp>0) then
          isegstart = lagmat%istsegline(lagmat%isfvctr+1)
          isegend = lagmat%istsegline(lagmat%isfvctr+lagmat%nfvctrp) + &
                    lagmat%nsegline(lagmat%isfvctr+lagmat%nfvctrp)-1
      else
          isegstart = 1
          isegend = 0
      end if
      if (data_strategy==GLOBAL_MATRIX) then
          stop 'symmetrize_matrix: option GLOBAL_MATRIX is deprecated'
          !!matrix_local = f_malloc_ptr(max(lagmat%nvctrp,1),id='matrix_local')
          !!tmp_mat_compr = sparsematrix_malloc(lagmat,iaction=SPARSE_FULL,id='tmp_mat_compr')
          !!call vcopy(lagmat%nvctr*lagmat%nspin, lagmat_%matrix_compr(1), 1, tmp_mat_compr(1), 1)
          !!do ispin=1,lagmat%nspin
          !!    ishift=(ispin-1)*lagmat%nvctr
          !!    if (isegend>=isegstart) then
          !!        !$omp parallel default(none) &
          !!        !$omp shared(isegstart,isegend,ishift,lagmat,matrix_local,tmp_mat_compr) &
          !!        !$omp private(iseg,ii,i,irowcol,ii_trans)
          !!        !$omp do
          !!        do iseg=isegstart,isegend
          !!            ii=lagmat%keyv(iseg)
          !!            ! A segment is always on one line, therefore no double loop
          !!            do i=lagmat%keyg(1,1,iseg),lagmat%keyg(2,1,iseg)
          !!               irowcol = orb_from_index(lagmat, i)
          !!               ii_trans = matrixindex_in_compressed(lagmat,lagmat%keyg(1,2,iseg),i)
          !!               matrix_local(ii-lagmat%isvctr) = -0.5d0*tmp_mat_compr(ii+ishift)-0.5d0*tmp_mat_compr(ii_trans+ishift)
          !!               ii=ii+1
          !!            end do
          !!        end do
          !!        !$omp end do
          !!        !$omp end parallel
          !!    end if
          !!    if (nproc>1) then
          !!        !!call mpi_allgatherv(matrix_local(1), lagmat%nvctrp, mpi_double_precision, &
          !!        !!     lagmat_%matrix_compr(ishift+1), lagmat%nvctr_par, lagmat%isvctr_par, mpi_double_precision, &
          !!        !!     bigdft_mpi%mpi_comm, ierr)
          !!        if (comm_strategy==ALLGATHERV) then
          !!            call mpi_allgatherv(matrix_local(1), lagmat%nvctrp, mpi_double_precision, &
          !!                 lagmat_%matrix_compr(ishift+1), lagmat%nvctr_par, lagmat%isvctr_par, mpi_double_precision, &
          !!                 bigdft_mpi%mpi_comm, ierr)
          !!            call f_free_ptr(matrix_local)
          !!        else if (comm_strategy==GET) then
          !!            !!call mpiget(iproc, nproc, bigdft_mpi%mpi_comm, lagmat%nvctrp, matrix_local, &
          !!            !!     lagmat%nvctr_par, lagmat%isvctr_par, lagmat%nvctr, lagmat_%matrix_compr(ishift+1:ishift+lagmat%nvctr))
          !!            call mpi_get_to_allgatherv(matrix_local(1), lagmat%nvctrp, &
          !!                 lagmat_%matrix_compr(ishift+1), &
          !!                 lagmat%nvctr_par, lagmat%isvctr_par, bigdft_mpi%mpi_comm)
          !!        else
          !!            stop 'symmetrize_matrix: wrong communication strategy'
          !!        end if
          !!    else
          !!        call vcopy(lagmat%nvctr, matrix_local(1), 1, lagmat_%matrix_compr(ishift+1), 1)
          !!    end if
          !!    if (ispin==lagmat%nspin) call f_free_ptr(matrix_local)
          !!end do
      else if (data_strategy==SUBMATRIX) then
          ! Directly use the large sparsity pattern as this one is used later
          ! for the matrix vector multiplication
          tmp_mat_compr = sparsematrix_malloc(linmat%l,iaction=SPARSE_TASKGROUP,id='tmp_mat_compr')
          call transform_sparse_matrix_local(linmat%m, linmat%l, lagmat_%matrix_compr, tmp_mat_compr, 'small_to_large')
          do ispin=1,lagmat%nspin
              ishift=(ispin-1)*linmat%l%nvctrp_tg
              !$omp parallel default(none) &
              !$omp shared(linmat,lagmat_large,tmp_mat_compr,ishift) &
              !$omp private(iseg,ii,i,ii_trans)
              !$omp do
              !do iseg=linmat%l%iseseg_tg(1),linmat%l%iseseg_tg(2)
              do iseg=linmat%l%istartendseg_local(1),linmat%l%istartendseg_local(2)
                  ii = linmat%l%keyv(iseg)
                  ! A segment is always on one line, therefore no double loop
                  do i=linmat%l%keyg(1,1,iseg),linmat%l%keyg(2,1,iseg) !this is too much, but for the moment ok
                      ii_trans = matrixindex_in_compressed(linmat%l,linmat%l%keyg(1,2,iseg),i)
                      lagmat_large(ii+ishift-linmat%l%isvctrp_tg) = &
                          - 0.5d0*tmp_mat_compr(ii+ishift-linmat%l%isvctrp_tg) &
                          - 0.5d0*tmp_mat_compr(ii_trans+ishift-linmat%l%isvctrp_tg)
                      ii=ii+1
                  end do
              end do
              !$omp end do
              !$omp end parallel
          end do
      else
          stop 'symmetrize_matrix: wrong data strategy'
      end if

      call f_free(tmp_mat_compr)
      call f_release_routine()

    end subroutine symmetrize_matrix


end subroutine orthoconstraintNonorthogonal


!!subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in) :: mpisource, mpidest, istsource, istdest, ncount, tag
!!  integer,dimension(6),intent(out) :: comarr
!!
!!
!!  ! From which MPI process shall the orbital be sent.
!!  comarr(1)=mpisource
!!
!!  ! Starting index on the sending process.
!!  comarr(2)=istsource
!!
!!  ! Amount of datat to be sent
!!  comarr(3)=ncount
!!
!!  ! To which MPI process shall the orbital be sent.
!!  comarr(4)=mpidest
!!
!!  ! Starting index on the receiving process.
!!  comarr(5)=istdest
!!
!!  ! Tag for the communication
!!  comarr(6)=tag
!!
!!  ! comarr(7): this entry is used as request for the mpi_isend.
!!
!!  ! comarr(8): this entry is used as request for the mpi_irecv.
!!
!!end subroutine setCommsParameters





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

subroutine overlap_power_minus_one_half_parallel(iproc, nproc, meth_overlap, orbs, ovrlp, ovrlp_mat, &
           inv_ovrlp_half, inv_ovrlp_half_)
  use module_base
  use module_types
  use module_interfaces
  use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, &
                               allocate_matrices, deallocate_matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  use matrix_operations, only: overlap_plus_minus_one_half_exact
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, meth_overlap
  type(orbitals_data),intent(in) :: orbs
  type(sparse_matrix),intent(inout) :: ovrlp
  type(matrices),intent(inout) :: ovrlp_mat
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half
  type(matrices),intent(inout) :: inv_ovrlp_half_

  ! Local variables
  integer(kind=8) :: ii, iend
  integer :: i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
  integer :: iiorb, ierr, iseg, ind, ishift_ovrlp, ishift_inv_ovrlp, ispin
  real(kind=8) :: error
  real(kind=8),dimension(:,:),pointer :: ovrlp_tmp, ovrlp_tmp_inv_half
  logical,dimension(:),allocatable :: in_neighborhood
  character(len=*),parameter :: subname='overlap_power_minus_one_half_parallel'
  !type(matrices) :: inv_ovrlp_half_
  !!integer :: itaskgroups, iitaskgroup, imin, imax

  !!imin=ovrlp%nvctr
  !!imax=0
  !!do itaskgroups=1,ovrlp%ntaskgroupp
  !!    iitaskgroup = ovrlp%taskgroupid(itaskgroups)
  !!    imin = min(imin,ovrlp%taskgroup_startend(1,1,iitaskgroup))
  !!    imax = max(imax,ovrlp%taskgroup_startend(2,1,iitaskgroup))
  !!end do


  call timing(iproc,'lovrlp^-1/2par','ON')
  call f_routine('overlap_power_minus_one_half_parallel')

  in_neighborhood = f_malloc(ovrlp%nfvctr,id='in_neighborhood')

  !inv_ovrlp_half_ = matrices_null()
  !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_)
  call f_zero(inv_ovrlp_half%nvctrp_tg*inv_ovrlp_half%nspin, inv_ovrlp_half_%matrix_compr(1))

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

  !!ii=0
  !!do ispin=1,ovrlp%nspin
  !!    do i=1,ovrlp%nvctr
  !!        ii=ii+1
  !!        write(950+iproc,*) 'ii, i, val', ii, i, ovrlp_mat%matrix_compr(ii)
  !!    end do
  !!end do


  spin_loop: do ispin=1,ovrlp%nspin

      ishift_ovrlp=(ispin-1)*ovrlp%nvctr
      ishift_inv_ovrlp=(ispin-1)*inv_ovrlp_half%nvctr

      do iorb=1,ovrlp%nfvctrp
         iiorb=ovrlp%isfvctr+iorb
         !ilr=orbs%inwhichlocreg(iiorb)
         ! We are at the start of a new atom
         ! Count all orbitals that are in the neighborhood

         iseg=ovrlp%istsegline(iiorb)
         iend=int(iiorb,kind=8)*int(ovrlp%nfvctr,kind=8)
         n=0
         in_neighborhood(:)=.false.
         do 
            do i=ovrlp%keyg(1,1,iseg),ovrlp%keyg(2,1,iseg)
               in_neighborhood(i)=.true.
               !if (iproc==0) write(*,*) 'iiorb, iseg, i, n', iiorb, iseg, i, n
               n=n+1
            end do
            iseg=iseg+1
            if (iseg>ovrlp%nseg) exit
            ii = int((ovrlp%keyg(1,2,iseg)-1),kind=8)*int(ovrlp%nfvctr,kind=8) + &
                 int(ovrlp%keyg(1,1,iseg),kind=8)
            if (ii>iend) exit
         end do
         !if (iproc==0) write(*,*) 'iiorb, n', iiorb, n

         ovrlp_tmp = f_malloc0_ptr((/n,n/),id='ovrlp_tmp')

         jjorb=0
         do jorb=1,ovrlp%nfvctr
            if (.not.in_neighborhood(jorb)) cycle
            jjorb=jjorb+1
            kkorb=0
            do korb=1,ovrlp%nfvctr
               if (.not.in_neighborhood(korb)) cycle
               kkorb=kkorb+1
               ind = matrixindex_in_compressed(ovrlp,korb,jorb)
               if (ind>0) then
                  !!if (ind<imin) then
                  !!    write(*,*) 'ind,imin',ind,imin
                  !!    stop 'ind<imin'
                  !!end if
                  !!if (ind>imax) then
                  !!    write(*,*) 'ind,imax',ind,imax
                  !!    stop 'ind>imax'
                  !!end if
                  ind=ind+ishift_ovrlp
                  ovrlp_tmp(kkorb,jjorb)=ovrlp_mat%matrix_compr(ind-ovrlp%isvctrp_tg)
               else
                  ovrlp_tmp(kkorb,jjorb)=0.d0
               end if
               !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
            end do
         end do
              
         ovrlp_tmp_inv_half = f_malloc_ptr((/n,n/),id='ovrlp_tmp_inv_half')
         call vcopy(n*n, ovrlp_tmp(1,1), 1, ovrlp_tmp_inv_half(1,1), 1)
         !!do jorb=1,n
         !!    do korb=1,n
         !!        write(900,'(a,2i8,es14.5)') 'jorb, korb, ovrlp_tmp(korb,jorb)', jorb, korb, ovrlp_tmp(korb,jorb)
         !!    end do
         !!end do

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
         call overlap_plus_minus_one_half_exact(1, n, -8, .false., ovrlp_tmp_inv_half,inv_ovrlp_half)


         !if (iiorb==orbs%norb) then
         !print*,''
         !print*,'inv_ovrlp_tmp',n,iiorb,error
         !do jorb=1,n
         !print*,jorb,ovrlp_tmp_inv_half(:,jorb)
         !end do
         !end if

         jjorb=0
         do jorb=1,ovrlp%nfvctr
            if (.not.in_neighborhood(jorb)) cycle
            jjorb=jjorb+1
            kkorb=0
            if (jorb==iiorb) then
               do korb=1,ovrlp%nfvctr
                  if (.not.in_neighborhood(korb)) cycle
                  kkorb=kkorb+1
                  ind = matrixindex_in_compressed(inv_ovrlp_half,korb,jorb)
                  if (ind>0) then
                     ind=ind+ishift_inv_ovrlp
                     inv_ovrlp_half_%matrix_compr(ind-inv_ovrlp_half%isvctrp_tg)=ovrlp_tmp_inv_half(kkorb,jjorb)
                     !if (iproc==0) write(*,'(a,6i8,es16.8)') 'ind, inv_ovrlp_half%isvctrp_tg, jorb, korb, jjorb, kkorb, val', &
                     !                  ind, inv_ovrlp_half%isvctrp_tg, jorb, korb, jjorb, kkorb, ovrlp_tmp_inv_half(kkorb,jjorb)
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

      !!if (nproc>1)then
      !!    call mpiallred(inv_ovrlp_half_%matrix_compr(1), inv_ovrlp_half%nvctr*inv_ovrlp_half%nspin, mpi_sum, bigdft_mpi%mpi_comm)
      !!end if
      call synchronize_matrix_taskgroups(iproc, nproc, inv_ovrlp_half, inv_ovrlp_half_)

  end do spin_loop

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


  call f_release_routine
  call timing(iproc,'lovrlp^-1/2par','OF')

end subroutine overlap_power_minus_one_half_parallel

!> Orthonormalize a subset of orbitals
subroutine orthonormalize_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalize_subset
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices, &
                               sparsematrix_malloc, assignment(=), SPARSE_FULL
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
  use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed, &
                                   normalize_transposed
  use matrix_operations, only: overlapPowerGeneral
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
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm, tmparr
  !type(sparse_matrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='orthonormalize_subset'
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: max_error, mean_error
  type(matrices) :: ovrlp_
  type(matrices),dimension(1) :: inv_ovrlp_half_

  call f_routine(id='orthonormalize_subset')


  !call nullify_sparse_matrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)

  inv_ovrlp_half_(1) = matrices_null()
  call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_(1))


  if(.not.can_use_transposed) then
      !!if(associated(psit_c)) then
      !!    call f_free_ptr(psit_c)
      !!end if
      !!if(associated(psit_f)) then
      !!    call f_free_ptr(psit_f)
      !!end if
      !!psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
      !!psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')

      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_FULL, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.

  end if

  ovrlp_ = matrices_null()
  call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, ovrlp, ovrlp_)
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
           orbs, ovrlp, ovrlp_, inv_ovrlp_half, inv_ovrlp_half_(1))
      !call deallocate_matrices(ovrlp_)
  else
      nullify(inv_ovrlp_null)
      ! do sparse.. check later
      !ovrlp%matrix_compr=ovrlp_%matrix_compr
      tmparr = sparsematrix_malloc(ovrlp,iaction=SPARSE_FULL,id='tmparr')
      call vcopy(ovrlp%nvctr*ovrlp%nspin, ovrlp_%matrix_compr(1), 1, tmparr(1), 1)
      call extract_taskgroup_inplace(ovrlp, ovrlp_)
      call overlapPowerGeneral(iproc, nproc, methTransformOverlap, 1, (/-2/), &
           orthpar%blocksize_pdsyev, &
           imode=1, check_accur=.true., &
           ovrlp_mat=ovrlp_, inv_ovrlp_mat=inv_ovrlp_half_, &
           ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
           max_error=max_error, mean_error=mean_error)
      call vcopy(ovrlp%nvctr*ovrlp%nspin, tmparr(1), 1, ovrlp_%matrix_compr(1), 1)
      call f_free(tmparr)
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
  call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_(1), &
       psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)



  call deallocate_matrices(ovrlp_)

  call deallocate_matrices(inv_ovrlp_half_(1))


  norm = f_malloc(orbs%norb,id='norm')
  call normalize_transposed(iproc, nproc, orbs, ovrlp%nspin, collcom, psit_c, psit_f, norm)
  call f_free(norm)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, & 
       TRANSPOSE_FULL, psit_c, psit_f, lphi, lzd)

  call f_free(psittemp_c)
  call f_free(psittemp_f)
  !!call f_free_ptr(inv_ovrlp_half%matrix_compr)

  call f_release_routine()

end subroutine orthonormalize_subset



subroutine gramschmidt_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => gramschmidt_subset
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace
  use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed
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

  call f_routine('gramschmidt_subset')


  !call nullify_sparse_matrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
  !!inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')


  if(.not.can_use_transposed) then
      !!if(associated(psit_c)) then
      !!    call f_free_ptr(psit_c)
      !!end if
      !!if(associated(psit_f)) then
      !!    call f_free_ptr(psit_f)
      !!end if
      !!psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
      !!psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')

      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_FULL, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.

  end if


  ovrlp_ = matrices_null()
  call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, ovrlp, ovrlp_)
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
      call f_zero(jcount_norb)
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



  !inv_ovrlp_ = matrices_null()
  !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
  !@WARNING CHECK THIS
  !!ovrlp_%matrix_compr = inv_ovrlp_half%matrix_compr
  call build_linear_combination_transposed(collcom, ovrlp, ovrlp_, &
       psittemp_c, psittemp_f, .false., psit_c, psit_f, iproc)
  !call deallocate_matrices(ovrlp_)


  call deallocate_matrices(ovrlp_)


  norm = f_malloc(orbs%norb,id='norm')
  !call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)

  call f_free(norm)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
       TRANSPOSE_FULL, psit_c, psit_f, lphi, lzd)

  call f_free(psittemp_c)
  call f_free(psittemp_f)
  !!call f_free_ptr(inv_ovrlp_half%matrix_compr)

  call f_release_routine()

end subroutine gramschmidt_subset

subroutine gramschmidt_coeff(iproc,nproc,norb,basis_orbs,basis_overlap,basis_overlap_mat,coeff)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix, matrices
  implicit none

  integer, intent(in) :: iproc, nproc, norb
  type(orbitals_data), intent(in) :: basis_orbs
  type(sparse_matrix),intent(inout) :: basis_overlap
  type(matrices),intent(inout) :: basis_overlap_mat
  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(inout) :: coeff

  integer :: iorb, jtmb, corb, ierr
  real(kind=8), dimension(:,:), allocatable :: ovrlp_coeff, coeff_tmp, coeff_trans
  real(kind=4) :: tr0, tr1
  real(kind=8) :: time0, time1, time2, time3, time4, time5

  call f_routine(id='gramschmidt_coeff')

  time0=0.0d0
  time1=0.0d0
  time2=0.0d0
  time3=0.0d0
  time4=0.0d0
  time5=0.0d0

  ! orthonormalizing all iorb<corb wrt corb (assume original vectors were normalized)
  do corb=norb,1,-1
     ovrlp_coeff=f_malloc((/corb,1/),id='ovrlp_coeff')
     coeff_tmp=f_malloc((/corb,basis_orbs%norbp/),id='coeff_tmp')
     ! calculate relevant part of cSc

     call cpu_time(tr0)
     if (basis_orbs%norbp>0) then
        call dgemm('t', 'n', corb, basis_orbs%norbp, basis_orbs%norb, 1.d0, &
             coeff(1,1), basis_orbs%norb, &
             basis_overlap_mat%matrix(1,basis_orbs%isorb+1,1), basis_orbs%norb, 0.d0, &
             coeff_tmp(1,1), corb)

        call dgemm('n', 'n', corb, 1, basis_orbs%norbp, 1.d0, &
             coeff_tmp(1,1), corb, &
             coeff(basis_orbs%isorb+1,corb), basis_orbs%norb, 0.d0, &
             ovrlp_coeff(1,1), corb)
     else
        call f_zero(ovrlp_coeff)
     end if

     if (nproc>1) then
        call mpiallred(ovrlp_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
     end if

     call cpu_time(tr1)
     time1=time1+real(tr1-tr0,kind=8)

     ! (c_corb S c_iorb) * c_corb
     if (basis_orbs%norbp>0) then
        call dgemm('n', 't', corb-1, basis_orbs%norbp, 1, 1.d0, &
             ovrlp_coeff(1,1), corb, &
             coeff(1+basis_orbs%isorb,corb), basis_orbs%norb, 0.d0, &
             coeff_tmp(1,1), corb)
     end if
     call cpu_time(tr0)
     time2=time2+real(tr0-tr1,kind=8)

     ! sum and transpose coeff for allgatherv
     !$omp parallel do default(private) shared(coeff,coeff_tmp,corb,basis_orbs,ovrlp_coeff)
     do iorb=1,corb-1
        do jtmb=1,basis_orbs%norbp
           coeff_tmp(iorb,jtmb) = coeff(jtmb+basis_orbs%isorb,iorb) - coeff_tmp(iorb,jtmb)/ovrlp_coeff(corb,1)
        end do
     end do
     !$omp end parallel do

     do jtmb=1,basis_orbs%norbp
        coeff_tmp(corb,jtmb) = coeff(jtmb+basis_orbs%isorb,corb)/sqrt(ovrlp_coeff(corb,1))
     end do
     call cpu_time(tr1)
     time3=time3+real(tr1-tr0,kind=8)

     call f_free(ovrlp_coeff)
     coeff_trans=f_malloc((/corb,basis_orbs%norb/),id='coeff_tmp')
     if(nproc > 1) then
        call mpi_allgatherv(coeff_tmp(1,1), basis_orbs%norbp*corb, mpi_double_precision, coeff_trans(1,1), &
           corb*basis_orbs%norb_par(:,0), corb*basis_orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call vcopy(basis_orbs%norbp*corb,coeff_tmp(1,1),1,coeff_trans(1,1),1)
     end if
     call f_free(coeff_tmp)

     call cpu_time(tr0)
     time4=time4+real(tr0-tr1,kind=8)

     ! untranspose coeff
     !$omp parallel do default(private) shared(coeff,coeff_trans,corb,basis_orbs)
     do jtmb=1,basis_orbs%norb
        do iorb=1,corb
           coeff(jtmb,iorb) = coeff_trans(iorb,jtmb)
        end do
     end do
     !$omp end parallel do

     call cpu_time(tr1)
     time5=time5+real(tr1-tr0,kind=8)

     call f_free(coeff_trans)
  end do

  !if (iproc==0) print*,'Time in gramschmidt_coeff',time0,time1,time2,time3,time4,time5,&
  !     time0+time1+time2+time3+time4+time5

  call f_release_routine()

end subroutine gramschmidt_coeff

subroutine gramschmidt_coeff_trans(iproc,nproc,norbu,norb,basis_orbs,basis_overlap,basis_overlap_mat,coeff)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix, matrices
  implicit none

  integer, intent(in) :: iproc, nproc, norbu, norb
  type(orbitals_data), intent(in) :: basis_orbs
  type(sparse_matrix),intent(inout) :: basis_overlap
  type(matrices),intent(inout) :: basis_overlap_mat
  real(kind=8),dimension(basis_overlap%nfvctr,basis_orbs%norb),intent(inout) :: coeff

  integer :: iorb, jtmb, corb, ierr, isorb, ieorb, ispin
  real(kind=8), dimension(:,:), allocatable :: ovrlp_coeff, coeff_tmp, coeff_trans, coeff_transp

  real(kind=4) :: tr0, tr1
  real(kind=8) :: time0, time1, time2, time3, time4, time5

  call f_routine(id='gramschmidt_coeff_trans')

  time0=0.0d0
  time1=0.0d0
  time2=0.0d0
  time3=0.0d0
  time4=0.0d0
  time5=0.0d0

  call cpu_time(tr0)

  spin_loop: do ispin=1,basis_overlap%nspin

      if (ispin==1) then
          isorb=1
          ieorb=norbu
      else if (ispin==2) then
          isorb=norbu+1
          ieorb=norb
      end if

      coeff_transp=f_malloc((/(ieorb-isorb+1),basis_overlap%nfvctrp/),id='coeff_transp')
      !$omp parallel do default(private) shared(coeff,coeff_transp,isorb,ieorb,basis_overlap)
      do iorb=isorb,ieorb
         do jtmb=1,basis_overlap%nfvctrp
            coeff_transp(iorb-isorb+1,jtmb) = coeff(jtmb+basis_overlap%isfvctr,iorb)
         end do
      end do
      !$omp end parallel do

      call cpu_time(tr1)
      time0=time0+real(tr1-tr0,kind=8)

      ! orthonormalizing all iorb<corb wrt corb (assume original vectors were normalized)
      do corb=ieorb-isorb+1,1,-1


         call cpu_time(tr0)

         coeff_tmp=f_malloc((/basis_overlap%nfvctr,1/),id='coeff_tmp')
         ! calculate relevant part of cSc
         if (basis_overlap%nfvctrp>0) then
            call dgemm('n', 't', basis_overlap%nfvctr, 1, basis_overlap%nfvctrp, 1.d0, &
                 basis_overlap_mat%matrix(1,basis_overlap%isfvctr+1,1), basis_overlap%nfvctr, &
                 coeff_transp(corb,1), ieorb-isorb+1, 0.d0, &
                 coeff_tmp(1,1), basis_overlap%nfvctr)
         else
            !!call to_zero(corb,coeff_tmp(1,1)) !!!LG: is this a typo?
            call f_zero(coeff_tmp)
         end if

         if (nproc>1) then
            call mpiallred(coeff_tmp, mpi_sum, comm=bigdft_mpi%mpi_comm)
         end if

         call cpu_time(tr1)
         time1=time1+real(tr1-tr0,kind=8)

         ovrlp_coeff=f_malloc((/corb,1/),id='ovrlp_coeff')
         if (basis_overlap%nfvctrp>0) then
            call dgemm('n', 'n', corb, 1, basis_overlap%nfvctrp, 1.d0, &
                 coeff_transp(1,1), ieorb-isorb+1, &
                 coeff_tmp(1+basis_overlap%isfvctr,1), basis_overlap%nfvctr, 0.d0, &
                 ovrlp_coeff(1,1), corb)
         else
            call f_zero(ovrlp_coeff)
         end if

         if (nproc>1) then
            call mpiallred(ovrlp_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
         end if
         call f_free(coeff_tmp)

         call cpu_time(tr0)
         time2=time2+real(tr0-tr1,kind=8)

         ! (c_corb S c_iorb) * c_corb
         coeff_tmp=f_malloc((/corb,basis_overlap%nfvctrp/),id='coeff_tmp')
         if (basis_overlap%nfvctrp>0) then
            call dgemm('n', 'n', corb-1, basis_overlap%nfvctrp, 1, 1.d0, &
                 ovrlp_coeff(1,1), corb, &
                 coeff_transp(corb,1), ieorb-isorb+1, 0.d0, &
                 coeff_tmp(1,1), corb)
         end if

         call cpu_time(tr1)
         time3=time3+real(tr1-tr0,kind=8)

         ! sum and transpose coeff for allgatherv
         !$omp parallel do default(private) shared(coeff_transp,coeff_tmp,isorb,corb,basis_overlap,ovrlp_coeff)
         do iorb=isorb,corb-1
            do jtmb=1,basis_overlap%nfvctrp
               coeff_transp(iorb,jtmb) = coeff_transp(iorb,jtmb) - coeff_tmp(iorb,jtmb)/ovrlp_coeff(corb,1)
            end do
         end do
         !$omp end parallel do

         do jtmb=1,basis_overlap%nfvctrp
            coeff_transp(corb,jtmb) = coeff_transp(corb,jtmb)/sqrt(ovrlp_coeff(corb,1))
         end do

         call cpu_time(tr0)
         time4=time4+real(tr0-tr1,kind=8)

         call f_free(ovrlp_coeff)
         call f_free(coeff_tmp)
      end do

      call cpu_time(tr0)

      coeff_trans=f_malloc((/(ieorb-isorb+1),basis_overlap%nfvctr/),id='coeff_tmp')
      if(nproc > 1) then
         call mpi_allgatherv(coeff_transp(1,1), basis_overlap%nfvctrp*(ieorb-isorb+1), mpi_double_precision, coeff_trans(1,1), &
              (ieorb-isorb+1)*basis_overlap%nfvctr_par(:), (ieorb-isorb+1)*basis_overlap%isfvctr_par, &
              mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call vcopy(basis_overlap%nfvctrp*(ieorb-isorb+1),coeff_transp(1,1),1,coeff_trans(1,1),1)
      end if
      call f_free(coeff_transp)

      call cpu_time(tr1)
      time5=time5+real(tr1-tr0,kind=8)

      ! untranspose coeff
      !$omp parallel do default(private) shared(coeff,coeff_trans,isorb,ieorb,basis_overlap)
      do jtmb=1,basis_overlap%nfvctr
         do iorb=isorb,ieorb
            coeff(jtmb,iorb) = coeff_trans(iorb-isorb+1,jtmb)
         end do
      end do
      !$omp end parallel do

      call cpu_time(tr0)
      time0=time0+real(tr0-tr1,kind=8)

      call f_free(coeff_trans)

      !if (iproc==0) print*,'Time in gramschmidt_coeff',time0,time1,time2,time3,time4,time5,&
      !     time0+time1+time2+time3+time4+time5

  end do spin_loop

  call f_release_routine()


end subroutine gramschmidt_coeff_trans









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



subroutine check_taylor_order(error, max_error, order_taylor)
  use module_base
  use yaml_output
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: error, max_error
  integer,intent(inout) :: order_taylor

  ! Local variables
  character(len=12) :: act
  integer,parameter :: max_order_positive=50
  integer,parameter :: max_order_negative=-20
  logical :: is_ice

  if (order_taylor>=1000) then
      order_taylor=order_taylor-1000
      is_ice=.true.
  else
      is_ice=.false.
  end if

  if (order_taylor/=0) then
      ! only do this if approximations (Taylor or "negative thing") are actually used
      if (error<=1.d-1*max_error) then
          !! error is very small, so decrease the order of the polynomial
          !if (order_taylor>20) then
          !    ! always keep a minimum of 20
          !    act=' (decreased)'
          !    if (order_taylor>0) then
          !        order_taylor = floor(0.9d0*real(order_taylor,kind=8))
          !    else
          !        order_taylor = ceiling(0.9d0*real(order_taylor,kind=8))
          !    end if
          !end if
      else if (error>max_error) then
          ! error is too big, increase the order of the Taylor series by 10%
          act=' (increased)'
          if (order_taylor>0) then
              order_taylor = ceiling(max(1.1d0*real(order_taylor,kind=8),real(order_taylor+5,kind=8)))
          else
              order_taylor = floor(min(1.1d0*real(order_taylor,kind=8),real(order_taylor-5,kind=8)))
          end if
      else
          ! error is small enough, so do nothing
          act=' (unchanged)'
      end if
      !if (bigdft_mpi%iproc==0) call yaml_map('new Taylor order',trim(yaml_toa(order_taylor,fmt='(i0)'))//act)
  end if

  if (order_taylor>0) then
      if (order_taylor>max_order_positive) then
          order_taylor=max_order_positive
          if (bigdft_mpi%iproc==0) call yaml_warning('Taylor order reached maximum')
      end if
  else
      if (order_taylor<max_order_negative) then
          order_taylor=max_order_negative
          if (bigdft_mpi%iproc==0) call yaml_warning('Taylor order reached maximum')
      end if
  end if

  if (is_ice) then
      order_taylor=order_taylor+1000
  end if

end subroutine check_taylor_order
