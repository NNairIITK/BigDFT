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
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix
  use foe_base, only: foe_data
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



  if(.not.can_use_transposed) then
      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
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

  call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_(1), &
       psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)


  norm = f_malloc(orbs%norb,id='norm')
  call normalize_transposed(iproc, nproc, orbs, ovrlp%nspin, collcom, psit_c, psit_f, norm)

  call f_free(norm)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

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
           npsidim_orbs_small, lzd_small, hpsi_noprecond)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  use yaml_output
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: matrices_null, allocate_matrices, deallocate_matrices, sparsematrix_malloc, &
                               sparsematrix_malloc_ptr, DENSE_FULL, DENSE_MATMUL, SPARSE_FULL, SPARSEMM_SEQ, &
                               assignment(=)
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: uncompress_matrix, uncompress_matrix_distributed, compress_matrix_distributed, &
                          sequential_acces_matrix_fast, sparsemm, transform_sparse_matrix, orb_from_index
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

  ! Local variables
  integer :: iorb, jorb, ii, ii_trans, irow, jcol, info, lwork, jj, ispin, iseg, i
  integer :: isegstart, isegend, ierr
  real(kind=8) :: max_error, mean_error
  real(kind=8),dimension(:),allocatable :: tmp_mat_compr, hpsit_tmp_c, hpsit_tmp_f, hphi_nococontra
  integer,dimension(:),allocatable :: ipiv
  type(matrices),dimension(1) :: inv_ovrlp_
  real(8),dimension(:),allocatable :: inv_ovrlp_seq, lagmat_large
  real(8),dimension(:,:),allocatable :: lagmatp, inv_lagmatp
  integer,dimension(2) :: irowcol
  real(kind=8),dimension(:),pointer :: matrix_local
  integer,parameter :: GLOBAL_MATRIX=101, SUBMATRIX=102
  integer,parameter :: data_strategy_main=SUBMATRIX!GLOBAL_MATRIX

  call f_routine(id='orthoconstraintNonorthogonal')


  if(.not. can_use_transposed) then
      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.
  end if

  inv_ovrlp_seq = sparsematrix_malloc(linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_seq')
  inv_lagmatp = sparsematrix_malloc(linmat%l, iaction=DENSE_MATMUL, id='inv_lagmatp')
  lagmatp = sparsematrix_malloc(linmat%l, iaction=DENSE_MATMUL, id='lagmatp')
  inv_ovrlp_(1) = matrices_null()
  inv_ovrlp_(1)%matrix_compr = sparsematrix_malloc_ptr(linmat%l,iaction=SPARSE_FULL,id='inv_ovrlp_(1)%matrix_compr')

  !!if (calculate_inverse) then
      ! Invert the overlap matrix
      if (iproc==0) call yaml_map('calculation of S^-1','direct calculation')
      call overlapPowerGeneral(iproc, nproc, norder_taylor, 1, (/1/), -1, &
           imode=1, ovrlp_smat=linmat%s, inv_ovrlp_smat=linmat%l, &
           ovrlp_mat=linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp_, &
           check_accur=.true., max_error=max_error, mean_error=mean_error)
      !if (iproc==0) call yaml_map('max error',max_error)
      !if (iproc==0) call yaml_map('mean error',mean_error)
      !if (iproc==0) call yaml_scalar('no check taylor')
      call check_taylor_order(mean_error, max_inversion_error, norder_taylor)

  !!else

  !!    if (iproc==0) call yaml_map('calculation of S^-1','square of S^-1/2')
  !!    !@NEW instead of calculating S^-1, take S^-1/2 from memory and square it
  !!    call sequential_acces_matrix_fast(linmat%l, linmat%ovrlp_minusonehalf_%matrix_compr, inv_ovrlp_seq)
  !!    call uncompress_matrix_distributed(iproc, linmat%l, DENSE_MATMUL, linmat%ovrlp_minusonehalf_%matrix_compr, lagmatp)
  !!    call sparsemm(linmat%l, inv_ovrlp_seq, lagmatp, inv_lagmatp)
  !!    call compress_matrix_distributed(iproc, nproc, linmat%l, DENSE_MATMUL, &
  !!         inv_lagmatp, inv_ovrlp_%matrix_compr(linmat%l%isvctrp_tg+1:))
  !!end if


  ! Calculate <phi_alpha|g_beta>
  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat, lagmat_)

  lagmat_large = sparsematrix_malloc(linmat%l, iaction=SPARSE_FULL, id='lagmat_large')

  call symmetrize_matrix()


  ! Apply S^-1
  call sequential_acces_matrix_fast(linmat%l, inv_ovrlp_(1)%matrix_compr, inv_ovrlp_seq)
  ! Transform the matrix to the large sparsity pattern (necessary for the following uncompress_matrix_distributed)
  if (correction_orthoconstraint==0) then
      if (data_strategy_main==GLOBAL_MATRIX) then
          call transform_sparse_matrix(linmat%m, linmat%l, lagmat_%matrix_compr, lagmat_large, 'small_to_large')
      end if
      if (iproc==0) call yaml_map('correction orthoconstraint',.true.)
      call uncompress_matrix_distributed(iproc, linmat%l, DENSE_MATMUL, lagmat_large, lagmatp)
      call sparsemm(linmat%l, inv_ovrlp_seq, lagmatp, inv_lagmatp)
      call compress_matrix_distributed(iproc, nproc, linmat%l, DENSE_MATMUL, &
           inv_lagmatp, lagmat_large(linmat%l%isvctrp_tg+1:))
  end if
  if (data_strategy_main==SUBMATRIX) then
      call transform_sparse_matrix(linmat%m, linmat%l, lagmat_%matrix_compr, lagmat_large, 'large_to_small')
  end if
  call f_free(lagmat_large)
  call f_free(inv_ovrlp_seq)
  call f_free(lagmatp)
  call f_free(inv_lagmatp)
  call build_linear_combination_transposed(collcom, lagmat, lagmat_, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)

  ! The symmetrized Lagrange multiplier matrix has now the wrong sign
  call dscal(lagmat%nvctr*lagmat%nspin, -1.d0, lagmat_%matrix_compr(1), 1)


  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, hpsit_c, hpsit_f, lhphi, lzd)

  overlap_calculated=.false.
  

  ! @NEW apply S^-1 to the gradient
  hpsit_tmp_c = f_malloc(collcom%ndimind_c,id='psit_tmp_c')
  hpsit_tmp_f = f_malloc(7*collcom%ndimind_f,id='psit_tmp_f')
  hphi_nococontra = f_malloc(npsidim_orbs,id='hphi_nococontra')
  call build_linear_combination_transposed(collcom, linmat%l, inv_ovrlp_(1), hpsit_c, hpsit_f, .true., hpsit_tmp_c, hpsit_tmp_f, iproc)
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, hpsit_tmp_c, hpsit_tmp_f, hphi_nococontra, lzd)
  call large_to_small_locreg(iproc, npsidim_orbs_small, npsidim_orbs, lzd_small, lzd, &
       orbs, hphi_nococontra, hpsi_noprecond)
  ! END @NEW

  call deallocate_matrices(inv_ovrlp_(1))
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
          matrix_local = f_malloc_ptr(max(lagmat%nvctrp,1),id='matrix_local')
          tmp_mat_compr = sparsematrix_malloc(lagmat,iaction=SPARSE_FULL,id='tmp_mat_compr')
          call vcopy(lagmat%nvctr*lagmat%nspin, lagmat_%matrix_compr(1), 1, tmp_mat_compr(1), 1)
          do ispin=1,lagmat%nspin
              ishift=(ispin-1)*lagmat%nvctr
              if (isegend>=isegstart) then
                  !$omp parallel default(none) &
                  !$omp shared(isegstart,isegend,ishift,lagmat,matrix_local,tmp_mat_compr) &
                  !$omp private(iseg,ii,i,irowcol,ii_trans)
                  !$omp do
                  do iseg=isegstart,isegend
                      ii=lagmat%keyv(iseg)
                      ! A segment is always on one line, therefore no double loop
                      do i=lagmat%keyg(1,1,iseg),lagmat%keyg(2,1,iseg)
                         irowcol = orb_from_index(lagmat, i)
                         ii_trans = matrixindex_in_compressed(lagmat,lagmat%keyg(1,2,iseg),i)
                         matrix_local(ii-lagmat%isvctr) = -0.5d0*tmp_mat_compr(ii+ishift)-0.5d0*tmp_mat_compr(ii_trans+ishift)
                         ii=ii+1
                      end do
                  end do
                  !$omp end do
                  !$omp end parallel
              end if
              if (nproc>1) then
                  !!call mpi_allgatherv(matrix_local(1), lagmat%nvctrp, mpi_double_precision, &
                  !!     lagmat_%matrix_compr(ishift+1), lagmat%nvctr_par, lagmat%isvctr_par, mpi_double_precision, &
                  !!     bigdft_mpi%mpi_comm, ierr)
                  if (comm_strategy==ALLGATHERV) then
                      call mpi_allgatherv(matrix_local(1), lagmat%nvctrp, mpi_double_precision, &
                           lagmat_%matrix_compr(ishift+1), lagmat%nvctr_par, lagmat%isvctr_par, mpi_double_precision, &
                           bigdft_mpi%mpi_comm, ierr)
                      call f_free_ptr(matrix_local)
                  else if (comm_strategy==GET) then
                      !!call mpiget(iproc, nproc, bigdft_mpi%mpi_comm, lagmat%nvctrp, matrix_local, &
                      !!     lagmat%nvctr_par, lagmat%isvctr_par, lagmat%nvctr, lagmat_%matrix_compr(ishift+1:ishift+lagmat%nvctr))
                      call mpi_get_to_allgatherv(matrix_local(1), lagmat%nvctrp, &
                           lagmat_%matrix_compr(ishift+1), &
                           lagmat%nvctr_par, lagmat%isvctr_par, bigdft_mpi%mpi_comm)
                  else
                      stop 'symmetrize_matrix: wrong communication strategy'
                  end if
              else
                  call vcopy(lagmat%nvctr, matrix_local(1), 1, lagmat_%matrix_compr(ishift+1), 1)
              end if
              if (ispin==lagmat%nspin) call f_free_ptr(matrix_local)
          end do
      else if (data_strategy==SUBMATRIX) then
          ! Directly use the large sparsity pattern as this one is used later
          ! for the matrix vector multiplication
          tmp_mat_compr = sparsematrix_malloc(linmat%l,iaction=SPARSE_FULL,id='tmp_mat_compr')
          call transform_sparse_matrix(linmat%m, linmat%l, lagmat_%matrix_compr, tmp_mat_compr, 'small_to_large')
          do ispin=1,lagmat%nspin
              ishift=(ispin-1)*linmat%l%nvctr
              !!!!$omp parallel default(none) &
              !!!!$omp shared(linmat,lagmat_large,tmp_mat_compr,ishift) &
              !!!!$omp private(iseg,ii,i,ii_trans)
              !!!!$omp do
              !!!do iseg=linmat%l%istartendseg_t(1),linmat%l%istartendseg_t(2)
              !!!    ii=linmat%l%keyv(iseg)
              !!!    ! A segment is always on one line, therefore no double loop
              !!!    do i=linmat%l%keyg(1,1,iseg),linmat%l%keyg(2,1,iseg) !this is too much, but for the moment ok
              !!!        ii_trans = matrixindex_in_compressed(linmat%l,linmat%l%keyg(1,2,iseg),i)
              !!!        lagmat_large(ii+ishift) = -0.5d0*tmp_mat_compr(ii+ishift)-0.5d0*tmp_mat_compr(ii_trans+ishift)
              !!!        ii=ii+1
              !!!    end do
              !!!end do
              !!!!$omp end do
              !!!!$omp end parallel
              !!do itg=1,linmat%l%ntaskgroupp
              !!    iitg = linmat%l%inwhichtaskgroup(itg)
              !!    do iseg=1,linmat%l%nseg
              !!        ii=linmat%l%keyv(iseg)
              !!        if (ii+linmat%l%keyg(2,1,iseg)-linmat%l%keyg(1,1,iseg)<linmat%l%taskgroup_startend(1,1,iitg)) cycle
              !!        if (ii>linmat%l%taskgroup_startend(2,1,iitg)) exit
              !!        ! A segment is always on one line, therefore no double loop
              !!        do i=linmat%l%keyg(1,1,iseg),linmat%l%keyg(2,1,iseg) !this is too much, but for the moment ok
              !!            if (ii>=linmat%l%taskgroup_startend(1,1,iitg) .and.  ii<=linmat%l%taskgroup_startend(2,1,iitg)) then
              !!                ii_trans = matrixindex_in_compressed(linmat%l,linmat%l%keyg(1,2,iseg),i)
              !!                lagmat_large(ii+ishift) = -0.5d0*tmp_mat_compr(ii+ishift)-0.5d0*tmp_mat_compr(ii_trans+ishift)
              !!            end if
              !!            ii=ii+1
              !!        end do
              !!    end do
              !!end do

              !$omp parallel default(none) &
              !$omp shared(linmat,lagmat_large,tmp_mat_compr,ishift) &
              !$omp private(iseg,ii,i,ii_trans)
              !$omp do
              do iseg=linmat%l%iseseg_tg(1),linmat%l%iseseg_tg(2)
                  ii = linmat%l%keyv(iseg)
                  ! A segment is always on one line, therefore no double loop
                  do i=linmat%l%keyg(1,1,iseg),linmat%l%keyg(2,1,iseg) !this is too much, but for the moment ok
                      ii_trans = matrixindex_in_compressed(linmat%l,linmat%l%keyg(1,2,iseg),i)
                      lagmat_large(ii+ishift) = -0.5d0*tmp_mat_compr(ii+ishift)-0.5d0*tmp_mat_compr(ii_trans+ishift)
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


!> S^-1 exact only works for symmetric matrices
!! BOTH sparse matrices must be present together and inv_ovrlp should be nullified pointer, NOT inv_ovrlp_smat%matrix
!! when sparse matrices present, check is performed to see whether %matrix is allocated so that its allocated status remains unchanged
!! contents of %matrix not guaranteed to be correct though - inv_ovrlp_smat%can_use_dense set accordingly
!! power: -2 -> S^-1/2, 2 -> S^1/2, 1 -> S^-1
subroutine overlapPowerGeneral(iproc, nproc, iorder, ncalc, power, blocksize, imode, &
           ovrlp_smat, inv_ovrlp_smat, ovrlp_mat, inv_ovrlp_mat, check_accur, &
           max_error, mean_error, nspinx)
     !!foe_nseg, foe_kernel_nsegline, foe_istsegline, foe_keyg)
  use module_base
  use module_types
  use module_interfaces, except_this_one => overlapPowerGeneral
  use sparsematrix_base, only: sparse_matrix, &
                          sparsematrix_malloc_ptr, sparsematrix_malloc, sparsematrix_malloc0, sparsematrix_malloc0_ptr, &
                          assignment(=), &
                          SPARSE_FULL, DENSE_PARALLEL, DENSE_MATMUL, DENSE_FULL, SPARSEMM_SEQ
  use sparsematrix, only: compress_matrix, uncompress_matrix, transform_sparse_matrix, &
                          compress_matrix_distributed, uncompress_matrix_distributed, &
                          sequential_acces_matrix_fast, sparsemm
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, ncalc
  integer,dimension(ncalc),intent(in) :: power
  integer,intent(in) :: imode
  type(sparse_matrix),intent(inout) :: ovrlp_smat, inv_ovrlp_smat
  type(matrices),intent(inout) :: ovrlp_mat
  type(matrices),dimension(ncalc),intent(inout) :: inv_ovrlp_mat
  logical,intent(in) :: check_accur
  real(kind=8),intent(out),optional :: max_error, mean_error
  integer,intent(in),optional :: nspinx !< overwrite the default spin value
  
  ! Local variables
  integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i, its, maxits
  integer :: matrixindex_in_compressed, nmaxvalk
  real(kind=8), dimension(:,:), pointer :: inv_ovrlpp, ovrlppowerp
  real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp
  real(kind=8), dimension(:,:,:), pointer :: ovrlpminone, ovrlp_local, inv_ovrlp_local, ovrlppoweroldp, ovrlpminonep
  real(kind=8) :: factor, newfactor
  logical :: ovrlp_allocated, inv_ovrlp_allocated

  ! new for sparse taylor
  integer :: nout, nseq, ispin, ishift, isshift, ilshift, nspin
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
  real(kind=8) :: ex, max_error_p, mean_error_p


  !!write(*,*) 'iorder',iorder


  call f_routine(id='overlapPowerGeneral')
  call timing(iproc,'lovrlp^-1     ','ON')

  ! Overwrite the default spin value is present. Usefull to manipulate spinless
  ! matrices even in a polarized calculation
  if (present(nspinx)) then
      nspin=nspinx
  else
      nspin=ovrlp_smat%nspin
  end if


  if (iproc==0) then
      call yaml_newline()
      call yaml_mapping_open('calculate S^x')
      if (imode==SPARSE) then
          call yaml_map('mode','sparse')
      else if (imode==DENSE) then
          call yaml_map('mode','dense')
      end if
      !call yaml_map('power(1)',power(1))
      select case (power(1))
      case (-2)
          call yaml_map('x','-1/2')
      case (2)
          call yaml_map('x','1/2')
      case (1)
          call yaml_map('x','-1')
      case default
          stop 'wrong power(1)'
      end select
      call yaml_map('order',iorder)
  end if


  ! Perform a check of the arguments

  if (imode/=SPARSE .and. imode/=DENSE) stop 'wrong imode'

  if (imode==DENSE) then
      if (.not.associated(ovrlp_mat%matrix)) stop 'ovrlp_mat%matrix not associated'
      if (.not.associated(inv_ovrlp_mat(1)%matrix)) stop 'inv_ovrlp_mat(1)%matrix not associated'
  end if
  
  if (check_accur) then
      if (.not.present(max_error)) stop 'max_error not present'
      if (.not.present(mean_error)) stop 'mean_error not present'
  end if

  if (power(1)/=-2 .and. power(1)/=1 .and. power(1)/=2) stop 'wrong value of power(1)'

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
          call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr*nspin,ovrlp_mat%matrix(1,1,1),1,inv_ovrlp_mat(1)%matrix(1,1,1),1)
          if (power(1)==1) then
             if (blocksize<0) then
                 do ispin=1,nspin
                     call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_mat(1)%matrix(1,1,ispin))
                 end do
             else
                stop 'check if working - upper half may not be filled'
                call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                     ovrlp_smat%nfvctr, inv_ovrlp_mat(1)%matrix(1,1,1), ovrlp_smat%nfvctr)
                call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                     ovrlp_smat%nfvctr, inv_ovrlp_mat(1)%matrix(1,1,1), ovrlp_smat%nfvctr)
             end if
          else if (power(1)==2) then
              do ispin=1,nspin
                  call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                       blocksize,.true.,inv_ovrlp_mat(1)%matrix(1,1,ispin),inv_ovrlp_smat)
              end do
          else if (power(1)==-2) then
              do ispin=1,nspin 
                  call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                       blocksize,.false.,inv_ovrlp_mat(1)%matrix(1,1,ispin),inv_ovrlp_smat)
              end do
          end if
      else if (iorder<0) then
          ! sign approach as used in CP2K
          ! use 4 submatrices
          !if (nproc>1) then
              !Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat12p')
              !Amat21p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat21p')
              !Amat11p = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat11p')
          !else
              Amat12p = sparsematrix_malloc(ovrlp_smat,iaction=DENSE_PARALLEL,id='Amat12p')
              Amat21p = sparsematrix_malloc(ovrlp_smat,iaction=DENSE_PARALLEL,id='Amat21p')
              Amat11p = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='Amat11p')
          !end if
          ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
          Amat22p=>Amat11p
          !if (nproc>1) then
          !    Amat21=sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat21')
          !else
              Amat21=sparsematrix_malloc0(ovrlp_smat,iaction=DENSE_FULL,id='Amat21')
          !end if

          do ispin=1,nspin

              Amat12=>inv_ovrlp_mat(1)%matrix(:,:,ispin)

              call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr,ovrlp_mat%matrix(1,1,ispin),1,Amat12(1,1),1)
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
                  if(nproc > 1) then
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call timing(iproc,'lovrlp_comm   ','ON')
                      call mpi_allgatherv(Amat12p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, Amat12, &
                           ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                           mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                      call mpi_allgatherv(Amat21p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, Amat21, &
                           ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                           mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                      call timing(iproc,'lovrlp_comm   ','OF')
                      call timing(iproc,'lovrlp^-1     ','ON')
                  else
                      call vcopy(ovrlp_smat%nfvctr**2,Amat12p(1,1),1,Amat12(1,1),1)
                      call vcopy(ovrlp_smat%nfvctr**2,Amat21p(1,1),1,Amat21(1,1),1)
                  end if
              end do

              nullify(Amat22p)
              call f_free_ptr(Amat11p)

              if (power(1)==1) then
                  if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.0d0, Amat21(1,1), &
                       ovrlp_smat%nfvctr, Amat21p(1,1), ovrlp_smat%nfvctr, 0.0d0, Amat12p(1,1), ovrlp_smat%nfvctr)
                  if (nproc>1) then
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call timing(iproc,'lovrlp_comm   ','ON')
                      call mpi_allgatherv(Amat12p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, inv_ovrlp_mat(1)%matrix(1,1,ispin), &
                           ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                           mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                      call timing(iproc,'lovrlp_comm   ','OF')
                      call timing(iproc,'lovrlp^-1     ','ON')
                  else
                      call vcopy(ovrlp_smat%nfvctr**2, Amat12p(1,1), 1, inv_ovrlp_mat(1)%matrix(1,1,ispin), 1)
                  end if
              !else if (power(1)==2) then
              !   call vcopy(ovrlp_smat%nfvctr**2,Amat12(1,1),1,inv_ovrlp_mat(1)%matrix(1,1),1)
              else if (power(1)==-2) then
                  call vcopy(ovrlp_smat%nfvctr**2,Amat21(1,1),1,inv_ovrlp_mat(1)%matrix(1,1,ispin),1)
              end if

          end do

          call f_free(Amat12p)
          call f_free(Amat21p)
          nullify(Amat12)
          call f_free(Amat21)

      else
          if (iorder>1) then
              if (nproc>1) then
                  ovrlpminone => inv_ovrlp_mat(1)%matrix(:,:,:)
                  !ovrlpminonep = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='ovrlpminonep')
                  ovrlpminonep = f_malloc_ptr((/ovrlp_smat%nfvctr,max(ovrlp_smat%nfvctrp,1),nspin/),id='ovrlpminonep')
              else
                  ovrlpminone = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_FULL,id='ovrlpminone')
                  ovrlpminonep => ovrlpminone
              end if

              do ispin=1,nspin
                  if (ovrlp_smat%nfvctrp>0) call matrix_minus_identity_dense(ovrlp_smat%nfvctr,&
                                    ovrlp_smat%isfvctr,ovrlp_smat%nfvctrp, &
                                    ovrlp_mat%matrix(1,ovrlp_smat%isfvctr+1,ispin),ovrlpminonep(1,1,ispin))


                  !!if (iproc==0) write(*,*) 'isorb, ovrlp_mat%matrix(1,isorb+1,ispin)',isorb, ovrlp_mat%matrix(1,isorb+1,ispin)
                  !!do iorb=1,norbp
                  !!    do jorb=1,ovrlp_smat%nfvctr
                  !!        write(2800+10*iproc+ispin,'(a,3i8,3es14.6)') 'ispin, iorb, jorb, vals', &
                  !!             ispin, iorb, jorb, ovrlpminonep(jorb,iorb,ispin), ovrlp_mat%matrix(jorb,isorb+iorb,ispin)
                  !!    end do
                  !!end do


                  if(nproc > 1) then
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call timing(iproc,'lovrlp_comm   ','ON')
                      call mpi_allgatherv(ovrlpminonep(1,1,ispin), ovrlp_smat%nfvctr*norbp, &
                           mpi_double_precision, ovrlpminone(1,1,ispin), &
                           ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                           mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                      call timing(iproc,'lovrlp_comm   ','OF')
                      call timing(iproc,'lovrlp^-1     ','ON')
                  end if

              end do

              ovrlppoweroldp = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='ovrlppoweroldp')
              if (norbp>0) call vcopy(ovrlp_smat%nfvctr*norbp*nspin,ovrlpminonep(1,1,1),1,ovrlppoweroldp(1,1,1),1)

              if (nproc>1) then
                  call f_free_ptr(ovrlpminonep)
              else
                  nullify(ovrlpminonep)
              end if

          end if

          ovrlppowerp = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='ovrlppowerp')


          do ispin=1,nspin

              if (power(1)==1) then
                  factor=-1.0d0
              else if (power(1)==2) then
                  factor=0.5d0
              else if (power(1)==-2) then
                  factor=-0.5d0
              end if

              if (nproc>1) then
                  inv_ovrlpp = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='inv_ovrlpp')
              else
                  inv_ovrlpp => inv_ovrlp_mat(1)%matrix(:,:,ispin)
              end if

              if (norbp>0) call first_order_taylor_dense(ovrlp_smat%nfvctr,isorb,norbp,power(1),&
                  ovrlp_mat%matrix(1,isorb+1,ispin),inv_ovrlpp)
              !!do iorb=1,norbp
              !!    do jorb=1,ovrlp_smat%nfvctr
              !!        write(2900+10*iproc+ispin,'(a,3i8,3es14.6)') 'ispin, iorb, jorb, vals', &
              !!             ispin, iorb, jorb, inv_ovrlpp(jorb,iorb), ovrlp_mat%matrix(jorb,isorb+iorb,ispin)
              !!    end do
              !!end do

              do i=2,iorder
                  if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, &
                                    1.d0, ovrlpminone(1,1,ispin), &
                                    ovrlp_smat%nfvctr, ovrlppoweroldp(1,1,ispin), ovrlp_smat%nfvctr, &
                                    0.d0, ovrlppowerp(1,1), ovrlp_smat%nfvctr)
                 factor=newfactor(power(1),i,factor)
                  !!do iorb=1,norbp
                  !!    do jorb=1,ovrlp_smat%nfvctr
                  !!        write(3000+10*iproc+ispin,'(a,3i8,3es14.6)') 'ispin, iorb, jorb, vals', &
                  !!             ispin, iorb, jorb, ovrlppowerp(jorb,iorb), inv_ovrlpp(jorb,iorb), ovrlp_mat%matrix(jorb,isorb+iorb,ispin)
                  !!    end do
                  !!end do
                  call daxpy(ovrlp_smat%nfvctr*norbp,factor,ovrlppowerp,1,inv_ovrlpp,1)
                  !!do iorb=1,norbp
                  !!    do jorb=1,ovrlp_smat%nfvctr
                  !!        write(3100+10*iproc+ispin,'(a,4i8,3es14.6)') 'ispin, i, iorb, jorb, vals', &
                  !!             ispin, i, iorb, jorb, factor, ovrlppowerp(jorb,iorb), inv_ovrlpp(jorb,iorb)
                  !!    end do
                  !!end do
                  !!if (iproc==0) write(*,'(a,2i8,es16.9)') 'ispin, i, sum(inv_ovrlpp)', ispin, i, sum(inv_ovrlpp)
                  if (i/=iorder.and.norbp>0) call vcopy(ovrlp_smat%nfvctr*norbp,ovrlppowerp(1,1),1,ovrlppoweroldp(1,1,ispin),1)
              end do


              !!write(*,'(a,2i8,es15.6)') 'iproc, ispin, sum(inv_ovrlpp)', iproc, ispin, sum(inv_ovrlpp)
              if(nproc > 1) then
                  call timing(iproc,'lovrlp^-1     ','OF')
                  call timing(iproc,'lovrlp_comm   ','ON')
                  call mpi_allgatherv(inv_ovrlpp, ovrlp_smat%nfvctr*norbp, mpi_double_precision, &
                       inv_ovrlp_mat(1)%matrix(1,1,ispin), &
                       ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                       mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                  call timing(iproc,'lovrlp_comm   ','OF')
                  call timing(iproc,'lovrlp^-1     ','ON')
                  call f_free_ptr(inv_ovrlpp)
              end if
              !!if (iproc==0) write(*,'(a,2i8,es15.6)') 'iproc, ispin, sum(inv_ovrlp_mat(1)%matrix(:,:,ispin))', iproc, ispin, sum(inv_ovrlp_mat(1)%matrix(:,:,ispin))
          end do


          call f_free_ptr(ovrlppowerp)

          if (iorder>1) then
              if(nproc > 1) then
                  nullify(ovrlpminone)
              else
                  call f_free_ptr(ovrlpminone)
              end if
          end if

          if (iorder>1) then
              call f_free_ptr(ovrlppoweroldp)
          else
              nullify(inv_ovrlpp)
          end if
      end if

      if (check_accur) then
          do ispin=1,nspin
              call check_accur_overlap_minus_one(iproc,nproc,ovrlp_smat%nfvctr,ovrlp_smat%nfvctrp,ovrlp_smat%isfvctr,power(1),&
                   ovrlp_mat%matrix(:,:,ispin),inv_ovrlp_mat(1)%matrix(:,:,ispin),ovrlp_smat,max_error,mean_error)
              if (iproc==0) then
                  call yaml_newline()
                  if (nspin==1) then
                      call yaml_map('max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                  else
                      if (ispin==1) then
                          call yaml_map('spin up, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                      else if (ispin==2) then
                          call yaml_map('spin down, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                      end if
                  end if
              end if
          end do
      end if
  else if (imode==SPARSE) then
      if (iorder==0) then
          ovrlp_local = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='ovrlp_local')
          inv_ovrlp_local = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='inv_ovrlp_local')
          call timing(iproc,'lovrlp^-1     ','OF')
          call uncompress_matrix(iproc, ovrlp_smat, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_local)
          call timing(iproc,'lovrlp^-1     ','ON')
          do ispin=1,nspin
              !!write(*,*) 'sum(ovrlp_local(:,:,ispin))',sum(ovrlp_local(:,:,ispin))
              call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr,ovrlp_local(1,1,ispin),1,inv_ovrlp_local(1,1,ispin),1)
              if (power(1)==1) then
                 if (blocksize<0) then
                    call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_local(1,1,ispin))
                 else
                    stop 'check if working - upper half may not be filled'
                    call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                         ovrlp_smat%nfvctr, inv_ovrlp_local(1,1,ispin), ovrlp_smat%nfvctr)
                    call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                         ovrlp_smat%nfvctr, inv_ovrlp_local(1,1,ispin), ovrlp_smat%nfvctr)
                 end if
              else if (power(1)==2) then
                  call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                       blocksize,.true.,inv_ovrlp_local(1,1,ispin),inv_ovrlp_smat)
              else if (power(1)==-2) then
                  call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                       blocksize,.false.,inv_ovrlp_local(1,1,ispin),inv_ovrlp_smat)
              end if
              call timing(iproc,'lovrlp^-1     ','OF')
              call timing(iproc,'lovrlp^-1     ','ON')
              !!write(*,*) 'sum(inv_ovrlp_local(:,:,ispin))',sum(inv_ovrlp_local(:,:,ispin))
          end do
          call compress_matrix(iproc, inv_ovrlp_smat, inmat=inv_ovrlp_local, outmat=inv_ovrlp_mat(1)%matrix_compr)
          call f_free_ptr(ovrlp_local)
          call f_free_ptr(inv_ovrlp_local)
          ! #############################################################################
      else if (iorder<0) then ! could improve timing for checking, but for now just making sure it works
          ! use 4 submatrices
          if (nproc>0) then
              Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='Amat12p')
              Amat21p = sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='Amat21p')
              Amat11p = sparsematrix_malloc0_ptr(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='Amat11p')
          else
              Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat12p')
              Amat21p = sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat21p')
              Amat11p = sparsematrix_malloc0_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat11p')
          end if
          ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
          Amat22p=>Amat11p
          Amat12_compr=>inv_ovrlp_mat(1)%matrix_compr

          call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
               ovrlp_mat%matrix_compr, Amat12_compr, 'small_to_large')
          Amat12_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='Amat12_seq')

          do ispin=1,nspin

              ishift=(ispin-1)*inv_ovrlp_smat%nvctr

              call sequential_acces_matrix_fast(inv_ovrlp_smat, &
                   Amat12_compr(ishift+1:ishift+inv_ovrlp_smat%nvctr), Amat12_seq)
              call timing(iproc,'lovrlp^-1     ','OF')
              call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                   Amat12_compr(ishift+1:ishift+inv_ovrlp_smat%nvctr), Amat12p)
              call timing(iproc,'lovrlp^-1     ','ON')

              do iorb=1,inv_ovrlp_smat%smmm%nfvctrp
                  Amat21p(iorb+inv_ovrlp_smat%smmm%isfvctr,iorb)=1.0d0
              end do
              Amat21_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='Amat21_compr')
              call timing(iproc,'lovrlp^-1     ','OF')
              call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, &
                   Amat21p, Amat21_compr(inv_ovrlp_smat%isvctrp_tg+1:))
              call timing(iproc,'lovrlp^-1     ','ON')
              Amat21_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='Amat21_seq')
              call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat21_compr, Amat21_seq)

              ! calculate Xn+1=0.5*Xn*(3I-Xn**2)
              do its=1,abs(iorder)
                  call timing(iproc,'lovrlp^-1     ','OF')
                  call sparsemm(inv_ovrlp_smat, Amat12_seq, Amat21p, Amat11p)
                  call timing(iproc,'lovrlp^-1     ','ON')

                  if (inv_ovrlp_smat%smmm%nfvctrp>0) then
                      call vscal(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp,-0.5d0,Amat11p(1,1),1)
                  end if
                  !call vscal(ovrlp_smat%nfvctr*norbp,-0.5d0,Amat22p(1,1),1)
                  do iorb=1,inv_ovrlp_smat%smmm%nfvctrp
                      Amat11p(iorb+inv_ovrlp_smat%smmm%isfvctr,iorb)=Amat11p(iorb+inv_ovrlp_smat%smmm%isfvctr,iorb)+1.5d0
                  !    Amat22p(iorb+isorb,iorb)=Amat22p(iorb+isorb,iorb)+1.5d0
                  end do

                  call timing(iproc,'lovrlp^-1     ','OF')
                  call sparsemm(inv_ovrlp_smat, Amat12_seq, Amat22p, Amat12p)
                  call sparsemm(inv_ovrlp_smat, Amat21_seq, Amat11p, Amat21p)
                  call timing(iproc,'lovrlp^-1     ','ON')

                  if (its/=abs(iorder).or.power(1)/=2) then
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, &
                           Amat21p, Amat21_compr(inv_ovrlp_smat%isvctrp_tg+1:))
                      call timing(iproc,'lovrlp^-1     ','ON')
                  end if
                  if (its/=abs(iorder).or.power(1)==1) then
                      call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat21_compr, Amat21_seq)
                  end if
                  if (its/=abs(iorder).or.power(1)==2) then
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, &
                           Amat12p, Amat12_compr(inv_ovrlp_smat%isvctrp_tg+1:))
                      call timing(iproc,'lovrlp^-1     ','ON')
                  end if
                  if (its/=abs(iorder)) then
                      call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat12_compr, Amat12_seq)
                  end if
              end do

              call f_free(Amat12_seq)
              nullify(Amat22p)
              call f_free_ptr(Amat11p)

              if (power(1)==1) then
                  call timing(iproc,'lovrlp^-1     ','OF')
                  call sparsemm(inv_ovrlp_smat, Amat21_seq, Amat21p, Amat12p)
                  call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, Amat12p, &
                       inv_ovrlp_mat(1)%matrix_compr(ishift+inv_ovrlp_smat%isvctrp_tg+1:))
                  call timing(iproc,'lovrlp^-1     ','ON')
              !else if (power(1)==2) then
              !    call vcopy(inv_ovrlp_smat%nvctr,Amat12_compr(1),1,inv_ovrlp_smat%matrix_compr(1),1)
              else if (power(1)==-2) then
                  call vcopy(inv_ovrlp_smat%nvctr,Amat21_compr(1),1,inv_ovrlp_mat(1)%matrix_compr(ishift+1),1)
              end if

          end do

          nullify(Amat12_compr)
          call f_free(Amat21_compr)
          call f_free(Amat12p)
          call f_free(Amat21p)
          call f_free(Amat21_seq)

      else
          if (iorder<1000) then
              ovrlp_large_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='ovrlp_large_compr')
              call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
                   ovrlp_mat%matrix_compr, ovrlp_large_compr, 'small_to_large')

              ovrlpminonep = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='ovrlpminonep')
              invovrlpp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='invovrlpp')

              if (iorder>1) then
                  ovrlpminone_sparse_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, &
                       id='ovrlpminone_sparse_seq')
                  ovrlpminone_sparse = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=SPARSE_FULL, &
                       id='ovrlpminone_sparse')
                  ovrlpminoneoldp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='ovrlpminoneoldp')
              end if

              do ispin=1,nspin

                  isshift=(ispin-1)*ovrlp_smat%nvctr
                  ilshift=(ispin-1)*inv_ovrlp_smat%nvctr

                  if (iorder>1) then
                      !!ovrlpminone_sparse_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, &
                      !!     id='ovrlpminone_sparse_seq')
                      !!ovrlpminone_sparse = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=SPARSE_FULL, &
                      !!     id='ovrlpminone_sparse')
                      !!ovrlpminoneoldp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlpminoneoldp')

                      call matrix_minus_identity_sparse(ovrlp_smat%nfvctr, inv_ovrlp_smat, &
                           ovrlp_large_compr(isshift+1), ovrlpminone_sparse)
                      call sequential_acces_matrix_fast(inv_ovrlp_smat, ovrlpminone_sparse, ovrlpminone_sparse_seq)
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, ovrlpminone_sparse, ovrlpminoneoldp)
                      call timing(iproc,'lovrlp^-1     ','ON')

                      !!call f_free_ptr(ovrlpminone_sparse)

                      if (power(1)==1) then
                          factor=-1.0d0
                      else if (power(1)==2) then
                          factor=0.5d0
                      else if (power(1)==-2) then
                          factor=-0.5d0
                      end if
                  end if


                  if (inv_ovrlp_smat%smmm%nfvctrp>0) then
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                           ovrlp_large_compr, ovrlpminonep(:,:,1))
                      call timing(iproc,'lovrlp^-1     ','ON')
                      if (.not.check_accur) call f_free(ovrlp_large_compr)
                      call first_order_taylor_dense(inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%isfvctr, &
                           inv_ovrlp_smat%smmm%nfvctrp,power(1),ovrlpminonep,invovrlpp)
                  end if

                  do i=2,iorder
                      call timing(iproc,'lovrlp^-1     ','OF')
                      call sparsemm(inv_ovrlp_smat, ovrlpminone_sparse_seq, ovrlpminoneoldp, ovrlpminonep)
                      call timing(iproc,'lovrlp^-1     ','ON')
                      factor=newfactor(power(1),i,factor)
                      call daxpy(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp,factor,ovrlpminonep,1,invovrlpp,1)
                      if (i/=iorder.and.inv_ovrlp_smat%smmm%nfvctrp>0) then
                          call vcopy(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp,&
                          ovrlpminonep(1,1,1),1,ovrlpminoneoldp(1,1),1)
                      end if
                  end do
                  !!call to_zero(inv_ovrlp_smat%nvctr, inv_ovrlp_smat%matrix_compr(1))
                  call timing(iproc,'lovrlp^-1     ','OF')
                  call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, invovrlpp, &
                       inv_ovrlp_mat(1)%matrix_compr(ilshift+inv_ovrlp_smat%isvctrp_tg+1:))
                  call timing(iproc,'lovrlp^-1     ','ON')

              end do

              if (iorder>1) then
                  call f_free(ovrlpminone_sparse_seq)
                  call f_free(ovrlpminoneoldp)
                  call f_free_ptr(ovrlpminone_sparse)
              end if

              if (.not.check_accur) call f_free(invovrlpp)
              call f_free_ptr(ovrlpminonep)

          else

              ! @ NEW: ICE ##########################
              !!select case (power(1))
              !!case (-2)
              !!    ex=-0.5d0
              !!case (2)
              !!    ex=0.5d0
              !!case (1)
              !!    ex=-1.d0
              !!case default
              !!    stop 'wrong power(1)'
              !!end select
              call ice(iproc, nproc, iorder-1000, ovrlp_smat, inv_ovrlp_smat, power(1), ovrlp_mat, inv_ovrlp_mat(1))
              ! #####################################
          end if
      end if

      if (check_accur) then
          ! HERE STARTS LINEAR CHECK ##########################
          if (iorder<1 .or. iorder>=1000) then
              invovrlpp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='invovrlpp')
              ovrlp_large_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='ovrlp_large_compr')
              call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
                   ovrlp_mat%matrix_compr, ovrlp_large_compr, 'small_to_large')
          end if
          invovrlp_compr_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlp_large_compr_seq')
          ovrlp_largep = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='ovrlp_largep')

          do ispin=1,nspin
              isshift=(ispin-1)*ovrlp_smat%nvctr
              ilshift=(ispin-1)*inv_ovrlp_smat%nvctr
              call timing(iproc,'lovrlp^-1     ','OF')
              call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, ovrlp_large_compr(ilshift+1), ovrlp_largep)
              call timing(iproc,'lovrlp^-1     ','ON')
              call sequential_acces_matrix_fast(inv_ovrlp_smat, &
                   inv_ovrlp_mat(1)%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr), invovrlp_compr_seq)
              !!write(*,*) 'sum(inv_ovrlp_mat(1)%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr)', sum(inv_ovrlp_mat(1)%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr))

              if (power(1)==1) then
                  call check_accur_overlap_minus_one_sparse(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, &
                       inv_ovrlp_smat%smmm%nfvctrp, inv_ovrlp_smat%smmm%isfvctr, &
                       inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                       inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                       invovrlp_compr_seq, ovrlp_largep, power(1), &
                       max_error, mean_error)
              else if (power(1)==2) then
                  call timing(iproc,'lovrlp^-1     ','OF')
                  call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                       inv_ovrlp_mat(1)%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr), invovrlpp)
                  call timing(iproc,'lovrlp^-1     ','ON')
                  call check_accur_overlap_minus_one_sparse(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, &
                       inv_ovrlp_smat%smmm%nfvctrp, inv_ovrlp_smat%smmm%isfvctr, &
                       inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                       inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                       invovrlp_compr_seq, invovrlpp, power(1), &
                       max_error, mean_error, cmatp=ovrlp_largep)
              else if (power(1)==-2) then
                  ovrlp_compr_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlp_compr_seq') 
                  call sequential_acces_matrix_fast(inv_ovrlp_smat, ovrlp_large_compr(ilshift+1), ovrlp_compr_seq)
                  call timing(iproc,'lovrlp^-1     ','OF')
                  call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                       inv_ovrlp_mat(1)%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr), invovrlpp)
                  call timing(iproc,'lovrlp^-1     ','ON')
                  call check_accur_overlap_minus_one_sparse(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, &
                       inv_ovrlp_smat%smmm%nfvctrp, inv_ovrlp_smat%smmm%isfvctr, &
                       inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                       inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                       invovrlp_compr_seq, invovrlpp, power(1), &
                       max_error, mean_error, &
                       ovrlp_compr_seq)
                  call f_free(ovrlp_compr_seq)
              else
                  stop 'wrong power(1)'
              end if
              if (iproc==0) then
                  call yaml_newline()
                  if (nspin==1) then
                      call yaml_map('max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                  else
                      if (ispin==1) then
                          call yaml_map('spin up, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                      else if (ispin==2) then
                          call yaml_map('spin down, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                      end if
                  end if
              end if
          end do
          call f_free(invovrlp_compr_seq)
          call f_free(ovrlp_largep)
          call f_free(invovrlpp)
          call f_free(ovrlp_large_compr)
          !HERE ENDS LINEAR CHECK #############################
      end if
  end if sparse_dense

  if (iproc==0) then
      call yaml_mapping_close()
      call yaml_newline()
  end if

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
  real(kind=8),dimension(norb,norb) :: inv_ovrlp_half
  logical, intent(in) :: plusminus
  type(sparse_matrix),intent(in) :: smat

  integer :: info, iorb, jorb, ierr, iiorb, isorb, norbp, lwork, jjorb
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:),allocatable :: tempArr, orig_ovrlp
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_halfp
  real(kind=8),dimension(:,:), allocatable :: vr,vl ! for non-symmetric LAPACK
  real(kind=8),dimension(:),allocatable:: eval1 ! for non-symmetric LAPACK
  real(dp) :: temp, max_error, mean_error
  real(dp), allocatable, dimension(:) :: temp_vec
  logical, parameter :: symmetric=.true.
  logical, parameter :: check_lapack=.true.
  integer :: korb


  call f_routine(id='overlap_plus_minus_one_half_exact')

  !!!if (nproc>1) then
  !!    if (.not.present(smat)) then 
  !!        call f_err_throw('overlap_plus_minus_one_half_exact: for nproc>1, smat must be present!', &
  !!             err_name='BIGDFT_RUNTIME_ERROR')
  !!    end if
  !!!end if
           

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
        !!do jorb=1,norb
        !!    do korb=1,norb
        !!        write(910,'(a,2i8,es14.5)') 'jorb, korb, inv_ovrlp_half(korb,jorb)', jorb, korb, inv_ovrlp_half(korb,jorb)
        !!    end do
        !!end do
        call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, lwork, info)
        !!do jorb=1,norb
        !!    do korb=1,norb
        !!        write(920,'(a,2i8,es14.5)') 'jorb, korb, inv_ovrlp_half(korb,jorb)', jorb, korb, inv_ovrlp_half(korb,jorb)
        !!    end do
        !!end do
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
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, smat, max_error, mean_error)
           if (bigdft_mpi%iproc==0.and.abs(max_error)>1.0d-8) then
              print*,'LAPACK error for dsyev in overlap_plus_minus_one_half_exact',max_error
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
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, smat, max_error, mean_error)
           if (bigdft_mpi%iproc==0.and.abs(max_error)>1.0d-8) then
              print*,'LAPACK error for dgeev in overlap_plus_minus_one_half_exact',max_error
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

  call f_release_routine()

end subroutine overlap_plus_minus_one_half_exact


subroutine check_accur_overlap_minus_one_sparse(iproc, nproc, smat, norb, norbp, isorb, nseq, nout, &
           ivectorindex, onedimindices, amat_seq, bmatp, power, &
           max_error, mean_error, dmat_seq, cmatp)
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
  real(kind=8),intent(out) :: max_error, mean_error
  real(kind=8),dimension(nseq),intent(in),optional :: dmat_seq
  real(kind=8),dimension(norb,norbp),intent(in),optional :: cmatp

  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
  real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
  integer :: ierr, i,j

  call f_routine(id='check_accur_overlap_minus_one_sparse')

  tmpp=f_malloc0((/norb,norbp/),id='tmpp')
  if (power==1) then
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(smat, amat_seq, bmatp, tmpp)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, smat, max_error, mean_error)
  else if (power==2) then
      if (.not.present(cmatp)) stop 'cmatp not present'
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     !write(*,*) 'iproc, sum(amat_seq), sum(bmatp)', iproc, sum(amat_seq), sum(bmatp)
     call sparsemm(smat, amat_seq, bmatp, tmpp)
     call max_matrix_diff_parallel(iproc, norb, norbp, isorb, tmpp, cmatp, smat, max_error, mean_error)
     !max_error=0.5d0*max_error
     !mean_error=0.5d0*mean_error
  else if (power==-2) then
     if (.not.present(dmat_seq)) stop 'dmat_seq not present'
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(smat, amat_seq, bmatp, tmpp)
     tmp2p=f_malloc0((/norb,norbp/),id='tmp2p')
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
     !!     norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
     call sparsemm(smat, dmat_seq, tmpp, tmp2p)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, smat, max_error, mean_error)
     !max_error=0.5d0*max_error
     !mean_error=0.5d0*mean_error
     call f_free(tmp2p)
  else
     stop 'Error in check_accur_overlap_minus_one_sparse'
  end if
  call f_free(tmpp)

  call f_release_routine()

end subroutine check_accur_overlap_minus_one_sparse


subroutine check_accur_overlap_minus_one(iproc,nproc,norb,norbp,isorb,power,ovrlp,inv_ovrlp,&
           smat,max_error,mean_error)
  use module_base
  use sparsematrix_base, only: sparse_matrix
  implicit none
  integer,intent(in) :: iproc, nproc, norb, norbp, isorb, power
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp, inv_ovrlp
  type(sparse_matrix),intent(in) :: smat
  real(kind=8),intent(out) :: max_error, mean_error

  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
  real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
  integer :: ierr, i,j

  call f_routine(id='check_accur_overlap_minus_one')

  tmpp=f_malloc((/norb,norbp/),id='tmpp')
  if (power==1) then
     if (norbp>0) then
        call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
             norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     end if
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, smat, max_error, mean_error)
  else if (power==2) then
     if (norbp>0) then
        call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
             norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
         call max_matrix_diff_parallel(iproc, norb, norbp, isorb, tmpp, ovrlp(1,isorb+1), smat, max_error, mean_error)
         max_error=0.5d0*max_error
         mean_error=0.5d0*mean_error
     else
         max_error=0.d0
         mean_error=0.d0
     end if
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
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, smat, max_error, mean_error)
     max_error=0.5d0*max_error
     mean_error=0.5d0*mean_error
     call f_free(tmp2p)
  else
     stop 'Error in check_accur_overlap_minus_one'
  end if
  call f_free(tmpp)

  call f_release_routine()

end subroutine check_accur_overlap_minus_one


subroutine max_matrix_diff(iproc, norb, mat1, mat2, smat, max_deviation, mean_deviation)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: mat1, mat2
  type(sparse_matrix),intent(in) :: smat
  real(8),intent(out):: max_deviation, mean_deviation

  ! Local variables
  integer:: iorb, jorb, ind
  real(8):: error, num

  call timing(iproc,'dev_from_unity','ON') 
  max_deviation=0.d0
  mean_deviation=0.d0
  num=0.d0
  do iorb=1,norb
     do jorb=1,norb
        ind=matrixindex_in_compressed(smat,jorb,iorb)
        if (ind>0) then
            error=(mat1(jorb,iorb)-mat2(jorb,iorb))**2
            max_deviation=max(error,max_deviation)
            mean_deviation=mean_deviation+error
            num=num+1.d0
        end if
     end do
  end do
  mean_deviation=mean_deviation/num
  mean_deviation=sqrt(mean_deviation)
  max_deviation=sqrt(max_deviation)
  call timing(iproc,'dev_from_unity','OF') 

end subroutine max_matrix_diff


subroutine max_matrix_diff_parallel(iproc, norb, norbp, isorb, mat1, mat2, &
           smat, max_deviation, mean_deviation)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb, norbp, isorb
  real(8),dimension(norb,norbp),intent(in):: mat1, mat2
  type(sparse_matrix),intent(in) :: smat
  real(8),intent(out):: max_deviation, mean_deviation

  ! Local variables
  integer:: iorb, iiorb, jorb, ierr, ind
  real(8):: error, num
  real(kind=8),dimension(2) :: reducearr

  call timing(iproc,'dev_from_unity','ON') 
  max_deviation=0.d0
  mean_deviation=0.d0
  num=0.d0
  do iorb=1,norbp
     iiorb=iorb+isorb
     !$omp parallel default(private) shared(iorb, iiorb, norb, mat1, mat2, max_deviation, mean_deviation, num, smat)
     !$omp do reduction(max:max_deviation) reduction(+:mean_deviation,num)
     do jorb=1,norb
        ind=matrixindex_in_compressed(smat,jorb,iiorb)
        if (ind>0) then
            ! This entry is within the sparsity pattern, i.e. it matters for the error.
            error=(mat1(jorb,iorb)-mat2(jorb,iorb))**2
            max_deviation=max(error,max_deviation)
            mean_deviation=mean_deviation+error
            num=num+1.d0
        end if
     end do
     !$omp end do
     !$omp end parallel
  end do

  if (bigdft_mpi%nproc>1) then
      reducearr(1)=mean_deviation
      reducearr(2)=num
      call mpiallred(max_deviation, 1, mpi_max, bigdft_mpi%mpi_comm)
      call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
      mean_deviation=reducearr(1)
      num=reducearr(2)
  end if

  mean_deviation=mean_deviation/num
  mean_deviation=sqrt(mean_deviation)
  max_deviation=sqrt(max_deviation)

  call timing(iproc,'dev_from_unity','OF') 

end subroutine max_matrix_diff_parallel


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


subroutine deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, ovrlp, smat, max_deviation, mean_deviation)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb, norbp, isorb
  real(8),dimension(norb,norbp),intent(in):: ovrlp
  type(sparse_matrix),intent(in) :: smat
  real(8),intent(out):: max_deviation, mean_deviation

  ! Local variables
  integer:: iorb, iiorb, jorb, ierr, ind
  real(8):: error, num
  real(kind=8),dimension(2) :: reducearr


  call timing(iproc,'dev_from_unity','ON') 
  max_deviation=0.d0
  mean_deviation=0.d0
  num=0.d0
  do iorb=1,norbp
     iiorb=iorb+isorb
     !$omp parallel default(private) shared(norb, iiorb, ovrlp, iorb, max_deviation, mean_deviation, num, smat)
     !$omp do reduction(max:max_deviation) reduction(+:mean_deviation,num)
     do jorb=1,norb
        ind=matrixindex_in_compressed(smat,jorb,iiorb)
        if (ind>0) then
            ! This entry is within the sparsity pattern, i.e. it matters for the error.
            if(iiorb==jorb) then
               error=(ovrlp(jorb,iorb)-1.d0)**2
            else
               error=ovrlp(jorb,iorb)**2
            end if
            max_deviation=max(error,max_deviation)
            mean_deviation=mean_deviation+error
            num=num+1.d0
        end if
     end do
     !$omp end do
     !$omp end parallel
  end do
  if (nproc>1) then
      reducearr(1)=mean_deviation
      reducearr(2)=num
      call mpiallred(max_deviation, 1, mpi_max, bigdft_mpi%mpi_comm)
      call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
      mean_deviation=reducearr(1)
      num=reducearr(2)
  end if

  mean_deviation=mean_deviation/num
  mean_deviation=sqrt(mean_deviation)
  max_deviation=sqrt(max_deviation)

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
  integer(kind=8) :: ii, iend
  integer :: i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
  integer :: iiorb, ierr, iseg, ind, ishift, ispin
  real(kind=8) :: error
  real(kind=8),dimension(:,:),pointer :: ovrlp_tmp, ovrlp_tmp_inv_half
  logical,dimension(:),allocatable :: in_neighborhood
  character(len=*),parameter :: subname='overlap_power_minus_one_half_parallel'
  !type(matrices) :: inv_ovrlp_half_
  !!integer :: itaskgroups, iitaskgroup, imin, imax

  !!imin=ovrlp%nvctr
  !!imax=0
  !!do itaskgroups=1,ovrlp%ntaskgroupp
  !!    iitaskgroup = ovrlp%inwhichtaskgroup(itaskgroups)
  !!    imin = min(imin,ovrlp%taskgroup_startend(1,1,iitaskgroup))
  !!    imax = max(imax,ovrlp%taskgroup_startend(2,1,iitaskgroup))
  !!end do


  call timing(iproc,'lovrlp^-1/2par','ON')
  call f_routine('overlap_power_minus_one_half_parallel')

  in_neighborhood = f_malloc(ovrlp%nfvctr,id='in_neighborhood')

  !inv_ovrlp_half_ = matrices_null()
  !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_)
  call to_zero(inv_ovrlp_half%nvctr*inv_ovrlp_half%nspin, inv_ovrlp_half_%matrix_compr(1))

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

      ishift=(ispin-1)*ovrlp%nvctr

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
               n=n+1
            end do
            iseg=iseg+1
            if (iseg>ovrlp%nseg) exit
            ii = int((ovrlp%keyg(1,2,iseg)-1),kind=8)*int(ovrlp%nfvctr,kind=8) + &
                 int(ovrlp%keyg(1,1,iseg),kind=8)
            if (ii>iend) exit
         end do

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
                  ind=ind+ishift
                  ovrlp_tmp(kkorb,jjorb)=ovrlp_mat%matrix_compr(ind)
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
                     ind=ind+ishift
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
          call mpiallred(inv_ovrlp_half_%matrix_compr(1), inv_ovrlp_half%nvctr*inv_ovrlp_half%nspin, mpi_sum, bigdft_mpi%mpi_comm)
      end if

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
           orbs, ovrlp, ovrlp_, inv_ovrlp_half, inv_ovrlp_half_(1))
      !call deallocate_matrices(ovrlp_)
  else
      nullify(inv_ovrlp_null)
      ! do sparse.. check later
      !ovrlp%matrix_compr=ovrlp_%matrix_compr
      call overlapPowerGeneral(iproc, nproc, methTransformOverlap, 1, (/-2/), &
           orthpar%blocksize_pdsyev, &
           imode=1, check_accur=.true., &
           ovrlp_mat=ovrlp_, inv_ovrlp_mat=inv_ovrlp_half_, &
           ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
           max_error=max_error, mean_error=mean_error)
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
  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

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
  !@WARNING CHECK THIS
  !!ovrlp_%matrix_compr = inv_ovrlp_half%matrix_compr
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
        call to_zero(corb,ovrlp_coeff(1,1))
     end if

     if (nproc>1) then
        call mpiallred(ovrlp_coeff(1,1), corb, mpi_sum, bigdft_mpi%mpi_comm)
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

subroutine gramschmidt_coeff_trans(iproc,nproc,norb,basis_orbs,basis_overlap,basis_overlap_mat,coeff)
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

  coeff_transp=f_malloc((/norb,basis_orbs%norbp/),id='coeff_transp')
  !$omp parallel do default(private) shared(coeff,coeff_transp,norb,basis_orbs)
  do iorb=1,norb
     do jtmb=1,basis_orbs%norbp
        coeff_transp(iorb,jtmb) = coeff(jtmb+basis_orbs%isorb,iorb)
     end do
  end do
  !$omp end parallel do

  call cpu_time(tr1)
  time0=time0+real(tr1-tr0,kind=8)

  ! orthonormalizing all iorb<corb wrt corb (assume original vectors were normalized)
  do corb=norb,1,-1

     call cpu_time(tr0)

     coeff_tmp=f_malloc((/basis_orbs%norb,1/),id='coeff_tmp')
     ! calculate relevant part of cSc
     if (basis_orbs%norbp>0) then
        call dgemm('n', 't', basis_orbs%norb, 1, basis_orbs%norbp, 1.d0, &
             basis_overlap_mat%matrix(1,basis_orbs%isorb+1,1), basis_orbs%norb, &
             coeff_transp(corb,1), norb, 0.d0, &
             coeff_tmp(1,1), basis_orbs%norb)
     else
        call to_zero(corb,coeff_tmp(1,1))
     end if

     if (nproc>1) then
        call mpiallred(coeff_tmp(1,1), basis_orbs%norb, mpi_sum, bigdft_mpi%mpi_comm)
     end if

     call cpu_time(tr1)
     time1=time1+real(tr1-tr0,kind=8)

     ovrlp_coeff=f_malloc((/corb,1/),id='ovrlp_coeff')
     if (basis_orbs%norbp>0) then
        call dgemm('n', 'n', corb, 1, basis_orbs%norbp, 1.d0, &
             coeff_transp(1,1), norb, &
             coeff_tmp(1+basis_orbs%isorb,1), basis_orbs%norb, 0.d0, &
             ovrlp_coeff(1,1), corb)
     else
        call to_zero(corb,ovrlp_coeff(1,1))
     end if

     if (nproc>1) then
        call mpiallred(ovrlp_coeff(1,1), corb, mpi_sum, bigdft_mpi%mpi_comm)
     end if
     call f_free(coeff_tmp)

     call cpu_time(tr0)
     time2=time2+real(tr0-tr1,kind=8)

     ! (c_corb S c_iorb) * c_corb
     coeff_tmp=f_malloc((/corb,basis_orbs%norbp/),id='coeff_tmp')
     if (basis_orbs%norbp>0) then
        call dgemm('n', 'n', corb-1, basis_orbs%norbp, 1, 1.d0, &
             ovrlp_coeff(1,1), corb, &
             coeff_transp(corb,1), norb, 0.d0, &
             coeff_tmp(1,1), corb)
     end if

     call cpu_time(tr1)
     time3=time3+real(tr1-tr0,kind=8)

     ! sum and transpose coeff for allgatherv
     !$omp parallel do default(private) shared(coeff_transp,coeff_tmp,corb,basis_orbs,ovrlp_coeff)
     do iorb=1,corb-1
        do jtmb=1,basis_orbs%norbp
           coeff_transp(iorb,jtmb) = coeff_transp(iorb,jtmb) - coeff_tmp(iorb,jtmb)/ovrlp_coeff(corb,1)
        end do
     end do
     !$omp end parallel do

     do jtmb=1,basis_orbs%norbp
        coeff_transp(corb,jtmb) = coeff_transp(corb,jtmb)/sqrt(ovrlp_coeff(corb,1))
     end do

     call cpu_time(tr0)
     time4=time4+real(tr0-tr1,kind=8)

     call f_free(ovrlp_coeff)
     call f_free(coeff_tmp)
  end do

  call cpu_time(tr0)

  coeff_trans=f_malloc((/norb,basis_orbs%norb/),id='coeff_tmp')
  if(nproc > 1) then
     call mpi_allgatherv(coeff_transp(1,1), basis_orbs%norbp*norb, mpi_double_precision, coeff_trans(1,1), &
        norb*basis_orbs%norb_par(:,0), norb*basis_orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call vcopy(basis_orbs%norbp*norb,coeff_transp(1,1),1,coeff_trans(1,1),1)
  end if
  call f_free(coeff_transp)

  call cpu_time(tr1)
  time5=time5+real(tr1-tr0,kind=8)

  ! untranspose coeff
  !$omp parallel do default(private) shared(coeff,coeff_trans,norb,basis_orbs)
  do jtmb=1,basis_orbs%norb
     do iorb=1,norb
        coeff(jtmb,iorb) = coeff_trans(iorb,jtmb)
     end do
  end do
  !$omp end parallel do

  call cpu_time(tr0)
  time0=time0+real(tr0-tr1,kind=8)

  call f_free(coeff_trans)

  !if (iproc==0) print*,'Time in gramschmidt_coeff',time0,time1,time2,time3,time4,time5,&
  !     time0+time1+time2+time3+time4+time5

  call f_release_routine()


end subroutine gramschmidt_coeff_trans





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
          ! A segment is always on one line, therefore no double loop
         do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
             iiorb = smat%keyg(1,2,iseg)
             jjorb = jorb
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
           smat, max_error, mean_error)
  use module_base
  use module_types
  use module_interfaces, except_this_one => overlap_minus_one_half_serial
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, power, norb
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp_matrix
  real(kind=8),dimension(:,:),pointer,intent(out) :: inv_ovrlp_matrix
  type(sparse_matrix),intent(in) :: smat
  logical,intent(in) :: check_accur
  real(kind=8),intent(out),optional :: max_error, mean_error
  
  ! Local variables
  integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i, its, maxits
  integer :: matrixindex_in_compressed, nmaxvalk
  real(kind=8), dimension(:,:), pointer :: ovrlpminonep, ovrlpminone, inv_ovrlpp, ovrlppowerp, ovrlppoweroldp
  real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp, ovrlp_local, inv_ovrlp_local
  real(kind=8) :: factor, newfactor
  logical :: ovrlp_allocated, inv_ovrlp_allocated

  ! new for sparse taylor
  integer :: nout, nseq
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


  call f_routine(id='overlap_minus_one_half_serial')
  call timing(iproc,'lovrlp^-1     ','ON')

  if (nproc>1) then
      stop 'this routine only works in serial'
  end if

  
  if (check_accur) then
      if (.not.present(max_error)) stop 'max_error not present'
      if (.not.present(mean_error)) stop 'mean_error not present'
  end if

  if (power/=-2 .and. power/=1 .and. power/=2) stop 'wrong value of power'


      if (iorder==0) then
          call vcopy(norb*norb,ovrlp_matrix(1,1),1,inv_ovrlp_matrix(1,1),1)
          if (power==1) then
             call overlap_minus_one_exact_serial(norb,inv_ovrlp_matrix)
          else if (power==2) then
             call overlap_plus_minus_one_half_exact(1,norb,blocksize,.true.,inv_ovrlp_matrix,smat)
          else if (power==-2) then
             call overlap_plus_minus_one_half_exact(1,norb,blocksize,.false.,inv_ovrlp_matrix,smat)
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
          call check_accur_overlap_minus_one(iproc,nproc,norb,norb,0,power,ovrlp_matrix,inv_ovrlp_matrix,&
               smat, max_error,mean_error)
      end if


  call f_release_routine()
  call timing(iproc,'lovrlp^-1     ','OF')


end subroutine overlap_minus_one_half_serial


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
!!  call to_zero(orbs%norb, orbs%eval(1))
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
