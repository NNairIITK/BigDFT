subroutine calculate_residue_ks(iproc, nproc, num_extra, ksorbs, tmb, hpsit_c, hpsit_f)
  use module_base
  use module_types
  !use module_interfaces, fake_name => calculate_residue_ks,fake_B=>calculate_energy_and_gradient_linear
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix, &
                               matrices_null, allocate_matrices, deallocate_matrices,copy_sparse_matrix
  use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                          extract_taskgroup_inplace, uncompress_matrix2
  use transposed_operations, only: calculate_overlap_transposed
  use coeffs, only: calculate_kernel_and_energy
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, num_extra
  type(dft_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: ksorbs
  real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f

  integer :: iorb, istat, iall, ierr
  real(kind=8) :: ksres_sum
  real(kind=8), dimension(:), allocatable :: ksres
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, grad_coeff
  type(sparse_matrix) :: grad_ovrlp
  type(matrices) :: grad_ovrlp_
  character(len=256) :: subname='calculate_residue_ks'


  call f_routine(id='calculate_residue_ks')

  ! want to calculate the residue of the KS states here, not just the tmbs
  ! for now just occupied, eventually would want occupied+num_extra
  ! probably better to calculate c_i^a|g_a> first but copying and pasting for now (INEFFICIENT but just testing)
  !!if(associated(hpsit_c)) then
  !!    iall=-product(shape(hpsit_c))*kind(hpsit_c)
  !!    deallocate(hpsit_c, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_c', subname)
  !!end if
  !!if(associated(hpsit_f)) then
  !!    iall=-product(shape(hpsit_f))*kind(hpsit_f)
  !!    deallocate(hpsit_f, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_f', subname)
  !!end if
  !!allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
  !!call memocc(istat, hpsit_c, 'hpsit_c', subname)
  !!allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
  !!call memocc(istat, hpsit_f, 'hpsit_f', subname)

  ! should already be done in orthoconstraintnonorthogonal so can use directly
  !call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !     tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
  !!can_use_transposed=.true.

  !call nullify_sparse_matrix(grad_ovrlp)
  grad_ovrlp=sparse_matrix_null()
  !call sparse_copy_pattern(tmb%linmat%smat(2), grad_ovrlp, iproc, subname)
  call copy_sparse_matrix(tmb%linmat%smat(2), grad_ovrlp)
  !grad_ovrlp%matrix_compr=f_malloc_ptr(grad_ovrlp%nvctr,id='grad_ovrlp%matrix_compr')
  grad_ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%smat(2), allocate_full=.false., &
       matname='grad_ovrlp_', mat=grad_ovrlp_)

  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, hpsit_c, hpsit_c, &
       hpsit_f, hpsit_f, tmb%linmat%smat(2), tmb%linmat%auxm, grad_ovrlp_)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), grad_ovrlp_)
  !! This can then be deleted if the transition to the new type has been completed.
  !grad_ovrlp%matrix_compr=grad_ovrlp_%matrix_compr


  grad_coeff = f_malloc((/ tmb%orbs%norb, tmb%orbs%norb /),id='grad_coeff')
  coeff_tmp = f_malloc((/ tmb%orbs%norbp, max(tmb%orbs%norb, 1) /),id='coeff_tmp')

  !grad_ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='grad_ovrlp%matrix')
  call uncompress_matrix2(iproc,nproc,bigdft_mpi%mpi_comm, &
       grad_ovrlp,grad_ovrlp_%matrix_compr,grad_ovrlp_%matrix)

  ! can change this so only go up to ksorbs%norb...
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, 1.d0, grad_ovrlp_%matrix(tmb%orbs%isorb+1,1,1), &
          tmb%orbs%norb, tmb%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp, tmb%orbs%norbp)
     call dgemm('t', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norbp, 1.d0, tmb%coeff(tmb%orbs%isorb+1,1), &
          tmb%orbs%norb, coeff_tmp, tmb%orbs%norbp, 0.d0, grad_coeff, tmb%orbs%norb)
  else
     call f_zero(grad_coeff)
  end if

  call f_free(coeff_tmp)
  call deallocate_matrices(grad_ovrlp_)

  if (nproc>1) then
      call mpiallred(grad_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if

  ! now calculate sqrt(<g_i|g_i>) and mean value
  ksres = f_malloc(ksorbs%norb+num_extra,id='ksres')
  
  ksres_sum=0.0d0
  do iorb=1,ksorbs%norb+num_extra
    ksres(iorb)=dsqrt(grad_coeff(iorb,iorb))
    !ksres_sum=ksres_sum+ksres(iorb)
    ksres_sum=ksres_sum+grad_coeff(iorb,iorb)
    if (iproc==0) write(*,*) 'KS residue',iorb,ksres(iorb)!,tmb%orbs%occup(iorb)
  end do
  if (iproc==0) write(*,*) 'Average KS residue',sqrt(ksres_sum/real(ksorbs%norb+num_extra,gp))


  !call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
  !     tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ham)
  ! calculate Tr[Kg]  (not recalculating kernel as don't have the correct occs here)
  !call calculate_kernel_and_energy(iproc,nproc,denskern,grad_coeff,ksres_sum,tmb%coeff,orbs,tmb%orbs,.true.)
  grad_ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%smat(2), allocate_full=.false., matname='grad_ovrlp_', mat=grad_ovrlp_)
  !grad_ovrlp_%matrix_compr=grad_ovrlp%matrix_compr
  !!call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
  call extract_taskgroup_inplace(grad_ovrlp, grad_ovrlp_)
  call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),grad_ovrlp,&
       tmb%linmat%kernel_, grad_ovrlp_, &
       ksres_sum,tmb%coeff,tmb%orbs%norbp, tmb%orbs%isorb, tmb%orbs%norbu, tmb%orbs%norb, tmb%orbs%occup, .false.)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, grad_ovrlp, grad_ovrlp_)
  call deallocate_matrices(grad_ovrlp_)
  if (iproc==0) write(*,*) 'KS residue from trace',&
       dsqrt(ksres_sum)/real(tmb%orbs%norb,gp) ! should update normalization as would only be occ here not extra?

  call deallocate_sparse_matrix(grad_ovrlp)

  call f_free(grad_coeff)
  call f_free(ksres)

  call f_release_routine()

end subroutine calculate_residue_ks


