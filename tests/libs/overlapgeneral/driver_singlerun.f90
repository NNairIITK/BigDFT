program driver_singlerun
  use bigdft_run
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               matrices_null, deallocate_sparse_matrix, deallocate_matrices, &
                               assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
  use sparsematrix_init, only: read_ccs_format, ccs_to_sparsebigdft, ccs_values_to_bigdft, &
                               read_bigdft_format, bigdft_to_sparsebigdft
  use sparsematrix, only: write_matrix_compressed, check_symmetry, &
                          write_sparsematrix_CCS, write_sparsematrix
  use module_interfaces, only: overlapPowerGeneral
  implicit none

  ! Variables
  integer :: iproc, nproc, ncol, nnonzero, nseg, ncolp, iscol, ierr
  character(len=*),parameter :: filename='matrix.dat'
  integer,dimension(:),pointer :: col_ptr, row_ind, keyv
  integer,dimension(:,:,:),pointer :: keyg
  real(kind=8),dimension(:),pointer :: val
  type(sparse_matrix) :: smat
  type(matrices) :: matA
  type(matrices),dimension(1) :: matB
  real(kind=8) :: max_error, mean_error
  logical :: symmetric
  real(kind=8) :: time_start, time_end

  ! Initialize
  call f_lib_initialize()
  call bigdft_init()!mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  iproc=bigdft_mpi%iproc!mpi_info(1)
  nproc=bigdft_mpi%nproc!mpi_info(2)

  ! Nullify all pointers
  nullify(col_ptr)
  nullify(row_ind)
  nullify(keyv)
  nullify(keyg)
  nullify(val)

  
  ! Read in a file in the CCS format
  !call read_ccs_format(filename, ncol, nnonzero, col_ptr, row_ind, val)

  ! Read in a file in the BigDFT format
  call read_bigdft_format(filename, ncol, nnonzero, nseg, keyv, keyg, val)

  ! Create the corresponding BigDFT sparsity pattern
  !call ccs_to_sparsebigdft(iproc, nproc, ncol, ncol, 0, nnonzero, row_ind, col_ptr, smat)
  call distribute_columns_on_processes(iproc, nproc, ncol, ncolp, iscol)
  write(*,'(a,4i9)') 'iproc, ncol, ncolp, iscol', iproc, ncol, ncolp, iscol
  call bigdft_to_sparsebigdft(iproc, nproc, ncol, ncolp, iscol, nnonzero, nseg, keyg, smat)

  matA = matrices_null()
  matB(1) = matrices_null()

  ! Assign the values
  matA%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matA%matrix_compr')
  !call ccs_values_to_bigdft(ncol, nnonzero, row_ind, col_ptr, smat, val, matA)
  matA%matrix_compr = val

  ! Check the symmetry
  symmetric = check_symmetry(ncol, smat)
  if (.not.symmetric) stop 'ERROR not symmetric'
  
  ! Write the original matrix
  !if (iproc==0) call write_matrix_compressed('Original matrix', smat, matA)
  if (iproc==0) call write_sparsematrix_CCS('original_css.dat', smat, matA)
  if (iproc==0) call write_sparsematrix('original_bigdft.dat', smat, matA)

  ! Calculate the inverse
  call mpibarrier(bigdft_mpi%mpi_comm)
  time_start = mpi_wtime()
  matB(1)%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matB(1)%matrix_compr')
  call overlapPowerGeneral(iproc, nproc, 1020, 1, (/1/), -8, &
       1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
       check_accur=.true., max_error=max_error, mean_error=mean_error)
  !!call overlapPowerGeneral(iproc, nproc, 0, 1, (/1/), -8, &
  !!     1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
  !!     check_accur=.true., max_error=max_error, mean_error=mean_error)
  call mpibarrier(bigdft_mpi%mpi_comm)
  time_end = mpi_wtime()
  if (iproc==0) write(*,*) 'walltime',time_end-time_start

  ! Write the results
  !if (iproc==0) call write_matrix_compressed('Result', smat, matB(1))

  if (iproc==0) call write_sparsematrix_CCS('result_css.dat', smat, matB(1))
  if (iproc==0) call write_sparsematrix('result_bigdft.dat', smat, matB(1))

  ! Deallocations
  call deallocate_sparse_matrix(smat)
  call deallocate_matrices(matA)
  call deallocate_matrices(matB(1))
  call f_free_ptr(col_ptr)
  call f_free_ptr(row_ind)
  call f_free_ptr(keyv)
  call f_free_ptr(keyg)
  call f_free_ptr(val)

  call bigdft_finalize(ierr)

  call f_lib_finalize()

end program driver_singlerun


subroutine distribute_columns_on_processes(iproc, nproc, ncol, ncolp, iscol)
  implicit none
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, ncol
  integer,intent(out) :: ncolp, iscol

  ! Local variables
  integer :: ncolpx, ii, i, jproc
  real(kind=8) :: tt

  ! Determine the number of columns per process
  tt = real(ncol,kind=8)/real(nproc,kind=8)
  ncolpx = floor(tt)
  ii = ncol - nproc*ncolpx
  if (iproc<ii) then
      ncolp = ncolpx + 1
  else
      ncolp = ncolpx
  end if
  
  ! Determine the first column of each process
  i = 0
  do jproc=0,nproc-1
      if (iproc==jproc) iscol = i
      if (jproc<ii) then
          i = i + ncolpx + 1
      else
          i = i + ncolpx
      end if
  end do
end subroutine distribute_columns_on_processes
