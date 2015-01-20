program driver_singlerun
  use bigdft_run
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
  use sparsematrix_init, only: read_ccs_format, ccs_to_bigdft, ccs_values_to_bigdft
  use sparsematrix, only: write_matrix_compressed, check_symmetry
  use module_interfaces, only: overlapPowerGeneral
  implicit none

  ! Variables
  integer :: iproc, nproc, ncol, nnonzero
  character(len=*),parameter :: filename='matrix.dat'
  integer,dimension(:),pointer :: col_ptr, row_ind
  real(kind=8),dimension(:),pointer :: val
  type(sparse_matrix) :: smat
  type(matrices) :: matA
  type(matrices),dimension(1) :: matB
  real(kind=8) :: max_error, mean_error
  logical :: symmetric

  ! Initialize
  call f_lib_initialize()
  call bigdft_init()!mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  iproc=bigdft_mpi%iproc!mpi_info(1)
  nproc=bigdft_mpi%nproc!mpi_info(2)

  
  ! Read in a file in the CCS format
  call read_ccs_format(filename, ncol, nnonzero, col_ptr, row_ind, val)

  ! Create the corresponding BigDFT sparsity pattern
  call ccs_to_bigdft(iproc, nproc, ncol, ncol, 0, nnonzero, row_ind, col_ptr, smat)

  ! Assign the values
  matA%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matA%matrix_compr')
  call ccs_values_to_bigdft(ncol, nnonzero, row_ind, col_ptr, smat, val, matA)

  ! Check the symmetry
  symmetric = check_symmetry(ncol, smat)
  if (.not.symmetric) stop 'ERROR not symmetric'
  
  ! Write the original matrix
  call write_matrix_compressed('Original matrix', smat, matA)

  ! Calculate the inverse
  matB(1)%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matB(1)%matrix_compr')
  !!call overlapPowerGeneral(iproc, nproc, 1020, 1, (/1/), -8, &
  !!     1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
  !!     check_accur=.true., max_error=max_error, mean_error=mean_error)
  call overlapPowerGeneral(iproc, nproc, 0, 1, (/1/), -8, &
       1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
       check_accur=.true., max_error=max_error, mean_error=mean_error)

  ! Write the results
  call write_matrix_compressed('Result', smat, matB(1))

end program driver_singlerun
