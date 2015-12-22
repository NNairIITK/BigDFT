program driver
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               deallocate_sparse_matrix, deallocate_matrices
  use bigdft_run, only: bigdft_init
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_ccs, &
                                    sparse_matrix_init_from_file_ccs, matrices_init, &
                                    matrices_set, sparse_matrix_init_from_data
  use ice, only: inverse_chebyshev_expansion
  use utilities, only: get_ccs_data_from_file
  implicit none

  ! Variables
  integer :: i
  integer :: norder_polynomial
  type(sparse_matrix) :: smat_s, smat_l, smat3
  type(matrices) :: overlap
  type(matrices),dimension(1) :: inv_overlap, mat2
  type(matrices),dimension(2) :: mat3
  integer :: nfvctr, nvctr
  integer,dimension(:),pointer :: row_ind, col_ptr


  call f_lib_initialize()

  call bigdft_init()

  ! Creates the type which has the sparse matrix descriptors (smat_s) as well as
  ! the type which contains the matrix data (overlap). The matrix element are stored in overlap%matrix_compr.
  call sparse_matrix_and_matrices_init_from_file_ccs('overlap_ccs.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s, overlap)

  ! Creates the type which has the sparse matrix descriptors (smat_l).
  call sparse_matrix_init_from_file_ccs('kernel_ccs.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, smat_l)

  ! This routine prepares the type containing the matrix data.
  call matrices_init(smat_l, inv_overlap(1))

  ! Calculate the square root of the matrix
  norder_polynomial = 30
  call inverse_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, norder_polynomial, &
       smat_s, smat_l, 1, (/0.5d0/), overlap, inv_overlap)

  ! Create another matrix type, this time directly with the CCS format descriptors. Get these descriptors
  ! from an auxiliary routine.
  call get_ccs_data_from_file('kernel_ccs.dat', nfvctr, nvctr, row_ind, col_ptr)
  call sparse_matrix_init_from_data(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr, nvctr, row_ind, col_ptr, smat3)
  call f_free_ptr(row_ind)
  call f_free_ptr(col_ptr)

  ! Prepare a new matrix data type.
  call matrices_init(smat_l, mat2(1))
  ! Set the contents of this new matrix type
  call matrices_set(smat_l, inv_overlap(1)%matrix_compr, mat2(1))

  ! Prepare two new matrices
  do i=1,2
      call matrices_init(smat3, mat3(i))
  end do

  ! Calculate the square root and the inverse square of the matrix
  norder_polynomial = 40
  call inverse_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, norder_polynomial, &
       smat_l, smat3, 2, (/0.5d0,-0.5d0/), mat2(1), mat3)

  call deallocate_sparse_matrix(smat_s)
  call deallocate_sparse_matrix(smat_l)
  call deallocate_sparse_matrix(smat3)
  call deallocate_matrices(overlap)
  call deallocate_matrices(inv_overlap(1))
  call deallocate_matrices(mat2(1))
  do i=1,2
      call deallocate_matrices(mat3(i))
  end do

  call f_lib_finalize()


end program driver
