program driver_css
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               deallocate_sparse_matrix, deallocate_matrices, &
                               assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
  use foe_base, only: foe_data, foe_data_deallocate
  use foe_common, only: init_foe
  use bigdft_run, only: bigdft_init
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_ccs, &
                                    sparse_matrix_init_from_file_ccs, matrices_init, &
                                    matrices_get_values, matrices_set_values, &
                                    sparse_matrix_init_from_data_ccs, &
                                    ccs_data_from_sparse_matrix, ccs_matrix_write, &
                                    matrix_matrix_multiplication, matrix_fermi_operator_expansion, &
                                    trace_AB
  use utilities, only: get_ccs_data_from_file
  use sparsematrix, only: write_matrix_compressed
  use yaml_output
  implicit none

  ! Variables
  integer :: i
  type(sparse_matrix) :: smat_s, smat_h, smat_k
  type(matrices) :: mat_s, mat_h, mat_k
  type(matrices),dimension(1) :: mat_ovrlpminusonehalf
  integer :: nfvctr, nvctr, ierr
  integer,dimension(:),pointer :: row_ind, col_ptr
  real(kind=8),dimension(:),pointer :: val
  real(kind=8),dimension(:),allocatable :: charge
  real(kind=8) :: energy, tr_KS
  type(foe_data) :: foe_obj

  ! Initialize flib
  call f_lib_initialize()

  ! General initialization, including MPI. Here we have:
  ! bigdft_mpi%iproc is the task ID
  ! bigdft_mpi%nproc is the total number of tasks
  ! PROBLEM: This is in src/modules...
  call bigdft_init()

  ! Read from matrix1.dat and create the type containing the sparse matrix descriptors (smat_s) as well as
  ! the type which contains the matrix data (overlap). The matrix element are stored in mat_s%matrix_compr.
  ! Do the same also for matrix2.dat
  call sparse_matrix_and_matrices_init_from_file_ccs('overlap_ccs.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s, mat_s)
  call sparse_matrix_and_matrices_init_from_file_ccs('hamiltonian_ccs.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, smat_h, mat_h)

  ! Create another matrix type, this time directly with the CCS format descriptors.
  ! Get these descriptors from an auxiliary routine using matrix3.dat
  call get_ccs_data_from_file('density_kernel_ccs.dat', nfvctr, nvctr, row_ind, col_ptr)
  call sparse_matrix_init_from_data_ccs(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr, nvctr, row_ind, col_ptr, smat_k)

  ! Prepares the type containing the matrix data.
  call matrices_init(smat_k, mat_k)
  call matrices_init(smat_k, mat_ovrlpminusonehalf(1))

  ! Initialize the opaque object holding the parameters required for the Fermi Operator Expansion.
  ! Only provide the mandatory values and take for the optional values the default ones.
  charge = f_malloc(smat_s%nspin,id='charge')
  charge(:) = 722.d0
  call init_foe(smat_s%nspin, charge, foe_obj)

  ! Calculate the density kernel for the system described by the pair smat_s/mat_s and smat_h/mat_h and 
  ! store the result in smat_k/mat_k. Attention: The sparsity pattern of smat_s must be contained within that of smat_h
  ! and the one of smat_h within that of smat_k. It is your responsabilty to assure this, 
  ! the routine does only some minimal checks.
  ! The final result is contained in mat_k(1)%matrix_compr.
  call matrix_fermi_operator_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, foe_obj, smat_s, smat_h, smat_k, &
       mat_s, mat_h, mat_ovrlpminusonehalf, mat_k, energy, &
       calculate_minusonehalf=.true., foe_verbosity=0)

  ! Write the result in YAML format to the standard output (required for non-regression tests).
  if (bigdft_mpi%iproc==0) call write_matrix_compressed('Result of FOE', smat_k, mat_k)

  ! Calculate trace(KS)
  !tr_KS = trace_sparse(bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s, smat_k, mat_s%matrix_compr, mat_k%matrix_compr, 1)
  tr_KS = trace_AB(smat_s, smat_k, mat_s, mat_k, 1)

  ! Write the result
  if (bigdft_mpi%iproc==0) call yaml_map('trace(KS)',tr_KS)



  ! Extract the compressed matrix from the data type. The first routine allocates an array with the correct size,
  ! the second one extracts the result.
  val = sparsematrix_malloc_ptr(smat_k, iaction=SPARSE_FULL, id='val')
  call matrices_get_values(smat_k, mat_k, val)

!!  ! Create another matrix type, this time directly with the CCS format descriptors.
!!  ! Get these descriptors from an auxiliary routine using again matrix2.dat
!!  call get_ccs_data_from_file('matrix2.dat', nfvctr, nvctr, row_ind, col_ptr)
!!  call sparse_matrix_init_from_data_ccs(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr, nvctr, row_ind, col_ptr, smat_k)

!!  ! Extract the compressed matrix from the data type. The first routine allocates an array with the correct size,
!!  ! the second one extracts the result.
!!  val = sparsematrix_malloc_ptr(smat_h, iaction=SPARSE_FULL, id='val')
!!  call matrices_get_values(smat_h, mat_h, val)

!!  ! Prepare a new matrix data type.
!!  call matrices_init(smat_h, mat_h(2))
!!  ! Set the contents of this new matrix type, using the above result.
!!  call matrices_set_values(smat_h, val, mat_h(2))

  ! Deallocate the array
  call f_free_ptr(val)

!!  ! Prepare two new matrix data types
!!  do i=1,2
!!      call matrices_init(smat_k, mat_k(i))
!!  end do
!!
!!  ! Calculate at the same time the square root and the inverse square of the matrix described  by the pair smat_h/mat_h
!!  ! and store the result in smat_k/mat_k. The final results are thus in mat_k(1)%matrix_compr and mat_k(2)%matrix_compr.
!!  ! The same wraning as above applies.
!!  call matrix_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, 2, (/0.5d0,-0.5d0/), &
!!       smat_h, smat_k, mat_h(2), mat_k)
!!
!!  ! Write the result in YAML format to the standard output (required for non-regression tests).
!!  if (bigdft_mpi%iproc==0) call write_matrix_compressed('Result of second matrix', smat_k, mat_k(1))
!!  if (bigdft_mpi%iproc==0) call write_matrix_compressed('Result of third matrix', smat_k, mat_k(2))
!!
!!  ! Calculate the CCS descriptors from the sparse_matrix type.
!!  call ccs_data_from_sparse_matrix(smat_k, row_ind, col_ptr)
!!
!!  ! Write the two matrices to disk, using the CCS format
!!  if (bigdft_mpi%iproc==0) call ccs_matrix_write('squareroot.dat', smat_k, row_ind, col_ptr, mat_k(1))
!!  if (bigdft_mpi%iproc==0) call ccs_matrix_write('invsquareroot.dat', smat_k, row_ind, col_ptr, mat_k(2))
!!
!!  ! Multiply the two matrices calculated above (i.e. the square root and the inverse square root) and
!!  ! store the result in mat_h. The final result is thus contained in mat_h%matrix_compr.
!!  call matrix_matrix_multiplication(bigdft_mpi%iproc, bigdft_mpi%nproc, smat_k, &
!!       mat_k(1), mat_k(2), mat_h(1))
!!
!!  ! Write the result of the above multiplication to a file. Since we multiply the square root times the 
!!  ! inverse square root, the result will be the unity matrix.
!!  call ccs_matrix_write('unity.dat', smat_k, row_ind, col_ptr, mat_h(1))
!!
!!  ! Write the result also in YAML format to the standard output (required for non-regression tests).
!!  if (bigdft_mpi%iproc==0) call write_matrix_compressed('Result of fourth matrix', smat_k, mat_h(1))
!!

  call foe_data_deallocate(foe_obj)

  ! Deallocate all the sparse matrix descriptrs types
  call deallocate_sparse_matrix(smat_s)
  call deallocate_sparse_matrix(smat_h)
  call deallocate_sparse_matrix(smat_k)
!!
!!  ! Deallocate all the matrix data types
  call deallocate_matrices(mat_s)
  call deallocate_matrices(mat_h)
  call deallocate_matrices(mat_k)
  call deallocate_matrices(mat_ovrlpminusonehalf(1))
  call f_free(charge)
!!  do i=1,2
!!      call deallocate_matrices(mat_h(i))
!!      call deallocate_matrices(mat_k(i))
!!  end do
!!
!!  ! Deallocate the CCS format descriptors
  call f_free_ptr(row_ind)
  call f_free_ptr(col_ptr)

  ! Finalize MPI
  call bigdft_finalize(ierr)

  ! Finalize flib
  call f_lib_finalize()


end program driver_css
