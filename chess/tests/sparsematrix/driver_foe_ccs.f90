program driver_css
  ! The following module are part of the sparsematrix library
  use sparsematrix_base
  use foe_base, only: foe_data, foe_data_deallocate
  use foe_common, only: init_foe
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_ccs, &
                                    sparse_matrix_init_from_file_ccs, matrices_init, &
                                    matrices_get_values, matrices_set_values, &
                                    sparse_matrix_init_from_data_ccs, &
                                    ccs_data_from_sparse_matrix, ccs_matrix_write, &
                                    matrix_matrix_multiplication, matrix_fermi_operator_expansion, &
                                    trace_A, trace_AB
  use sparsematrix, only: write_matrix_compressed, transform_sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  ! The following module is an auxiliary module for this test
  use utilities, only: get_ccs_data_from_file
  use futile
  implicit none

  ! Variables
  type(sparse_matrix) :: smat_s, smat_h, smat_k
  type(matrices) :: mat_s, mat_h, mat_k, mat_ek
  type(matrices),dimension(1) :: mat_ovrlpminusonehalf
  integer :: nfvctr, nvctr, ierr, iproc, nproc
  integer,dimension(:),pointer :: row_ind, col_ptr
  real(mp),dimension(:),pointer :: kernel, overlap, overlap_large
  real(mp),dimension(:),allocatable :: charge
  real(mp) :: energy, tr_KS, tr_KS_check
  type(foe_data) :: foe_obj, ice_obj
  real(mp) :: tr

  ! Initialize flib
  call f_lib_initialize()

  ! MPI initialization; we have:
  ! iproc is the task ID
  ! nproc is the total number of tasks
  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()

  ! Initialize the sparsematrix error handling and timing.
  call sparsematrix_init_errors()
  call sparsematrix_initialize_timing_categories()

  ! Read from matrix1.dat and create the type containing the sparse matrix descriptors (smat_s) as well as
  ! the type which contains the matrix data (overlap). The matrix element are stored in mat_s%matrix_compr.
  ! Do the same also for matrix2.dat
  call sparse_matrix_and_matrices_init_from_file_ccs('overlap_ccs.dat', &
       iproc, nproc, mpi_comm_world, smat_s, mat_s)
  call sparse_matrix_and_matrices_init_from_file_ccs('hamiltonian_ccs.dat', &
       iproc, nproc, mpi_comm_world, smat_h, mat_h)

  ! Create another matrix type, this time directly with the CCS format descriptors.
  ! Get these descriptors from an auxiliary routine using matrix3.dat
  call get_ccs_data_from_file('density_kernel_ccs.dat', nfvctr, nvctr, row_ind, col_ptr)
  call sparse_matrix_init_from_data_ccs(iproc, nproc, mpi_comm_world, &
       nfvctr, nvctr, row_ind, col_ptr, smat_k, &
       init_matmul=.true., nvctr_mult=nvctr, row_ind_mult=row_ind, col_ptr_mult=col_ptr)

  ! Prepares the type containing the matrix data.
  call matrices_init(smat_k, mat_k)
  call matrices_init(smat_k, mat_ek)
  call matrices_init(smat_k, mat_ovrlpminusonehalf(1))

  ! Initialize the opaque object holding the parameters required for the Fermi Operator Expansion.
  ! Only provide the mandatory values and take for the optional values the default ones.
  charge = f_malloc(smat_s%nspin,id='charge')
  charge(:) = 722.d0
  call init_foe(iproc, nproc, smat_s%nspin, charge, foe_obj)
  ! Initialize the same object for the calculation of the inverse. Charge does not really make sense here...
  call init_foe(iproc, nproc, smat_s%nspin, charge, ice_obj, evlow=0.5_mp, evhigh=1.5_mp)

  ! Calculate the density kernel for the system described by the pair smat_s/mat_s and smat_h/mat_h and 
  ! store the result in smat_k/mat_k.
  ! Attention: The sparsity pattern of smat_s must be contained within that of smat_h
  ! and the one of smat_h within that of smat_k. It is your responsabilty to assure this, 
  ! the routine does only some minimal checks.
  ! The final result will be contained in mat_k%matrix_compr.
  call matrix_fermi_operator_expansion(iproc, nproc, mpi_comm_world, &
       foe_obj, ice_obj, smat_s, smat_h, smat_k, &
       mat_s, mat_h, mat_ovrlpminusonehalf, mat_k, energy, &
       calculate_minusonehalf=.true., foe_verbosity=1, symmetrize_kernel=.true., &
       calculate_energy_density_kernel=.true., energy_kernel=mat_ek)

  if (iproc==0) then
      call yaml_map('Energy from FOE',energy)
      tr = trace_A(iproc, nproc, mpi_comm_world, smat_k, mat_ek, 1)
      call yaml_map('Trace of energy density kernel', tr)
      call yaml_map('Difference',abs(energy-tr))
  end if

  ! Write the result in YAML format to the standard output (required for non-regression tests).
  if (iproc==0) call write_matrix_compressed('Result of FOE', smat_k, mat_k)

  ! Calculate trace(KS)
  !tr_KS = trace_sparse(iproc, nproc, smat_s, smat_k, mat_s%matrix_compr, mat_k%matrix_compr, 1)
  tr_KS = trace_AB(iproc, nproc, mpi_comm_world, smat_s, smat_k, mat_s, mat_k, 1)

  ! Write the result
  if (iproc==0) call yaml_map('trace(KS)',tr_KS)

  ! Extract the compressed kernel matrix from the data type.
  ! The first routine allocates an array with the correct size, the second one extracts the result.
  kernel = sparsematrix_malloc_ptr(smat_k, iaction=SPARSE_FULL, id='kernel')
  call matrices_get_values(smat_k, mat_k, kernel)

  ! Do the same also for the overlap matrix
  overlap = sparsematrix_malloc_ptr(smat_s, iaction=SPARSE_FULL, id='overlap')
  call matrices_get_values(smat_s, mat_s, overlap)

  ! Transform the overlap matrix to the sparsity pattern of the kernel
  overlap_large = sparsematrix_malloc_ptr(smat_k, iaction=SPARSE_FULL, id='overlap_large')
  call transform_sparse_matrix(iproc, smat_s, smat_k, SPARSE_FULL, 'small_to_large', &
       smat_in=overlap, lmat_out=overlap_large)

  ! Again calculate trace(KS), this time directly with the array holding the data.
  ! Since both matrices are symmetric and have now the same sparsity pattern, this is a simple ddot.
  tr_KS_check = dot(smat_k%nvctr, kernel(1), 1, overlap_large(1), 1)

  ! Write the result
  if (iproc==0) call yaml_map('trace(KS) check',tr_KS_check)

  ! Write the difference to the previous result
  if (iproc==0) call yaml_map('difference',tr_KS-tr_KS_check)

  ! Deallocate the object holding the FOE parameters
  call foe_data_deallocate(foe_obj)
  call foe_data_deallocate(ice_obj)

  ! Deallocate all the sparse matrix descriptors types
  call deallocate_sparse_matrix(smat_s)
  call deallocate_sparse_matrix(smat_h)
  call deallocate_sparse_matrix(smat_k)

  ! Deallocate all the matrix data types
  call deallocate_matrices(mat_s)
  call deallocate_matrices(mat_h)
  call deallocate_matrices(mat_k)
  call deallocate_matrices(mat_ovrlpminusonehalf(1))

  ! Deallocate all the remaining arrays
  call f_free(charge)
  call f_free_ptr(row_ind)
  call f_free_ptr(col_ptr)
  call f_free_ptr(kernel)
  call f_free_ptr(overlap)
  call f_free_ptr(overlap_large)

  ! Finalize MPI
  call mpifinalize()

  ! Finalize flib
  ! SM: I have the impression that every task should call this routine, but if I do so
  ! some things are printed nproc times instead of once.
  if (iproc==0) call f_lib_finalize()

end program driver_css
