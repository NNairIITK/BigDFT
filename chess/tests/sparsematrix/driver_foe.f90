!> @file
!!   Test of FOE using the CCS or the BigDFT format
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


program driver_foe
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
                                    trace_A, trace_AB, sparse_matrix_metadata_init_from_file, &
                                    sparse_matrix_and_matrices_init_from_file_bigdft
  use sparsematrix, only: write_matrix_compressed, transform_sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed, write_sparsematrix_info, &
                               get_number_of_electrons
  ! The following module is an auxiliary module for this test
  use utilities, only: get_ccs_data_from_file
  use futile
  use wrapper_MPI
  use wrapper_linalg

  implicit none

  ! Variables
  type(sparse_matrix) :: smat_s, smat_h, smat_k
  type(matrices) :: mat_s, mat_h, mat_k, mat_ek
  type(matrices),dimension(1) :: mat_ovrlpminusonehalf
  type(sparse_matrix_metadata) :: smmd
  integer :: nfvctr, nvctr, ierr, iproc, nproc, nthread, ncharge, nfvctr_mult, nvctr_mult
  integer,dimension(:),pointer :: row_ind, col_ptr, row_ind_mult, col_ptr_mult
  real(mp),dimension(:),pointer :: kernel, overlap, overlap_large
  real(mp),dimension(:),allocatable :: charge
  real(mp) :: energy, tr_KS, tr_KS_check
  type(foe_data) :: foe_obj, ice_obj
  real(mp) :: tr
  type(dictionary),pointer :: dict_timing_info, options
  type(yaml_cl_parse) :: parser !< command line parser
  character(len=1024) :: metadata_file, overlap_file, hamiltonian_file, kernel_file, kernel_matmul_file
  character(len=1024) :: sparsity_format, matrix_format
  external :: gather_timings
  !$ integer :: omp_get_max_threads

  ! Initialize flib
  call f_lib_initialize()

  ! MPI initialization; we have:
  ! iproc is the task ID
  ! nproc is the total number of tasks
  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()

  call f_malloc_set_status(memory_limit=0.e0,iproc=iproc)

  ! Initialize the sparsematrix error handling and timing.
  call sparsematrix_init_errors()
  call sparsematrix_initialize_timing_categories()


  if (iproc==0) then
      call yaml_new_document()
      !call print_logo()
  end if

  !Time initialization
  call f_timing_reset(filename='time.yaml',master=(iproc==0),verbose_mode=.false.)

  if (iproc==0) then
      call yaml_scalar('',hfill='~')
      call yaml_scalar('CHESS FOE TEST DRIVER',hfill='~')
  end if

  if (iproc==0) then
      call yaml_map('Timestamp of the run',yaml_date_and_time_toa())
      call yaml_mapping_open('Parallel environment')
      call yaml_map('MPI tasks',nproc)
      nthread = 1
      !$ nthread = omp_get_max_threads()
      call yaml_map('OpenMP threads',nthread)
      call yaml_mapping_close()
  end if

  ! Read in the parameters for the run and print them.
  if (iproc==0) then
      parser=yaml_cl_parse_null()
      call commandline_options(parser)
      call yaml_cl_parse_cmd_line(parser,args=options)
      call yaml_cl_parse_free(parser)

      metadata_file = options//'metadata_file'
      overlap_file = options//'overlap_file'
      hamiltonian_file = options//'hamiltonian_file'
      kernel_file = options//'kernel_file'
      kernel_matmul_file = options//'kernel_matmul_file'
      sparsity_format = options//'sparsity_format'
      matrix_format = options//'matrix_format'
     
      call dict_free(options)

      call yaml_mapping_open('Input parameters')
      call yaml_map('Sparsity format',trim(sparsity_format))
      call yaml_map('Matrix format',trim(matrix_format))
      call yaml_map('Metadata file',trim(metadata_file))
      call yaml_map('Overlap matrix file',trim(overlap_file))
      call yaml_map('Hamiltonian matrix file',trim(hamiltonian_file))
      call yaml_map('Density kernel matrix file',trim(kernel_file))
      call yaml_map('Density kernel matrix multiplication file',trim(kernel_matmul_file))
      call yaml_mapping_close()
  end if

  ! Send the input parameters to all MPI tasks
  call mpibcast(sparsity_format, root=0, comm=mpi_comm_world)
  call mpibcast(matrix_format, root=0, comm=mpi_comm_world)
  call mpibcast(metadata_file, root=0, comm=mpi_comm_world)
  call mpibcast(overlap_file, root=0, comm=mpi_comm_world)
  call mpibcast(hamiltonian_file, root=0, comm=mpi_comm_world)
  call mpibcast(kernel_file, root=0, comm=mpi_comm_world)
  call mpibcast(kernel_matmul_file, root=0, comm=mpi_comm_world)


  ! Read in the overlap matrix and create the type containing the sparse matrix descriptors (smat_s) as well as
  ! the type which contains the matrix data (overlap). The matrix element are stored in mat_s%matrix_compr.
  ! Do the same also for the Hamiltonian matrix.
  if (trim(sparsity_format)=='ccs') then
      if (trim(matrix_format)/='serial_text') then
          call f_err_throw('Wrong matrix format')
      end if
      !if (iproc==0) then
      !    call yaml_scalar('Initializing overlap matrix',hfill='-')
      !    call yaml_map('Reading from file','overlap_ccs.txt')
      !end if
      call sparse_matrix_and_matrices_init_from_file_ccs(trim(overlap_file), &
           iproc, nproc, mpi_comm_world, smat_s, mat_s)

      !if (iproc==0) then
      !    call yaml_scalar('Initializing Hamiltonian matrix',hfill='-')
      !    call yaml_map('Reading from file','hamiltonian_ccs.txt')
      !end if
      call sparse_matrix_and_matrices_init_from_file_ccs(trim(hamiltonian_file), &
           iproc, nproc, mpi_comm_world, smat_h, mat_h)

      ! Create another matrix type, this time directly with the CCS format descriptors.
      ! Get these descriptors from an auxiliary routine using matrix3.dat
      !if (iproc==0) then
      !    call yaml_scalar('Initializing density kernel matrix',hfill='-')
      !    call yaml_map('Reading from file','density_kernel_ccs.txt')
      !end if
      call get_ccs_data_from_file(trim(kernel_file), nfvctr, nvctr, row_ind, col_ptr)
      !call get_ccs_data_from_file('density_kernel_matmul_ccs.dat', nfvctr_mult, nvctr_mult, row_ind_mult, col_ptr_mult)
      call get_ccs_data_from_file(trim(kernel_matmul_file), nfvctr_mult, nvctr_mult, row_ind_mult, col_ptr_mult)
      if (nfvctr_mult/=nfvctr) then
          call f_err_throw('nfvctr_mult/=nfvctr',err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
      end if
      call sparse_matrix_init_from_data_ccs(iproc, nproc, mpi_comm_world, &
           nfvctr, nvctr, row_ind, col_ptr, smat_k, &
           init_matmul=.true., nvctr_mult=nvctr_mult, row_ind_mult=row_ind_mult, col_ptr_mult=col_ptr_mult)
      call f_free_ptr(row_ind)
      call f_free_ptr(col_ptr)
      call f_free_ptr(row_ind_mult)
      call f_free_ptr(col_ptr_mult)
  else if (trim(sparsity_format)=='bigdft') then
      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, overlap_file, &
           iproc, nproc, mpi_comm_world, smat_s, mat_s, init_matmul=.false.)
      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, hamiltonian_file, &
           iproc, nproc, mpi_comm_world, smat_h, mat_h, init_matmul=.false.)
      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, kernel_file, &
           iproc, nproc, mpi_comm_world, smat_k, mat_k, init_matmul=.true., filename_mult=trim(kernel_matmul_file))
  else
      call f_err_throw('Wrong sparsity format')
  end if

  if (iproc==0) then
      call yaml_mapping_open('Matrix properties')
      call write_sparsematrix_info(smat_s, 'Overlap matrix')
      call write_sparsematrix_info(smat_h, 'Hamiltonian matrix')
      call write_sparsematrix_info(smat_k, 'Density kernel')
      call yaml_mapping_close()
  end if

  call sparse_matrix_metadata_init_from_file('sparsematrix_metadata.dat', smmd)
  call get_number_of_electrons(smmd, ncharge)
  if (iproc==0) then
      call yaml_map('Number of electrons',ncharge)
  end if

  ! Prepares the type containing the matrix data.
  if (trim(sparsity_format)=='ccs') then
      ! Otherwise already done above
      call matrices_init(smat_k, mat_k)
  end if
  call matrices_init(smat_k, mat_ek)
  call matrices_init(smat_k, mat_ovrlpminusonehalf(1))

  ! Initialize the opaque object holding the parameters required for the Fermi Operator Expansion.
  ! Only provide the mandatory values and take for the optional values the default ones.
  charge = f_malloc(smat_s%nspin,id='charge')
  charge(:) = real(ncharge,kind=mp)
  call init_foe(iproc, nproc, smat_s%nspin, charge, foe_obj)
  ! Initialize the same object for the calculation of the inverse. Charge does not really make sense here...
  call init_foe(iproc, nproc, smat_s%nspin, charge, ice_obj, evlow=0.5_mp, evhigh=1.5_mp)

  call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(), &
       gather_routine=gather_timings)

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

  call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(), &
       gather_routine=gather_timings)

  !tr = trace_A(iproc, nproc, mpi_comm_world, smat_k, mat_ek, 1)
  tr = trace_AB(iproc, nproc, mpi_comm_world, smat_s, smat_k, mat_s, mat_ek, 1)
  if (iproc==0) then
      call yaml_map('Energy from FOE',energy)
      call yaml_map('Trace of energy density kernel', tr)
      call yaml_map('Difference',abs(energy-tr))
  end if

  !! Write the result in YAML format to the standard output (required for non-regression tests).
  !if (iproc==0) call write_matrix_compressed('Result of FOE', smat_k, mat_k)

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
  call deallocate_sparse_matrix_metadata(smmd)

  ! Deallocate all the matrix data types
  call deallocate_matrices(mat_s)
  call deallocate_matrices(mat_h)
  call deallocate_matrices(mat_k)
  call deallocate_matrices(mat_ek)
  call deallocate_matrices(mat_ovrlpminusonehalf(1))

  ! Deallocate all the remaining arrays
  call f_free(charge)
  call f_free_ptr(kernel)
  call f_free_ptr(overlap)
  call f_free_ptr(overlap_large)

  call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(), &
       gather_routine=gather_timings)

  call build_dict_info(iproc, nproc, dict_timing_info)
  call f_timing_stop(mpi_comm=mpi_comm_world,nproc=nproc,&
       gather_routine=gather_timings,dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  if (iproc==0) then
      call yaml_release_document()
  end if

  ! Finalize MPI
  call mpifinalize()

  ! Finalize flib
  call f_lib_finalize()


  !!contains

  !!  !> construct the dictionary needed for the timing information
  !!  !! SM: This routine should go to a module
  !!  subroutine build_dict_info(dict_info)
  !!    use wrapper_MPI
  !!    use dynamic_memory
  !!    use dictionaries
  !!    implicit none

  !!    type(dictionary), pointer :: dict_info
  !!    !local variables
  !!    integer :: ierr,namelen,nthreads
  !!    character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  !!    character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
  !!    type(dictionary), pointer :: dict_tmp
  !!    !$ integer :: omp_get_max_threads

  !!    call dict_init(dict_info)
! !! bastian: comment out 4 followinf lines for debug purposes (7.12.2014)
  !!    !if (DoLastRunThings) then
  !!       call f_malloc_dump_status(dict_summary=dict_tmp)
  !!       call set(dict_info//'Routines timing and number of calls',dict_tmp)
  !!    !end if
  !!    nthreads = 0
  !!    !$  nthreads=omp_get_max_threads()
  !!    call set(dict_info//'CPU parallelism'//'MPI tasks',nproc)
  !!    if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
  !!         nthreads)

  !!    nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0.to.nproc-1,id='nodename')
  !!    if (nproc>1) then
  !!       call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
  !!       !gather the result between all the process
  !!       call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
  !!            nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
  !!            mpi_comm_world,ierr)
  !!       if (iproc==0) call set(dict_info//'Hostnames',&
  !!               list_new(.item. nodename))
  !!    end if
  !!    call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)

  !!  end subroutine build_dict_info

end program driver_foe


subroutine commandline_options(parser)
  use yaml_parse
  use dictionaries, only: dict_new,operator(.is.)
  implicit none
  type(yaml_cl_parse),intent(inout) :: parser

  call yaml_cl_parse_option(parser,'metadata_file','sparsematrix_metadata.dat',&
       'input file with the matrix metadata',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the matrix metadata',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'overlap_file','overlap_sparse.txt',&
       'input file with the overlap matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the overlap matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'hamiltonian_file','hamiltonian_sparse.txt',&
       'input file with the Hamiltonian matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the Hamiltonian matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'kernel_file','kernel_sparse.txt',&
       'input file with the density kernel matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the density kernel matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'kernel_matmul_file','kernel_matmul_sparse.txt',&
       'input file with the density kernel matrix multiplication matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the density kernel matrix multiplication matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'sparsity_format','bigdft',&
       'indicate the sparsity format',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the sparsity format of the input matrices',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'matrix_format','serial_text',&
       'indicate the matrix format',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the format of the input matrices',&
       'Allowed values' .is. &
       'String'))

end subroutine commandline_options
