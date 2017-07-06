!> @file
!!   Flexible test of the matrix power expansion
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


program driver_random
  use futile
  use sparsematrix_base
  use sparsematrix_types, only: sparse_matrix
  use sparsematrix_init, only: generate_random_symmetric_sparsity_pattern, &
                               matrixindex_in_compressed, write_sparsematrix_info, &
                               init_matrix_taskgroups_wrapper
  use sparsematrix, only: symmetrize_matrix, check_deviation_from_unity_sparse, &
                          matrix_power_dense_lapack, get_minmax_eigenvalues, &
                          uncompress_matrix, uncompress_matrix2, transform_sparse_matrix, &
                          resize_matrix_to_taskgroup
  use sparsematrix_highlevel, only: matrix_chebyshev_expansion, matrices_init, &
                                    matrix_matrix_multiplication, &
                                    sparse_matrix_init_from_file_bigdft, &
                                    sparse_matrix_and_matrices_init_from_file_bigdft, &
                                    sparse_matrix_metadata_init_from_file
  use sparsematrix_io, only: write_dense_matrix, write_sparse_matrix
  !use random, only: builtin_rand
  use f_random, only: f_random_number
  use foe_base, only: foe_data, foe_data_deallocate
  use foe_common, only: init_foe
  use wrapper_MPI
  use selinv, only: selinv_wrapper
  use utilities, only: calculate_error, median

  implicit none

  ! Variables
  integer :: iproc, nproc, iseg, ierr, idum, ii, i, nthread, nit, it, icons
  integer :: nfvctr, nvctr, nbuf_large, nbuf_mult, iwrite, blocksize_diag, blocksize_matmul
  integer :: ithreshold, icheck, j, pexsi_np_sym_fact, output_level, profiling_depth
  type(sparse_matrix),dimension(2) :: smat
  type(sparse_matrix_metadata) :: smmd
  real(kind=4) :: tt_real
  real(mp) :: tt, tt_rel, t1, t2
  real(mp),dimension(1) :: eval_min, eval_max
  type(matrices) :: mat1, mat2
  type(matrices),dimension(3) :: mat3
  real(mp) :: condition_number, expo, max_error, mean_error, betax
  real(mp) :: max_error_rel, mean_error_rel, evlow, evhigh, eval_multiplicator, accuracy_ice, accuracy_penalty
  real(mp),dimension(:),allocatable :: charge_fake, times, sumarr
  type(foe_data),dimension(:),allocatable :: ice_obj
  character(len=1024) :: infile, outfile, outmatmulfile, sparsegen_method, matgen_method, diag_algorithm
  character(len=1024) :: metadata_file, solution_method
  logical :: write_matrices, do_cubic_check, init_matmul, do_consistency_checks
  type(dictionary), pointer :: options, dict_timing_info
  type(yaml_cl_parse) :: parser !< command line parser
  external :: gather_timings
  !$ integer :: omp_get_max_threads
  integer,parameter :: nthreshold = 10 !< number of checks with threshold
  real(mp),dimension(nthreshold),parameter :: threshold = (/ 1.e-3_mp, &
                                                             1.e-4_mp, &
                                                             1.e-5_mp, &
                                                             1.e-6_mp, &
                                                             1.e-7_mp, &
                                                             1.e-8_mp, &
                                                             1.e-9_mp, &
                                                             1.e-10_mp,&
                                                             1.e-11_mp,&
                                                             1.e-12_mp /) !< threshold for the relative errror
  !!integer,dimension(nthreshold) :: nrel_threshold
  !!real(mp),dimension(nthreshold) :: max_error_rel_threshold, mean_error_rel_threshold

  ! Initialize flib
  call f_lib_initialize()

  ! MPI initialization; we have:
  ! iproc is the task ID
  ! nproc is the total number of tasks
  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()

  ! Read in the parameters for the run.
  call read_and_communicate_input_variables()

  call f_malloc_set_status(iproc=iproc,output_level=output_level,&
       logfile_name='mem.log')

  ! Initialize the sparse matrix errors and timings.
  call sparsematrix_init_errors()
  call sparsematrix_initialize_timing_categories()

  if (iproc==0) then
      call yaml_new_document()
      !call print_logo()
  end if

  ! Timing initialization
  call mpibarrier()
  call f_timing_reset(filename='time.yaml', master=(iproc==0), verbose_mode=.false.)

  if (iproc==0) then
      call yaml_scalar('',hfill='~')
      call yaml_scalar('DRIVER FOR THE CHEBYSHEV MATRIX EXPANSION',hfill='~')
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



  if (iproc==0) then
      call yaml_mapping_open('Input parameters')
      call yaml_map('Metadata file',trim(metadata_file))
      call yaml_map('Sparsity pattern generation',trim(sparsegen_method))
      call yaml_map('Matrix content generation',trim(matgen_method))
      if (trim(sparsegen_method)=='random') then
          if (trim(matgen_method)=='file') then
              call f_err_throw("Inconsistency between 'sparsegen_method' and 'matgen_method'")
          end if
          call yaml_map('Matrix dimension',nfvctr)
          call yaml_map('Number of non-zero entries',nvctr)
          call yaml_map('Buffer for large matrix',nbuf_large)
          call yaml_map('Buffer for sparse multiplications',nbuf_mult)
      else if (trim(sparsegen_method)=='file') then
          !!if (trim(matgen_method)=='random') then
          !!    call f_err_throw("Inconsistency between 'sparsegen_method' and 'matgen_method'")
          !!end if
          call yaml_map('File with input sparsity pattern',trim(infile))
          call yaml_map('File with output sparsity pattern',trim(outfile))
          call yaml_map('File with output matrix multiplication sparsity pattern',trim(outmatmulfile))
      else
          call f_err_throw("Wrong value for 'sparsegen_method'")
      end if
      if (trim(matgen_method)=='random') then
          call yaml_map('Condition number',condition_number)
      else if (trim(matgen_method)=='file') then
      else
          call f_err_throw("Wrong value for 'matgen_method'")
      end if
      call yaml_map('Exponent for the matrix power calculation',expo)
      call yaml_map('Solution method',trim(solution_method))
      call yaml_map('Write the matrices',write_matrices)
      call yaml_map('betax',betax,fmt='(f9.1)')
      call yaml_map('Initial minimal eigenvalue',evlow)
      call yaml_map('Initial maximal eigenvalue',evhigh)
      call yaml_map('blocksize_diag',blocksize_diag)
      call yaml_map('blocksize_matmul',blocksize_matmul)
      call yaml_map('ScaLAPACK diagonalization algorithm',diag_algorithm)
      call yaml_map('ICE multiplication factor',eval_multiplicator)
      call yaml_map('PEXSI number of procs for symbolic factorization',pexsi_np_sym_fact)
      call yaml_map('Do a check with cubic scaling (Sca)LAPACK',do_cubic_check)
      call yaml_mapping_close()
      call yaml_map('Accuracy of Chebyshev fit for ICE',accuracy_ice)
      call yaml_map('Accuracy of Chebyshev fit for the penalty function',accuracy_penalty)
      call yaml_map('Number of iterations',nit)
      call yaml_map('Do consistency checks',do_consistency_checks)
      call yaml_map('Routine profiling output level',output_level)
      call yaml_map('Routine timing profiling depth',profiling_depth)
  end if


  if (trim(sparsegen_method)=='random') then
      ! Generate a random sparsity pattern
      call generate_random_symmetric_sparsity_pattern(iproc, nproc, mpi_comm_world, &
           nfvctr, nvctr, nbuf_mult, .false., smat(1), &
           1, (/nbuf_large/), (/.true./), smat(2:2))
  else
      if (trim(matgen_method)=='random') then
          call sparse_matrix_init_from_file_bigdft('serial_text', trim(infile), &
              iproc, nproc, mpi_comm_world, smat(1), &
              init_matmul=.false.)
      else if (trim(matgen_method)=='file') then
          call sparse_matrix_and_matrices_init_from_file_bigdft('serial_text', trim(infile), &
              iproc, nproc, mpi_comm_world, smat(1), mat2,&
              init_matmul=.false.)
      end if
      select case(trim(solution_method))
      case ('ICE', 'SelInv')
          init_matmul = .true.
      case('LAPACK')
          init_matmul = .false.
      case default
          call f_err_throw('wrong value for kernel_method')
      end select
      if (do_consistency_checks) then
          ! Since these checks are done with CheSS, init_matmul must be true
          init_matmul = .true.
      end if
      call sparse_matrix_init_from_file_bigdft('serial_text', trim(outfile), &
          iproc, nproc, mpi_comm_world, smat(2), &
          init_matmul=init_matmul, filename_mult=trim(outmatmulfile))
  end if


  !!call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)

  call init_matrix_taskgroups_wrapper(iproc, nproc, mpi_comm_world, .true., 2, smat)


  ! Write a summary of the sparse matrix layout 
  if (iproc==0) then
      call yaml_mapping_open('Matrix properties')
      call write_sparsematrix_info(smat(1), 'Random matrix')
      call write_sparsematrix_info(smat(2), 'Large random matrix')
      call yaml_mapping_close()
  end if


  times = f_malloc(nit,id='times')
  sumarr = f_malloc(nit,id='sumarr')
  allocate(ice_obj(nit))


  ! Initialize an object which holds some parameters for the Chebyshev expansion.
  ! It would also work without (then always the default values are taken), but when this
  ! object is provided some informations are stored between calls to the Chebyshev expansion,
  ! in this way improving the performance.
  ! Should maybe go to a wrapper.
  charge_fake = f_malloc0(1,id='charge_fake')
  do it=1,nit
      call init_foe(iproc, nproc, 1, charge_fake, ice_obj(it), evlow=evlow, evhigh=evhigh, &
           betax=betax, eval_multiplicator=eval_multiplicator, &
           accuracy_function=accuracy_ice, accuracy_penalty=accuracy_penalty)
  end do
  call f_free(charge_fake)


  ! Allocate the matrices
  mat3(:) = matrices_null()
  call matrices_init(smat(2), mat3(1), matsize=SPARSE_TASKGROUP)
  if (do_consistency_checks) then
      call matrices_init(smat(2), mat3(2), matsize=SPARSE_TASKGROUP)
  end if
  if (do_consistency_checks .or. do_cubic_check) then
      call matrices_init(smat(2), mat3(3), matsize=SPARSE_TASKGROUP)
  end if

  if (trim(matgen_method)=='random') then

      call matrices_init(smat(1), mat1)
      call matrices_init(smat(1), mat2, matsize=SPARSE_TASKGROUP)


    
      ! Fill the matrix with random entries
      idum = 0
      !tt_real = builtin_rand(idum, reset=.true.)
      call f_random_number(tt_real, reset=.true.)
      do i=1,smat(1)%nvctr
          !tt_real = builtin_rand(idum)
          call f_random_number(tt_real)
          mat1%matrix_compr(i) = real(tt_real,kind=8)
      end do

    
      ! By construction, all elements lie between 0 and 1. If we thus scale the
      ! matrix by the inverse of nfvctr (i.e. its dimension), then the sum of each line
      ! is between 0 and 1.
      call dscal(smat(1)%nvctr, 1.0_mp/real(smat(1)%nfvctr,kind=8), mat1%matrix_compr(1), 1)
    
      ! By construction, the sum of each line is between 0 and 1. If we thus set the diagonal
      ! to 1, we get a diagonally dominant matrix, which is positive definite.
      ! Additionally, we add to the diagonal elements a random number between 0 and the condition number.
      do i=1,smat(1)%nfvctr
          ii = matrixindex_in_compressed(smat(1), i, i)
          !tt_real = builtin_rand(idum)
          call f_random_number(tt_real)
          tt = real(tt_real,kind=8)*condition_number
          mat1%matrix_compr(ii) = 1.0_mp + tt
      end do

    
      ! Scale the matrix by the condition number, which should move down the largest
      ! eigenvalue of the order of 1.
      call dscal(smat(1)%nvctr, 1.0_mp/condition_number, mat1%matrix_compr(1), 1)


      ! Resize the matrix to the task group
      call resize_matrix_to_taskgroup(smat(1), mat1)

      ! Symmetrize the matrix
      call symmetrize_matrix(smat(1), 'plus', mat1%matrix_compr, mat2%matrix_compr)

      call deallocate_matrices(mat1)
  else
      ! Resize the matrix to the rask group
      call resize_matrix_to_taskgroup(smat(1), mat2)
  end if


  ! Initialization part done
  !call timing(mpi_comm_world,'INIT','PR')
  call mpibarrier()
  call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)

  ! Calculate the minimal and maximal eigenvalue, to determine the condition number
  call get_minmax_eigenvalues(iproc, nproc, mpiworld(), 'standard', blocksize_diag, &
       smat(1), mat2, eval_min, eval_max, &
       algorithm=diag_algorithm, quiet=.true.)
  if (iproc==0) then
      call yaml_mapping_open('Eigenvalue properties')
      call yaml_map('Minimal',eval_min)
      call yaml_map('Maximal',eval_max)
      call yaml_map('Condition number',eval_max/eval_min)
      call yaml_mapping_close()
  end if

  !call write_dense_matrix(iproc, nproc, mpi_comm_world, smat(1), mat2, 'randommatrix.dat', binary=.false.)
  if (write_matrices) then
      call write_sparse_matrix('serial_text', iproc, nproc, mpi_comm_world, smat(1), mat2, 'randommatrix_sparse')
  end if

  call mpibarrier()
  call f_timing_checkpoint(ctr_name='INFO',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)


  if (iproc==0) then
      call yaml_sequence_open('Matrix power calculations')
  end if
  it_loop: do it=1,nit
      t1 = mpi_wtime()
      if (iproc==0) then
          call yaml_comment('Kernel iteration number'//trim(yaml_toa(it)),hfill='=')
          call yaml_sequence(advance='no')
      end if
      ! Calculate the desired matrix power
      if (iproc==0) then
          call yaml_comment('Calculating mat^x',hfill='~')
          call yaml_mapping_open('Calculating mat^x')
      end if
      if (trim(solution_method)=='ICE') then
          call matrix_chebyshev_expansion(iproc, nproc, mpi_comm_world, &
               1, (/expo/), smat(1), smat(2), mat2, mat3(1), ice_obj=ice_obj(it))
      else if (trim(solution_method)=='SelInv') then
          if (expo/=-1.0_mp) then
              call f_err_throw('Selecetd Inversion is only possible for the calculation of the inverse')
          end if
          call selinv_wrapper(iproc, nproc, mpi_comm_world, smat(1), smat(2), mat2, pexsi_np_sym_fact, mat3(1))
      else if (trim(solutioN_method)=='LAPACK') then
          mat2%matrix = sparsematrix_malloc_ptr(smat(1), iaction=DENSE_FULL, id='mat2%matrix')
          mat3(1)%matrix = sparsematrix_malloc_ptr(smat(1), iaction=DENSE_FULL, id='mat3(3)%matrix')
          call matrix_power_dense_lapack(iproc, nproc, mpiworld(), blocksize_diag, blocksize_matmul, .false., &
                expo, smat(1), smat(2), mat2, mat3(1), algorithm=diag_algorithm, overwrite=.true.)
          call f_free_ptr(mat2%matrix)
          call f_free_ptr(mat3(1)%matrix)
      else
          call f_err_throw("wrong value for 'solution_method'; possible values are 'ICE', 'SelInv' or 'LAPACK'")
      end if
      if (iproc==0) then
          call yaml_mapping_close()
      end if
      ! Calculation part done

      call mpibarrier()
      t2 = mpi_wtime()
      times(it) = t2-t1
      sumarr(it) = sum(mat3(1)%matrix_compr)
  end do it_loop

  if (iproc==0) then
      call yaml_comment('Timings summary',hfill='=')
      call yaml_sequence_close()
      call yaml_sequence_open('Timinig summary')
      do it=1,nit
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
           call yaml_map('iteration',it)
           call yaml_map('time',times(it))
           call yaml_map('sum(M^a)',sumarr(it))
          call yaml_mapping_close()
      end do
      call yaml_sequence_close()
      call yaml_mapping_open('Timing analysis')
       call yaml_mapping_open('All runs')
        call yaml_map('Minimal',minval(times))
        call yaml_map('Maximal',maxval(times))
        call yaml_map('Average',sum(times)/real(nit,kind=mp))
        call yaml_map('Median',median(nit, times))
       call yamL_mapping_close()
       if (nit>1) then
           call yaml_mapping_open('Excluding first run')
            call yaml_map('Minimal',minval(times(2:nit)))
            call yaml_map('Maximal',maxval(times(2:nit)))
            call yaml_map('Average',sum(times(2:nit))/real(nit-1,kind=mp))
            call yaml_map('Median',median(nit-1, times(2:nit)))
           call yamL_mapping_close()
       end if
      call yaml_mapping_close()
      call yaml_scalar('',hfill='=')
  end if

  !call timing(mpi_comm_world,'CALC','PR')
  call mpibarrier()
  call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)

  if (write_matrices) then
      call write_sparse_matrix('serial_text', iproc, nproc, mpi_comm_world, smat(2), mat3(1), 'solutionmatrix_sparse')
  end if



  consistency_checks: if (do_consistency_checks) then

      ! Calculate the inverse of the desired matrix power
      if (iproc==0) then
          call yaml_comment('Calculating mat^-x',hfill='~')
          call yaml_mapping_open('Calculating mat^-x')
      end if
      call matrix_chebyshev_expansion(iproc, nproc, mpi_comm_world, &
           1, (/-expo/), smat(1), smat(2), mat2, mat3(2), ice_obj=ice_obj(1))
    
      if (iproc==0) then
          call yaml_mapping_close()
      end if
    
      ! Multiply the two resulting matrices.
      if (iproc==0) then
          call yaml_comment('Calculating mat^x*mat^-x',hfill='~')
      end if
      call matrix_matrix_multiplication(iproc, nproc, smat(2), mat3(1), mat3(2), mat3(3))
    
      ! Check the result
      call check_deviation_from_unity_sparse(iproc, smat(2), mat3(3), max_error, mean_error)
      if (iproc==0) then
          call yaml_mapping_open('Check the deviation from unity of the operation mat^x*mat^-x')
          call yaml_map('max_error',max_error,fmt='(es10.3)')
          call yaml_map('mean_error',mean_error,fmt='(es10.3)')
          call yaml_mapping_close()
      end if
    
    
      ! Calculate the inverse operation, applied to the operation itsel, i.e. (mat^x)^(1/x)
      if (iproc==0) then
          call yaml_comment('Calculating (mat^x)^(1/x)',hfill='~')
          call yaml_mapping_open('Calculating (mat^x)^(1/x)')
      end if
    
      ! Reset the ICE object
      call foe_data_deallocate(ice_obj(1))
      charge_fake = f_malloc0(1,id='charge_fake')
      call init_foe(iproc, nproc, 1, charge_fake, ice_obj(1), evlow=evlow, evhigh=evhigh, &
           betax=betax, eval_multiplicator=eval_multiplicator, &
           accuracy_function=accuracy_ice, accuracy_penalty=accuracy_penalty)
      call f_free(charge_fake)
      call matrix_chebyshev_expansion(iproc, nproc, mpi_comm_world, &
           1, (/1.0_mp/expo/), smat(2), smat(2), mat3(1), mat3(3), ice_obj=ice_obj(1))
    
      if (iproc==0) then
          call yaml_mapping_close()
      end if
    
      ! Calculate the errors
      call transform_sparse_matrix(iproc, smat(1), smat(2), SPARSE_TASKGROUP, 'small_to_large', &
                   smat_in=mat2%matrix_compr, lmat_out=mat3(2)%matrix_compr)
    
      call calculate_error(iproc, nproc, mpiworld(), smat(2), mat3(3), mat3(2), nthreshold, threshold, .false., &
           'Check the deviation from the original matrix')
    
      !call timing(mpi_comm_world,'CHECK_LINEAR','PR')
      call mpibarrier()
      call f_timing_checkpoint(ctr_name='CHECK_LINEAR',mpi_comm=mpiworld(),nproc=mpisize(),&
           gather_routine=gather_timings)

  end if consistency_checks


  cubic_check:if (do_cubic_check) then
      ! Do the operation using exact LAPACK and the dense matrices
      if (iproc==0) then
          call yaml_comment('Do the same calculation using dense LAPACK',hfill='~')
      end if
      !call operation_using_dense_lapack(iproc, nproc, smat(1)_in, mat_in)
      mat2%matrix = sparsematrix_malloc_ptr(smat(1), iaction=DENSE_FULL, id='mat2%matrix')
      mat3(3)%matrix = sparsematrix_malloc_ptr(smat(1), iaction=DENSE_FULL, id='mat3(3)%matrix')
      call matrix_power_dense_lapack(iproc, nproc, mpiworld(), blocksize_diag, blocksize_matmul, .true., &
            expo, smat(1), smat(2), mat2, mat3(3), algorithm=diag_algorithm)
      call mpibarrier()
      call f_timing_checkpoint(ctr_name='CALC_CUBIC',mpi_comm=mpiworld(),nproc=mpisize(),&
           gather_routine=gather_timings)
      if (write_matrices) then
          call write_dense_matrix(iproc, nproc, mpiworld(), smat(2), mat3(3), &
               uncompress=.false., filename='solutionmatrix_dense', binary=.false.)
          call write_dense_matrix(iproc, nproc, mpiworld(), smat(2), mat2, &
               uncompress=.false., filename='randommatrix_dense', binary=.false.)
      end if
      !call write_dense_matrix(iproc, nproc, mpi_comm_world, smat(2), mat3(1), 'resultchebyshev.dat', binary=.false.)
      !call write_dense_matrix(iproc, nproc, mpi_comm_world, smat(2), mat3(3), 'resultlapack.dat', binary=.false.)

      mat3(1)%matrix = sparsematrix_malloc0_ptr(smat(2), iaction=DENSE_FULL,id=' mat3(1)%matrix')
      call uncompress_matrix2(iproc, nproc, mpiworld(), smat(2), mat3(1)%matrix_compr, mat3(1)%matrix)

      ! Calculate the errors
      call calculate_error(iproc, nproc, mpiworld(), smat(2), mat3(1), mat3(3), nthreshold, threshold, .true., &
           'Check the deviation from the exact result using BLAS')

      !!! Sparse matrices
      !!max_error = 0.0_mp
      !!mean_error = 0.0_mp
      !!max_error_rel = 0.0_mp
      !!mean_error_rel = 0.0_mp
      !!max_error_rel_threshold(:) = 0.0_mp
      !!mean_error_rel_threshold(:) = 0.0_mp
      !!nrel_threshold(:) = 0
      !!do i=1,smat(2)%nvctr
      !!    tt = abs(mat3(1)%matrix_compr(i)-mat3(3)%matrix_compr(i))
      !!    tt_rel = tt/abs(mat3(3)%matrix_compr(i))
      !!    mean_error = mean_error + tt
      !!    max_error = max(max_error,tt)
      !!    mean_error_rel = mean_error_rel + tt_rel
      !!    max_error_rel = max(max_error_rel,tt_rel)
      !!    do ithreshold=1,nthreshold
      !!        if (abs(mat3(3)%matrix_compr(i))>threshold(ithreshold)) then
      !!            nrel_threshold(ithreshold) = nrel_threshold(ithreshold) + 1
      !!            mean_error_rel_threshold(ithreshold) = mean_error_rel_threshold(ithreshold) + tt_rel
      !!            max_error_rel_threshold(ithreshold) = max(max_error_rel_threshold(ithreshold),tt_rel)
      !!        end if
      !!    end do
      !!end do
      !!mean_error = mean_error/real(smat(2)%nvctr,kind=8)
      !!mean_error_rel = mean_error_rel/real(smat(2)%nvctr,kind=8)
      !!do ithreshold=1,nthreshold
      !!    mean_error_rel_threshold(ithreshold) = mean_error_rel_threshold(ithreshold)/real(nrel_threshold(ithreshold),kind=8)
      !!end do
      !!if (iproc==0) then
      !!    call yaml_mapping_open('Check the deviation from the exact result using BLAS (only within the sparsity pattern)')
      !!    call yaml_mapping_open('absolute error')
      !!    call yaml_map('max error',max_error,fmt='(es10.3)')
      !!    call yaml_map('mean error',mean_error,fmt='(es10.3)')
      !!    call yaml_mapping_close()
      !!    call yaml_mapping_open('relative error')
      !!    call yaml_map('max error relative',max_error_rel,fmt='(es10.3)')
      !!    call yaml_map('mean error relative',mean_error_rel,fmt='(es10.3)')
      !!    call yaml_mapping_close()
      !!    !call yaml_mapping_open('relative error with threshold')
      !!    call yaml_sequence_open('relative error with threshold')
      !!    do ithreshold=1,nthreshold
      !!        call yaml_sequence(advance='no')
      !!        call yaml_mapping_open(flow=.true.)
      !!        call yaml_map('threshold value',threshold(ithreshold),fmt='(es8.1)')
      !!        call yaml_map('max error relative',max_error_rel_threshold(ithreshold),fmt='(es10.3)')
      !!        call yaml_map('mean error relative',mean_error_rel_threshold(ithreshold),fmt='(es10.3)')
      !!        call yaml_mapping_close()
      !!    end do
      !!    call yaml_sequence_close()
      !!    call yaml_mapping_close()
      !!    call yaml_mapping_close()
      !!end if

      !!! Full matrices
      !!max_error = 0.0_mp
      !!mean_error = 0.0_mp
      !!max_error_rel = 0.0_mp
      !!mean_error_rel = 0.0_mp
      !!max_error_rel_threshold(:) = 0.0_mp
      !!mean_error_rel_threshold(:) = 0.0_mp
      !!nrel_threshold(:) = 0
      !!do i=1,smat(2)%nfvctr
      !!    do j=1,smat(2)%nfvctr
      !!        tt = abs(mat3(1)%matrix(j,i,1)-mat3(3)%matrix(j,i,1))
      !!        tt_rel = tt/abs(mat3(3)%matrix(j,i,1))
      !!        mean_error = mean_error + tt
      !!        max_error = max(max_error,tt)
      !!        mean_error_rel = mean_error_rel + tt_rel
      !!        max_error_rel = max(max_error_rel,tt_rel)
      !!        do ithreshold=1,nthreshold
      !!            if (abs(mat3(3)%matrix(j,i,1))>threshold(ithreshold)) then
      !!                nrel_threshold(ithreshold) = nrel_threshold(ithreshold) + 1
      !!                mean_error_rel_threshold(ithreshold) = mean_error_rel_threshold(ithreshold) + tt_rel
      !!                max_error_rel_threshold(ithreshold) = max(max_error_rel_threshold(ithreshold),tt_rel)
      !!            end if
      !!        end do
      !!    end do
      !!end do
      !!mean_error = mean_error/real(smat(2)%nvctr,kind=8)
      !!mean_error_rel = mean_error_rel/real(smat(2)%nvctr,kind=8)
      !!do ithreshold=1,nthreshold
      !!    mean_error_rel_threshold(ithreshold) = mean_error_rel_threshold(ithreshold)/real(nrel_threshold(ithreshold),kind=8)
      !!end do

      !!if (iproc==0) then
      !!    call yaml_mapping_open('Check the deviation from the exact result using BLAS (for the entire matrix)')
      !!    call yaml_mapping_open('absolute error')
      !!    call yaml_map('max error',max_error,fmt='(es10.3)')
      !!    call yaml_map('mean error',mean_error,fmt='(es10.3)')
      !!    call yaml_mapping_close()
      !!    call yaml_mapping_open('relative error')
      !!    call yaml_map('max error relative',max_error_rel,fmt='(es10.3)')
      !!    call yaml_map('mean error relative',mean_error_rel,fmt='(es10.3)')
      !!    call yaml_mapping_close()
      !!    !call yaml_mapping_open('relative error with threshold')
      !!    call yaml_sequence_open('relative error with threshold')
      !!    do ithreshold=1,nthreshold
      !!        call yaml_sequence(advance='no')
      !!        call yaml_mapping_open(flow=.true.)
      !!        call yaml_map('threshold value',threshold(ithreshold),fmt='(es8.1)')
      !!        call yaml_map('max error relative',max_error_rel_threshold(ithreshold),fmt='(es10.3)')
      !!        call yaml_map('mean error relative',mean_error_rel_threshold(ithreshold),fmt='(es10.3)')
      !!        call yaml_mapping_close()
      !!    end do
      !!    call yaml_sequence_close()
      !!    call yaml_mapping_close()
      !!    call yaml_mapping_close()
      !!end if

      call f_free_ptr(mat3(1)%matrix)

      !call timing(mpi_comm_world,'CHECK_CUBIC','PR')
      call mpibarrier()
      call f_timing_checkpoint(ctr_name='CHECK_CUBIC',mpi_comm=mpiworld(),nproc=mpisize(),&
           gather_routine=gather_timings)
  end if cubic_check


  ! Deallocate the sparse matrix descriptors type
  call deallocate_sparse_matrix(smat(1))
  call deallocate_sparse_matrix(smat(2))

  !! Deallocate the meta data
  !call deallocate_sparse_matrix_metadata(smmd)

  ! Deallocat the sparse matrices
  call deallocate_matrices(mat2)
  call deallocate_matrices(mat3(1))
  call deallocate_matrices(mat3(2))
  call deallocate_matrices(mat3(3))

  ! Deallocate the object holding the parameters
  call f_free(times)
  call f_free(sumarr)
  do it=1,nit
      call foe_data_deallocate(ice_obj(it))
  end do

  ! Gather the timings
  call mpibarrier()
  call build_dict_info(iproc, nproc, dict_timing_info)
  call f_timing_stop(mpi_comm=mpi_comm_world, nproc=nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  if (iproc==0) then
      call yaml_release_document()
  end if

  ! Finalize MPI
  call mpifinalize()

  ! Finalize flib
  call f_lib_finalize()



  !!!SM: This routine should go to a module
  !!contains
  !! !> construct the dictionary needed for the timing information
  !!  subroutine build_dict_info(dict_info)
  !!    !use module_base
  !!    use dynamic_memory
  !!    use dictionaries
  !!    implicit none
  !!    include 'mpif.h'
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

  contains 

    subroutine read_and_communicate_input_variables()

      ! Read in the parameters
      if (iproc==0) then
          !!call yaml_comment('Required input: nfvctr, nvctr, nbuf_large, nbuf_mult, condition_number, expo')

          parser=yaml_cl_parse_null()
          call commandline_options(parser)
          call yaml_cl_parse_cmd_line(parser,args=options)
          call yaml_cl_parse_free(parser)

          metadata_file = options//'metadata_file'
          nfvctr = options//'nfvctr'
          nvctr = options//'nvctr'
          nbuf_large = options//'nbuf_large'
          nbuf_mult = options//'nbuf_mult'
          condition_number = options//'condition_number'
          expo = options//'expo'
          infile = options//'infile'
          outfile = options//'outfile'
          outmatmulfile = options//'outmatmulfile'
          sparsegen_method = options//'sparsegen_method'
          matgen_method = options//'matgen_method'
          write_matrices = options//'write_matrices'
          betax = options//'betax'
          blocksize_diag = options//'blocksize_diag'
          blocksize_matmul = options//'blocksize_matmul'
          evlow = options//'evlow'
          evhigh = options//'evhigh'
          do_cubic_check = options//'do_cubic_check'
          diag_algorithm = options//'diag_algorithm'
          eval_multiplicator = options//'eval_multiplicator'
          solution_method = options//'solution_method'
          pexsi_np_sym_fact = options//'pexsi_np_sym_fact'
          accuracy_ice = options//'accuracy_ice'
          accuracy_penalty = options//'accuracy_penalty'
          nit = options//'nit'
          do_consistency_checks = options//'do_consistency_checks'
          output_level = options//'output_level'
          profiling_depth = options//'profiling_depth'

          call dict_free(options)
      end if

      ! Send the input parameters to all MPI tasks
      call mpibcast(metadata_file, root=0, comm=mpi_comm_world)
      call mpibcast(nfvctr, root=0, comm=mpi_comm_world)
      call mpibcast(nvctr, root=0, comm=mpi_comm_world)
      call mpibcast(nbuf_large, root=0, comm=mpi_comm_world)
      call mpibcast(nbuf_mult, root=0, comm=mpi_comm_world)
      call mpibcast(condition_number, root=0, comm=mpi_comm_world)
      call mpibcast(expo, root=0, comm=mpi_comm_world)
      call mpibcast(infile, root=0, comm=mpi_comm_world)
      call mpibcast(outfile, root=0, comm=mpi_comm_world)
      call mpibcast(outmatmulfile, root=0, comm=mpi_comm_world)
      call mpibcast(sparsegen_method, root=0, comm=mpi_comm_world)
      call mpibcast(matgen_method, root=0, comm=mpi_comm_world)
      call mpibcast(solution_method, root=0, comm=mpi_comm_world)
      call mpibcast(betax, root=0, comm=mpi_comm_world)
      call mpibcast(blocksize_diag, root=0, comm=mpi_comm_world)
      call mpibcast(blocksize_matmul, root=0, comm=mpi_comm_world)
      call mpibcast(evlow, root=0, comm=mpi_comm_world)
      call mpibcast(evhigh, root=0, comm=mpi_comm_world)
      call mpibcast(diag_algorithm, root=0, comm=mpi_comm_world)
      call mpibcast(eval_multiplicator, root=0, comm=mpi_comm_world)
      call mpibcast(pexsi_np_sym_fact, root=0, comm=mpi_comm_world)
      call mpibcast(accuracy_ice, root=0, comm=mpi_comm_world)
      call mpibcast(accuracy_penalty, root=0, comm=mpi_comm_world)
      call mpibcast(nit, root=0, comm=mpi_comm_world)

      ! Since there is no wrapper for logicals...
      if (iproc==0) then
          if (write_matrices) then
              iwrite = 1
          else
              iwrite = 0
          end if
          if (do_cubic_check) then
              icheck = 1
          else
              icheck = 0
          end if
          if (do_consistency_checks) then
              icons = 1
          else
              icons = 0
          end if
      end if
      call mpibcast(iwrite, root=0, comm=mpi_comm_world)
      call mpibcast(icheck, root=0, comm=mpi_comm_world)
      call mpibcast(icons, root=0, comm=mpi_comm_world)
      if (iwrite==1) then
          write_matrices = .true.
      else
          write_matrices = .false.
      end if
      if (icheck==1) then
          do_cubic_check = .true.
      else
          do_cubic_check = .false.
      end if
      if (icons==1) then
          do_consistency_checks = .true.
      else
          do_consistency_checks = .false.
      end if

    end subroutine read_and_communicate_input_variables

end program driver_random



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

  call yaml_cl_parse_option(parser,'nfvctr','0',&
       'matrix size',&
       help_dict=dict_new('Usage' .is. &
       'Size of the matrix (number of rows/columns)',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'nvctr','0',&
       'nonzero entries',&
       help_dict=dict_new('Usage' .is. &
       'Number of nonzero entries of the matrix',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'nbuf_large','0',&
       'buffer for large matrix',&
       help_dict=dict_new('Usage' .is. &
       'Number of buffer elements around the sparisity pattern to create the large sparsity pattern',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'nbuf_mult','0',&
       'buffer for matrix multiplications',&
       help_dict=dict_new('Usage' .is. &
       'Number of buffer elements around the sparisity pattern to create the matrix multiplication sparsity pattern',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'condition_number','1.0',&
       'condition number',&
       help_dict=dict_new('Usage' .is. &
       'Target condition number of the random matrix',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'expo','1.0',&
       'exponent',&
       help_dict=dict_new('Usage' .is. &
       'Exponent for the matrix function to be calculated (M^expo)',&
       'Allowed values' .is. &
       'Double'))
   
  call yaml_cl_parse_option(parser,'infile','infile.dat',&
       'input file',&
       help_dict=dict_new('Usage' .is. &
       'File containing the input matrix descriptors',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'outfile','outfile.dat',&
       'output file',&
       help_dict=dict_new('Usage' .is. &
       'File containing the output matrix descriptors',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'outmatmulfile','outmatmulfile.dat',&
       'output matrix multiplication file',&
       help_dict=dict_new('Usage' .is. &
       'File containing the output matrix multiplication descriptors',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'sparsegen_method','file',&
       'sparsity pattern generation',&
       help_dict=dict_new('Usage' .is. &
       'Indicate whether the sparsity patterns should be created randomly or read from files',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'matgen_method','file',&
       'matrix content generation',&
       help_dict=dict_new('Usage' .is. &
       'Indicate whether the matrix contents should be created randomly or read from files',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'write_matrices','.false.',&
       'write the matrices to disk',&
       help_dict=dict_new('Usage' .is. &
       'Indicate whether the sparse matrices shall be written to disk',&
       'Allowed values' .is. &
       'Logical'))

  call yaml_cl_parse_option(parser,'betax','-500.0',&
       'betax for the penalty function',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the betax value, which is used in the exponential of the penalty function',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'blocksize_diag','-1',&
       'blocksize for ScaLAPACK diagonalization (negative for standard LAPACK)',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the blocksize to be used by ScaLAPACK for the diagonalization.'&
       &' If negative, then the standard LAPACK routines will be used',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'blocksize_matmul','-1',&
       'blocksize for ScaLAPACK matrix multiplication (negative for standard LAPACK)',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the blocksize to be used by ScaLAPACK for the matrix multiplication.'&
       &' If negative, then the standard LAPACK routines will be used',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'evlow','0.5',&
       'guess for the lowest matrix eigenvalue',&
       help_dict=dict_new('Usage' .is. &
       'Indicate a guess for the lowest eigenvalue of the matrix',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'evhigh','1.5',&
       'guess for the highest matrix eigenvalue',&
       help_dict=dict_new('Usage' .is. &
       'Indicate a guess for the highest eigenvalue of the matrix',&
       'Allowed values' .is. &
       'Double'))

   call yaml_cl_parse_option(parser,'do_cubic_check','.true.',&
       'perform a check using cubic scaling dense (Sca)LAPACK',&
       help_dict=dict_new('Usage' .is. &
       'Indicate whether a cubic scaling check using dense (Sca)LAPACK should be performed',&
       'Allowed values' .is. &
       'Logical'))

   call yaml_cl_parse_option(parser,'diag_algorithm','pdsyevx',&
       'ScaLAPACK algorithm to be used for the diagonalization (pdsyevx, pdsyevd)',&
       help_dict=dict_new('Usage' .is. &
       'ScaLAPACK algorithm to be used for the diagonalization: pdsyevx or pdsyevd',&
       'Allowed values' .is. &
       'String'))

   call yaml_cl_parse_option(parser,'eval_multiplicator','1.0',&
       'scale the matrix by this factor',&
       help_dict=dict_new('Usage' .is. &
       'scale the matrix by this factor to get a spectrum which is asier representable using the Chebyshe polynomials',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'solution_method','ICE',&
       'Indicate which solution method should be used (ICE or SelInv)',&
       help_dict=dict_new('Usage' .is. &
       'Indicate which solution method should be used (ICE or SelInv)',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'pexsi_np_sym_fact','16',&
       'Indicate the number of tasks used for the symbolic factorization within PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the number of tasks used for the symbolic factorization within PEXSI',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'accuracy_ice','1.e-8',&
       'Required accuracy for the Chebyshev fit for ICE',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the required accuracy for the Chebyshev fit for ICE',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'accuracy_penalty','1.e-5',&
       'Required accuracy for the Chebyshev fit for the penalty function',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the required accuracy for the Chebyshev fit for teh penalty function',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'nit','1',&
       'Number of iterations for the kernel caluculations (can be used for benchmarking)',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the number of iterations for the kernel caluculations (can be used for benchmarking)',&
       'Allowed values' .is. &
       'Integer'))

   call yaml_cl_parse_option(parser,'do_consistency_checks','.false.',&
       'Perform consistency checks of the result (done with CheSS): mat^x*mat^-x and (mat^x)^(1/x)',&
       help_dict=dict_new('Usage' .is. &
       'Perform consistency checks of the result (done with CheSS): mat^x*mat^-x and (mat^x)^(1/x)',&
       'Allowed values' .is. &
       'Logical'))

  call yaml_cl_parse_option(parser,'output_level','0',&
       'Output level of the routine profiling',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the output level of the routine profiling',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'profiling_depth','-1',&
       'Depth of the individual routine timing profiling',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the depth of the individual routine timing profiling',&
       'Allowed values' .is. &
       'Integer'))

end subroutine commandline_options
