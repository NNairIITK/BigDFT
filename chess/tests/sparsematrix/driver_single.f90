!> @file
!! Test of the sparsematrix library
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


program driver_single
  ! The following module are part of the sparsematrix library
  use sparsematrix_base
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft, &
                                    sparse_matrix_init_from_file_bigdft, &
                                    matrices_init, &
                                    matrix_chebyshev_expansion, &
                                    matrix_matrix_multiplication
  use sparsematrix, only: check_deviation_from_unity_sparse, operation_using_dense_lapack
  implicit none

  ! External routines
  external :: gather_timings

  ! Variables
  type(sparse_matrix) :: smat_in, smat_out
  type(matrices) :: mat_in
  type(matrices),dimension(1) :: mat_out
  type(matrices),dimension(2) :: mat_check_accur
  integer :: norder_polynomial, ierr, nthread, blocksize, iproc, nproc
  real(mp) :: exp_power, max_error, mean_error
  character(len=200) :: filename_in, filename_out
  type(dictionary), pointer :: dict_timing_info
  !$ integer :: omp_get_max_threads

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

  if (iproc==0) then
      call yaml_new_document()
  end if
  call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=.true. .and. nproc>1)


  if (iproc==0) then
      call yaml_scalar('',hfill='~')
      call yaml_scalar('DRIVER FOR MATRIX CHEBYSHEV EXPANSION',hfill='~')
  end if

  if (iproc==0) then
      call yaml_mapping_open('Parallel environment')
      call yaml_map('MPI tasks',nproc)
      nthread = 1
      !$ nthread = omp_get_max_threads()
      call yaml_map('OpenMP threads',nthread)
      call yaml_mapping_close()
  end if



  ! Read in the parameters for the run:
  ! 1) name of the file with the input matrix
  ! 2) name of the file with the descriptors for the output matrix
  ! 3) exponent for the operation (mat**exponent)
  ! 4) the degree of the polynomial that shall be used
  if (iproc==0) then
       read(*,*) filename_in, filename_out, exp_power, norder_polynomial

      call yaml_mapping_open('Input parameters')
      call yaml_map('File with input matrix',trim(filename_in))
      call yaml_map('File with output matrix descriptors',trim(filename_out))
      call yaml_map('Exponent',exp_power)
      call yaml_map('Polynomial degree',norder_polynomial)
      call yaml_mapping_close()
  end if

  ! Send the input parameters to all MPI tasks
  call mpibcast(filename_in, root=0, comm=mpi_comm_world)
  call mpibcast(filename_out, root=0, comm=mpi_comm_world)
  call mpibcast(exp_power, root=0, comm=mpi_comm_world)
  call mpibcast(norder_polynomial, root=0, comm=mpi_comm_world)

  ! Read the input matrix descriptors and the matrix itself, and create the correpsonding structures
  if (iproc==0) then
      call yaml_comment('Input matrix',hfill='-')
      call yaml_mapping_open('Input matrix structure')
  end if
  call sparse_matrix_and_matrices_init_from_file_bigdft(filename_in, &
       iproc, nproc, mpi_comm_world, smat_in, mat_in)
  if (iproc==0) then
      call yaml_mapping_close()
  end if
  mat_in%matrix_compr(:) = 100.d0*mat_in%matrix_compr(:)

  ! Read the output matrix descriptors, and create the correpsonding structures
  if (iproc==0) then
      call yaml_comment('Output matrix',hfill='-')
      call yaml_mapping_open('Output matrix structure')
  end if
  call sparse_matrix_init_from_file_bigdft(filename_out, &
       iproc, nproc, mpi_comm_world, smat_out)
  if (iproc==0) then
      call yaml_mapping_close()
  end if

  ! Prepare the output matrix structure
  call matrices_init(smat_out, mat_out(1))

  if (iproc==0) then
      call yaml_comment('Performing Matrix Chebyshev Expansion',hfill='-')
  end if

  !call timing(mpi_comm_world,'INIT','PR')
  call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)

  ! Perform the operation mat_out = mat_in**exp_power
  call matrix_chebyshev_expansion(iproc, nproc, mpi_comm_world, &
       1, (/exp_power/), &
       smat_in, smat_out, mat_in, mat_out, npl_auto=.true.)

  !call timing(mpi_comm_world,'CALC_LINEAR','PR')
  call f_timing_checkpoint(ctr_name='CALC_LINEAR',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)


  ! Now perform a check of the accuracy by calculating the matrix with the inverse power. The previously
  ! calculated matrix times this result should the give the identity.
  ! Prepare the structures needed for the check of the accuracy
  call matrices_init(smat_out, mat_check_accur(1))
  call matrices_init(smat_out, mat_check_accur(2))

  ! Perform the operation mat_out = mat_in**exp_power
  if (iproc==0) then
      call yaml_comment('Performing Matrix Chebyshev Expansion',hfill='-')
  end if
  call matrix_chebyshev_expansion(iproc, nproc, mpi_comm_world, &
       1, (/-exp_power/), &
       smat_in, smat_out, mat_in, mat_check_accur(1), npl_auto=.true.)
  ! Multiply the previously calculated one with this new none. The result should be the identity.
  call matrix_matrix_multiplication(iproc, nproc, smat_out, &
       mat_out(1), mat_check_accur(1), mat_check_accur(2))

  ! Check the deviation from unity
  if (iproc==0) then
      call yaml_comment('Check deviation from Unity')
  end if
  call check_deviation_from_unity_sparse(iproc, smat_out, mat_check_accur(2), &
       max_error, mean_error)
  if (iproc==0) then
      call yaml_mapping_open('Check the deviation from unity of the operation mat^x*mat^-x')
      call yaml_map('max_error',max_error,fmt='(es10.3)')
      call yaml_map('mean_error',mean_error,fmt='(es10.3)')
      call yaml_mapping_close()
  end if


  !call timing(mpi_comm_world,'CHECK_LINEAR','PR')
  call f_timing_checkpoint(ctr_name='CHECK_LINEAR',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)


  ! Do the operation using exact LAPACK and the dense matrices
  if (iproc==0) then
      call yaml_comment('Do the same calculation using dense LAPACK',hfill='-')
  end if
  call operation_using_dense_lapack(iproc, nproc, exp_power, smat_in, mat_in)

  ! Deallocate all structures
  call deallocate_sparse_matrix(smat_in)
  call deallocate_sparse_matrix(smat_out)
  call deallocate_matrices(mat_in)
  call deallocate_matrices(mat_out(1))
  call deallocate_matrices(mat_check_accur(1))
  call deallocate_matrices(mat_check_accur(2))

  !call timing(mpi_comm_world,'FINISH','PR')
  call f_timing_checkpoint(ctr_name='FINISH',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)


  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=mpi_comm_world, nproc=nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  ! Finalize MPI
  call mpifinalize()

  ! Finalize flib
  ! SM: I have the impression that every task should call this routine, but if I do so
  ! some things are printed nproc times instead of once.
  if (iproc==0) call f_lib_finalize()


  contains

   !> Construct the dictionary needed for the timing information
    subroutine build_dict_info(dict_info)
      !use module_base
      use dynamic_memory
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict_info
      !local variables
      integer :: ierr,namelen,nthreads
      character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
      character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
      type(dictionary), pointer :: dict_tmp
      !$ integer :: omp_get_max_threads

      call dict_init(dict_info)
!  bastian: comment out 4 followinf lines for debug purposes (7.12.2014)
      !if (DoLastRunThings) then
         call f_malloc_dump_status(dict_summary=dict_tmp)
         call set(dict_info//'Routines timing and number of calls',dict_tmp)
      !end if
      nthreads = 0
      !$  nthreads=omp_get_max_threads()
      call set(dict_info//'CPU parallelism'//'MPI tasks',nproc)
      if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
           nthreads)

      nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0.to.nproc-1,id='nodename')
      if (nproc>1) then
         call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
         !gather the result between all the process
         call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
              nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
              mpi_comm_world,ierr)
         if (iproc==0) call set(dict_info//'Hostnames',&
                 list_new(.item. nodename))
      end if
      call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)

    end subroutine build_dict_info





end program driver_single
