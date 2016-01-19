program driver_single
  use module_base
  use yaml_output
  use bigdft_run, only: bigdft_init
  use sparsematrix_base, only: sparse_matrix, &
                               matrices, &
                               deallocate_sparse_matrix, &
                               deallocate_matrices
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft, &
                                    sparse_matrix_init_from_file_bigdft, &
                                    matrices_init, &
                                    matrix_chebyshev_expansion, &
                                    matrix_matrix_multiplication
  implicit none

  ! External routines
  external :: gather_timings

  ! Variables
  type(sparse_matrix) :: smat_in, smat_out
  type(matrices) :: mat_in
  type(matrices),dimension(1) :: mat_out
  type(matrices),dimension(2) :: mat_check_accur
  integer :: norder_polynomial, ierr, nthread
  real(kind=8) :: exp_power
  character(len=200) :: filename_in, filename_out
  type(dictionary), pointer :: dict_timing_info
  !$ integer :: omp_get_max_threads

  ! Initialize flib
  call f_lib_initialize()

  ! General initialization, including MPI. Here we have:
  ! bigdft_mpi%iproc is the task ID
  ! bigdft_mpi%nproc is the total number of tasks
  ! PROBLEM: This is in src/modules...
  call bigdft_init()

  if (bigdft_mpi%iproc==0) then
      call yaml_new_document()
  end if
  call f_timing_reset(filename='time.yaml',master=bigdft_mpi%iproc==0,verbose_mode=.true. .and. bigdft_mpi%nproc>1)


  if (bigdft_mpi%iproc==0) then
      call yaml_scalar('',hfill='~')
      call yaml_scalar('DRIVER FOR MATRIX CHEBYSHEV EXPANSION',hfill='~')
  end if

  if (bigdft_mpi%iproc==0) then
      call yaml_mapping_open('Parallel environment')
      call yaml_map('MPI tasks',bigdft_mpi%nproc)
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
  read(*,*) filename_in, filename_out, exp_power, norder_polynomial

  if (bigdft_mpi%iproc==0) then
      call yaml_mapping_open('Input parameters')
      call yaml_map('File with input matrix',trim(filename_in))
      call yaml_map('File with output matrix descriptors',trim(filename_out))
      call yaml_map('Exponent',exp_power)
      call yaml_map('Polynomial degree',norder_polynomial)
      call yaml_mapping_close()
  end if

  ! Read the input matrix descriptors and the matrix itself, and create the correpsonding structures
  if (bigdft_mpi%iproc==0) then
      call yaml_comment('Input matrix',hfill='-')
      call yaml_mapping_open('Input matrix structure')
  end if
  call sparse_matrix_and_matrices_init_from_file_bigdft(filename_in, bigdft_mpi%iproc, bigdft_mpi%nproc, smat_in, mat_in)
  if (bigdft_mpi%iproc==0) then
      call yaml_mapping_close()
  end if
  mat_in%matrix_compr(:) = 100.d0*mat_in%matrix_compr(:)

  ! Read the output matrix descriptors, and create the correpsonding structures
  if (bigdft_mpi%iproc==0) then
      call yaml_comment('Output matrix',hfill='-')
      call yaml_mapping_open('Output matrix structure')
  end if
  call sparse_matrix_init_from_file_bigdft(filename_out, bigdft_mpi%iproc, bigdft_mpi%nproc, smat_out)
  if (bigdft_mpi%iproc==0) then
      call yaml_mapping_close()
  end if

  ! Prepare the output matrix structure
  call matrices_init(smat_out, mat_out(1))

  if (bigdft_mpi%iproc==0) then
      call yaml_comment('Performing Matrix Chebyshev Expansion',hfill='-')
  end if

  call timing(bigdft_mpi%mpi_comm,'INIT','PR')

  ! Perform the operation mat_out = mat_in**exp_power
  call matrix_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, norder_polynomial, 1, (/exp_power/), &
       smat_in, smat_out, mat_in, mat_out, npl_auto=.true.)

  call timing(bigdft_mpi%mpi_comm,'CALC','PR')

  ! Now perform a check of the accuracy by calculating the matrix with the inverse power. The previously
  ! calculated matrix times this result should the give the identity.
  ! Prepare the structures needed for the check of the accuracy
  call matrices_init(smat_out, mat_check_accur(1))
  call matrices_init(smat_out, mat_check_accur(2))

  ! Perform the operation mat_out = mat_in**exp_power
  if (bigdft_mpi%iproc==0) then
      call yaml_comment('Performing Matrix Chebyshev Expansion',hfill='-')
  end if
  call matrix_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, norder_polynomial, 1, (/-exp_power/), &
       smat_in, smat_out, mat_in, mat_check_accur(1), npl_auto=.true.)
  ! Multiply the previously calculated one with this new none. The result should be the identity.
  call matrix_matrix_multiplication(bigdft_mpi%iproc, bigdft_mpi%nproc, smat_out, &
       mat_out(1), mat_check_accur(1), mat_check_accur(2))

  ! Check the deviation from unity
  if (bigdft_mpi%iproc==0) then
      call yaml_comment('Check deviation from Unity')
  end if
  call check_deviation_from_unity(smat_out, mat_check_accur(2))

  ! Deallocate all structures
  call deallocate_sparse_matrix(smat_in)
  call deallocate_sparse_matrix(smat_out)
  call deallocate_matrices(mat_in)
  call deallocate_matrices(mat_out(1))
  call deallocate_matrices(mat_check_accur(1))
  call deallocate_matrices(mat_check_accur(2))

  call timing(bigdft_mpi%mpi_comm,'FINISH','PR')

  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm, nproc=bigdft_mpi%nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)


  ! Finalize MPI
  call bigdft_finalize(ierr)

  ! Finalize flib
  call f_lib_finalize()


  contains
   !> construct the dictionary needed for the timing information
    subroutine build_dict_info(dict_info)
      !use module_base
      use dynamic_memory
      use dictionaries
      implicit none
      include 'mpif.h'
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
      call set(dict_info//'CPU parallelism'//'MPI tasks',bigdft_mpi%nproc)
      if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
           nthreads)

      nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0.to.bigdft_mpi%nproc-1,id='nodename')
      if (bigdft_mpi%nproc>1) then
         call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
         !gather the result between all the process
         call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
              nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
              bigdft_mpi%mpi_comm,ierr)
         if (bigdft_mpi%iproc==0) call set(dict_info//'Hostnames',&
                 list_new(.item. nodename))
      end if
      call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)

    end subroutine build_dict_info


    subroutine check_deviation_from_unity(smat, mat)
      use sparsematrix_base, only: sparse_matrix, &
                                   matrices
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat

      ! Local variables
      integer :: iseg, ii, i, irow, icolumn
      real(kind=8) :: sum_error, max_error, error

      sum_error = 0.d0
      max_error = 0.d0
      do iseg=1,smat%nseg
          ii=smat%keyv(iseg)
          ! A segment is always on one line, therefore no double loop
          do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
             irow = i
             icolumn = smat%keyg(1,2,iseg)
             if (irow==icolumn) then
                 error = abs(mat%matrix_compr(ii)-1.d0)
             else
                 error = abs(mat%matrix_compr(ii))
             end if
             sum_error = sum_error + error
             max_error = max(max_error,error)
             ii=ii+1
         end do
      end do

      call yaml_mapping_open('Check the deviation from unity of the operation S^x*S^-x')
      call yaml_map('max_error',max_error,fmt='(es10.3)')
      call yaml_map('sum_error',sum_error/real(smat%nvctr,kind=8),fmt='(es10.3)')
      call yaml_mapping_close()

    end subroutine check_deviation_from_unity



end program driver_single
