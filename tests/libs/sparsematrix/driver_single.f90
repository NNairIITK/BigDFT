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
                                    matrix_chebyshev_expansion
  implicit none

  ! External routines
  external :: gather_timings

  ! Variables
  type(sparse_matrix) :: smat_in, smat_out
  type(matrices) :: mat_in
  type(matrices),dimension(1) :: mat_out
  integer :: norder_polynomial, ierr
  real(kind=8) :: exp_power
  character(len=200) :: filename_in, filename_out
  type(dictionary), pointer :: dict_timing_info

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

  ! Read in the parameters for the run:
  ! 1) name of the file with the input matrix
  ! 2) name of the file with the descriptors for the output matrix
  ! 3) exponent for the operation (mat**exponent)
  ! 4) the degree of the polynomial that shall be used
  read(*,*) filename_in, filename_out, exp_power, norder_polynomial

  write(*,*) 'filename_in',filename_in

  ! Read the input matrix descriptors and the matrix itself, and create the correpsonding structures
  call sparse_matrix_and_matrices_init_from_file_bigdft(filename_in, bigdft_mpi%iproc, bigdft_mpi%nproc, smat_in, mat_in)

  ! Read the output matrix descriptors, and create the correpsonding structures
  call sparse_matrix_init_from_file_bigdft(filename_out, bigdft_mpi%iproc, bigdft_mpi%nproc, smat_out)

  ! Prepare the output matrix structure
  call matrices_init(smat_out, mat_out(1))

  call timing(bigdft_mpi%mpi_comm,'INIT','PR')

  ! Perform the operation mat_out = mat_in**exp_power
  call matrix_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, norder_polynomial, 1, (/exp_power/), &
       smat_in, smat_out, mat_in, mat_out)

  call timing(bigdft_mpi%mpi_comm,'CALC','PR')

  ! Deallocate all structures
  call deallocate_sparse_matrix(smat_in)
  call deallocate_sparse_matrix(smat_out)
  call deallocate_matrices(mat_in)
  call deallocate_matrices(mat_out(1))

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


end program driver_single
