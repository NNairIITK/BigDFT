program smatmul
  use bigdft_run, only: bigdft_init
  use module_base
  use yaml_parse, only: yaml_cl_parse, yaml_cl_parse_null, yaml_cl_parse_free, yaml_cl_parse_cmd_line
  !use dictionaries
  use yaml_output
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               matrices_null, deallocate_sparse_matrix, deallocate_matrices, &
                               assignment(=), sparsematrix_malloc_ptr, sparsematrix_malloc, SPARSE_FULL, SPARSEMM_SEQ, &
                               SPARSE_MATMUL_SMALL
  !use sparsematrix_init, only: read_ccs_format, ccs_to_sparsebigdft, ccs_values_to_bigdft, &
  !                             read_bigdft_format, bigdft_to_sparsebigdft
  use sparsematrix_init, only: bigdft_to_sparsebigdft
  use sparsematrix, only: write_matrix_compressed, check_symmetry, &
                          write_sparsematrix_CCS, write_sparsematrix, &
                          sparsemm_new, sequential_acces_matrix_fast2, &
                          compress_matrix_distributed_wrapper
  !use matrix_operations, only: overlapPowerGeneral
  use io, only: read_sparse_matrix
  implicit none

  external :: gather_timings

  ! Variables
  integer :: iproc, nproc, ncol, nnonzero, nseg, ncolp, iscol, ierr, nspin, nfvctr, nvctr, isfvctr, nfvctrp, nit, it, verbosity
  !character(len=*),parameter :: filename='matrix.dat'
  character(len=1024) :: filename
  integer,dimension(:),pointer :: col_ptr, row_ind, keyv, on_which_atom
  integer,dimension(:,:,:),pointer :: keyg
  type(sparse_matrix) :: smat
  type(matrices) :: matA
  type(matrices),dimension(1) :: matB
  real(kind=8) :: max_error, mean_error
  logical :: symmetric
  real(kind=8) :: time_start, time_end
  real(kind=8),dimension(:),pointer :: mat_compr
  real(kind=8),dimension(:),allocatable :: mat_seq, vector_in, vector_out
  type(dictionary), pointer :: dict_timing_info
  type(dictionary), pointer :: options
  type(yaml_cl_parse) :: parser !< command line parser

  ! Initialize
  call f_lib_initialize()


  parser=yaml_cl_parse_null()
  call commandline_options(parser)
  call yaml_cl_parse_cmd_line(parser,args=options)
  call yaml_cl_parse_free(parser)

  filename = options//'filename'
  nit = options//'nit'
  verbosity = options//'verbosity'
  call dict_free(options)

  if (verbosity<1 .or. verbosity>2) then
      call f_err_throw('wrong value for the verbosity, only 1 or 2 is allowed', &
          err_name='GENERIC_ERROR')
  end if

  call bigdft_init()!mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  iproc=bigdft_mpi%iproc!mpi_info(1)
  nproc=bigdft_mpi%nproc!mpi_info(2)

  if (iproc==0) then
      call yaml_new_document()
  end if

  call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=.true. .and. nproc>1)

  ! Nullify all pointers
  nullify(keyv)
  nullify(keyg)
  nullify(mat_compr)
  nullify(on_which_atom)

  
  ! Read in a file in the sparse BigDFT format
  call read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, mat_compr, on_which_atom=on_which_atom)

  ! Create the corresponding BigDFT sparsity pattern
  call distribute_columns_on_processes(iproc, nproc, nfvctr, nfvctrp, isfvctr)
  call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nfvctrp, isfvctr, on_which_atom, nvctr, nseg, keyg, smat)

  matA = matrices_null()

  ! Assign the values
  matA%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matA%matrix_compr')
  matA%matrix_compr = mat_compr

  ! Check the symmetry
  symmetric = check_symmetry(nfvctr, smat)
  if (.not.symmetric) stop 'ERROR not symmetric'
  
  ! Write the original matrix
  if (iproc==0) call write_sparsematrix('original_bigdft.dat', smat, matA)

  call timing(bigdft_mpi%mpi_comm,'INIT','PR')

  ! Calculate the inverse
  call mpibarrier(bigdft_mpi%mpi_comm)
  time_start = mpi_wtime()

  mat_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='mat_seq')
  vector_in = f_malloc0(smat%smmm%nvctrp,id='vector_in')
  vector_out = f_malloc0(smat%smmm%nvctrp,id='vector_out')
  call sequential_acces_matrix_fast2(smat, mat_compr, mat_seq)

  call vcopy(smat%smmm%nvctrp, mat_compr(smat%smmm%isvctr_mm_par(iproc)+1), 1, vector_in(1), 1)
  do it=1,nit
      call sparsemm_new(smat, mat_seq, vector_in, vector_out)
      call vcopy(smat%smmm%nvctrp, vector_out(1), 1, vector_in(1), 1)
  end do


  call compress_matrix_distributed_wrapper(iproc, nproc, smat, SPARSE_MATMUL_SMALL, &
       vector_out, matA%matrix_compr)
  if (iproc==0 .and. verbosity==2) then
      call write_matrix_compressed('final result', smat, matA)
  end if

  call mpibarrier(bigdft_mpi%mpi_comm)
  time_end = mpi_wtime()
  call timing(bigdft_mpi%mpi_comm,'CALC','PR')

  ! Deallocations
  call deallocate_sparse_matrix(smat)
  call deallocate_matrices(matA)
  call f_free_ptr(keyv)
  call f_free_ptr(keyg)
  call f_free(mat_seq)
  call f_free_ptr(mat_compr)
  call f_free(vector_in)
  call f_free(vector_out)
  call f_free_ptr(on_which_atom)

  call timing(bigdft_mpi%mpi_comm,'FINISH','PR')

  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm, nproc=bigdft_mpi%nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  if (iproc==0) then
      call yaml_release_document()
  end if

  call bigdft_finalize(ierr)

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

end program smatmul


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



subroutine commandline_options(parser)
  use yaml_parse
  use dictionaries, only: dict_new,operator(.is.)
  implicit none
  type(yaml_cl_parse),intent(inout) :: parser

  call yaml_cl_parse_option(parser,'filename','matrix.dat',&
       'input file name','f',&
       dict_new('Usage' .is. &
       'File name from which the sparse matrix is read',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'nit','1',&
       'number of iterations','n',&
       dict_new('Usage' .is. &
       'Number of matrix matrix multiplications to be performed',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'verbosity','2',&
       'verbosity of the output','v',&
       dict_new('Usage' .is. &
       'If the verbosity is high, the final result will be printed to the scree',&
       'Allowed values' .is. &
       'Integer. Only 1 or 2 is possible'))

end subroutine commandline_options
