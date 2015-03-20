program smatmul
  use bigdft_run, only: bigdft_init
  use module_base
  use yaml_parse, only: yaml_cl_parse, yaml_cl_parse_null, yaml_cl_parse_free, yaml_cl_parse_cmd_line
  !use dictionaries
  use yaml_output
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               matrices_null, deallocate_sparse_matrix, deallocate_matrices, &
                               assignment(=), sparsematrix_malloc_ptr, sparsematrix_malloc, SPARSE_FULL, SPARSEMM_SEQ
  !use sparsematrix_init, only: read_ccs_format, ccs_to_sparsebigdft, ccs_values_to_bigdft, &
  !                             read_bigdft_format, bigdft_to_sparsebigdft
  use sparsematrix_init, only: bigdft_to_sparsebigdft
  use sparsematrix, only: write_matrix_compressed, check_symmetry, &
                          write_sparsematrix_CCS, write_sparsematrix, &
                          sparsemm_new, sequential_acces_matrix_fast2
  !use matrix_operations, only: overlapPowerGeneral
  use io, only: read_sparse_matrix
  implicit none

  external :: gather_timings

  ! Variables
  integer :: iproc, nproc, ncol, nnonzero, nseg, ncolp, iscol, ierr, nspin, nfvctr, nvctr, isfvctr, nfvctrp
  !character(len=*),parameter :: filename='matrix.dat'
  character(len=1024) :: filename
  integer,dimension(:),pointer :: col_ptr, row_ind, keyv
  integer,dimension(:,:,:),pointer :: keyg
  real(kind=8),dimension(:),pointer :: val
  type(sparse_matrix) :: smat
  type(matrices) :: matA
  type(matrices),dimension(1) :: matB
  real(kind=8) :: max_error, mean_error
  logical :: symmetric
  real(kind=8) :: time_start, time_end
  real(kind=8),dimension(:),pointer :: mat_compr
  real(kind=8),dimension(:),allocatable :: mat_seq, vectors
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
  call dict_free(options)


  call bigdft_init()!mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  iproc=bigdft_mpi%iproc!mpi_info(1)
  nproc=bigdft_mpi%nproc!mpi_info(2)

  call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=.true. .and. nproc>1)

  ! Nullify all pointers
  !!nullify(col_ptr)
  !!nullify(row_ind)
  nullify(keyv)
  nullify(keyg)
  nullify(mat_compr)
  !nullify(val)

  
  ! Read in a file in the sparse BigDFT format
  !call read_bigdft_format(filename, ncol, nnonzero, nseg, keyv, keyg, val)
  call read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, mat_compr)

  ! Create the corresponding BigDFT sparsity pattern
  !call ccs_to_sparsebigdft(iproc, nproc, ncol, ncol, 0, nnonzero, row_ind, col_ptr, smat)
  !call distribute_columns_on_processes(iproc, nproc, ncol, ncolp, iscol)
  call distribute_columns_on_processes(iproc, nproc, nfvctr, nfvctrp, isfvctr)
  !write(*,'(a,4i9)') 'iproc, nfvctr, nfvctrp, isfvctr', iproc, nfvctr, nfvctrp, isfvctr
  !call bigdft_to_sparsebigdft(iproc, nproc, ncol, ncolp, iscol, nnonzero, nseg, keyg, smat)
  call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nfvctrp, isfvctr, nvctr, nseg, keyg, smat)

  matA = matrices_null()
  !matB(1) = matrices_null()

  ! Assign the values
  matA%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matA%matrix_compr')
  !call ccs_values_to_bigdft(ncol, nnonzero, row_ind, col_ptr, smat, val, matA)
  matA%matrix_compr = mat_compr

  ! Check the symmetry
  symmetric = check_symmetry(nfvctr, smat)
  if (.not.symmetric) stop 'ERROR not symmetric'
  
  ! Write the original matrix
  !if (iproc==0) call write_matrix_compressed('Original matrix', smat, matA)
  !if (iproc==0) call write_sparsematrix_CCS('original_css.dat', smat, matA)
  if (iproc==0) call write_sparsematrix('original_bigdft.dat', smat, matA)

  call timing(bigdft_mpi%mpi_comm,'INIT','PR')

  ! Calculate the inverse
  !matB(1)%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matB(1)%matrix_compr')
  call mpibarrier(bigdft_mpi%mpi_comm)
  time_start = mpi_wtime()
  !!call overlapPowerGeneral(iproc, nproc, 1020, 1, (/1/), -8, &
  !!     1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
  !!     check_accur=.true., max_error=max_error, mean_error=mean_error)

  mat_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='mat_seq')
  vectors = f_malloc0(smat%smmm%nvctrp,id='vectors')
  call sequential_acces_matrix_fast2(smat, mat_compr, mat_seq)
  call sparsemm_new(smat, mat_seq, mat_compr, vectors)


  !!call overlapPowerGeneral(iproc, nproc, 0, 1, (/1/), -8, &
  !!     1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
  !!     check_accur=.true., max_error=max_error, mean_error=mean_error)
  call mpibarrier(bigdft_mpi%mpi_comm)
  time_end = mpi_wtime()
  call timing(bigdft_mpi%mpi_comm,'CALC','PR')
  !if (iproc==0) write(*,*) 'walltime',time_end-time_start

  ! Write the results
  !if (iproc==0) call write_matrix_compressed('Result', smat, matB(1))

  !if (iproc==0) call write_sparsematrix_CCS('result_css.dat', smat, matB(1))
  !if (iproc==0) call write_sparsematrix('result_bigdft.dat', smat, matB(1))

  ! Deallocations
  call deallocate_sparse_matrix(smat)
  call deallocate_matrices(matA)
  !call deallocate_matrices(matB(1))
  call f_free_ptr(col_ptr)
  call f_free_ptr(row_ind)
  call f_free_ptr(keyv)
  call f_free_ptr(keyg)
  call f_free_ptr(val)
  call f_free(mat_seq)
  call f_free_ptr(mat_compr)
  call f_free(vectors)

  call timing(bigdft_mpi%mpi_comm,'FINISH','PR')

  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm, nproc=bigdft_mpi%nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  !if (iproc == 0) then
  !   !call yaml_comment('Timing for root process',hfill='-')
  !   !call yaml_mapping_open('Timings for root process')
  !   !call yaml_map('CPU time (s)',tcpu1-tcpu0,fmt='(f12.2)')
  !   !call yaml_map('Elapsed time (s)',tel,fmt='(f12.2)')
  !   call yaml_mapping_close()
  !   call yaml_flush_document()
  !end if


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

end subroutine commandline_options
