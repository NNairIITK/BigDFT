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
  integer :: norder_polynomial, ierr, nthread, blocksize
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

  call mpiinit()

  !iproc=mpirank()

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
  if (bigdft_mpi%iproc==0) then
       read(*,*) filename_in, filename_out, exp_power, norder_polynomial

      call yaml_mapping_open('Input parameters')
      call yaml_map('File with input matrix',trim(filename_in))
      call yaml_map('File with output matrix descriptors',trim(filename_out))
      call yaml_map('Exponent',exp_power)
      call yaml_map('Polynomial degree',norder_polynomial)
      call yaml_mapping_close()
  end if

  ! Send the input parameters to all MPI tasks
  call mpibcast(filename_in, root=0, comm=bigdft_mpi%mpi_comm)
  call mpibcast(filename_out, root=0, comm=bigdft_mpi%mpi_comm)
  call mpibcast(exp_power, root=0, comm=bigdft_mpi%mpi_comm)
  call mpibcast(norder_polynomial, root=0, comm=bigdft_mpi%mpi_comm)

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

  call timing(bigdft_mpi%mpi_comm,'CALC_LINEAR','PR')

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

  call timing(bigdft_mpi%mpi_comm,'CHECK_LINEAR','PR')

  ! Do the operation using exact LAPACK and the dense matrices
  if (bigdft_mpi%iproc==0) then
      call yaml_comment('Do the same calculation using dense LAPACK',hfill='-')
  end if
  call operation_using_dense_lapack(bigdft_mpi%iproc, bigdft_mpi%nproc, smat_in, mat_in)

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
  call mpifinalize()
  !call bigdft_finalize(ierr)

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

      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_open('Check the deviation from unity of the operation S^x*S^-x')
          call yaml_map('max_error',max_error,fmt='(es10.3)')
          call yaml_map('sum_error',sum_error/real(smat%nvctr,kind=8),fmt='(es10.3)')
          call yaml_mapping_close()
      end if

    end subroutine check_deviation_from_unity


    subroutine check_deviation_from_unity_dense(n, mat)
      implicit none

      ! Calling arguments
      integer,intent(in) :: n
      real(kind=8),dimension(n,n),intent(in) :: mat

      ! Local variables
      integer :: i, j
      real(kind=8) :: sum_error, max_error, error

      sum_error = 0.d0
      max_error = 0.d0
      do i=1,n
          do j=1,n
              if (j==i) then
                  error = abs(mat(j,i)-1.d0)
              else
                  error = abs(mat(j,i))
              end if
              sum_error = sum_error + error
              max_error = max(max_error,error)
          end do
      end do


      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_open('Check the deviation from unity of the operation S^x*S^-x')
          call yaml_map('max_error',max_error,fmt='(es10.3)')
          call yaml_map('sum_error',sum_error/real(n**2,kind=8),fmt='(es10.3)')
          call yaml_mapping_close()
      end if

    end subroutine check_deviation_from_unity_dense



    !> Calculate matrix**power, using the dense matrix and exact LAPACK operations
    subroutine matrix_power_dense(iproc, nproc, blocksize, n, mat_in, ex, mat_out)
      use module_base
      use parallel_linalg, only: dgemm_parallel, dsyev_parallel
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, blocksize, n
      real(kind=8),dimension(n,n),intent(in) :: mat_in
      real(kind=8),intent(in) :: ex
      real(kind=8),dimension(n,n),intent(out) :: mat_out

      ! Local variables
      integer :: i, j, info
      real(kind=8) :: tt
      real(kind=8),dimension(:,:,:),allocatable :: mat_tmp
      real(kind=8),dimension(:),allocatable :: eval

      call f_routine(id='matrix_power_dense')


      ! Diagonalize the matrix
      mat_tmp = f_malloc((/n,n,2/),id='mat_tmp')
      eval = f_malloc(n,id='mat_tmp')
      call f_memcpy(src=mat_in, dest=mat_tmp)
      call dsyev_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'v', 'l', n, mat_tmp, n, eval, info)
      if (info /= 0) then
          call f_err_throw('wrong infocode, value ='//trim(yaml_toa(info)))
      end if

      ! Multiply a diagonal matrix containing the eigenvalues to the power ex with the diagonalized matrix
      do i=1,n
          tt = eval(i)**ex
          do j=1,n
              mat_tmp(j,i,2) = mat_tmp(j,i,1)*tt
          end do
      end do

      ! Apply the diagonalized overlap matrix to the matrix constructed above
      call dgemm_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'n', 't', n, n, n, 1.d0, mat_tmp(1,1,1), n, &
           mat_tmp(1,1,2), n, 0.d0, mat_out, n)

      call f_free(mat_tmp)
      call f_free(eval)

      call f_release_routine()


    end subroutine matrix_power_dense


    subroutine operation_using_dense_lapack(iproc, nproc, smat_in, mat_in)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix, only: uncompress_matrix
      use parallel_linalg, only: dgemm_parallel
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: smat_in
      type(matrices),intent(in) :: mat_in

      ! Local variables
      integer :: blocksize
      real(kind=8),dimension(:,:),allocatable :: mat_in_dense, mat_out_dense
      real(kind=8),dimension(:,:,:),allocatable :: mat_check_accur_dense

      call f_routine(id='operation_using_dense_lapack')

      blocksize = -100
      mat_in_dense = f_malloc((/smat_in%nfvctr,smat_in%nfvctr/),id='mat_in_dense')
      mat_out_dense = f_malloc((/smat_in%nfvctr,smat_in%nfvctr/),id='mat_out_dense')
      mat_check_accur_dense = f_malloc((/smat_in%nfvctr,smat_in%nfvctr,2/),id='mat_check_accur_dense')
      call uncompress_matrix(bigdft_mpi%iproc, smat_in, mat_in%matrix_compr, mat_in_dense)
      call timing(bigdft_mpi%mpi_comm,'INIT_CUBIC','PR')
      call matrix_power_dense(bigdft_mpi%iproc, bigdft_mpi%nproc, blocksize, smat_in%nfvctr, &
           mat_in_dense, exp_power, mat_out_dense)
      call timing(bigdft_mpi%mpi_comm,'CALC_CUBIC','PR')
      call matrix_power_dense(bigdft_mpi%iproc, bigdft_mpi%nproc, blocksize, smat_in%nfvctr, &
           mat_in_dense, -exp_power, mat_check_accur_dense)
      call dgemm_parallel(bigdft_mpi%iproc, bigdft_mpi%nproc, blocksize, bigdft_mpi%mpi_comm, 'n', 'n', &
           smat_in%nfvctr, smat_in%nfvctr, smat_in%nfvctr, &
           1.d0, mat_out_dense(1,1), smat_in%nfvctr, &
           mat_check_accur_dense(1,1,1), smat_in%nfvctr, 0.d0, mat_check_accur_dense(1,1,2), smat_in%nfvctr)
      call check_deviation_from_unity_dense(smat_in%nfvctr, mat_check_accur_dense(1,1,2))
      call timing(bigdft_mpi%mpi_comm,'CHECK_CUBIC','PR')
      call f_free(mat_check_accur_dense)
      call f_free(mat_in_dense)
      call f_free(mat_out_dense)

      call f_release_routine()

    end subroutine operation_using_dense_lapack

end program driver_single
