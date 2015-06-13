program driver_singlerun
  use bigdft_run
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               matrices_null, deallocate_sparse_matrix, deallocate_matrices, &
                               assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL, &
                               sparsematrix_malloc, DENSE_FULL
  use sparsematrix_init, only: read_ccs_format, ccs_to_sparsebigdft, ccs_values_to_bigdft, &
                               read_bigdft_format, bigdft_to_sparsebigdft, distribute_columns_on_processes_simple
  use sparsematrix, only: write_matrix_compressed, check_symmetry, &
                          write_sparsematrix_CCS, write_sparsematrix, &
                          uncompress_matrix
  use matrix_operations, only: overlapPowerGeneral
  use io, only: read_sparse_matrix, write_sparse_matrix
  use yaml_output
  implicit none

  external :: gather_timings

  ! Variables
  integer :: iproc, nproc, ncol, nnonzero, nseg, ncolp, iscol, ierr, nspin, nat, ntypes, lwork, info, ispin, i, j
  character(len=*),parameter :: filename1='matrix1.dat', filename2='matrix2.dat'
  integer,dimension(:),pointer :: col_ptr, row_ind, keyv, nzatom, nelpsp, iatype, on_which_atom
  character(len=20),dimension(:),pointer :: atomnames
  integer,dimension(:,:,:),pointer :: keyg
  real(kind=8),dimension(:,:,:),allocatable :: tempmat, matA_full, matB_full, matC_full
  real(kind=8),dimension(:),pointer :: val
  real(kind=8),dimension(:),allocatable :: work, eval
  real(kind=8),dimension(:,:),pointer :: rxyz
  type(sparse_matrix) :: smatA, smatB
  type(matrices) :: matA
  type(matrices),dimension(1) :: matB
  real(kind=8) :: max_error, mean_error, max_dev, mean_dev, tt
  logical :: symmetric
  real(kind=8) :: time_start, time_end
  type(dictionary), pointer :: dict_timing_info

  ! Initialize
  call f_lib_initialize()
  call bigdft_init()!mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  iproc=bigdft_mpi%iproc!mpi_info(1)
  nproc=bigdft_mpi%nproc!mpi_info(2)

  call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=.true. .and. nproc>1)

  ! Nullify all pointers
  nullify(col_ptr)
  nullify(row_ind)
  nullify(keyv)
  nullify(keyg)
  nullify(val)

  
  ! Read in a file in the CCS format
  !call read_ccs_format(filename, ncol, nnonzero, col_ptr, row_ind, val)

  ! Read in a file in the BigDFT format
  !call read_bigdft_format(filename, ncol, nnonzero, nseg, keyv, keyg, val, on_which_atom=on_which_atom)
  call read_sparse_matrix(filename1, nspin, ncol, nseg, nnonzero, keyv, keyg, val, &
       nat=nat, ntypes=ntypes, nzatom=nzatom, nelpsp=nelpsp, &
       atomnames=atomnames, iatype=iatype, rxyz=rxyz, on_which_atom=on_which_atom)

  ! Create the corresponding BigDFT sparsity pattern
  !call ccs_to_sparsebigdft(iproc, nproc, ncol, ncol, 0, nnonzero, row_ind, col_ptr, smat)
  call distribute_columns_on_processes_simple(iproc, nproc, ncol, ncolp, iscol)
  call bigdft_to_sparsebigdft(iproc, nproc, nspin, ncol, ncolp, iscol, on_which_atom, nnonzero, nseg, keyg, smatA)


  ! Check the symmetry
  symmetric = check_symmetry(ncol, smatA)
  if (.not.symmetric) stop 'ERROR not symmetric'

  matA = matrices_null()
  ! Assign the values
  matA%matrix_compr = sparsematrix_malloc_ptr(smatA, iaction=SPARSE_FULL, id='matA%matrix_compr')
  !call ccs_values_to_bigdft(ncol, nnonzero, row_ind, col_ptr, smat, val, matA)
  matA%matrix_compr = val

  call f_free_ptr(keyv)
  call f_free_ptr(keyg)
  call f_free_ptr(on_which_atom)
  call f_free_ptr(val)

  call read_sparse_matrix(filename2, nspin, ncol, nseg, nnonzero, keyv, keyg, val, on_which_atom=on_which_atom)
       !nat=nat, ntypes=ntypes, nzatom=nzatom, nelpsp=nelpsp, &
       !atomnames=atomnames, iatype=iatype, rxyz=rxyz, on_which_atom=on_which_atom)
  call distribute_columns_on_processes_simple(iproc, nproc, ncol, ncolp, iscol)
  call bigdft_to_sparsebigdft(iproc, nproc, nspin, ncol, ncolp, iscol, on_which_atom, nnonzero, nseg, keyg, smatB)

  matB(1) = matrices_null()


  ! Check the symmetry
  symmetric = check_symmetry(ncol, smatB)
  if (.not.symmetric) stop 'ERROR not symmetric'
  
  ! Write the original matrix
  !if (iproc==0) call write_matrix_compressed('Original matrix', smat, matA)
  !if (iproc==0) call write_sparsematrix_CCS('original_css.dat', smatA, matA)
  !if (iproc==0) call write_sparsematrix('original_bigdft.dat', smatA, matA)


  tempmat = sparsematrix_malloc(smatA, iaction=DENSE_FULL, id='tempmat')
  call uncompress_matrix(iproc, smatA, matA%matrix_compr, tempmat)
  eval = f_malloc(smatA%nfvctr,id='eval')
  lwork=100*smatA%nfvctr
  work = f_malloc(lwork,id='work')
  do ispin=1,smatA%nspin
      max_dev = 0.d0
      mean_dev = 0.d0
      do i=1,smatA%nfvctr
          do j=1,smatA%nfvctr
              if (i==j) then
                  tt = abs(tempmat(j,i,ispin)-1.d0)
              else
                  tt = abs(tempmat(j,i,ispin))
              end if
              max_dev = max(max_dev,tt)
              mean_dev = mean_dev + tt
          end do
      end do
      mean_dev = mean_dev/real(smatA%nvctr,kind=8)
      if (iproc==0) call yaml_map('max/mean dev from unity',(/max_dev,mean_dev/),fmt='(es16.6)')
      call dsyev('n','l', smatA%nfvctr, tempmat(1,1,ispin), smatA%nfvctr, eval, work, lwork, info)
      if (iproc==0) call yaml_map('eval max/min',(/eval(1),eval(smatA%nfvctr)/),fmt='(es16.6)')
  end do
  call f_free(eval)
  call f_free(work)
  call f_free(tempmat)



  call timing(bigdft_mpi%mpi_comm,'INIT','PR')

  ! Calculate the inverse
  !write(*,*) 'smatA%nvctr, smatB%nvctr', smatA%nvctr, smatB%nvctr
  matB(1)%matrix_compr = sparsematrix_malloc_ptr(smatB, iaction=SPARSE_FULL, id='matB(1)%matrix_compr')
  !matB(2)%matrix_compr = sparsematrix_malloc_ptr(smatB, iaction=SPARSE_FULL, id='matB(2)%matrix_compr')
  !matB(3)%matrix_compr = sparsematrix_malloc_ptr(smatB, iaction=SPARSE_FULL, id='matB(3)%matrix_compr')
  call mpibarrier(bigdft_mpi%mpi_comm)
  time_start = mpi_wtime()
  !call overlapPowerGeneral(iproc, nproc, 1020, 1, (/1/), -8, &
  !     1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
  !     check_accur=.true., max_error=max_error, mean_error=mean_error)
  call overlapPowerGeneral(iproc, nproc, 50, 1, (/1/), -8, &
       imode=1, ovrlp_smat=smatA, inv_ovrlp_smat=smatB, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
       check_accur=.true., max_error=max_error, mean_error=mean_error)
  !!call overlapPowerGeneral(iproc, nproc, 0, 1, (/1/), -8, &
  !!     1, ovrlp_smat=smat, inv_ovrlp_smat=smat, ovrlp_mat=matA, inv_ovrlp_mat=matB, &
  !!     check_accur=.true., max_error=max_error, mean_error=mean_error)
  call mpibarrier(bigdft_mpi%mpi_comm)
  time_end = mpi_wtime()
  call timing(bigdft_mpi%mpi_comm,'CALC','PR')
  if (iproc==0) write(*,*) 'walltime',time_end-time_start


  ! Write the original and the final matrix in dense format ####################
  if (iproc==0) then
      matA_full = sparsematrix_malloc(smatA, iaction=DENSE_FULL, id='matA_full')
      matB_full = sparsematrix_malloc(smatB, iaction=DENSE_FULL, id='matB_full')
      matC_full = sparsematrix_malloc(smatB, iaction=DENSE_FULL, id='matC_full')
      call uncompress_matrix(iproc, smatA, matA%matrix_compr, matA_full)
      !do ispin=1,smatA%nspin
      !   do i=1,smatA%nfvctr
      !      do j=1,smatA%nfvctr
      !         write(200,'(2(i6,1x),es19.12)') i,j,matA_full(i,j,ispin)
      !      end do
      !   end do
      !end do
      call uncompress_matrix(iproc, smatB, matB(1)%matrix_compr, matB_full)
      !do ispin=1,smatB%nspin
      !   do i=1,smatB%nfvctr
      !      do j=1,smatB%nfvctr
      !         write(300,'(2(i6,1x),es19.12)') i,j,matB_full(i,j,ispin)
      !      end do
      !   end do
      !end do

      call dgemm('n', 'n', smatB%nfvctr, smatB%nfvctr, smatB%nfvctr, 1.d0, &
           matA_full, smatB%nfvctr, matB_full, smatB%nfvctr, 0.d0, matC_full, smatB%nfvctr)

      max_dev = 0.d0
      mean_dev = 0.d0
      do ispin=1,smatB%nspin
          do i=1,smatB%nfvctr
              do j=1,smatB%nfvctr
                  !write(400,*) i, j, matC_full(j,i,ispin)
                  if (i==j) then
                      max_dev = max(max_dev,(matC_full(j,i,ispin)-1.d0)**2)
                      mean_dev = mean_dev + (matC_full(j,i,ispin)-1.d0)**2
                  else
                      max_dev = max(max_dev,matC_full(j,i,ispin)**2)
                      mean_dev = mean_dev + matC_full(j,i,ispin)**2
                  end if
              end do
          end do
      end do
      mean_dev = mean_dev/real(smatB%nfvctr,kind=8)**2
      mean_dev=sqrt(mean_dev)
      max_dev=sqrt(max_dev)

      call yaml_map('max / mean dev for full matrix',(/max_dev,mean_dev/),fmt='(es10.3)')
      call f_free(matA_full)
      call f_free(matB_full)
      call f_free(matC_full)

  end if
  ! ############################################################################

  ! Write the results
  !if (iproc==0) call write_matrix_compressed('Result', smat, matB(1))

  !if (iproc==0) call write_sparsematrix_CCS('result_css.dat', smatB, matB(1))
  !if (iproc==0) call write_sparsematrix('result_bigdft.dat', smatB, matB(1))
  ! write_sparse_matrix must be called by all MPI tasks
  !call write_sparse_matrix(nat, ntypes, iatype, rxyz, nzatom, nelpsp, atomnames, smatB, matB(1), 'result_bigdft_new.dat')

  ! Deallocations
  call deallocate_sparse_matrix(smatA)
  call deallocate_sparse_matrix(smatB)
  call deallocate_matrices(matA)
  call deallocate_matrices(matB(1))
  call f_free_ptr(col_ptr)
  call f_free_ptr(row_ind)
  call f_free_ptr(keyv)
  call f_free_ptr(keyg)
  call f_free_ptr(val)
  call f_free_ptr(on_which_atom)
  call f_free_ptr(nzatom)
  call f_free_ptr(nelpsp)
  call f_free_ptr(iatype)
  call f_free_str_ptr(len(atomnames),atomnames)
  call f_free_ptr(rxyz)

  call timing(bigdft_mpi%mpi_comm,'FINISH','PR')

  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm, nproc=bigdft_mpi%nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  if (iproc == 0) then
     !call yaml_comment('Timing for root process',hfill='-')
     !call yaml_mapping_open('Timings for root process')
     !call yaml_map('CPU time (s)',tcpu1-tcpu0,fmt='(f12.2)')
     !call yaml_map('Elapsed time (s)',tel,fmt='(f12.2)')
     call yaml_mapping_close()
     call yaml_flush_document()
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

end program driver_singlerun


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



