!> @file
!!   Test of the sparse matrix multiplication
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


program smatmul
  use wrapper_MPI
  use wrapper_linalg
  use futile
  !use module_base
  !use yaml_parse, only: yaml_cl_parse, yaml_cl_parse_null, yaml_cl_parse_free, yaml_cl_parse_cmd_line
  !use dictionaries
  !use yaml_output
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               matrices_null, deallocate_sparse_matrix, deallocate_matrices, &
                               assignment(=), sparsematrix_malloc_ptr, sparsematrix_malloc, SPARSE_FULL, SPARSEMM_SEQ, &
                               SPARSE_MATMUL_SMALL
  use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple, check_symmetry
  use sparsematrix, only: write_matrix_compressed, &
                          write_sparsematrix_CCS, write_sparsematrix, &
                          sparsemm_new, sequential_acces_matrix_fast2, &
                          compress_matrix_distributed_wrapper
  use sparsematrix_io, only: read_sparse_matrix
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft
  implicit none

  external :: gather_timings

  ! Variables
  integer :: iproc, nproc, ncol, nnonzero, nseg, ncolp, iscol, ierr, nspin,comm
  integer :: nfvctr, nvctr, isfvctr, nfvctrp, nit, it, verbosity, nat, ntypes
  !character(len=*),parameter :: filename='matrix.dat'
  character(len=1024) :: filename
  character(len=1) :: geocode
  integer,dimension(:),pointer :: col_ptr, row_ind, keyv, on_which_atom, nzatom, nelpsp, iatype
  integer,dimension(:,:,:),pointer :: keyg
  character(len=20),dimension(:),pointer :: atomnames
  real(kind=8),dimension(:,:),pointer :: rxyz
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
  real(kind=8),dimension(3) :: cell_dim


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

  !call bigdft_init()!mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()
  comm=mpiworld()

  if (iproc==0) then
      call yaml_new_document()
  end if

  call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=.true. .and. nproc>1)

  ! Nullify all pointers
  !!nullify(keyv)
  !!nullify(keyg)
  !!nullify(mat_compr)
  nullify(on_which_atom)

  
  ! Read in a file in the sparse BigDFT format
  !!call read_sparse_matrix(filename, nspin, geocode, cell_dim, nfvctr, nseg, nvctr, keyv, keyg, &
  !!     mat_compr, nat=nat, ntypes=ntypes, nzatom=nzatom, nelpsp=nelpsp, &
  !!     atomnames=atomnames, iatype=iatype, rxyz=rxyz, on_which_atom=on_which_atom)

  !!! Create the corresponding BigDFT sparsity pattern
  !!!call distribute_columns_on_processes_simple(iproc, nproc, nfvctr, nfvctrp, isfvctr)
  !!call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nvctr, nseg, keyg, smat, &
  !!     init_matmul=.true., nspin=nspin, geocode=geocode, cell_dim=cell_dim, on_which_atom=on_which_atom)

  !!matA = matrices_null()

  !!! Assign the values
  !!matA%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='matA%matrix_compr')
  !!matA%matrix_compr = mat_compr

  call sparse_matrix_and_matrices_init_from_file_bigdft('serial', 'serial', filename, iproc, nproc,comm, smat, matA, &
       init_matmul=.true.)!, nat=nat, ntypes=ntypes, nzatom=nzatom, nelpsp=nelpsp, &
       !atomnames=atomnames, iatype=iatype, rxyz=rxyz, on_which_atom=on_which_atom)

  ! Check the symmetry
  symmetric = check_symmetry(smat)
  if (.not.symmetric) stop 'ERROR not symmetric'
  
  ! Write the original matrix
  if (iproc==0) call write_sparsematrix('original_bigdft.dat', smat, matA)

  call f_timing_checkpoint(ctr_name='INIT',mpi_comm=comm,nproc=nproc,&
       gather_routine=gather_timings)

  ! Calculate the inverse
  call mpibarrier(comm)
  time_start = mpi_wtime()

  mat_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='mat_seq')
  vector_in = f_malloc0(smat%smmm%nvctrp,id='vector_in')
  vector_out = f_malloc0(smat%smmm%nvctrp,id='vector_out')
  call sequential_acces_matrix_fast2(smat, matA%matrix_compr, mat_seq)

  call vcopy(smat%smmm%nvctrp, matA%matrix_compr(smat%smmm%isvctr_mm_par(iproc)+1), 1, vector_in(1), 1)
  do it=1,nit
      call sparsemm_new(iproc, smat, mat_seq, vector_in, vector_out)
      call vcopy(smat%smmm%nvctrp, vector_out(1), 1, vector_in(1), 1)
  end do


  call compress_matrix_distributed_wrapper(iproc, nproc, smat, SPARSE_MATMUL_SMALL, &
       vector_out, matA%matrix_compr)
  if (iproc==0 .and. verbosity==2) then
      call write_matrix_compressed('final result', smat, matA)
  end if

  call mpibarrier(comm)
  time_end = mpi_wtime()
  call f_timing_checkpoint(ctr_name='CALC',mpi_comm=comm,nproc=nproc,&
       gather_routine=gather_timings)

  ! Deallocations
  call deallocate_sparse_matrix(smat)
  call deallocate_matrices(matA)
  !call f_free_ptr(keyv)
  !call f_free_ptr(keyg)
  call f_free(mat_seq)
  !call f_free_ptr(mat_compr)
  call f_free(vector_in)
  call f_free(vector_out)
  !call f_free_ptr(on_which_atom)
  !call f_free_ptr(nzatom)
  !call f_free_ptr(nelpsp)
  !call f_free_ptr(iatype)
  !call f_free_str_ptr(int(len(atomnames),kind=4),atomnames)
  !call f_free_ptr(rxyz)

  call f_timing_checkpoint(ctr_name='FINISH',mpi_comm=comm,nproc=nproc,&
       gather_routine=gather_timings)

  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=comm, nproc=nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  if (iproc==0) then
      call yaml_release_document()
  end if
  
  call mpifinalize()

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
      type(dictionary), pointer :: dict_tmp,list
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

      if (nproc>1) then
         call mpihostnames_list(comm,list)
         if (iproc==0) then
            call set(dict_info//'Hostnames',list)
         else
            call dict_free(list)
         end if
      end if

    end subroutine build_dict_info

end program smatmul



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
