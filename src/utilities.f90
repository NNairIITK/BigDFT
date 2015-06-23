!> @file
!!   Program to perform pre-/post-processing.
!! @author
!!   Copyright (C) 2015 BigDFT group (SM)
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Test the input files and estimates the memory occupation versus the number
!! of processors
program utilities

   use module_base
   use yaml_output
   use module_types, only: bigdft_init_errors, bigdft_init_timing_categories
   use module_atoms, only: atoms_data, atoms_data_null, deallocate_atoms_data
   use io, only: read_sparse_matrix
   use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, assignment(=), SPARSE_FULL, &
                                sparsematrix_malloc_ptr, deallocate_sparse_matrix, deallocate_matrices
   use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple
   use postprocessing_linear, only: loewdin_charge_analysis_core
   use bigdft_run, only: bigdft_init
   implicit none
   character(len=*), parameter :: subname='utilities'
   character(len=30) :: tatonam, radical
   character(len=128) :: overlap_file, kernel_file
   logical :: charge_analysis = .false.
   type(atoms_data) :: at
   integer :: istat, i_arg, ierr, nspin, icount, nthread
   integer :: nfvctr_s, nseg_s, nvctr_s, nfvctrp_s, isfvctr_s
   integer :: nfvctr_l, nseg_l, nvctr_l, nfvctrp_l, isfvctr_l
   integer,dimension(:),allocatable :: on_which_atom
   integer,dimension(:),pointer :: keyv_s, keyv_l, on_which_atom_s, on_which_atom_l
   integer,dimension(:,:,:),pointer :: keyg_s, keyg_l
   real(kind=8),dimension(:),pointer :: matrix_compr
   type(matrices) :: ovrlp_mat, kernel_mat
   logical :: mpi_init
   type(sparse_matrix) :: smat_s, smat_l
   !$ integer :: omp_get_max_threads

   call f_lib_initialize()

   ! Initialize MPI
   !call bigdft_mpi_init(ierr)
   call bigdft_init()

   if (bigdft_mpi%iproc==0) then
       call yaml_scalar('',hfill='~')
       call yaml_scalar('BIGDFT PRE-/POST-PROCESSING',hfill='~')
   end if

   if (bigdft_mpi%iproc==0) then
       call yaml_mapping_open('Parallel environment')
       call yaml_map('MPI tasks',bigdft_mpi%nproc)
       nthread = 1
       !$ nthread = omp_get_max_threads()
       call yaml_map('OpenMP threads',nthread)
       call yaml_mapping_close()
   end if

   ! Get arguments
   call get_command_argument(1, value = tatonam, status = istat)

   write(radical, "(A)") "input"
   if(trim(tatonam)=='' .or. istat>0) then
      write(*,'(1x,a)')&
         &   'Usage: ./utilities  [option]'
      write(*,'(1x,a)')&
         &   '[option] can be the following: '
      write(*,'(1x,a)')&
           &   '"charge-analysis"" ' 
      write(*,'(1x,a)')&
           & 'perform a Loewdin charge analysis'

      stop
   else
      i_arg = 1
      loop_getargs: do
         call get_command_argument(i_arg, value = tatonam, status = istat)
         !call getarg(i_arg,tatonam)
         if(trim(tatonam)=='' .or. istat > 0) then
            exit loop_getargs
         else if (trim(tatonam)=='charge-analysis') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            !write(*,'(1x,2a)')&
            !   &   'perform a Loewdin charge analysis'
            charge_analysis = .true.
            exit loop_getargs
         end if
         i_arg = i_arg + 1
      end do loop_getargs
   end if


   if (charge_analysis) then
       if (bigdft_mpi%iproc==0) then
           call yaml_comment('Loewdin charge analysis',hfill='-')
       end if
       
       !call set_astruct_from_file(trim(posinp_file),0,at%astruct,fcomment,energy,fxyz)

       at = atoms_data_null()

       call read_sparse_matrix(trim(overlap_file), nspin, nfvctr_s, nseg_s, nvctr_s, keyv_s, keyg_s, &
            matrix_compr, at%astruct%nat, at%astruct%ntypes, at%nzatom, at%nelpsp, &
            at%astruct%atomnames, at%astruct%iatype, at%astruct%rxyz,  on_which_atom=on_which_atom_s)
       at%refcnt=f_ref_new('atoms')
       call distribute_columns_on_processes_simple(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr_s, nfvctrp_s, isfvctr_s)
       call bigdft_to_sparsebigdft(bigdft_mpi%iproc, bigdft_mpi%nproc, nspin, nfvctr_s, nfvctrp_s, isfvctr_s, &
            on_which_atom_s, nvctr_s, nseg_s, keyg_s, smat_s)
       call f_free_ptr(keyv_s)
       call f_free_ptr(keyg_s)
       call f_free_ptr(on_which_atom_s)
       ovrlp_mat = matrices_null()
       ovrlp_mat%matrix_compr = sparsematrix_malloc_ptr(smat_s, iaction=SPARSE_FULL, id='ovrlp%matrix_compr')
       call vcopy(smat_s%nvctr*smat_s%nspin, matrix_compr(1), 1, ovrlp_mat%matrix_compr(1), 1)
       call f_free_ptr(matrix_compr)

       call read_sparse_matrix(trim(kernel_file), nspin, nfvctr_l, nseg_l, nvctr_l, keyv_l, keyg_l, &
            matrix_compr, on_which_atom=on_which_atom_l)
       call distribute_columns_on_processes_simple(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr_l, nfvctrp_l, isfvctr_l)
       call bigdft_to_sparsebigdft(bigdft_mpi%iproc, bigdft_mpi%nproc, nspin, nfvctr_l, nfvctrp_l, isfvctr_l, &
            on_which_atom_l, nvctr_l, nseg_l, keyg_l, smat_l)
       call f_free_ptr(keyv_l)
       call f_free_ptr(keyg_l)
       call f_free_ptr(on_which_atom_l)
       kernel_mat = matrices_null()
       kernel_mat%matrix_compr = sparsematrix_malloc_ptr(smat_l, iaction=SPARSE_FULL, id='kernel_mat%matrix_compr')
       call vcopy(smat_l%nvctr*smat_l%nspin, matrix_compr(1), 1, kernel_mat%matrix_compr(1), 1)
       call f_free_ptr(matrix_compr)

       call loewdin_charge_analysis_core(bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s%nfvctr, smat_s%nfvctrp, smat_s%isfvctr, &
            smat_s%nfvctr_par, smat_s%isfvctr_par, meth_overlap=1020, &
            smats=smat_s, smatl=smat_l, atoms=at, kernel=kernel_mat, ovrlp=ovrlp_mat)

       call deallocate_atoms_data(at)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_l)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)

       if (bigdft_mpi%iproc==0) then
           call yaml_comment('done',hfill='-')
       end if

   end if

   call bigdft_finalize(ierr)
   call f_lib_finalize()


end program utilities
