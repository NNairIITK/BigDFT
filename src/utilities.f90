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
   use postprocessing_linear, only: CHARGE_ANALYSIS_LOEWDIN, CHARGE_ANALYSIS_MULLIKEN, &
                                    loewdin_charge_analysis_core
   use bigdft_run, only: bigdft_init
   implicit none
   external :: gather_timings
   character(len=*), parameter :: subname='utilities'
   character(len=30) :: tatonam, radical
   character(len=128) :: method_name, overlap_file, kernel_file
   logical :: charge_analysis = .false.
   type(atoms_data) :: at
   integer :: istat, i_arg, ierr, nspin, icount, nthread, method
   integer :: nfvctr_s, nseg_s, nvctr_s, nfvctrp_s, isfvctr_s
   integer :: nfvctr_l, nseg_l, nvctr_l, nfvctrp_l, isfvctr_l
   integer,dimension(:),allocatable :: on_which_atom
   integer,dimension(:),pointer :: keyv_s, keyv_l, on_which_atom_s, on_which_atom_l
   integer,dimension(:,:,:),pointer :: keyg_s, keyg_l
   real(kind=8),dimension(:),pointer :: matrix_compr
   type(matrices) :: ovrlp_mat, kernel_mat
   logical :: mpi_init
   type(sparse_matrix) :: smat_s, smat_l
   type(dictionary), pointer :: dict_timing_info
   !$ integer :: omp_get_max_threads

   call f_lib_initialize()

   ! Initialize MPI
   !call bigdft_mpi_init(ierr)
   call bigdft_init()

    if (bigdft_mpi%iproc==0) then
        call yaml_new_document()
    end if


   !Time initialization
   call f_timing_reset(filename='time.yaml',master=(bigdft_mpi%iproc==0),verbose_mode=(.true..and.bigdft_mpi%nproc>1))


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
         &   'Usage: ./utilities -a [option]'
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
            call get_command_argument(i_arg, value = method_name)
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

       ! Determine the method
       select case(trim(method_name))
       case ('loewdin','LOEWDIN')
           method = CHARGE_ANALYSIS_LOEWDIN
       case ('mulliken','MULLIKEN')
           method = CHARGE_ANALYSIS_MULLIKEN
       case default
           call f_err_throw('Unknown Method for the charge analysis',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       end select

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
       call timing(bigdft_mpi%mpi_comm,'INIT','PR')

       call loewdin_charge_analysis_core(method, bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s%nfvctr, smat_s%nfvctrp, smat_s%isfvctr, &
            smat_s%nfvctr_par, smat_s%isfvctr_par, meth_overlap=1020, blocksize=-8, &
            smats=smat_s, smatl=smat_l, atoms=at, kernel=kernel_mat, ovrlp=ovrlp_mat)
       call timing(bigdft_mpi%mpi_comm,'CALC','PR')

       call deallocate_atoms_data(at)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_l)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)

       if (bigdft_mpi%iproc==0) then
           call yaml_comment('done',hfill='-')
       end if

   end if

   call build_dict_info(dict_timing_info)
   call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm,nproc=bigdft_mpi%nproc,&
        gather_routine=gather_timings,dict_info=dict_timing_info)
   call dict_free(dict_timing_info)

   if (bigdft_mpi%iproc==0) then
       call yaml_release_document()
   end if


   call bigdft_finalize(ierr)
   call f_lib_finalize()


  !SM: This routine should go to a module
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


end program utilities
