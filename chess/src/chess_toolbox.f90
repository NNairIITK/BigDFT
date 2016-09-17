!> @file
!!   Auxiliary program to perform pre-/post-processing
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


program chess_toolbox

   !use module_base
   use yaml_output
   !!use module_types, only: bigdft_init_errors, bigdft_init_timing_categories
   !!use module_atoms, only: atoms_data, atoms_data_null, deallocate_atoms_data
   use sparsematrix_base
   use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple, &
                                write_sparsematrix_info
   use sparsematrix_io, only: read_sparse_matrix, write_sparse_matrix, write_dense_matrix
   use sparsematrix, only: uncompress_matrix, uncompress_matrix_distributed2, diagonalizeHamiltonian2, &
                           transform_sparse_matrix
   use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft, &
                                     sparse_matrix_and_matrices_init_from_file_ccs, &
                                     sparse_matrix_metadata_init_from_file, &
                                     matrices_init_from_file_bigdft, &
                                     ccs_data_from_sparse_matrix, &
                                     ccs_matrix_write, &
                                     matrices_init, &
                                     get_selected_eigenvalues_from_FOE, &
                                     matrix_matrix_multiplication
   !!use postprocessing_linear, only: CHARGE_ANALYSIS_LOEWDIN, CHARGE_ANALYSIS_MULLIKEN
   !!use multipole, only: multipole_analysis_driver_new
   !!use multipole_base, only: lmax
   use sparsematrix_io, only: write_linear_coefficients, read_linear_coefficients
   !!use bigdft_run, only: bigdft_init
   implicit none
   external :: gather_timings
   character(len=*), parameter :: subname='utilities'
   character(len=1) :: geocode
   character(len=3) :: do_ortho
   character(len=30) :: tatonam, radical, colorname, linestart, lineend, cname, methodc
   character(len=128) :: method_name, overlap_file, hamiltonian_file, kernel_file, coeff_file, pdos_file, metadata_file
   character(len=128) :: line, cc, output_pdos, conversion, infile, outfile, iev_min_, iev_max_, fscale_, matrix_basis
   !!character(len=128),dimension(-lmax:lmax,0:lmax) :: multipoles_files
   character(len=128) :: kernel_matmul_file, fragment_file
   logical :: multipole_analysis = .false.
   logical :: solve_eigensystem = .false.
   logical :: calculate_pdos = .false.
   logical :: convert_matrix_format = .false.
   logical :: calculate_selected_eigenvalues = .false.
   logical :: kernel_purity = .false.
   !!type(atoms_data) :: at
   type(sparse_matrix_metadata) :: smmd
   integer :: iproc, nproc
   integer :: istat, i_arg, ierr, nspin, icount, nthread, method, ntypes
   integer :: nfvctr_m, nseg_m, nvctr_m
   integer :: nfvctr_l, nseg_l, nvctr_l
   !integer :: nfvctrp_l, isfvctr_l, nfvctrp_m, isfvctr_m, nfvctrp_s, isfvctr_s
   integer :: iconv, iev, iev_min, iev_max, l, ll, m, jat, jat_start, jtype
   integer,dimension(:),pointer :: on_which_atom
   integer,dimension(:),pointer :: keyv_s, keyv_m, keyv_l, on_which_atom_s, on_which_atom_m, on_which_atom_l
   integer,dimension(:),pointer :: iatype, nzatom, nelpsp
   integer,dimension(:),pointer :: col_ptr, row_ind
   integer,dimension(:,:,:),pointer :: keyg_s, keyg_m, keyg_l
   integer,dimension(:),allocatable :: fragment_atom_id, fragment_supfun_id
   !!logical,dimension(-lmax:lmax) :: file_present
   real(kind=8),dimension(:),pointer :: matrix_compr, eval_ptr
   real(kind=8),dimension(:,:),pointer :: rxyz, coeff_ptr
   real(kind=8),dimension(:),allocatable :: eval, energy_arr, occups
   real(kind=8),dimension(:,:),allocatable :: denskernel, pdos, occup_arr
   real(kind=8),dimension(:,:),allocatable :: kernel_fragment, overlap_fragment, ksk_fragment, tmpmat
   logical,dimension(:,:),allocatable :: calc_array
   logical :: file_exists, found, found_a_fragment, found_icol, found_irow
   type(matrices) :: ovrlp_mat, hamiltonian_mat, kernel_mat, mat, ovrlp_large, KS_large
   type(matrices),dimension(1) :: ovrlp_minus_one_half
   type(matrices),dimension(:,:),allocatable :: multipoles_matrices
   type(sparse_matrix) :: smat_s, smat_m, smat_l, smat
   type(dictionary), pointer :: dict_timing_info
   integer :: iunit, nat, iat, iat_prev, ii, iitype, iorb, itmb, itype, ival, ios, ipdos, ispin
   integer :: jtmb, norbks, npdos, npt, ntmb, jjtmb, nat_frag, nfvctr_frag, i, iiat
   integer :: icol, irow, icol_atom, irow_atom, iseg, iirow, iicol, j, ifrag
   character(len=20),dimension(:),pointer :: atomnames
   character(len=128),dimension(:),allocatable :: pdos_name, fragment_atomnames
   real(kind=8),dimension(3) :: cell_dim
   character(len=2) :: backslash, num
   real(kind=8) :: energy, occup, occup_pdos, total_occup, fscale, factor, maxdiff, meandiff, tt
   type(f_progress_bar) :: bar
   integer,parameter :: ncolors = 12
   ! Presumably well suited colorschemes from colorbrewer2.org
   character(len=20),dimension(ncolors),parameter :: colors=(/'#a6cee3', &
                                                              '#1f78b4', &
                                                              '#b2df8a', &
                                                              '#33a02c', &
                                                              '#fb9a99', &
                                                              '#e31a1c', &
                                                              '#fdbf6f', &
                                                              '#ff7f00', &
                                                              '#cab2d6', &
                                                              '#6a3d9a', &
                                                              '#ffff99', &
                                                              '#b15928'/)
   !$ integer :: omp_get_max_threads

   call f_lib_initialize()

   ! Initialize MPI
   !call bigdft_mpi_init(ierr)
   !!call bigdft_init()
   ! MPI initialization; we have:
   ! iproc is the task ID
   ! nproc is the total number of tasks
   call mpiinit()
   iproc=mpirank()
   nproc=mpisize()

   call f_malloc_set_status(memory_limit=0.e0,iproc=iproc)

   ! Initialize the sparse matrix errors and timings.
   call sparsematrix_init_errors
   call sparsematrix_initialize_timing_categories()


    if (iproc==0) then
        call yaml_new_document()
        !call print_logo()
    end if

   !Time initialization
   call f_timing_reset(filename='time.yaml',master=(iproc==0),verbose_mode=.false.)

   if (iproc==0) then
       call yaml_scalar('',hfill='~')
       call yaml_scalar('CHESS TOOLBOX',hfill='~')
   end if

   if (iproc==0) then
       call yaml_mapping_open('Parallel environment')
       call yaml_map('MPI tasks',nproc)
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
           &   '"multipoles-analysis"" ' 
      write(*,'(1x,a)')&
           & 'perform a charge analysis (Loewdin or Mulliken)'

      stop
   else
      i_arg = 1
      loop_getargs: do
         call get_command_argument(i_arg, value = tatonam, status = istat)
         !call getarg(i_arg,tatonam)
         if(trim(tatonam)=='' .or. istat > 0) then
            exit loop_getargs
         !!else if (trim(tatonam)=='multipole-analysis') then
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = method_name)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = matrix_basis)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = metadata_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = overlap_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = kernel_file)
         !!   i_arg = i_arg + 1
         !!   call get_command_argument(i_arg, value = hamiltonian_file)
         !!   do l=0,lmax
         !!       do m=-l,l
         !!           i_arg = i_arg + 1
         !!           call get_command_argument(i_arg, value = multipoles_files(m,l))
         !!       end do
         !!   end do
         !!   !write(*,'(1x,2a)')&
         !!   !   &   'perform a Loewdin charge analysis'
         !!   multipole_analysis = .true.
         !!   exit loop_getargs
         else if (trim(tatonam)=='solve-eigensystem') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            !write(*,'(1x,2a)')&
            !   &   'perform a Loewdin charge analysis'
            solve_eigensystem = .true.
            exit loop_getargs
         else if (trim(tatonam)=='pdos') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = pdos_file)
            calculate_pdos = .true.
        else if (trim(tatonam)=='convert-matrix-format') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = conversion)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = infile)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = outfile)
            convert_matrix_format = .true.
        else if (trim(tatonam)=='calculate-selected-eigenvalues') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_matmul_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = iev_min_)
            read(iev_min_,fmt=*,iostat=ierr) iev_min
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = iev_max_)
            read(iev_max_,fmt=*,iostat=ierr) iev_max
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fscale_)
            read(fscale_,fmt=*,iostat=ierr) fscale
            calculate_selected_eigenvalues = .true.
        else if (trim(tatonam)=='kernel-purity') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = metadata_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_matmul_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fragment_file)
            kernel_purity = .true.
         end if
         i_arg = i_arg + 1
      end do loop_getargs
   end if


   !!if (multipole_analysis) then
   !!    if (iproc==0) then
   !!        call yaml_comment('Multipole analysis',hfill='-')
   !!    end if

   !!    ! Determine the method
   !!    select case(trim(method_name))
   !!    case ('loewdin','LOEWDIN')
   !!        method = CHARGE_ANALYSIS_LOEWDIN
   !!    case ('mulliken','MULLIKEN')
   !!        method = CHARGE_ANALYSIS_MULLIKEN
   !!    !!case ('projector','PROJECTOR')
   !!    !!    method = CHARGE_ANALYSIS_PROJECTOR
   !!    case default
   !!        call f_err_throw('Unknown Method for the multipole analysis',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
   !!    end select

   !!    call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
   !!    if (iproc==0) then
   !!        call yaml_mapping_open('Atomic System Properties')
   !!        call yaml_map('Types of atoms',smmd%atomnames)
   !!        call yaml_mapping_close()
   !!    end if

   !!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
   !!         iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
   !!         init_matmul=.true.)

   !!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(kernel_file), &
   !!         iproc, nproc, mpiworld(), smat_l, kernel_mat, &
   !!         init_matmul=.true.)

   !!    call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
   !!         iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
   !!         init_matmul=.true.)

   !!    ! Check which multipole matrices are present
   !!    do l=0,lmax
   !!        file_present(:) = .false.
   !!        do m=-l,l
   !!            inquire(file=trim(multipoles_files(m,l)), exist=file_exists)
   !!            if (file_exists) then
   !!                file_present(m) = .true.
   !!            end if
   !!        end do
   !!        if (any(file_present(-l:l))) then
   !!            if (.not.all(file_present(-l:l))) then
   !!                call f_err_throw('for a given shell all matrices must be present', &
   !!                     err_name='BIGDFT_RUNTIME_ERROR')
   !!            end if
   !!            ll = l
   !!        end if
   !!    end do

   !!    allocate(multipoles_matrices(-ll:ll,0:ll))
   !!    do l=0,ll
   !!        do m=-l,l
   !!            call matrices_init_from_file_bigdft(trim(multipoles_files(m,l)), &
   !!                 iproc, nproc, mpiworld(), smat_s, multipoles_matrices(m,l))
   !!        end do
   !!    end do

   !!    call timing(mpiworld(),'INIT','PR')

   !!    select case(method)
   !!    case (CHARGE_ANALYSIS_MULLIKEN)
   !!        methodc='loewdin'
   !!        do_ortho='no'
   !!    case (CHARGE_ANALYSIS_LOEWDIN)
   !!        methodc='loewdin'
   !!        do_ortho='yes'
   !!    !!case (CHARGE_ANALYSIS_PROJECTOR)
   !!    !!    methodc='projector'
   !!    !!    do_ortho='no'
   !!    case default
   !!        call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
   !!    end select

   !!    ! Determine whether the multipoles matrices to be used are calculated analytically or on the grid.
   !!    ! For the analytic case, only the overlap matrix is possible.
   !!    select case(trim(matrix_basis))
   !!    case ('wavelet','WAVELET')
   !!        call multipole_analysis_driver_new(iproc, nproc, 0, 11, &
   !!             smmd, smat_s, smat_m, smat_l, &
   !!             ovrlp_mat, hamiltonian_mat, kernel_mat, smmd%rxyz, &
   !!             methodc, do_ortho=trim(do_ortho), projectormode='simple', &
   !!             calculate_multipole_matrices=.false., do_check=.false., &
   !!             write_multipole_matrices=.false., &
   !!             multipole_matrix_in=(/(/ovrlp_mat/)/))
   !!    case ('realspace','REALSPACE')
   !!        call multipole_analysis_driver_new(iproc, nproc, ll, 11, &
   !!             smmd, smat_s, smat_m, smat_l, &
   !!             ovrlp_mat, hamiltonian_mat, kernel_mat, smmd%rxyz, &
   !!             methodc, do_ortho=trim(do_ortho), projectormode='simple', &
   !!             calculate_multipole_matrices=.false., do_check=.false., &
   !!             write_multipole_matrices=.false., &
   !!             multipole_matrix_in=multipoles_matrices)
   !!    case default
   !!        call f_err_throw('wrong value for matrix_basis',err_name='BIGDFT_RUNTIME_ERROR')
   !!    end select

   !!    call timing(mpiworld(),'CALC','PR')

   !!    call deallocate_sparse_matrix_metadata(smmd)
   !!    call deallocate_sparse_matrix(smat_s)
   !!    call deallocate_sparse_matrix(smat_l)
   !!    call deallocate_matrices(ovrlp_mat)
   !!    call deallocate_matrices(kernel_mat)
   !!    call deallocate_sparse_matrix(smat_m)
   !!    call deallocate_matrices(hamiltonian_mat)
   !!    do l=0,ll
   !!        do m=-l,l
   !!            call deallocate_matrices(multipoles_matrices(m,l))
   !!        end do
   !!    end do
   !!    deallocate(multipoles_matrices)

   !!    if (iproc==0) then
   !!        call yaml_comment('done',hfill='-')
   !!    end if

   !!end if


   if (solve_eigensystem) then

       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
            iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
            init_matmul=.false.)!, nat=nat, rxyz=rxyz, iatype=iatype, ntypes=ntypes, &
            !nzatom=nzatom, nelpsp=nelpsp, atomnames=atomnames)
       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
            iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
            init_matmul=.false.)

       ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='ovrlp_mat%matrix')
       call uncompress_matrix(iproc, nproc, &
            smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
       hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='hamiltonian_mat%matrix')
       call uncompress_matrix(iproc, nproc, &
            smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
       eval = f_malloc(smat_s%nfvctr,id='eval')

       if (iproc==0) then
           call yaml_comment('Diagonalizing the matrix',hfill='~')
       end if
       call diagonalizeHamiltonian2(iproc, smat_s%nfvctr, hamiltonian_mat%matrix, ovrlp_mat%matrix, eval)
       if (iproc==0) then
           call yaml_comment('Matrix succesfully diagonalized',hfill='~')
       end if
       iunit=99
       call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
       !call writeLinearCoefficients(iunit, .true., nat, rxyz, smat_s%nfvctr, smat_s%nfvctr, &
       !     smat_s%nfvctr, hamiltonian_mat%matrix, eval)
       call write_linear_coefficients(iproc, 0, trim(coeff_file), 2, smmd%nat, smmd%rxyz, &
            smmd%iatype, smmd%ntypes, smmd%nzatom, &
            smmd%nelpsp, smmd%atomnames, smat_s%nfvctr, &
            smat_s%nfvctr, smat_s%nspin, hamiltonian_mat%matrix, eval)
       call f_close(iunit)

       call f_free(eval)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(hamiltonian_mat)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_m)
       call deallocate_sparse_matrix_metadata(smmd)
       !call f_free_ptr(rxyz)
       !call f_free_ptr(iatype)
       !call f_free_ptr(nzatom)
       !call f_free_ptr(nelpsp)
       !call f_free_str_ptr(len(atomnames),atomnames)
   end if


   if (calculate_pdos) then
       iunit = 99
       if (iproc==0) call yaml_comment('Reading from file '//trim(coeff_file),hfill='~')
       call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
       call read_linear_coefficients(trim(coeff_file), nspin, ntmb, norbks, coeff_ptr, &
            eval=eval_ptr)
       call f_close(iunit)

       if (iproc==0) call yaml_comment('Reading from file '//trim(overlap_file),hfill='~')
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
            iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
            init_matmul=.false.)!, iatype=iatype, ntypes=ntypes, atomnames=atomnames, &
            !on_which_atom=on_which_atom)
       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       if (ntmb/=smat_s%nfvctr) call f_err_throw('ntmb/=smat_s%nfvctr')
       ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_PARALLEL, id='ovrlp_mat%matrix')
       !call uncompress_matrix(iproc, smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
       call uncompress_matrix_distributed2(iproc, smat_s, DENSE_PARALLEL, &
            ovrlp_mat%matrix_compr, ovrlp_mat%matrix(1:,1:,1))

       if (iproc==0) call yaml_comment('Reading from file '//trim(hamiltonian_file),hfill='~')
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
            iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
            init_matmul=.false.)
       hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_PARALLEL, id='hamiltonian_mat%matrix')
       !call uncompress_matrix(iproc, smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
       call uncompress_matrix_distributed2(iproc, smat_m, DENSE_PARALLEL, &
            hamiltonian_mat%matrix_compr, hamiltonian_mat%matrix(1:,1:,1))

       iunit=99
       call f_open_file(iunit, file=pdos_file, binary=.false.)
       read(iunit,*) npdos

       calc_array = f_malloc((/ntmb,npdos/),id='calc_array')
       pdos_name = f_malloc_str(len(pdos_name),npdos,id='pdos_name')


       do ipdos=1,npdos
           do itmb=1,ntmb
               calc_array(itmb,ipdos) = .false.
           end do
       end do
       ipdos = 0
       !npdos_loop: do !ipdos=1,npdos
       if (iproc==0) then
           call yaml_sequence_open('Atoms and support functions to be taken into account for each partial density of states')
       end if
           do 
               !read(iunit01,*,iostat=ios) cc, ival
               read(iunit,'(a128)',iostat=ios) line
               if (ios/=0) exit
               !write(*,*) 'line',line
               read(line,*,iostat=ios) cc, cname
               if (cc=='#') then
                   ipdos = ipdos + 1
                   pdos_name(ipdos) = trim(cname)
                   if (iproc==0) then
                       if (ipdos>1) then
                           call yaml_mapping_close()
                       end if
                       call yaml_sequence(advance='no')
                       call yaml_mapping_open(trim(pdos_name(ipdos)))
                   end if
                   cycle 
               end if
               read(line,*,iostat=ios) cc, ival
               if (iproc==0) then
                   call yaml_map(trim(cc),ival)
               end if
               found = .false.
               do itype=1,smmd%ntypes
                   if (trim(smmd%atomnames(itype))==trim(cc)) then
                       iitype = itype
                       found = .true.
                       exit
                   end if
               end do
               if (.not.found) then
                   call f_err_throw("Atom type '"//trim(cc)//"' is unknown", &
                        err_name='SPARSEMATRIX_RUNTIME_ERROR')
               end if
               iat_prev = -1
               do itmb=1,ntmb
                   iat = smmd%on_which_atom(itmb)
                   if (iat/=iat_prev) then
                       ii = 0
                   end if
                   iat_prev = iat
                   itype = smmd%iatype(iat)
                   ii = ii + 1
                   if (itype==iitype .and. ii==ival) then
                       if (calc_array(itmb,ipdos)) then
                           call f_err_throw('calc_array(itmb) must not be true here', &
                                err_name='SPARSEMATRIX_RUNTIME_ERROR')
                       end if
                       calc_array(itmb,ipdos) = .true.
                   end if
               end do
           end do
       if (iproc==0) then
           call yaml_mapping_close()
           call yaml_sequence_close()
       end if

       !end do npdos_loop
       call f_close(iunit)


       energy_arr = f_malloc0(norbks,id='energy_arr')
       occup_arr = f_malloc0((/npdos,norbks/),id='occup_arr')
       occups = f_malloc0(npdos,id='occups')

       denskernel = f_malloc((/ntmb,smat_s%nfvctrp/),id='denskernel')
       pdos = f_malloc0((/npt,npdos/),id='pdos')
       if (iproc==0) then
           call yaml_comment('PDoS calculation',hfill='~')
           call yaml_mapping_open('Calculating PDoS')
       end if
       ! Calculate a partial kernel for each KS orbital
       !do ipdos=1,npdos
           !if (iproc==0) call yaml_map('PDoS number',ipdos)
           !if (iproc==0) call yaml_map('start, increment',(/ipdos,npdos/))
           !write(num,fmt='(i2.2)') ipdos
           if (iproc==0) bar=f_progress_bar_new(nstep=norbks)
           do iorb=1,norbks
               if (iproc==0) call yaml_map('orbital being processed',iorb)
               call gemm('n', 't', ntmb, smat_s%nfvctrp, 1, 1.d0, coeff_ptr(1,iorb), ntmb, &
                    coeff_ptr(smat_s%isfvctr+1,iorb), ntmb, 0.d0, denskernel(1,1), ntmb)
               energy = 0.d0
               call f_zero(occups)
               do ispin=1,nspin
                   !$omp parallel default(none) &
                   !$omp shared(ispin,smat_s,ntmb,denskernel,hamiltonian_mat,energy) &
                   !$omp private(itmb,jtmb,jjtmb)
                   !$omp do reduction(+:energy)
                   do jtmb=1,smat_s%nfvctrp
                       jjtmb = smat_s%isfvctr + jtmb
                       do itmb=1,ntmb
                           energy = energy + denskernel(itmb,jtmb)*hamiltonian_mat%matrix(itmb,jtmb,ispin)
                       end do
                   end do
                   !$omp end do
                   !$omp end parallel
                   energy_arr(iorb) = energy
               end do
               !call mpiallred(energy, 1, mpi_sum, comm=mpiworld())
               do ispin=1,nspin
                   !$omp parallel default(none) &
                   !$omp shared(ispin,smat_s,ntmb,denskernel,ovrlp_mat,calc_array,npdos,occups) &
                   !$omp private(itmb,jtmb,jjtmb,occup,ipdos)
                   !$omp do reduction(+:occups)
                   do jtmb=1,smat_s%nfvctrp
                       jjtmb = smat_s%isfvctr + jtmb
                       !if (calc_array(jjtmb,ipdos)) then
                           do itmb=1,ntmb!ipdos,ntmb,npdos
                               !if (calc_array(itmb,ipdos)) then
                                   occup = denskernel(itmb,jtmb)*ovrlp_mat%matrix(itmb,jtmb,ispin)
                                   do ipdos=1,npdos
                                       if (calc_array(itmb,ipdos) .and. calc_array(jjtmb,ipdos)) then
                                           occups(ipdos) = occups(ipdos) + occup
                                       end if
                                   end do
                                   !occup = occup + denskernel(itmb,jtmb)*ovrlp_mat%matrix(itmb,jtmb,ispin)
                               !end if
                           end do
                       !end if
                   end do
                   !$omp end do
                   !$omp end parallel
                   do ipdos=1,npdos
                       occup_arr(ipdos,iorb) = occups(ipdos)
                   end do
               end do
               !occup_pdos = occup_pdos + occup
               !!if (iorb<norbks) then
               !!    write(iunit,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2)) + '//trim(backslash)
               !!else
               !!    write(iunit,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2))'
               !!end if
               !!!write(*,'(a,i6,3es16.8)')'iorb, eval(iorb), energy, occup', iorb, eval(iorb), energy, occup
               if (iproc==0) call dump_progress_bar(bar,step=iorb)
           end do
           !total_occup = total_occup + occup_pdos
           !!if (ipdos==1) then
           !!    write(iunit,'(a,i0,a)') "plot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           !!else
           !!    write(iunit,'(a,i0,a)') "replot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           !!end if
           !!if (iproc==0) call yaml_map('sum of PDoS',occup_pdos)
           !!output_pdos='PDoS_'//num//'.dat'
           !!call yaml_map('output file',trim(output_pdos))
           !!call f_open_file(iunit01, file=trim(output_pdos), binary=.false.)
           !!write(iunit01,'(a)') '#             energy                pdos'
           !!do ipt=1,npt
           !!    write(iunit01,'(2es20.12)') energy_bins(1,ipt), pdos(ipt,ipdos)
           !!end do
           !!call f_close(iunit01)
       !end do
       call mpiallred(occup_arr, mpi_sum, comm=mpiworld())
       call mpiallred(energy_arr, mpi_sum, comm=mpiworld())
       if (iproc==0) call yaml_mapping_close()

       if (iproc==0) then
           call yaml_comment('Calculation complete',hfill='=')
           output_pdos='PDoS.gp'
           call yaml_map('output file',trim(output_pdos))
           iunit = 99
           call f_open_file(iunit, file=trim(output_pdos), binary=.false.)
           write(iunit,'(a)') '# plot the DOS as a sum of Gaussians'
           write(iunit,'(a)') 'set samples 1000'
           write(iunit,'(a,2(es12.5,a))') 'set xrange[',eval_ptr(1),':',eval_ptr(ntmb),']'
           write(iunit,'(a)') 'sigma=0.01'
           write(backslash,'(a)') '\ '
           total_occup = 0.d0
           do ipdos=1,npdos
               occup_pdos = 0.d0
               call yaml_map('PDoS number',ipdos)
               call yaml_map('start, increment',(/ipdos,npdos/))
               write(num,fmt='(i2.2)') ipdos
               write(iunit,'(a,i0,a)') 'f',ipdos,'(x) = '//trim(backslash)
               do iorb=1,norbks
                   if (iorb<norbks) then
                       write(iunit,'(2(a,es16.9),a)') '  ',occup_arr(ipdos,iorb),&
                           &'*exp(-(x-',energy_arr(iorb),')**2/(2*sigma**2)) + '//trim(backslash)
                   else
                       write(iunit,'(2(a,es16.9),a)') '  ',occup_arr(ipdos,iorb),&
                           &'*exp(-(x-',energy_arr(iorb),')**2/(2*sigma**2))'
                   end if
                   occup_pdos = occup_pdos + occup_arr(ipdos,iorb)
               end do
               total_occup = total_occup + occup_pdos
               call yaml_map('sum of PDoS',occup_pdos)
           end do
           do ipdos=1,npdos
               if (ipdos<ncolors) then
                   colorname = colors(ipdos)
               else
                   colorname = 'color'
               end if
               if (ipdos<npdos) then
                   lineend = ' ,'//trim(backslash)
               else
                   lineend = ''
               end if
               if (ipdos==1) then
                   linestart = 'plot'
               else
                   linestart = '    '
               end if
               write(iunit,'(a,i0,a)') trim(linestart)//" f",ipdos,"(x) lc rgb '"//trim(colorname)//&
                   &"' lt 1 lw 2 w l title '"//trim(pdos_name(ipdos))//"'"//trim(lineend)
           end do
           call f_close(iunit)
       end if
       if (iproc==0) call yaml_map('sum of total DoS',total_occup)
       call mpibarrier(mpiworld())

       call f_free(pdos)
       call f_free(denskernel)
       call f_free_ptr(eval_ptr)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(hamiltonian_mat)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_m)
       call deallocate_sparse_matrix_metadata(smmd)
       !call f_free_ptr(iatype)
       !call f_free_str_ptr(len(atomnames),atomnames)
       call f_free_str(len(pdos_name),pdos_name)
       call f_free(calc_array)
       !call f_free_ptr(on_which_atom)
       call f_free_ptr(coeff_ptr)
       call f_free(energy_arr)
       call f_free(occup_arr)
       call f_free(occups)
   end if

   if (convert_matrix_format) then
       select case (trim(conversion))
       case ('bigdft_to_ccs')
           iconv = 1
       case ('ccs_to_bigdft')
           iconv = 2
       case ('bigdft_to_dense')
           iconv = 3
       case default
           call f_err_throw("wrong value for conversion; possible are 'bigdft_to_ccs' and 'ccs_to_bigdft'")
       end select

       select case (iconv)
       case (1,3)
           call sparse_matrix_and_matrices_init_from_file_bigdft(trim(infile), &
                iproc, nproc, mpiworld(), &
                smat, mat, init_matmul=.false.)
       case (2)
           call sparse_matrix_and_matrices_init_from_file_ccs(trim(infile), iproc, nproc, &
                mpiworld(), smat, mat, init_matmul=.false.)
       end select
       select case (iconv)
       case (1)
           row_ind = f_malloc_ptr(smat%nvctr,id='row_ind')
           col_ptr = f_malloc_ptr(smat%nfvctr,id='col_ptr')
           call ccs_data_from_sparse_matrix(smat, row_ind, col_ptr)
           if (iproc==0) call ccs_matrix_write(trim(outfile), smat, row_ind, col_ptr, mat)
           call f_free_ptr(row_ind)
           call f_free_ptr(col_ptr)
       case (2)
           call sparse_matrix_and_matrices_init_from_file_ccs(trim(infile), iproc, nproc, &
                mpiworld(), smat, mat, init_matmul=.false.)
           call write_sparse_matrix(iproc, nproc, mpiworld(), &
                smat, mat, trim(outfile))
       case (3)
           call write_dense_matrix(iproc, nproc, mpiworld(), smat, mat, trim(outfile), .false.)

       end select

       call deallocate_sparse_matrix(smat)
       call deallocate_matrices(mat)
   end if


   if (calculate_selected_eigenvalues) then
       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
            iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
            iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(kernel_file), &
            iproc, nproc, mpiworld(), smat_l, kernel_mat, &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file))
       call matrices_init(smat_l, ovrlp_minus_one_half(1))

       !call timing(mpiworld(),'INIT','PR')
       call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)
       if (iev_min<1 .or. iev_min>smat_s%nfvctr .or. iev_max>smat_s%nfvctr .or. iev_max<1) then
           if (iproc==0) then
               call yaml_warning('The required eigenvalues are outside of the possible range, automatic ajustment')
           end if
       end if
       iev_min = max(iev_min,1)
       iev_min = min(iev_min,smat_s%nfvctr)
       iev_max = min(iev_max,smat_s%nfvctr)
       iev_max = max(iev_max,1)
       eval = f_malloc(iev_min.to.iev_max,id='eval')
       if (iproc==0) then
           call yaml_mapping_open('Calculating eigenvalues using FOE')
       end if
       call get_selected_eigenvalues_from_FOE(iproc, nproc, mpiworld(), &
            iev_min, iev_max, smat_s, smat_m, smat_l, ovrlp_mat, hamiltonian_mat, &
            ovrlp_minus_one_half, eval, fscale)

       !!call timing(mpiworld(),'CALC','PR')
       call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)

       if (iproc==0) then
           call yaml_sequence_open('values')
           do iev=iev_min,iev_max
               call yaml_sequence(advance='no')
               call yaml_mapping_open(flow=.true.)
               call yaml_map('ID',iev,fmt='(i6.6)')
               call yaml_map('eval',eval(iev),fmt='(es12.5)')
               call yaml_mapping_close()
           end do
           call yaml_sequence_close()
           call yaml_mapping_close()
       end if


       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_m)
       call deallocate_sparse_matrix(smat_l)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(hamiltonian_mat)
       call deallocate_matrices(kernel_mat)
       call deallocate_matrices(ovrlp_minus_one_half(1))
       call deallocate_sparse_matrix_metadata(smmd)
       call f_free(eval)

       !!call timing(mpiworld(),'LAST','PR')
       call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)

   end if



   if (kernel_purity) then
       call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
            iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
            init_matmul=.false.)
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(kernel_file), &
            iproc, nproc, mpiworld(), smat_l, kernel_mat, &
            init_matmul=.true., filename_mult=trim(kernel_matmul_file))

       if (iproc==0) then
           call yaml_mapping_open('Matrix properties')
           call write_sparsematrix_info(smat_s, 'Overlap matrix')
           call write_sparsematrix_info(smat_l, 'Density kernel')
           call yaml_mapping_close()
       end if

       iunit=99
       call f_open_file(iunit, file=fragment_file, binary=.false.)
       read(iunit,*) nat_frag
       fragment_atomnames = f_malloc_str(len(fragment_atomnames),nat_frag,id='nat_frag')
       do iat=1,nat_frag
           read(iunit,*) fragment_atomnames(iat)
       end do
       call f_close(iunit)


       ovrlp_large = matrices_null()
       ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smat_l, iaction=SPARSE_FULL, id='ovrlp_large%matrix_compr')
       KS_large = matrices_null()
       KS_large%matrix_compr = sparsematrix_malloc_ptr(smat_l, iaction=SPARSE_FULL, id='KS_large%matrix_compr')
       fragment_atom_id = f_malloc(nat_frag,id='fragment_atom_id')
       fragment_supfun_id = f_malloc(smat_l%nfvctr,id='fragment_supfun_id')
       jat_start = 1
       found_a_fragment = .false.
       ifrag = 0
       if (iproc==0) then
           call yaml_map('Fragment composition',fragment_atomnames)
           call yaml_comment('Starting fragment purity analysis',hfill='~')
           call yaml_sequence_open('Purity analysis of fragment')
       end if
       search_fragments: do
           ifrag = ifrag + 1
           fragment_loop: do iat=1,nat_frag
               found = .false.
               search_loop: do jat=jat_start,smmd%nat
                   jtype = smmd%iatype(jat)
                   if (trim(fragment_atomnames(iat))==trim(smmd%atomnames(jtype))) then
                       ! Found an atom beloning to the fragment
                       fragment_atom_id(iat) = jat
                       jat_start = jat + 1
                       found = .true.
                       exit search_loop
                   end if
               end do search_loop
               if (.not.found) then
                   jat_start = smmd%nat + 1
                   exit fragment_loop
               end if
           end do fragment_loop
           if (found) then
               found_a_fragment = .true.
           else
               exit search_fragments
           end if
           
           ! Count how many support functions belong to the fragment
           nfvctr_frag = 0
           do i=1,smat_l%nfvctr
               iiat = smmd%on_which_atom(i)
               do iat=1,nat_frag
                   if (iiat==fragment_atom_id(iat)) then
                       nfvctr_frag = nfvctr_frag + 1
                       exit
                   end if
               end do
               fragment_supfun_id(i) = nfvctr_frag
           end do



           kernel_fragment = f_malloc0((/nfvctr_frag,nfvctr_frag/),id='kernel_fragment')
           overlap_fragment = f_malloc0((/nfvctr_frag,nfvctr_frag/),id='overlap_fragment')
           ksk_fragment = f_malloc((/nfvctr_frag,nfvctr_frag/),id='ksk_fragment')
           tmpmat = f_malloc((/nfvctr_frag,nfvctr_frag/),id='tmpmat')

           if (smat_l%nspin==1) then
               factor = 0.5_mp
           else if (smat_l%nspin==2) then
               factor = 1.0_mp
           end if

           call transform_sparse_matrix(iproc, smat_s, smat_l, SPARSE_FULL, 'small_to_large', &
                                        smat_in=ovrlp_mat%matrix_compr, lmat_out=ovrlp_large%matrix_compr)

           if (iproc==0) then
               call yaml_comment('Fragment number'//trim(yaml_toa(ifrag)),hfill='-')
               call yaml_sequence(advance='no')
               call yaml_map('Atoms ID',fragment_atom_id)
               call yaml_map('Submatrix size',nfvctr_frag)
           end if

           ! Calculate KSK-K for the submatrices
           call extract_fragment_submatrix(smmd, smat_l, ovrlp_large, nat_frag, nfvctr_frag, &
                fragment_atom_id, fragment_supfun_id, overlap_fragment)
           call extract_fragment_submatrix(smmd, smat_l, kernel_mat, nat_frag, nfvctr_frag, &
                fragment_atom_id, fragment_supfun_id, kernel_fragment)
           call gemm('n', 'n', nfvctr_frag, nfvctr_frag, nfvctr_frag, factor, kernel_fragment(1,1), nfvctr_frag, &
                overlap_fragment(1,1), nfvctr_frag, 0.d0, tmpmat(1,1), nfvctr_frag)
           call gemm('n', 'n', nfvctr_frag, nfvctr_frag, nfvctr_frag, 1.d0, tmpmat(1,1), nfvctr_frag, &
                kernel_fragment(1,1), nfvctr_frag, 0.d0, ksk_fragment(1,1), nfvctr_frag)
           maxdiff = 0.0_mp
           meandiff = 0.0_mp
           do i=1,nfvctr_frag
               do j=1,nfvctr_frag
                   tt = abs(kernel_fragment(j,i)-ksk_fragment(j,i))
                   maxdiff = max(maxdiff,tt)
                   meandiff = meandiff + tt
               end do
           end do
           meandiff = meandiff/real(nfvctr_frag,kind=mp)**2
           if (iproc==0) then
               call yaml_mapping_open('analyzing KSK-K')
               call yaml_map('maximal deviation',maxdiff,fmt='(es10.3)')
               call yaml_map('average deviation',meandiff,fmt='(es10.3)')
               call yaml_mapping_close()
           end if

           ! Calculate KS, extract the submatrices and calculate (KS)^2-(KS)
           call matrix_matrix_multiplication(iproc, nproc, smat_l, kernel_mat, ovrlp_large, KS_large)
           call extract_fragment_submatrix(smmd, smat_l, KS_large, nat_frag, nfvctr_frag, &
                fragment_atom_id, fragment_supfun_id, kernel_fragment)
           call gemm('n', 'n', nfvctr_frag, nfvctr_frag, nfvctr_frag, factor, kernel_fragment(1,1), nfvctr_frag, &
                kernel_fragment(1,1), nfvctr_frag, 0.d0, tmpmat(1,1), nfvctr_frag)

           maxdiff = 0.0_mp
           meandiff = 0.0_mp
           do i=1,nfvctr_frag
               do j=1,nfvctr_frag
                   tt = abs(kernel_fragment(j,i)-tmpmat(j,i))
                   maxdiff = max(maxdiff,tt)
                   meandiff = meandiff + tt
               end do
           end do
           meandiff = meandiff/real(nfvctr_frag,kind=mp)**2

           if (iproc==0) then
               call yaml_mapping_open('analyzing (KS)^2-(KS)')
               call yaml_map('maximal deviation',maxdiff,fmt='(es10.3)')
               call yaml_map('average deviation',meandiff,fmt='(es10.3)')
               call yaml_mapping_close()
           end if

           call f_free(kernel_fragment)
           call f_free(overlap_fragment)
           call f_free(ksk_fragment)
           call f_free(tmpmat)

       end do search_fragments
       if (iproc==0) then
           call yaml_sequence_close()
           call yaml_comment('Fragment purity analysis completed',hfill='~')
       end if
       if (.not.found_a_fragment) then
           call f_err_throw('The specified fragment is not present')
       end if

       call deallocate_sparse_matrix_metadata(smmd)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_l)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)
       call deallocate_matrices(ovrlp_large)
       call deallocate_matrices(KS_large)
       call f_free_str(len(fragment_atomnames),fragment_atomnames)
       call f_free(fragment_atom_id)
       call f_free(fragment_supfun_id)

       call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
                    gather_routine=gather_timings)


   end if



   call build_dict_info(iproc, nproc, dict_timing_info)
   call f_timing_stop(mpi_comm=mpiworld(),nproc=nproc,&
        gather_routine=gather_timings,dict_info=dict_timing_info)
   call dict_free(dict_timing_info)

   if (iproc==0) then
       call yaml_release_document()
   end if


   !!call bigdft_finalize(ierr)
   call mpifinalize()

   call f_lib_finalize()




end program chess_toolbox

!!$!> extract the different wavefunctions to verify if the completeness relation is satisfied
!!$subroutine completeness_relation
!!$
!!$  call wfn_filename(filename_out,radical,binary,ikpt,nspinor,nspin,ispinor,spin,iorb)
!!$
!!$  !loop that has to be done for each of the wavefunctions
!!$  unitwf=99
!!$  call f_open_file(unit=unitwf,file=filename_out)
!!$  call readonewave(unitwf,.not. binary,iorb,iproc,&
!!$       it%lr%d%n1,it%lr%d%n2,it%lr%d%n3, &
!!$       Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
!!$       at,it%lr%wfd,rxyz_old,rxyz,&
!!$       psi_ptr,eval,psifscf)
!!$  call f_close(unitwf)
!!$
!!$end subroutine completeness_relation






subroutine extract_fragment_submatrix(smmd, smat, mat, nat_frag, nfvctr_frag, &
           fragment_atom_id, fragment_supfun_id, mat_frag)
  use sparsematrix_base
  implicit none

  ! Calling arguments
  type(sparse_matrix_metadata),intent(in) :: smmd
  type(sparse_matrix),intent(in) :: smat
  type(matrices),intent(in) :: mat
  integer,intent(in) :: nat_frag, nfvctr_frag
  integer,dimension(nat_frag),intent(in) :: fragment_atom_id
  integer,dimension(smat%nfvctr),intent(in) :: fragment_supfun_id
  real(mp),dimension(nfvctr_frag,nfvctr_frag),intent(out) :: mat_frag

  ! Local variables
  integer :: iicol, icol, icol_atom, iat, iiat, ii, iirow, i, irow, irow_atom, iseg, icol_old
  logical :: found_icol, found_irow

  iicol = 0
  icol_old = -1
  seg_loop: do iseg=1,smat%nseg
      icol = smat%keyg(1,2,iseg)
      icol_atom = smmd%on_which_atom(icol)
      ! Search whether this column belongs to the fragment
      found_icol = .false.
      do iat=1,nat_frag
          iiat = fragment_atom_id(iat)
          if (icol_atom==iiat) then
              found_icol = .true.
          end if
      end do
      if (found_icol) then
          if (icol/=icol_old) then
              iicol = iicol + 1
              icol_old = icol
          end if
      else
          cycle seg_loop
      end if
      ii=smat%keyv(iseg)
      iirow = 0
      do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg) 
          irow = i
          irow_atom = smmd%on_which_atom(irow)
          ! Search whether this column belongs to the fragment
          found_irow = .false.
          do iat=1,nat_frag
              iiat = fragment_atom_id(iat)
              if (irow_atom==iiat) then
                  found_irow = .true.
              end if
          end do
          if (found_irow .and. found_icol) then
              iirow = iirow + 1
              !mat_frag(iirow,iicol) = mat%matrix_compr(ii)
              iirow = fragment_supfun_id(irow)
              iicol = fragment_supfun_id(icol)
              mat_frag(iirow,iicol) = mat%matrix_compr(ii)
          end if
          !write(*,*) 'iirow, iicol, ii, mat_frag, mat_compr', iirow, iicol, ii, mat_frag(iirow,iicol), mat%matrix_compr(ii)
          ii = ii + 1
      end do
  end do seg_loop

end subroutine extract_fragment_submatrix
