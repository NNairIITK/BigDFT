!> @file
!!   Program to guess the used memory by BigDFT
!! @author
!!   Copyright (C) 2007-2013 BigDFT group (LG)
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Test the input files and estimates the memory occupation versus the number
!! of processors
program memguess

   use module_base
   use module_types
   use module_interfaces
   use module_xc
   use m_ab6_symmetry
   use module_fragments
   use yaml_output
   use bigdft_run
   use module_atoms, only: set_astruct_from_file,astruct_dump_to_file
   use internal_coordinates
   use gaussians, only: gaussian_basis, deallocate_gwf
   use communications_base, only: deallocate_comms
   use psp_projectors, only: free_DFT_PSP_projectors
   use io, only: read_linear_matrix_dense, read_coeff_minbasis, writeLinearCoefficients, &
                 read_sparse_matrix, read_linear_coefficients
   use sparsematrix_base, only: sparse_matrix, matrices_null, assignment(=), SPARSE_FULL, &
                                sparsematrix_malloc_ptr, sparsematrix_malloc0_ptr, DENSE_FULL
   use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple
   use sparsematrix, only: uncompress_matrix
   use postprocessing_linear, only: loewdin_charge_analysis_core
                                
   implicit none
   character(len=*), parameter :: subname='memguess'
   character(len=30) :: tatonam, radical
   character(len=2) :: num
   character(len=40) :: comment
   character(len=1024) :: fcomment
   character(len=128) :: fileFrom, fileTo,filename_wfn, coeff_file, ham_file, overlap_file, kernel_file, matrix_file
   character(len=128) :: ntmb_, norbks_, interval_, npdos_, nat_, nsubmatrices_, ncategories_, cutoff_, power_
   character(len=128) :: output_pdos, amatrix_file, bmatrix_file, cmatrix_file, inmatrix_file, outmatrix_file, wf_file
   character(len=128) :: posinp_file, pdos_file, cc
   logical :: optimise,GPUtest,atwf,convert=.false.,exportwf=.false.,logfile=.false.
   logical :: disable_deprecation = .false.,convertpos=.false.,transform_coordinates=.false.
   logical :: calculate_pdos = .false., kernel_analysis = .false., extract_submatrix = .false.
   logical :: solve_eigensystem = .false., analyze_coeffs = .false., peel_matrix = .false.
   logical :: multiply_matrices = .false., matrixpower = .false., plot_wavefunction = .false.
   logical :: suggest_cutoff = .false., charge_analysis = .false.
   integer :: ntimes,nproc,output_grid, i_arg,istat
   integer :: nspin,iorb,norbu,norbd,nspinor,norb,iorbp,iorb_out,lwork
   integer :: norbgpu,ng, nsubmatrices, ncategories
   integer :: export_wf_iband, export_wf_ispin, export_wf_ikpt, export_wf_ispinor,irad
   real(gp) :: hx,hy,hz,energy,occup,interval,tt,cutoff,power,d,occup_pdos, total_occup
   type(memory_estimation) :: mem
   type(run_objects) :: runObj
   type(orbitals_data) :: orbstst
   type(DFT_PSP_projectors) :: nlpsp
   type(gaussian_basis) :: G !basis for davidson IG
   type(atoms_data) :: at
   type(denspot_distribution) :: dpbox
   real(gp), dimension(3) :: shift
   real(gp), dimension(:,:), pointer :: fxyz
   real(wp), dimension(:), allocatable :: rhoexpo
   real(wp), dimension(:,:,:,:), pointer :: rhocoeff
   real(gp), dimension(:), pointer :: gbd_occ
   type(system_fragment), dimension(:), pointer :: ref_frags
   character(len=3) :: in_name !lr408
   character(len=128) :: line
   integer :: i, inputpsi, input_wf_format, nneighbor_min, nneighbor_max, nneighbor, ntypes
   integer,parameter :: nconfig=1
   type(dictionary), pointer :: run
   integer,dimension(:),pointer :: nzatom, nelpsp, iatype
   character(len=20),dimension(:),pointer :: atomnames
   !character(len=60),dimension(nconfig) :: arr_radical,arr_posinp
   !character(len=60) :: run_id, infile, outfile
   !integer, dimension(4) :: mpi_info
   !real(gp) :: tcpu0,tcpu1,tel
   !integer :: ncount0,ncount1,ncount_max,ncount_rate
   !! By Ali
   integer :: ierror, iat, itmb, jtmb, iitmb, jjtmb, ntmb, norbks, npdos, iunit01, iunit02, norb_dummy, ipt, npt, ipdos, nat
   integer :: iproc, isub, jat, icat, info, itype, iiat, jjat, jtype, ios, ival, iat_prev, ii, iitype, ispin
   integer :: nfvctr_s, nseg_s, nvctr_s, nfvctrp_s, isfvctr_s
   integer :: nfvctr_m, nseg_m, nvctr_m, nfvctrp_m, isfvctr_m
   integer :: nfvctr_l, nseg_l, nvctr_l, nfvctrp_l, isfvctr_l
   integer,dimension(:),allocatable :: na, nb, nc, on_which_atom
   integer,dimension(:),pointer :: keyv_s, keyv_m, keyv_l, on_which_atom_s, on_which_atom_m, on_which_atom_l
   integer,dimension(:,:),allocatable :: atoms_ref, imin_list
   integer,dimension(:,:,:),pointer :: keyg_s, keyg_m, keyg_l
   real(kind=8),dimension(:,:),allocatable :: rxyz, rxyz_int, denskernel, ham, overlap, coeff, pdos, energy_bins, matrix
   real(kind=8),dimension(:,:),pointer :: coeff_ptr
   real(kind=8),dimension(:,:),allocatable :: amatrix, bmatrix, cmatrix, temparr, d1min_list
   real(kind=8),dimension(:),allocatable :: eval, coeff_cat, work, d2min_list, dtype, rcov
   real(kind=8),dimension(:),pointer :: eval_ptr
   real(kind=8),dimension(:),pointer :: matrix_compr
   type(matrices) :: ovrlp_mat, ham_mat, kernel_mat
   !real(kind=8),parameter :: degree=57.295779513d0
   real(kind=8),parameter :: degree=1.d0
   character(len=6) :: direction
   character(len=2) :: backslash
   logical :: file_exists, found_bin, mpi_init
   logical,dimension(:,:),allocatable :: calc_array
   real(kind=8),parameter :: eps_roundoff=1.d-5
   type(sparse_matrix) :: smat_s, smat_m, smat_l

   call f_lib_initialize()
   !initialize errors and timings as bigdft routines are called
   call bigdft_init_errors()
   call bigdft_init_timing_categories()
   ! Get arguments
   !call getarg(1,tatonam)
   call get_command_argument(1, value = tatonam, status = istat)

   write(radical, "(A)") "input"
   optimise=.false.
   GPUtest=.false.
   atwf=.false.
   if(trim(tatonam)=='' .or. istat>0) then
      write(*,'(1x,a)')&
         &   'Usage: ./memguess <nproc> [option]'
      write(*,'(1x,a)')&
         &   'Indicate the number of processes after the executable'
      write(*,'(1x,a)')&
         &   '[option] can be the following: '
      write(*,'(1x,a)')&
         &   '"y": grid to be plotted with V_Sim'
      write(*,'(1x,a)')&
         &   '"o" rotate the molecule such that the volume of the simulation box is optimised'
      write(*,'(1x,a)')&
         &   '"GPUtest <nrep>" case of a CUDAGPU calculation, to test the speed of 3d operators'
      write(*,'(1x,a)')&
         &   '         <nrep> is the number of repeats'
      write(*,'(1x,a)')&
         &   '"upgrade" upgrades input files older than 1.2 into actual format'
      write(*,'(1x,a)')&
         &   '"convert" <from.[cube,etsf]> <to.[cube,etsf]>" converts "from" to file "to" using the given formats'
      write(*,'(1x,a)')&
         &   '"exportwf" <n>[u,d] <from.[bin,formatted,etsf]> "'//&
      ' converts n-th wavefunction of file "from" to cube using BigDFT uncompression'
      write(*,'(1x,a)')&
         &   '"atwf" <ng> calculates the atomic wavefunctions of the first atom in the gatom basis and write their expression '
      write(*,'(1x,a)')&
         &   '            in the "gatom-wfn.dat" file '
      write(*,'(1x,a)')&
         &   '           <ng> is the number of gaussians used for the gatom calculation'
      write(*,'(1x,a)')&
           &   '"convert-positions" <from.[xyz,ascii,yaml]> <to.[xyz,ascii,yaml]>" ' 
      write(*,'(1x,a)')&
           & 'converts input positions file "from" to file "to" using the given formats'
      write(*,'(1x,a)')&
           &   '"pdos" <ntmb> <norb> <coeffs.bin> <npdos>" ' 
      write(*,'(1x,a)')&
           & 'reads in the expansion coefficients "coeffs.bin" of dimension (nmtb x norb) &
           &and calculate "npdos" partial density of states'
      write(*,'(1x,a)')&
           &   '"kernel-analysis" <coeffs.bin> <kernel.bin>" ' 
      write(*,'(1x,a)')&
           & 'calculates a full kernel from the expansion coefficients "coeffs.bin" of dimension (nmtb x norb) &
           &and compare it with the sparse kernel in "kernel.bin"'
      write(*,'(1x,a)')&
           &   '"solve-eigensystem" <ham.bin> <overlap.bin> <coeffs.bin>" ' 
      write(*,'(1x,a)')&
           & 'solve the eigensystem Hc = lSc and write the coeffs c to disk'
      write(*,'(1x,a)')&
           &   '"analyse-coeffs" <coeff.bin>" ' 
      write(*,'(1x,a)')&
           & 'analyse the coefficients by assiging them in to ncategories categories'
      write(*,'(1x,a)')&
           &   '"peel-matrix" <matrix.bin>" ' 
      write(*,'(1x,a)')&
           & 'peel a matrix by stripping off elements which are outside of a cutoff'
      write(*,'(1x,a)')&
           &   '"multiply-matrices" <matrix.bin>" ' 
      write(*,'(1x,a)')&
           & 'multiply two matrices'
      write(*,'(1x,a)')&
           &   '"matrixpower" <matrix.bin>" ' 
      write(*,'(1x,a)')&
           & 'caluclate the power of a matrix'
      write(*,'(1x,a)')&
           &   '"suggest-cutoff" <posinp.xyz>" ' 
      write(*,'(1x,a)')&
           & 'suggest cutoff radii for the linear scaling version'
      write(*,'(1x,a)')&
           &   '"charge-analysis"" ' 
      write(*,'(1x,a)')&
           & 'perform a Loewdin charge analysis'

      stop
   else
      read(unit=tatonam,fmt=*,iostat=ierror) nproc
      if (ierror /= 0) then
         call deprecation_message()
         stop
      end if
      i_arg = 2
      output_grid=0
      loop_getargs: do
         call get_command_argument(i_arg, value = tatonam, status = istat)
         !call getarg(i_arg,tatonam)
         if(trim(tatonam)=='' .or. istat > 0) then
            exit loop_getargs
         else if (trim(tatonam)=='y') then
            output_grid=1
            write(*,'(1x,a)') 'The system grid will be displayed in the "grid.xyz" file'
            exit loop_getargs
         else if (trim(tatonam)=='o') then
            optimise=.true.
            output_grid=1
            write(*,'(1x,a)')&
               &   'The optimised system grid will be displayed in the "grid.xyz" file'
            exit loop_getargs
         else if (trim(tatonam)=='GPUtest') then
            GPUtest=.true.
            write(*,'(1x,a)')&
               &   'Perform the test with GPU, if present.'
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            !call getarg(i_arg,tatonam)
            ntimes=1
            norbgpu=0
            read(unit=tatonam,fmt=*,iostat=ierror) ntimes
            if (ierror == 0) then
               write(*,'(1x,a,i0,a)')&
                  &   'Repeat each calculation ',ntimes,' times.'
               i_arg = i_arg + 1
               call get_command_argument(i_arg, value = tatonam)
               !call getarg(i_arg,tatonam)
               read(unit=tatonam,fmt=*,iostat=ierror) norbgpu
            end if
            exit loop_getargs
         else if (trim(tatonam)=='convert') then
            convert=.true.
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fileFrom)
            !call getarg(i_arg,fileFrom)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fileTo)
            !call getarg(i_arg,fileTo)
            write(*,'(1x,5a)')&
               &   'convert "', trim(fileFrom),'" file to "', trim(fileTo),'"'
            exit loop_getargs

         else if (trim(tatonam)=='exportwf') then
            !Export wavefunctions (cube format)
            exportwf=.true.
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = filename_wfn)
            !call getarg(i_arg,filename_wfn)
            ! Read optional additional arguments with the iband, up/down and ikpt
            export_wf_iband = 1
            export_wf_ispin = 1
            export_wf_ikpt  = 1
            export_wf_ispinor = 1
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            if(trim(tatonam)=='' .or. istat > 0) then
               exit loop_getargs
            else
               read(unit=tatonam,fmt=*) export_wf_iband
            end if
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            if(trim(tatonam)=='' .or. istat > 0) then
               exit loop_getargs
            else
               read(unit=tatonam,fmt=*) export_wf_ispin
            end if
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            if(trim(tatonam)=='' .or. istat > 0) then
               exit loop_getargs
            else
               read(unit=tatonam,fmt=*) export_wf_ikpt
            end if
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            if(trim(tatonam)=='' .or. istat > 0) then
               exit loop_getargs
            else
               read(unit=tatonam,fmt=*) export_wf_ispinor
            end if
            exit loop_getargs
         else if (trim(tatonam)=='atwf') then
            atwf=.true.
            write(*,'(1x,a)')&
               &   'Perform the calculation of atomic wavefunction of the first atom'
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam)
            !call getarg(i_arg,tatonam)
            read(unit=tatonam,fmt=*,iostat=ierror) ng
            write(*,'(1x,a,i0,a)')&
               &   'Use gaussian basis of',ng,' elements.'
            exit loop_getargs
         else if (trim(tatonam)=='convert-positions') then
            convertpos=.true.
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fileFrom)
            !call getarg(i_arg,fileFrom)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fileTo)
            !call getarg(i_arg,fileTo)
            write(*,'(1x,5a)')&
               &   'convert input file "', trim(fileFrom),'" file to "', trim(fileTo),'"'
            exit loop_getargs
         else if (trim(tatonam)=='transform-coordinates') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = direction)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fileFrom)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = fileTo)
            write(*,'(1x,5a)')&
               &   'convert input file "', trim(fileFrom),'" file to "', trim(fileTo),'"'
            transform_coordinates=.true.
            exit loop_getargs
         else if (trim(tatonam)=='pdos') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ham_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            !i_arg = i_arg + 1
            !call get_command_argument(i_arg, value = ntmb_)
            !read(ntmb_,fmt=*,iostat=ierror) ntmb
            !i_arg = i_arg + 1
            !call get_command_argument(i_arg, value = norbks_)
            !read(norbks_,fmt=*,iostat=ierror) norbks
            !i_arg = i_arg + 1
            !call get_command_argument(i_arg, value = nat_)
            !read(nat_,fmt=*,iostat=ierror) nat
            !i_arg = i_arg + 1
            !call get_command_argument(i_arg, value = interval_)
            !read(interval_,fmt=*,iostat=ierror) interval
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = pdos_file)
            !i_arg = i_arg + 1
            !call get_command_argument(i_arg, value = posinp_file)
            npdos = 1
            write(*,'(1x,a,i0,3a)')&
               &   'calculate ', npdos,' PDOS based on the coeffs in the file "', trim(coeff_file),'"'
               !&   'calculate ', npdos,' PDOS based on the coeffs (', ntmb, 'x', norbks, ') in the file "', trim(coeff_file),'"'
            calculate_pdos=.true.
            exit loop_getargs
         else if (trim(tatonam)=='kernel-analysis') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ntmb_)
            read(ntmb_,fmt=*,iostat=ierror) ntmb
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = norbks_)
            read(norbks_,fmt=*,iostat=ierror) norbks
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nat_)
            read(nat_,fmt=*,iostat=ierror) nat
            write(*,'(1x,5a)')&
               &   'calculate a full kernel from the coeffs in "', trim(coeff_file), &
               &'" and compres it to the sparse kernel in "', trim(kernel_file),'"'
            kernel_analysis = .true.
            exit loop_getargs
         else if (trim(tatonam)=='extract-submatrix') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ntmb_)
            read(ntmb_,fmt=*,iostat=ierror) ntmb
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nat_)
            read(nat_,fmt=*,iostat=ierror) nat
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nsubmatrices_)
            read(nsubmatrices_,fmt=*,iostat=ierror) nsubmatrices
            write(*,'(1x,a,i0,3a,2(i0,a))')&
               &   'extract ',nsubmatrices,' submatrices from the matrix in "', trim(matrix_file),'" (size ',ntmb,'x',ntmb,')'
            extract_submatrix = .true.
            exit loop_getargs
         else if (trim(tatonam)=='solve-eigensystem') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ham_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ntmb_)
            read(ntmb_,fmt=*,iostat=ierror) ntmb
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nat_)
            read(nat_,fmt=*,iostat=ierror) nat
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            read(nsubmatrices_,fmt=*,iostat=ierror) nsubmatrices
            write(*,'(1x,2(a,i0))')&
               &   'solve the eigensystem Hc=lSc of size ',ntmb,'x',ntmb
            solve_eigensystem = .true.
            exit loop_getargs
         else if (trim(tatonam)=='analyze-coeffs') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = coeff_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ntmb_)
            read(ntmb_,fmt=*,iostat=ierror) ntmb
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = norbks_)
            read(norbks_,fmt=*,iostat=ierror) norbks
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nat_)
            read(nat_,fmt=*,iostat=ierror) nat
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ncategories_)
            read(ncategories_,fmt=*,iostat=ierror) ncategories
            write(*,'(1x,a)')&
               &   'analyze the coeffs'
            analyze_coeffs = .true.
            exit loop_getargs
         else if (trim(tatonam)=='peel-matrix') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = matrix_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ntmb_)
            read(ntmb_,fmt=*,iostat=ierror) ntmb
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nat_)
            read(nat_,fmt=*,iostat=ierror) nat
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = cutoff_)
            read(cutoff_,fmt=*,iostat=ierror) cutoff
            write(*,'(1x,a)')&
               &   'peel the matrix'
            peel_matrix = .true.
            exit loop_getargs
         else if (trim(tatonam)=='multiply-matrices') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = amatrix_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = bmatrix_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = cmatrix_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ntmb_)
            read(ntmb_,fmt=*,iostat=ierror) ntmb
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nat_)
            read(nat_,fmt=*,iostat=ierror) nat
            i_arg = i_arg + 1
            write(*,'(1x,a,2(i0,a))')&
               &   'multiply the matrices (size ',ntmb,'x',ntmb,')'
            multiply_matrices = .true.
            exit loop_getargs
         else if (trim(tatonam)=='matrixpower') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = inmatrix_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = outmatrix_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = ntmb_)
            read(ntmb_,fmt=*,iostat=ierror) ntmb
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = nat_)
            read(nat_,fmt=*,iostat=ierror) nat
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = power_)
            read(power_,fmt=*,iostat=ierror) power
            i_arg = i_arg + 1
            write(*,'(1x,a,2(i0,a))')&
               &   'calculate the power of a matrix'
            matrixpower = .true.
            exit loop_getargs
         else if (trim(tatonam)=='plot-wavefunction') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = wf_file)
            write(*,'(1x,a,2(i0,a))')&
               &   'plot the wave function from file ',trim(wf_file)
            plot_wavefunction = .true.
            exit loop_getargs
         else if (trim(tatonam)=='suggest-cutoff') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = posinp_file)
            write(*,'(1x,2a)')&
               &   'suggest cutoff radii based on the atomic positions in ',trim(posinp_file)
            suggest_cutoff = .true.
            exit loop_getargs
         else if (trim(tatonam)=='charge-analysis') then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = overlap_file)
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = kernel_file)
            write(*,'(1x,2a)')&
               &   'perform a Loewdin charge analysis'
            charge_analysis = .true.
            exit loop_getargs
         else if (trim(tatonam) == 'dd') then
            ! dd: disable deprecation message
            disable_deprecation = .true.
         else if (trim(tatonam) == 'l') then
            ! l: log to disk
            logfile = .true.
         else
            ! Use value as radical for input files.
            write(radical, "(A)") trim(tatonam)
         end if
         i_arg = i_arg + 1
      end do loop_getargs
   end if

   !!!  open(unit=1,file='input.memguess',status='old')
   !!!  
   !!!  !line number, to control the input values
   !!!  iline=0
   !!!  
   !!!  !number of MPI proccessors
   !!!  read(1,*) nproc
   !!!  write(*,*) 'Number of mpi processes is: ',nproc
   !!!  
   !!!  read(1,*) optimise
   !!!  if (optimise) write(*,*) 'Molecule will be rotated to minimize simulation box size and workarrays in BigDFT'
   !!!  
   !!!  !    "T"  If the system grid is to be displayed in the "grid.xyz" file
   !!!  read(1,*) output_grid
   !!!  write(*,*)  'output_grid= ',output_grid
   !!!  
   !!!  !    "T"   'Perform the test with GPU, if present.'   
   !!!  read(1,*) GPUtest
   !!!  if (GPUtest) write(*,*) 'Perform the test with GPU'
   !!!!!! END of By Ali
   if (.not. disable_deprecation) then
      call deprecation_message()
   end if

   !welcome screen
   !call print_logo()

   if (convert) then
      at%astruct%geocode = "P"
      write(*,*) "Read density file..."
      call read_density(trim(fileFrom), at%astruct%geocode, &
           & dpbox%ndims(1), dpbox%ndims(2), dpbox%ndims(3), &
           & nspin, dpbox%hgrids(1), dpbox%hgrids(2), dpbox%hgrids(3), &
           & rhocoeff, at%astruct%nat, at%astruct%rxyz, at%astruct%iatype, at%nzatom)
      at%astruct%ntypes = size(at%nzatom)
      write(*,*) "Write new density file..."
      dpbox%ngatherarr = f_malloc_ptr((/ 0.to.0, 1.to.2 /),id='dpbox%ngatherarr')

      call plot_density(0,1,trim(fileTo),at,at%astruct%rxyz,dpbox,nspin,rhocoeff)
      write(*,*) "Done"
      stop
   end if

   if (convertpos) then
      call set_astruct_from_file(trim(fileFrom),0,at%astruct,fcomment,energy,fxyz)
      
      !find the format of the output file
      if (index(fileTo,'.xyz') > 0) then
         irad=index(fileTo,'.xyz')
         at%astruct%inputfile_format='xyz  '
      else if (index(fileTo,'.ascii') > 0) then
         irad=index(fileTo,'.ascii')
         at%astruct%inputfile_format='ascii'
      else if (index(fileTo,'.yaml') > 0) then
         irad=index(fileTo,'.yaml')
         at%astruct%inputfile_format='yaml '
      else
         irad = len(trim(fileTo)) + 1
      end if
      
      if (associated(fxyz)) then
         call astruct_dump_to_file(at%astruct,fileTo(1:irad-1),&
              trim(fcomment) // ' (converted from '//trim(fileFrom)//")",&
              energy,forces=fxyz)
!!$         call write_atomic_file(fileTo(1:irad-1),energy,at%astruct%rxyz,at%astruct%ixyz_int,at,&
!!$              trim(fcomment) // ' (converted from '//trim(fileFrom)//")", forces=fxyz)

         call f_free_ptr(fxyz)
      else
         call astruct_dump_to_file(at%astruct,fileTo(1:irad-1),&
              trim(fcomment) // ' (converted from '//trim(fileFrom)//")",&
              energy)
!!$
!!$         call write_atomic_file(fileTo(1:irad-1),energy,at%astruct%rxyz,at%astruct%ixyz_int,at,&
!!$              trim(fcomment) // ' (converted from '//trim(fileFrom)//")")
      end if
      stop
   end if

   if (transform_coordinates) then
       if (direction=='carint') then
           write(*,*) 'Converting cartesian coordinates to internal coordinates.'
       else if (direction=='intcar') then
           write(*,*) 'Converting internal coordinates to cartesian coordinates.'
       else if (direction=='carcar') then
           write(*,*) 'Converting cartesian coordinates to cartesian coordinates.'
       else
           call f_err_throw("wrong switch for coordinate transforms", err_name='BIGDFT_RUNTIME_ERROR')
       end if
       call set_astruct_from_file(trim(fileFrom),0,at%astruct,fcomment,energy,fxyz)

       !find the format of the output file
       if (index(fileTo,'.xyz') > 0) then
          irad=index(fileTo,'.xyz')
          at%astruct%inputfile_format='xyz  '
       else if (index(fileTo,'.ascii') > 0) then
          irad=index(fileTo,'.ascii')
          at%astruct%inputfile_format='ascii'
       else if (index(fileTo,'.int') > 0) then
          irad=index(fileTo,'.int')
          at%astruct%inputfile_format='int  '
       else if (index(fileTo,'.yaml') > 0) then
          irad=index(fileTo,'.yaml')
          at%astruct%inputfile_format='yaml '
       else
          irad = len(trim(fileTo)) + 1
       end if

       ! Check whether the output file format is correct
       if (direction=='carint') then
           ! output file must be .int
           if (at%astruct%inputfile_format/='int  ')then
               call f_err_throw("wrong output format for the coordinate transform", err_name='BIGDFT_RUNTIME_ERROR')
           end if
       else if (direction=='intcar' .or. direction=='carcar') then
           ! output file must be .xyz, .ascii or. .yaml
           if (at%astruct%inputfile_format/='xyz  ' .and. &
               at%astruct%inputfile_format/='ascii' .and. &
               at%astruct%inputfile_format/='yaml ')then
               call f_err_throw("wrong output format for the coordinate transform", err_name='BIGDFT_RUNTIME_ERROR')
           end if
       end if


       !!call bigdft_get_run_ids(nconfig,trim(run_id),arr_radical,arr_posinp,ierror)
       !!call run_objects_init_from_files(runObj, radical, posinp)
       na = f_malloc(at%astruct%nat,id='na')
       nb = f_malloc(at%astruct%nat,id='nb')
       nc = f_malloc(at%astruct%nat,id='nc')
       rxyz_int = f_malloc((/ 3, at%astruct%nat /),id='rxyz_int')

       if (direction=='carint') then
           inquire(file='posinp.fix',exist=file_exists)
           if (file_exists) then
               atoms_ref = f_malloc((/3,at%astruct%nat/),id='atoms_ref')
               open(unit=123,file='posinp.fix')
               do iat=1,at%astruct%nat
                   read(123,*) atoms_ref(1:3,iat)
               end do
               call get_neighbors(at%astruct%rxyz, at%astruct%nat, &
                    at%astruct%ixyz_int(1,:), at%astruct%ixyz_int(2,:),at%astruct%ixyz_int(3,:), &
                    atoms_ref)
               call f_free(atoms_ref)
           else
               call get_neighbors(at%astruct%rxyz, at%astruct%nat, &
                    at%astruct%ixyz_int(1,:), at%astruct%ixyz_int(2,:),at%astruct%ixyz_int(3,:))
           end if
           !call xyzint(at%astruct%rxyz, at%astruct%nat, &
           !     at%astruct%ixyz_int(1,:), at%astruct%ixyz_int(2,:),at%astruct%ixyz_int(3,:), &
           !     degree, at%astruct%rxyz_int)
           !! The bond angle must be modified (take 180 degrees minus the angle)
           !at%astruct%rxyz_int(2:2,1:at%astruct%nat) = pi_param - at%astruct%rxyz_int(2:2,1:at%astruct%nat)
           call astruct_dump_to_file(at%astruct,fileTo(1:irad-1),&
                trim(fcomment) // ' (converted from '//trim(fileFrom)//")")

!!$           call write_atomic_file(trim(fileTo(1:irad-1)),UNINITIALIZED(123.d0),at%astruct%rxyz_int,&
!!$                at%astruct%ixyz_int,at,trim(fcomment) // ' (converted from '//trim(fileFrom)//")")
       else if (direction=='intcar' .or. direction=='carcar') then
          call astruct_dump_to_file(at%astruct,fileTo(1:irad-1),&
               trim(fcomment) // ' (converted from '//trim(fileFrom)//")")

!!$           call write_atomic_file(trim(fileTo(1:irad-1)),UNINITIALIZED(123.d0),at%astruct%rxyz,&
!!$                at%astruct%ixyz_int,at,trim(fcomment) // ' (converted from '//trim(fileFrom)//")")
       end if

       write(*,*) 'Done.'

       call f_free(na)
       call f_free(nb)
       call f_free(nc)
       call f_free(rxyz_int)
       stop
   end if

   if (calculate_pdos) then
       call mpi_initialized(mpi_init, ierror)
       if (mpi_init) then
           call mpi_comm_rank(mpi_comm_world, iproc, ierror)
           call mpi_comm_size(mpi_comm_world, nproc, ierror)
       else
           iproc = 0
           nproc = 1
       end if
       !coeff = f_malloc((/ntmb,norbks/),id='coeff')
       !eval = f_malloc(norbks,id='eval')
       !denskernel = f_malloc((/ntmb,ntmb/),id='denskernel')
       !ham = f_malloc((/ntmb,ntmb/),id='ham')
       !overlap = f_malloc((/ntmb,ntmb/),id='overlap')
       !on_which_atom = f_malloc(ntmb,id='on_which_atom')
       !calc_array = f_malloc(ntmb,id='calc_array')
       call f_open_file(iunit01, file=trim(coeff_file), binary=.false.)
       !call read_coeff_minbasis(iunit01, .true., iproc, norbks, norb_dummy, ntmb, coeff, eval, nat)
       !write(*,*) 'trim(coeff_file)',trim(coeff_file)
       call read_linear_coefficients(trim(coeff_file), nspin, ntmb, norbks, coeff_ptr, &
            eval=eval_ptr)
       call f_close(iunit01)
       !write(*,*) 'ntmb',ntmb
       !write(*,*) 'coeff_ptr',coeff_ptr

       !call f_open_file(iunit01, file=ham_file, binary=.false.)
       !call read_linear_matrix_dense(iunit01, ntmb, nat, ham, on_which_atom=on_which_atom)
       !call read_sparse_matrix(trim(ham_file), nspin, nfvctr_m, nseg_m, nvctr_m, keyv_m, keyg_m, &
       !     matrix_compr, at%astruct%nat, at%astruct%ntypes, at%nzatom, at%nelpsp, &
       !     at%astruct%atomnames, at%astruct%iatype, at%astruct%rxyz,  on_which_atom=on_which_atom_m)
       call read_sparse_matrix(trim(ham_file), nspin, nfvctr_m, nseg_m, nvctr_m, keyv_m, keyg_m, &
            matrix_compr, on_which_atom=on_which_atom_m)
       call distribute_columns_on_processes_simple(iproc, nproc, nfvctr_m, nfvctrp_m, isfvctr_m)
       call bigdft_to_sparsebigdft(iproc, nproc, nfvctr_m, nfvctrp_m, isfvctr_m, &
            on_which_atom_m, nvctr_m, nseg_m, keyg_m, smat_m)
       ham_mat = matrices_null()
       ham_mat%matrix = sparsematrix_malloc0_ptr(smat_m,iaction=DENSE_FULL,id='smat_m%matrix')
       call uncompress_matrix(iproc, smat_m, matrix_compr, ham_mat%matrix)
       call f_free_ptr(matrix_compr)
       !call f_close(iunit01)

       !call f_open_file(iunit01, file=overlap_file, binary=.false.)
       !call read_linear_matrix_dense(iunit01, ntmb, nat, overlap)
       call read_sparse_matrix(trim(overlap_file), nspin, nfvctr_s, nseg_s, nvctr_s, keyv_s, keyg_s, &
            matrix_compr, at%astruct%nat, at%astruct%ntypes, at%nzatom, at%nelpsp, &
            at%astruct%atomnames, at%astruct%iatype, at%astruct%rxyz,  on_which_atom=on_which_atom_s)
       !!call read_sparse_matrix(trim(ham_file), nspin, nfvctr_s, nseg_s, nvctr_s, keyv_s, keyg_s, &
       !!     matrix_compr, on_which_atom=on_which_atom_s)
       call distribute_columns_on_processes_simple(iproc, nproc, nfvctr_s, nfvctrp_s, isfvctr_s)
       call bigdft_to_sparsebigdft(iproc, nproc, nfvctr_s, nfvctrp_s, isfvctr_s, &
            on_which_atom_s, nvctr_s, nseg_s, keyg_s, smat_s)
       ovrlp_mat = matrices_null()
       ovrlp_mat%matrix = sparsematrix_malloc0_ptr(smat_s,iaction=DENSE_FULL,id='smat_s%matrix')
       call uncompress_matrix(iproc, smat_s, matrix_compr, ovrlp_mat%matrix)
       call f_free_ptr(matrix_compr)
       !call f_close(iunit01)

       !!!write(*,*) 'trim(kernel_file)',trim(kernel_file)
       !!call read_sparse_matrix(trim(kernel_file), nspin, nfvctr_l, nseg_l, nvctr_l, keyv_l, keyg_l, &
       !!     matrix_compr, at%astruct%nat, at%astruct%ntypes, at%nzatom, at%nelpsp, &
       !!     at%astruct%atomnames, at%astruct%iatype, at%astruct%rxyz,  on_which_atom=on_which_atom_l)
       !!call distribute_columns_on_processes_simple(iproc, nproc, nfvctr_l, nfvctrp_l, isfvctr_l)
       !!call bigdft_to_sparsebigdft(iproc, nproc, nfvctr_l, nfvctrp_l, isfvctr_l, &
       !!     on_which_atom_l, nvctr_l, nseg_l, keyg_l, smat_l)
       !!kernel_mat = matrices_null()
       !!kernel_mat%matrix = sparsematrix_malloc0_ptr(smat_l,iaction=DENSE_FULL,id='smat_s%matrix')
       !!call uncompress_matrix(iproc, smat_l, matrix_compr, kernel_mat%matrix)
       !!call f_free_ptr(matrix_compr)

       !coeff = f_malloc((/ntmb,norbks/),id='coeff')
       !eval = f_malloc(norbks,id='eval')
       denskernel = f_malloc((/ntmb,ntmb/),id='denskernel')
       !ham = f_malloc((/ntmb,ntmb/),id='ham')
       !overlap = f_malloc((/ntmb,ntmb/),id='overlap')
       !on_which_atom = f_malloc(ntmb,id='on_which_atom')

       !call set_astruct_from_file(trim(posinp_file),0,at%astruct,fcomment,energy,fxyz)
       !write(*,*) 'trim(pdos_file)', trim(pdos_file)
       call f_open_file(iunit01, file=pdos_file, binary=.false.)

       read(iunit01,*) npdos

       calc_array = f_malloc((/ntmb,npdos/),id='calc_array')

       calc_array = .false.
       npdos_loop: do ipdos=1,npdos
           do 
               !read(iunit01,*,iostat=ios) cc, ival
               read(iunit01,'(a128)',iostat=ios) line
               if (ios/=0) exit
               write(*,*) 'line',line
               read(line,*,iostat=ios) cc, ival
               if (cc=='#') cycle npdos_loop
               do itype=1,at%astruct%ntypes
                   if (trim(at%astruct%atomnames(itype))==trim(cc)) then
                       iitype = itype
                       exit
                   end if
               end do
               !write(*,'(a,i7,a,2i7)') 'ipdos, cc, ival, iitype', ipdos, trim(cc), ival, iitype
               iat_prev = -1
               do itmb=1,ntmb
                   iat = on_which_atom_s(itmb)
                   if (iat/=iat_prev) then
                       ii = 0
                   end if
                   iat_prev = iat
                   itype = at%astruct%iatype(iat)
                   ii = ii + 1
                   if (itype==iitype .and. ii==ival) then
                       if (calc_array(itmb,ipdos)) stop 'calc_array(itmb)'
                       calc_array(itmb,ipdos) = .true.
                   end if
               end do
           end do

       end do npdos_loop

       !do ipdos=1,npdos
       !    do itmb=1,ntmb
       !        iat = on_which_atom_m(itmb)
       !        itype = at%astruct%iatype(iat)
       !        write(*,'(a,4i8,l5)') 'ipdos, itmb, iat, itype, calc_array(itmb,ipdos)', &
       !            ipdos, itmb, iat, itype, calc_array(itmb,ipdos)
       !    end do
       !end do

       call f_close(iunit01)

       !!npt = ceiling((eval_ptr(ntmb)-eval_ptr(1))/interval)
       pdos = f_malloc0((/npt,npdos/),id='pdos')
       !!energy_bins = f_malloc((/2,npt/),id='energy_bins')
       !!! Determine the energy bins
       !!do ipt=1,npt
       !!    energy_bins(1,ipt) = eval_ptr(1) + real(ipt-1,kind=8)*interval - eps_roundoff
       !!    energy_bins(2,ipt) = energy_bins(1,ipt) + interval
       !!end do
       output_pdos='PDoS.gp'
       call yaml_map('output file',trim(output_pdos))
       iunit02 = 99
       call f_open_file(iunit02, file=trim(output_pdos), binary=.false.)
       write(iunit02,'(a)') '# plot the DOS as a sum of Gaussians'
       write(iunit02,'(a)') 'set samples 1000'
       write(iunit02,'(a,2(es12.5,a))') 'set xrange[',eval_ptr(1),':',eval_ptr(ntmb),']'
       write(iunit02,'(a)') 'sigma=0.01'
       write(backslash,'(a)') '\ '
       ! Calculate a partial kernel for each KS orbital
       total_occup = 0.d0
       do ipdos=1,npdos
           call yaml_map('PDoS number',ipdos)
           call yaml_map('start, increment',(/ipdos,npdos/))
           write(num,fmt='(i2.2)') ipdos
           write(iunit02,'(a,i0,a)') 'f',ipdos,'(x) = '//trim(backslash)
           occup_pdos = 0.d0
           do iorb=1,norbks
               call yaml_map('orbital being processed',iorb)
               call gemm('n', 't', ntmb, ntmb, 1, 1.d0, coeff_ptr(1,iorb), ntmb, &
                    coeff_ptr(1,iorb), ntmb, 0.d0, denskernel(1,1), ntmb)
                !write(*,*) 'denskernel',denskernel
                !write(*,*) 'ovrlp_mat%matrix',ovrlp_mat%matrix
               energy = 0.d0
               occup = 0.d0
               do ispin=1,nspin
                   !$omp parallel default(none) &
                   !$omp shared(ispin,ntmb,denskernel,ham_mat,energy) &
                   !$omp private(itmb,jtmb)
                   !$omp do reduction(+:energy)
                   do itmb=1,ntmb
                       do jtmb=1,ntmb
                           energy = energy + denskernel(itmb,jtmb)*ham_mat%matrix(jtmb,itmb,ispin)
                       end do
                   end do
                   !$omp end do
                   !$omp end parallel
               end do
               do ispin=1,nspin
                   !!$omp do reduction(+:occup)
                   do itmb=1,ntmb!ipdos,ntmb,npdos
                       if (.not.calc_array(itmb,ipdos)) cycle
                       do jtmb=1,ntmb!ipdos,ntmb,npdos
                           if (.not.calc_array(jtmb,ipdos)) cycle
                           occup = occup + denskernel(itmb,jtmb)*ovrlp_mat%matrix(jtmb,itmb,ispin)
                       end do
                   end do
                   !!$omp end do
               end do
               !write(*,*) 'OCCUP',occup, energy, eval_ptr(iorb)
               !!$omp end parallel
               !!found_bin = .false.
               !!do ipt=1,npt
               !!    if (energy>=energy_bins(1,ipt) .and. energy<energy_bins(2,ipt)) then
               !!        pdos(ipt,ipdos) = pdos(ipt,ipdos) + occup
               !!        found_bin = .true.
               !!        exit
               !!    end if
               !!end do
               !!if (.not.found_bin) then
               !!    call f_err_throw('could not determine energy bin, energy='//yaml_toa(energy),err_name='BIGDFT_RUNTIME_ERROR')
               !!end if
               occup_pdos = occup_pdos + occup
               if (iorb<norbks) then
                   write(iunit02,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2)) + '//trim(backslash)
               else
                   write(iunit02,'(2(a,es16.9),a)') '  ',occup,'*exp(-(x-',energy,')**2/(2*sigma**2))'
               end if
               !write(*,'(a,i6,3es16.8)')'iorb, eval(iorb), energy, occup', iorb, eval(iorb), energy, occup
           end do
           total_occup = total_occup + occup_pdos
           if (ipdos==1) then
               write(iunit02,'(a,i0,a)') "plot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           else
               write(iunit02,'(a,i0,a)') "replot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           end if
           call yaml_map('sum of PDoS',occup_pdos)
           !!output_pdos='PDoS_'//num//'.dat'
           !!call yaml_map('output file',trim(output_pdos))
           !!call f_open_file(iunit01, file=trim(output_pdos), binary=.false.)
           !!write(iunit01,'(a)') '#             energy                pdos'
           !!do ipt=1,npt
           !!    write(iunit01,'(2es20.12)') energy_bins(1,ipt), pdos(ipt,ipdos)
           !!end do
           !!call f_close(iunit01)
       end do
       call f_close(iunit02)
       call yaml_map('sum of total DoS',total_occup)

       stop
   end if

   if (kernel_analysis) then
       call mpi_initialized(mpi_init, ierror)
       if (mpi_init) then
           call mpi_comm_rank(mpi_comm_world, iproc, ierror)
           call mpi_comm_size(mpi_comm_world, nproc, ierror)
       else
           iproc = 0
           nproc = 1
       end if
       coeff = f_malloc((/ntmb,norbks/),id='coeff')
       eval = f_malloc(norbks,id='eval')
       denskernel = f_malloc((/ntmb,ntmb/),id='denskernel')
       rxyz = f_malloc((/3,nat/),id='rxyz')
       on_which_atom = f_malloc(ntmb,id='on_which_atom')
       call f_open_file(iunit01, file=trim(coeff_file), binary=.false.)
       call read_coeff_minbasis(iunit01, .true., iproc, norbks, norb_dummy, ntmb, coeff, eval, nat, rxyz)
       call f_close(iunit01)
       call f_open_file(iunit01, file=trim(kernel_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, denskernel, on_which_atom=on_which_atom)
       call f_close(iunit01)
       call analyze_kernel(ntmb, norbks, nat, coeff, denskernel, rxyz, on_which_atom)
       stop
   end if

   if (extract_submatrix) then
       if (mod(ntmb,nsubmatrices)/=0) then
           call f_err_throw('nsubmatrices must be a divisor of the number of basis function',&
                err_name='BIGDFT_RUNETIME_ERROR')
       end if
       matrix = f_malloc((/ntmb,ntmb/),id='matrix')
       on_which_atom = f_malloc(ntmb,id='on_which_atom')
       rxyz = f_malloc((/3,nat/),id='rxyz')
       call f_open_file(iunit01, file=trim(matrix_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, matrix, rxyz=rxyz, on_which_atom=on_which_atom)
       call f_close(iunit01)
       do isub=1,nsubmatrices
           write(num,fmt='(i2.2)') isub
           call f_open_file(iunit01, file=trim(matrix_file)//'_sub'//num, binary=.false.)
           write(iunit01,'(a,2i10,a)') '#  ',ntmb/nsubmatrices, nat, &
               '    number of basis functions, number of atoms'
           do iat=1,nat
                   write(iunit01,'(a,3es24.16)') '#  ',rxyz(1:3,iat)
           end do
           iitmb = 0
           do itmb=isub,ntmb,nsubmatrices
               iitmb = iitmb + 1
               iat = on_which_atom(itmb)
               jjtmb = 0
               do jtmb=isub,ntmb,nsubmatrices
                   jjtmb = jjtmb + 1
                   jat = on_which_atom(jtmb)
                   write(iunit01,'(2(i6,1x),e19.12,2(1x,i6))') iitmb,jjtmb,matrix(itmb,jtmb),iat,jat
               end do
           end do
           call f_close(iunit01)
       end do
       stop
   end if

   if (solve_eigensystem) then
       call mpi_initialized(mpi_init, ierror)
       if (mpi_init) then
           call mpi_comm_rank(mpi_comm_world, iproc, ierror)
           call mpi_comm_size(mpi_comm_world, nproc, ierror)
       else
           iproc = 0
           nproc = 1
       end if
       ham = f_malloc((/ntmb,ntmb/),id='ham')
       overlap = f_malloc((/ntmb,ntmb/),id='overlap')
       eval = f_malloc(ntmb,id='eval')
       rxyz = f_malloc((/3,nat/),id='rxyz')
       call f_open_file(iunit01, file=trim(ham_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, ham, rxyz=rxyz)
       call f_close(iunit01)
       call f_open_file(iunit01, file=trim(overlap_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, overlap)
       call f_close(iunit01)
       call diagonalizeHamiltonian2(iproc, ntmb, ham, overlap, eval)
       call f_open_file(iunit01, file=trim(coeff_file), binary=.false.)
       call writeLinearCoefficients(iunit01, .true., nat, rxyz, &
            ntmb, ntmb, ntmb, ham, eval)
       stop
   end if

   if (analyze_coeffs) then
       call mpi_initialized(mpi_init, ierror)
       if (mpi_init) then
           call mpi_comm_rank(mpi_comm_world, iproc, ierror)
           call mpi_comm_size(mpi_comm_world, nproc, ierror)
       else
           iproc = 0
           nproc = 1
       end if
       coeff = f_malloc((/ntmb,norbks/),id='coeff')
       coeff_cat = f_malloc(ncategories,id='coeff_cat')
       eval = f_malloc(ntmb,id='eval')
       call f_open_file(iunit01, file=trim(coeff_file), binary=.false.)
       call read_coeff_minbasis(iunit01, .true., iproc, norbks, norb_dummy, ntmb, coeff, eval, nat)
       call f_close(iunit01)
       do iorb=1,norbks
           do icat=1,ncategories
               tt = 0.d0
               do itmb=icat,ntmb,ncategories
                   tt = tt + coeff(itmb,iorb)**2
               end do
               coeff_cat(icat) = tt
           end do
           write(*,'(a,i8,4es12.4)') 'iorb, vals', iorb, coeff_cat(1:ncategories)
       end do
       stop
   end if

   if (peel_matrix) then
       matrix = f_malloc((/ntmb,ntmb/),id='matrix')
       on_which_atom = f_malloc(ntmb,id='on_which_atom')
       rxyz = f_malloc((/3,nat/),id='rxyz')
       call f_open_file(iunit01, file=trim(matrix_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, matrix, rxyz=rxyz, on_which_atom=on_which_atom)
       call f_close(iunit01)
       nneighbor_min = huge(nneighbor_min)
       nneighbor_max = -huge(nneighbor_max)
       do itmb=1,ntmb
           iat = on_which_atom(itmb)
           nneighbor = 0
           do jtmb=1,ntmb
               jat = on_which_atom(jtmb)
               tt = sqrt((rxyz(1,jat)-rxyz(1,iat))**2 + &
                         (rxyz(2,jat)-rxyz(2,iat))**2 + &
                         (rxyz(3,jat)-rxyz(3,iat))**2)
               if (tt>cutoff) then
                   matrix(jtmb,itmb) = 0.d0
               else
                   nneighbor = nneighbor + 1
               end if
           end do
           write(*,*) 'itmb, nneighbor', itmb, nneighbor
           nneighbor_min = min (nneighbor_min,nneighbor)
           nneighbor_max = max (nneighbor_max,nneighbor)
       end do
       call yaml_map('min number of neighbors',nneighbor_min)
       call yaml_map('max number of neighbors',nneighbor_max)
       call f_open_file(iunit01, file=trim(matrix_file)//'_peeled', binary=.false.)
       write(iunit01,'(a,2i10,a)') '#  ',ntmb, nat, &
           '    number of basis functions, number of atoms'
       do iat=1,nat
               write(iunit01,'(a,3es24.16)') '#  ',rxyz(1:3,iat)
       end do
       do itmb=1,ntmb
           iat = on_which_atom(itmb)
           do jtmb=1,ntmb
               jat = on_which_atom(jtmb)
               write(iunit01,'(2(i6,1x),e19.12,2(1x,i6))') itmb,jtmb,matrix(itmb,jtmb),iat,jat
           end do
       end do
       call f_close(iunit01)
       stop
   end if

   if (multiply_matrices) then
       amatrix = f_malloc((/ntmb,ntmb/),id='amatrix')
       bmatrix = f_malloc((/ntmb,ntmb/),id='bmatrix')
       cmatrix = f_malloc((/ntmb,ntmb/),id='cmatrix')
       on_which_atom = f_malloc(ntmb,id='on_which_atom')
       rxyz = f_malloc((/3,nat/),id='rxyz')
       call f_open_file(iunit01, file=trim(amatrix_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, amatrix, rxyz=rxyz, on_which_atom=on_which_atom)
       call f_close(iunit01)
       call f_open_file(iunit01, file=trim(bmatrix_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, bmatrix, rxyz=rxyz, on_which_atom=on_which_atom)
       call f_close(iunit01)
       call gemm('n', 'n', ntmb, ntmb, ntmb, 1.d0, amatrix(1,1), ntmb, &
            bmatrix(1,1), ntmb, 0.d0, cmatrix(1,1), ntmb)
       call f_open_file(iunit01, file=trim(cmatrix_file), binary=.false.)
       write(iunit01,'(a,2i10,a)') '#  ',ntmb, nat, &
           '    number of basis functions, number of atoms'
       do iat=1,nat
               write(iunit01,'(a,3es24.16)') '#  ',rxyz(1:3,iat)
       end do
!cmatrix=0.d0
       do itmb=1,ntmb
!cmatrix(itmb,itmb)=1.d0
           iat = on_which_atom(itmb)
           do jtmb=1,ntmb
               jat = on_which_atom(jtmb)
               write(iunit01,'(2(i6,1x),e19.12,2(1x,i6))') itmb,jtmb,cmatrix(itmb,jtmb),iat,jat
           end do
       end do
       call f_close(iunit01)
       stop
   end if

   if (matrixpower) then
       amatrix = f_malloc((/ntmb,ntmb/),id='amatrix')
       bmatrix = f_malloc((/ntmb,ntmb/),id='bmatrix')
       on_which_atom = f_malloc(ntmb,id='on_which_atom')
       rxyz = f_malloc((/3,nat/),id='rxyz')
       eval = f_malloc(ntmb,id='eval')
       call f_open_file(iunit01, file=trim(inmatrix_file), binary=.false.)
       call read_linear_matrix_dense(iunit01, ntmb, nat, amatrix, rxyz=rxyz, on_which_atom=on_which_atom)
       call f_close(iunit01)


       lwork = 10*ntmb
       work = f_malloc(lwork,id='work')
       tempArr=f_malloc((/ntmb,ntmb/), id='tempArr')
       call dsyev('v', 'l', ntmb, amatrix(1,1), ntmb, eval, work, lwork, info)
       do itmb=1,ntmb
          do jtmb=1,ntmb
             tempArr(jtmb,itmb)=amatrix(jtmb,itmb)*eval(itmb)**power
          end do
       end do
       call gemm('n', 't', ntmb, ntmb, ntmb, 1.d0, amatrix(1,1), &
            ntmb, tempArr(1,1), ntmb, 0.d0, bmatrix(1,1), ntmb)
       call f_free(tempArr)
       call f_open_file(iunit01, file=trim(outmatrix_file), binary=.false.)
       write(iunit01,'(a,2i10,a)') '#  ',ntmb, nat, &
           '    number of basis functions, number of atoms'
       do iat=1,nat
               write(iunit01,'(a,3es24.16)') '#  ',rxyz(1:3,iat)
       end do
!cmatrix=0.d0
       do itmb=1,ntmb
!cmatrix(itmb,itmb)=1.d0
           iat = on_which_atom(itmb)
           do jtmb=1,ntmb
               jat = on_which_atom(jtmb)
               write(iunit01,'(2(i6,1x),e19.12,2(1x,i6))') itmb,jtmb,bmatrix(itmb,jtmb),iat,jat
           end do
       end do
       call f_close(iunit01)
       stop
   end if

   if (plot_wavefunction) then
       stop
   end if

   if (suggest_cutoff) then
       call set_astruct_from_file(trim(posinp_file),0,at%astruct,fcomment,energy,fxyz)
       rcov = f_malloc(at%astruct%ntypes,id='rcov')
       rcov(1) = 3.30d0
       rcov(2) = 2.50d0
       rcov(3) = 1.45d0
       rcov(4) = 1.42d0
       rcov(5) = 0.75d0
       !rcov(1) = 1.45d0

       d1min_list = f_malloc((/2,at%astruct%nat/),id='d2min_list')
       d2min_list = f_malloc(at%astruct%nat,id='d1min_list')
       imin_list = f_malloc((/2,at%astruct%nat/),id='imin_list')
       dtype = f_malloc(at%astruct%ntypes,id='dtype')

       d1min_list = huge(d1min_list)
       imin_list = 0
       do iat=1,at%astruct%nat
           itype = at%astruct%iatype(iat)
           do jat=1,at%astruct%nat
               jtype = at%astruct%iatype(jat)
               if (jat/=iat) then
                   !if (rcov(jtype)<=rcov(itype)) then
                       d = (at%astruct%rxyz(1,jat)-at%astruct%rxyz(1,iat))**2 + &
                           (at%astruct%rxyz(2,jat)-at%astruct%rxyz(2,iat))**2 + &
                           (at%astruct%rxyz(3,jat)-at%astruct%rxyz(3,iat))**2
                       d = sqrt(d)
                   !else
                   !    d = 3.d0*rcov(itype)
                   !end if
                   if (d<d1min_list(1,iat)) then
                       d1min_list(2,iat) = d1min_list(1,iat)
                       d1min_list(1,iat) = d
                       imin_list(2,iat) = imin_list(1,iat)
                       imin_list(1,iat) = jat
                   else if (d<d1min_list(2,iat)) then
                       d1min_list(2,iat) = d
                       imin_list(2,iat) = jat
                   end if
               end if
           end do
       end do

       d2min_list = huge(d2min_list)
       do iat=1,at%astruct%nat
           iiat = imin_list(1,iat)
           itype = at%astruct%iatype(iat)
           iitype = at%astruct%iatype(iiat)
           write(*,'(a,i5,2es12.4)') 'iat, d1min_list(1:2,iat)', iat, d1min_list(1:2,iat)
           do jat=1,2
               jjat = imin_list(jat,iiat)
               jtype = at%astruct%iatype(jjat)
               if (jjat/=iat) then
                   !if (rcov(jtype)<=rcov(itype)) then
                   if (rcov(iitype)<=rcov(itype)) then
                       d = (at%astruct%rxyz(1,iat)-at%astruct%rxyz(1,jjat))**2 + &
                           (at%astruct%rxyz(2,iat)-at%astruct%rxyz(2,jjat))**2 + &
                           (at%astruct%rxyz(3,iat)-at%astruct%rxyz(3,jjat))**2
                       d = sqrt(d)
                   else
                       d = 3.d0*rcov(itype)
                   end if
                   write(*,'(a,5i8,es12.4)') 'itype, iat, iiat, jat, jjat, d', itype, iat, iiat, jat, jjat, d
                   d2min_list(iat) = d
               end if
           end do
       end do

       dtype = -huge(dtype)
       do iat=1,at%astruct%nat
           itype = at%astruct%iatype(iat)
           d = d2min_list(iat)
           dtype(itype) = max(d,dtype(itype))
       end do

       do itype=1,at%astruct%ntypes
           write(*,'(a,i7,a,f7.2,es12.4)') 'itype, name, rcov dtype(itype)', &
               itype, at%astruct%atomnames(itype), rcov(itype), dtype(itype)
       end do

       stop
   end if

   if (charge_analysis) then
       call mpi_initialized(mpi_init, ierror)
       if (mpi_init) then
           call mpi_comm_rank(mpi_comm_world, iproc, ierror)
           call mpi_comm_size(mpi_comm_world, nproc, ierror)
       else
           iproc = 0
           nproc = 1
       end if
       !call set_astruct_from_file(trim(posinp_file),0,at%astruct,fcomment,energy,fxyz)

       call read_sparse_matrix(trim(overlap_file), nspin, nfvctr_s, nseg_s, nvctr_s, keyv_s, keyg_s, &
            matrix_compr, at%astruct%nat, at%astruct%ntypes, at%nzatom, at%nelpsp, &
            at%astruct%atomnames, at%astruct%iatype, at%astruct%rxyz,  on_which_atom=on_which_atom_s)
       call distribute_columns_on_processes_simple(iproc, nproc, nfvctr_s, nfvctrp_s, isfvctr_s)
       call bigdft_to_sparsebigdft(iproc, nproc, nfvctr_s, nfvctrp_s, isfvctr_s, &
            on_which_atom_s, nvctr_s, nseg_s, keyg_s, smat_s)
       ovrlp_mat = matrices_null()
       ovrlp_mat%matrix_compr = sparsematrix_malloc_ptr(smat_s, iaction=SPARSE_FULL, id='ovrlp%matrix_compr')
       call vcopy(smat_s%nvctr, matrix_compr(1), 1, ovrlp_mat%matrix_compr(1), 1)
       call f_free_ptr(matrix_compr)

       call read_sparse_matrix(trim(kernel_file), nspin, nfvctr_l, nseg_l, nvctr_l, keyv_l, keyg_l, &
            matrix_compr, on_which_atom=on_which_atom_l)
       call distribute_columns_on_processes_simple(iproc, nproc, nfvctr_l, nfvctrp_l, isfvctr_l)
       call bigdft_to_sparsebigdft(iproc, nproc, nfvctr_l, nfvctrp_l, isfvctr_l, &
            on_which_atom_l, nvctr_l, nseg_l, keyg_l, smat_l)
       kernel_mat = matrices_null()
       kernel_mat%matrix_compr = sparsematrix_malloc_ptr(smat_l, iaction=SPARSE_FULL, id='kernel_mat%matrix_compr')
       call vcopy(smat_l%nvctr, matrix_compr(1), 1, kernel_mat%matrix_compr(1), 1)
       call f_free_ptr(matrix_compr)

       call loewdin_charge_analysis_core(iproc, nproc, smat_s%nfvctr, smat_s%nfvctrp, smat_s%isfvctr, &
            smat_s%nfvctr_par, smat_s%isfvctr_par, meth_overlap=0, &
            smats=smat_s, smatl=smat_l, atoms=at, kernel=kernel_mat, ovrlp=ovrlp_mat)
       stop
   end if

   nullify(run)
   call bigdft_set_run_properties(run, run_id = trim(radical), run_from_files = .true., log_to_disk = logfile)

   call run_objects_init(runObj, run)
   call dict_free(run)

   if (optimise) then
      if (runObj%atoms%astruct%geocode =='F') then
         call optimise_volume(runObj%atoms,&
              & runObj%inputs%crmult,runObj%inputs%frmult,&
              & runObj%inputs%hx,runObj%inputs%hy,runObj%inputs%hz,&
              & runObj%atoms%astruct%rxyz)
      else
         call shift_periodic_directions(runObj%atoms,runObj%atoms%astruct%rxyz,runObj%atoms%radii_cf)
      end if
      write(*,'(1x,a)')'Writing optimised positions in file posopt.[xyz,ascii]...'
      write(comment,'(a)')'POSITIONS IN OPTIMIZED CELL '

      call astruct_dump_to_file(at%astruct,'posopt',trim(comment))

!!$      call write_atomic_file('posopt',0.d0,runObj%atoms%astruct%rxyz,&
!!$           runObj%atoms%astruct%ixyz_int,runObj%atoms,trim(comment))
      !call wtxyz('posopt',0.d0,rxyz,atoms,trim(comment))
   end if


   call print_dft_parameters(runObj%inputs,runObj%atoms)

   !Time initialization
   !call cpu_time(tcpu0)
   !call system_clock(ncount0,ncount_rate,ncount_max)

   inputpsi = runObj%inputs%inputPsiId
   call system_initialization(0, nproc, .true.,inputpsi, input_wf_format, .true., &
        & runObj%inputs, runObj%atoms, runObj%atoms%astruct%rxyz, runObj%rst%GPU%OCLconv, &
        & runObj%rst%KSwfn%orbs, runObj%rst%tmb%npsidim_orbs, runObj%rst%tmb%npsidim_comp, &
        & runObj%rst%tmb%orbs, runObj%rst%KSwfn%Lzd, runObj%rst%tmb%Lzd, nlpsp, runObj%rst%KSwfn%comms, &
        & shift, ref_frags, output_grid = (output_grid > 0))
   call MemoryEstimator(nproc,runObj%inputs%idsx,runObj%rst%KSwfn%Lzd%Glr,&
        & runObj%rst%KSwfn%orbs%norb,runObj%rst%KSwfn%orbs%nspinor,&
        & runObj%rst%KSwfn%orbs%nkpts,nlpsp%nprojel,&
        runObj%inputs%nspin,runObj%inputs%itrpmax,runObj%inputs%iscf,mem)
   
   if (.not. exportwf) then
      call print_memory_estimation(mem)
   else
      runObj%rst%KSwfn%psi = f_malloc_ptr((runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_c+&
           & 7*runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_f)*runObj%rst%KSwfn%orbs%nspinor,&
           id='runObj%rst%KSwfn%psi')

      ! Optionally compute iorbp from arguments in case of ETSF.
      if (export_wf_ikpt < 1 .or. export_wf_ikpt > runObj%rst%KSwfn%orbs%nkpts) stop "Wrong k-point"
      if (export_wf_ispin < 1 .or. export_wf_ispin > runObj%rst%KSwfn%orbs%nspin) stop "Wrong spin"
      if ((export_wf_ispin == 1 .and. &
           & (export_wf_iband < 1 .or. export_wf_iband > runObj%rst%KSwfn%orbs%norbu)) .or. &
           & (export_wf_ispin == 0 .and. &
           & (export_wf_iband < 1 .or. export_wf_iband > runObj%rst%KSwfn%orbs%norbd))) stop "Wrong orbital"
      iorbp = (export_wf_ikpt - 1) * runObj%rst%KSwfn%orbs%norb + &
           & (export_wf_ispin - 1) * runObj%rst%KSwfn%orbs%norbu + export_wf_iband

      ! ref_frags to be allocated here
      i = index(filename_wfn, "/",back=.true.)+1
      read(filename_wfn(i:i+3),*) in_name ! lr408
      if (in_name == 'min') then
         stop 'Ref fragment not initialized, linear reading currently nonfunctional, to be fixed'
      end if

      call yaml_map("Export wavefunction from file", trim(filename_wfn))
      ! @todo Very ugly patch for ref_frags that is nullified by system_initialization
      ! but used as an allocated array by take_psi_from_file().
      ! TO BE CORRECTED !!!!!
      if (.not.associated(ref_frags)) allocate(ref_frags(runObj%inputs%frag%nfrag_ref))
      call take_psi_from_file(filename_wfn,runObj%inputs%frag, &
           & runObj%inputs%hx,runObj%inputs%hy,runObj%inputs%hz,runObj%rst%KSwfn%Lzd%Glr, &
           & runObj%atoms,runObj%atoms%astruct%rxyz,runObj%rst%KSwfn%orbs,runObj%rst%KSwfn%psi,&
           & iorbp,export_wf_ispinor,ref_frags)
      call filename_of_iorb(.false.,"wavefunction",runObj%rst%KSwfn%orbs,iorbp, &
           & export_wf_ispinor,filename_wfn,iorb_out)

      call plot_wf(.false.,filename_wfn,1,runObj%atoms,1.0_wp,runObj%rst%KSwfn%Lzd%Glr, &
           & runObj%inputs%hx,runObj%inputs%hy,runObj%inputs%hz,runObj%atoms%astruct%rxyz, &
           & runObj%rst%KSwfn%psi((runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_c+&
           & 7*runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_f) * (export_wf_ispinor - 1) + 1:))
      deallocate(ref_frags)
      nullify(ref_frags)
   end if

   if (GPUtest) then
      !test the hamiltonian in CPU or GPU
      !create the orbitals data structure for one orbital
      !test orbitals
      nspin=1
      if (norbgpu == 0) then
         norb=runObj%rst%KSwfn%orbs%norb
      else
         norb=norbgpu
      end if
      norbu=norb
      norbd=0
      nspinor=1

      call orbitals_descriptors(0,nproc,norb,norbu,norbd,runObj%inputs%nspin,nspinor, &
           runObj%inputs%gen_nkpt,runObj%inputs%gen_kpt,runObj%inputs%gen_wkpt,orbstst,LINEAR_PARTITION_NONE)
      orbstst%eval = f_malloc_ptr(orbstst%norbp,id='orbstst%eval')
      do iorb=1,orbstst%norbp
         orbstst%eval(iorb)=-0.5_gp
      end do

      do iorb=1,orbstst%norb
         orbstst%occup(iorb)=1.0_gp
         orbstst%spinsgn(iorb)=1.0_gp
      end do

      call check_linear_and_create_Lzd(0,1,runObj%inputs%linear,runObj%rst%KSwfn%Lzd,&
           & runObj%atoms,orbstst,runObj%inputs%nspin,runObj%atoms%astruct%rxyz)

      !for the given processor (this is only the cubic strategy)
      orbstst%npsidim_orbs=(runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_c+&
           & 7*runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_f)*orbstst%norbp*orbstst%nspinor
      orbstst%npsidim_comp=1


      call compare_cpu_gpu_hamiltonian(0,1,runObj%inputs%matacc,runObj%atoms,&
           orbstst,nspin,runObj%inputs%ncong,runObj%inputs%ixc,&
           runObj%rst%KSwfn%Lzd,hx,hy,hz,runObj%atoms%astruct%rxyz,ntimes)

      call deallocate_orbs(orbstst)


      call f_free_ptr(orbstst%eval)

   end if

   if (atwf) then
      !here the treatment of the AE Core charge density
      !number of gaussians defined in the input of memguess
      !ng=31
      !plot the wavefunctions for the pseudo atom
      nullify(G%rxyz)
      call gaussian_pswf_basis(ng,.false.,0,runObj%inputs%nspin,runObj%atoms,runObj%atoms%astruct%rxyz,G,gbd_occ)
      !for the moment multiply the number of coefficients for each channel
      rhocoeff = f_malloc_ptr((/ (ng*(ng+1))/2, 4, 1, 1 /),id='rhocoeff')
      rhoexpo = f_malloc((ng*(ng+1))/2,id='rhoexpo')

      call plot_gatom_basis('gatom',1,ng,G,gbd_occ,rhocoeff,rhoexpo)

      if (associated(gbd_occ)) then
         call f_free_ptr(gbd_occ)
         nullify(gbd_occ)
      end if
      !deallocate the gaussian basis descriptors
      call deallocate_gwf(G)

      !!$  !plot the wavefunctions for the AE atom
      !!$  !not possible, the code should recognize the AE eleconf
      !!$  call to_zero(35,runObj%atoms%psppar(0,0,runObj%atoms%astruct%iatype(1)))
      !!$  runObj%atoms%psppar(0,0,runObj%atoms%astruct%iatype(1))=0.01_gp
      !!$  nullify(G%rxyz)
      !!$  call gaussian_pswf_basis(ng,.false.,0,runObj%inputs%nspin,atoms,runObj%atoms%astruct%rxyz,G,gbd_occ)
      !!$  !for the moment multiply the number of coefficients for each channel
      !!$  allocate(rhocoeff((ng*(ng+1))/2,4+ndebug),stat=i_stat)
      !!$  call memocc(i_stat,rhocoeff,'rhocoeff',subname)
      !!$  allocate(rhoexpo((ng*(ng+1))/2+ndebug),stat=i_stat)
      !!$  call memocc(i_stat,rhoexpo,'rhoexpo',subname)
      !!$  
      !!$  call plot_gatom_basis('all-elec',1,ng,G,gbd_occ,rhocoeff,rhoexpo)
      !!$
      !!$  if (associated(gbd_occ)) then
      !!$     i_all=-product(shape(gbd_occ))*kind(gbd_occ)
      !!$     deallocate(gbd_occ,stat=i_stat)
      !!$     call memocc(i_stat,i_all,'gbd_occ',subname)
      !!$     nullify(gbd_occ)
      !!$  end if
      !!$  !deallocate the gaussian basis descriptors
      !!$  call deallocate_gwf(G)

      call f_free(rhoexpo)
      call f_free_ptr(rhocoeff)

   end if

   ! Add the comparison between cuda hamiltonian and normal one if it is the case

   ! De-allocations
   call deallocate_Lzd_except_Glr(runObj%rst%KSwfn%Lzd)
   call deallocate_comms(runObj%rst%KSwfn%comms)
   call deallocate_orbs(runObj%rst%KSwfn%orbs)
   call free_DFT_PSP_projectors(nlpsp)

   !remove the directory which has been created if it is possible
   call deldir(runObj%inputs%dir_output,len(trim(runObj%inputs%dir_output)),ierror)

   call free_run_objects(runObj)
!   !finalize memory counting
!   call memocc(0,0,'count','stop')

   !Elapsed time
   !call cpu_time(tcpu1)
   !call system_clock(ncount1,ncount_rate,ncount_max)
   !tel=dble(ncount1-ncount0)/dble(ncount_rate)
   !write( *,'(1x,a,2(1x,f12.2))') 'CPU time/ELAPSED time ', tel, tcpu1-tcpu0

   if (.not. disable_deprecation) then
      call deprecation_message()
   end if

   call f_lib_finalize()

END PROGRAM memguess


!> Rotate the molecule via an orthogonal matrix in order to minimise the
!! volume of the cubic cell
subroutine optimise_volume(atoms,crmult,frmult,hx,hy,hz,rxyz)
   use module_base
   use module_types
   use module_interfaces, only: system_size
   implicit none
   type(atoms_data), intent(inout) :: atoms
   real(gp), intent(in) :: crmult,frmult
   real(gp), intent(inout) :: hx,hy,hz
   real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
   !local variables
   character(len=*), parameter :: subname='optimise_volume'
   integer :: iat,it,i
   real(gp) :: x,y,z,vol,tx,ty,tz,tvol,s,diag,dmax
   type(locreg_descriptors) :: Glr
   real(gp), dimension(3) :: shift
   real(gp), dimension(3,3) :: urot
   real(gp), dimension(:,:), allocatable :: txyz

   txyz = f_malloc((/ 3, atoms%astruct%nat /),id='txyz')
   call system_size(atoms,rxyz,crmult,frmult,hx,hy,hz,.false.,Glr,shift)
   !call volume(nat,rxyz,vol)
   vol=atoms%astruct%cell_dim(1)*atoms%astruct%cell_dim(2)*atoms%astruct%cell_dim(3)
   write(*,'(1x,a,1pe16.8)')'Initial volume (Bohr^3)',vol

   it=0
   diag=1.d-2 ! initial small diagonal element allows for search over all angles
   loop_rotations: do  ! loop over all trial rotations
      diag=diag*1.0001_gp ! increase diag to search over smaller angles
      it=it+1
      if (diag > 100._gp) exit loop_rotations ! smaller angle rotations do not make sense

      ! create a random orthogonal (rotation) matrix
      call random_number(urot)
      urot(:,:)=urot(:,:)-.5_gp
      do i=1,3
         urot(i,i)=urot(i,i)+diag
      enddo

      s=urot(1,1)**2+urot(2,1)**2+urot(3,1)**2
      s=1._gp/sqrt(s)
      urot(:,1)=s*urot(:,1) 

      s=urot(1,1)*urot(1,2)+urot(2,1)*urot(2,2)+urot(3,1)*urot(3,2)
      urot(:,2)=urot(:,2)-s*urot(:,1)
      s=urot(1,2)**2+urot(2,2)**2+urot(3,2)**2
      s=1._gp/sqrt(s)
      urot(:,2)=s*urot(:,2) 

      s=urot(1,1)*urot(1,3)+urot(2,1)*urot(2,3)+urot(3,1)*urot(3,3)
      urot(:,3)=urot(:,3)-s*urot(:,1)
      s=urot(1,2)*urot(1,3)+urot(2,2)*urot(2,3)+urot(3,2)*urot(3,3)
      urot(:,3)=urot(:,3)-s*urot(:,2)
      s=urot(1,3)**2+urot(2,3)**2+urot(3,3)**2
      s=1._gp/sqrt(s)
      urot(:,3)=s*urot(:,3) 

      ! eliminate reflections
      if (urot(1,1) <= 0._gp) urot(:,1)=-urot(:,1)
      if (urot(2,2) <= 0._gp) urot(:,2)=-urot(:,2)
      if (urot(3,3) <= 0._gp) urot(:,3)=-urot(:,3)

      ! apply the rotation to all atomic positions! 
      do iat=1,atoms%astruct%nat
         x=rxyz(1,iat) 
         y=rxyz(2,iat) 
         z=rxyz(3,iat)

         txyz(:,iat)=x*urot(:,1)+y*urot(:,2)+z*urot(:,3)
      enddo

      call system_size(atoms,txyz,crmult,frmult,hx,hy,hz,.false.,Glr,shift)
      tvol=atoms%astruct%cell_dim(1)*atoms%astruct%cell_dim(2)*atoms%astruct%cell_dim(3)
      !call volume(nat,txyz,tvol)
      if (tvol < vol) then
         write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol,it,diag
         rxyz(:,:)=txyz(:,:)
         vol=tvol
         dmax=max(atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3))
         ! if box longest along x switch x and z
         if (atoms%astruct%cell_dim(1) == dmax)  then
            do  iat=1,atoms%astruct%nat
               tx=rxyz(1,iat)
               tz=rxyz(3,iat)

               rxyz(1,iat)=tz
               rxyz(3,iat)=tx
            enddo
            ! if box longest along y switch y and z
         else if (atoms%astruct%cell_dim(2) == dmax .and. atoms%astruct%cell_dim(1) /= dmax)  then
            do  iat=1,atoms%astruct%nat
               ty=rxyz(2,iat) 
               tz=rxyz(3,iat)

               rxyz(2,iat)=tz 
               rxyz(3,iat)=ty
            enddo
         endif
      endif
   end do loop_rotations

   call f_free(txyz)

END SUBROUTINE optimise_volume


!> Add a shift in the periodic directions such that the system
!! uses as less as possible the modulo operation
subroutine shift_periodic_directions(at,rxyz,radii_cf)
   use module_base
   use module_types
   implicit none
   type(atoms_data), intent(inout) :: at
   real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
   real(gp), dimension(3,at%astruct%nat), intent(inout) :: rxyz
   !local variables
   character(len=*), parameter :: subname='shift_periodic_directions'
   integer :: iat,i,ityp
   real(gp) :: vol,tvol,maxsh,shiftx,shifty,shiftz
   real(gp), dimension(:,:), allocatable :: txyz

   !calculate maximum shift between these values
   !this is taken as five times the coarse radius around atoms
   maxsh=0.0_gp
   do ityp=1,at%astruct%ntypes
      maxsh=max(maxsh,5_gp*radii_cf(ityp,1))
   end do

   txyz = f_malloc((/ 3, at%astruct%nat /),id='txyz')

   call calc_vol(at%astruct%geocode,at%astruct%nat,rxyz,vol)

   if (at%astruct%geocode /= 'F') then
      loop_shiftx: do i=1,5000 ! loop over all trial rotations
         ! create a random orthogonal (rotation) matrix
         call random_number(shiftx)

         !apply the shift to all atomic positions taking into account the modulo operation
         do iat=1,at%astruct%nat
            txyz(1,iat)=modulo(rxyz(1,iat)+shiftx*maxsh,at%astruct%cell_dim(1))
            txyz(2,iat)=rxyz(2,iat)
            txyz(3,iat)=rxyz(3,iat)
         end do

         call calc_vol(at%astruct%geocode,at%astruct%nat,txyz,tvol)
         !print *,'vol',tvol

         if (tvol < vol) then
            write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
            rxyz(:,:)=txyz(:,:)
            vol=tvol
         endif
      end do loop_shiftx
   end if

   if (at%astruct%geocode == 'P') then
      loop_shifty: do i=1,5000 ! loop over all trial rotations
         ! create a random orthogonal (rotation) matrix
         call random_number(shifty)

         !apply the shift to all atomic positions taking into account the modulo operation
         do iat=1,at%astruct%nat
            txyz(1,iat)=rxyz(1,iat)
            txyz(2,iat)=modulo(rxyz(2,iat)+shifty*maxsh,at%astruct%cell_dim(2))
            txyz(3,iat)=rxyz(3,iat)
         end do

         call calc_vol(at%astruct%geocode,at%astruct%nat,txyz,tvol)

         if (tvol < vol) then
            write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
            rxyz(:,:)=txyz(:,:)
            vol=tvol
         endif
      end do loop_shifty
   end if

   if (at%astruct%geocode /= 'F') then
      loop_shiftz: do i=1,5000 ! loop over all trial rotations
         ! create a random orthogonal (rotation) matrix
         call random_number(shiftz)

         !apply the shift to all atomic positions taking into account the modulo operation
         do iat=1,at%astruct%nat
            txyz(1,iat)=rxyz(1,iat)
            txyz(2,iat)=rxyz(2,iat)
            txyz(3,iat)=modulo(rxyz(3,iat)+shiftz*maxsh,at%astruct%cell_dim(3))
         end do

         call calc_vol(at%astruct%geocode,at%astruct%nat,txyz,tvol)

         if (tvol < vol) then
            write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
            rxyz(:,:)=txyz(:,:)
            vol=tvol
         endif
      end do loop_shiftz
   end if

   call f_free(txyz)

END SUBROUTINE shift_periodic_directions


!> Calculate the extremes of the boxes taking into account the spheres around the atoms
subroutine calc_vol(geocode,nat,rxyz,vol)
   use module_base
   implicit none
   character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
   integer, intent(in) :: nat
   real(gp), dimension(3,nat), intent(in) :: rxyz
   real(gp), intent(out) :: vol
   !local variables
   integer :: iat
   real(gp) :: cxmin,cxmax,cymin,cymax,czmin,czmax

   cxmax=-1.e10_gp 
   cxmin=1.e10_gp

   cymax=-1.e10_gp 
   cymin=1.e10_gp

   czmax=-1.e10_gp 
   czmin=1.e10_gp

   do iat=1,nat
      cxmax=max(cxmax,rxyz(1,iat)) 
      cxmin=min(cxmin,rxyz(1,iat))

      cymax=max(cymax,rxyz(2,iat)) 
      cymin=min(cymin,rxyz(2,iat))

      czmax=max(czmax,rxyz(3,iat)) 
      czmin=min(czmin,rxyz(3,iat))
   enddo
   !print *,cxmax,cxmin,cymax,cymin,czmax,czmin
   !now calculate the volume for the periodic part
   if (geocode == 'P') then
      vol=(cxmax-cxmin)*(cymax-cymin)*(czmax-czmin)
   else if (geocode == 'S') then
      vol=(cxmax-cxmin)*(czmax-czmin)
   end if

END SUBROUTINE calc_vol


subroutine compare_cpu_gpu_hamiltonian(iproc,nproc,matacc,at,orbs,&
     nspin,ixc,ncong,Lzd,hx,hy,hz,rxyz,ntimes)
   use module_base
   use module_types
   use module_interfaces
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
   use gaussians, only: gaussian_basis, deallocate_gwf
   use module_xc

   implicit none
   integer, intent(in) :: iproc,nproc,nspin,ncong,ixc,ntimes
   real(gp), intent(in) :: hx,hy,hz
   type(material_acceleration), intent(in) :: matacc
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(inout) :: orbs
   type(local_zone_descriptors), intent(inout) :: Lzd
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   !local variables
   character(len=*), parameter :: subname='compare_cpu_gpu_hamiltonian'
   logical :: rsflag
   integer :: icoeff,i1,i2,i3,ispin,j
   integer :: iorb,n3d,n3p,n3pi,i3xcsh,i3s,jproc,nrhotot,nspinn,nvctrp
   integer(kind=8) :: itsc0,itsc1
   real(kind=4) :: tt
   real(gp) :: ttd,x,y,z,r2,arg,sigma2,ekin_sum,epot_sum,ekinGPU,epotGPU,gnrm,gnrm_zero,gnrmGPU
   real(gp) :: Rden,Rham,Rgemm,Rsyrk,Rprec,eSIC_DC
   real(kind=8) :: CPUtime,GPUtime
   type(gaussian_basis) :: G
   type(GPU_pointers) :: GPU
   type(xc_info) :: xc
   integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
   real(wp), dimension(:,:,:,:), allocatable :: pot,rho
   real(wp), dimension(:), pointer:: pottmp
   real(wp), dimension(:,:), allocatable :: gaucoeffs,psi,hpsi
   real(wp), dimension(:,:,:), allocatable :: overlap
   real(wp), dimension(:), pointer :: gbd_occ
   type(coulomb_operator) :: fake_pkernelSIC
   type(confpot_data), dimension(orbs%norbp) :: confdatarr

   call default_confinement_data(confdatarr,orbs%norbp)

   !nullify pkernelSIC pointer
   nullify(fake_pkernelSIC%kernel)

   !nullify the G%rxyz pointer
   nullify(G%rxyz)
   !extract the gaussian basis from the pseudowavefunctions
   call gaussian_pswf_basis(21,.false.,iproc,nspin,at,rxyz,G,gbd_occ)

   gaucoeffs = f_malloc((/ G%ncoeff, orbs%norbp*orbs%nspinor /),id='gaucoeffs')

   !fill randomly the gaussian coefficients for the orbitals considered
   do iorb=1,orbs%norbp*orbs%nspinor
      do icoeff=1,G%ncoeff
         call random_number(tt)
         gaucoeffs(icoeff,iorb)=real(tt,wp)
      end do
   end do

   !allocate the wavefunctions
   psi = f_malloc0((/ Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f, orbs%nspinor*orbs%norbp /),id='psi')
   hpsi = f_malloc0((/ Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f , orbs%nspinor*orbs%norbp /),id='hpsi')

   !call to_zero(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,psi)
   !call to_zero(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,hpsi)

   !convert the gaussians in wavelets
   call gaussians_to_wavelets(iproc,nproc,at%astruct%geocode,orbs,Lzd%Glr%d,&
           hx,hy,hz,Lzd%Glr%wfd,G,gaucoeffs,psi)

   call f_free(gaucoeffs)
   call f_free_ptr(gbd_occ)

   !deallocate the gaussian basis descriptors
   call deallocate_gwf(G)

   !allocate and initialise the potential and the density
   pot = f_malloc((/ Lzd%Glr%d%n1i, Lzd%Glr%d%n2i, Lzd%Glr%d%n3i, nspin /),id='pot')
   rho = f_malloc((/ Lzd%Glr%d%n1i, Lzd%Glr%d%n2i, Lzd%Glr%d%n3i, nspin /),id='rho')

   !here the potential can be used for building the density
   nscatterarr = f_malloc((/ 0.to.nproc-1, 1.to.4 /),id='nscatterarr')
   ngatherarr = f_malloc((/ 0.to.nproc-1, 1.to.2 /),id='ngatherarr')

   if (ixc < 0) then
      call xc_init(xc, ixc, XC_MIXED, nspin)
   else
      call xc_init(xc, ixc, XC_ABINIT, nspin)
   end if

   !normally nproc=1
   do jproc=0,nproc-1
      call PS_dim4allocation(at%astruct%geocode,'D',jproc,nproc,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,xc_isgga(xc),(ixc/=13),&
         &   n3d,n3p,n3pi,i3xcsh,i3s)
      nscatterarr(jproc,1)=n3d
      nscatterarr(jproc,2)=n3p
      nscatterarr(jproc,3)=i3s+i3xcsh-1
      nscatterarr(jproc,4)=i3xcsh
   end do

   ngatherarr(:,1)=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(:,2)
   ngatherarr(:,2)=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(:,3)

   !components of the charge density
   if (orbs%nspinor ==4) then
      nspinn=4
   else
      nspinn=nspin
   end if

   !flag for toggling the REDUCE_SCATTER stategy
   rsflag = .not.xc_isgga(xc)

   !calculate dimensions of the complete array to be allocated before the reduction procedure
   if (rsflag) then
      nrhotot=0
      do jproc=0,nproc-1
         nrhotot=nrhotot+nscatterarr(jproc,1)
      end do
   else
      nrhotot=Lzd%Glr%d%n3i
   end if

   call local_potential_dimensions(iproc,Lzd,orbs,xc,ngatherarr(0,1))

   !allocate the necessary objects on the GPU
   !set initialisation of GPU part 
   !initialise the acceleration strategy if required
   call init_material_acceleration(iproc,matacc,GPU)

   !allocate arrays for the GPU if a card is present
   if (GPU%OCLconv) then
      !the same with OpenCL, but they cannot exist at same time
      call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,Lzd%Glr%geocode,&
           nspin,Lzd%Glr%wfd,orbs,GPU)
   end if
   if (iproc == 0) write(*,*)&
      &   'GPU data allocated'

   write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Density calculation'

   !for each of the orbitals treated by the processor build the partial densities
   !call cpu_time(t0)
   call nanosec(itsc0)
   do j=1,ntimes
      call tenminustwenty(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot*nspinn,pot,nproc)
      call local_partial_density(nproc,rsflag,nscatterarr,&
           nrhotot,Lzd%Glr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,orbs,&
           psi,pot)
   end do
   call nanosec(itsc1)
   !call cpu_time(t1)
   !CPUtime=real(t1-t0,kind=8)
   CPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   !now the GPU part
   !for each of the orbitals treated by the processor build the partial densities
   !call cpu_time(t0)
   call nanosec(itsc0)
   do j=1,ntimes
      !switch between GPU/CPU treatment of the density
      if (GPU%OCLconv) then
         call local_partial_density_OCL(orbs,nrhotot,Lzd%Glr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,psi,rho,GPU)
      end if
   end do
   call nanosec(itsc1)
   !call cpu_time(t1)
   !GPUtime=real(t1-t0,kind=8)
   GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   call f_free(nscatterarr)
   call f_free(ngatherarr)


   !compare the results between the different actions of the hamiltonian
   !check the differences between the results
   call compare_data_and_gflops(CPUtime,GPUtime,&
        & real(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,kind=8)*192.d0,pot,rho,&
        Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,ntimes*orbs%norbp,.false.,Rden)

   call f_free(rho)


   !here the grid spacings are the small ones
   sigma2=0.125_gp*((Lzd%Glr%d%n1i*hx)**2+(Lzd%Glr%d%n2i*hy)**2+(Lzd%Glr%d%n3i*hz)**2)
   do ispin=1,nspin
      do i3=1,Lzd%Glr%d%n3i
         z=hz*real(i3-Lzd%Glr%d%n3i/2-1,gp)
         do i2=1,Lzd%Glr%d%n2i
            y=hy*real(i2-Lzd%Glr%d%n2i/2-1,gp)
            do i1=1,Lzd%Glr%d%n1i
               x=hx*real(i1-Lzd%Glr%d%n1i/2-1,gp)
               !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
               r2=x**2+y**2+z**2
               arg=0.5d0*r2/sigma2
               ttd=dexp(-arg)

               pot(i1,i2,i3,ispin)=ttd
            end do
         end do
      end do
   end do

   write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Local Hamiltonian calculation'

   !warm-up
   !call local_hamiltonian(iproc,orbs,Lzd%Glr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum) 

   !apply the CPU hamiltonian
   !take timings
   call nanosec(itsc0)
   xc%ixc = 0
   do j=1,ntimes
      pottmp = f_malloc_ptr(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*(nspin),id='pottmp')
      call vcopy(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*(nspin),pot(1,1,1,1),1,pottmp(1),1)
      call local_hamiltonian(iproc,nproc,orbs%npsidim_orbs,orbs,Lzd,hx,hy,hz,0,confdatarr,pottmp,psi,hpsi, &
           fake_pkernelSIC,xc,0.0_gp,ekin_sum,epot_sum,eSIC_DC)
      call f_free_ptr(pottmp)
   end do
   xc%ixc = ixc
   call nanosec(itsc1)
   CPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   print *,'ekin,epot=',ekin_sum,epot_sum

   !WARNING: local hamiltonian overwrites the psis
   !warm-up
   !call gpu_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,hx,hy,hz,orbs,GPU,ekinGPU,epotGPU)

   !apply the GPU hamiltonian and put the results in the hpsi_GPU array
   GPU%hpsi_ASYNC = f_malloc_ptr((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp,id='GPU%hpsi_ASYNC')

   !take timings
   call nanosec(itsc0)
   do j=1,ntimes
      if (GPU%OCLconv) then
         call local_hamiltonian_OCL(orbs,Lzd%Glr,hx,hy,hz,orbs%nspin,pot,psi,GPU%hpsi_ASYNC,ekinGPU,epotGPU,GPU)
      end if
   end do
   if(ASYNCconv .and. GPU%OCLconv) call finish_hamiltonian_OCL(orbs,ekinGPU,epotGPU,GPU)
   call nanosec(itsc1)
   GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   print *,'ekinGPU,epotGPU',ekinGPU,epotGPU

   !compare the results between the different actions of the hamiltonian
   !check the differences between the results
   call compare_data_and_gflops(CPUtime,GPUtime,&
      &   real(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,kind=8)*real(192+46*3+192+2,kind=8),hpsi,GPU%hpsi_ASYNC,&
   orbs%norbp*orbs%nspinor*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),ntimes*orbs%norbp,.false.,Rham)

   call f_free(pot)

   write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Linear Algebra (Blas)'

   !perform the scalar product between the hpsi wavefunctions
   !actually this is <hpsi|hpsi> it has no meaning.
   !this works only if nspinor==1
   overlap = f_malloc((/ orbs%norbp, orbs%norbp, 2 /),id='overlap')

   call nanosec(itsc0)
   do j=1,ntimes
      call DGEMM('T','N',orbs%norbp,orbs%norbp,(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),1.0_wp,&
         &   psi(1,1),(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),&
      hpsi(1,1),(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),0.0_wp,&
         &   overlap(1,1,1),orbs%norbp)
   end do
   call nanosec(itsc1)
   CPUtime=real(itsc1-itsc0,kind=8)*1.d-9


   call nanosec(itsc0)
   do j=1,ntimes
      call GEMMSY('T','N',orbs%norbp,orbs%norbp,(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),1.0_wp,&
         &   psi(1,1),(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),&
      hpsi(1,1),(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),0.0_wp,&
         &   overlap(1,1,2),orbs%norbp)
   end do
   call nanosec(itsc1)
   GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   !comparison between the results
   call compare_data_and_gflops(CPUtime,GPUtime,&
      &   real(orbs%norbp**2,kind=8)*real((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*2,kind=8),overlap(1,1,1),overlap(1,1,2),&
   orbs%norbp**2,ntimes,.false.,Rgemm)


   nvctrp=Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f

   call nanosec(itsc0)
   do j=1,ntimes
      call dsyrk('L','T',orbs%norbp,nvctrp,1.0_wp,psi(1,1),nvctrp,0.0_wp,&
         &   overlap(1,1,1),orbs%norbp)
   end do
   call nanosec(itsc1)
   CPUtime=real(itsc1-itsc0,kind=8)*1.d-9


   call nanosec(itsc0)
   do j=1,ntimes
      call syrk('L','T',orbs%norbp,nvctrp,1.0_wp,psi(1,1),nvctrp,0.0_wp,&
         &   overlap(1,1,2),orbs%norbp)
   end do
   call nanosec(itsc1)
   GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   call compare_data_and_gflops(CPUtime,GPUtime,&
        real(orbs%norbp*(orbs%norbp+1),kind=8)*&
        real(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,kind=8),overlap(1,1,1),overlap(1,1,2),&
        orbs%norbp**2,ntimes,.false.,Rsyrk)

   call f_free(overlap)


   !-------------------now the same for preconditioning
   write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Preconditioner'

   !the input function is psi
   call nanosec(itsc0)
   do j=1,ntimes
      call preconditionall(orbs,Lzd%Glr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
   end do
   call nanosec(itsc1)

   CPUtime=real(itsc1-itsc0,kind=8)*1.d-9
   print *,'gnrm',gnrm


   !GPU data are aLzd%Glready on the card, must be only copied back
   !the input function is GPU%hpsi in that case
   call nanosec(itsc0)
   do j=1,ntimes
      !Preconditions all orbitals belonging to iproc
      !and calculate the partial norm of the residue
      !switch between CPU and GPU treatment
      if (GPU%OCLconv) then
         call preconditionall_OCL(orbs,Lzd%Glr,hx,hy,hz,ncong,&
            &   GPU%hpsi_ASYNC,gnrmGPU,gnrm_zero,GPU)
      end if
   end do
   call nanosec(itsc1)

   GPUtime=real(itsc1-itsc0,kind=8)*1.d-9
   print *,'gnrmGPU',gnrmGPU

   call compare_data_and_gflops(CPUtime,GPUtime,&
      &   real(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,kind=8)*real((192+46*3+192+2-1+12)*(ncong+1),kind=8),hpsi,GPU%hpsi_ASYNC,&
   orbs%norbp*orbs%nspinor*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),ntimes*orbs%norbp,.false.,Rprec)

   call f_free_ptr(GPU%hpsi_ASYNC)
   call f_free(psi)
   call f_free(hpsi)


   !free the card at the end
   if (GPU%OCLconv) then
      call free_gpu_OCL(GPU,orbs,nspin)
   end if

   call xc_end(xc)

   !finalise the material accelearion usage
   call release_material_acceleration(GPU)


   write(*,'(1x,a,5(1x,f7.3))')'Ratios:',Rden,Rham,Rgemm,Rsyrk,Rprec

END SUBROUTINE compare_cpu_gpu_hamiltonian


subroutine compare_data_and_gflops(CPUtime,GPUtime,GFlopsfactor,&
      &   CPUdata,GPUdata,n,ntimes,dowrite,ratio)
   use module_base
   implicit none
   logical, intent(in) :: dowrite
   integer, intent(in) :: n,ntimes
   real(gp), intent(in) :: CPUtime,GPUtime,GFlopsfactor
   real(gp), intent(out) :: ratio
   real(wp), dimension(n), intent(in) :: CPUdata,GPUdata
   !local variables
   integer :: i
   real(gp) :: CPUGflops,GPUGflops,maxdiff,comp,threshold

   threshold=1.d-12
   !un-initialize valies which might suffer from fpe
   GPUGflops=-1.0_gp
   CPUGflops=-1.0_gp
   ratio=-1.0_gp

   if (CPUtime > 0.0_gp) CPUGflops=GFlopsfactor*real(ntimes,gp)/(CPUtime*1.d9)
   if (GPUtime > 0.0_gp) GPUGflops=GFlopsfactor*real(ntimes,gp)/(GPUtime*1.d9)

   maxdiff=0.0_gp

   rewind(17)

   do i=1,n
      if (dowrite) write(17,'(i6,2(1pe24.17))')i,CPUdata(i),GPUdata(i)
      comp=abs(CPUdata(i)-GPUdata(i))
      maxdiff=max(maxdiff,comp)
   end do
   if (GPUtime > 0.0_gp) ratio=CPUtime/GPUtime
   write(*,'(1x,a)')'| CPU: ms  |  Gflops  || GPU:  ms |  GFlops  || Ratio  | No. Elements | Max. Diff. |'

   write(*,'(1x,2(2(a,f10.2),a),a,f8.3,a,i14,a,1pe12.4,a)',advance='no')&
      &   '|',CPUtime*1.d3/real(ntimes,kind=8),'|',& ! Time CPU (ms)
      &   CPUGflops,'|',& !Gflops CPU (ms)
      &   '|',GPUtime*1.d3/real(ntimes,kind=8),'|',& ! Time GPU (ms)
      &   GPUGflops,'|',&!Gflops GPU (ms)
      &   '|',ratio,'|',& ! ratio
      &   n,'|',& !No. elements
      &   maxdiff,'|' ! maxdiff
   if (maxdiff <= threshold) then
      write(*,'(a)')''
   else
      write(*,'(a)')'<<<< WARNING' 
   end if

END SUBROUTINE compare_data_and_gflops


!> Extract the compressed wavefunction from the given file 
subroutine take_psi_from_file(filename,in_frag,hx,hy,hz,lr,at,rxyz,orbs,psi,iorbp,ispinor,ref_frags)
   use module_base
   use module_types
   use module_interfaces
   use module_fragments
   use locreg_operations, only: lpsi_to_global2
   implicit none
   integer, intent(inout) :: iorbp, ispinor
   real(gp), intent(in) :: hx,hy,hz
   character(len=*), intent(in) :: filename
   type(locreg_descriptors), intent(in) :: lr
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs
   type(fragmentInputParameters), intent(in) :: in_frag
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor), intent(out) :: psi
   type(system_fragment), dimension(in_frag%nfrag_ref), intent(inout) :: ref_frags
   !local variables
   character(len=*), parameter :: subname='take_psi_form_file'
   logical :: perx,pery,perz
   integer :: nb1,nb2,nb3,ikpt, ispin, i
   integer :: wave_format_from_filename,iformat
   real(gp) :: eval_fake
   real(wp), dimension(:,:,:), allocatable :: psifscf
   real(gp), dimension(:,:), allocatable :: rxyz_file
   character(len = 1) :: code

   integer :: confPotOrder !lr408
   real(gp) :: locrad, confPotprefac !lr408
   real(gp), dimension(3) :: locregCenter !lr408
   character(len=3) :: in_name !lr408
   type(local_zone_descriptors) :: Lzd 
   integer, dimension(1) :: orblist
   character(len=100) :: filename_start
   real(wp), allocatable, dimension(:) :: lpsi
   type(orbitals_data) :: lin_orbs

   rxyz_file = f_malloc((/ at%astruct%nat, 3 /),id='rxyz_file')

   iformat = wave_format_from_filename(0, filename)
   if (iformat == WF_FORMAT_PLAIN .or. iformat == WF_FORMAT_BINARY) then
      !conditions for periodicity in the three directions
      perx=(at%astruct%geocode /= 'F')
      pery=(at%astruct%geocode == 'P')
      perz=(at%astruct%geocode /= 'F')

      !buffers related to periodicity
      !WARNING: the boundary conditions are not assumed to change between new and old
      call ext_buffers_coarse(perx,nb1)
      call ext_buffers_coarse(pery,nb2)
      call ext_buffers_coarse(perz,nb3)

      psifscf = f_malloc((/ -nb1.to.2*lr%d%n1+1+nb1, -nb2.to.2*lr%d%n2+1+nb2, -nb3.to.2*lr%d%n3+1+nb3 /),id='psifscf')

      !find the value of iorbp
      read(filename(index(filename, ".", back = .true.)+2:len(filename)),*) iorbp
      i = index(filename, "-k", back = .true.)+2
      read(filename(i:i+2),*) ikpt
      i = index(filename, "-", back = .true.)+1
      read(filename(i:i),*) code
      if (code == "U" .or. code == "N") ispin = 1
      if (code == "D") ispin = 2
      read(filename(i+1:i+1),*) code
      if (code == "R") ispinor = 1
      if (code == "I") ispinor = 2



      i = index(filename, "/",back=.true.)+1
      read(filename(i:i+3),*) in_name ! lr408

      if (in_name == 'min') then
         ! Create orbs data structure.
         call nullify_orbitals_data(lin_orbs)
         call copy_orbitals_data(orbs, lin_orbs, subname)

         lin_orbs%norb = 1
         lin_orbs%norbp = 1

         ! need to change the lr info so relates to locregs not global
         ! need to copy Glr and hgrids into Lzd
         Lzd%Glr = lr
         Lzd%hgrids(1) = hx
         Lzd%hgrids(2) = hy
         Lzd%hgrids(3) = hz
         orblist = iorbp

         i = index(filename, "-",back=.true.)+1
         read(filename(1:i),*) filename_start

         print*,'Initialize linear'
         call initialize_linear_from_file(0,1,in_frag,at%astruct,rxyz,lin_orbs,Lzd,&
              WF_FORMAT_BINARY,filename_start//"/","minBasis",ref_frags,orblist)

         filename_start = trim(filename_start)//"/minBasis"

         lpsi = f_malloc(1.to.Lzd%llr(1)%wfd%nvctr_c+7*Lzd%llr(1)%wfd%nvctr_f,id='lpsi')
      end if

      if (iformat == WF_FORMAT_BINARY) then
         open(unit=99,file=trim(filename),status='unknown',form="unformatted")
      else
         open(unit=99,file=trim(filename),status='unknown')
      end if

      !@todo geocode should be passed in the localisation regions descriptors
      if (in_name /= 'min') then
         call readonewave(99, (iformat == WF_FORMAT_PLAIN),iorbp,0,lr%d%n1,lr%d%n2,lr%d%n3, &
              & hx,hy,hz,at,lr%wfd,rxyz_file,rxyz,psi(1,ispinor),eval_fake,psifscf)
      else
         call readonewave_linear(99, (iformat == WF_FORMAT_PLAIN),iorbp,0,&
              lr%d%n1,lr%d%n2,lr%d%n3,hx,hy,hz,at,Lzd%llr(1),rxyz_file,rxyz,&
              locrad,locregCenter,confPotOrder,confPotPrefac,&
              lpsi(1),eval_fake,psifscf)

         call f_zero(psi)

         call Lpsi_to_global2(0,Lzd%llr(1)%wfd%nvctr_c+7*Lzd%llr(1)%wfd%nvctr_f, &
              lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,1,1,1,lr,Lzd%Llr(1),lpsi,psi)

         call f_free(lpsi)
      end if

      ! Update iorbp
      iorbp = (ikpt - 1) * orbs%norb + (ispin - 1) * orbs%norbu + iorbp

      close(99)

      call f_free(psifscf)

   else if (iformat == WF_FORMAT_ETSF) then
      call read_one_wave_etsf(0,filename,iorbp,0,orbs%nspinor,lr%d%n1,lr%d%n2,lr%d%n3,&
           & hx,hy,hz,at,rxyz_file,rxyz,lr%wfd,psi,eval_fake)
   end if
   call f_free(rxyz_file)
END SUBROUTINE take_psi_from_file


!> Deprecated message for memguess (do not use directly!!)
subroutine deprecation_message()
   implicit none
   write(*, "(15x,A)") "+--------------------------------------------+"
   write(*, "(15x,A)") "|                                            |"
   write(*, "(15x,A)") "| /!\ memguess is deprecated since 1.6.0 /!\ |"
   write(*, "(15x,A)") "|                                            |"
   write(*, "(15x,A)") "|     Use bigdft-tool  instead,  located     |"
   write(*, "(15x,A)") "|     in the  build directory or in  the     |"
   write(*, "(15x,A)") "|     bin directory of the install path.     |"
   write(*, "(15x,A)") "|       $ bigdft-tool -h for help            |"
   write(*, "(15x,A)") "|                                            |"
   write(*, "(15x,A)") "+--------------------------------------------+"
END SUBROUTINE deprecation_message
