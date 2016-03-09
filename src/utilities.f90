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
   use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, assignment(=), &
                                SPARSE_FULL, DENSE_FULL, DENSE_PARALLEL, &
                                sparsematrix_malloc_ptr, deallocate_sparse_matrix, deallocate_matrices, &
                                sparse_matrix_metadata, deallocate_sparse_matrix_metadata
   use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple
   use sparsematrix_io, only: read_sparse_matrix, write_sparse_matrix
   use sparsematrix, only: uncompress_matrix, uncompress_matrix_distributed2
   use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft, &
                                     sparse_matrix_and_matrices_init_from_file_ccs, &
                                     sparse_matrix_metadata_init_from_file, &
                                     ccs_data_from_sparse_matrix, &
                                     ccs_matrix_write
   use postprocessing_linear, only: CHARGE_ANALYSIS_LOEWDIN, CHARGE_ANALYSIS_MULLIKEN, &
                                    CHARGE_ANALYSIS_PROJECTOR, &
                                    loewdin_charge_analysis_core
   use multipole, only: projector_for_charge_analysis, multipole_analysis_driver
   use io, only: write_linear_coefficients, read_linear_coefficients
   use bigdft_run, only: bigdft_init
   implicit none
   external :: gather_timings
   character(len=*), parameter :: subname='utilities'
   character(len=1) :: geocode
   character(len=3) :: do_ortho
   character(len=30) :: tatonam, radical, colorname, linestart, lineend, cname, methodc
   character(len=128) :: method_name, overlap_file, hamiltonian_file, kernel_file, coeff_file, pdos_file
   character(len=128) :: line, cc, output_pdos, conversion, infile, outfile
   logical :: charge_analysis = .false.
   logical :: solve_eigensystem = .false.
   logical :: calculate_pdos = .false.
   logical :: convert_matrix_format = .false.
   type(atoms_data) :: at
   type(sparse_matrix_metadata) :: smmd
   integer :: istat, i_arg, ierr, nspin, icount, nthread, method, ntypes
   integer :: nfvctr_s, nseg_s, nvctr_s, nfvctrp_s, isfvctr_s
   integer :: nfvctr_m, nseg_m, nvctr_m, nfvctrp_m, isfvctr_m
   integer :: nfvctr_l, nseg_l, nvctr_l, nfvctrp_l, isfvctr_l
   integer :: iconv
   integer,dimension(:),pointer :: on_which_atom
   integer,dimension(:),pointer :: keyv_s, keyv_m, keyv_l, on_which_atom_s, on_which_atom_m, on_which_atom_l
   integer,dimension(:),pointer :: iatype, nzatom, nelpsp
   integer,dimension(:),pointer :: col_ptr, row_ind
   integer,dimension(:,:,:),pointer :: keyg_s, keyg_m, keyg_l
   real(kind=8),dimension(:),pointer :: matrix_compr, eval_ptr
   real(kind=8),dimension(:,:),pointer :: rxyz, coeff_ptr
   real(kind=8),dimension(:),allocatable :: eval, energy_arr, occups
   real(kind=8),dimension(:,:),allocatable :: denskernel, pdos, occup_arr
   logical,dimension(:,:),allocatable :: calc_array
   type(matrices) :: ovrlp_mat, hamiltonian_mat, kernel_mat, mat
   type(sparse_matrix) :: smat_s, smat_m, smat_l, smat
   type(dictionary), pointer :: dict_timing_info
   integer :: iunit, nat, iat, iat_prev, ii, iitype, iorb, itmb, itype, ival, ios, ipdos, ispin
   integer :: jtmb, norbks, npdos, npt, ntmb, jjtmb
   character(len=20),dimension(:),pointer :: atomnames
   character(len=30),dimension(:),allocatable :: pdos_name
   real(kind=8),dimension(3) :: cell_dim
   character(len=2) :: backslash, num
   real(kind=8) :: energy, occup, occup_pdos, total_occup
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
   call bigdft_init()

    if (bigdft_mpi%iproc==0) then
        call yaml_new_document()
        call print_logo()
    end if


   !Time initialization
   call f_timing_reset(filename='time.yaml',master=(bigdft_mpi%iproc==0),verbose_mode=.false.)


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
           & 'perform a charge analysis (Loewdin or Mulliken)'

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
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = hamiltonian_file)
            !write(*,'(1x,2a)')&
            !   &   'perform a Loewdin charge analysis'
            charge_analysis = .true.
            exit loop_getargs
         else if (trim(tatonam)=='solve-eigensystem') then
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
         end if
         i_arg = i_arg + 1
      end do loop_getargs
   end if


   if (charge_analysis) then
       if (bigdft_mpi%iproc==0) then
           call yaml_comment('Charge analysis',hfill='-')
       end if
       
       ! Determine the method
       select case(trim(method_name))
       case ('loewdin','LOEWDIN')
           method = CHARGE_ANALYSIS_LOEWDIN
       case ('mulliken','MULLIKEN')
           method = CHARGE_ANALYSIS_MULLIKEN
       case ('projector','PROJECTOR')
           method = CHARGE_ANALYSIS_PROJECTOR
       case default
           call f_err_throw('Unknown Method for the charge analysis',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       end select

       call sparse_matrix_metadata_init_from_file('sparsematrix_metadata.bin', smmd)
       if (bigdft_mpi%iproc==0) then
           call yaml_mapping_open('Atomic System Properties')
           call yaml_map('Types of atoms',smmd%atomnames)
           call yaml_mapping_close()
       end if

       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s, ovrlp_mat, &
            init_matmul=.true.)

       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(kernel_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, smat_l, kernel_mat, &
            init_matmul=.true.)

       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, smat_m, hamiltonian_mat, &
            init_matmul=.true.)

       call timing(bigdft_mpi%mpi_comm,'INIT','PR')

       select case(method)
       case (CHARGE_ANALYSIS_MULLIKEN)
           methodc='loewdin'
           do_ortho='no'
       case (CHARGE_ANALYSIS_LOEWDIN)
           methodc='loewdin'
           do_ortho='yes'
       case (CHARGE_ANALYSIS_PROJECTOR)
           methodc='projector'
           do_ortho='no'
       case default
           call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
       end select

       call multipole_analysis_driver(bigdft_mpi%iproc, bigdft_mpi%nproc, 0, 11, &
            smmd, smat_s, smat_m, smat_l, &
            ovrlp_mat, hamiltonian_mat, kernel_mat, smmd%rxyz, &
            methodc, do_ortho=trim(do_ortho), projectormode='simple', &
            calculate_multipole_matrices=.false., do_check=.false., &
            multipole_matrix_in=(/(/ovrlp_mat/)/))

       call timing(bigdft_mpi%mpi_comm,'CALC','PR')

       call deallocate_sparse_matrix_metadata(smmd)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_l)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(kernel_mat)
       call deallocate_sparse_matrix(smat_m)
       call deallocate_matrices(hamiltonian_mat)

       if (bigdft_mpi%iproc==0) then
           call yaml_comment('done',hfill='-')
       end if

   end if


   if (solve_eigensystem) then

       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s, ovrlp_mat, &
            init_matmul=.false.)!, nat=nat, rxyz=rxyz, iatype=iatype, ntypes=ntypes, &
            !nzatom=nzatom, nelpsp=nelpsp, atomnames=atomnames)
       call sparse_matrix_metadata_init_from_file('sparsematrix_metadata.bin', smmd)
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, smat_m, hamiltonian_mat, &
            init_matmul=.false.)

       ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='ovrlp_mat%matrix')
       call uncompress_matrix(bigdft_mpi%iproc, smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
       hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='hamiltonian_mat%matrix')
       call uncompress_matrix(bigdft_mpi%iproc, smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
       eval = f_malloc(smat_s%nfvctr,id='eval')

       if (bigdft_mpi%iproc==0) then
           call yaml_comment('Diagonalizing the matrix',hfill='~')
       end if
       call diagonalizeHamiltonian2(bigdft_mpi%iproc, smat_s%nfvctr, hamiltonian_mat%matrix, ovrlp_mat%matrix, eval)
       if (bigdft_mpi%iproc==0) then
           call yaml_comment('Matrix succesfully diagonalized',hfill='~')
       end if
       iunit=99
       call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
       !call writeLinearCoefficients(iunit, .true., nat, rxyz, smat_s%nfvctr, smat_s%nfvctr, &
       !     smat_s%nfvctr, hamiltonian_mat%matrix, eval)
       call write_linear_coefficients(0, trim(coeff_file), nat, rxyz, iatype, ntypes, nzatom, &
            nelpsp, atomnames, smat_s%nfvctr, smat_s%nfvctr, smat_s%nspin, hamiltonian_mat%matrix, eval)
       call f_close(iunit)

       call f_free(eval)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(hamiltonian_mat)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_m)
       call f_free_ptr(rxyz)
       call f_free_ptr(iatype)
       call f_free_ptr(nzatom)
       call f_free_ptr(nelpsp)
       call f_free_str_ptr(len(atomnames),atomnames)
   end if


   if (calculate_pdos) then
       iunit = 99
       if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(coeff_file),hfill='~')
       call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
       call read_linear_coefficients(trim(coeff_file), nspin, ntmb, norbks, coeff_ptr, &
            eval=eval_ptr)
       call f_close(iunit)

       if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(overlap_file),hfill='~')
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(overlap_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s, ovrlp_mat, &
            init_matmul=.false.)!, iatype=iatype, ntypes=ntypes, atomnames=atomnames, &
            !on_which_atom=on_which_atom)
       call sparse_matrix_metadata_init_from_file('sparsematrix_metadata.bin', smmd)
       if (ntmb/=smat_s%nfvctr) call f_err_throw('ntmb/=smat_s%nfvctr')
       ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_PARALLEL, id='ovrlp_mat%matrix')
       !call uncompress_matrix(bigdft_mpi%iproc, smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
       call uncompress_matrix_distributed2(bigdft_mpi%iproc, smat_s, DENSE_PARALLEL, &
            ovrlp_mat%matrix_compr, ovrlp_mat%matrix(1:,1:,1))

       if (bigdft_mpi%iproc==0) call yaml_comment('Reading from file '//trim(hamiltonian_file),hfill='~')
       call sparse_matrix_and_matrices_init_from_file_bigdft(trim(hamiltonian_file), &
            bigdft_mpi%iproc, bigdft_mpi%nproc, smat_m, hamiltonian_mat, &
            init_matmul=.false.)
       hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_PARALLEL, id='hamiltonian_mat%matrix')
       !call uncompress_matrix(bigdft_mpi%iproc, smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
       call uncompress_matrix_distributed2(bigdft_mpi%iproc, smat_m, DENSE_PARALLEL, &
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
           do 
               !read(iunit01,*,iostat=ios) cc, ival
               read(iunit,'(a128)',iostat=ios) line
               if (ios/=0) exit
               !write(*,*) 'line',line
               read(line,*,iostat=ios) cc, cname
               if (cc=='#') then
                   ipdos = ipdos + 1
                   pdos_name(ipdos) = trim(cname)
                   cycle 
               end if
               write(*,*) 'ipdos, line', ipdos, line
               read(line,*,iostat=ios) cc, ival
               do itype=1,ntypes
                   if (trim(atomnames(itype))==trim(cc)) then
                       iitype = itype
                       exit
                   end if
               end do
               iat_prev = -1
               do itmb=1,ntmb
                   iat = on_which_atom(itmb)
                   if (iat/=iat_prev) then
                       ii = 0
                   end if
                   iat_prev = iat
                   itype = iatype(iat)
                   ii = ii + 1
                   if (itype==iitype .and. ii==ival) then
                       if (calc_array(itmb,ipdos)) stop 'calc_array(itmb)'
                       calc_array(itmb,ipdos) = .true.
                   end if
               end do
           end do

       !end do npdos_loop
       call f_close(iunit)

       energy_arr = f_malloc0(norbks,id='energy_arr')
       occup_arr = f_malloc0((/npdos,norbks/),id='occup_arr')
       occups = f_malloc0(npdos,id='occups')

       denskernel = f_malloc((/ntmb,smat_s%nfvctrp/),id='denskernel')
       pdos = f_malloc0((/npt,npdos/),id='pdos')
       ! Calculate a partial kernel for each KS orbital
       !do ipdos=1,npdos
           !if (bigdft_mpi%iproc==0) call yaml_map('PDoS number',ipdos)
           !if (bigdft_mpi%iproc==0) call yaml_map('start, increment',(/ipdos,npdos/))
           !write(num,fmt='(i2.2)') ipdos
           if (bigdft_mpi%iproc==0) bar=f_progress_bar_new(nstep=norbks)
           do iorb=1,norbks
               if (bigdft_mpi%iproc==0) call yaml_map('orbital being processed',iorb)
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
               !call mpiallred(energy, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
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
               if (bigdft_mpi%iproc==0) call dump_progress_bar(bar,step=iorb)
           end do
           !total_occup = total_occup + occup_pdos
           !!if (ipdos==1) then
           !!    write(iunit,'(a,i0,a)') "plot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           !!else
           !!    write(iunit,'(a,i0,a)') "replot f",ipdos,"(x) lc rgb 'color' lt 1 lw 2 w l title 'name'"
           !!end if
           !!if (bigdft_mpi%iproc==0) call yaml_map('sum of PDoS',occup_pdos)
           !!output_pdos='PDoS_'//num//'.dat'
           !!call yaml_map('output file',trim(output_pdos))
           !!call f_open_file(iunit01, file=trim(output_pdos), binary=.false.)
           !!write(iunit01,'(a)') '#             energy                pdos'
           !!do ipt=1,npt
           !!    write(iunit01,'(2es20.12)') energy_bins(1,ipt), pdos(ipt,ipdos)
           !!end do
           !!call f_close(iunit01)
       !end do
       call mpiallred(occup_arr, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(energy_arr, mpi_sum, comm=bigdft_mpi%mpi_comm)

       if (bigdft_mpi%iproc==0) then
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
                   lineend = ' ,\'
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
       if (bigdft_mpi%iproc==0) call yaml_map('sum of total DoS',total_occup)
       call mpibarrier(bigdft_mpi%mpi_comm)

       call f_free(pdos)
       call f_free(denskernel)
       call f_free_ptr(eval_ptr)
       call deallocate_matrices(ovrlp_mat)
       call deallocate_matrices(hamiltonian_mat)
       call deallocate_sparse_matrix(smat_s)
       call deallocate_sparse_matrix(smat_m)
       call f_free_ptr(iatype)
       call f_free_str_ptr(len(atomnames),atomnames)
       call f_free_str(len(pdos_name),pdos_name)
       call f_free(calc_array)
       call f_free_ptr(on_which_atom)
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
       case default
           call f_err_throw("wrong value for conversion; possible are 'bigdft_to_ccs' and 'ccs_to_bigdft'")
       end select

       if (iconv==1) then
           call sparse_matrix_and_matrices_init_from_file_bigdft(trim(infile), bigdft_mpi%iproc, bigdft_mpi%nproc, &
                smat, mat, init_matmul=.false.)
           row_ind = f_malloc_ptr(smat%nvctr,id='row_ind')
           col_ptr = f_malloc_ptr(smat%nfvctr,id='col_ptr')
           call ccs_data_from_sparse_matrix(smat, row_ind, col_ptr)
           if (bigdft_mpi%iproc==0) call ccs_matrix_write(trim(outfile), smat, row_ind, col_ptr, mat)
           call f_free_ptr(row_ind)
           call f_free_ptr(col_ptr)
       else if (iconv==2) then
           call sparse_matrix_and_matrices_init_from_file_ccs(trim(infile), bigdft_mpi%iproc, bigdft_mpi%nproc, &
                bigdft_mpi%mpi_comm, smat, mat, init_matmul=.false.)
           call write_sparse_matrix(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
                smat, mat, trim(outfile))
       end if

       call deallocate_sparse_matrix(smat)
       call deallocate_matrices(mat)
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
