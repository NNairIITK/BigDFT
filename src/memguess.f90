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
   use module_atoms, only: set_astruct_from_file
   implicit none
   character(len=*), parameter :: subname='memguess'
   character(len=20) :: tatonam, radical
   character(len=40) :: comment
   character(len=1024) :: fcomment
   character(len=128) :: fileFrom, fileTo,filename_wfn
   character(len=50) :: posinp
   logical :: optimise,GPUtest,atwf,convert=.false.,exportwf=.false.
   logical :: disable_deprecation = .false.,convertpos=.false.
   integer :: ntimes,nproc,i_stat,i_all,output_grid, i_arg,istat
   integer :: nspin,iorb,norbu,norbd,nspinor,norb,iorbp,iorb_out
   integer :: norbgpu,ng
   integer :: export_wf_iband, export_wf_ispin, export_wf_ikpt, export_wf_ispinor,irad
   real(gp) :: hx,hy,hz,energy
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
   integer :: i, inputpsi, input_wf_format
   !real(gp) :: tcpu0,tcpu1,tel
   !integer :: ncount0,ncount1,ncount_max,ncount_rate
   !! By Ali
   integer :: ierror

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
         else if (trim(tatonam) == 'dd') then
            ! dd: disable deprecation message
            disable_deprecation = .true.
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
      at%astruct%ntypes = size(at%nzatom) - ndebug
      write(*,*) "Write new density file..."
      dpbox%ngatherarr = f_malloc_ptr((/ 0.to.0, 1.to.2 /),id='dpbox%ngatherarr')

      call plot_density(0,1,trim(fileTo),at,at%astruct%rxyz,dpbox,nspin,rhocoeff)
      write(*,*) "Done"
      stop
   end if

   if (convertpos) then
      call set_astruct_from_file(trim(fileFrom),0,at%astruct,i_stat,fcomment,energy,fxyz)
      if (i_stat /=0) stop 'error on input file parsing' 
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
         call write_atomic_file(fileTo(1:irad-1),energy,at%astruct%rxyz,at,&
              trim(fcomment) // ' (converted from '//trim(fileFrom)//")", fxyz)

         call f_free_ptr(fxyz)
      else
         call write_atomic_file(fileTo(1:irad-1),energy,at%astruct%rxyz,at,&
              trim(fcomment) // ' (converted from '//trim(fileFrom)//")")
      end if
      stop
   end if

   if (trim(radical) == "input") then
      posinp='posinp'
   else
      posinp=trim(radical)
   end if

   call run_objects_init_from_files(runObj, radical, posinp)

   if (optimise) then
      if (runObj%atoms%astruct%geocode =='F') then
         call optimise_volume(runObj%atoms,&
              & runObj%inputs%crmult,runObj%inputs%frmult,&
              & runObj%inputs%hx,runObj%inputs%hy,runObj%inputs%hz,&
              & runObj%atoms%astruct%rxyz,runObj%radii_cf)
      else
         call shift_periodic_directions(runObj%atoms,runObj%atoms%astruct%rxyz,runObj%radii_cf)
      end if
      write(*,'(1x,a)')'Writing optimised positions in file posopt.[xyz,ascii]...'
      write(comment,'(a)')'POSITIONS IN OPTIMIZED CELL '
      call write_atomic_file('posopt',0.d0,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment))
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
        & shift,runObj%radii_cf, ref_frags, output_grid = (output_grid > 0))
   call MemoryEstimator(nproc,runObj%inputs%idsx,runObj%rst%KSwfn%Lzd%Glr,&
        & runObj%rst%KSwfn%orbs%norb,runObj%rst%KSwfn%orbs%nspinor,&
        & runObj%rst%KSwfn%orbs%nkpts,nlpsp%nprojel,&
        runObj%inputs%nspin,runObj%inputs%itrpmax,runObj%inputs%iscf,mem)
   
   if (.not. exportwf) then
      call print_memory_estimation(mem)
   else
      runObj%rst%KSwfn%psi = f_malloc_ptr((runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_c+&
           & 7*runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_f)*runObj%rst%KSwfn%orbs%nspinor+ndebug,&
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

      call plot_wf(filename_wfn,1,runObj%atoms,1.0_wp,runObj%rst%KSwfn%Lzd%Glr, &
           & runObj%inputs%hx,runObj%inputs%hy,runObj%inputs%hz,runObj%atoms%astruct%rxyz, &
           & runObj%rst%KSwfn%psi((runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_c+&
           & 7*runObj%rst%KSwfn%Lzd%Glr%wfd%nvctr_f) * (export_wf_ispinor - 1) + 1))
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
           runObj%inputs%gen_nkpt,runObj%inputs%gen_kpt,runObj%inputs%gen_wkpt,orbstst,.false.)
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

      call deallocate_orbs(orbstst,subname)


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
      call deallocate_gwf(G,subname)

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
      !!$  call deallocate_gwf(G,subname)

      call f_free(rhoexpo)
      call f_free_ptr(rhocoeff)

   end if

   ! Add the comparison between cuda hamiltonian and normal one if it is the case

   ! De-allocations
   call deallocate_Lzd_except_Glr(runObj%rst%KSwfn%Lzd, subname)
   call deallocate_comms(runObj%rst%KSwfn%comms,subname)
   call deallocate_orbs(runObj%rst%KSwfn%orbs,subname)
   call free_DFT_PSP_projectors(nlpsp)

   !remove the directory which has been created if it is possible
   call deldir(runObj%inputs%dir_output,len(trim(runObj%inputs%dir_output)),ierror)

   call run_objects_free(runObj, subname)
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


!>  Rotate the molecule via an orthogonal matrix in order to minimise the
!!  volume of the cubic cell
subroutine optimise_volume(atoms,crmult,frmult,hx,hy,hz,rxyz,radii_cf)
   use module_base
   use module_types
   implicit none
   type(atoms_data), intent(inout) :: atoms
   real(gp), intent(in) :: crmult,frmult,hx,hy,hz
   real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
   real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
   !local variables
   character(len=*), parameter :: subname='optimise_volume'
   integer :: iat,i_all,i_stat,it,i
   real(gp) :: x,y,z,vol,tx,ty,tz,tvol,s,diag,dmax
   type(locreg_descriptors) :: Glr
   real(gp), dimension(3) :: shift
   real(gp), dimension(3,3) :: urot
   real(gp), dimension(:,:), allocatable :: txyz

   txyz = f_malloc((/ 3, atoms%astruct%nat /),id='txyz')
   call system_size(atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,.false.,Glr,shift)
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

      call system_size(atoms,txyz,radii_cf,crmult,frmult,hx,hy,hz,.false.,Glr,shift)
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


!>  Add a shift in the periodic directions such that the system
!!  uses as less as possible the modulo operation
subroutine shift_periodic_directions(at,rxyz,radii_cf)
   use module_base
   use module_types
   implicit none
   type(atoms_data), intent(inout) :: at
   real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
   real(gp), dimension(3,at%astruct%nat), intent(inout) :: rxyz
   !local variables
   character(len=*), parameter :: subname='shift_periodic_directions'
   integer :: iat,i_all,i_stat,i,ityp
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
   integer :: icoeff,i_stat,i_all,i1,i2,i3,ispin,j
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
   psi = f_malloc((/ Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f, orbs%nspinor*orbs%norbp /),id='psi')
   hpsi = f_malloc((/ Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f , orbs%nspinor*orbs%norbp+ndebug /),id='hpsi')

   call to_zero(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,psi)
   call to_zero(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,hpsi)

   !convert the gaussians in wavelets
   call gaussians_to_wavelets(iproc,nproc,at%astruct%geocode,orbs,Lzd%Glr%d,&
           hx,hy,hz,Lzd%Glr%wfd,G,gaucoeffs,psi)

   call f_free(gaucoeffs)
   call f_free_ptr(gbd_occ)

   !deallocate the gaussian basis descriptors
   call deallocate_gwf(G,subname)

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

   if (GPUconv .eqv. GPU%OCLconv) stop 'ERROR: One (and only one) acceleration should be present with GPUtest'

   !allocate arrays for the GPU if a card is present
   if (GPUconv) then
      call prepare_gpu_for_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,nspin,&
           hx,hy,hz,Lzd%Glr%wfd,orbs,GPU)
   else if (GPU%OCLconv) then
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
      if (GPUconv) then
         call local_partial_density_GPU(orbs,nrhotot,Lzd%Glr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,psi,rho,GPU)
      else if (GPU%OCLconv) then
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
      pottmp = f_malloc_ptr(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*(nspin+ndebug),id='pottmp')
      call vcopy(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*(nspin+ndebug),pot(1,1,1,1),1,pottmp(1),1)
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
      if (GPUconv) then
         call local_hamiltonian_GPU(orbs,Lzd%Glr,hx,hy,hz,orbs%nspin,pot,psi,GPU%hpsi_ASYNC,ekinGPU,epotGPU,GPU)
      else if (GPU%OCLconv) then
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
      if (GPUconv) then
         call preconditionall_GPU(orbs,Lzd%Glr,hx,hy,hz,ncong,&
            &   GPU%hpsi_ASYNC,gnrmGPU,gnrm_zero,GPU)
      else if (GPU%OCLconv) then
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
   if (GPUconv) then
      call free_gpu(GPU,orbs%norbp)
   else if (GPU%OCLconv) then
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
   integer :: nb1,nb2,nb3,i_stat,i_all, ikpt, ispin, i
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

         call to_zero(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,psi(1,1))

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
