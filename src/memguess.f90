!> @file
!!   Program to guess the used memory by BigDFT
!! @author
!!   Copyright (C) 2007-2011 BigDFT group (LG)
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
   integer :: norbe,norbsc,nspin,iorb,norbu,norbd,nspinor,norb,iorbp,iorb_out
   integer :: norbgpu,nspin_ig,ng,ncount0,ncount1,ncount_max,ncount_rate
   integer :: export_wf_iband, export_wf_ispin, export_wf_ikpt, export_wf_ispinor,irad
   real(gp) :: peakmem,hx,hy,hz,tcpu0,tcpu1,tel,energy
   type(input_variables) :: in
   type(atoms_data) :: atoms
   type(orbitals_data) :: orbs,orbstst
   type(communications_arrays) :: comms
   type(local_zone_descriptors) :: Lzd
   type(nonlocal_psp_descriptors) :: nlpspd
   type(gaussian_basis) :: G !basis for davidson IG
   type(denspot_distribution) :: dpbox
   real(gp), dimension(3) :: shift
   logical, dimension(:,:,:), allocatable :: logrid
   integer, dimension(:,:), allocatable :: norbsc_arr
   real(gp), dimension(:,:), pointer :: rxyz, fxyz
   real(wp), dimension(:), allocatable :: rhoexpo,psi
   real(wp), dimension(:,:,:,:), pointer :: rhocoeff
   real(kind=8), dimension(:,:), allocatable :: radii_cf
   logical, dimension(:,:,:), allocatable :: scorb
   real(kind=8), dimension(:), allocatable :: locrad
   real(gp), dimension(:), pointer :: gbd_occ
   type(system_fragment), dimension(:), pointer :: ref_frags
   character(len=3) :: in_name !lr408
   integer :: i
   !! By Ali
   integer :: ierror

   call f_malloc_set_status(memory_limit=0.e0)

   ! Get arguments
   !call getarg(1,tatonam)
   call get_command_argument(1, value = tatonam, status = istat)

   write(radical, "(A)") ""
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
      read(unit=tatonam,fmt=*) nproc
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
            read(tatonam,*,iostat=ierror)ntimes
            if (ierror==0) then
               write(*,'(1x,a,i0,a)')&
                  &   'Repeat each calculation ',ntimes,' times.'
               i_arg = i_arg + 1
               call get_command_argument(i_arg, value = tatonam)
               !call getarg(i_arg,tatonam)
               read(tatonam,*,iostat=ierror)norbgpu
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
            exportwf=.true.
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = filename_wfn)
            !call getarg(i_arg,filename_wfn)
            write(*,'(1x,3a)')&
               &   'export wavefunction file: "', trim(filename_wfn),'" in .cube format'
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
               read(tatonam, *) export_wf_iband
            end if
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            if(trim(tatonam)=='' .or. istat > 0) then
               exit loop_getargs
            else
               read(tatonam, *) export_wf_ispin
            end if
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            if(trim(tatonam)=='' .or. istat > 0) then
               exit loop_getargs
            else
               read(tatonam, *) export_wf_ikpt
            end if
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam, status = istat)
            if(trim(tatonam)=='' .or. istat > 0) then
               exit loop_getargs
            else
               read(tatonam, *) export_wf_ispinor
            end if
            exit loop_getargs
         else if (trim(tatonam)=='atwf') then
            atwf=.true.
            write(*,'(1x,a)')&
               &   'Perform the calculation of atomic wavefunction of the first atom'
            i_arg = i_arg + 1
            call get_command_argument(i_arg, value = tatonam)
            !call getarg(i_arg,tatonam)
            read(tatonam,*,iostat=ierror)ng
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
            if (trim(radical) /= "") then
               write(*,'(1x,a)')&
                  &   'Usage: ./memguess <nproc> [y]'
               write(*,'(1x,a)')&
                  &   'Indicate the number of processes after the executable'
               write(*,'(1x,a)')&
                  &   'ERROR: The only second argument which is accepted is "y", "o","convert", "GPUtest" or "atwf" ' 
               write(*,'(1x,a)')&
                  &   '       (type "memguess" without arguments to have an help)'
               stop
            end if
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
      atoms%astruct%geocode = "P"
      write(*,*) "Read density file..."
      call read_density(trim(fileFrom), atoms%astruct%geocode, Lzd%Glr%d%n1i, Lzd%Glr%d%n2i, Lzd%Glr%d%n3i, &
         &   nspin, hx, hy, hz, rhocoeff, atoms%astruct%nat, rxyz, atoms%astruct%iatype, atoms%nzatom)
      atoms%astruct%ntypes = size(atoms%nzatom) - ndebug
      write(*,*) "Write new density file..."
      dpbox%ndims(1)=Lzd%Glr%d%n1i
      dpbox%ndims(2)=Lzd%Glr%d%n2i
      dpbox%ndims(3)=Lzd%Glr%d%n3i
      dpbox%hgrids(1)=hx
      dpbox%hgrids(2)=hy
      dpbox%hgrids(3)=hz
      allocate(dpbox%ngatherarr(0:0,2+ndebug),stat=i_stat)
      call memocc(i_stat,dpbox%ngatherarr,'ngatherarr',subname)

      call plot_density(0,1,trim(fileTo),atoms,rxyz,dpbox,nspin,rhocoeff)
      write(*,*) "Done"
      stop
   end if
   if (convertpos) then
      call read_atomic_file(trim(fileFrom),0,atoms%astruct,i_stat,fcomment,energy,fxyz)
      call allocate_atoms_nat(atoms, subname)
      call allocate_atoms_ntypes(atoms, subname)
      rxyz=>atoms%astruct%rxyz
      if (i_stat /=0) stop 'error on input file parsing' 
      !find the format of the output file
      if (index(fileTo,'.xyz') > 0) then
         irad=index(fileTo,'.xyz')
         atoms%astruct%inputfile_format='xyz  '
      else if (index(fileTo,'.ascii') > 0) then
         irad=index(fileTo,'.ascii')
         atoms%astruct%inputfile_format='ascii'
      else if (index(fileTo,'.yaml') > 0) then
         irad=index(fileTo,'.yaml')
         atoms%astruct%inputfile_format='yaml '
      else
         irad = len(trim(fileTo)) + 1
      end if
      
      if (associated(fxyz)) then
         call write_atomic_file(fileTo(1:irad-1),energy,rxyz,atoms,&
              trim(fcomment) // ' (converted from '//trim(fileFrom)//")", fxyz)

         i_all=-product(shape(fxyz))*kind(fxyz)
         deallocate(fxyz,stat=i_stat)
         call memocc(i_stat,i_all,'fxyz',subname)
      else
         call write_atomic_file(fileTo(1:irad-1),energy,rxyz,atoms,&
              trim(fcomment) // ' (converted from '//trim(fileFrom)//")")
      end if
      stop
   end if

   !standard names
   call standard_inputfile_names(in, radical, 1)


   if (trim(radical) == "") then
      posinp='posinp'
   else
      posinp=trim(radical)
   end if
   call bigdft_set_input(radical,posinp,rxyz,in,atoms)
   !initialize memory counting
   !call memocc(0,0,'count','start')

   if (in%ixc < 0) then
      call xc_init(in%ixc, XC_MIXED, nspin)
   else
      call xc_init(in%ixc, XC_ABINIT, nspin)
   end if

   call print_general_parameters(nproc,in,atoms)
   call print_dft_parameters(in,atoms)
   call xc_dump()

   !Time initialization
   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max)

   ! store PSP parameters
   allocate(radii_cf(atoms%astruct%ntypes,3+ndebug),stat=i_stat)
   call memocc(i_stat,radii_cf,'radii_cf',subname)

   call system_properties(0,nproc,in,atoms,orbs,radii_cf)

   if (optimise) then
      if (atoms%astruct%geocode =='F') then
         call optimise_volume(atoms,in%crmult,in%frmult,in%hx,in%hy,in%hz,rxyz,radii_cf)
      else
         call shift_periodic_directions(atoms,rxyz,radii_cf)
      end if
      write(*,'(1x,a)')'Writing optimised positions in file posopt.[xyz,ascii]...'
      write(comment,'(a)')'POSITIONS IN OPTIMIZED CELL '
      call write_atomic_file('posopt',0.d0,rxyz,atoms,trim(comment))
      !call wtxyz('posopt',0.d0,rxyz,atoms,trim(comment))
   end if

   !In the case in which the number of orbitals is not "trivial" check whether they are too many
   !Always True! (TD)
   !if ( max(orbs%norbu,orbs%norbd) /= ceiling(real(nelec,kind=4)/2.0) .or. .true.) then
      ! Allocations for readAtomicOrbitals (check inguess.dat and psppar files + give norbe)
      allocate(scorb(4,2,atoms%natsc+ndebug),stat=i_stat)
      call memocc(i_stat,scorb,'scorb',subname)
      allocate(norbsc_arr(atoms%natsc+1,in%nspin+ndebug),stat=i_stat)
      call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
      allocate(locrad(atoms%astruct%nat+ndebug),stat=i_stat)
      call memocc(i_stat,locrad,'locrad',subname)

      !calculate the inputguess orbitals
      !spin for inputguess orbitals
      if (in%nspin==4) then
         nspin_ig=1
      else
         nspin_ig=in%nspin
      end if

      ! Read the inguess.dat file or generate the input guess via the inguess_generator
      call readAtomicOrbitals(atoms,norbe,norbsc,nspin_ig,orbs%nspinor,&
         &   scorb,norbsc_arr,locrad)

      if (in%nspin==4) then
         !in that case the number of orbitals doubles
         norbe=2*norbe
      end if

      ! De-allocations
      i_all=-product(shape(locrad))*kind(locrad)
      deallocate(locrad,stat=i_stat)
      call memocc(i_stat,i_all,'locrad',subname)
      i_all=-product(shape(scorb))*kind(scorb)
      deallocate(scorb,stat=i_stat)
      call memocc(i_stat,i_all,'scorb',subname)
      i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
      deallocate(norbsc_arr,stat=i_stat)
      call memocc(i_stat,i_all,'norbsc_arr',subname)

      if (in%inputpsiId /= INPUT_PSI_RANDOM) then
         ! Check the maximum number of orbitals
         if (in%nspin==1 .or. in%nspin==4) then
            if (orbs%norb>norbe) then
               write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals (',orbs%norb,&
                  &   ') must not be greater than the number of orbitals (',norbe,&
                  &   ') generated from the input guess.'
               stop
            end if
         else if (in%nspin == 2) then
            if (orbs%norbu > norbe) then
               write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals up (',orbs%norbu,&
                  &   ') must not be greater than the number of orbitals (',norbe,&
                  &   ') generated from the input guess.'
               stop
            end if
            if (orbs%norbd > norbe) then
               write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals down (',orbs%norbd,&
                  &   ') must not be greater than the number of orbitals (',norbe,&
                  &   ') generated from the input guess.'
               stop
            end if
         end if
      end if

   !end if

   ! Determine size alat of overall simulation cell and shift atom positions
   ! then calculate the size in units of the grid space
   hx=in%hx
   hy=in%hy
   hz=in%hz

   call system_size(0,atoms,rxyz,radii_cf,in%crmult,in%frmult,hx,hy,hz,Lzd%Glr,shift)

   ! Build and print the communicator scheme.
   call createWavefunctionsDescriptors(0,hx,hy,hz,&
      &   atoms,rxyz,radii_cf,in%crmult,in%frmult,Lzd%Glr, output_denspot = (output_grid > 0))
   call orbitals_communicators(0,nproc,Lzd%Glr,orbs,comms)  

   if (exportwf) then

      allocate(psi((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor+ndebug),stat=i_stat)
      call memocc(i_stat,psi,'psi',subname)

      ! Optionally compute iorbp from arguments in case of ETSF.
      if (export_wf_ikpt < 1 .or. export_wf_ikpt > orbs%nkpts) stop "Wrong k-point"
      if (export_wf_ispin < 1 .or. export_wf_ispin > orbs%nspin) stop "Wrong spin"
      if ((export_wf_ispin == 1 .and. &
           & (export_wf_iband < 1 .or. export_wf_iband > orbs%norbu)) .or. &
           & (export_wf_ispin == 0 .and. &
           & (export_wf_iband < 1 .or. export_wf_iband > orbs%norbd))) stop "Wrong orbital"
      iorbp = (export_wf_ikpt - 1) * orbs%norb + (export_wf_ispin - 1) * orbs%norbu + export_wf_iband

      ! ref_frags to be allocated here
      i = index(filename_wfn, "/",back=.true.)+1
      read(filename_wfn(i:i+3),*) in_name ! lr408
      if (in_name == 'min') then
         stop 'Ref fragment not initialized, linear reading currently nonfunctional, to be fixed'
      end if

      call take_psi_from_file(filename_wfn,in%frag,hx,hy,hz,Lzd%Glr, &
           & atoms,rxyz,orbs,psi,iorbp,export_wf_ispinor,ref_frags)
      call filename_of_iorb(.false.,"wavefunction",orbs,iorbp, &
           & export_wf_ispinor,filename_wfn,iorb_out)
      call plot_wf(filename_wfn,1,atoms,1.0_wp,Lzd%Glr,hx,hy,hz,rxyz, &
           & psi((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f) * (export_wf_ispinor - 1) + 1))

      i_all=-product(shape(psi))*kind(psi)
      deallocate(psi,stat=i_stat)
      call memocc(i_stat,i_all,'psi',subname)

   end if

   if (GPUtest) then
      !test the hamiltonian in CPU or GPU
      !create the orbitals data structure for one orbital
      !test orbitals
      nspin=1
      if (norbgpu == 0) then
         norb=orbs%norb
      else
         norb=norbgpu
      end if
      norbu=norb
      norbd=0
      nspinor=1

      call orbitals_descriptors(0,nproc,norb,norbu,norbd,in%nspin,nspinor, &
           in%nkpt,in%kpt,in%wkpt,orbstst,.false.)
      allocate(orbstst%eval(orbstst%norbp+ndebug),stat=i_stat)
      call memocc(i_stat,orbstst%eval,'orbstst%eval',subname)
      do iorb=1,orbstst%norbp
         orbstst%eval(iorb)=-0.5_gp
      end do

      do iorb=1,orbstst%norb
         orbstst%occup(iorb)=1.0_gp
         orbstst%spinsgn(iorb)=1.0_gp
      end do

      call check_linear_and_create_Lzd(0,1,in%linear,Lzd,atoms,orbstst,in%nspin,rxyz)

      !for the given processor (this is only the cubic strategy)
      orbstst%npsidim_orbs=(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbstst%norbp*orbstst%nspinor
      orbstst%npsidim_comp=1


      call compare_cpu_gpu_hamiltonian(0,1,in%matacc,atoms,&
           orbstst,nspin,in%ncong,in%ixc,&
           Lzd,hx,hy,hz,rxyz,ntimes)

      call deallocate_orbs(orbstst,subname)


      i_all=-product(shape(orbstst%eval))*kind(orbstst%eval)
      deallocate(orbstst%eval,stat=i_stat)
      call memocc(i_stat,i_all,'orbstst%eval',subname)

   end if

   call deallocate_comms(comms,subname)

   ! determine localization region for all projectors, but do not yet fill the descriptor arrays
   allocate(nlpspd%plr(atoms%astruct%nat))
!!$   allocate(nlpspd%nseg_p(0:2*atoms%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,nlpspd%nseg_p,'nseg_p',subname)
!!$   allocate(nlpspd%nvctr_p(0:2*atoms%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,nlpspd%nvctr_p,'nvctr_p',subname)
!!$   allocate(nlpspd%nboxp_c(2,3,atoms%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,nlpspd%nboxp_c,'nboxp_c',subname)
!!$   allocate(nlpspd%nboxp_f(2,3,atoms%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,nlpspd%nboxp_f,'nboxp_f',subname)

   allocate(logrid(0:Lzd%Glr%d%n1,0:Lzd%Glr%d%n2,0:Lzd%Glr%d%n3+ndebug),stat=i_stat)
   call memocc(i_stat,logrid,'logrid',subname)

   call localize_projectors(0,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,hx,hy,hz,&
      &   in%frmult,in%frmult,rxyz,radii_cf,logrid,atoms,orbs,nlpspd)
   deallocate(nlpspd%plr)
   !allocations for arrays holding the data descriptors
!!$   !just for modularity
!!$   allocate(nlpspd%keyg_p(2,nlpspd%nseg_p(2*atoms%astruct%nat)+ndebug),stat=i_stat)
!!$   call memocc(i_stat,nlpspd%keyg_p,'nlpspd%keyg_p',subname)
!!$   allocate(nlpspd%keyv_p(nlpspd%nseg_p(2*atoms%astruct%nat)+ndebug),stat=i_stat)
!!$   call memocc(i_stat,nlpspd%keyv_p,'nlpspd%keyv_p',subname)

   if (atwf) then
      !here the treatment of the AE Core charge density
      !number of gaussians defined in the input of memguess
      !ng=31
      !plot the wavefunctions for the pseudo atom
      nullify(G%rxyz)
      call gaussian_pswf_basis(ng,.false.,0,in%nspin,atoms,rxyz,G,gbd_occ)
      !for the moment multiply the number of coefficients for each channel
      allocate(rhocoeff((ng*(ng+1))/2,4,1,1+ndebug),stat=i_stat)
      call memocc(i_stat,rhocoeff,'rhocoeff',subname)
      allocate(rhoexpo((ng*(ng+1))/2+ndebug),stat=i_stat)
      call memocc(i_stat,rhoexpo,'rhoexpo',subname)

      call plot_gatom_basis('gatom',1,ng,G,gbd_occ,rhocoeff,rhoexpo)

      if (associated(gbd_occ)) then
         i_all=-product(shape(gbd_occ))*kind(gbd_occ)
         deallocate(gbd_occ,stat=i_stat)
         call memocc(i_stat,i_all,'gbd_occ',subname)
         nullify(gbd_occ)
      end if
      !deallocate the gaussian basis descriptors
      call deallocate_gwf(G,subname)

      !!$  !plot the wavefunctions for the AE atom
      !!$  !not possible, the code should recognize the AE eleconf
      !!$  call razero(35,atoms%psppar(0,0,atoms%astruct%iatype(1)))
      !!$  atoms%psppar(0,0,atoms%astruct%iatype(1))=0.01_gp
      !!$  nullify(G%rxyz)
      !!$  call gaussian_pswf_basis(ng,.false.,0,in%nspin,atoms,rxyz,G,gbd_occ)
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

      i_all=-product(shape(rhoexpo))*kind(rhoexpo)
      deallocate(rhoexpo,stat=i_stat)
      call memocc(i_stat,i_all,'rhoexpo',subname)
      i_all=-product(shape(rhocoeff))*kind(rhocoeff)
      deallocate(rhocoeff,stat=i_stat)
      call memocc(i_stat,i_all,'rhocoeff',subname)

   end if

   i_all=-product(shape(logrid))*kind(logrid)
   deallocate(logrid,stat=i_stat)
   call memocc(i_stat,i_all,'logrid',subname)

   !call deallocate_proj_descr(nlpspd,subname)

   call MemoryEstimator(nproc,in%idsx,Lzd%Glr,&
        atoms%astruct%nat,orbs%norb,orbs%nspinor,orbs%nkpts,nlpspd%nprojel,&
        in%nspin,in%itrpmax,in%iscf,peakmem)

   !add the comparison between cuda hamiltonian and normal one if it is the case

   call deallocate_atoms(atoms,subname)

   call deallocate_lr(Lzd%Glr,subname)

   call xc_end()

   i_all=-product(shape(radii_cf))*kind(radii_cf)
   deallocate(radii_cf,stat=i_stat)
   call memocc(i_stat,i_all,'radii_cf',subname)
   i_all=-product(shape(rxyz))*kind(rxyz)
   deallocate(rxyz,stat=i_stat)
   call memocc(i_stat,i_all,'rxyz',subname)

   ! De-allocations
   call deallocate_orbs(orbs,subname)

   !remove the directory which has been created if it is possible
   call deldir(in%dir_output,len(trim(in%dir_output)),ierror)
   call free_input_variables(in)  

   !finalize memory counting
   call memocc(0,0,'count','stop')

   !Elapsed time
   call cpu_time(tcpu1)
   call system_clock(ncount1,ncount_rate,ncount_max)
   tel=dble(ncount1-ncount0)/dble(ncount_rate)
   write( *,'(1x,a,2(1x,f12.2))') 'CPU time/ELAPSED time ', tel, tcpu1-tcpu0

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

   allocate(txyz(3,atoms%astruct%nat+ndebug),stat=i_stat)
   call memocc(i_stat,txyz,'txyz',subname)
   call system_size(1,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
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

      call system_size(1,atoms,txyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
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

   i_all=-product(shape(txyz))*kind(txyz)
   deallocate(txyz,stat=i_stat)
   call memocc(i_stat,i_all,'txyz',subname)

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

   allocate(txyz(3,at%astruct%nat+ndebug),stat=i_stat)
   call memocc(i_stat,txyz,'txyz',subname)

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

   i_all=-product(shape(txyz))*kind(txyz)
   deallocate(txyz,stat=i_stat)
   call memocc(i_stat,i_all,'txyz',subname)

END SUBROUTINE shift_periodic_directions


subroutine calc_vol(geocode,nat,rxyz,vol)
   use module_base
   implicit none
   character(len=1), intent(in) :: geocode
   integer, intent(in) :: nat
   real(gp), dimension(3,nat), intent(in) :: rxyz
   real(gp), intent(out) :: vol
   !local variables
   integer :: iat
   real(gp) :: cxmin,cxmax,cymin,cymax,czmin,czmax

   !calculate the extremes of the boxes taking into account the spheres around the atoms
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

   allocate(gaucoeffs(G%ncoeff,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
   call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)

   !fill randomly the gaussian coefficients for the orbitals considered
   do iorb=1,orbs%norbp*orbs%nspinor
      do icoeff=1,G%ncoeff
         call random_number(tt)
         gaucoeffs(icoeff,iorb)=real(tt,wp)
      end do
   end do

   !allocate the wavefunctions
   allocate(psi(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
   call memocc(i_stat,psi,'psi',subname)
   allocate(hpsi(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
   call memocc(i_stat,hpsi,'hpsi',subname)

   call razero(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,psi)
   call razero(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,hpsi)

   !convert the gaussians in wavelets
   call gaussians_to_wavelets(iproc,nproc,at%astruct%geocode,orbs,Lzd%Glr%d,&
           hx,hy,hz,Lzd%Glr%wfd,G,gaucoeffs,psi)

   i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
   deallocate(gaucoeffs,stat=i_stat)
   call memocc(i_stat,i_all,'gaucoeffs',subname)

   i_all=-product(shape(gbd_occ))*kind(gbd_occ)
   deallocate(gbd_occ,stat=i_stat)
   call memocc(i_stat,i_all,'gbd_occ',subname)

   !deallocate the gaussian basis descriptors
   call deallocate_gwf(G,subname)

   !allocate and initialise the potential and the density
   allocate(pot(Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,nspin+ndebug),stat=i_stat)
   call memocc(i_stat,pot,'pot',subname)
   allocate(rho(Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,nspin+ndebug),stat=i_stat)
   call memocc(i_stat,rho,'rho',subname)

   !here the potential can be used for building the density
   allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
   call memocc(i_stat,nscatterarr,'nscatterarr',subname)
   allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
   call memocc(i_stat,nscatterarr,'nscatterarr',subname)

   !normally nproc=1
   do jproc=0,nproc-1
      call PS_dim4allocation(at%astruct%geocode,'D',jproc,nproc,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,xc_isgga(),(ixc/=13),&
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
   rsflag = .not.xc_isgga()

   !calculate dimensions of the complete array to be allocated before the reduction procedure
   if (rsflag) then
      nrhotot=0
      do jproc=0,nproc-1
         nrhotot=nrhotot+nscatterarr(jproc,1)
      end do
   else
      nrhotot=Lzd%Glr%d%n3i
   end if

   call local_potential_dimensions(Lzd,orbs,ngatherarr(0,1))

   !allocate the necessary objects on the GPU
   !set initialisation of GPU part 
   !initialise the acceleration strategy if required
   call init_material_acceleration(iproc,matacc,GPU)

   if (GPUconv .eqv. OCLconv) stop 'ERROR: One (and only one) acceleration should be present with GPUtest'

   !allocate arrays for the GPU if a card is present
   if (GPUconv) then
      call prepare_gpu_for_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,nspin,&
           hx,hy,hz,Lzd%Glr%wfd,orbs,GPU)
   else if (OCLconv) then
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
      else if (OCLconv) then
         call local_partial_density_OCL(orbs,nrhotot,Lzd%Glr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,psi,rho,GPU)
      end if
   end do
   call nanosec(itsc1)
   !call cpu_time(t1)
   !GPUtime=real(t1-t0,kind=8)
   GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   i_all=-product(shape(nscatterarr))*kind(nscatterarr)
   deallocate(nscatterarr,stat=i_stat)
   call memocc(i_stat,i_all,'nscatterarr',subname)
   i_all=-product(shape(ngatherarr))*kind(ngatherarr)
   deallocate(ngatherarr,stat=i_stat)
   call memocc(i_stat,i_all,'ngatherarr',subname)



   !compare the results between the different actions of the hamiltonian
   !check the differences between the results
   call compare_data_and_gflops(CPUtime,GPUtime,&
        & real(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,kind=8)*192.d0,pot,rho,&
        Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,ntimes*orbs%norbp,.false.,Rden)

   i_all=-product(shape(rho))*kind(rho)
   deallocate(rho,stat=i_stat)
   call memocc(i_stat,i_all,'rho',subname)


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
   do j=1,ntimes
      allocate(pottmp(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*(nspin+ndebug)),stat=i_stat)
      call memocc(i_stat,pottmp,'pottmp',subname)
      call dcopy(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*(nspin+ndebug),pot(1,1,1,1),1,pottmp(1),1)
      call local_hamiltonian(iproc,nproc,orbs%npsidim_orbs,orbs,Lzd,hx,hy,hz,0,confdatarr,pottmp,psi,hpsi, &
           fake_pkernelSIC,0,0.0_gp,ekin_sum,epot_sum,eSIC_DC)
      i_all=-product(shape(pottmp))*kind(pottmp)
      deallocate(pottmp,stat=i_stat)
      call memocc(i_stat,i_all,'pottmp',subname)
   end do
   call nanosec(itsc1)
   CPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   print *,'ekin,epot=',ekin_sum,epot_sum

   !WARNING: local hamiltonian overwrites the psis
   !warm-up
   !call gpu_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,hx,hy,hz,orbs,GPU,ekinGPU,epotGPU)

   !apply the GPU hamiltonian and put the results in the hpsi_GPU array
   allocate(GPU%hpsi_ASYNC((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),stat=i_stat)
   call memocc(i_stat,GPU%hpsi_ASYNC,'GPU%hpsi_ASYNC',subname)

   !take timings
   call nanosec(itsc0)
   do j=1,ntimes
      if (GPUconv) then
         call local_hamiltonian_GPU(orbs,Lzd%Glr,hx,hy,hz,orbs%nspin,pot,psi,GPU%hpsi_ASYNC,ekinGPU,epotGPU,GPU)
      else if (OCLconv) then
         call local_hamiltonian_OCL(orbs,Lzd%Glr,hx,hy,hz,orbs%nspin,pot,psi,GPU%hpsi_ASYNC,ekinGPU,epotGPU,GPU)
      end if
   end do
   if(ASYNCconv .and. OCLconv) call finish_hamiltonian_OCL(orbs,ekinGPU,epotGPU,GPU)
   call nanosec(itsc1)
   GPUtime=real(itsc1-itsc0,kind=8)*1.d-9

   print *,'ekinGPU,epotGPU',ekinGPU,epotGPU

   !compare the results between the different actions of the hamiltonian
   !check the differences between the results
   call compare_data_and_gflops(CPUtime,GPUtime,&
      &   real(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,kind=8)*real(192+46*3+192+2,kind=8),hpsi,GPU%hpsi_ASYNC,&
   orbs%norbp*orbs%nspinor*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f),ntimes*orbs%norbp,.false.,Rham)

   i_all=-product(shape(pot))*kind(pot)
   deallocate(pot,stat=i_stat)
   call memocc(i_stat,i_all,'pot',subname)

   write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Linear Algebra (Blas)'

   !perform the scalar product between the hpsi wavefunctions
   !actually this is <hpsi|hpsi> it has no meaning.
   !this works only if nspinor==1
   allocate(overlap(orbs%norbp,orbs%norbp,2+ndebug),stat=i_stat)
   call memocc(i_stat,overlap,'overlap',subname)

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

   i_all=-product(shape(overlap))*kind(overlap)
   deallocate(overlap,stat=i_stat)
   call memocc(i_stat,i_all,'overlap',subname)


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
      else if (OCLconv) then
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

   i_all=-product(shape(GPU%hpsi_ASYNC))*kind(GPU%hpsi_ASYNC)
   deallocate(GPU%hpsi_ASYNC,stat=i_stat)
   call memocc(i_stat,i_all,'GPU%hpsi_ASYNC',subname)
   i_all=-product(shape(psi))*kind(psi)
   deallocate(psi,stat=i_stat)
   call memocc(i_stat,i_all,'psi',subname)
   i_all=-product(shape(hpsi))*kind(hpsi)
   deallocate(hpsi,stat=i_stat)
   call memocc(i_stat,i_all,'hpsi',subname)


   !free the card at the end
   if (GPUconv) then
      call free_gpu(GPU,orbs%norbp)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,orbs,nspin)
   end if

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

   allocate(rxyz_file(at%astruct%nat,3+ndebug),stat=i_stat)
   call memocc(i_stat,rxyz_file,'rxyz_file',subname)

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

      allocate(psifscf(-nb1:2*lr%d%n1+1+nb1,-nb2:2*lr%d%n2+1+nb2, &
           & -nb3:2*lr%d%n3+1+nb3+ndebug),stat=i_stat)
      call memocc(i_stat,psifscf,'psifscf',subname)

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

         allocate(lpsi(1:Lzd%llr(1)%wfd%nvctr_c+7*Lzd%llr(1)%wfd%nvctr_f))
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
              lr%d%n1,lr%d%n2,lr%d%n3,hx,hy,hz,at,Lzd%llr(1)%wfd,rxyz_file,rxyz,&
              locrad,locregCenter,confPotOrder,confPotPrefac,&
              lpsi(1),eval_fake,psifscf)

         call to_zero(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,psi(1,1))

         call Lpsi_to_global2(0,Lzd%llr(1)%wfd%nvctr_c+7*Lzd%llr(1)%wfd%nvctr_f, &
              lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,1,1,1,lr,Lzd%Llr(1),lpsi,psi)

         deallocate(lpsi)
      end if

      ! Update iorbp
      iorbp = (ikpt - 1) * orbs%norb + (ispin - 1) * orbs%norbu + iorbp

      close(99)

      i_all=-product(shape(psifscf))*kind(psifscf)
      deallocate(psifscf,stat=i_stat)
      call memocc(i_stat,i_all,'psifscf',subname)

   else if (iformat == WF_FORMAT_ETSF) then
      call read_one_wave_etsf(0,filename,iorbp,0,orbs%nspinor,lr%d%n1,lr%d%n2,lr%d%n3,&
           & hx,hy,hz,at,rxyz_file,rxyz,lr%wfd,psi,eval_fake)
   end if
   i_all=-product(shape(rxyz_file))*kind(rxyz_file)
   deallocate(rxyz_file,stat=i_stat)
   call memocc(i_stat,i_all,'rxyz_file',subname)
END SUBROUTINE take_psi_from_file

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
