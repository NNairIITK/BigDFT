!!****p* BigDFT/BigDFT
!! FUNCTION
!!  Main program to calculate electronic structures
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA, UNIBAS
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
program abscalc_main

  use module_base
  use module_types
  use module_interfaces
  use ab6_symmetry
!  use minimization, only: parameterminimization 

  !implicit real(kind=8) (a-h,o-z)
  !as a general policy, I will put "implicit none" by assuming the same
  !name convention as "implicit real(kind=8) (a-h,o-z)"
  !such that the implicit statement can be commented at will

  implicit none
  character(len=*), parameter :: subname='BigDFT'
  character(len=20) :: units
  integer :: iproc,nproc,iat,ityp,j,i_stat,i_all,ierr,infocode
  integer ::  ncount_bigdft
  real(gp) :: etot,sumx,sumy,sumz
  logical :: exist_list
  !input variables
  type(atoms_data) :: atoms
  type(input_variables) :: inputs
  type(restart_objects) :: rst
  character(len=20), dimension(:), allocatable :: atomnames
  character(len=50), dimension(:), allocatable :: arr_posinp
  character(len=60)  :: filename
  ! atomic coordinates, forces
  real(gp), dimension(:,:), allocatable :: fxyz
  real(gp), dimension(:,:), pointer :: rxyz
  integer :: npr,iam,iconfig,nconfig
  integer  :: nfluct
  real(gp) :: fluctsum
  logical :: exists

 
  ! Start MPI in parallel version
  !in the case of MPIfake libraries the number of processors is automatically adjusted
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Find out which input files will be used
  inquire(file="list_posinp",exist=exist_list)
  if (exist_list) then
     open(unit=54,file="list_posinp")
     read(54,*) nconfig
     if (nconfig > 0) then 
        !allocation not referenced since memocc count not initialised
        allocate(arr_posinp(1:nconfig))

        do iconfig=1,nconfig
           read(54,*) arr_posinp(iconfig)
        enddo
     else
        nconfig=1
        allocate(arr_posinp(1:1))
        arr_posinp(1)='posinp'
     endif
     close(unit=54)
  else
     nconfig=1
     allocate(arr_posinp(1:1))
     arr_posinp(1)='posinp'
  end if

  do iconfig=1,nconfig

     !Initialize memory counting
     call memocc(0,iproc,'count','start')

     !Welcome screen
     if (iproc==0) call print_logo()

     ! Read all input files.
     call read_input_variables(iproc,trim(arr_posinp(iconfig)), &
          & "input.dft", "input.kpt", "input.geopt", inputs, atoms, rxyz)

     
     !Read absorption-calculation input variables
     !inquire for the needed file 
     !if not present, set default (no absorption calculation)
          
     inquire(file="input.abscalc",exist=exists)
     if (.not. exists) then
        if (iproc == 0) write(*,*)'ERROR: need file input.abscalc for x-ray absorber treatment.'
        if(nproc/=0)   call MPI_FINALIZE(ierr)
        stop
     end if
     call abscalc_input_variables(iproc,'input.abscalc',inputs)
     if( inputs%iat_absorber <1 .or. inputs%iat_absorber > atoms%nat) then
        if (iproc == 0) write(*,*)'ERROR: inputs%iat_absorber  must .ge. 1 and .le. number_of_atoms '
        if(nproc/=0)   call MPI_FINALIZE(ierr)
        stop
     endif


     allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fxyz,'fxyz',subname)

     call init_restart_objects(atoms,rst,subname)

     call call_abscalc(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,rst,infocode)

     if (iproc.eq.0) then
        sumx=0.d0
        sumy=0.d0
        sumz=0.d0
        write(*,'(1x,a,19x,a)') 'Final values of the Forces for each atom'
        do iat=1,atoms%nat
           write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
                iat,trim(atoms%atomnames(atoms%iatype(iat))),(fxyz(j,iat),j=1,3)
           sumx=sumx+fxyz(1,iat)
           sumy=sumy+fxyz(2,iat)
           sumz=sumz+fxyz(3,iat)
        enddo
        if (.not. inputs%gaussian_help .or. .true.) then !zero of the forces calculated
           write(*,'(1x,a)')'the sum of the forces is'
           write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
           write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
           write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
        end if
     endif

     call deallocate_atoms(atoms,subname) 

     call free_restart_objects(rst,subname)

     i_all=-product(shape(rxyz))*kind(rxyz)
     deallocate(rxyz,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz',subname)
     i_all=-product(shape(fxyz))*kind(fxyz)
     deallocate(fxyz,stat=i_stat)
     call memocc(i_stat,i_all,'fxyz',subname)

     call free_input_variables(inputs)

     !finalize memory counting
     call memocc(0,0,'count','stop')

!     call sg_end()

  enddo !loop over iconfig

  deallocate(arr_posinp)

  call MPI_FINALIZE(ierr)

end program abscalc_main
!!***


!!****f* BigDFT/call_abscalc
!! FUNCTION
!!   Routines to use abscalc as a blackbox
!! COPYRIGHT
!!   Copyright (C) 2005-2010 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
 subroutine call_abscalc(nproc,iproc,atoms,rxyz,in,energy,fxyz,rst,infocode)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(input_variables),intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  type(restart_objects), intent(inout) :: rst
  integer, intent(inout) :: infocode
  real(gp), intent(out) :: energy
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
  !local variables
  character(len=*), parameter :: subname='call_bigdft'
  character(len=40) :: comment
  integer :: i_stat,i_all,ierr,inputPsiId_orig,icycle

  !temporary interface
  interface
     subroutine abscalc(nproc,iproc,atoms,rxyz,energy,&
          psi,Glr,gaucoeffs,gbd,orbs,rxyz_old,hx_old,hy_old,hz_old,in,infocode)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nproc,iproc
       integer, intent(out) :: infocode
       real(gp), intent(inout) :: hx_old,hy_old,hz_old
       type(input_variables), intent(in) :: in
       type(locreg_descriptors), intent(inout) :: Glr
       type(atoms_data), intent(inout) :: atoms
       type(gaussian_basis), intent(inout) :: gbd
       type(orbitals_data), intent(inout) :: orbs
       real(gp), intent(out) :: energy
       real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz_old
       real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
       real(wp), dimension(:), pointer :: psi
       real(wp), dimension(:,:), pointer :: gaucoeffs
     end subroutine abscalc 
  end interface

  !put a barrier for all the processes
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !assign the verbosity of the output
  !the verbose variables is defined in module_base
  verbose=in%verbosity

  inputPsiId_orig=in%inputPsiId

  loop_cluster: do icycle=1,in%nrepmax

     if (in%inputPsiId == 0 .and. associated(rst%psi)) then
        i_all=-product(shape(rst%psi))*kind(rst%psi)
        deallocate(rst%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(rst%orbs%eval))*kind(rst%orbs%eval)
        deallocate(rst%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(rst%Glr%wfd,subname)
     end if

     if(.not. in%c_absorbtion) then 

        stop 'ERROR'
     else
        call abscalc(nproc,iproc,atoms,rxyz,energy,&
             rst%psi,rst%Glr,rst%gaucoeffs,rst%gbd,rst%orbs,&
             rst%rxyz_old,rst%hx_old,rst%hy_old,rst%hz_old,in,infocode)
        fxyz(:,:) = 0.d0
     endif

     if (in%inputPsiId==1 .and. infocode==2) then
        if (in%gaussian_help) then
           in%inputPsiId=11
        else
           in%inputPsiId=0
        end if
     else if ((in%inputPsiId==1 .or. in%inputPsiId==0) .and. infocode==1) then
        !in%inputPsiId=0 !better to diagonalise that to restart an input guess
        in%inputPsiId=1
        if(iproc==0) then
           write(*,*)&
             ' WARNING: Wavefunctions not converged after cycle',icycle
           write(*,*)' restart after diagonalisation'
        end if
        
     else if (in%inputPsiId == 0 .and. infocode==3) then
        if (iproc.eq.0) then
           write( *,'(1x,a)')'Convergence error, cannot proceed.'
           write( *,'(1x,a)')' writing positions in file posfail.xyz then exiting'
           write(comment,'(a)')'UNCONVERGED WF '
           !call wtxyz('posfail',energy,rxyz,atoms,trim(comment))

           call write_atomic_file("posfail",energy,rxyz,atoms,trim(comment))

        end if 

        i_all=-product(shape(rst%psi))*kind(rst%psi)
        deallocate(rst%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(rst%orbs%eval))*kind(rst%orbs%eval)
        deallocate(rst%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(rst%Glr%wfd,subname)
        !finalize memory counting (there are still the positions and the forces allocated)
        call memocc(0,0,'count','stop')

        if (nproc > 1) call MPI_FINALIZE(ierr)

        stop 'normal end'
     else
        exit loop_cluster
     end if

  end do loop_cluster

  !preserve the previous value
  in%inputPsiId=inputPsiId_orig

  !put a barrier for all the processes
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end subroutine call_abscalc
!!***


subroutine abscalc(nproc,iproc,atoms,rxyz,energy,&
     psi,Glr,gaucoeffs,gbd,orbs,rxyz_old,hx_old,hy_old,hz_old,in,infocode)
  ! inputPsiId = 0 : compute input guess for Psi by subspace diagonalization of atomic orbitals
  ! inputPsiId = 1 : read waves from argument psi, using n1, n2, n3, hgrid and rxyz_old
  !                  as definition of the previous system.
  ! inputPsiId = 2 : read waves from disk
  ! does an electronic structure calculation. Output is the total energy and the forces 
  ! psi, keyg, keyv and eval should be freed after use outside of the routine.
  ! infocode -> encloses some information about the status of the run
  !          =0 run succesfully succeded
  !          =1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
  !             forces may be meaningless   
  !          =2 (present only for inputPsiId=1) gnrm of the first iteration > 1 AND growing in
  !             the second iteration OR grnm 1st >2.
  !             Input wavefunctions need to be recalculated. Routine exits.
  !          =3 (present only for inputPsiId=0) gnrm > 4. SCF error. Routine exits.
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  use libxc_functionals
  use vdwcorrection, only: vdwcorrection_calculate_energy, vdwcorrection_calculate_forces, vdwcorrection_warnings
  use esatto
  implicit none
  integer, intent(in) :: nproc,iproc
  real(gp), intent(inout) :: hx_old,hy_old,hz_old
  type(input_variables), intent(in) :: in
  type(locreg_descriptors), intent(inout) :: Glr
  type(atoms_data), intent(inout) :: atoms
  type(gaussian_basis), intent(inout) :: gbd
  type(orbitals_data), intent(inout) :: orbs
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz_old
  real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
  integer, intent(out) :: infocode
  real(gp), intent(out) :: energy
  real(wp), dimension(:), pointer :: psi
  real(wp), dimension(:,:), pointer :: gaucoeffs
  !local variables
  character(len=*), parameter :: subname='abscalc'
  character(len=3) :: PSquiet
  character(len=4) :: f4
  character(len=50) :: filename
  logical :: endloop,allfiles
  integer :: ixc,ncong,idsx,ncongt,nspin,itermax,idsx_actual,idsx_actual_before
  integer :: nvirt
  integer :: nelec,ndegree_ip,j,i,k, iorb
  integer :: n1_old,n2_old,n3_old,n3d,n3p,n3pi,i3xcsh,i3s,n1,n2,n3
  integer :: ncount0,ncount1,ncount_rate,ncount_max,n1i,n2i,n3i,i03,i04
  integer :: i1,i2,i3,ind,iat,i_all,i_stat,iter,ierr,jproc,ispin,inputpsi
  real :: tcpu0,tcpu1
  real(gp), dimension(3) :: shift
  real(kind=8) :: crmult,frmult,cpmult,fpmult,gnrm_cv,rbuf,hxh,hyh,hzh,hx,hy,hz
  real(kind=8) :: peakmem,energy_old
  real(kind=8) :: eion,epot_sum,ekin_sum,eproj_sum,ehart,eexcu,alpha,gnrm
  real(kind=8) :: energybs,tel,eexcu_fake,ehart_fake,energy_min,psoffset
  real(gp) :: edisp ! Dispersion energy
  type(nonlocal_psp_descriptors) :: nlpspd
  type(communications_arrays) :: comms
  type(orbitals_data) :: orbsv
  type(gaussian_basis) :: Gvirt
  type(GPU_pointers) :: GPU

  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(kind=8), dimension(:,:), allocatable :: radii_cf,gxyz,fion
  real(gp), dimension(:,:),allocatable :: fdisp
  ! Charge density/potential,ionic potential, pkernel
  real(kind=8), dimension(:), allocatable :: pot_ion
  real(kind=8), dimension(:,:,:,:), allocatable :: rhopot,pot
  real(kind=8), dimension(:), pointer :: pkernel
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(kind=8), dimension(:), pointer :: hpsi,psit,psivirt,psidst,hpsidst
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: proj
  ! arrays for DIIS convergence accelerator
  real(kind=8), dimension(:,:,:), pointer :: ads
  ! Arrays for the symmetrisation, not used here...
  integer, dimension(:,:,:), allocatable :: irrzon
  real(dp), dimension(:,:,:), allocatable :: phnons
  
  !for xabsorber
  logical in_refinement
  integer lpot_a, ix, iy, iz
  real(gp) rpot_a, spot_a, hpot_a, espo, harmo, r, rx, ry, rz, minrx, maxrx,   minry, maxry,   minrz, maxrz, minr  
  real(gp), pointer :: radpot(:,:)
  integer radpotcount, igrid

  type(atoms_data) :: atoms_b2B
  real(gp), dimension(:,:), pointer :: rxyz_b2B
  real(gp) :: shift_b2B(3)
  integer itype, nd
  integer n1i_bB,n2i_bB,n3i_bB
  real(gp), dimension(:), pointer :: pot_bB
  real(gp) alat1_bB, alat2_bB, alat3_bB
  real(gp), dimension(:), pointer ::  intfunc_x, intfunc_y
  real(gp) factx, facty, factz
  character(len=80) :: comment
  integer idelta
  integer ix_bB, iy_bB, iz_bB
  integer maxX_B, maxY_B, maxZ_B
  integer minX_B, minY_B, minZ_B
  real(gp)  rx_bB, ry_bB, rz_bB
  integer nrange
  real(gp), pointer :: auxint(:)
  logical exists 
  integer nat_b2B

  ! ----------------------------------
  ! per monitorare il minimo del pot letto in b2B attorno al 1 atomo
  real(gp)  potx(21,3), potcoors(3,21,3)
  

  !------------------------------------------
  !copying the input variables for readability
  !this section is of course not needed
  !note that this procedure is convenient ONLY in the case of scalar variables
  !an array would have been copied, thus occupying more memory space
  !Hence WARNING: these variables are copied, in case of an update the new value should be 
  !reassigned inside the structure





  

  if (  in%potshortcut==0) then
     if(nproc>1) call MPI_Finalize(ierr)
     stop '   in%potshortcut==0 calculating spectra. Use rather box2Box option      '
  endif




  crmult=in%crmult
  frmult=in%frmult
  cpmult=in%frmult
  fpmult=in%frmult
  ixc=in%ixc
  gnrm_cv=in%gnrm_cv
  itermax=in%itermax
  ncong=in%ncong
  idsx=in%idsx
  rbuf=in%rbuf
  ncongt=in%ncongt
  nspin=in%nspin

  nvirt=in%nvirt

  hx=in%hx
  hy=in%hy
  hz=in%hz

  if (ixc < 0) then
     call libxc_functionals_init(ixc, nspin)
  end if

  !character string for quieting the Poisson solver
  if (verbose >1) then
     PSquiet='NO'
  else
     PSquiet='YES'
  end if

  if (iproc == 0) then
     write( *,'(1x,a,1x,i0)') &
          '===================== BigDFT XANE calculation =============== inputPsiId=',&
          in%inputPsiId
     call print_input_parameters(in,atoms)
  end if
  if (nproc > 1) then
     call timing(iproc,'parallel     ','IN')
  else
     call timing(iproc,'             ','IN')
  end if
  call cpu_time(tcpu0)
  call system_clock(ncount0,ncount_rate,ncount_max)


 
  if(nspin/=1 .and. nspin/=2 .and. nspin/=4) nspin=1

  ! grid spacing (same in x,y and z direction)

  if (iproc==0) then
     write( *,'(1x,a)')&
          '------------------------------------------------------------------ System Properties'
  end if

  !these routines can be regrouped in one

  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)


  call system_properties(iproc,nproc,in,atoms,orbs,radii_cf,nelec)

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space

  call system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)

  

  !variables substitution for the PSolver part
  hxh=0.5d0*hx
  hyh=0.5d0*hy
  hzh=0.5d0*hz
  n1i=Glr%d%n1i
  n2i=Glr%d%n2i
  n3i=Glr%d%n3i

  n1=Glr%d%n1
  n2=Glr%d%n2
  n3=Glr%d%n3

  ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
  call timing(iproc,'CrtDescriptors','ON')
  call createWavefunctionsDescriptors(iproc,hx,hy,hz,&
       atoms,rxyz,radii_cf,crmult,frmult,Glr,orbs)
  call timing(iproc,'CrtDescriptors','OF')
  ! Calculate all projectors, or allocate array for on-the-fly calculation
  call timing(iproc,'CrtProjectors ','ON')
  call createProjectorsArrays(iproc,n1,n2,n3,rxyz,atoms,orbs,&
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
  call timing(iproc,'CrtProjectors ','OF')


 

  !calculate the partitioning of the orbitals between the different processors
  !memory estimation
  if (iproc==0 .and. verbose > 0) then
     call MemoryEstimator(atoms%geocode,nproc,idsx,n1,n2,n3,&
          atoms%alat1,atoms%alat2,atoms%alat3,&
          hx,hy,hz,atoms%nat,atoms%ntypes,atoms%iatype,rxyz,radii_cf,crmult,frmult,&
          orbs%norb,nlpspd%nprojel,atoms%atomnames,0,in%nspin,peakmem)
  end if


  !allocate communications arrays
  !call allocate_comms(nproc,orbs,comms,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbs,comms)  

  !these arrays should be included in the comms descriptor
  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  !allocate array for the communications of the potential
  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ngatherarr,'ngatherarr',subname)
  !create the descriptors for the density and the potential
  !these descriptors should take into account the localisation regions
  call createDensPotDescriptors(iproc,nproc,atoms%geocode,'D',n1i,n2i,n3i,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  !allocate ionic potential
  if (n3pi > 0) then
     allocate(pot_ion(n1i*n2i*n3pi+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  else
     allocate(pot_ion(1+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  end if

  allocate(fion(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fion,'fion',subname)

  ! A message about dispersion forces.
  if (iproc == 0) call vdwcorrection_warnings(atoms, in)

  !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
  ndegree_ip=16 !default value 
  call createKernel(iproc,nproc,atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,pkernel,&
       quiet=PSquiet)


  call IonicEnergyandForces(iproc,nproc,atoms,hxh,hyh,hzh,in%elecfield,rxyz,eion,fion,&
       psoffset,in%nvacancy,n1,n2,n3,n1i,n2i,n3i,i3s+i3xcsh,n3pi,pot_ion,pkernel)

  call createIonicPotential(atoms%geocode,iproc,nproc,atoms,rxyz,hxh,hyh,hzh,&
       in%elecfield,n1,n2,n3,n3pi,i3s+i3xcsh,n1i,n2i,n3i,pkernel,pot_ion,psoffset,in%nvacancy,&
       in%correct_offset)

  !this can be inserted inside the IonicEnergyandForces routine
  !(after insertion of the non-regression test)
  call vdwcorrection_calculate_energy(edisp,rxyz,atoms,in,iproc)

  allocate(fdisp(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fdisp,'fdisp',subname)
  !this can be inserted inside the IonicEnergyandForces routine
  call vdwcorrection_calculate_forces(fdisp,rxyz,atoms,in) 

  !Allocate Charge density, Potential in real space
  if (n3d >0) then
     allocate(rhopot(n1i,n2i,n3d,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  else
     allocate(rhopot(1,1,1,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  end if


  !check the communication distribution
  !call check_communications(iproc,nproc,orbs,Glr,comms)


  if(in%potshortcut/=2) then



     inputpsi=in%inputPsiId


     nspin=in%nspin
     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(iproc,nproc,atoms,&
          orbs,orbsv,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,pot_ion,&
          nlpspd,proj,pkernel,ixc,psi,hpsi,psit,psivirt,Gvirt,&
          nscatterarr,ngatherarr,nspin, in%potshortcut, -1, irrzon, phnons)


     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)




  end if
   
  if (nproc > 1  ) then
     i_all=-product(shape(hpsi))*kind(hpsi)
     deallocate(hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi',subname)
  endif
  
  if (nproc > 1  ) then
     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  else
     nullify(psit)
  end if
  


  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)
       
  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'kernel',subname)




  ! needs something to let to  bigdft to deallocate
  allocate(psi(2+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)


  allocate(orbs%eval(2+ndebug),stat=i_stat)
  call memocc(i_stat, orbs%eval,'eval',subname)


  if (in%c_absorbtion ) then

!!$
!!$     rhopot(10,9,8+i3xcsh,1)=100.0

     if (abs(in%output_grid)==2) then
        if (in%output_grid==2) then
          if (iproc == 0) write(*,*) 'writing local_potential.pot'
           call plot_density(atoms%geocode,'local_potentialb2B.pot',iproc,nproc,&
                n1,n2,n3,n1i,n2i,n3i,n3p,&
                atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rhopot(1,1,1+i3xcsh,1))
        else
           call plot_density_cube(atoms%geocode,'local_potentialb2B',iproc,nproc,&
                n1,n2,n3,n1i,n2i,n3i,n3p,&
                in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rhopot(1,1,1+i3xcsh,1))
        endif
     end if

!!$     call  read_potfile4b2B("local_potential.pot",n1i_bB,n2i_bB,n3i_bB, pot_bB, alat1_bB, alat2_bB, alat3_bB)
!!$     print *, pot_bB(10  + (9-1)*n1i_bB  + (8-1)*n1i_bB*n2i_bB)
!!$     stop
     
     if(iproc==0) print *, " going to calculate spectra "

     
     if(in%potshortcut==2) then

        inquire(file="b2B_xanes.cube",exist=exists)
        if(exists) then

           call read_density_cube("b2B_xanes", n1i_bB,n2i_bB,n3i_bB, 1 , hx_old ,hy_old ,hz_old , nat_b2B, rxyz_b2B, pot_bB )
           hx_old=hx_old*2
           hy_old=hy_old*2
           hz_old=hz_old*2


           if( atoms%nat/= nat_b2B ) then
              if(iproc==0) write(*,*)  "   b2B_xanes cube  is not compatible with actual positions" 
              if(nproc>1) call MPI_Finalize(ierr)
              stop '      b2B_xanes cube  is not compatible with actual positions          '
           end if

        else
           print  *, " reading atomic positions from file ","b2B_xanes.xyz"
           call read_atomic_file("b2B_xanes.xyz",iproc, atoms_b2B, rxyz_b2B )
           if( atoms%nat/=atoms_b2B%nat) then
              if(iproc==0) write(*,*)  "   b2B_xanes.xyz  is not compatible with actual positions" 
              if(nproc>1) call MPI_Finalize(ierr)
              stop '      b2B_xanes.xyz  is not compatible with actual positions          '
           end if
           
           print  *, " reading potential from file ","b2B_xanes.pot"
           call  read_potfile4b2B("b2B_xanes.pot",n1i_bB,n2i_bB,n3i_bB, pot_bB, alat1_bB, alat2_bB, alat3_bB)
           print  *, " reading OK "
           
           
           if( atoms_b2B%geocode/='F') then
              hx_old = 2*alat1_bB / (n1i_bB)
              hy_old = 2*alat2_bB / (n2i_bB)
              hz_old = 2*alat3_bB / (n3i_bB)
           else
              hx_old = 2*alat1_bB / (n1i_bB-2)
              hy_old = 2*alat2_bB / (n2i_bB-2)
              hz_old = 2*alat3_bB / (n3i_bB-2)
           endif

           call deallocate_atoms(atoms_b2B,subname) 

        endif

        if(  n1i_bB > n1i .or.  n2i_bB > n2i .or.   n3i_bB > n3i    ) then
           if(nproc>1) call MPI_Finalize(ierr)
           stop '  b2B potential must be defined on a smaller ( in number of points  ) source grid then the target grid    '
        endif

        do j=1,3
           shift_b2B(j) = rxyz(j,1) - rxyz_b2B(j,1) 
           do iat=2,atoms%nat
              if( abs(shift_b2B(j) - (rxyz(j,iat) - rxyz_b2B(j,iat) )  )>1.0e-4 ) then
                 if(iproc==0) write(*,*)  "   b2B_xanes.xyz  positions are not compatible with actual positions" 
                 if(nproc>1) call MPI_Finalize(ierr)
                 stop '      b2B_xanes.xyz positions are not compatible with actual positions          '
              end if
           enddo
        enddo

        print *, "SHIFT " , shift_b2B
        
        ! passing from an old grid  x  to a new grid one
        !    ((x-1)*hx_old/2.0+shift_x)/(hx/2.0) +1 
        ! inverse change
        !    ((x-1)*hx/2.0-shift_x)/(hx_old/2.0) +1 
        

        itype=16
        nd=2**20

        allocate(  intfunc_x(0:nd+ndebug),stat=i_stat )
        call memocc(i_stat,intfunc_x,'intfunc_x',subname)
        allocate( intfunc_y(0:nd+ndebug) ,stat=i_stat )
        call memocc(i_stat,intfunc_y,'intfunc_y',subname)

        print *, " scaling function for interpolation "

        call scaling_function4b2B(itype,nd,nrange,intfunc_x,intfunc_y)  ! intervallo di 32 con 2**20 punti
        if( abs(intfunc_y(nd/2)-1)>1.0e-10 ) then
           stop " wrong scaling function 4b2B: not a centered one "
        endif

        i_all=-product(shape(intfunc_x))*kind(intfunc_x)
        deallocate(intfunc_x,stat=i_stat)
        call memocc(i_stat,i_all,'intfunc_x',subname)
        

        ! how to acceed to a  function using a  
        ! x variable given in units of the old  grids intfunc_y(( x+16) *2**15)

        rhopot=0.0_gp


        do j=-10,10
           do k=1,3
              potcoors(:,j+11,k)= rxyz_b2B(:,1)
              potcoors(k,j+11,k) = potcoors(k,j+11,k)+j*0.1
           enddo
        enddo
        potx=0.0

        do iz_bB = 1,n3i_bB
           do iy_bB=1,n2i_bB
              do ix_bB=1,n1i_bB
                 rhopot(ix_bB,iy_bB,iz_bB +i3xcsh,1) =  pot_bB(ix_bB  + (iy_bB-1)*n1i_bB  + (iz_bB-1)*n1i_bB*n2i_bB)
                 
                 do j=-10,10
                    do k=1,3
                       rx_bB = hx_old*(ix_bB-1)           /2.0   - potcoors(1,j+11,k)
                       ry_bB = hy_old*(iy_bB-1)           /2.0   - potcoors(2,j+11,k)
                       rz_bB = hz_old*(iz_bB-1)           /2.0   - potcoors(3,j+11,k)

                       idelta = (NINT((rx_bB)*2**15/(hx_old/2)))  + nd/2
                       if(idelta<nd .and. idelta>0 ) then 
                          factx = intfunc_y(idelta)
                          
                          idelta =  abs(NINT((ry_bB)*2**15/(hy_old/2)) ) + nd/2
                          if(idelta<nd) then 
                             facty = intfunc_y(idelta)
                             
                             idelta = abs( NINT((rz_bB)*2**15/(hz_old/2)) ) + nd/2
                             if(idelta<nd) then 
                                factz = intfunc_y(idelta)
                                
                                potx(j+11 ,k) =potx(j+11,k)+factx*facty*factz*rhopot(ix_bB,iy_bB,iz_bB +i3xcsh,1)
                             end if
                          end if
                       end if
                    end do
                 end do
              enddo
           enddo
        enddo
        

        open(unit=22,file='potx.dat', status='unknown')
        do j=1,21
           write(22,*)  (potx(j,k), k=1,3)
        enddo
        close(unit=22)
        

        

        allocate(auxint(n1i+n2i+n3i+ndebug),stat=i_stat)
        call memocc(i_stat,auxint,'auxint',subname)

        do iz_bB = 1,n3i_bB-1
           do iy_bB=1,n2i_bB
              auxint = 0.0_gp
              do ix_bB=1,n1i_bB
                 rx_bB = hx_old*(ix_bB-1)           /2.0   +  shift_b2B(1)
                 minX_B  =  max(1,NINT((rx_bB -8*hx_old/2)/(hx/2.0)))
                 maxX_B  =  min(n1i,NINT((rx_bB +8*hx_old/2)/(hx/2.0)))
                 do ix= minX_B , maxX_B 
                    rx = hx*(ix-1  )/2.0  
                    idelta = NINT((rx-rx_bB)*2**15/(hx_old/2))  
                    factx = intfunc_y(nd/2+idelta)
!!$                    print *, rx, rx_bB, ix, ix_bB      , factx
                    auxint(ix) = auxint(ix) + &
                         factx * rhopot(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                 enddo
!!$                 print *, auxint(ix_bB) ,  rhopot(ix_bB,iy_bB,iz_bB+i3xcsh,1)
              enddo
              rhopot(:,iy_bB,iz_bB+i3xcsh,1)=auxint(1:n1i)
           enddo
        enddo

        do iz_bB = 1,n3i_bB-1
           do ix_bB=1,n1i
              auxint = 0.0_gp
              do iy_bB=1,n2i_bB
                 ry_bB = hy_old*(iy_bB-1)/2.0   +  shift_b2B(2)
                 minY_B  =  max(1  ,NINT((ry_bB -8*hy_old/2)/(hy/2.0)))
                 maxY_B  =  min(n2i,NINT((ry_bB +8*hy_old/2)/(hy/2.0)))
                 do iy= minY_B , maxY_B 
                    ry = hy*(iy-1  )/2.0  
                    idelta = NINT((ry-ry_bB)*2**15/(hy_old/2))
                    facty = intfunc_y(nd/2+idelta)
                    auxint(iy) = auxint(iy) + &
                         facty * rhopot(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                 enddo
              enddo
              rhopot(ix_bB ,:,iz_bB+i3xcsh,1)=auxint(1:n2i)
           enddo
        enddo
 
        do ix_bB=1,n1i
           do iy_bB=1,n2i
              auxint = 0.0_gp
              do iz_bB = 1,n3i_bB-1
                 rz_bB = hz_old*(iz_bB-1)           /2.0   +  shift_b2B(3)

                 minZ_B  =  max(1  ,  NINT((rz_bB -8*hz_old/2)/(hz/2.0)))
                 maxZ_B  =  min(n3i-i3xcsh , NINT((rz_bB +8*hz_old/2)/(hz/2.0)))

                 do iz= minZ_B , maxZ_B 
                    rz = hz*(iz-1  )/2.0  
                    idelta = NINT((rz-rz_bB)*2**15/(hz_old/2.0))     
                    factz = intfunc_y(nd/2+idelta)
                    auxint(iz+i3xcsh) = auxint(iz+i3xcsh) + &
                         factz * rhopot(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                 enddo
              enddo
              rhopot(ix_bB ,iy_bB, : ,1)=auxint(1:n3i)
           enddo
        enddo

        if (iproc == 0) write(*,*) 'writing NEW local_potential.pot'



        call plot_density(atoms%geocode,'local_potentialb2BNEW.pot',iproc,nproc,&
             n1,n2,n3,n1i,n2i,n3i,n3p,&
             atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rhopot(1,1,1+i3xcsh,1))





        i_all=-product(shape(auxint))*kind(auxint)
        deallocate(auxint,stat=i_stat)
        call memocc(i_stat,i_all,'auxint',subname)
 
        i_all=-product(shape(pot_bB))*kind(pot_bB)
        deallocate(pot_bB,stat=i_stat)
        call memocc(i_stat,i_all,'rho',subname)
        
        i_all=-product(shape(intfunc_y))*kind(intfunc_y)
        deallocate(intfunc_y,stat=i_stat)
        call memocc(i_stat,i_all,'intfunc_y',subname)
        print *," exiting b2B"
  

        i_all=-product(shape(rxyz_b2B))*kind(rxyz_b2B)
        deallocate(rxyz_b2B,stat=i_stat)
        call memocc(i_stat,i_all,'rxyz',subname)
              
     endif
     
     if(in%abscalc_alterpot) then
        ! Attention :  modification of the  potential for the  
        ! exactly resolvable case 

        lpot_a=1
        rpot_a = 6.0d0
        spot_a = 1.0d0
        hpot_a = 3.0d0


        allocate(radpot(30000 ,2+ndebug ))
        radpotcount=30000

        open(unit=22,file='pot.dat', status='old')
        do igrid=1, radpotcount
           read(22,*)  radpot(igrid ,1 ),  radpot(igrid , 2 )
        enddo
        close(unit=22)
        
        minr=1000.0
        do ix=1,n1i
           do iy=1,n2i
              do iz = 1,n3p
                 rx = hx*(ix-1)           /2.0  -  rxyz(1,in%iat_absorber )
                 ry = hy*(iy-1)           /2.0  -  rxyz(2,in%iat_absorber )
                 rz = hz*(iz-1 +i3xcsh + i3s -1 )/2.0  -  rxyz(3,in%iat_absorber )
                 
                 r  = sqrt( rx*rx+ry*ry+rz*rz)
                 
                 if(r>2.5) then
                    if( r>29) then
                       rhopot(ix,iy,iz+i3xcsh,1)=0.0
                    else
                       igrid = binary_search( r, radpot, radpotcount )
                       rhopot(ix,iy,iz+i3xcsh,1) = &
                            ( radpot(igrid,2)*(radpot(igrid+1,1)-R) + radpot(igrid+1,2)*(R-radpot(igrid,1)) )/&
                            ( radpot(igrid+1,1) -radpot(igrid,1) )
                    endif
                 endif
                 
                 if(r<minr) minr=r
                 
                 if( r.ne.0.0) then
                    harmo = rz/r *sqrt(3.0/4.0/3.1415926535)
                 else
                    harmo=0.0_gp
                 endif
                 
                 espo  = ((r-rpot_a)**2)/spot_a/spot_a/2.0
                 if(espo<100) then
                    rhopot(ix,iy,iz+i3xcsh,1) = rhopot(ix,iy,iz+i3xcsh,1) +  hpot_a * exp(-espo) *harmo
                 endif
              enddo
           enddo
        enddo
        
     end if
     infocode=0
 
     if(in%iabscalc_type==2) then
        call xabs_lanczos(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
             rhopot(1,1,1+i3xcsh,1) ,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU &
             , in%iat_absorber  , .false., orbs%norb,   psit , orbs%eval , in )
        
     else
        call xabs_chebychev(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
             rhopot(1,1,1+i3xcsh,1) ,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU &
             , in%iat_absorber, in)
        
     endif
     
  end if

  !    No tail calculation
  if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call deallocate_before_exiting


contains


  !routine which deallocate the pointers and the arrays before exiting 
  subroutine deallocate_before_exiting
    
    !when this condition is verified we are in the middle of the SCF cycle

    !! if (infocode /=0 .and. infocode /=1) then
    if (.true.) then
       
!!$       if (idsx_actual > 0) then
!!$          i_all=-product(shape(psidst))*kind(psidst)
!!$          deallocate(psidst,stat=i_stat)
!!$          call memocc(i_stat,i_all,'psidst',subname)
!!$          i_all=-product(shape(hpsidst))*kind(hpsidst)
!!$          deallocate(hpsidst,stat=i_stat)
!!$          call memocc(i_stat,i_all,'hpsidst',subname)
!!$          i_all=-product(shape(ads))*kind(ads)
!!$          deallocate(ads,stat=i_stat)
!!$          call memocc(i_stat,i_all,'ads',subname)
!!$       end if
       
!!$       if (nproc > 1) then
!!$          i_all=-product(shape(psit))*kind(psit)
!!$          deallocate(psit,stat=i_stat)
!!$          call memocc(i_stat,i_all,'psit',subname)
!!$       end if
!!$       
!!$       i_all=-product(shape(hpsi))*kind(hpsi)
!!$       deallocate(hpsi,stat=i_stat)
!!$       call memocc(i_stat,i_all,'hpsi',subname)
       

       
       
!!$       i_all=-product(shape(pot_ion))*kind(pot_ion)
!!$       deallocate(pot_ion,stat=i_stat)
!!$       call memocc(i_stat,i_all,'pot_ion',subname)
!!$       
!!$       i_all=-product(shape(pkernel))*kind(pkernel)
!!$       deallocate(pkernel,stat=i_stat)
!!$       call memocc(i_stat,i_all,'pkernel',subname)
!!$

!!$       if (in%read_ref_den) then
!!$          i_all=-product(shape(pkernel_ref))*kind(pkernel_ref)
!!$          deallocate(pkernel_ref,stat=i_stat)
!!$          call memocc(i_stat,i_all,'pkernel_ref',subname)
!!$       end if
       
       ! calc_tail false
       i_all=-product(shape(rhopot))*kind(rhopot)
       deallocate(rhopot,stat=i_stat)
       call memocc(i_stat,i_all,'rhopot',subname)

!!$       if (in%read_ref_den) then
!!$          i_all=-product(shape(rhoref))*kind(rhoref)
!!$          deallocate(rhoref,stat=i_stat)
!!$          call memocc(i_stat,i_all,'rhoref',subname)
!!$       end if
       
       i_all=-product(shape(nscatterarr))*kind(nscatterarr)
       deallocate(nscatterarr,stat=i_stat)
       call memocc(i_stat,i_all,'nscatterarr',subname)


       i_all=-product(shape(ngatherarr))*kind(ngatherarr)
       deallocate(ngatherarr,stat=i_stat)
       call memocc(i_stat,i_all,'ngatherarr',subname)
       
       i_all=-product(shape(fion))*kind(fion)
       deallocate(fion,stat=i_stat)
       call memocc(i_stat,i_all,'fion',subname)

       i_all=-product(shape(fdisp))*kind(fdisp)
       deallocate(fdisp,stat=i_stat)
       call memocc(i_stat,i_all,'fdisp',subname)
       
       
    end if

    !deallocate wavefunction for virtual orbitals
    !if it is the case
    if (in%nvirt > 0) then
       !call deallocate_gwf(Gvirt,subname)
       i_all=-product(shape(psivirt))*kind(psivirt)
       deallocate(psivirt,stat=i_stat)
       call memocc(i_stat,i_all,'psivirt',subname)
    end if
    
    if (atoms%geocode == 'F') then
       call deallocate_bounds(Glr%bounds,subname)
    end if
    
    if (atoms%geocode == 'P' .and. Glr%hybrid_on) then 
       
       i_all=-product(shape(Glr%bounds%kb%ibxy_f))*kind(Glr%bounds%kb%ibxy_f)
       deallocate(Glr%bounds%kb%ibxy_f,stat=i_stat)
       call memocc(i_stat,i_all,'Glr%bounds%kb%ibxy_f',subname)
       
       i_all=-product(shape(Glr%bounds%kb%ibxz_f))*kind(Glr%bounds%kb%ibxz_f)
       deallocate(Glr%bounds%kb%ibxz_f,stat=i_stat)
       call memocc(i_stat,i_all,'Glr%bounds%kb%ibxz_f',subname)
       
       i_all=-product(shape(Glr%bounds%kb%ibyz_f))*kind(Glr%bounds%kb%ibyz_f)
       deallocate(Glr%bounds%kb%ibyz_f,stat=i_stat)
       call memocc(i_stat,i_all,'Glr%bounds%kb%ibyz_f',subname)
       
       i_all=-product(shape(Glr%bounds%sb%ibxy_ff))*kind(Glr%bounds%sb%ibxy_ff)
       deallocate(Glr%bounds%sb%ibxy_ff,stat=i_stat)
       call memocc(i_stat,i_all,'ibxy_ff',subname)
       i_all=-product(shape(Glr%bounds%sb%ibzzx_f))*kind(Glr%bounds%sb%ibzzx_f)
       deallocate(Glr%bounds%sb%ibzzx_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibzzx_f',subname)
       i_all=-product(shape(Glr%bounds%sb%ibyyzz_f))*kind(Glr%bounds%sb%ibyyzz_f)
       deallocate(Glr%bounds%sb%ibyyzz_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibyyzz_f',subname)
       
       i_all=-product(shape(Glr%bounds%gb%ibyz_ff))*kind(Glr%bounds%gb%ibyz_ff)
       deallocate(Glr%bounds%gb%ibyz_ff,stat=i_stat)
       call memocc(i_stat,i_all,'ibyz_ff',subname)
       
       i_all=-product(shape(Glr%bounds%gb%ibzxx_f))*kind(Glr%bounds%gb%ibzxx_f)
       deallocate(Glr%bounds%gb%ibzxx_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibzxx_f',subname)
       
       i_all=-product(shape(Glr%bounds%gb%ibxxyy_f))*kind(Glr%bounds%gb%ibxxyy_f)
       deallocate(Glr%bounds%gb%ibxxyy_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibxxyy_f',subname)
    endif
    
    
    call deallocate_comms(comms,subname)

  
    call deallocate_orbs(orbs,subname)
    

    !semicores useful only for the input guess
    i_all=-product(shape(atoms%iasctype))*kind(atoms%iasctype)
    deallocate(atoms%iasctype,stat=i_stat)
    call memocc(i_stat,i_all,'iasctype',subname)
    i_all=-product(shape(atoms%aocc))*kind(atoms%aocc)
    deallocate(atoms%aocc,stat=i_stat)
    call memocc(i_stat,i_all,'aocc',subname)
    i_all=-product(shape(atoms%nzatom))*kind(atoms%nzatom)
    deallocate(atoms%nzatom,stat=i_stat)
    call memocc(i_stat,i_all,'nzatom',subname)
    i_all=-product(shape(nlpspd%nboxp_c))*kind(nlpspd%nboxp_c)
    deallocate(nlpspd%nboxp_c,stat=i_stat)
    call memocc(i_stat,i_all,'nboxp_c',subname)
    i_all=-product(shape(nlpspd%nboxp_f))*kind(nlpspd%nboxp_f)
    deallocate(nlpspd%nboxp_f,stat=i_stat)
    call memocc(i_stat,i_all,'nboxp_f',subname)
    i_all=-product(shape(nlpspd%keyg_p))*kind(nlpspd%keyg_p)
    deallocate(nlpspd%keyg_p,stat=i_stat)
    call memocc(i_stat,i_all,'keyg_p',subname)
    i_all=-product(shape(nlpspd%keyv_p))*kind(nlpspd%keyv_p)
    deallocate(nlpspd%keyv_p,stat=i_stat)
    call memocc(i_stat,i_all,'keyv_p',subname)
    i_all=-product(shape(nlpspd%nvctr_p))*kind(nlpspd%nvctr_p)
    deallocate(nlpspd%nvctr_p,stat=i_stat)
    call memocc(i_stat,i_all,'nvctr_p',subname)
    i_all=-product(shape(nlpspd%nseg_p))*kind(nlpspd%nseg_p)
    deallocate(nlpspd%nseg_p,stat=i_stat)
    call memocc(i_stat,i_all,'nseg_p',subname)

    i_all=-product(shape(proj))*kind(proj)
    deallocate(proj,stat=i_stat)
    call memocc(i_stat,i_all,'proj',subname)

    i_all=-product(shape(atoms%psppar))*kind(atoms%psppar)
    deallocate(atoms%psppar,stat=i_stat)
    call memocc(i_stat,i_all,'psppar',subname)
    i_all=-product(shape(atoms%nelpsp))*kind(atoms%nelpsp)
    deallocate(atoms%nelpsp,stat=i_stat)
    call memocc(i_stat,i_all,'nelpsp',subname)
    i_all=-product(shape(radii_cf))*kind(radii_cf)
    deallocate(radii_cf,stat=i_stat)
    call memocc(i_stat,i_all,'radii_cf',subname)
    i_all=-product(shape(atoms%npspcode))*kind(atoms%npspcode)
    deallocate(atoms%npspcode,stat=i_stat)
    call memocc(i_stat,i_all,'npspcode',subname)

    ! Free the libXC stuff if necessary.
    if (ixc < 0) then
       call libxc_functionals_end()
    end if

    !end of wavefunction minimisation
    call timing(iproc,'LAST','PR')
    call timing(iproc,'              ','RE')
    call cpu_time(tcpu1)
    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)
    if (iproc == 0) &
         write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', iproc,tel,tcpu1-tcpu0

 

  end subroutine deallocate_before_exiting

END SUBROUTINE abscalc
!!***

