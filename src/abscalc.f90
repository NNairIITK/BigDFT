!> @file
!! Main file for XANES calculation
!! @author Copyright (C) 2009-2011 ESRF
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!>  Main program for XANES calculation (absorption calculation)
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
  integer :: iproc,nproc,iat,j,i_stat,i_all,ierr,infocode
  real(gp) :: etot,sumx,sumy,sumz
  logical :: exist_list
  !input variables
  type(atoms_data) :: atoms
  type(input_variables) :: inputs
  type(restart_objects) :: rst
  character(len=50), dimension(:), allocatable :: arr_posinp
  ! atomic coordinates, forces
  real(gp), dimension(:,:), allocatable :: fxyz
  real(gp), dimension(:,:), pointer :: rxyz
  integer :: iconfig,nconfig
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
        allocate(arr_posinp(1:nconfig),stat=i_stat)
        do iconfig=1,nconfig
           read(54,*) arr_posinp(iconfig)
        enddo
     else
        nconfig=1
        allocate(arr_posinp(1:1),stat=i_stat)
        arr_posinp(1)='posinp'
     endif
     close(unit=54)
  else
     nconfig=1
     allocate(arr_posinp(1:1),stat=i_stat)
     arr_posinp(1)='posinp'
  end if

  do iconfig=1,nconfig

     !Welcome screen
     if (iproc==0) call print_logo()

     ! Read all input files.
     call read_input_variables(iproc,trim(arr_posinp(iconfig)), &
          & "input.dft", "input.kpt","input.mix","input.geopt", "input.perf", inputs, atoms, rxyz)

     !Initialize memory counting
     !call memocc(0,iproc,'count','start')
     
     !Read absorption-calculation input variables
     !inquire for the needed file 
     !if not present, set default (no absorption calculation)
          
     inquire(file="input.abscalc",exist=exists)
     if (.not. exists) then
        if (iproc == 0) write(*,*) 'ERROR: need file input.abscalc for x-ray absorber treatment.'
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

     call init_restart_objects(iproc,inputs%iacceleration,atoms,rst,subname)

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

  i_all=-product(shape(arr_posinp))*kind(arr_posinp)
  deallocate(arr_posinp)
  call memocc(i_stat,i_all,'arr_posinp',subname)

  call MPI_FINALIZE(ierr)

end program abscalc_main


!>   Routines to use abscalc as a blackbox
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
     subroutine abscalc(nproc,iproc,atoms,rxyz,&
          psi,Glr,orbs,hx_old,hy_old,hz_old,in,GPU,infocode)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nproc,iproc
       integer, intent(out) :: infocode
       real(gp), intent(inout) :: hx_old,hy_old,hz_old
       type(input_variables), intent(in) :: in
       type(locreg_descriptors), intent(inout) :: Glr
       type(atoms_data), intent(inout) :: atoms
       type(orbitals_data), intent(inout) :: orbs
       type(GPU_pointers), intent(inout) :: GPU
       real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
       real(wp), dimension(:), pointer :: psi
     END SUBROUTINE abscalc 
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
        nullify(rst%orbs%eval)

        call deallocate_wfd(rst%Glr%wfd,subname)
     end if

     if(.not. in%c_absorbtion) then 

        stop 'ERROR'
     else
        call abscalc(nproc,iproc,atoms,rxyz,&
             rst%psi,rst%Glr,rst%orbs,&
             rst%hx_old,rst%hy_old,rst%hz_old,in,rst%GPU,infocode)
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
        nullify(rst%orbs%eval)

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

END SUBROUTINE call_abscalc


!>   Absorption (XANES) calculation
!!   @param inputPsiId = 
!!          - 0 : compute input guess for Psi by subspace diagonalization of atomic orbitals
!!          - 1 : read waves from argument psi, using n1, n2, n3, hgrid
!!                as definition of the previous system.
!!          - 2 : read waves from disk
!!                does an electronic structure calculation. Output is the total energy and the forces 
!!   @param psi, keyg, keyv and eval should be freed after use outside of the routine.
!!   @param infocode -> encloses some information about the status of the run
!!          - 0 run succesfully succeded
!!          - 1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
!!               forces may be meaningless   
!!          - 2 (present only for inputPsiId=1) gnrm of the first iteration > 1 AND growing in
!!               the second iteration OR grnm 1st >2.
!!               Input wavefunctions need to be recalculated. Routine exits.
!!          - 3 (present only for inputPsiId=0) gnrm > 4. SCF error. Routine exits.
!!
subroutine abscalc(nproc,iproc,atoms,rxyz,&
     psi,Glr,orbs,hx_old,hy_old,hz_old,in,GPU,infocode)
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
  type(orbitals_data), intent(inout) :: orbs
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
  integer, intent(out) :: infocode
  real(wp), dimension(:), pointer :: psi
  !local variables
  character(len=*), parameter :: subname='abscalc'
  character(len=3) :: PSquiet
  integer :: ixc,ncong,idsx,ncongt,nspin,itermax
  integer :: nvirt
  integer :: nelec,ndegree_ip,j
  integer :: n3d,n3p,n3pi,i3xcsh,i3s,n1,n2,n3
  integer :: ncount0,ncount1,ncount_rate,ncount_max,n1i,n2i,n3i
  integer :: iat,i_all,i_stat,ierr,inputpsi
  real :: tcpu0,tcpu1
  real(gp), dimension(3) :: shift
  real(kind=8) :: crmult,frmult,cpmult,fpmult,gnrm_cv,rbuf,hxh,hyh,hzh,hx,hy,hz
  real(kind=8) :: peakmem
  real(kind=8) :: eion,epot_sum,ekin_sum,eproj_sum
  real(kind=8) :: tel,psoffset
  real(gp) :: edisp ! Dispersion energy
  type(nonlocal_psp_descriptors) :: nlpspd
  type(communications_arrays) :: comms
  type(gaussian_basis) :: Gvirt

  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(kind=8), dimension(:,:), allocatable :: radii_cf,fion
  !real(kind=8), dimension(:,:), allocatable :: gxyz
  real(gp), dimension(:,:),allocatable :: fdisp
  ! Charge density/potential,ionic potential, pkernel
  real(kind=8), dimension(:), allocatable :: pot_ion

  real(kind=8), dimension(:,:,:,:), allocatable, target :: rhopot
  real(kind=8), dimension(:,:,:,:), pointer ::  rhopottmp, rhopotExtra, rhoXanes, rhotarget
  integer b2Bcounter, b2BN
  character(len=100) filename
  real(kind=8), dimension(:), pointer :: pkernel

  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(kind=8), dimension(:), pointer :: hpsi,psit,psivirt,rhocore
  !real(kind=8), dimension(:), pointer :: psidst,hpsidst
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: proj
  ! arrays for DIIS convergence accelerator
  !real(kind=8), dimension(:,:,:), pointer :: ads
  ! Arrays for the symmetrisation, not used here...
  integer, dimension(:,:,:), allocatable :: irrzon
  real(dp), dimension(:,:,:), allocatable :: phnons
  
  !for xabsorber
  integer :: lpot_a, ix, iy, iz , ixnl, iynl, iznl
  real(gp) :: rpot_a,spot_a,hpot_a,espo,harmo,r,rx,ry,rz,minr  
  real(gp), pointer :: radpot(:,:)
  integer :: radpotcount, igrid

  type(atoms_data) :: atoms_b2B
  real(gp), dimension(:,:), pointer :: rxyz_b2B
  integer, dimension(:), pointer :: iatype_b2B, znucl_b2B
  real(gp) :: shift_b2B(3)
  integer itype, nd
  integer n1i_bB,n2i_bB,n3i_bB
  real(gp), dimension(:,:), pointer :: pot_bB
  real(gp) alat1_bB, alat2_bB, alat3_bB
  real(gp), dimension(:), pointer ::  intfunc_x, intfunc_y
  real(gp) factx, facty, factz
  integer :: idelta
  integer :: ix_bB, iy_bB, iz_bB
  integer :: maxX_B, maxY_B, maxZ_B
  integer :: minX_B, minY_B, minZ_B
  real(gp) :: rx_bB, ry_bB, rz_bB
  integer :: nrange
  real(gp), pointer :: auxint(:)

  logical exists 
  integer nat_b2B
  integer Nreplicas, ireplica, replicaoffset
  real(gp) dumvect3D(3)
  real(gp) shiftdiff
  real(gp) potmodified_maxr, potmodified_shift


  type(atoms_data) :: atoms_clone
  integer :: nsp, nspinor, noncoll
  integer, parameter :: nelecmax=32,nmax=6,lmax=4
  integer, parameter :: noccmax=2


  if (in%potshortcut==0) then
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
     call print_dft_parameters(in,atoms)
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
       atoms,rxyz,radii_cf,crmult,frmult,Glr)
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
          orbs%norb,orbs%nspinor,orbs%nkpts,nlpspd%nprojel,atoms%atomnames,0,&
          in%nspin,in%itrpmax,in%iscf,peakmem)
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
  print *, " allocate ionic potential " 
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

  print *, " IonicEnergyandForces  " 

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

  if( iand( in%potshortcut,16)>0) then
     allocate(rhoXanes(n1i,n2i,n3i,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhoXanes,'rhoXanes',subname)
  else
     allocate(rhoXanes(1,1,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,rhoXanes,'rhoXanes',subname)     
  endif
  

  nullify(rhocore)
  !check the communication distribution
  !call check_communications(iproc,nproc,orbs,Glr,comms)


  if( iand( in%potshortcut,4)  .gt. 0 ) then
     if (n3d >0) then
        allocate(rhopotExtra(n1i,n2i,n3d,in%nspin+ndebug),stat=i_stat)
        call memocc(i_stat,rhopotExtra,'rhopotExtra',subname)
     else
        allocate(rhopotExtra(1,1,1,in%nspin+ndebug),stat=i_stat)
        call memocc(i_stat,rhopotExtra,'rhopotExtra',subname)
     end if

     atoms_clone = atoms
     nullify(atoms_clone%aocc)
     nullify(atoms_clone%iasctype)
     
     
     allocate(atoms_clone%aocc(lbound(atoms%aocc,1 ):ubound(atoms%aocc,1),&
          lbound(atoms%aocc,2):ubound(atoms%aocc,2)),stat=i_stat)
     call memocc(i_stat,atoms%aocc,'atoms_clone%aocc',subname)

     allocate(atoms_clone%iasctype(lbound(atoms%iasctype,1 ):ubound(atoms%iasctype,1)),stat=i_stat)
     call memocc(i_stat,atoms%iasctype,'atoms_clone%iasctype',subname)

  
     atoms_clone%aocc=0.0_gp
     atoms_clone%iasctype=0


     read(in%extraOrbital,*,iostat=ierr)iat
     !control the spin
     select case(in%nspin)
     case(1)
        nsp=1
        nspinor=1
        noncoll=1
     case(2)
        nsp=2
        nspinor=1
        noncoll=1
     case(4)
        nsp=1
        nspinor=4
        noncoll=2
     case default
        write(*,*)' ERROR: nspin not valid:',nspin
        stop
     end select

     
     print *, " Going to create extra potential for orbital "
     print *, in%extraOrbital
     print *, "using hard-coded parameters "
     print *, "noccmax, nelecmax,lmax ", noccmax, nelecmax,lmax


     call read_eleconf(in%extraOrbital ,nsp,nspinor,noccmax, nelecmax,lmax, &
     atoms_clone%aocc(1,iat), atoms_clone%iasctype(iat))

     
     nspin=in%nspin
     
     !Fake allocations
     allocate(irrzon(1,2,1+ndebug),stat=i_stat)
     call memocc(i_stat,irrzon,'irrzon',subname)
     allocate(phnons(2,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,phnons,'phnons',subname)
     
     !-- calculate input guess from non-diagonalised of LCAO basis (written in wavelets)
     !-- if spectra calculation is energy dependent  input_wf_diag will write
     !    the density to the file electronic_density.cube
     !   To tell  input_wf_diag to do this ne must set the 5th bit on in in%potshortcut
     if( iand( in%potshortcut,16)>0) then
        if( in%iabscalc_type<3) then
           if(iproc==0) write(*,*)  ' Energy dependent potential has been asked in abscalc but  iabscalc_type<3 '
           if(nproc>1) call MPI_Finalize(ierr)
           stop '      Energy dependent potential has been asked in abscalc but  iabscalc_type<3    '
        endif
     endif


     call input_wf_diag(iproc,nproc,atoms_clone,&
          orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopotExtra,rhocore,pot_ion,&
          nlpspd,proj,pkernel,pkernel,ixc,psi,hpsi,psit,Gvirt,&
          nscatterarr,ngatherarr,nspin, in%potshortcut, -1, irrzon, phnons, GPU,in)
     
     if( iand( in%potshortcut,16)>0) then
        if(iproc==0) write(*,*) "re-reading electronic_density for Xanes energy dependent potential "
        call read_density_cube_old("electronic_density", n1i,n2i,n3i,1, hx ,hy ,hz, atoms%nat, rxyz_b2B, pot_bB )
        rhoXanes=0.0_gp
        do iz = 1,n3i
           do iy=1,n2i
              do ix=1,n1i
                 rhopottmp(ix,iy,iz +i3xcsh,1) =  pot_bB(ix  + (iy-1)*n1i  + (iz-1)*n1i*n2i,1)  
              enddo
           enddo
        enddo
        
        i_all=-product(shape(pot_bB))*kind(pot_bB)
        deallocate(pot_bB,stat=i_stat)
        call memocc(i_stat,i_all,'rho',subname)
        i_all=-product(shape(rxyz_b2B))*kind(rxyz_b2B)
        deallocate(rxyz_b2B,stat=i_stat)
        call memocc(i_stat,i_all,'rxyz',subname)


     endif

 

     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)
     
     i_all=-product(shape(irrzon))*kind(irrzon)
     deallocate(irrzon,stat=i_stat)
     call memocc(i_stat,i_all,'irrzon',subname)
     
     i_all=-product(shape(phnons))*kind(phnons)
     deallocate(phnons,stat=i_stat)
     call memocc(i_stat,i_all,'phnons',subname)
     
     i_all=-product(shape(atoms_clone%aocc))*kind(atoms_clone%aocc)
     deallocate(atoms_clone%aocc,stat=i_stat)
     call memocc(i_stat,i_all,'atoms_clone%aocc',subname)
     i_all=-product(shape(atoms_clone%iasctype))*kind(atoms_clone%iasctype)
     deallocate(atoms_clone%iasctype,stat=i_stat)
     call memocc(i_stat,i_all,'atoms_clone%iasctype',subname)

  endif


  if( iand( in%potshortcut,1)  .gt. 0 ) then

     inputpsi=in%inputPsiId

     nspin=in%nspin

     !Fake allocations
     allocate(irrzon(1,2,1+ndebug),stat=i_stat)
     call memocc(i_stat,irrzon,'irrzon',subname)
     allocate(phnons(2,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,phnons,'phnons',subname)

     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(iproc,nproc,atoms,&
          orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
          nlpspd,proj,pkernel,pkernel,ixc,psi,hpsi,psit,Gvirt,&
          nscatterarr,ngatherarr,nspin, in%potshortcut, -1, irrzon, phnons, GPU, in)

     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)

     i_all=-product(shape(irrzon))*kind(irrzon)
     deallocate(irrzon,stat=i_stat)
     call memocc(i_stat,i_all,'irrzon',subname)

     i_all=-product(shape(phnons))*kind(phnons)
     deallocate(phnons,stat=i_stat)
     call memocc(i_stat,i_all,'phnons',subname)

  end if

  nullify(psit)

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

     !put i3xcsh=0 for the moment, should be eliminated from the potential
     i3xcsh=0

!!$
!!$     rhopot(10,9,8+i3xcsh,1)=100.0

     if (in%output_grid == OUTPUT_GRID_DENSPOT) then
        if (in%output_grid_format == OUTPUT_GRID_FORMAT_TEXT) then
          if (iproc == 0) write(*,*) 'writing local_potential.pot'
           call plot_density_old(atoms%geocode,'local_potentialb2B.pot',iproc,nproc,&
                n1,n2,n3,n1i,n2i,n3i,n3p,&
                atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rhopot(1,1,1,1))
        else
           call plot_density_cube_old(atoms%geocode,'local_potentialb2B',iproc,nproc,&
                n1,n2,n3,n1i,n2i,n3i,n3p,&
                in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rhopot(1,1,1,1))
        endif
     end if

!!$     call  read_potfile4b2B("local_potential.pot",n1i_bB,n2i_bB,n3i_bB, pot_bB, alat1_bB, alat2_bB, alat3_bB)
!!$     print *, pot_bB(10  + (9-1)*n1i_bB  + (8-1)*n1i_bB*n2i_bB)
!!$     stop
     
     if(iproc==0) print *, " going to calculate spectra "
     
     if(  iand( in%potshortcut, 2)  > 0 ) then

        if(  iand( in%potshortcut, 16)  > 0 ) then
           b2BN=2
        else
           b2BN=1
        endif
        

        do b2Bcounter=1,b2BN

           if(b2Bcounter==1) then
              write(filename,'(A)' ) 'b2B_xanes'
              rhotarget=>rhopot
           else
              write(filename,'(A)') 'b2B_rho'
              rhotarget=>rhoXanes
           endif


        inquire(file=trim(trim(filename)//'.cube'),exist=exists)
        print *, "controllo ",  trim(filename)//'.cube', exists
        if(exists) then

           call read_cube(trim(filename),atoms%geocode,n1i_bB,n2i_bB,n3i_bB, &
                & nspin , hx_old ,hy_old ,hz_old ,pot_bB, nat_b2B, rxyz_b2B, iatype_b2B, znucl_b2B)
           !call read_density_cube_old(trim(filename), n1i_bB,n2i_bB,n3i_bB, 1 , hx_old ,hy_old ,hz_old , nat_b2B, rxyz_b2B, pot_bB )
           hx_old=hx_old*2
           hy_old=hy_old*2
           hz_old=hz_old*2


           if( (atoms%nat/nat_b2B)*nat_b2B /=  atoms%nat ) then
              if(iproc==0) write(*,*)  "   b2B_xanes cube  is not compatible with actual positions" 
              if(nproc>1) call MPI_Finalize(ierr)
              stop '      b2B_xanes cube  is not compatible with actual positions          '
           end if

        else

           if(b2BN>1) then
              if(iproc==0) write(*,*)  " b2B must be read only from *.cube when potential is energy dependent  " 
              if(nproc>1) call MPI_Finalize(ierr)
              stop '   b2B must be read only from *.cube when potential is energy dependent         '
           endif
           print  *, " reading atomic positions from file ","b2B_xanes.xyz"
           call read_atomic_file("b2B_xanes.xyz",iproc, atoms_b2B, rxyz_b2B )
           print *, "OK ", shape( rxyz_b2B )

           nat_b2B= (   Ubound(rxyz_b2B,2)  - Lbound(rxyz_b2B,2)  +  1 ) - ndebug


           if( (atoms%nat/nat_b2B)*nat_b2B /= atoms%nat ) then
              if(iproc==0) write(*,*)  "   b2B_xanes.xyz  is not compatible with actual positions" 
              if(nproc>1) call MPI_Finalize(ierr)
              stop '      b2B_xanes.xyz  is not compatible with actual positions          '
           end if
           
           print  *, " reading potential from file ","b2B_xanes.pot"
           !call  read_potfile4b2B("b2B_xanes.pot",n1i_bB,n2i_bB,n3i_bB, pot_bB, alat1_bB, alat2_bB, alat3_bB)
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


        allocate(rhopottmp( max(n1i_bB,n1i),max(n2i_bB,n2i),max(n3i_bB,n3d),in%nspin+ndebug),stat=i_stat)
        call memocc(i_stat,rhopottmp,'rhopottmp',subname)


        rhotarget=0.0_gp



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
        
        allocate(auxint(n1i+n2i+n3i+ndebug),stat=i_stat)
        call memocc(i_stat,auxint,'auxint',subname)

        Nreplicas = atoms%nat / nat_b2B
        dumvect3d(1)=atoms%alat1
        dumvect3d(2)=atoms%alat2
        dumvect3d(3)=atoms%alat3

        do ireplica=0, Nreplicas-1
           
           replicaoffset = (ireplica)*(   nat_b2B      )

           do j=1,3
              shift_b2B(j) = rxyz(j,1+replicaoffset) - rxyz_b2B(j,1) 

              do iat=2+ireplica*nat_b2B, (ireplica+1)*nat_b2B

                 shiftdiff = shift_b2B(j) - (rxyz(j,iat) -rxyz_b2B(j,iat-replicaoffset)   )
          
                 if( abs( shiftdiff )>1.0e-4 .and.  abs(abs(shiftdiff)-dumvect3d(j))  >1.0e-4) then
                    if(iproc==0) write(*,*)  "   b2B_xanes  positions are not compatible with actual positions" 
                    if(nproc>1) call MPI_Finalize(ierr)
                    stop '      b2B_xanes positions are not compatible with actual positions          '
                 end if
              enddo
           enddo

           print '(a,i6,a,3(1x,f18.14))',"for replica ", ireplica,  "SHIFT " , shift_b2B

           rhopottmp=0.0_gp
           do iz_bB = 1,n3i_bB
              do iy_bB=1,n2i_bB
                 do ix_bB=1,n1i_bB
                    rhopottmp(ix_bB,iy_bB,iz_bB +i3xcsh,1) =  pot_bB(ix_bB  + (iy_bB-1)*n1i_bB  + (iz_bB-1)*n1i_bB*n2i_bB,1)
                 enddo
              enddo
           enddo

           i_all=-product(shape(pot_bB))*kind(pot_bB)
           deallocate(pot_bB,stat=i_stat)
           call memocc(i_stat,i_all,'rho',subname)



           do iz_bB = 1,n3i_bB-1
              do iy_bB=1,n2i_bB
                 auxint = 0.0_gp
                 do ix_bB=1,n1i_bB
                    rx_bB = hx_old*(ix_bB-1)           /2.0   +  shift_b2B(1)
                    minX_B  =  max(1,NINT((rx_bB -8*hx_old/2)/(hx/2.0)))
                    maxX_B  =  min(n1i,NINT((rx_bB +8*hx_old/2)/(hx/2.0)))

                    minX_B  =  NINT((rx_bB -8*hx_old/2)/(hx/2.0))
                    maxX_B  =  NINT((rx_bB +8*hx_old/2)/(hx/2.0))

                    do ixnl= minX_B , maxX_B 
                       ix = mod(ixnl-1 + n1i , n1i) +1

                       rx = hx*(ix-1  )/2.0  

                       shiftdiff = (rx-rx_bB)
                       if ( abs(shiftdiff -atoms%alat1) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat1
                       if ( abs(shiftdiff +atoms%alat1) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat1

                       idelta = NINT( shiftdiff *2**15/(hx_old/2))  
                       factx = intfunc_y(nd/2+idelta)
!!$                    print *, rx, rx_bB, ix, ix_bB      , factx
                       auxint(ix) = auxint(ix) + &
                            factx * rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                    enddo
!!$                 print *, auxint(ix_bB) ,  rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                 enddo
                 rhopottmp(:,iy_bB,iz_bB+i3xcsh,1)=auxint(1:n1i)
              enddo
           enddo

           do iz_bB = 1,n3i_bB-1
              do ix_bB=1,n1i
                 auxint = 0.0_gp
                 do iy_bB=1,n2i_bB
                    ry_bB = hy_old*(iy_bB-1)/2.0   +  shift_b2B(2)
                    minY_B  =  max(1  ,NINT((ry_bB -8*hy_old/2)/(hy/2.0)))
                    maxY_B  =  min(n2i,NINT((ry_bB +8*hy_old/2)/(hy/2.0)))

                    minY_B  =  NINT((ry_bB -8*hy_old/2)/(hy/2.0))
                    maxY_B  =  NINT((ry_bB +8*hy_old/2)/(hy/2.0))

                    do iynl= minY_B , maxY_B 
                       iy = mod(iynl-1 + n2i , n2i) +1
                       
                       ry = hy*(iy-1  )/2.0  

                       shiftdiff = (ry-ry_bB)
                       if ( abs(shiftdiff -atoms%alat2) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat2
                       if ( abs(shiftdiff +atoms%alat2) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat2


                       idelta = NINT(shiftdiff *2**15/(hy_old/2))
                       facty = intfunc_y(nd/2+idelta)
                       auxint(iy) = auxint(iy) + &
                            facty * rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                    enddo
                 enddo
                 rhopottmp(ix_bB ,:,iz_bB+i3xcsh,1)=auxint(1:n2i)
              enddo
           enddo

           do ix_bB=1,n1i
              do iy_bB=1,n2i
                 auxint = 0.0_gp
                 do iz_bB = 1,n3i_bB-1
                    rz_bB = hz_old*(iz_bB-1)           /2.0   +  shift_b2B(3)

                    minZ_B  =  max(1  ,  NINT((rz_bB -8*hz_old/2)/(hz/2.0)))
                    maxZ_B  =  min(n3i-i3xcsh , NINT((rz_bB +8*hz_old/2)/(hz/2.0)))

                    minZ_B  =    NINT((rz_bB -8*hz_old/2)/(hz/2.0))
                    maxZ_B  =    NINT((rz_bB +8*hz_old/2)/(hz/2.0))

                    do iznl= minZ_B , maxZ_B 

                       iz = mod(iznl-1 + n3i , n3i) +1

                       rz = hz*(iz-1  )/2.0  

                       shiftdiff = (rz-rz_bB)
                       if ( abs(shiftdiff -atoms%alat3) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat3
                       if ( abs(shiftdiff +atoms%alat3) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat3

                       idelta = NINT( shiftdiff *2**15/(hz_old/2.0))     
                       factz = intfunc_y(nd/2+idelta)
                       auxint(iz+i3xcsh) = auxint(iz+i3xcsh) + &
                            factz * rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                    enddo
                 enddo
                 rhotarget(ix_bB ,iy_bB, : ,1)= rhotarget(ix_bB ,iy_bB, : ,1)+auxint(1:n3i)
              enddo
           enddo
        enddo
        i_all=-product(shape(rxyz_b2B))*kind(rxyz_b2B)
        deallocate(rxyz_b2B,stat=i_stat)
        call memocc(i_stat,i_all,'rxyz',subname)
        i_all=-product(shape(iatype_b2B))*kind(rxyz_b2B)
        deallocate(iatype_b2B,stat=i_stat)
        call memocc(i_stat,i_all,'iatype',subname)
        i_all=-product(shape(znucl_b2B))*kind(rxyz_b2B)
        deallocate(znucl_b2B,stat=i_stat)
        call memocc(i_stat,i_all,'znucl',subname)
        i_all=-product(shape(rhopottmp))*kind(rhopottmp)
        deallocate(rhopottmp,stat=i_stat)
        call memocc(i_stat,i_all,'rhopottmp',subname)
        i_all=-product(shape(auxint))*kind(auxint)
        deallocate(auxint,stat=i_stat)
        call memocc(i_stat,i_all,'auxint',subname)
        i_all=-product(shape(intfunc_y))*kind(intfunc_y)
        deallocate(intfunc_y,stat=i_stat)
        call memocc(i_stat,i_all,'intfunc_y',subname)
     enddo
        


        if (iproc == 0) write(*,*) 'writing NEW local_potential.pot'

        call plot_density_old(atoms%geocode,'local_potentialb2BNEW.pot',iproc,nproc,&
             n1,n2,n3,n1i,n2i,n3i,n3p,&
             atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rhopot(1,1,1+i3xcsh,1))


        print *," exiting b2B"



     endif

     if( iand( in%potshortcut,4)  .gt. 0 ) then
        do ix=1,n1i
           do iy=1,n2i
              do iz = 1,n3i
                 rhopot(ix ,iy, iz ,1)= rhopot(ix ,iy, iz ,1)+rhopotExtra(ix ,iy, iz ,1)
              enddo
           enddo
        enddo
        i_all=-product(shape(rhopotExtra))*kind(rhopotExtra)
        deallocate(rhopotExtra,stat=i_stat)
        call memocc(i_stat,i_all,'rhopotExtra',subname)
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
        potmodified_maxr=0
        potmodified_shift=0

        

        do ix=1,n1i
           do iy=1,n2i
              do iz = 1,n3p
                 rx = hx*(ix-1)           /2.0  -  rxyz(1,in%iat_absorber )
                 ry = hy*(iy-1)           /2.0  -  rxyz(2,in%iat_absorber )
                 rz = hz*(iz-1 +i3xcsh + i3s -1 )/2.0  -  rxyz(3,in%iat_absorber )

                 r  = sqrt( rx*rx+ry*ry+rz*rz)

                 if(r>2.0) then
                    
                    if( r>29) then
                       rhopot(ix,iy,iz+i3xcsh,1)=0.0
                    else
                       igrid = binary_search( r, radpot, radpotcount )
                       rhopot(ix,iy,iz+i3xcsh,1) = &
                            ( radpot(igrid,2)*(radpot(igrid+1,1)-R) + radpot(igrid+1,2)*(R-radpot(igrid,1)) )/&
                            ( radpot(igrid+1,1) -radpot(igrid,1) )
                    endif
                 else
                    if(potmodified_maxr<r) then 
                       potmodified_maxr=r
                       igrid = binary_search( r, radpot, radpotcount )
                       potmodified_shift =&
                            ( radpot(igrid,2)*(radpot(igrid+1,1)-R) + radpot(igrid+1,2)*(R-radpot(igrid,1)) )/&
                            ( radpot(igrid+1,1) -radpot(igrid,1) ) &
                            -rhopot(ix,iy,iz+i3xcsh,1) 
                    endif
                 endif
                 
                 if(r<minr) minr=r
                 
                 if( r.ge.2.0) then
                    harmo = (rx+ry+rz)/sqrt(3.0)/r *sqrt( 3.0/4.0/3.1415926535)
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
        do ix=1,n1i
           do iy=1,n2i
              do iz = 1,n3p
                 rx = hx*(ix-1)           /2.0  -  rxyz(1,in%iat_absorber )
                 ry = hy*(iy-1)           /2.0  -  rxyz(2,in%iat_absorber )
                 rz = hz*(iz-1 +i3xcsh + i3s -1 )/2.0  -  rxyz(3,in%iat_absorber )

                 r  = sqrt( rx*rx+ry*ry+rz*rz)

                 if(r<=2.0) then
                    print *, " potmodified_shift ", potmodified_shift
                    rhopot(ix,iy,iz+i3xcsh,1)=rhopot(ix,iy,iz+i3xcsh,1)+potmodified_shift
                 endif
              enddo
           enddo
        enddo
     end if
     infocode=0

     if (in%iabscalc_type==2) then
        call xabs_lanczos(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
             rhopot(1,1,1+i3xcsh,1) ,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU &
             , in%iat_absorber  , in )
        
     else if (in%iabscalc_type==1) then
        call xabs_chebychev(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
             rhopot(1,1,1+i3xcsh,1) ,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU &
             , in%iat_absorber, in)
     else if (in%iabscalc_type==3) then
        call xabs_cg(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
             rhopot(1,1,1+i3xcsh,1) ,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU &
             , in%iat_absorber, in, rhoXanes(1,1,1,1))
     else
        if (iproc == 0) write(*,*)' iabscalc_type not known, does not perform calculation'
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

       if(associated(rhoXanes)) then
          i_all=-product(shape(rhoXanes))*kind(rhoXanes)
          deallocate(rhoXanes,stat=i_stat)
          call memocc(i_stat,i_all,'rhoXanes',subname)
       endif




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
    

    call deallocate_bounds(atoms%geocode,Glr%hybrid_on,  Glr%bounds,subname)

!!$    if (atoms%geocode == 'F') then
!!$       call deallocate_bounds(Glr%bounds,subname)
!!$    end if
!!$    
!!$    if (atoms%geocode == 'P' .and. Glr%hybrid_on) then 
!!$       
!!$       i_all=-product(shape(Glr%bounds%kb%ibxy_f))*kind(Glr%bounds%kb%ibxy_f)
!!$       deallocate(Glr%bounds%kb%ibxy_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'Glr%bounds%kb%ibxy_f',subname)
!!$       
!!$       i_all=-product(shape(Glr%bounds%kb%ibxz_f))*kind(Glr%bounds%kb%ibxz_f)
!!$       deallocate(Glr%bounds%kb%ibxz_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'Glr%bounds%kb%ibxz_f',subname)
!!$       
!!$       i_all=-product(shape(Glr%bounds%kb%ibyz_f))*kind(Glr%bounds%kb%ibyz_f)
!!$       deallocate(Glr%bounds%kb%ibyz_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'Glr%bounds%kb%ibyz_f',subname)
!!$       
!!$       i_all=-product(shape(Glr%bounds%sb%ibxy_ff))*kind(Glr%bounds%sb%ibxy_ff)
!!$       deallocate(Glr%bounds%sb%ibxy_ff,stat=i_stat)
!!$       call memocc(i_stat,i_all,'ibxy_ff',subname)
!!$       i_all=-product(shape(Glr%bounds%sb%ibzzx_f))*kind(Glr%bounds%sb%ibzzx_f)
!!$       deallocate(Glr%bounds%sb%ibzzx_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'ibzzx_f',subname)
!!$       i_all=-product(shape(Glr%bounds%sb%ibyyzz_f))*kind(Glr%bounds%sb%ibyyzz_f)
!!$       deallocate(Glr%bounds%sb%ibyyzz_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'ibyyzz_f',subname)
!!$       
!!$       i_all=-product(shape(Glr%bounds%gb%ibyz_ff))*kind(Glr%bounds%gb%ibyz_ff)
!!$       deallocate(Glr%bounds%gb%ibyz_ff,stat=i_stat)
!!$       call memocc(i_stat,i_all,'ibyz_ff',subname)
!!$       
!!$       i_all=-product(shape(Glr%bounds%gb%ibzxx_f))*kind(Glr%bounds%gb%ibzxx_f)
!!$       deallocate(Glr%bounds%gb%ibzxx_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'ibzxx_f',subname)
!!$       
!!$       i_all=-product(shape(Glr%bounds%gb%ibxxyy_f))*kind(Glr%bounds%gb%ibxxyy_f)
!!$       deallocate(Glr%bounds%gb%ibxxyy_f,stat=i_stat)
!!$       call memocc(i_stat,i_all,'ibxxyy_f',subname)
!!$    endif

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

  END SUBROUTINE deallocate_before_exiting

END SUBROUTINE abscalc
