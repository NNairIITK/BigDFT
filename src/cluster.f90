!>   Routines to use bigdft as a blackbox
!!
!! @author
!!   Copyright (C) 2005-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!!
!!
subroutine call_bigdft(nproc,iproc,atoms,rxyz0,in,energy,fxyz,fnoise,rst,infocode)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(input_variables),intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  type(restart_objects), intent(inout) :: rst
  integer, intent(inout) :: infocode
  real(gp), intent(out) :: energy,fnoise
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz0
  real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
  !local variables
  character(len=*), parameter :: subname='call_bigdft'
  character(len=40) :: comment
  integer :: i_stat,i_all,ierr,inputPsiId_orig,iat

  !temporary interface
  interface
     subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,fnoise,&
          psi,Glr,gaucoeffs,gbd,orbs,rxyz_old,hx_old,hy_old,hz_old,in,GPU,infocode)
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
       type(GPU_pointers), intent(inout) :: GPU
       real(gp), intent(out) :: energy,fnoise
       real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz_old
       real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
       real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
       real(wp), dimension(:), pointer :: psi
       real(wp), dimension(:,:), pointer :: gaucoeffs
     END SUBROUTINE cluster
  end interface

  !put a barrier for all the processes
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !fill the rxyz array with the positions
  !wrap the atoms in the periodic directions when needed
  do iat=1,atoms%nat
     if (atoms%geocode == 'P') then
        rst%rxyz_new(1,iat)=modulo(rxyz0(1,iat),atoms%alat1)
        rst%rxyz_new(2,iat)=modulo(rxyz0(2,iat),atoms%alat2)
        rst%rxyz_new(3,iat)=modulo(rxyz0(3,iat),atoms%alat3)
     else if (atoms%geocode == 'S') then
        rst%rxyz_new(1,iat)=modulo(rxyz0(1,iat),atoms%alat1)
        rst%rxyz_new(2,iat)=rxyz0(2,iat)
        rst%rxyz_new(3,iat)=modulo(rxyz0(3,iat),atoms%alat3)
     else if (atoms%geocode == 'F') then
        rst%rxyz_new(1,iat)=rxyz0(1,iat)
        rst%rxyz_new(2,iat)=rxyz0(2,iat)
        rst%rxyz_new(3,iat)=rxyz0(3,iat)
     end if
  end do

  !assign the verbosity of the output
  !the verbose variables is defined in module_base
  verbose=in%verbosity

  inputPsiId_orig=in%inputPsiId

  loop_cluster: do

     if (in%inputPsiId == 0 .and. associated(rst%psi)) then
        i_all=-product(shape(rst%psi))*kind(rst%psi)
        deallocate(rst%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(rst%orbs%eval))*kind(rst%orbs%eval)
        deallocate(rst%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(rst%Glr%wfd,subname)
     end if

     call cluster(nproc,iproc,atoms,rst%rxyz_new,energy,fxyz,fnoise,&
          rst%psi,rst%Glr,rst%gaucoeffs,rst%gbd,rst%orbs,&
          rst%rxyz_old,rst%hx_old,rst%hy_old,rst%hz_old,in,rst%GPU,infocode)

     if (in%inputPsiId==1 .and. infocode==2) then
        if (in%gaussian_help) then
           in%inputPsiId=11
        else
           in%inputPsiId=0
        end if
     else if ((in%inputPsiId==1 .or. in%inputPsiId==0) .and. infocode==1) then
        !in%inputPsiId=0 !better to diagonalise than to restart an input guess
        in%inputPsiId=1
        if(iproc==0) then
           write(*,*)&
                ' WARNING: Self-consistent cycle did not met convergence criteria'
        end if
        exit loop_cluster
     else if (in%inputPsiId == 0 .and. infocode==3) then
        if (iproc == 0) then
           write( *,'(1x,a)')'Convergence error, cannot proceed.'
           write( *,'(1x,a)')' writing positions in file posfail.xyz then exiting'
           write(comment,'(a)')'UNCONVERGED WF '
           !call wtxyz('posfail',energy,rxyz,atoms,trim(comment))

           call write_atomic_file("posfail",energy,rst%rxyz_new,atoms,trim(comment))

        end if

        i_all=-product(shape(rst%psi))*kind(rst%psi)
        deallocate(rst%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(rst%orbs%eval))*kind(rst%orbs%eval)
        deallocate(rst%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(rst%Glr%wfd,subname)

        !finalize memory counting (there are still at least positions and the forces allocated)
        call memocc(0,0,'count','stop')

        if (nproc > 1) call MPI_FINALIZE(ierr)

        stop 'unnormal end'
     else
        exit loop_cluster
     end if

  end do loop_cluster

  !preserve the previous value
  in%inputPsiId=inputPsiId_orig

  !put a barrier for all the processes
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE call_bigdft



!>  Main routine which does self-consistent loop.
!!  Does not parse input file and no geometry optimization.
!!
!!   inputPsiId = 0 : compute input guess for Psi by subspace diagonalization of atomic orbitals
!!   inputPsiId = 1 : read waves from argument psi, using n1, n2, n3, hgrid and rxyz_old
!!                    as definition of the previous system.
!!   inputPsiId = 2 : read waves from disk
!!   does an electronic structure calculation. Output is the total energy and the forces 
!!   psi, keyg, keyv and eval should be freed after use outside of the routine.
!!   infocode -> encloses some information about the status of the run
!!            =0 run succesfully succeded
!!            =1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
!!               forces may be meaningless   
!!            =2 (present only for inputPsiId=1) gnrm of the first iteration > 1 AND growing in
!!               the second iteration OR grnm 1st >2.
!!               Input wavefunctions need to be recalculated. Routine exits.
!!            =3 (present only for inputPsiId=0) gnrm > 4. SCF error. Routine exits.
!!
!!
subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,fnoise,&
     psi,Glr,gaucoeffs,gbd,orbs,rxyz_old,hx_old,hy_old,hz_old,in,GPU,infocode)
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  use libxc_functionals
  use vdwcorrection, only: vdwcorrection_calculate_energy, vdwcorrection_calculate_forces, vdwcorrection_warnings
  use esatto
  use ab6_symmetry
  use m_ab6_mixing
  implicit none
  integer, intent(in) :: nproc,iproc
  real(gp), intent(inout) :: hx_old,hy_old,hz_old
  type(input_variables), intent(in) :: in
  type(locreg_descriptors), intent(inout) :: Glr
  type(atoms_data), intent(inout) :: atoms
  type(gaussian_basis), intent(inout) :: gbd
  type(orbitals_data), intent(inout) :: orbs
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz_old
  real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
  integer, intent(out) :: infocode
  real(gp), intent(out) :: energy,fnoise
  real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
  real(wp), dimension(:), pointer :: psi
  real(wp), dimension(:,:), pointer :: gaucoeffs
  !local variables
  character(len=*), parameter :: subname='cluster'
  character(len=3) :: PSquiet
  character(len=4) :: f4
  character(len=5) :: gridformat, wfformat, final_out
  character(len=50) :: filename
  character(len=500) :: errmess
  logical :: endloop,endlooprp,allfiles,onefile,refill_proj
  logical :: DoDavidson,counterions,DoLastRunThings=.false.,lcs,scpot
  integer :: ixc,ncong,idsx,ncongt,nspin,nsym,icycle,potden
  integer :: nvirt,ndiis_sd_sw,norbv,idsx_actual_before
  integer :: nelec,ndegree_ip,j,i,iorb,npoints
  integer :: n1_old,n2_old,n3_old,n3d,n3p,n3pi,i3xcsh,i3s,n1,n2,n3
  integer :: ncount0,ncount1,ncount_rate,ncount_max,n1i,n2i,n3i
  integer :: iat,i_all,i_stat,iter,itrp,ierr,jproc,inputpsi,igroup,ikpt,jkpt,ispin
  real :: tcpu0,tcpu1
  real(kind=8) :: crmult,frmult,cpmult,fpmult,gnrm_cv,rbuf,hxh,hyh,hzh,hx,hy,hz
  real(gp) :: peakmem,evsum
  real(gp) :: eion,epot_sum,ekin_sum,eproj_sum,eexctX,ehart,eexcu,vexcu,rpnrm,gnrm,gnrm_zero
  real(gp) :: trH,energybs,tt,tel,ehart_fake,psoffset
  real(kind=8) :: ttsum
  real(gp) :: edisp ! Dispersion energy
  type(wavefunctions_descriptors) :: wfd_old
  type(nonlocal_psp_descriptors) :: nlpspd
  type(communications_arrays) :: comms, commsv
  type(orbitals_data) :: orbsv
  type(gaussian_basis) :: Gvirt
  type(diis_objects) :: diis
  real(gp), dimension(3) :: shift,chargec
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(kind=8), dimension(:), allocatable :: rho,psirocc,psirvirt
  real(gp), dimension(:,:), allocatable :: radii_cf,gxyz,fion,thetaphi,band_structure_eval
  real(gp), dimension(:,:),allocatable :: fdisp
  ! Charge density/potential,ionic potential, pkernel
  type(ab6_mixing_object) :: mix
  real(dp), dimension(:), allocatable :: pot_ion,rhopot,counter_ions
  real(kind=8), dimension(:,:,:,:), allocatable :: pot,potxc,dvxcdrho
  real(wp), dimension(:), pointer :: potential
  real(wp), dimension(:,:), pointer :: pot_from_disk
  real(dp), dimension(:), pointer :: pkernel,pkernelseq
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(kind=8), dimension(:), pointer :: hpsi,psit,psi_old,psivirt
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: proj,gbd_occ,rhocore
  ! Arrays for the symmetrisation.
  integer, dimension(:,:,:), allocatable :: irrzon
  real(dp), dimension(:,:,:), allocatable :: phnons
  ! Variables for the virtual orbitals and band diagram.
  integer :: nkptv, nvirtu, nvirtd
  real(gp), allocatable :: wkptv(:)


  ! ----------------------------------

  !copying the input variables for readability
  !this section is of course not needed
  !note that this procedure is convenient ONLY in the case of scalar variables
  !an array would have been copied, thus occupying more memory space
  !Hence WARNING: these variables are copied, in case of an update the new value should be 
  !reassigned inside the structure


  crmult=in%crmult
  frmult=in%frmult
  cpmult=in%frmult
  fpmult=in%frmult
  ixc=in%ixc
  gnrm_cv=in%gnrm_cv
  ncong=in%ncong
  idsx=in%idsx
  rbuf=in%rbuf
  ncongt=in%ncongt
  nspin=in%nspin
  write(gridformat, "(A)") ""
  select case (in%output_grid_format)
     case (OUTPUT_GRID_FORMAT_ETSF)
        write(gridformat, "(A)") ".etsf"
     case (OUTPUT_GRID_FORMAT_CUBE)
        write(gridformat, "(A)") ".bin"
  end select
  write(wfformat, "(A)") ""
  select case (in%output_wf_format)
     case (WF_FORMAT_ETSF)
        write(wfformat, "(A)") ".etsf"
     case (WF_FORMAT_BINARY)
        write(wfformat, "(A)") ".bin"
  end select
     
  norbv=abs(in%norbv)
  nvirt=in%nvirt

  hx=in%hx
  hy=in%hy
  hz=in%hz

  call libxc_functionals_init(ixc, nspin)

  !character string for quieting the Poisson solver
  if (verbose >1) then
     PSquiet='NO'
  else
     PSquiet='YES'
  end if

  if (iproc == 0) then
     write( *,'(1x,a,1x,i0)') &
          '===================== BigDFT Wavefunction Optimization =============== inputPsiId=',&
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

  ! We save the variables that defined the previous psi if the restart is active
  if (in%inputPsiId == 1) then
     !regenerate grid spacings
     if (atoms%geocode == 'P') then
        call correct_grid(atoms%alat1,hx_old,Glr%d%n1)
        call correct_grid(atoms%alat2,hy_old,Glr%d%n2)
        call correct_grid(atoms%alat3,hz_old,Glr%d%n3)
     else if (atoms%geocode == 'S') then 
        call correct_grid(atoms%alat1,hx_old,Glr%d%n1)
        call correct_grid(atoms%alat3,hz_old,Glr%d%n3)
     end if
     call copy_old_wavefunctions(nproc,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
          Glr%wfd,psi,n1_old,n2_old,n3_old,wfd_old,psi_old)
  else if (in%inputPsiId == 11) then
     !deallocate wavefunction and descriptors for placing the gaussians

     call deallocate_wfd(Glr%wfd,subname)

     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)

  end if

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

  ! A message about dispersion forces.
  if (iproc == 0) call vdwcorrection_warnings(atoms, in)

  !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
  ndegree_ip=16 !default value 
  call createKernel(iproc,nproc,atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,pkernel,&
       quiet=PSquiet)

  !create the sequential kernel if the exctX parallelisation scheme requires it
  if (libxc_functionals_exctXfac() /= 0.0_gp .and. in%exctxpar=='OP2P' .and. nproc > 1) then
     call createKernel(0,1,atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,&
          pkernelseq,quiet='YES')
  else
     pkernelseq => pkernel
  end if


  ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
  call timing(iproc,'CrtDescriptors','ON')
  call createWavefunctionsDescriptors(iproc,hx,hy,hz,&
       atoms,rxyz,radii_cf,crmult,frmult,Glr)
  call timing(iproc,'CrtDescriptors','OF')

  !allocate communications arrays (allocate it before Projectors because of the definition
  !of iskpts and nkptsp)
  call orbitals_communicators(iproc,nproc,Glr,orbs,comms)  

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
          orbs%norb,orbs%nkpts,nlpspd%nprojel,atoms%atomnames,0,in%nspin,peakmem)
  end if

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

  !calculate the irreductible zone, if necessary.
  if (atoms%symObj >= 0) then
     call ab6_symmetry_get_n_sym(atoms%symObj, nsym, i_stat)
     if (nsym > 1) then
        ! Current third dimension is set to 1 always
        ! since nspin == nsppol always in BigDFT
        allocate(irrzon(n1i*n2i*n3i,2,1+ndebug),stat=i_stat)
        call memocc(i_stat,irrzon,'irrzon',subname)
        allocate(phnons(2,n1i*n2i*n3i,1+ndebug),stat=i_stat)
        call memocc(i_stat,phnons,'phnons',subname)
        call ab6_symmetry_get_irreductible_zone(atoms%symObj, irrzon, phnons, &
             & n1i, n2i, n3i, in%nspin, in%nspin, i_stat)
     end if
  end if
  if (.not. allocated(irrzon)) then
     ! Allocate anyway to small size other size the bounds check does not pass.
     allocate(irrzon(1,2,1+ndebug),stat=i_stat)
     call memocc(i_stat,irrzon,'irrzon',subname)
     allocate(phnons(2,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,phnons,'phnons',subname)
  end if

  !allocate ionic potential
  if (n3pi > 0) then
     allocate(pot_ion(n1i*n2i*n3pi+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  else
     allocate(pot_ion(1+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  end if

  !here calculate the ionic energy and forces accordingly
  allocate(fion(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fion,'fion',subname)

  call IonicEnergyandForces(iproc,nproc,atoms,hxh,hyh,hzh,in%elecfield,rxyz,eion,fion,&
       psoffset,in%nvacancy,n1,n2,n3,n1i,n2i,n3i,i3s+i3xcsh,n3pi,pot_ion,pkernel)

  call createIonicPotential(atoms%geocode,iproc,nproc,atoms,rxyz,hxh,hyh,hzh,&
       in%elecfield,n1,n2,n3,n3pi,i3s+i3xcsh,n1i,n2i,n3i,pkernel,pot_ion,psoffset,&
       in%nvacancy,in%correct_offset)

  !inquire for the counter_ion potential calculation (for the moment only xyz format)
  inquire(file='posinp_ci.xyz',exist=counterions)
  if (counterions) then
     if (n3pi > 0) then
        allocate(counter_ions(n1i*n2i*n3pi+ndebug),stat=i_stat)
        call memocc(i_stat,counter_ions,'counter_ions',subname)
     else
        allocate(counter_ions(1+ndebug),stat=i_stat)
        call memocc(i_stat,counter_ions,'counter_ions',subname)
     end if

     call CounterIonPotential(atoms%geocode,iproc,nproc,in,shift,&
          hxh,hyh,hzh,Glr%d,n3pi,i3s,pkernel,counter_ions)

     !sum that to the ionic potential
     call axpy(Glr%d%n1i*Glr%d%n2i*n3p,1.0_dp,counter_ions(1),1,&
          pot_ion(1),1)

  end if

  !this can be inserted inside the IonicEnergyandForces routine
  !(after insertion of the non-regression test)
  call vdwcorrection_calculate_energy(edisp,rxyz,atoms,in,iproc)

  allocate(fdisp(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fdisp,'fdisp',subname)
  !this can be inserted inside the IonicEnergyandForces routine
  call vdwcorrection_calculate_forces(fdisp,rxyz,atoms,in) 

  !Allocate Charge density, Potential in real space
  if (n3d >0) then
     allocate(rhopot(n1i*n2i*n3d*in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  else
     allocate(rhopot(in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  end if
  !Allocate XC potential
  if (n3p >0) then
     allocate(potxc(n1i,n2i,n3p,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,potxc,'potxc',subname)
  else
     allocate(potxc(1,1,1,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,potxc,'potxc',subname)
  end if

  !check if non-linear core correction should be applied, and allocate the 
  !pointer if it is the case
  call calculate_rhocore(iproc,atoms,Glr%d,rxyz,hxh,hyh,hzh,i3s,i3xcsh,n3d,n3p,rhocore)

  !check the communication distribution
  call check_communications(iproc,nproc,orbs,Glr,comms)

  !avoid allocation of the eigenvalues array in case of restart
  if (in%inputPsiId /= 1 .and. in%inputPsiId /= 11) then
     allocate(orbs%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
     call memocc(i_stat,orbs%eval,'eval',subname)
  end if

  inputpsi=in%inputPsiId

  !for the inputPsiId==2 case, check 
  !if the wavefunctions are all present
  !otherwise switch to normal input guess
  if (in%inputPsiId ==2) then
     allfiles=.true.
     ! Test ETSF file.
     inquire(file="wavefunction.etsf",exist=onefile)
     if (.not. onefile) then
        ! Test bin files
        do iorb=1,orbs%norb*orbs%nspinor
           write(f4,'(i4.4)')  iorb
           filename = 'wavefunction.bin.'//f4
           inquire(file=filename,exist=onefile)
           allfiles=allfiles .and. onefile
           if (.not. allfiles) then
              exit
           end if
        end do
        if (.not. allfiles) then           
           allfiles = .true.
           ! Test plain files.
           do iorb=1,orbs%norb*orbs%nspinor
              write(f4,'(i4.4)')  iorb
              filename = 'wavefunction.'//f4
              inquire(file=filename,exist=onefile)
              allfiles=allfiles .and. onefile
              if (.not. allfiles) then
                 exit
              end if
           end do
        end if
     end if
     if (.not. allfiles) then
        if (iproc == 0) write(*,*)' WARNING: Missing wavefunction files, switch to normal input guess'
        inputpsi = 0
     end if
  end if

  !all the input formats need to allocate psi except the LCAO input_guess
  if (inputpsi /= 0) then
     allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)
  end if

  ! INPUT WAVEFUNCTIONS, added also random input guess
  select case(inputpsi)
  case(INPUT_PSI_EMPTY)
     !allocate fake psit and hpsi
     allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hpsi,'hpsi',subname)
     if (nproc > 1) then
        allocate(psit(orbs%npsidim+ndebug),stat=i_stat)
        call memocc(i_stat,psit,'psit',subname)
     else
        psit => psi
     end if
     !fill the rhopot array with the read potential if needed
     if (trim(in%band_structure_filename) /= '') then
        !only the first processor should read this
        if (iproc == 0) then
           write(*,'(1x,a)')'Reading local potential from file:'//trim(in%band_structure_filename)
           call read_density(trim(in%band_structure_filename),atoms%geocode,&
                n1i,n2i,n3i,nspin,hxh,hyh,hzh,pot_from_disk)
           if (nspin /= in%nspin) stop
        else
           allocate(pot_from_disk(1,in%nspin+ndebug),stat=i_stat)
           call memocc(i_stat,pot_from_disk,'pot_from_disk',subname)
        end if

        if (nproc > 1) then
           do ispin=1,in%nspin
              call MPI_SCATTERV(pot_from_disk(1,ispin),&
                   ngatherarr(0,1),ngatherarr(0,2),mpidtypw, &
                   rhopot((ispin-1)*Glr%d%n1i*Glr%d%n2i*n3p+1),&
                   Glr%d%n1i*Glr%d%n2i*n3p,mpidtypw,0,MPI_COMM_WORLD,ierr)
           end do
        else
           call dcopy(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*in%nspin,pot_from_disk,1,rhopot,1)
        end if

        i_all=-product(shape(pot_from_disk))*kind(pot_from_disk)
        deallocate(pot_from_disk,stat=i_stat)
        call memocc(i_stat,i_all,'pot_from_disk',subname)

        !add pot_ion potential to the local_potential
        !do ispin=1,in%nspin
        !   !spin up and down together with the XC part
        !   call axpy(Glr%d%n1i*Glr%d%n2i*n3p,1.0_dp,pot_ion(1),1,&
        !        rhopot((ispin-1)*Glr%d%n1i*Glr%d%n2i*n3p+1),1)
        !end do
     end if

  case(INPUT_PSI_RANDOM)

     if (iproc == 0) then
        write( *,'(1x,a)')&
             '------------------------------------------------ Random wavefunctions initialization'
     end if

     !random initialisation of the wavefunctions

     psi=0.0d0
     ttsum=0.0d0
     do i=1,orbs%npsidim
        do j=0,iproc-1
           call random_number(tt)
        end do
        call random_number(tt)
        psi(i)=real(tt,wp)*0.01_wp
        ttsum=ttsum+psi(i)
        do j=iproc+1,nproc
           call random_number(tt)
        end do
     end do

     orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0

  case(INPUT_PSI_CP2K)

     !import gaussians form CP2K (data in files gaubasis.dat and gaucoeff.dat)
     !and calculate eigenvalues
     if (iproc == 0) then
        write(*,'(1x,a)')&
             '--------------------------------------------------------- Import Gaussians from CP2K'
     end if

     if (in%nspin /= 1) then
        if (iproc==0) then
           write(*,'(1x,a)')&
                'Gaussian importing is possible only for non-spin polarised calculations'
           write(*,'(1x,a)')&
                'The reading rules of CP2K files for spin-polarised orbitals are not implemented'
        end if
        stop
     end if

     !eliminate the old_import_gaussians routine
     call parse_cp2k_files(iproc,'gaubasis.dat','gaucoeff.dat',&
          atoms%nat,atoms%ntypes,orbs,atoms%iatype,rxyz,gbd,gaucoeffs)
     call gaussians_to_wavelets_new(iproc,nproc,Glr,orbs,hx,hy,hz,gbd,gaucoeffs,psi)
     !deallocate gaussian structure and coefficients
     call deallocate_gwf(gbd,subname)
     i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
     deallocate(gaucoeffs,stat=i_stat)
     call memocc(i_stat,i_all,'gaucoeffs',subname)
     nullify(gbd%rxyz)

     !call dual_gaussian_coefficients(orbs%norbp,gbd,gaucoeffs)
     orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0

  case(INPUT_PSI_LCAO)
     nspin=in%nspin
     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(iproc,nproc, atoms,&
          orbs,norbv,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
          nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,Gvirt,&
          nscatterarr,ngatherarr,nspin,0,atoms%symObj,irrzon,phnons,GPU,in)
     if (nvirt > norbv) then
        nvirt = norbv
     end if
  case(INPUT_PSI_MEMORY_WVL)
     !these parts should be reworked for the non-collinear spin case

     !restart from previously calculated wavefunctions, in memory
     if (iproc == 0) then
        write( *,'(1x,a)')&
             '-------------------------------------------------------------- Wavefunctions Restart'
     end if

     call reformatmywaves(iproc,orbs,atoms,hx_old,hy_old,hz_old,&
          n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,hx,hy,hz,n1,n2,n3,rxyz,Glr%wfd,psi)

     call deallocate_wfd(wfd_old,subname)

     i_all=-product(shape(psi_old))*kind(psi_old)
     deallocate(psi_old,stat=i_stat)
     call memocc(i_stat,i_all,'psi_old',subname)

  case(INPUT_PSI_DISK_WVL)
     !restart from previously calculated wavefunctions, on disk
     if (iproc == 0) then
        write( *,'(1x,a)')&
             '---------------------------------------------------- Reading Wavefunctions from disk'
     end if

     call readmywaves(iproc,"wavefunction", &
          & orbs,n1,n2,n3,hx,hy,hz,atoms,rxyz_old,rxyz,Glr%wfd,psi)

     if (in%itrpmax > 1) then
        !recalculate orbitals occupation numbers
        call evaltoocc(iproc,nproc,.false.,in%Tel,orbs)
     end if

  case(INPUT_PSI_MEMORY_GAUSS)
     !restart from previously calculated gaussian coefficients
     if (iproc == 0) then
        write( *,'(1x,a)')&
             '--------------------------------------- Quick Wavefunctions Restart (Gaussian basis)'
     end if

     call restart_from_gaussians(iproc,nproc,orbs,Glr,hx,hy,hz,psi,gbd,gaucoeffs)

  case(INPUT_PSI_DISK_GAUSS)
     !reading wavefunctions from gaussian file
     if (iproc == 0) then
        write( *,'(1x,a)')&
             '------------------------------------------- Reading Wavefunctions from gaussian file'
     end if

     call read_gaussian_information(orbs,gbd,gaucoeffs,'wavefunctions.gau')
     !associate the new positions, provided that the atom number is good
     if (gbd%nat == atoms%nat) then
        gbd%rxyz=>rxyz
     else
        !        if (iproc == 0) then
        write( *,*)&
             ' ERROR: the atom number does not coincide with the number of gaussian centers'
        !        end if
        stop
     end if

     call restart_from_gaussians(iproc,nproc,orbs,Glr,hx,hy,hz,psi,gbd,gaucoeffs)

  case default

     !     if (iproc == 0) then
     write( *,'(1x,a,I0,a)')'ERROR: illegal value of inputPsiId (', in%inputPsiId, ').'
     call input_psi_help()
     stop
     !     end if

  end select

  !all the input format need first_orthon except the LCAO input_guess
  if (inputpsi /= 0 .and. inputpsi /=-1000) then
     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,orbs,Glr%wfd,comms,psi,hpsi,psit,in)
  end if

  !save the new atomic positions in the rxyz_old array
  do iat=1,atoms%nat
     rxyz_old(1,iat)=rxyz(1,iat)
     rxyz_old(2,iat)=rxyz(2,iat)
     rxyz_old(3,iat)=rxyz(3,iat)
  enddo
  !save the new grid spacing into the hgrid_old value
  hx_old=in%hx
  hy_old=in%hy
  hz_old=in%hz

  ! allocate arrays necessary for DIIS convergence acceleration
  !the allocation with npsidim is not necessary here since DIIS arrays
  !are always calculated in the transpsed form
  if (idsx > 0) then
     call allocate_diis_objects(idsx,sum(comms%ncntt(0:nproc-1)),orbs%nkptsp,diis,subname)
  endif

  !allocate arrays for the GPU if a card is present
  if (GPUconv) then
     call prepare_gpu_for_locham(Glr%d%n1,Glr%d%n2,Glr%d%n3,in%nspin,&
          hx,hy,hz,Glr%wfd,orbs,GPU)
  end if
  !the same with OpenCL, but they cannot exist at same time
  if (OCLconv) then
     call allocate_data_OCL(Glr%d%n1,Glr%d%n2,Glr%d%n3,atoms%geocode,&
          in%nspin,hx,hy,hz,Glr%wfd,orbs,GPU)
     if (iproc == 0) write(*,*)&
          'GPU data allocated'
  end if

  energy=1.d10
  energybs=1.d10
  gnrm=1.d10
  rpnrm=1.d10
  gnrm_zero=0.0d0
  ekin_sum=0.d0 
  epot_sum=0.d0 
  eproj_sum=0.d0
  !diis initialisation variables
  diis%alpha=in%alphadiis
  diis%alpha_max=in%alphadiis
  diis%energy=1.d10
  !minimum value of the energy during the minimisation procedure
  diis%energy_min=1.d10
  !local variable for the diis history
  diis%idsx=idsx
  !number of switching betweed DIIS and SD during self-consistent loop
  ndiis_sd_sw=0
  !previous value of idsx_actual to control if switching has appeared
  idsx_actual_before=diis%idsx
  !logical control variable for switch DIIS-SD
  diis%switchSD=.false.

  !end of the initialization part
  call timing(iproc,'INIT','PR')

  !Davidson is set to false first because used in deallocate_before_exiting
  DoDavidson= .false.

  !allocate the rhopot_old array needed for mixing
  if (in%iscf < 10) then
     potden = AB6_MIXING_POTENTIAL
     npoints = n1i*n2i*n3p
  else
     potden = AB6_MIXING_DENSITY
     npoints = n1i*n2i*n3d
  end if
  if (n3d >0 .and. in%itrpmax>1) then
     call ab6_mixing_new(mix, modulo(in%iscf, 10), potden, &
          & AB6_MIXING_REAL_SPACE, npoints, in%nspin, 0, &
          & ierr, errmess, useprec = .false.)
  else if (in%itrpmax >1) then
     call ab6_mixing_new(mix, modulo(in%iscf, 10), potden, &
          & AB6_MIXING_REAL_SPACE, 1, in%nspin, 0, &
          & ierr, errmess, useprec = .false.)
  end if
  if (in%itrpmax >1) call ab6_mixing_eval_allocate(mix)
  endlooprp=.false.

  !if we are in the last_run case, validate the last_run only for the last cycle
  DoLastRunThings=(in%last_run == 1 .and. in%nrepmax == 0) !do the last_run things regardless of infocode

  rhopot_loop: do itrp=1,in%itrpmax
     !set the infocode to the value it would have in the case of no convergence
     infocode=1
     subd_loop : do icycle=1,in%nrepmax
        !if we are in the last_run case, validate the last_run only for the last cycle
        DoLastRunThings=(in%last_run == 1 .and. icycle == in%nrepmax) !do the last_run things regardless of infocode

        wfn_loop: do iter=1,in%itermax

           if (iproc == 0 .and. verbose > 0) then 
              write( *,'(1x,a,i0)') &
                   & repeat('-',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
           endif
           !control whether the minimisation iterations ended
           endloop= gnrm <= gnrm_cv .or. iter == in%itermax

           !control how many times the DIIS has switched into SD
           if (diis%idsx /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1

           !leave SD if the DIIS did not work the second time
           if (ndiis_sd_sw > 1) then
              diis%switchSD=.false.
           end if

           !stop the partial timing counter if necessary
           if (endloop .and. in%itrpmax==1) call timing(iproc,'WFN_OPT','PR')
           !logical flac for the self-consistent potential
           scpot=(in%itrpmax /= 1 .and. iter==1 .and. icycle==1) .or. & !mixing to be done
                (in%itrpmax == 1) .or. & !direct minimisation
                (itrp==1 .and. in%itrpmax/=1 .and. gnrm > in%gnrm_startmix) !startmix condition

           !calculate the self-consistent potential
           if (scpot) then
              ! Potential from electronic charge density
              call sumrho(iproc,nproc,orbs,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
                   n1i*n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)

              !here the density can be mixed
              if (mix%kind == AB6_MIXING_DENSITY .and. in%itrpmax>1) then
                 call mix_rhopot(iproc,nproc,mix%nfft*mix%nspden,in%alphamix,mix,&
                      & rhopot,itrp,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,hx*hy*hz,rpnrm)
                 if (iproc == 0 .and. itrp > 1) write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
                      'DENSITY iteration,Delta P (Norm 2/Volume)',itrp,rpnrm
                 endlooprp= (itrp > 1 .and. rpnrm <= in%rpnrm_cv) .or. itrp == in%itrpmax
                 ! DEBUGGING !!!!!!
                 rhopot = abs(rhopot) + 1.0d-20
              end if

              if(orbs%nspinor==4) then
                 !this wrapper can be inserted inside the poisson solver 
                 call PSolverNC(atoms%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,n3d,&
                      ixc,hxh,hyh,hzh,&
                      rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
              else
                 call XC_potential(atoms%geocode,'D',iproc,nproc,&
                      Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
                      rhopot,eexcu,vexcu,in%nspin,rhocore,potxc)

                 call H_potential(atoms%geocode,'D',iproc,nproc,&
                      Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
                      rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.,&
                      quiet=PSquiet) !optional argument

                 !sum the two potentials in rhopot array
                 !fill the other part, for spin, polarised
                 if (in%nspin == 2) then
                    call dcopy(Glr%d%n1i*Glr%d%n2i*n3p,rhopot(1),1,&
                         rhopot(1+n1i*n2i*n3p),1)
                 end if
                 !spin up and down together with the XC part
                 call axpy(Glr%d%n1i*Glr%d%n2i*n3p*in%nspin,1.0_dp,potxc(1,1,1,1),1,&
                      rhopot(1),1)

              end if

              !here the potential can be mixed
              if (mix%kind == AB6_MIXING_POTENTIAL .and. in%itrpmax>1) then
                 call mix_rhopot(iproc,nproc,mix%nfft*mix%nspden,in%alphamix,mix,&
                      & rhopot,itrp,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,hx*hy*hz,rpnrm)
                 if (iproc == 0 .and. itrp > 1) write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
                      'POTENTIAL iteration,Delta P (Norm 2/Volume)',itrp,rpnrm
                 endlooprp= (itrp > 1 .and. rpnrm <= in%rpnrm_cv) .or. itrp == in%itrpmax
              end if

           end if

           !temporary, to be corrected with comms structure
           if (in%exctxpar == 'OP2P') eexctX = -99.0_gp

           !allocate the potential in the full box
           call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,in%nspin,&
                orbs%norb,orbs%norbp,ngatherarr,rhopot,potential)

           call HamiltonianApplication(iproc,nproc,atoms,orbs,hx,hy,hz,rxyz,&
                nlpspd,proj,Glr,ngatherarr,potential,psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,&
                in%nspin,GPU,pkernel=pkernelseq)

           !deallocate potential
           call free_full_potential(nproc,potential,subname)

           energybs=ekin_sum+epot_sum+eproj_sum !the potential energy contains also exctX
           energy=energybs-ehart+eexcu-vexcu-eexctX+eion+edisp

           !check for convergence or whether max. numb. of iterations exceeded
           if (endloop) then
              if (gnrm < gnrm_cv) infocode=0
              exit wfn_loop 
           endif

           !control the previous value of idsx_actual
           idsx_actual_before=diis%idsx
           !previous value
           diis%energy_old=diis%energy
           !new value without the trace, to be added in hpsitopsi
           if (in%itrpmax >1) then
              diis%energy=0.0_gp
           else
              diis%energy=-ehart+eexcu-vexcu-eexctX+eion+edisp
           end if

           call hpsitopsi(iproc,nproc,orbs,hx,hy,hz,Glr,comms,ncong,&
                iter,diis,idsx,gnrm,gnrm_zero,trH,psi,psit,hpsi,in%nspin,GPU,in)

           tt=(energybs-trH)/trH
           if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
                (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
              !write this warning only if the system is closed shell
              call check_closed_shell(in%nspin,orbs,lcs)
              if (lcs) then
                 write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
                      'ERROR: inconsistency between gradient and energy',tt,energybs,trH
              end if
           endif
           if (iproc == 0) then
              if (verbose > 0 .and. in%itrpmax==1) then
                 write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
                      ekin_sum,epot_sum,eproj_sum
                 write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
              end if
              if (.not. scpot) then
                 if (gnrm_zero == 0.0_gp) then
                    write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter, tr(H),gnrm',iter,trH,gnrm
                 else
                    write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter, tr(H),gnrm,gnrm_zero',iter,trH,gnrm,gnrm_zero
                 end if
              else
                 if (gnrm_zero == 0.0_gp) then
                    write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
                 else
                    write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter,total energy,gnrm,gnrm_zero',iter,energy,gnrm,gnrm_zero
                 end if
              end if
           endif

           if (in%inputPsiId == 0) then
              if ((gnrm > 4.d0 .and. orbs%norbu /= orbs%norbd) .or. &
                   (orbs%norbu == orbs%norbd .and. gnrm > 10.d0)) then
                 if (iproc == 0) then
                    write( *,'(1x,a)')&
                         'ERROR: the norm of the residue is too large also with input wavefunctions.'
                 end if
                 infocode=3
                 call deallocate_before_exiting
                 return
              end if
           else if (in%inputPsiId == 1) then
              if (gnrm > 1.d0) then
                 if (iproc == 0) then
                    write( *,'(1x,a)')&
                         'The norm of the residue is too large, need to recalculate input wavefunctions'
                 end if
                 infocode=2
                 if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                 call deallocate_before_exiting
                 return
              end if
           end if
           !flush all writings on standart output
           call flush(6)
        end do wfn_loop

        if (iproc == 0) then 
           if (verbose > 1) write( *,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write( *,'(1x,a)') &
                '--------------------------------------------------- End of Wavefunction Optimisation'
           write( *,'(1x,a,3(1x,1pe18.11))') &
                'final  ekin,  epot,  eproj ',ekin_sum,epot_sum,eproj_sum
           write( *,'(1x,a,3(1x,1pe18.11))') &
                'final ehart, eexcu,  vexcu ',ehart,eexcu,vexcu
           if ((in%itrpmax >1 .and. endlooprp) .or. in%itrpmax == 1) then
              write(final_out, "(A5)") "FINAL"
           else
              write(final_out, "(A5)") "final"
           end if
           if (gnrm_zero == 0.0_gp) then
              write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') &
                   final_out // ' iter,total energy,gnrm',iter,energy,gnrm
           else
              write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') &
                   final_out // ' iter,total energy,gnrm,gnrm_zero',iter,energy,gnrm,gnrm_zero

           end if
           !write(61,*)hx,hy,hz,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
           if (in%itrpmax >1) then
              if ( diis%energy > diis%energy_min) write( *,'(1x,a,2(1pe9.2))')&
                   & 'WARNING: Found an energy value lower than the ' // final_out // &
                   & ' energy, delta:',diis%energy-diis%energy_min
           else
              !write this warning only if the system is closed shell
              call check_closed_shell(in%nspin,orbs,lcs)
              if (lcs) then
                 if ( energy > diis%energy_min) write( *,'(1x,a,2(1pe19.12))')&
                      'WARNING: Found an energy value lower than the FINAL energy, delta:',energy,diis%energy_min
              end if
           end if
        end if

        if (iter == in%itermax .and. iproc == 0 .and. infocode/=0) &
             write( *,'(1x,a)')'No convergence within the allowed number of minimization steps'

        call last_orthon(iproc,nproc,orbs,Glr%wfd,in%nspin,&
             comms,psi,hpsi,psit,evsum,.true.) !never deallocate psit and hpsi

        !exit if the infocode is correct
        if (infocode == 0) then
           exit subd_loop
        else
           if(iproc==0) then
              write(*,*)&
                   ' WARNING: Wavefunctions not converged after cycle',icycle
              if (icycle < in%nrepmax) write(*,*)' restart after diagonalisation'
           end if
           gnrm=1.d10
        end if

        if (in%itrpmax == 1 .and. in%norbsempty > 0) then
           !recalculate orbitals occupation numbers
           call evaltoocc(iproc,nproc,.false.,in%Tel,orbs)

           gnrm =1.d10
           diis%energy_min=1.d10
           diis%alpha=2.d0
        end if

     end do subd_loop

     if (in%itrpmax > 1) then
        !stop the partial timing counter if necessary
        if (endlooprp .and. in%itrpmax >1) then
           call timing(iproc,'WFN_OPT','PR')
           exit rhopot_loop
        end if

        !recalculate orbitals occupation numbers
        call evaltoocc(iproc,nproc,.false.,in%Tel,orbs)

        gnrm =1.d10
        diis%energy_min=1.d10
        diis%alpha=2.d0
     end if

  end do rhopot_loop

  !deallocate psit and hpsi if it was not already been done
  if (nproc > 1) then
     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  else
     nullify(psit)
  end if
  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)
  if (in%itrpmax > 1) then
     call ab6_mixing_deallocate(mix)
  end if

  if (in%inputPsiId /=-1000) then
     if (abs(evsum-energybs) > 1.d-8 .and. iproc==0) write( *,'(1x,a,2(1x,1pe20.13))')&
          'Difference:evsum,energybs',evsum,energybs
  end if



  if (in%idsx > 0) then
     call deallocate_diis_objects(diis,subname)
  end if

  !last run things has to be done:
  !if it is the last run and the infocode is zero
  !if infocode is not zero but the last run has been done for nrepmax times
  DoLastRunThings= (in%last_run == 1 .and. infocode == 0) .or. DoLastRunThings

  !analyse the possiblity to calculate Davidson treatment
  !(nvirt > 0 .and. in%inputPsiId == 0)
  DoDavidson= abs(in%norbv) > 0 .and. DoLastRunThings

  !project the wavefunctions on a gaussian basis and keep in memory
  if (in%gaussian_help) then
     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '---------------------------------------------------------- Gaussian Basis Projection'
     end if

     !extract the gaussian basis from the pseudowavefunctions
!!!     if (in%inputPsiId == 11) then
!!!        !extract the gaussian basis from the pseudowavefunctions
!!!        call gaussian_pswf_basis(21,.false.,iproc,atoms,rxyz,gbd)
!!!     else if (in%inputPsiId == 12) then
!!!        !extract the gaussian basis from the pseudopotential
!!!        call gaussian_psp_basis(atoms,rxyz,gbd)
!!!     end if

     !extract the gaussian basis from the pseudowavefunctions
     call gaussian_pswf_basis(21,.false.,iproc,in%nspin,atoms,rxyz,gbd,gbd_occ)

     if (associated(gbd_occ)) then
        i_all=-product(shape(gbd_occ))*kind(gbd_occ)
        deallocate(gbd_occ,stat=i_stat)
        call memocc(i_stat,i_all,'gbd_occ',subname)
        nullify(gbd_occ)
     end if


     if (.not. associated(gaucoeffs)) then
        allocate(gaucoeffs(gbd%ncoeff,orbs%norbp+ndebug),stat=i_stat)
        call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)
     end if

     allocate(thetaphi(2,gbd%nat+ndebug),stat=i_stat)
     call memocc(i_stat,thetaphi,'thetaphi',subname)
     thetaphi=0.0_gp

     call wavelets_to_gaussians(atoms%geocode,orbs%norbp,orbs%nspinor,&
          n1,n2,n3,gbd,thetaphi,&
          hx,hy,hz,Glr%wfd,psi,gaucoeffs)

     i_all=-product(shape(thetaphi))*kind(thetaphi)
     deallocate(thetaphi,stat=i_stat)
     call memocc(i_stat,i_all,'thetaphi',subname)

  end if

  !  write all the wavefunctions into files
  if (in%output_wf .and. DoLastRunThings) then
     !add flag for writing waves in the gaussian basis form
     if (in%gaussian_help) then

!!!        call gaussian_orthogonality(iproc,nproc,norb,norbp,gbd,gaucoeffs)
!!!
!!!        call gaussian_orthogonality(iproc,nproc,norb,norbp,gbd,gaucoeffs)
        !write the coefficients and the basis on a file
        call write_gaussian_information(iproc,nproc,orbs,gbd,gaucoeffs,'wavefunctions.gau')

        !build dual coefficients
        call dual_gaussian_coefficients(orbs%norbp,gbd,gaucoeffs)

        !control the accuracy of the expansion
        call check_gaussian_expansion(iproc,nproc,orbs,Glr,hx,hy,hz,psi,gbd,gaucoeffs)

        call deallocate_gwf(gbd,subname)
        i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
        deallocate(gaucoeffs,stat=i_stat)
        call memocc(i_stat,i_all,'gaucoeffs',subname)
        nullify(gbd%rxyz)

     else
        call  writemywaves(iproc,"wavefunction" // trim(wfformat), &
             & orbs,n1,n2,n3,hx,hy,hz,atoms,rxyz,Glr%wfd,psi)
     end if
  end if

  !plot the ionic potential, if required by output_grid
  if (in%output_grid == OUTPUT_GRID_DENSPOT .and. DoLastRunThings) then
     if (iproc == 0) write(*,*) 'writing ionic_potential' // gridformat
     call plot_density('ionic_potential' // gridformat,iproc,nproc,&
          n1,n2,n3,n1i,n2i,n3i,n3p,&
          1,hxh,hyh,hzh,atoms,rxyz,ngatherarr,pot_ion)
     if (iproc == 0) write(*,*) 'writing local_potential' // gridformat
     call plot_density('local_potential' // gridformat,iproc,nproc,&
          n1,n2,n3,n1i,n2i,n3i,n3p,&
          1,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rhopot)
  end if

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)
  if (counterions) then
     i_all=-product(shape(counter_ions))*kind(counter_ions)
     deallocate(counter_ions,stat=i_stat)
     call memocc(i_stat,i_all,'counter_ions',subname)
  end if

  if (inputpsi /= -1000) then
     !------------------------------------------------------------------------
     ! here we start the calculation of the forces
     if (iproc == 0) then
        write( *,'(1x,a)')&
             '----------------------------------------------------------------- Forces Calculation'
     end if

     ! Selfconsistent potential is saved in rhopot, 
     ! new arrays rho,pot for calculation of forces ground state electronic density

     ! Potential from electronic charge density

     !manipulate scatter array for avoiding the GGA shift
     do jproc=0,nproc-1
        !n3d=n3p
        nscatterarr(jproc,1)=nscatterarr(jproc,2)
        !i3xcsh=0
        nscatterarr(jproc,4)=0
     end do

     if (n3p>0) then
        allocate(rho(n1i*n2i*n3p*in%nspin+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     else
        allocate(rho(1+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     end if
     call sumrho(iproc,nproc,orbs,Glr,0,hxh,hyh,hzh,psi,rho,n1i*n2i*n3p,&
          nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)

     !plot the density on the density.pot file
     if ((in%output_grid >= OUTPUT_GRID_DENSITY .or. in%nvacancy /=0) .and. DoLastRunThings) then
        if (iproc == 0) write(*,*) 'writing electronic_density' // gridformat

        call plot_density('electronic_density' // gridformat,&
             iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,  & 
             in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rho)
        if (associated(rhocore)) then
           if (iproc == 0) write(*,*) 'writing grid core_density' // gridformat
           call plot_density('core_density' // gridformat,&
                iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,  & 
                1,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rhocore(1+n1i*n2i*i3xcsh:))
        end if
     end if
     !calculate the total density in the case of nspin==2
     if (in%nspin==2) then
        call axpy(n1i*n2i*n3p,1.0_dp,rho(1+n1i*n2i*n3p),1,rho(1),1)
     end if
     if (n3p>0) then
        allocate(pot(n1i,n2i,n3p,1+ndebug),stat=i_stat)
        call memocc(i_stat,pot,'pot',subname)
     else
        allocate(pot(1,1,1,1+ndebug),stat=i_stat)
        call memocc(i_stat,pot,'pot',subname)
     end if

     !calculate electrostatic potential
     call dcopy(n1i*n2i*n3p,rho,1,pot,1) 
     call H_potential(atoms%geocode,'D',iproc,nproc,&
          n1i,n2i,n3i,hxh,hyh,hzh,pot,pkernel,pot,ehart_fake,0.0_dp,.false.)

     !plot also the electrostatic potential
     if (in%output_grid == OUTPUT_GRID_DENSPOT .and. DoLastRunThings) then
        if (iproc == 0) write(*,*) 'writing hartree_potential' // gridformat
        call plot_density('hartree_potential' // gridformat, &
             & iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
             & in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,pot)
     end if

     allocate(gxyz(3,atoms%nat+ndebug),stat=i_stat)
     call memocc(i_stat,gxyz,'gxyz',subname)

     call timing(iproc,'Forces        ','ON')
     ! calculate local part of the forces gxyz
     call local_forces(iproc,atoms,rxyz,hxh,hyh,hzh,&
          n1,n2,n3,n3p,i3s+i3xcsh,n1i,n2i,n3i,rho,pot,gxyz)

     i_all=-product(shape(rho))*kind(rho)
     deallocate(rho,stat=i_stat)
     call memocc(i_stat,i_all,'rho',subname)
     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot',subname)

     if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)',advance='no')'Calculate nonlocal forces...'

     !refill projectors for tails, davidson
     refill_proj=(in%calc_tail .or. DoDavidson) .and. DoLastRunThings

     call nonlocal_forces(iproc,n1,n2,n3,hx,hy,hz,atoms,rxyz,&
          orbs,nlpspd,proj,Glr%wfd,psi,gxyz,refill_proj)

     if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)')'done.'

     ! Add up all the force contributions
     if (nproc > 1) then
        call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        do iat=1,atoms%nat
           fxyz(1,iat)=gxyz(1,iat)
           fxyz(2,iat)=gxyz(2,iat)
           fxyz(3,iat)=gxyz(3,iat)
        enddo
     end if

!!$  if (iproc == 0) then
!!$     sumx=0.d0 ; sumy=0.d0 ; sumz=0.d0
!!$     fumx=0.d0 ; fumy=0.d0 ; fumz=0.d0
!!$     do iat=1,atoms%nat
!!$        sumx=sumx+fxyz(1,iat) ; sumy=sumy+fxyz(2,iat) ; sumz=sumz+fxyz(3,iat)
!!$        fumx=fumx+fion(1,iat) ; fumy=fumy+fion(2,iat) ; fumz=fumz+fion(3,iat)
!!$     enddo
!!$     write(77,'(a30,3(1x,e10.3))') 'translat. force total pot ',sumx,sumy,sumz
!!$     write(77,'(a30,3(1x,e10.3))') 'translat. force ionic pot ',fumx,fumy,fumz
!!$  endif

     !add to the forces the ionic and dispersion contribution 
     do iat=1,atoms%nat
        fxyz(1,iat)=fxyz(1,iat)+fion(1,iat)+fdisp(1,iat)
        fxyz(2,iat)=fxyz(2,iat)+fion(2,iat)+fdisp(2,iat)
        fxyz(3,iat)=fxyz(3,iat)+fion(3,iat)+fdisp(3,iat)
     enddo

     i_all=-product(shape(gxyz))*kind(gxyz)
     deallocate(gxyz,stat=i_stat)
     call memocc(i_stat,i_all,'gxyz',subname)

     !clean the center mass shift and the torque in isolated directions
     call clean_forces(iproc,atoms,rxyz,fxyz,fnoise)

     call timing(iproc,'Forces        ','OF')
  end if

  i_all=-product(shape(fion))*kind(fion)
  deallocate(fion,stat=i_stat)
  call memocc(i_stat,i_all,'fion',subname)
  i_all=-product(shape(fdisp))*kind(fdisp)
  deallocate(fdisp,stat=i_stat)
  call memocc(i_stat,i_all,'fdisp',subname)

  !if (nvirt > 0 .and. in%inputPsiId == 0) then
  if (DoDavidson) then

     !for a band structure calculation allocate the array in which to put the eigenvalues
     if (associated(in%kptv)) then
        allocate(band_structure_eval(orbs%norbu+orbs%norbd+in%nspin*norbv,in%nkptv+ndebug),stat=i_stat)
        call memocc(i_stat,band_structure_eval,'band_structure_eval',subname)
     end if

     !calculate davidson procedure for all the groups of k-points which are chosen
     ikpt=1
     do igroup=1,in%ngroups_kptv

        ! Set-up number of states and shifting values.
        nvirtu = norbv
        nvirtd = 0
        if (in%nspin==2) nvirtd=nvirtu
        ! Create the orbitals.
        if (associated(in%kptv)) then
           nvirtu = nvirtu + orbs%norbu
           nvirtd = nvirtd + orbs%norbd
           nvirt  = nvirtu+nvirtd

           !number of k-points for this group
           nkptv = in%nkptsv_group(igroup) !size(in%kptv, 2)

           allocate(wkptv(nkptv+ndebug),stat=i_stat)
           call memocc(i_stat,wkptv,'wkptv',subname)
           wkptv(:) = real(1.0, gp) / real(nkptv, gp)

           call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
                & orbs%nspinor,nkptv,in%kptv(1,ikpt),wkptv,orbsv)
           !allocate communications arrays for virtual orbitals
           call orbitals_communicators(iproc,nproc,Glr,orbsv,commsv)  

           i_all=-product(shape(wkptv))*kind(wkptv)
           deallocate(wkptv,stat=i_stat)
           call memocc(i_stat,i_all,'wkptv',subname)

           !recreate the memory space for the projectors 
           call deallocate_proj_descr(nlpspd,subname)  
           i_all=-product(shape(proj))*kind(proj)
           deallocate(proj,stat=i_stat)
           call memocc(i_stat,i_all,'proj',subname)

           ! Calculate all projectors, or allocate array for on-the-fly calculation
           call timing(iproc,'CrtProjectors ','ON')
           call createProjectorsArrays(iproc,n1,n2,n3,rxyz,atoms,orbsv,&
                radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj) 
           call timing(iproc,'CrtProjectors ','OF') 

        else
           call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
                & orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsv)
           !allocate communications arrays for virtual orbitals
           call orbitals_communicators(iproc,nproc,Glr,orbsv,commsv)  

        end if

        !allocate psivirt pointer (note the orbs dimension)
        allocate(psivirt(orbsv%npsidim+ndebug),stat=i_stat)
        call memocc(i_stat,psivirt,'psivirt',subname)

        if (in%norbv < 0) then

           call direct_minimization(iproc,nproc,n1i,n2i,in,atoms,&
                orbs,orbsv,nvirt,Glr,comms,commsv,&
                hx,hy,hz,rxyz,rhopot,n3p,nlpspd,proj, &
                pkernelseq,psi,psivirt,ngatherarr,GPU)

        else if (in%norbv > 0) then
           call davidson(iproc,nproc,n1i,n2i,in,atoms,&
                orbs,orbsv,nvirt,Glr,comms,commsv,&
                hx,hy,hz,rxyz,rhopot,n3p,nlpspd,proj, &
                pkernelseq,psi,psivirt,ngatherarr,GPU)

        end if

        if (atoms%geocode == 'F' .and. .false.) then
           ! Potential from electronic charge density
           call sumrho(iproc,nproc,orbs,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
                n1i*n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)

           !Allocate second Exc derivative
           if (n3p >0) then
              allocate(dvxcdrho(n1i,n2i,n3p,in%nspin+1+ndebug),stat=i_stat)
              call memocc(i_stat,dvxcdrho,'dvxcdrho',subname)
           else
              allocate(dvxcdrho(1,1,1,in%nspin+1+ndebug),stat=i_stat)
              call memocc(i_stat,dvxcdrho,'dvxcdrho',subname)
           end if

           call XC_potential(atoms%geocode,'D',iproc,nproc,&
                Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
                rhopot,eexcu,vexcu,in%nspin,rhocore,potxc,dvxcdrho)

           !temporary call to the coupling matrix calculation
           allocate(psirocc(max(max(n1i*n2i*n3i*orbs%norbp,&
                ngatherarr(0,1)*orbs%norb),1)+ndebug),stat=i_stat)
           call memocc(i_stat,psirocc,'psirocc',subname)

           allocate(psirvirt(max(max(n1i*n2i*n3i*orbsv%norbp,&
                ngatherarr(0,1)*orbsv%norb),1)+ndebug),stat=i_stat)
           call memocc(i_stat,psirvirt,'psirvirt',subname)

           call prepare_psirocc(iproc,nproc,Glr,orbs,n3p,ngatherarr(0,1),psi,psirocc)

           call prepare_psirocc(iproc,nproc,Glr,orbsv,n3p,ngatherarr(0,1),psivirt,psirvirt)

           call center_of_charge(atoms,rxyz,chargec)

           call coupling_matrix_prelim(iproc,nproc,atoms%geocode,in%nspin,Glr,orbs,orbsv,&
                i3s+i3xcsh,n3p,hxh,hyh,hzh,chargec,pkernelseq,dvxcdrho,psirocc,psirvirt)

           i_all=-product(shape(psirocc))*kind(psirocc)
           deallocate(psirocc,stat=i_stat)
           call memocc(i_stat,i_all,'psirocc',subname)

           i_all=-product(shape(psirvirt))*kind(psirvirt)
           deallocate(psirvirt,stat=i_stat)
           call memocc(i_stat,i_all,'psirvirt',subname)


           i_all=-product(shape(dvxcdrho))*kind(dvxcdrho)
           deallocate(dvxcdrho,stat=i_stat)
           call memocc(i_stat,i_all,'dvxcdrho',subname)

        end if

        call deallocate_comms(commsv,subname)
        call deallocate_orbs(orbsv,subname)

        !in the case of band structure calculation, copy the values of the eigenvectors
        !into a new array to write them afterwards
        if (associated(in%kptv)) then
           call dcopy(orbsv%norb*nkptv,orbsv%eval(1),1,band_structure_eval(1,ikpt),1)
           !increment the value of ikpt
           ikpt=ikpt+in%nkptsv_group(igroup)
        end if

        i_all=-product(shape(orbsv%eval))*kind(orbsv%eval)
        deallocate(orbsv%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        i_all=-product(shape(psivirt))*kind(psivirt)
        deallocate(psivirt,stat=i_stat)
        call memocc(i_stat,i_all,'psivirt',subname)

     end do

     if (associated(in%kptv)) then
        !dump the band structure eigenvalue on a file and deallocate it
        
        i_all=-product(shape(band_structure_eval))*kind(band_structure_eval)
        deallocate(band_structure_eval,stat=i_stat)
        call memocc(i_stat,i_all,'band_structure_eval',subname)
     end if

  end if


  !perform here the mulliken charge and density of states
  !localise them on the basis of gatom of a number of atoms
  if (in%gaussian_help .and. DoLastRunThings) then
     !here one must check if psivirt should have been kept allocated
     if (.not. DoDavidson) then
        orbsv%norb=0
        orbsv%norbp=0
     end if
     call local_analysis(iproc,nproc,hx,hy,hz,in,atoms,rxyz,shift,Glr,orbs,orbsv,psi,psivirt)
  else if (DoLastRunThings .and. verbose > 2) then
     ! Do a full DOS calculation.
     call global_analysis(iproc, nproc, orbs)
  end if

  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'kernel',subname)

  if (in%exctxpar == 'OP2P') then
     i_all=-product(shape(pkernelseq))*kind(pkernelseq)
     deallocate(pkernelseq,stat=i_stat)
     call memocc(i_stat,i_all,'kernelseq',subname)
  end if



  !------------------------------------------------------------------------
  if (in%calc_tail .and. atoms%geocode == 'F' .and. DoLastRunThings ) then
     call timing(iproc,'Tail          ','ON')
     !    Calculate energy correction due to finite size effects
     !    ---reformat potential
     allocate(pot(n1i,n2i,n3i,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)

     if (nproc > 1) then
        call MPI_ALLGATHERV(rhopot,n1i*n2i*n3p,&
             mpidtypd,pot(1,1,1,1),ngatherarr(0,1),ngatherarr(0,2), & 
             mpidtypd,MPI_COMM_WORLD,ierr)
        !print '(a,2f12.6)','RHOup',sum(abs(rhopot(:,:,:,1))),sum(abs(pot(:,:,:,1)))
        if(in%nspin==2) then
           !print '(a,2f12.6)','RHOdw',sum(abs(rhopot(:,:,:,2))),sum(abs(pot(:,:,:,2)))
           call MPI_ALLGATHERV(rhopot(1+n1i*n2i*n3p),n1i*n2i*n3p,&
                mpidtypd,pot(1,1,1,2),ngatherarr(0,1),ngatherarr(0,2), & 
                mpidtypd,MPI_COMM_WORLD,ierr)
        end if
     else
        call dcopy(n1i*n2i*n3i*in%nspin,rhopot,1,pot,1)
     end if
     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
     deallocate(nscatterarr,stat=i_stat)
     call memocc(i_stat,i_all,'nscatterarr',subname)
     i_all=-product(shape(ngatherarr))*kind(ngatherarr)
     deallocate(ngatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'ngatherarr',subname)
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot',subname)
     i_all=-product(shape(potxc))*kind(potxc)
     deallocate(potxc,stat=i_stat)
     call memocc(i_stat,i_all,'potxc',subname)


     !pass hx instead of hgrid since we are only in free BC
     call CalculateTailCorrection(iproc,nproc,atoms,rbuf,orbs,&
          Glr,nlpspd,ncongt,pot,hx,rxyz,radii_cf,crmult,frmult,in%nspin,&
          proj,psi,(in%output_grid /= 0),ekin_sum,epot_sum,eproj_sum)

     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot',subname)

     !if (iproc==0) then
     !   open(61)
     !   write(61,'(4(f9.3),1x,7(1pe19.11))',advance='no')&
     !        hgrid,alat1,alat2,alat3,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
     !end if

     energybs=ekin_sum+epot_sum+eproj_sum
     energy=energybs-ehart+eexcu-vexcu+eion+edisp

     !if (iproc==0) then
     !   write(61,'(1pe19.11)')energy
     !   close(61)
     !end if

     if (iproc == 0) then
        write( *,'(1x,a,3(1x,1pe18.11))')&
             '  Corrected ekin,epot,eproj',ekin_sum,epot_sum,eproj_sum
        write( *,'(1x,a,1x,1pe24.17)')&
             'Total energy with tail correction',energy
     endif

     call timing(iproc,'Tail          ','OF')
  else
     !    No tail calculation
     if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot',subname)
     i_all=-product(shape(potxc))*kind(potxc)
     deallocate(potxc,stat=i_stat)
     call memocc(i_stat,i_all,'potxc',subname)
     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
     deallocate(nscatterarr,stat=i_stat)
     call memocc(i_stat,i_all,'nscatterarr',subname)
     i_all=-product(shape(ngatherarr))*kind(ngatherarr)
     deallocate(ngatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'ngatherarr',subname)
  endif
  ! --- End if of tail calculation

  call deallocate_before_exiting

contains

  !routine which deallocate the pointers and the arrays before exiting 
  subroutine deallocate_before_exiting

    !when this condition is verified we are in the middle of the SCF cycle
    if (infocode /=0 .and. infocode /=1 .and. in%inputPsiId /=-1000) then

       if (in%idsx > 0) then
          call deallocate_diis_objects(diis,subname)
       end if

       if (nproc > 1) then
          i_all=-product(shape(psit))*kind(psit)
          deallocate(psit,stat=i_stat)
          call memocc(i_stat,i_all,'psit',subname)
       end if

       i_all=-product(shape(hpsi))*kind(hpsi)
       deallocate(hpsi,stat=i_stat)
       call memocc(i_stat,i_all,'hpsi',subname)

       !free GPU if it is the case
       if (GPUconv .and. .not.(DoDavidson)) then
          call free_gpu(GPU,orbs%norbp)
       else if (OCLconv .and. .not.(DoDavidson)) then
          call free_gpu_OCL(GPU,orbs,in%nspin)
       end if

       i_all=-product(shape(pot_ion))*kind(pot_ion)
       deallocate(pot_ion,stat=i_stat)
       call memocc(i_stat,i_all,'pot_ion',subname)
       if (counterions) then
          i_all=-product(shape(counter_ions))*kind(counter_ions)
          deallocate(counter_ions,stat=i_stat)
          call memocc(i_stat,i_all,'counter_ions',subname)
       end if

       if (in%exctxpar == 'OP2P') then
          i_all=-product(shape(pkernelseq))*kind(pkernelseq)
          deallocate(pkernelseq,stat=i_stat)
          call memocc(i_stat,i_all,'kernelseq',subname)
       end if


       i_all=-product(shape(pkernel))*kind(pkernel)
       deallocate(pkernel,stat=i_stat)
       call memocc(i_stat,i_all,'kernel',subname)

       ! calc_tail false
       i_all=-product(shape(rhopot))*kind(rhopot)
       deallocate(rhopot,stat=i_stat)
       call memocc(i_stat,i_all,'rhopot',subname)
       i_all=-product(shape(potxc))*kind(potxc)
       deallocate(potxc,stat=i_stat)
       call memocc(i_stat,i_all,'potxc',subname)

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

    i_all=-product(shape(irrzon))*kind(irrzon)
    deallocate(irrzon,stat=i_stat)
    call memocc(i_stat,i_all,'irrzon',subname)

    i_all=-product(shape(phnons))*kind(phnons)
    deallocate(phnons,stat=i_stat)
    call memocc(i_stat,i_all,'phnons',subname)

    call deallocate_bounds(Glr%geocode,Glr%hybrid_on,Glr%bounds,subname)

    !free GPU if it is the case
    if (GPUconv .and. .not.(DoDavidson)) then
       call free_gpu(GPU,orbs%norbp)
    else if (OCLconv .and. .not.(DoDavidson)) then
       call free_gpu_OCL(GPU,orbs,in%nspin)
    end if

    call deallocate_comms(comms,subname)

    call deallocate_orbs(orbs,subname)
    call deallocate_atoms_scf(atoms,subname) 

    i_all=-product(shape(radii_cf))*kind(radii_cf)
    deallocate(radii_cf,stat=i_stat)
    call memocc(i_stat,i_all,'radii_cf',subname)

    call deallocate_proj_descr(nlpspd,subname)

    i_all=-product(shape(proj))*kind(proj)
    deallocate(proj,stat=i_stat)
    call memocc(i_stat,i_all,'proj',subname)

    !deallocate the core density if it has been allocated
    if(associated(rhocore)) then
       i_all=-product(shape(rhocore))*kind(rhocore)
       deallocate(rhocore,stat=i_stat)
       call memocc(i_stat,i_all,'rhocore',subname)
    end if

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

END SUBROUTINE cluster

