!!****f* BigDFT/call_bigdft
!! FUNCTION
!!   Routines to use bigdft as a blackbox
!! COPYRIGHT
!!   Copyright (C) 2005-2010 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
 subroutine call_bigdft(nproc,iproc,atoms,rxyz,in,energy,fxyz,rst,infocode)
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
     subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
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
       real(gp), intent(out) :: energy
       real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz_old
       real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
       real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
       real(wp), dimension(:), pointer :: psi
       real(wp), dimension(:,:), pointer :: gaucoeffs
     END SUBROUTINE cluster 
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

     call cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
          rst%psi,rst%Glr,rst%gaucoeffs,rst%gbd,rst%orbs,&
          rst%rxyz_old,rst%hx_old,rst%hy_old,rst%hz_old,in,rst%GPU,infocode)

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

END SUBROUTINE call_bigdft
!!***


!!****f* BigDFT/cluster
!! FUNCTION
!!  Main routine which does self-consistent loop.
!!  Does not parse input file and no geometry optimization.
!! DESCRIPTION
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
!! SOURCE
!!
subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
     psi,Glr,gaucoeffs,gbd,orbs,rxyz_old,hx_old,hy_old,hz_old,in,GPU,infocode)
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  use libxc_functionals
  use vdwcorrection, only: vdwcorrection_calculate_energy, vdwcorrection_calculate_forces, vdwcorrection_warnings
  use esatto
  use ab6_symmetry
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
  real(gp), intent(out) :: energy
  real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
  real(wp), dimension(:), pointer :: psi
  real(wp), dimension(:,:), pointer :: gaucoeffs
  !local variables
  character(len=*), parameter :: subname='cluster'
  character(len=3) :: PSquiet
  character(len=4) :: f4
  character(len=50) :: filename
  logical :: endloop,potion_overwritten=.false.,allfiles,onefile,refill_proj,DoDavidson
  integer :: ixc,ncong,idsx,ncongt,nspin,itermax,idsx_actual,idsx_actual_before,nsym
  integer :: nvirt,ndiis_sd_sw,norbv
  integer :: nelec,ndegree_ip,j,i,iorb
  integer :: n1_old,n2_old,n3_old,n3d,n3p,n3pi,i3xcsh,i3s,n1,n2,n3
  integer :: ncount0,ncount1,ncount_rate,ncount_max,n1i,n2i,n3i,i03,i04
  integer :: i1,i2,i3,ind,iat,i_all,i_stat,iter,ierr,jproc,ispin,inputpsi
  real :: tcpu0,tcpu1
  real(kind=8) :: crmult,frmult,cpmult,fpmult,gnrm_cv,rbuf,hxh,hyh,hzh,hx,hy,hz
  real(gp) :: peakmem,energy_old,sumz,evsum,sumx,sumy
  real(gp) :: eion,epot_sum,ekin_sum,eproj_sum,eexctX,ehart,eexcu,vexcu,alpha,gnrm
  real(gp) :: scprsum,energybs,tt,tel,ehart_fake,energy_min,psoffset
  real(kind=8) :: ttsum
  real(gp) :: edisp ! Dispersion energy
  type(wavefunctions_descriptors) :: wfd_old
  type(nonlocal_psp_descriptors) :: nlpspd
  type(communications_arrays) :: comms
  type(orbitals_data) :: orbsv
  type(gaussian_basis) :: Gvirt
  real(gp), dimension(3) :: shift
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(kind=8), dimension(:), allocatable :: rho
  real(kind=8), dimension(:,:), allocatable :: radii_cf,gxyz,fion,thetaphi
  real(gp), dimension(:,:),allocatable :: fdisp
  ! Charge density/potential,ionic potential, pkernel
  real(dp), dimension(:), allocatable :: pot_ion,rhopot
  real(gp), dimension(:), allocatable :: atchgs,radii
  real(kind=8), dimension(:,:,:,:), allocatable :: pot,rhoref,potxc
  real(kind=8), dimension(:), pointer :: pkernel,pkernel_ref
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(kind=8), dimension(:), pointer :: hpsi,psit,psi_old,psivirt,psidst,hpsidst
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: proj,gbd_occ,rhocore
  ! arrays for DIIS convergence accelerator
  real(kind=8), dimension(:,:,:), pointer :: ads
  ! Arrays for the symmetrisation.
  integer, dimension(:,:,:), allocatable :: irrzon
  real(dp), dimension(:,:,:), allocatable :: phnons


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
  itermax=in%itermax
  ncong=in%ncong
  idsx=in%idsx
  rbuf=in%rbuf
  ncongt=in%ncongt
  nspin=in%nspin

  norbv=in%norbv
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
          orbs%norb,nlpspd%nprojel,atoms%atomnames,0,in%nspin,peakmem)
  end if

  !allocate communications arrays
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
     do iorb=1,orbs%norb*orbs%nspinor
        write(f4,'(i4.4)')  iorb
        filename = 'wavefunction.'//f4
        inquire(file=filename,exist=onefile)
        allfiles=allfiles .and. onefile
        if (.not. allfiles) then
           if (iproc == 0) write(*,*)' WARNING: The wavefunction ',filename,&
                'does not exist, switch to normal input guess'
           inputpsi = 0
           exit
        end if
     end do
  end if

  ! INPUT WAVEFUNCTIONS, added also random input guess
  select case(inputpsi)
  case(-2)

     if (iproc == 0) then
        write( *,'(1x,a)')&
             '------------------------------------------------ Random wavefunctions initialization'
     end if

     !random initialisation of the wavefunctions
     allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

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

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,orbs,Glr%wfd,comms,psi,hpsi,psit)

  case(-1)

     !import gaussians form CP2K (data in files gaubasis.dat and gaucoeff.dat)
     !and calculate eigenvalues
     call import_gaussians(iproc,nproc,atoms,orbs,comms,&
          & Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,nlpspd,proj, &
          & pkernel,ixc,psi,psit,hpsi,nscatterarr,ngatherarr,in%nspin,&
          & atoms%symObj,irrzon,phnons)

  case(0)
     nspin=in%nspin
     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(iproc,nproc, atoms,&
          orbs,orbsv,norbv,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
          nlpspd,proj,pkernel,ixc,psi,hpsi,psit,psivirt,Gvirt,&
          nscatterarr,ngatherarr,nspin,0,atoms%symObj,irrzon,phnons,GPU)

  case(1)
     !these parts should be reworked for the non-collinear spin case

     !restart from previously calculated wavefunctions, in memory

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     if (iproc == 0) then
        write( *,'(1x,a)')&
             '-------------------------------------------------------------- Wavefunctions Restart'
     end if

     call reformatmywaves(iproc,orbs,atoms,hx_old,hy_old,hz_old,&
          n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,hx,hy,hz,n1,n2,n3,rxyz,Glr%wfd,psi)

     call deallocate_wfd(wfd_old,'cluster')

     i_all=-product(shape(psi_old))*kind(psi_old)
     deallocate(psi_old,stat=i_stat)
     call memocc(i_stat,i_all,'psi_old',subname)

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,orbs,Glr%wfd,comms,psi,hpsi,psit)

  case(2)
     !restart from previously calculated wavefunctions, on disk

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition

     allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     if (iproc == 0) then
        write( *,'(1x,a)')&
             '---------------------------------------------------- Reading Wavefunctions from disk'
     end if

     call readmywaves(iproc,orbs,n1,n2,n3,hx,hy,hz,atoms,rxyz_old,rxyz,Glr%wfd,psi)

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,orbs,Glr%wfd,comms,psi,hpsi,psit)

  case(11)
     !restart from previously calculated gaussian coefficients
     if (iproc == 0) then
        write( *,'(1x,a)')&
             '--------------------------------------- Quick Wavefunctions Restart (Gaussian basis)'
     end if

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     call restart_from_gaussians(iproc,nproc,orbs,Glr,hx,hy,hz,psi,gbd,gaucoeffs)

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,orbs,Glr%wfd,comms,psi,hpsi,psit)

  case(12)
     !reading wavefunctions from gaussian file
     if (iproc == 0) then
        write( *,'(1x,a)')&
             '------------------------------------------- Reading Wavefunctions from gaussian file'
     end if

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition

     allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

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

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,orbs,Glr%wfd,comms,psi,hpsi,psit)

  case default

!     if (iproc == 0) then
        write( *,'(1x,a)')'ERROR:values of inputPsiId must be integers from -2 to  2'
        write( *,'(1x,a)')'                                         or from 10 to 12'
        write( *,'(1x,a,i0)')'                               while we found',in%inputPsiId
!     end if
     stop

  end select

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
     allocate(psidst(sum(comms%ncntt(0:nproc-1))*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,psidst,'psidst',subname)
     allocate(hpsidst(sum(comms%ncntt(0:nproc-1))*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,hpsidst,'hpsidst',subname)
     allocate(ads(idsx+1,idsx+1,orbs%nkptsp*3+ndebug),stat=i_stat)
     call memocc(i_stat,ads,'ads',subname)
     call razero(orbs%nkptsp*3*(idsx+1)**2,ads)
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
  end if

  alpha=2.d0
  energy=1.d10
  gnrm=1.d10
  ekin_sum=0.d0 
  epot_sum=0.d0 
  eproj_sum=0.d0
  !minimum value of the energy during the minimisation procedure
  energy_min=1.d10
  !set the infocode to the value it would have in the case of no convergence
  infocode=1
  !local variable for the diis history
  idsx_actual=idsx
  !number of switching betweed DIIS and SD during self-consistent loop
  ndiis_sd_sw=0
  !previous value of idsx_actual to control if switching has appeared
  idsx_actual_before=idsx_actual

  !control whether there is a reference density
  potion_overwritten=.false.
 
  if (in%read_ref_den ) then

     !allocate the kernel for the reference density case
     call createKernel(iproc,nproc,'F',n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,pkernel_ref,&
          quiet=PSquiet)

     allocate(rhoref(n1i,n2i,max(n3d,1),in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhoref,'rhoref',subname)

     call read_potfile(atoms%geocode,'density.pot',n1,n2,n3,n1i,n2i,n3i,n3d,i3s,rhoref)

     potion_overwritten=.false.

  end if

  !end of the initialization part
  call timing(iproc,'INIT','PR')

  !Davidson is set to false first because used in deallocate_before_exiting
  DoDavidson= .false.

  wfn_loop: do iter=1,itermax

     if (iproc == 0 .and. verbose > 0) then 
        write( *,'(1x,a,i0)') &
             & repeat('-',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
     endif
     !control whether the minimisation iterations ended
     endloop= gnrm <= gnrm_cv .or. iter == itermax

     
     !control how many times the DIIS has switched into SD
     if (idsx_actual /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1

     !terminate SCF loop if forced to switch more than once from DIIS to SD
     endloop=endloop .or. ndiis_sd_sw > 2

     !stop the partial timing counter if necessary
     if (endloop) call timing(iproc,'WFN_OPT','PR')

     ! Potential from electronic charge density
     call sumrho(iproc,nproc,orbs,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
          n1i*n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)
!!$
!!$     open(78)
!!$     do i_stat=1,n3d
!!$        if (i3s+i3xcsh+i_stat-1 == 160) then
!!$           do i_all=1,n1i
!!$              write(78,*)i_all,rhopot(i_all,160,i_stat,1)!,rhopot(i_all,160,i_stat,2)
!!$           end do
!!$        end if
!!$     end do
!!$     close(78)


     if(orbs%nspinor==4) then
        !this wrapper can be inserted inside the poisson solver 
        call PSolverNC(atoms%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,n3d,&
             ixc,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
     else

        if (in%read_ref_den .and. gnrm <= in%gnrm_sw .or. potion_overwritten) then
           if (.not. potion_overwritten) then
              !overwrite pot_ion with the potential previously created
              call read_potfile(atoms%geocode,'potion_corr.pot',n1,n2,n3,n1i,n2i,n3i,n3pi,&
                   i3s+i3xcsh,pot_ion)

              if (.not. in%correct_offset) then
                 !read the ionic energy from disk
                 open(unit=22,file='eion_corr.tmp',status='unknown')
                 read(22,*)eion,ehart_fake
                 close(unit=22)
              end if
              potion_overwritten=.true.
           end if
           call correct_hartree_potential(atoms,iproc,nproc,&
                Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
                n3p,n3pi,n3d,i3s,i3xcsh,hxh,hyh,hzh,pkernel,ngatherarr,&
                rhoref,pkernel_ref,pot_ion,rhopot,ixc,in%nspin,ehart,eexcu,vexcu,PSquiet,&
                in%correct_offset)

        else

           call XC_potential(atoms%geocode,'D',iproc,nproc,&
                Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
                rhopot,eexcu,vexcu,in%nspin,rhocore,potxc)

           call H_potential(atoms%geocode,'D',iproc,nproc,&
                Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
                rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.,&
                quiet=PSquiet) !optional argument

           !atomic charge calculation with York's method, to be commented out
           if (endloop .and. in%gaussian_help .and. in%last_run == 1 .and. .false.) then

              allocate(atchgs(atoms%nat+ndebug),stat=i_stat)
              call memocc(i_stat,atchgs,'atchgs',subname)
              allocate(radii(atoms%nat+ndebug),stat=i_stat)
              call memocc(i_stat,radii,'radii',subname)

              call assign_atomic_radii(atoms,radii)

              !calculation of the atomic charges to compare wrt Mulliken
              !radii to be defined
              call atomic_charges(iproc,nproc,atoms%geocode,rxyz,radii,&
                   atoms%alat1,atoms%alat2,atoms%alat3,nelec,atoms%nat,Glr%d,&
                   hxh,hyh,hzh,n3p,i3s+i3xcsh,rhopot,atchgs)

              i_all=-product(shape(radii))*kind(radii)
              deallocate(radii,stat=i_stat)
              call memocc(i_stat,i_all,'radii',subname)
              i_all=-product(shape(atchgs))*kind(atchgs)
              deallocate(atchgs,stat=i_stat)
              call memocc(i_stat,i_all,'atchgs',subname)

           end if

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

     end if

     call HamiltonianApplication(iproc,nproc,atoms,orbs,hx,hy,hz,rxyz,&
          nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
          rhopot,psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,&
          in%nspin,GPU,pkernel=pkernel)

     energybs=ekin_sum+epot_sum+eproj_sum
     energy_old=energy
     energy=energybs-ehart+eexcu-vexcu-eexctX+eion+edisp

     !check for convergence or whether max. numb. of iterations exceeded
     if (endloop) then 
        if (iproc.eq.0) then 
           if (verbose > 1) write( *,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write( *,'(1x,a)') &
                '--------------------------------------------------- End of Wavefunction Optimisation'
           write( *,'(1x,a,3(1x,1pe18.11))') &
                'final  ekin,  epot,  eproj ',ekin_sum,epot_sum,eproj_sum
           write( *,'(1x,a,3(1x,1pe18.11))') &
                'final ehart, eexcu,  vexcu ',ehart,eexcu,vexcu
           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') &
                'FINAL iter,total energy,gnrm',iter,energy,gnrm
           !write(61,*)hx,hy,hz,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
           if (energy > energy_min) write( *,'(1x,a,1pe9.2)')&
                'WARNING: Found an energy value lower than the FINAL energy, delta:',energy-energy_min
        end if
        if (gnrm <= gnrm_cv) infocode=0
        exit wfn_loop 
     endif

     !control the previous value of idsx_actual
     idsx_actual_before=idsx_actual

     call hpsitopsi(iproc,nproc,orbs,hx,hy,hz,Glr,comms,ncong,&
          iter,idsx,idsx_actual,ads,energy,energy_old,energy_min,&
          alpha,gnrm,scprsum,psi,psit,hpsi,psidst,hpsidst,in%nspin,GPU)

     tt=(energybs-scprsum)/scprsum
     if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
          (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
        write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
     endif
     if (iproc.eq.0) then
        if (verbose > 0) then
           write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
                ekin_sum,epot_sum,eproj_sum
           write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
        end if
        write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
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
 
  end do wfn_loop
  if (iter == itermax .and. iproc == 0 ) &
       write( *,'(1x,a)')'No convergence within the allowed number of minimization steps'

  if (idsx_actual > 0) then
     i_all=-product(shape(psidst))*kind(psidst)
     deallocate(psidst,stat=i_stat)
     call memocc(i_stat,i_all,'psidst',subname)
     i_all=-product(shape(hpsidst))*kind(hpsidst)
     deallocate(hpsidst,stat=i_stat)
     call memocc(i_stat,i_all,'hpsidst',subname)
     i_all=-product(shape(ads))*kind(ads)
     deallocate(ads,stat=i_stat)
     call memocc(i_stat,i_all,'ads',subname)
  end if

  !analyse the possiblity to calculate Davidson rteatment
  !(nvirt > 0 .and. in%inputPsiId == 0)
  DoDavidson= in%norbv > 0 .and. (infocode==0 .or. in%nrepmax == 1) .and. in%last_run == 1
  
  call last_orthon(iproc,nproc,orbs,Glr%wfd,in%nspin,&
       comms,psi,hpsi,psit,evsum)


  if (abs(evsum-energybs) > 1.d-8 .and. iproc==0) write( *,'(1x,a,2(1x,1pe20.13))')&
       'Difference:evsum,energybs',evsum,energybs

  !project the wavefunctions on a gaussian basis and keep in memory
  if (in%gaussian_help) then
     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '---------------------------------------------------------- Gaussian Basis Projection'
     end if

     !extract the gaussian basis from the pseudowavefunctions
!!!     if (in%inputPsiId == 11) then
!!!        !extract the gaussian basis from the pseudowavefunctions
!!!        call gaussian_pswf_basis(21,iproc,atoms,rxyz,gbd)
!!!     else if (in%inputPsiId == 12) then
!!!        !extract the gaussian basis from the pseudopotential
!!!        call gaussian_psp_basis(atoms,rxyz,gbd)
!!!     end if

     !extract the gaussian basis from the pseudowavefunctions
     call gaussian_pswf_basis(21,iproc,in%nspin,atoms,rxyz,gbd,gbd_occ)

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
  if (in%output_wf .and. in%last_run==1) then
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
        call  writemywaves(iproc,orbs,n1,n2,n3,hx,hy,hz,atoms%nat,rxyz,Glr%wfd,psi)
        if (verbose >0) write( *,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
     end if
  end if

  !plot the ionic potential, if required by output_grid
  if (abs(in%output_grid)==2 .and. in%last_run==1) then
     if (in%output_grid==2) then
        if (iproc == 0) write(*,*) 'writing ionic_potential.pot'
        call plot_density(atoms%geocode,'ionic_potential.pot',iproc,nproc,&
             n1,n2,n3,n1i,n2i,n3i,n3p,&
             atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,pot_ion)
        if (iproc == 0) write(*,*) 'writing local_potential.pot'
        call plot_density(atoms%geocode,'local_potential.pot',iproc,nproc,&
             n1,n2,n3,n1i,n2i,n3i,n3p,&
             atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rhopot)
     else
        if (iproc == 0) write(*,*) 'writing ionic_potential.cube'
        call plot_density_cube(atoms%geocode,'ionic_potential',iproc,nproc,&
             n1,n2,n3,n1i,n2i,n3i,n3p,&
             in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,pot_ion)
        if (iproc == 0) write(*,*) 'writing local_potential.cube'
        call plot_density_cube(atoms%geocode,'local_potential',iproc,nproc,&
             n1,n2,n3,n1i,n2i,n3i,n3p,&
             in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rhopot)
     endif
  end if


  if (in%output_grid==3 .and. in%last_run==1) then
!!$        call plot_density(atoms%geocode,'b2B_xanes.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,1,&
!!$             atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rhopot)
!!$        write(comment,'(a)')'this file to check the positions and calculate shift '
!!$        call write_atomic_file("b2B_xanes",energy,rxyz,atoms,trim(comment))
  endif

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)

  !------------------------------------------------------------------------
  ! here we start the calculation of the forces
  if (iproc.eq.0) then
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
          nscatterarr,in%nspin,GPU, atoms%symObj, irrzon, phnons)

  !plot the density on the density.pot file
  if ((abs(in%output_grid) >= 1 .or. in%nvacancy /=0) .and. in%last_run==1) then
     if (in%output_grid >= 0) then
        if (in%nspin == 2 ) then
           if(iproc==0) write(*,*) 'ERROR: density cannot be plotted in .pot format for a spin-polarised calculation'
        else
           if (iproc == 0) write(*,*) 'writing electronic_density.pot'
           call plot_density(atoms%geocode,'electronic_density.pot',&
                iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
                atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rho)
           if (associated(rhocore)) then
              if (iproc == 0) write(*,*) 'writing grid core_density.pot'
              call plot_density(atoms%geocode,'core_density.pot',&
                   iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
                   atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,&
                   rhocore(1+n1i*n2i*i3xcsh))
           end if

        end if
     else 
        if (iproc == 0) write(*,*) 'writing electronic_density.cube'
        call plot_density_cube(atoms%geocode,'electronic_density',&
             iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,  & 
             in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rho)
        if (associated(rhocore)) then
           if (iproc == 0) write(*,*) 'writing grid core_density.cube'
           call plot_density_cube(atoms%geocode,'core_density',&
                iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,  & 
                1,hxh,hyh,hzh,atoms,rxyz,ngatherarr,rhocore(1+n1i*n2i*i3xcsh))
        end if
     endif
  end if
  !calculate the total density in the case of nspin==2
  if (in%nspin==2) then
     call axpy(n1i*n2i*n3p,1.0_dp,rho(1+n1i*n2i*n3p),1,rho(1),1)
!!$     do i3=1,n3p
!!$        do i2=1,n2i
!!$           do i1=1,n1i
!!$              ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!$              rho(ind)=rho(ind)+rho(ind+n1i*n2i*n3p)
!!$           end do
!!$        end do
!!$     end do
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
  if (abs(in%output_grid) == 2 .and. in%last_run==1) then
     if (in%output_grid == 2) then
        if (iproc.eq.0) write(*,*) 'writing hartree_potential.pot'
        call plot_density(atoms%geocode,'hartree_potential.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
             atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,pot)
     else
        if (iproc.eq.0) write(*,*) 'writing hartree_potential.cube'
        call plot_density_cube(atoms%geocode,'hartree_potential',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
             in%nspin,hxh,hyh,hzh,atoms,rxyz,ngatherarr,pot)
     end if
  end if

  if (in%read_ref_den) then
     i_all=-product(shape(pkernel_ref))*kind(pkernel_ref)
     deallocate(pkernel_ref,stat=i_stat)
     call memocc(i_stat,i_all,'pkernel_ref',subname)
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
  refill_proj=(in%calc_tail .or. DoDavidson) .and. in%last_run==1

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

  !add to the forces the ionic and dispersion contribution 
  do iat=1,atoms%nat
     fxyz(1,iat)=fxyz(1,iat)+fion(1,iat)+fdisp(1,iat)
     fxyz(2,iat)=fxyz(2,iat)+fion(2,iat)+fdisp(2,iat)
     fxyz(3,iat)=fxyz(3,iat)+fion(3,iat)+fdisp(3,iat)
  enddo

  i_all=-product(shape(fion))*kind(fion)
  deallocate(fion,stat=i_stat)
  call memocc(i_stat,i_all,'fion',subname)
  i_all=-product(shape(fdisp))*kind(fdisp)
  deallocate(fdisp,stat=i_stat)
  call memocc(i_stat,i_all,'fdisp',subname)
  i_all=-product(shape(gxyz))*kind(gxyz)
  deallocate(gxyz,stat=i_stat)
  call memocc(i_stat,i_all,'gxyz',subname)

  !subtraction of zero of the forces, disabled for the moment
  !the zero of the forces depends on the atomic positions
  if (in%gaussian_help .and. .false.) then
     sumx=0.d0
     sumy=0.d0
     sumz=0.d0
     do iat=1,atoms%nat
        sumx=sumx+fxyz(1,iat)
        sumy=sumy+fxyz(2,iat)
        sumz=sumz+fxyz(3,iat)
     enddo
     sumx=sumx/real(atoms%nat,gp)
     sumy=sumy/real(atoms%nat,gp)
     sumz=sumz/real(atoms%nat,gp)
     if (iproc==0) write( *,'(1x,a,1x,3(1x,1pe9.2))') &
          'Subtracting center-mass shift of',sumx,sumy,sumz

     do iat=1,atoms%nat
        fxyz(1,iat)=fxyz(1,iat)-sumx
        fxyz(2,iat)=fxyz(2,iat)-sumy
        fxyz(3,iat)=fxyz(3,iat)-sumz
     enddo
  end if

  call timing(iproc,'Forces        ','OF')

  !if (nvirt > 0 .and. in%inputPsiId == 0) then
  if (DoDavidson) then

!!$     !allocate the virtual orbitals descriptors
!!$     !allocated previously
!!$     call deallocate_orbs(orbsv,subname)
!!$     
!!$     !create orbitals equal to the other
!!$     call orbitals_descriptors(iproc,nproc,orbs%norb,orbs%norbu,orbs%norbd, &
!!$          & orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsv)
!!$     nvirt=orbs%norb !temporary
!!$
!!$     !allocate psivirt pointer (note the orbs dimension)
!!$     allocate(psivirt(orbs%npsidim+ndebug),stat=i_stat)
!!$     call memocc(i_stat,psivirt,'psivirt',subname)

     !allocate psivirt pointer (note the orbs dimension)
     allocate(psivirt(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psivirt,'psivirt',subname)

     call davidson(iproc,nproc,n1i,n2i,in,atoms,&
          orbs,orbsv,nvirt,Glr,comms,&
          hx,hy,hz,rxyz,rhopot,i3xcsh,n3p,nlpspd,proj, &
          pkernel,psi,psivirt,ngatherarr,GPU)
  end if


  !perform here the mulliken charge and density of states
  !localise them on the basis of gatom of a number of atoms
  if (in%gaussian_help .and. in%last_run==1) then
     call local_analysis(iproc,nproc,hx,hy,hz,shift,Glr,orbs,orbsv,psi,psivirt)
  end if



  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'kernel',subname)


  !------------------------------------------------------------------------
  if (in%calc_tail .and. atoms%geocode == 'F' .and. in%last_run==1 ) then
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

     
     if (in%read_ref_den) then
        i_all=-product(shape(rhoref))*kind(rhoref)
        deallocate(rhoref,stat=i_stat)
        call memocc(i_stat,i_all,'rhoref',subname)
     end if
     
     !pass hx instead of hgrid since we are only in free BC
     call CalculateTailCorrection(iproc,nproc,atoms,rbuf,orbs,&
          Glr,nlpspd,ncongt,pot,hx,rxyz,radii_cf,crmult,frmult,in%nspin,&
          proj,psi,in%output_grid,ekin_sum,epot_sum,eproj_sum)
     
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
     if (in%read_ref_den) then
        i_all=-product(shape(rhoref))*kind(rhoref)
        deallocate(rhoref,stat=i_stat)
        call memocc(i_stat,i_all,'rhoref',subname)
     end if
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
    if (infocode /=0 .and. infocode /=1) then
       
       if (idsx_actual > 0) then
          i_all=-product(shape(psidst))*kind(psidst)
          deallocate(psidst,stat=i_stat)
          call memocc(i_stat,i_all,'psidst',subname)
          i_all=-product(shape(hpsidst))*kind(hpsidst)
          deallocate(hpsidst,stat=i_stat)
          call memocc(i_stat,i_all,'hpsidst',subname)
          i_all=-product(shape(ads))*kind(ads)
          deallocate(ads,stat=i_stat)
          call memocc(i_stat,i_all,'ads',subname)
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
          call free_gpu_OCL(GPU,orbs%norbp)
       end if
       
       i_all=-product(shape(pot_ion))*kind(pot_ion)
       deallocate(pot_ion,stat=i_stat)
       call memocc(i_stat,i_all,'pot_ion',subname)
       
       i_all=-product(shape(pkernel))*kind(pkernel)
       deallocate(pkernel,stat=i_stat)
       call memocc(i_stat,i_all,'kernel',subname)
       if (in%read_ref_den) then
          i_all=-product(shape(pkernel_ref))*kind(pkernel_ref)
          deallocate(pkernel_ref,stat=i_stat)
          call memocc(i_stat,i_all,'pkernel_ref',subname)
       end if
       
       ! calc_tail false
       i_all=-product(shape(rhopot))*kind(rhopot)
       deallocate(rhopot,stat=i_stat)
       call memocc(i_stat,i_all,'rhopot',subname)
       i_all=-product(shape(potxc))*kind(potxc)
       deallocate(potxc,stat=i_stat)
       call memocc(i_stat,i_all,'potxc',subname)
       if (in%read_ref_den) then
          i_all=-product(shape(rhoref))*kind(rhoref)
          deallocate(rhoref,stat=i_stat)
          call memocc(i_stat,i_all,'rhoref',subname)
       end if
       
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
    if (DoDavidson) then
       !call deallocate_gwf(Gvirt,subname)

       call deallocate_orbs(orbsv,subname)

       i_all=-product(shape(orbsv%eval))*kind(orbsv%eval)
       deallocate(orbsv%eval,stat=i_stat)
       call memocc(i_stat,i_all,'eval',subname)

       i_all=-product(shape(psivirt))*kind(psivirt)
       deallocate(psivirt,stat=i_stat)
       call memocc(i_stat,i_all,'psivirt',subname)
    end if
    
    i_all=-product(shape(irrzon))*kind(irrzon)
    deallocate(irrzon,stat=i_stat)
    call memocc(i_stat,i_all,'irrzon',subname)

    i_all=-product(shape(phnons))*kind(phnons)
    deallocate(phnons,stat=i_stat)
    call memocc(i_stat,i_all,'phnons',subname)

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
    
    !free GPU if it is the case
    if (GPUconv .and. .not.(DoDavidson)) then
       call free_gpu(GPU,orbs%norbp)
    else if (OCLconv .and. .not.(DoDavidson)) then
       call free_gpu_OCL(GPU,orbs%norbp)
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
!!***

