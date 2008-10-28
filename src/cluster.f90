!!****f* BigDFT/call_bigdft
!! NAME
!!  call_bigdft
!!
!! FUNCTION
!!  Routines to use bigdft as a blackbox
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 
!!
!!
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
  integer :: i_stat,i_all,ierr,inputPsiId_orig,icycle
  !temporary interface
  interface
     subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
          psi,wfd,gaucoeffs,gbd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nproc,iproc
       integer, intent(inout) :: n1,n2,n3,norbp,norb
       integer, intent(out) :: infocode
       type(input_variables), intent(in) :: in
       type(wavefunctions_descriptors), intent(inout) :: wfd
       type(atoms_data), intent(inout) :: atoms
       type(gaussian_basis), intent(inout) :: gbd
       real(kind=8), intent(out) :: energy
       real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz_old
       real(kind=8), dimension(3,atoms%nat), target, intent(inout) :: rxyz
       real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz
       real(kind=8), dimension(:), pointer :: eval
       real(kind=8), dimension(:), pointer :: psi
       real(kind=8), dimension(:,:), pointer :: gaucoeffs
     end subroutine cluster 
  end interface

  inputPsiId_orig=in%inputPsiId

  loop_cluster: do icycle=1,10

     if (in%inputPsiId == 0 .and. associated(rst%psi)) then
        i_all=-product(shape(rst%psi))*kind(rst%psi)
        deallocate(rst%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(rst%eval))*kind(rst%eval)
        deallocate(rst%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(rst%wfd,subname)
     end if

     call cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
          rst%psi,rst%wfd,rst%gaucoeffs,rst%gbd,rst%norbp,rst%norb,rst%eval,&
          rst%n1,rst%n2,rst%n3,rst%rxyz_old,in,infocode)

     if (in%inputPsiId==1 .and. infocode==2) then
        if (in%gaussian_help) then
           in%inputPsiId=11
        else
           in%inputPsiId=0
        end if
     else if (in%inputPsiId==1 .and. infocode==1) then
        !in%inputPsiId=0 !better to diagonalise that to restart an input guess
        if(iproc==0)write(*,*)' WARNING: Wavefunctions converged after cycle',icycle 
     else if (in%inputPsiId == 0 .and. infocode==3) then
        if (iproc.eq.0) then
           write( *,'(1x,a)')'Convergence error, cannot proceed.'
           write( *,'(1x,a)')' writing positions in file posout_999.xyz then exiting'

           call wtposout(999,energy,rxyz,atoms)

        end if

        i_all=-product(shape(rst%psi))*kind(rst%psi)
        deallocate(rst%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(rst%eval))*kind(rst%eval)
        deallocate(rst%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(rst%wfd,subname)
        !finalize memory counting (there are still the positions and the forces allocated)
        call memocc(0,0,'count','stop')

        if (nproc > 1) call MPI_FINALIZE(ierr)

        stop
     else
        exit loop_cluster
     end if

  end do loop_cluster

  !preserve the previous value
  in%inputPsiId=inputPsiId_orig

end subroutine call_bigdft
!!***


!!****f* BigDFT/cluster
!! NAME
!!  cluster
!!
!! FUNCTION
!!  Main routine which does self-consistent loop.
!!  Do not parse input file and no geometry optimization.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 
!!
!!
!! SOURCE
!!
subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
     psi,wfd,gaucoeffs,gbd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
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
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: n1,n2,n3,norbp,norb
  integer, intent(out) :: infocode
  type(input_variables), intent(in) :: in
  type(wavefunctions_descriptors), intent(inout) :: wfd
  type(atoms_data), intent(inout) :: atoms
  type(gaussian_basis), intent(inout) :: gbd
  real(kind=8), intent(out) :: energy
  real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz_old
  real(kind=8), dimension(3,atoms%nat), target, intent(inout) :: rxyz
  real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz
  real(kind=8), dimension(:), pointer :: eval
  real(kind=8), dimension(:), pointer :: psi
  real(kind=8), dimension(:,:), pointer :: gaucoeffs
  !local variables
  character(len=*), parameter :: subname='cluster'
  character(len=10) :: orbname
  logical :: calc_tail,switchSD
  integer :: ixc,ncharge,ncong,idsx,ncongt,nspin,mpol,itermax,idsx_actual,nvirte,nvirtep,nvirt
  integer :: nelec,norbu,norbd,ndegree_ip,nvctrp,mids,iorb,iounit,ids,idiistol,j,norbpeff
  integer :: n1_old,n2_old,n3_old,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n3d,n3p,n3pi,i3xcsh,i3s
  integer :: ncount0,ncount1,ncount_rate,ncount_max,iunit,n1i,n2i,n3i,nl1,nl2,nl3
  integer :: i1,i2,i3,ind,iat,ierror,i_all,i_stat,iter,ierr,i03,i04,jproc,ispin,nspinor,nplot
  real :: tcpu0,tcpu1
  real(kind=8) :: hgrid,crmult,frmult,cpmult,fpmult,elecfield,gnrm_cv,rbuf,hx,hy,hz,hxh,hyh,hzh
  real(kind=8) :: peakmem,gnrm_check,hgrid_old,energy_old,sumz
  real(kind=8) :: eion,epot_sum,ekin_sum,eproj_sum,ehart,eexcu,vexcu,alpha,gnrm,evsum,sumx,sumy
  real(kind=8) :: scprsum,energybs,tt,tel,eexcu_fake,vexcu_fake,ehart_fake,energy_min,psoffset
  real(kind=8) :: factor,rhon,rhos,ttsum,hx_old,hy_old,hz_old
  type(wavefunctions_descriptors) :: wfd_old
  type(convolutions_bounds) :: bounds
  type(nonlocal_psp_descriptors) :: nlpspd
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(kind=8), dimension(:), allocatable :: occup,spinsgn,spinsgn_foo,rho
  real(kind=8), dimension(:,:), allocatable :: radii_cf,gxyz,fion,thetaphi
  ! Charge density/potential,ionic potential, pkernel
  real(kind=8), dimension(:), allocatable :: pot_ion
  real(kind=8), dimension(:,:,:,:), allocatable :: rhopot,pot,rho_diag
  real(kind=8), dimension(:,:,:), allocatable :: m_norm
  real(kind=8), dimension(:), pointer :: pkernel
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(kind=8), dimension(:), pointer :: hpsi,psit,psi_old,psivirt,psidst,hpsidst
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: proj
  ! arrays for DIIS convergence accelerator
  real(kind=8), dimension(:,:,:), pointer :: ads
  ! tmp debug array
  real(kind=8), dimension(:,:), allocatable :: tmred
  
!*****added by Alexey**********************************************************************	   
  integer,parameter::lupfil=14
  logical::hybrid_on
!******************************************************************************************  

  !copying the input variables for readability
  !this section is of course not needed
  !note that this procedure is convenient ONLY in the case of scalar variables
  !an array would have been copied, thus occupying more memory space
  !Hence WARNING: these variables are copied, in case of an update the new value should be 
  !reassigned inside the structure

  hgrid=in%hgrid
  crmult=in%crmult
  frmult=in%frmult
  cpmult=in%frmult
  fpmult=in%frmult
  ixc=in%ixc
  ncharge=in%ncharge
  elecfield=in%elecfield
  gnrm_cv=in%gnrm_cv
  itermax=in%itermax
  ncong=in%ncong
  idsx=in%idsx
  calc_tail=in%calc_tail
  rbuf=in%rbuf
  ncongt=in%ncongt
  nspin=in%nspin
  if(nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if
  mpol=in%mpol

  nvirt=in%nvirt
  nplot=in%nplot

  hx=in%hgrid
  hy=in%hgrid
  hz=in%hgrid

!!$  geocode=atoms%geocode
!!$  alat1=atoms%alat1
!!$  alat2=atoms%alat2
!!$  alat3=atoms%alat3

  if (iproc.eq.0) then
     write( *,'(1x,a,1x,i0)') &
       '===================== BigDFT Wavefunction Optimization =============== inputPsiId=',&
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


  ! We save the variables that defined the previous psi if the restart is active
  if (in%inputPsiId == 1) then
     !regenerate grid spacings
     if (atoms%geocode == 'P') then
        call correct_grid(atoms%alat1,hx,n1)
        call correct_grid(atoms%alat2,hy,n2)
        call correct_grid(atoms%alat3,hz,n3)
     else if (atoms%geocode == 'S') then 
        call correct_grid(atoms%alat1,hx,n1)
        call correct_grid(atoms%alat3,hz,n3)
     end if
     call copy_old_wavefunctions(iproc,nproc,norb,norbp,nspinor,hx,hy,hz,n1,n2,n3,wfd,psi,&
          hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,wfd_old,psi_old)
  else if (in%inputPsiId == 11) then
     !deallocate wavefunction and descriptors for placing the gaussians
     
     call deallocate_wfd(wfd,subname)

     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)

  end if

  if(nspin/=1 .and. nspin/=2 .and. nspin/=4) nspin=1
  if(nspin==1) mpol=0

  ! grid spacing (same in x,y and z direction)

  if (iproc==0) then
     write( *,'(1x,a)')&
          '------------------------------------------------------------------ System Properties'
  end if

  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call read_system_variables(iproc,nproc,in,atoms,radii_cf,nelec,norb,norbu,norbd,norbp,iunit)

  allocate(occup(norb+ndebug),stat=i_stat)
  call memocc(i_stat,occup,'occup',subname)
  allocate(spinsgn(norb+ndebug),stat=i_stat)
  call memocc(i_stat,spinsgn,'spinsgn',subname)

  ! Occupation numbers
  call input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinsgn)

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,&
       n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)

!*****added by Alexey**********************************************************************	   
  hybrid_on=               (nfu1-nfl1+lupfil.lt.n1+1)
  hybrid_on=(hybrid_on.and.(nfu2-nfl2+lupfil.lt.n2+1))
  hybrid_on=(hybrid_on.and.(nfu3-nfl3+lupfil.lt.n3+1))
!******************************************************************************************  
		
  hxh=0.5d0*hx
  hyh=0.5d0*hy
  hzh=0.5d0*hz

  !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
  ndegree_ip=16 !default value to be put to 16 and update references for test
  call createKernel(atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,iproc,nproc,pkernel)

  ! Create wavefunctions descriptors and allocate them
  call timing(iproc,'CrtDescriptors','ON')
  call createWavefunctionsDescriptors(iproc,nproc,n1,n2,n3,in%output_grid,hx,hy,hz,&
       atoms,rxyz,radii_cf,crmult,frmult,wfd,&
       nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,nspinor,hybrid_on)
  call timing(iproc,'CrtDescriptors','OF')

  ! Calculate all projectors, or allocate array for on-the-fly calculation
  call timing(iproc,'CrtProjectors ','ON')
  call createProjectorsArrays(iproc,n1,n2,n3,rxyz,atoms,&
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
  call timing(iproc,'CrtProjectors ','OF')

  !memory estimation
  if (iproc==0) then
     call MemoryEstimator(atoms%geocode,nproc,idsx,n1,n2,n3,atoms%alat1,atoms%alat2,atoms%alat3,&
          hx,hy,hz,atoms%nat,atoms%ntypes,atoms%iatype,rxyz,radii_cf,crmult,frmult,norb,&
          nlpspd%nprojel,atoms%atomnames,.false.,nspin,peakmem)
  end if

  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  !allocate array for the communications of the potential
  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ngatherarr,'ngatherarr',subname)

  !create the descriptors for the density and the potential
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

  !here calculate the ionic energy and forces accordingly
  allocate(fion(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fion,'fion',subname)

  call IonicEnergyandForces(iproc,nproc,atoms,hxh,hyh,hzh,rxyz,eion,fion,&
       psoffset,n1,n2,n3,n1i,n2i,n3i,i3s+i3xcsh,n3pi,pot_ion,pkernel)

  !can pass the atoms data structure as argument
  call createIonicPotential(atoms%geocode,iproc,nproc,atoms%nat,atoms%ntypes,atoms%iatype,&
       atoms%psppar,atoms%nelpsp,rxyz,hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s+i3xcsh,&
       n1i,n2i,n3i,pkernel,pot_ion,eion,psoffset)

  !Allocate Charge density, Potential in real space
  if (n3d >0) then
     allocate(rhopot(n1i,n2i,n3d,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  else
     allocate(rhopot(1,1,1,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  end if

  !avoid allocation of the eigenvalues array in case of restart
  if (in%inputPsiId /= 1 .and. in%inputPsiId /= 11) then
     allocate(eval(norb+ndebug),stat=i_stat)
     call memocc(i_stat,eval,'eval',subname)
  end if

  ! INPUT WAVEFUNCTIONS, added also random input guess
  if (in%inputPsiId == -2) then

     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '------------------------------------------------ Random wavefunctions initialization'
     end if

     !random initialisation of the wavefunctions
     allocate(psi(nvctrp*nspinor*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     psi=0.0d0
     ttsum=0.0d0
     do iorb=1,norbp*nproc*max(1,nspinor)
        if(mod(iorb-1,nspinor)==0) then
           do i1=1,nvctrp
              do j=0,iproc-1
                 call random_number(tt)
              end do
              call random_number(tt)
              psi(i1+nvctrp*(iorb-1))=real(tt,kind=8)*0.01d0
              ttsum=ttsum+psi(i1+nvctrp*(iorb-1))
              do j=iproc+1,nproc
                 call random_number(tt)
              end do
           end do
        end if
     end do
     !write( *,'(a,30f10.4)') 'Rand Check',ttsum,(sum(psi(:,iorb)),iorb=1,norbp*nproc*nspinor)
 
     eval(:)=-0.5d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else if (in%inputPsiId == -1) then

     !import gaussians form CP2K (data in files gaubasis.dat and gaucoeff.dat)
     !and calculate eigenvalues
     call import_gaussians(iproc,nproc,cpmult,fpmult,radii_cf,atoms,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          norb,norbp,occup,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,wfd,bounds,nlpspd,proj, &
          pkernel,ixc,psi,psit,hpsi,eval,nscatterarr,ngatherarr,nspin,spinsgn,hybrid_on)

  else if (in%inputPsiId == 0) then 

     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(iproc,nproc,cpmult,fpmult,radii_cf,atoms,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          norb,norbp,nvirte,nvirtep,nvirt,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,&
          wfd,bounds,nlpspd,proj,pkernel,ixc,psi,hpsi,psit,psivirt,eval,&
          nscatterarr,ngatherarr,nspin,spinsgn,hybrid_on)
  
  else if (in%inputPsiId == 1) then 
     !these parts should be reworked for the non-collinear spin case

     !restart from previously calculated wavefunctions, in memory

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(psi(nvctrp*norbp*nproc*nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '-------------------------------------------------------------- Wavefunctions Restart'
     end if

     call reformatmywaves(iproc,norb*nspinor,norbp*nspinor,atoms%nat,hx_old,hy_old,hz_old,&
          n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)

     call deallocate_wfd(wfd_old,'cluster')

     i_all=-product(shape(psi_old))*kind(psi_old)
     deallocate(psi_old,stat=i_stat)
     call memocc(i_stat,i_all,'psi_old',subname)

     !initialise control value for gnrm in the case of a restart
     gnrm_check=0.d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else if (in%inputPsiId == 2 ) then 
     !restart from previously calculated wavefunctions, on disk

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(psi(nvctrp*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '---------------------------------------------------- Reading Wavefunctions from disk'
     end if

     call readmywaves(iproc,norb,norbp,n1,n2,n3,hx,hy,hz,atoms%nat,rxyz,wfd,psi,eval)

     !initialise control value for gnrm in the case of a restart
     gnrm_check=0.d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else if (in%inputPsiId == 11 ) then 
     !restart from previously calculated gaussian coefficients
     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '--------------------------------------- Quick Wavefunctions Restart (Gaussian basis)'
     end if

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(psi(nvctrp*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     call restart_from_gaussians(atoms%geocode,iproc,nproc,norb,norbp,&
     n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,hy,hz,wfd,psi,gbd,gaucoeffs)

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else if (in%inputPsiId == 12 ) then 
     !reading wavefunctions from gaussian file
     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '------------------------------------------- Reading Wavefunctions from gaussian file'
     end if
     
     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(psi(nvctrp*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     call read_gaussian_information(iproc,nproc,norb,norbp,gbd,gaucoeffs,eval,'wavefunctions.gau')
     !associate the new positions, provided that the atom number is good
     if (gbd%nat == atoms%nat) then
        gbd%rxyz=>rxyz
     else
        if (iproc == 0) then
           write( *,*)&
                ' ERROR: the atom number does not coincide with the number of gaussian centers'
        end if
        stop
     end if
 
     call restart_from_gaussians(atoms%geocode,iproc,nproc,norb,norbp,&
     n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,hy,hz,wfd,psi,gbd,gaucoeffs)

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else

     if (iproc == 0) then
        write( *,'(1x,a)')'ERROR:values of inputPsiId must be integers from -2 to  2'
        write( *,'(1x,a)')'                                         or from 10 to 12'
        write( *,'(1x,a,i0)')'                               while we found',in%inputPsiId
     end if
     stop

  end if

  !this rearrange the value of norbu, to be changed
  if(nspinor==4) then
     norbu=norb
     norbd=0
  end if

  !save the new atomic positions in the rxyz_old array
  do iat=1,atoms%nat
     rxyz_old(1,iat)=rxyz(1,iat)
     rxyz_old(2,iat)=rxyz(2,iat)
     rxyz_old(3,iat)=rxyz(3,iat)
  enddo

  ! allocate arrays necessary for DIIS convergence acceleration
  if (idsx > 0) then
     allocate(psidst(nvctrp*nspinor*norbp*nproc*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,psidst,'psidst',subname)
     allocate(hpsidst(nvctrp*nspinor*norbp*nproc*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,hpsidst,'hpsidst',subname)
     allocate(ads(idsx+1,idsx+1,3+ndebug),stat=i_stat)
     call memocc(i_stat,ads,'ads',subname)
     call razero(3*(idsx+1)**2,ads)
  endif

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

!!$  !logical control variable for switch DIIS-SD
!!$  switchSD=.false.
!!$  ids=0
!!$  idiistol=0

  !end of the initialization part
  call timing(iproc,'INIT','PR')

  ! loop for wavefunction minimization
  wfn_loop: do iter=1,itermax
!!$     if (idsx > 0) then
!!$        mids=mod(ids,idsx)+1
!!$        ids=ids+1
!!$     end if
     if (iproc == 0) then 
        write( *,'(1x,a,i0)')&
           '---------------------------------------------------------------------------- iter= ',&
           iter
     endif

     !control whether the minimisation iterations ended and stop the partial timing counter
     if (gnrm <= gnrm_cv .or. iter == itermax) call timing(iproc,'WFN_OPT','PR')

     ! Potential from electronic charge density
     call sumrho(atoms%geocode,iproc,nproc,norb,norbp,ixc,n1,n2,n3,hxh,hyh,hzh,occup,  & 
     wfd,psi,rhopot,n1i*n2i*n3d,nscatterarr,nspin,nspinor,spinsgn,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,hybrid_on)
     
     if(nspinor==4) then
        !this wrapper can be inserted inside the poisson solver 
         call PSolverNC(atoms%geocode,'D',iproc,nproc,n1i,n2i,n3i,n3d,ixc,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,nspin)
     else
              
        call PSolver(atoms%geocode,'D',iproc,nproc,n1i,n2i,n3i,ixc,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,nspin)
        
     end if

     call HamiltonianApplication(iproc,nproc,atoms,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
          norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          wfd,bounds,nlpspd,proj,ngatherarr,n1i*n2i*n3p,&
          rhopot(1,1,1+i3xcsh,1),psi,hpsi,ekin_sum,epot_sum,&
          eproj_sum,nspin,nspinor,spinsgn,hybrid_on)

     energybs=ekin_sum+epot_sum+eproj_sum
     energy_old=energy
     energy=energybs-ehart+eexcu-vexcu+eion

     !check for convergence or whether max. numb. of iterations exceeded
     if (gnrm <= gnrm_cv .or. iter == itermax) then 
        if (iproc.eq.0) then 
           write( *,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write( *,'(1x,a)') &
                '--------------------------------------------------- End of Wavefunction Optimisation'
           write( *,'(1x,a,3(1x,1pe18.11))') &
                'final  ekin,  epot,  eproj ',ekin_sum,epot_sum,eproj_sum
           write( *,'(1x,a,3(1x,1pe18.11))') &
                'final ehart, eexcu,  vexcu ',ehart,eexcu,vexcu
           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') &
                'FINAL iter,total energy,gnrm',iter,energy,gnrm
           !write(61,*)hgrid,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
           if (energy > energy_min) write( *,'(1x,a,1pe9.2)')&
                'WARNING: Found an energy value lower than the FINAL energy, delta:',energy-energy_min
        end if
        if (gnrm <= gnrm_cv) infocode=0
        exit wfn_loop 
     endif

     call hpsitopsi(atoms%geocode,iproc,nproc,norb,norbp,occup,hx,hy,hz,n1,n2,n3,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,bounds%kb,&
          eval,ncong,iter,idsx,idsx_actual,ads,energy,energy_old,energy_min,&
          alpha,gnrm,scprsum,psi,psit,hpsi,psidst,hpsidst,nspin,nspinor,spinsgn,hybrid_on)

     tt=(energybs-scprsum)/scprsum
     if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
          (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
        write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
     endif
     if (iproc.eq.0) then
        write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
             ekin_sum,epot_sum,eproj_sum
        write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
        write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
     endif

     if (in%inputPsiId == 0) then
        if ((gnrm > 4.d0 .and. norbu/=norbd) .or. (norbu==norbd .and. gnrm > 10.d0)) then
           if (iproc == 0) then
              write( *,'(1x,a)')&
                   'Error: the norm of the residue is too large also with input wavefunctions.'
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

  ! transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
  call last_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit,&
       occup,evsum,eval)

  if (abs(evsum-energybs) > 1.d-8 .and. iproc==0) write( *,'(1x,a,2(1x,1pe20.13))')&
       'Difference:evsum,energybs',evsum,energybs
 
  if (nvirt > 0 .and. in%inputPsiId == 0) then
     call davidson(iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i,atoms,&
          cpmult,fpmult,radii_cf,&
          norb,norbu,norbp,nvirte,nvirtep,nvirt,gnrm_cv,nplot,n1,n2,n3,nvctrp,&
          hx,hy,hz,rxyz,rhopot,occup,i3xcsh,n3p,itermax,wfd,bounds,nlpspd,proj,  &
          pkernel,ixc,psi,psivirt,eval,ncong,nscatterarr,ngatherarr,hybrid_on)
  end if
  
  !project the wavefunctions on a gaussian basis and keep in memory
  if (in%gaussian_help) then
     if (iproc.eq.0) then
        write( *,'(1x,a)')&
             '---------------------------------------------------------- Gaussian Basis Projection'
     end if

     !extract the gaussian basis from the pseudowavefunctions
!!$     if (in%inputPsiId == 11) then
!!$        !extract the gaussian basis from the pseudowavefunctions
!!$        call gaussian_pswf_basis(iproc,atoms,rxyz,gbd)
!!$     else if (in%inputPsiId == 12) then
!!$        !extract the gaussian basis from the pseudopotential
!!$        call gaussian_psp_basis(atoms,rxyz,gbd)
!!$     end if

     !extract the gaussian basis from the pseudowavefunctions
     call gaussian_pswf_basis(iproc,atoms,rxyz,gbd)

     if (.not. associated(gaucoeffs)) then
        allocate(gaucoeffs(gbd%ncoeff,norbp),stat=i_stat)
        call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)
     end if

     allocate(thetaphi(2,gbd%nat),stat=i_stat)
     call memocc(i_stat,thetaphi,'thetaphi',subname)
     thetaphi=0.0_gp

     call wavelets_to_gaussians(atoms%geocode,norbp,n1,n2,n3,gbd,thetaphi,hx,hy,hz,wfd,psi,gaucoeffs)

     i_all=-product(shape(thetaphi))*kind(thetaphi)
     deallocate(thetaphi,stat=i_stat)
     call memocc(i_stat,i_all,'thetaphi',subname)

  end if


  !  write all the wavefunctions into files
  if (in%output_wf) then
     !add flag for writing waves in the gaussian basis form
     if (in%gaussian_help) then

!!$        call gaussian_orthogonality(iproc,nproc,norb,norbp,gbd,gaucoeffs)
!!$
!!$        call gaussian_orthogonality(iproc,nproc,norb,norbp,gbd,gaucoeffs)
        !write the coefficients and the basis on a file
        call write_gaussian_information(iproc,nproc,norb,norbp,gbd,gaucoeffs,eval,&
             'wavefunctions.gau')

        !build dual coefficients
        call dual_gaussian_coefficients(norbp,gbd,gaucoeffs)
        !control the accuracy of the expansion
        call check_gaussian_expansion(atoms%geocode,iproc,nproc,norb,norbp,&
             n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,hy,hz,wfd,psi,gbd,gaucoeffs)

        call deallocate_gwf(gbd,subname)
        i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
        deallocate(gaucoeffs,stat=i_stat)
        call memocc(i_stat,i_all,'gaucoeffs',subname)
        nullify(gbd%rxyz)

     else
        call  writemywaves(iproc,norb,norbp,n1,n2,n3,hx,hy,hz,atoms%nat,rxyz,wfd,psi,eval)
        write( *,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
     end if
  end if


  !------------------------------------------------------------------------
  ! here we start the calculation of the forces
  if (iproc.eq.0) then
     write( *,'(1x,a)')&
          '----------------------------------------------------------------- Forces Calculation'
  end if

  ! Selfconsistent potential is saved in rhopot, 
  ! new arrays rho,pot for calculation of forces ground state electronic density

  allocate(spinsgn_foo(norb+ndebug),stat=i_stat)
  call memocc(i_stat,spinsgn_foo,'spinsgn_foo',subname)
  spinsgn_foo(:)=1.0d0
  ! Potential from electronic charge density

  !manipulate scatter array for avoiding the GGA shift
  do jproc=0,nproc-1
     !n3d=n3p
     nscatterarr(jproc,1)=nscatterarr(jproc,2)
     !i3xcsh=0
     nscatterarr(jproc,4)=0
  end do

  !here there are the spinor which must be taken into account
  if(nproc>1) then
     if (n3p>0) then
        allocate(rho(n1i*n2i*n3p+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     else
        allocate(rho(1+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     end if
  else
     if (n3p>0) then
        allocate(rho(n1i*n2i*n3p*nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     else
        allocate(rho(1+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     end if
  end if

  call sumrho(atoms%geocode,iproc,nproc,norb,norbp,0,n1,n2,n3,hxh,hyh,hzh,occup,  & 
       wfd,psi,rho,n1i*n2i*n3p,nscatterarr,1,nspinor,spinsgn_foo,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,hybrid_on)

  i_all=-product(shape(spinsgn_foo))*kind(spinsgn_foo)
  deallocate(spinsgn_foo,stat=i_stat)
  call memocc(i_stat,i_all,'spinsgn_foo',subname)

  !plot also the ionic potential
  if (in%output_grid) then
     call plot_density(atoms%geocode,'pot_ion.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,1,&
          atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,pot_ion)
  end if

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)

  !plot the density on the density.pot file
  if (in%output_grid) then
     call plot_density(atoms%geocode,'density.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nelec,&
     atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rho)
  end if


  if (n3p>0) then
     allocate(pot(n1i,n2i,n3p,1+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)
  else
     allocate(pot(1,1,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)
  end if

  !calculate electrostatic potential
  call DCOPY(n1i*n2i*n3p,rho,1,pot,1) 
  call PSolver(atoms%geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
       pot,pkernel,pot,ehart_fake,eexcu_fake,vexcu_fake,0.d0,.false.,1)
  !here nspin=1 since ixc=0

  !plot also the electrostatic potential
  if (in%output_grid) then
     call plot_density(atoms%geocode,'potential.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,1,&
          atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,pot)
  end if

  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'pkernel',subname)

  allocate(gxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,gxyz,'gxyz',subname)

  call timing(iproc,'Forces        ','ON')
  ! calculate local part of the forces gxyz
  call local_forces(iproc,nproc,atoms,rxyz,hxh,hyh,hzh,&
       n1,n2,n3,n3p,i3s+i3xcsh,n1i,n2i,n3i,rho,pot,gxyz)

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)
  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)

  if (iproc == 0) write( *,'(1x,a)',advance='no')'Calculate projectors derivatives...'

  if (iproc == 0) write( *,'(1x,a)',advance='no')'done, calculate nonlocal forces...'

  call nonlocal_forces(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,atoms,rxyz,radii_cf,&
     norb,norbp,nspinor,occup,nlpspd,proj,wfd,psi,gxyz,calc_tail) !refill projectors for tails

  if (iproc == 0) write( *,'(1x,a)')'done.'

  ! Add up all the force contributions
  if (nproc > 1) then
     call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     do iat=1,atoms%nat
        fxyz(1,iat)=gxyz(1,iat)
        fxyz(2,iat)=gxyz(2,iat)
        fxyz(3,iat)=gxyz(3,iat)
     enddo
  end if

  !add to the forces the ionic contribution 
  do iat=1,atoms%nat
     fxyz(1,iat)=fxyz(1,iat)+fion(1,iat)
     fxyz(2,iat)=fxyz(2,iat)+fion(2,iat)
     fxyz(3,iat)=fxyz(3,iat)+fion(3,iat)
  enddo

  i_all=-product(shape(fion))*kind(fion)
  deallocate(fion,stat=i_stat)
  call memocc(i_stat,i_all,'fion',subname)
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

  !------------------------------------------------------------------------
  if (calc_tail .and. atoms%geocode == 'F') then
     call timing(iproc,'Tail          ','ON')
     !    Calculate energy correction due to finite size effects
     !    ---reformat potential
     allocate(pot(n1i,n2i,n3i,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)

     if (nproc > 1) then
        call MPI_ALLGATHERV(rhopot(1,1,1+i3xcsh,1),n1i*n2i*n3p,&
             MPI_DOUBLE_PRECISION,pot(1,1,1,1),ngatherarr(0,1),ngatherarr(0,2), & 
             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        !print '(a,2f12.6)','RHOup',sum(abs(rhopot(:,:,:,1))),sum(abs(pot(:,:,:,1)))
        if(nspin==2) then
           !print '(a,2f12.6)','RHOdw',sum(abs(rhopot(:,:,:,2))),sum(abs(pot(:,:,:,2)))
           if (n3d /= n3p) then
              i03=1+i3xcsh+n3p
              i04=1
           else
              i03=1
              i04=2
           end if
           call MPI_ALLGATHERV(rhopot(1,1,i03,i04),n1i*n2i*n3p,&
                MPI_DOUBLE_PRECISION,pot(1,1,1,2),ngatherarr(0,1),ngatherarr(0,2), & 
                MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        end if
     else
        do ispin=1,nspin
           !here one could have not allocated pot and: call move_alloc(rhopot,pot) 
           !(but it is a Fortran 95/2003 spec)
           do i3=1,n3i
              do i2=1,n2i
                 do i1=1,n1i
                    pot(i1,i2,i3,ispin)=rhopot(i1,i2,i3,ispin)
                 enddo
              enddo
           enddo
        end do
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

     !pass hx instead of hgrid since we are only in free BC
     call CalculateTailCorrection(iproc,nproc,atoms,n1,n2,n3,rbuf,norb,norbp,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,nlpspd,ncongt,eval,&
          pot,hx,rxyz,radii_cf,crmult,frmult,cpmult,fpmult,nspin,spinsgn,&
          proj,psi,occup,in%output_grid,ekin_sum,epot_sum,eproj_sum)

     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot',subname)

     !if (iproc==0) then
     !   open(61)
     !   write(61,'(4(f9.3),1x,7(1pe19.11))',advance='no')&
     !        hgrid,alat1,alat2,alat3,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
     !end if

     energybs=ekin_sum+epot_sum+eproj_sum
     energy=energybs-ehart+eexcu-vexcu+eion

     !if (iproc==0) then
     !   write(61,'(1pe19.11)')energy
     !   close(61)
     !end if

     if (iproc.eq.0) then
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

       i_all=-product(shape(pot_ion))*kind(pot_ion)
       deallocate(pot_ion,stat=i_stat)
       call memocc(i_stat,i_all,'pot_ion',subname)

       i_all=-product(shape(pkernel))*kind(pkernel)
       deallocate(pkernel,stat=i_stat)
       call memocc(i_stat,i_all,'pkernel',subname)

       ! calc_tail false
       i_all=-product(shape(rhopot))*kind(rhopot)
       deallocate(rhopot,stat=i_stat)
       call memocc(i_stat,i_all,'rhopot',subname)
       i_all=-product(shape(nscatterarr))*kind(nscatterarr)
       deallocate(nscatterarr,stat=i_stat)
       call memocc(i_stat,i_all,'nscatterarr',subname)
       i_all=-product(shape(ngatherarr))*kind(ngatherarr)
       deallocate(ngatherarr,stat=i_stat)
       call memocc(i_stat,i_all,'ngatherarr',subname)

       i_all=-product(shape(fion))*kind(fion)
       deallocate(fion,stat=i_stat)
       call memocc(i_stat,i_all,'fion',subname)

    end if
    !deallocate wavefunction for virtual orbitals
    if (in%nvirt > 0) then
       i_all=-product(shape(psivirt))*kind(psivirt)
       deallocate(psivirt)
       call memocc(i_stat,i_all,'psivirt','cluster')
    end if

    if (atoms%geocode == 'F') then
       call deallocate_bounds(bounds,'cluster')
    end if
    !****************Added by Alexey***********************************************************	
    if (atoms%geocode == 'P' .and. hybrid_on) then 

       i_all=-product(shape(bounds%kb%ibxy_f))*kind(bounds%kb%ibxy_f)
       deallocate(bounds%kb%ibxy_f,stat=i_stat)
       call memocc(i_stat,i_all,'bounds%kb%ibxy_f',subname)

       i_all=-product(shape(bounds%kb%ibxz_f))*kind(bounds%kb%ibxz_f)
       deallocate(bounds%kb%ibxz_f,stat=i_stat)
       call memocc(i_stat,i_all,'bounds%kb%ibxz_f',subname)

       i_all=-product(shape(bounds%kb%ibyz_f))*kind(bounds%kb%ibyz_f)
       deallocate(bounds%kb%ibyz_f,stat=i_stat)
       call memocc(i_stat,i_all,'bounds%kb%ibyz_f',subname)

       i_all=-product(shape(bounds%sb%ibxy_ff))*kind(bounds%sb%ibxy_ff)
       deallocate(bounds%sb%ibxy_ff,stat=i_stat)
       call memocc(i_stat,i_all,'ibxy_ff',subname)
       i_all=-product(shape(bounds%sb%ibzzx_f))*kind(bounds%sb%ibzzx_f)
       deallocate(bounds%sb%ibzzx_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibzzx_f',subname)
       i_all=-product(shape(bounds%sb%ibyyzz_f))*kind(bounds%sb%ibyyzz_f)
       deallocate(bounds%sb%ibyyzz_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibyyzz_f',subname)

       i_all=-product(shape(bounds%gb%ibyz_ff))*kind(bounds%gb%ibyz_ff)
       deallocate(bounds%gb%ibyz_ff,stat=i_stat)
       call memocc(i_stat,i_all,'ibyz_ff',subname)

       i_all=-product(shape(bounds%gb%ibzxx_f))*kind(bounds%gb%ibzxx_f)
       deallocate(bounds%gb%ibzxx_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibzxx_f',subname)

       i_all=-product(shape(bounds%gb%ibxxyy_f))*kind(bounds%gb%ibxxyy_f)
       deallocate(bounds%gb%ibxxyy_f,stat=i_stat)
       call memocc(i_stat,i_all,'ibxxyy_f',subname)
    endif

    !***************************************************************************************	

    !semicores useful only for the input guess
    i_all=-product(shape(atoms%iasctype))*kind(atoms%iasctype)
    deallocate(atoms%iasctype,stat=i_stat)
    call memocc(i_stat,i_all,'iasctype',subname)
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

    i_all=-product(shape(occup))*kind(occup)
    deallocate(occup,stat=i_stat)
    call memocc(i_stat,i_all,'occup',subname)
    i_all=-product(shape(spinsgn))*kind(spinsgn)
    deallocate(spinsgn,stat=i_stat)
    call memocc(i_stat,i_all,'spinsgn',subname)
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

    !end of wavefunction minimisation
    call timing(iproc,'LAST','PR')
    call timing(iproc,'              ','RE')
    call cpu_time(tcpu1)
    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)
    if (iproc == 0) &
         write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time for root process ', iproc,tel,tcpu1-tcpu0

  end subroutine deallocate_before_exiting

END SUBROUTINE cluster
!!***

