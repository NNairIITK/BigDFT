!> @file
!!  Sandbox for test on linear version
!!
!! @author
!!    Copyright (C) 2011-2012 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
program sandbox
!!  use BigDFT_API
!!  use Poisson_Solver
  implicit none
!!  character(len=*), parameter :: subname='sandbox'
!!  logical :: dokernel=.false.
!!  logical :: endloop
!!  integer :: n1i,n2i,n3i,iproc,nproc,ifine,i_stat,i_all,nelec
!!  integer :: n3d,n3p,n3pi,i3xcsh,i3s,n1,n2,n3,ndegree_ip
!!  integer :: idsx_actual,ndiis_sd_sw,idsx_actual_before,iter,iorb,jorb
!!  real(gp) :: hxh,hyh,hzh
!!  real(gp) :: tt,gnrm,gnrm_zero,epot_sum,eexctX,ekin_sum,eproj_sum,alpha
!!  real(gp) :: energy,energy_min,energy_old,energybs,evsum,scprsum
!!  type(atoms_data) :: atoms
!!  type(input_variables) :: in
!!  type(orbitals_data) :: orbs
!!  type(locreg_descriptors) :: Glr
!!  type(nonlocal_psp_descriptors) :: nlpspd
!!  type(comms_cubic) :: comms
!!  type(GPU_pointers) :: GPU
!!  type(diis_objects) :: diis
!!  character(len=4) :: itername
!!  real(gp), dimension(3) :: shift
!!  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
!!  real(gp), dimension(:,:), allocatable :: radii_cf
!!  real(wp), dimension(:), pointer :: hpsi,psit,psi,psidst,hpsidst,proj
!!  real(dp), dimension(:), pointer :: pkernel,pot_ion
!!  real(gp), dimension(:,:), pointer :: rxyz
!!  ! arrays for DIIS convergence accelerator
!!  real(wp), dimension(:,:,:), pointer :: ads
!!  !##########################################
!!  ! My variables
!!  !##########################################
!!  integer :: ilr
!!  integer, parameter :: nlr=2,alr=1,blr=2
!!  integer :: ldim   ! dimension of lpsi
!!  integer :: fmin   ! min(1,nseg_f)
!!  integer :: isovrlp
!!  real(wp) :: scpr
!!  integer,dimension(:),allocatable :: projflg
!!  integer,dimension(3,nlr) :: outofzone
!!  real(gp), dimension(3,2) :: cxyz
!!!  real(gp), dimension(nlr),parameter :: locrad=(/16.0, 16.0 /)
!!  real(wp),dimension(nlr),parameter :: locrad=(/ 3.0, 3.0  /)
!!!  real(wp),dimension(nlr),parameter :: locrad=(/ 2.565213453, 2.565213453  /) 
!!  real(wp), dimension(:),allocatable :: lpsi    ! local projection of |Psi>
!!  real(wp), dimension(:),allocatable :: lhpsi   ! local projection of H|Psi>
!!  real(wp), dimension(:),allocatable :: lppsi   ! local projection of H|Psi>
!!  real(wp), dimension(:),allocatable :: ppsi   ! local projection of H|Psi>
!!  real(wp), dimension(:), pointer :: Lproj
!!  real(wp),dimension(:,:),allocatable :: overlap1
!!  real(wp),dimension(:,:),allocatable :: koverlap
!!  real(wp),dimension(:,:),allocatable :: poverlap
!!  type(locreg_descriptors), dimension(nlr) :: Llr 
!!  type(locreg_descriptors), allocatable :: Olr(:)
!!  real(wp), dimension(:), pointer :: potential
!!  real :: sum_pot !debug
!!  integer :: ii    ! debug  
!!  type(nonlocal_psp_descriptors) :: Lnlpspd
!!  type(atoms_data) :: Latoms
!!  real(gp) :: eproj
!!  !for the moment no need to have parallelism
!!  iproc=0
!!  nproc=1
!!
!!  !initalise the variables for the calculation
!!  !standard names
!!  call standard_inputfile_names(in)
!!  call read_input_variables(iproc,'posinp',in, atoms, rxyz)
!!
!!  if (iproc == 0) then
!!     call print_general_parameters(in,atoms)
!!  end if
!!       
!!  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
!!  call memocc(i_stat,radii_cf,'radii_cf',subname)
!!
!!  call system_properties(iproc,nproc,in,atoms,orbs,radii_cf,nelec)
!!
!!  ! Determine size alat of overall simulation cell and shift atom positions
!!  ! then calculate the size in units of the grid space
!!  call system_size(iproc,atoms,rxyz,radii_cf,in%crmult,in%frmult,in%hx,in%hy,in%hz,&
!!       Glr,shift)
!!
!!  !variables substitution for the PSolver part
!!  hxh=0.5d0*in%hx
!!  hyh=0.5d0*in%hy
!!  hzh=0.5d0*in%hz
!!  n1i=Glr%d%n1i
!!  n2i=Glr%d%n2i
!!  n3i=Glr%d%n3i
!!
!!  n1=Glr%d%n1
!!  n2=Glr%d%n2
!!  n3=Glr%d%n3
!!
!!  call timing(iproc,'CrtDescriptors','ON')
!!  call createWavefunctionsDescriptors(iproc,in%hx,in%hy,in%hz,&
!!       atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr)
!!  call timing(iproc,'CrtDescriptors','OF')
!!
!!  ! Calculate all projectors, or allocate array for on-the-fly calculation
!!  call timing(iproc,'CrtProjectors ','ON')
!!  call createProjectorsArrays(iproc,Glr,rxyz,atoms,orbs,&
!!       radii_cf,in%frmult,in%frmult,in%hx,in%hy,in%hz,nlpspd,proj)
!!  call timing(iproc,'CrtProjectors ','OF')
!!
!!  !allocate communications arrays
!!  call orbitals_communicators(iproc,nproc,Glr,orbs,comms)  
!!
!!  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
!!  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
!!  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
!!  call memocc(i_stat,ngatherarr,'ngatherarr',subname)
!!
!!  call createDensPotDescriptors(iproc,nproc,atoms%geocode,'D',&
!!       Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,&
!!       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)
!!
!!  !commented out, to be used in the future
!!  if (dokernel) then
!!     ndegree_ip=16 !default value 
!!     call createKernel(iproc,nproc,atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,&
!!          pkernel)
!!  else
!!     nullify(pkernel)
!!  end if
!!
!!  !allocate ionic potential
!!  if (n3pi > 0) then
!!     allocate(pot_ion(n1i*n2i*n3pi+ndebug),stat=i_stat)
!!     call memocc(i_stat,pot_ion,'pot_ion',subname)
!!  else
!!     allocate(pot_ion(1+ndebug),stat=i_stat)
!!     call memocc(i_stat,pot_ion,'pot_ion',subname)
!!  end if
!!
!!  allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
!!  call memocc(i_stat,psi,'psi',subname)
!!  allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
!!  call memocc(i_stat,hpsi,'hpsi',subname)
!!  if (nproc == 1) then
!!     nullify(psit)
!!  else
!!     stop 'for the moment only sequential runs'
!!  end if
!!
!!  ! allocate arrays necessary for DIIS convergence acceleration
!!  !the allocation with npsidim is not necessary here since DIIS arrays
!!  !are always calculated in the transposed form
!!  if (in%idsx > 0) then
!!     call allocate_diis_objects(in%idsx,sum(comms%ncntt(0:nproc-1)),orbs%nkptsp,diis)  
!!  endif
!!
!!  sum_pot = 0.0d0
!!  do ii=1,n1i*n2i*n3pi
!!     sum_pot = sum_pot + pot_ion(ii)
!!  end do
!!  print *,'Sum of the potential :',sum_pot
!!
!!  !create input guess wavefunction
!!  call psi_from_gaussians(iproc,nproc,atoms,orbs,Glr,rxyz,in%hx,in%hy,in%hz,in%nspin,psi)
!!
!!  !calculate scalar product between wavefunctions
!!  !this scheme works only in sequential
!!  write(*,'(A26)') 'Overlap matrix with ddot:'
!!  do iorb=1,orbs%norb
!!     do jorb=iorb,orbs%norb
!!        write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
!!             dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
!!             psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
!!     end do
!!  end do
!!
!!
!!! ########################################
!!!  Added by me to have kinetic part 
!!! ########################################
!!! Fake Eigenvalues for orbitals
!!     orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0
!!! #################
!!!  Kinetic part
!!! #################
!!
!!  !allocate the potential in the full box
!!     call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3pi,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,in%nspin,&
!!          orbs%norb,orbs%norbp,ngatherarr,pot_ion,potential)
!!
!!  !put in hpsi the kinetic operator (pot_ion = 0d.0)
!!     call HamiltonianApplication(iproc,nproc,atoms,orbs,in%hx,in%hy,in%hz,rxyz,&
!!         nlpspd,proj,Glr,ngatherarr,potential,psi,hpsi,ekin_sum,epot_sum,eexctX,&
!!         eproj_sum,in%nspin,GPU)
!!
!!  !calculate scalar product between wavefunctions
!!  !this scheme works only in sequential
!!  !   write(*,*)' '
!!  !   write(*,'(A34)') 'Kinetic overlap matrix with ddot:'
!!  !   do iorb=1,orbs%norb
!!  !      do jorb=iorb,orbs%norb
!!  !         write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
!!  !              dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
!!  !              hpsi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
!!  !      end do
!!  !   end do
!! 
!!!#######################################################################################################################
!!! Test the overlap calculation of wavefunctions defined in different localisation regions
!!!########################################################################################################################
!! write(*,'(A32)') 'Entering the localisation code.'
!!
!!  cxyz(:,1) = rxyz(:,1)
!!  cxyz(:,2) = rxyz(:,size(rxyz,2))
!!
!!! Write some physical information on the Glr
!!  write(*,'(a24,3i4)')'Global region n1,n2,n3:',Glr%d%n1,Glr%d%n2,Glr%d%n3
!!  write(*,'(a27,f6.2,f6.2,f6.2)')'Global dimension (x,y,z):',Glr%d%n1*in%hx,Glr%d%n2*in%hy,Glr%d%n3*in%hz
!!  write(*,'(a17,f12.2)')'Global volume: ',Glr%d%n1*in%hx*Glr%d%n2*in%hy*Glr%d%n3*in%hz
!!  print *,'Global statistics:',Glr%wfd%nseg_c,Glr%wfd%nseg_f,Glr%wfd%nvctr_c,Glr%wfd%nvctr_f
!!   
!!! First, determine the localisation regions
!!  call determine_locreg_periodic(iproc,nlr,rxyz,locrad,in%hx,in%hy,in%hz,Glr,Llr)
!!
!!  print *,'Outside determine_locreg2'
!!! Plot the localization regions
!!!  call draw_locregs(2,in%hx,in%hy,in%hz,Llr)
!!
!!! Second, calculate the number of overlap regions using periodicity
!!  call get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr)
!!
!!  write(*,*)'Outside get_number_of_overlap_region:',isovrlp
!!
!!  if(isovrlp > 0) then
!!
!!     ! allocate the overlap region descriptors (Olr)
!!     print *,'isovrlp:',isovrlp
!!     allocate(Olr(isovrlp),stat=i_stat)
!!!     call memocc(i_stat,Olr,'Olr',subname)
!!   
!!     ! Third, construct the overlap region descriptors
!!     call get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr)
!!  
!!!    Write some physical information on the overlap
!!     do ilr=1,isovrlp
!!        write(*,'(a25,3i4)')'Overlap region',ilr
!!        write(*,*)' n1,n2,n3:',Olr(ilr)%d%n1,Olr(ilr)%d%n2,Olr(ilr)%d%n3
!!        write(*,'(a27,f6.2,f6.2,f6.2)')'Overlap dimension (x,y,z):',Olr(ilr)%d%n1*in%hx,Olr(ilr)%d%n2*in%hy,Olr(ilr)%d%n3*in%hz
!!        write(*,'(a17,f12.2)')'Overlap volume: ',Olr(ilr)%d%n1*in%hx*Olr(ilr)%d%n2*in%hy*Olr(ilr)%d%n3*in%hz
!!        write(*,*)'Overlap starting point:',Olr(ilr)%ns1,Olr(ilr)%ns2,Olr(ilr)%ns3
!!     end do
!!
!!!    Plot an overlap region
!!     call draw_locregs(1,in%hx,in%hy,in%hz,Olr(1))
!!
!!! Allocate overlaps accumulators
!!     allocate(overlap1(orbs%norb,orbs%norb),stat=i_stat)
!!     call memocc(i_stat,overlap1,'overlap1',subname)
!!     call to_zero(orbs%norb*orbs%norb,overlap1)
!!
!!     allocate(koverlap(orbs%norb,orbs%norb),stat=i_stat)
!!     call memocc(i_stat,koverlap,'koverlap',subname)
!!     call to_zero(orbs%norb*orbs%norb,koverlap)
!!
!!     allocate(poverlap(orbs%norb,orbs%norb),stat=i_stat)
!!     call memocc(i_stat,poverlap,'poverlap',subname)
!!     call to_zero(orbs%norb*orbs%norb,poverlap)    
!!
!!! Third, transform the wavefunction to overlap regions
!!     do ilr=1,isovrlp
!!       print *,'Treating overlap region (ilr):',ilr
!!       ldim = (Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*orbs%norb*orbs%nspinor
!!
!!       ! Allocate the local wavefunction (in one overlap region)
!!       allocate(lpsi(ldim+ndebug), stat=i_stat)
!!       call memocc(i_stat,lpsi,'lpsi',subname)
!! 
!!       ! Project the wavefunction inside the overlap region
!!       call psi_to_locreg(Glr,ilr,ldim,Olr,lpsi,isovrlp,orbs,psi)
!!
!!! Calculate overlap in localization region using ddot
!!!     write(*,'(A26)') 'Overlap matrix with ddot:'
!!!     do iorb=1,orbs%norb
!!!        do jorb=iorb,orbs%norb
!!!           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
!!!                dot(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f,lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
!!!                lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
!!!        end do
!!!     end do
!!
!! !Calculate the overlap
!!     do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!         overlap1(iorb,jorb) = overlap1(iorb,jorb) + dot(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f,lpsi((Olr(ilr)%wfd%nvctr_c&
!!                  +7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
!!        end do
!!     end do
!!
!!! Calculate overlap in localization region using wpdot
!!!     write(*,'(A27)')'Overlap matrix with wpdot:'
!!!     fmin = min(Olr(ilr)%wfd%nseg_f,1)  ! checks if there is some fine_grid in the region, if not, do not shift keyv
!!!     do iorb=1,orbs%norb
!!!        do jorb=iorb,orbs%norb
!!!           call wpdot(Olr(ilr)%wfd%nvctr_c,Olr(ilr)%wfd%nvctr_f,Olr(ilr)%wfd%nseg_c,Olr(ilr)%wfd%nseg_f,&
!!!&                  Olr(ilr)%wfd%keyv(1),Olr(ilr)%wfd%keyv(Olr(ilr)%wfd%nseg_c+fmin),Olr(ilr)%wfd%keyg(1,1),&
!!!&                  Olr(ilr)%wfd%keyg(1,Olr(ilr)%wfd%nseg_c+fmin),lpsi(1+(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)),&
!!!&                  lpsi(Olr(ilr)%wfd%nvctr_c+fmin+(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)),&
!!!&                  Olr(ilr)%wfd%nvctr_c,Olr(ilr)%wfd%nvctr_f,Olr(ilr)%wfd%nseg_c,Olr(ilr)%wfd%nseg_f,&
!!!&                  Olr(ilr)%wfd%keyv(1),Olr(ilr)%wfd%keyv(Olr(ilr)%wfd%nseg_c+fmin),Olr(ilr)%wfd%keyg(1,1),&
!!!&                  Olr(ilr)%wfd%keyg(1,Olr(ilr)%wfd%nseg_c+fmin),lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),&
!!!&                  lpsi(Olr(ilr)%wfd%nvctr_c+(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+fmin),scpr)
!!!           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,scpr
!!!        end do
!!!     end do
!!
!!! Plot the orbitals that are inside the overlap region
!!! plot first orbital
!!!     call plot_wf_sandbox('orb1_ovrlp',1,atoms,Olr(ilr),hxh,hyh,hzh,rxyz,lpsi,'')
!!! plot second orbital
!!!     call plot_wf_sandbox('orb2_ovrlp',1,atoms,Olr(ilr),hxh,hyh,hzh,rxyz,&
!!!&                     lpsi(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f+1),'')
!!
!!! Plot first orbital
!!     call plot_wf_sandbox('iter0-1',1,atoms,Glr,hxh,hyh,hzh,rxyz,psi,'')
!!! Plot second orbital
!!     call plot_wf_sandbox('iter0-2',1,atoms,Glr,hxh,hyh,hzh,rxyz,psi(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f+1),'')
!!
!!     allocate(orbs%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
!!     call memocc(i_stat,orbs%eval,'eval',subname)
!!
!!! ########################
!!!  Kinetic Part
!!! ########################
!!! Now transform hpsi to overlap region
!!     allocate(lhpsi(ldim+ndebug), stat=i_stat)
!!     call memocc(i_stat,lhpsi,'lhpsi',subname)
!!
!!!     call psi_to_locreg(Glr,ilr,ldim,Olr,lhpsi,isovrlp,orbs,hpsi)
!!
!!!     write(*,'(A51)') 'Kinetic overlap matrix with ddot and localization:'
!!     do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!           koverlap(iorb,jorb) = koverlap(iorb,jorb) + dot(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f,&
!!           lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
!!           lhpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
!!        end do
!!     end do
!!
!!! ##################
!!! LOCAL: Non-local Projectors part
!!! ##################
!!
!!! allocate  projflg
!!  allocate(projflg(atoms%nat+ndebug),stat=i_stat)
!!  call memocc(i_stat,projflg,'projflg',subname)
!!
!!! Need to restrict non-local projectors to overlap region
!!! First make the descriptors
!!  call nlpspd_to_locreg(in,iproc,Glr,Olr(ilr),rxyz,atoms,orbs,&
!!&       radii_cf,in%frmult,in%frmult,in%hx,in%hy,in%hz,nlpspd,Lnlpspd,projflg)
!!
!!! Allocate Lproj and Lppsi
!!  allocate(Lproj(Lnlpspd%nprojel+ndebug),stat=i_stat)
!!  call memocc(i_stat,Lproj,'Lproj',subname)
!!  allocate(lppsi(ldim+ndebug), stat=i_stat)
!!  call memocc(i_stat,lppsi,'lppsi',subname)
!!  allocate(ppsi((Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f)*orbs%norb*orbs%nspinor+ndebug), stat=i_stat)
!!  call memocc(i_stat,ppsi,'ppsi',subname)
!!
!!! initialise lppsi and ppsi
!!  call to_zero((Olr(ilr)%wfd%nvctr_c + 7*Olr(ilr)%wfd%nvctr_f)*orbs%norb*orbs%nspinor,lppsi)
!!  call to_zero((Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f)*orbs%norb*orbs%nspinor,ppsi)
!!
!!! Apply without the localization (for comparison)
!!  eproj_sum = 0.0_gp
!!  call applyprojectorsonthefly(iproc,orbs,atoms,Glr,&
!!          rxyz,in%hx,in%hy,in%hz,Glr%wfd,nlpspd,proj,psi,ppsi,eproj_sum)
!!
!!! Calculate dotprod: <Psi_a|p><p|Psi_b>
!!   if(ilr == 1) then
!!     write(*,'(A37)') 'NL-Operator overlap matrix with ddot:'
!!     do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
!!               dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
!!               ppsi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
!!        end do
!!     end do
!!   end if
!!
!!
!!! Fill and Apply the projectors on the wavefunctions
!!  call apply_local_projectors(ilr,in%nspin,atoms,in%hx,in%hy,in%hz,Olr(ilr),Lnlpspd,orbs,projflg,lpsi,rxyz,lppsi,eproj)
!!
!!! Calculate dotprod: <Psi_a|p><p|Psi_b>
!!! write(*,'(A54)') 'NL-Operator overlap matrix with ddot and localization:'
!!     do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!           poverlap(iorb,jorb) = poverlap(iorb,jorb) + dot(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f,&
!!               lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
!!               lppsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
!!        end do
!!     end do
!!
!!! Deallocations
!!     deallocate(projflg,stat=i_stat)
!!     call memocc(i_stat,i_all,'projflg',subname) 
!!
!!     deallocate(lpsi,stat=i_stat)
!!     call memocc(i_stat,i_all,'lpsi',subname)
!!
!!     deallocate(lhpsi,stat=i_stat)
!!     call memocc(i_stat,i_all,'lhpsi',subname)
!!
!!     deallocate(lppsi,stat=i_stat)
!!     call memocc(i_stat,i_all,'lppsi',subname)
!!
!!     deallocate(ppsi,stat=i_stat)
!!     call memocc(i_stat,i_all,'ppsi',subname)
!!
!!   end do ! ilr
!!
!!! Write total overlap of wavefunction for localisation regions
!!     write(*,'(A26)') 'Overlap matrix with ddot:'
!!     do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,overlap1(iorb,jorb)
!!        end do
!!     end do
!!
!!     write(*,'(A34)') 'Kinetic overlap matrix with ddot:'
!!     do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
!!                dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
!!                hpsi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
!!        end do
!!     end do
!!
!!     write(*,'(A51)') 'Kinetic overlap matrix with ddot and localization:'
!!     do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,koverlap(iorb,jorb)
!!        end do
!!     end do
!!
!!     write(*,'(A54)') 'NL-Operator overlap matrix with ddot and localization:'
!!       do iorb=1,orbs%norb
!!        do jorb=iorb,orbs%norb
!!           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,poverlap(iorb,jorb)
!!        end do
!!     end do
!!
!!     deallocate(overlap1,stat=i_stat)
!!     call memocc(i_stat,orbs%norb*orbs%norb,'overlap1',subname)
!!
!!     deallocate(koverlap,stat=i_stat)
!!     call memocc(i_stat,orbs%norb*orbs%norb,'koverlap',subname)
!!  end if   !isovrlp
!!
!!  stop  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------------------------------------------
!!!#######################################################################################################################
!!! End of Localization part
!!!######################################################################################################################
!!
!!
!!
!!  !othogonalise them
!!  !transpose the psi wavefunction
!!  call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,&
!!       psi,work=hpsi)
!!  call orthogonalize(iproc,nproc,orbs,comms,Glr%wfd,psi,in)
!!  !untranspose psi
!!  call untranspose_v(iproc,nproc,orbs,Glr%wfd,comms,psi,work=hpsi)
!!
!!
!!  alpha=2.d0
!!  energy=1.d10
!!  gnrm=1.d10
!!  gnrm_zero=0.0_gp
!!  ekin_sum=0.d0 
!!  epot_sum=0.d0 
!!  eproj_sum=0.d0
!!  !minimum value of the energy during the minimisation procedure
!!  energy_min=1.d10
!!  !local variable for the diis history
!!  idsx_actual=in%idsx
!!  !number of switching betweed DIIS and SD during self-consistent loop
!!  ndiis_sd_sw=0
!!  !previous value of idsx_actual to control if switching has appeared
!!  idsx_actual_before=idsx_actual
!!
!!  wfn_loop: do iter=1,in%itermax
!!
!!     if (iproc == 0 .and. verbose > 0) then 
!!        write( *,'(1x,a,i0)') &
!!             & repeat('~',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
!!     endif
!!     !control whether the minimisation iterations ended
!!     endloop= gnrm <= in%gnrm_cv .or. iter == in%itermax
!!     
!!     !control how many times the DIIS has switched into SD
!!     if (idsx_actual /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1
!!
!!     !terminate SCF loop if forced to switch more than once from DIIS to SD
!!     endloop=endloop .or. ndiis_sd_sw > 2
!!
!!     call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3pi,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,in%nspin,&
!!          orbs%norb,orbs%norbp,ngatherarr,pot_ion,potential)
!!
!!     call HamiltonianApplication(iproc,nproc,atoms,orbs,in%hx,in%hy,in%hz,rxyz,&
!!          nlpspd,proj,Glr,ngatherarr,potential,psi,hpsi,ekin_sum,epot_sum,eexctX,&
!!          eproj_sum,in%nspin,GPU)
!!
!!     energybs=ekin_sum+epot_sum+eproj_sum
!!     energy_old=energy
!!     energy=energybs-eexctX
!!
!!     !check for convergence or whether max. numb. of iterations exceeded
!!     if (endloop) then 
!!        if (iproc == 0) then 
!!           if (verbose > 1) write( *,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
!!           write( *,'(1x,a)') &
!!                '------------------------------------------- End of Virtual Wavefunction Optimisation'
!!           write( *,'(1x,a,3(1x,1pe18.11))') &
!!                'final  ekin,  epot,  eproj ',ekin_sum,epot_sum,eproj_sum
!!           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') &
!!                'FINAL iter,total "energy",gnrm',iter,energy,gnrm
!!           !write(61,*)hx,hy,hz,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
!!           if (energy > energy_min) write( *,'(1x,a,1pe9.2)')&
!!                'WARNING: Found an "energy" value lower than the FINAL "energy", delta:',energy-energy_min
!!        end if
!!        exit wfn_loop 
!!     endif
!!
!!     !evaluate the functional of the wavefucntions and put it into the diis structure
!!     !the energy values is printed out here
!!     call calculate_energy_and_gradient(iter,iproc,nproc,orbs,comms,GPU,Glr,in%hx,in%hy,in%hz,in%ncong,in%iscf,&
!!          ekin_sum,epot_sum,eproj_sum,0.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp,&
!!          psi,psit,hpsi,gnrm,gnrm_zero,diis%energy)
!!
!!     !control the previous value of idsx_actual
!!     idsx_actual_before=diis%idsx
!!
!!     call hpsitopsi(iproc,nproc,orbs,Glr,comms,iter,diis,in%idsx,psi,psit,hpsi,in%nspin,in)
!!
!!     write(itername,'(i4.4)')iter
!!     call plot_wf_sandbox('iter'//itername,1,atoms,Glr,hxh,hyh,hzh,rxyz,psi,'')
!!
!!     tt=(energybs-scprsum)/scprsum
!!     if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
!!          (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
!!        write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
!!             'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
!!     endif
!!!!$     if (iproc.eq.0) then
!!!!$        if (verbose > 0) then
!!!!$           write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
!!!!$                ekin_sum,epot_sum,eproj_sum
!!!!$        end if
!!!!$        write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total "energy",gnrm',iter,energy,gnrm
!!!!$     endif
!!
!!  end do wfn_loop
!!  if (iter == in%itermax .and. iproc == 0 ) &
!!       write( *,'(1x,a)')'No convergence within the allowed number of minimization steps'
!!
!!  !this deallocates also hpsivirt and psitvirt
!!  call last_orthon(iproc,nproc,orbs,Glr%wfd,in%nspin,&
!!       comms,psi,hpsi,psit,evsum)
!!  
!!  if (in%idsx > 0) then
!!     call deallocate_diis_objects(diis)
!!  end if
!!
!!  if (nproc > 1) then
!!     i_all=-product(shape(psit))*kind(psit)
!!     deallocate(psit,stat=i_stat)
!!     call memocc(i_stat,i_all,'psit',subname)
!!  end if
!!
!!  i_all=-product(shape(psi))*kind(psi)
!!  deallocate(psi,stat=i_stat)
!!  call memocc(i_stat,i_all,'psi',subname)
!!
!!  i_all=-product(shape(proj))*kind(proj)
!!  deallocate(proj,stat=i_stat)
!!  call memocc(i_stat,i_all,'proj',subname)
!!
!!  i_all=-product(shape(orbs%eval))*kind(orbs%eval)
!!  deallocate(orbs%eval,stat=i_stat)
!!  call memocc(i_stat,i_all,'eval',subname)
!!
!!
!!  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
!!  deallocate(nscatterarr,stat=i_stat)
!!  call memocc(i_stat,i_all,'nscatterarr',subname)
!!  i_all=-product(shape(ngatherarr))*kind(ngatherarr)
!!  deallocate(ngatherarr,stat=i_stat)
!!  call memocc(i_stat,i_all,'ngatherarr',subname)
!!
!!  i_all=-product(shape(radii_cf))*kind(radii_cf)
!!  deallocate(radii_cf,stat=i_stat)
!!  call memocc(i_stat,i_all,'radii_cf',subname)
!!
!!
!!  call deallocate_lr(Glr)
!!  call deallocate_comms(comms)
!!  call deallocate_orbs(orbs)
!!  call deallocate_atoms_scf(atoms,subname) 
!!  call deallocate_proj_descr(nlpspd,subname)
!!
!!  if (dokernel) then
!!     i_all=-product(shape(pkernel))*kind(pkernel)
!!     deallocate(pkernel,stat=i_stat)
!!     call memocc(i_stat,i_all,'kernel',subname)
!!  end if
!!
!!  i_all=-product(shape(pot_ion))*kind(pot_ion)
!!  deallocate(pot_ion,stat=i_stat)
!!  call memocc(i_stat,i_all,'pot_ion',subname)
!!  
!!  call deallocate_atoms(atoms,subname) 
!!  i_all=-product(shape(rxyz))*kind(rxyz)
!!  deallocate(rxyz,stat=i_stat)
!!  call memocc(i_stat,i_all,'rxyz',subname)
!!  call free_input_variables(in)
!!
!!  !finalize memory counting
!!  call memocc(0,0,'count','stop')

end program sandbox


subroutine psi_from_gaussians(iproc,nproc,at,orbs,lr,rxyz,hx,hy,hz,nspin,psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(orbs%npsidim), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='psi_from_gaussians'
  logical ::  randinp
  integer :: iorb,icoeff,i_all,i_stat,jproc,nwork,info,jorb,i,j
  real(kind=4) :: tt
  real(wp), dimension(:,:), allocatable :: gaucoeffs
  real(gp), dimension(:), allocatable :: work,ev
  real(gp), dimension(:,:), allocatable :: ovrlp
  type(gaussian_basis) :: G
  real(wp), dimension(:), pointer :: gbd_occ


  !initialise some coefficients in the gaussian basis
  !nullify the G%rxyz pointer
  nullify(G%rxyz)
  !extract the gaussian basis from the pseudowavefunctions
  !use a better basis than the input guess
  call gaussian_pswf_basis(31,.false.,iproc,nspin,at,rxyz,G,gbd_occ)

  gaucoeffs = f_malloc((/ G%ncoeff, orbs%norbp*orbs%nspinor /),id='gaucoeffs')

  !in view of complete gaussian calculation
  ovrlp = f_malloc((/ G%ncoeff, G%ncoeff /),id='ovrlp')


  !the kinetic overlap is correctly calculated only with Free BC
  randinp = .true.!.false.!lr%geocode /= 'F'

  if (randinp) then
     !fill randomly the gaussian coefficients for the orbitals considered
     do iorb=1,orbs%norbp*orbs%nspinor
        do icoeff=1,G%ncoeff
           !be sure to call always a different random number
           do jproc=0,iproc-1
              call random_number(tt)
           end do
           call random_number(tt)
           !modification: initialisation coefficient equal to delta
           if (icoeff==iorb) then
              gaucoeffs(icoeff,iorb)=1.0_gp!real(tt,wp)
           else
              gaucoeffs(icoeff,iorb)=0.0_gp
           end if
           do jproc=iproc+1,nproc-1
              call random_number(tt)
           end do
        end do
     end do

     !othogonalise the gaussian basis (wrong with k-points)
     !call gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)
     
     !calculate the overlap matrix
     call gaussian_overlap(G,G,ovrlp)

     !print it 
     write(*,'(A37)') 'Overlap matrix with gaussian_overlap:'
     do i=1,G%ncoeff
        write(*,'(i4,10(1x,1pe19.12))')i,(ovrlp(i,j),j=1,G%ncoeff)
     end do

     !calculate the overlap matrix
     call kinetic_overlap(G,G,ovrlp)

     !print it 
     write(*,'(A37)') 'Overlap matrix with kinetic_overlap:'
     do i=1,G%ncoeff
        write(*,'(i4,10(1x,1pe19.12))')i,(ovrlp(i,j),j=1,G%ncoeff)
     end do


  else
     !as an alternative strategy we may take the eigenvectors of the kinetic+k hamiltonian


     !overlap calculation of the kinetic operator, upper triangular part
     !call kinetic_overlap(G,G,ovrlp)
     call gaussian_overlap(G,G,ovrlp)
     nwork=3*G%ncoeff+1
     work = f_malloc(nwork,id='work')
     ev = f_malloc(G%ncoeff,id='ev')

!!$  if (iproc == 0) then
!!$     do iat=1,G%ncoeff
!!$        write(*,'(a,i0,10(1pe15.7))')'T',iat,ovrlp(1:iat,iat)
!!$     end do
!!$  end if

     !print *,'nwork',nwork,3*nbasis-1
     call dsyev('V','U',G%ncoeff,ovrlp(1,1),G%ncoeff,ev(1),work(1),nwork,info)
     if (info /= 0) then
        if (iproc == 0) then
           write(*,*)'DSyev Error',info
        end if
        stop
     end if

!!$  if (iproc == 0) then
!!$     do iat=1,G%ncoeff
!!$        write(*,'(a,i0,10(1pe15.7))')'Ev',iat,ovrlp(:,iat)
!!$     end do
!!$     do iat=1,G%ncoeff
!!$        write(*,'(a,i0,10(1pe15.7))')'K',iat,ev(iat)
!!$     end do
!!$  end if

     !copy the eigenvectors to the matrix
     call to_zero(G%ncoeff*orbs%norbp*orbs%nspinor,gaucoeffs)
     if (orbs%norb > G%ncoeff) stop 'wrong gaussian basis'
     jorb=mod(orbs%isorb,orbs%norb)
     do iorb=1,orbs%norbp
        jorb=jorb+1
        if (jorb == orbs%norb+1) jorb=1 !for k-points calculation
        call vcopy(G%ncoeff,ovrlp(1,jorb),1,gaucoeffs(1,orbs%nspinor*(iorb-1)+1),orbs%nspinor)
     end do

     call f_free(work)
     call f_free(ev)

     !call MPI_BARRIER(MPI_COMM_WORLD,info)
     !stop

  end if

  call f_free(ovrlp)

!WARNING: not correct!
  call gaussians_to_wavelets_new(iproc,nproc,lr,orbs,G,&
       gaucoeffs,psi)
  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G)

  !deallocate gaussian array
  call f_free(gaucoeffs)
  call f_free_ptr(gbd_occ)

END SUBROUTINE psi_from_gaussians


subroutine plot_wf_sandbox(orbname,nexpo,at,lr,hxh,hyh,hzh,rxyz,psi,comment)
  use module_base
  use module_types
  implicit none
  character(len=10) :: comment
  character(len=*) :: orbname
  integer, intent(in) :: nexpo
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  !local variables
  character(len=*), parameter :: subname='plot_wf'
  integer :: i_stat,i_all,tmp
  integer :: nl1,nl2,nl3,n1i,n2i,n3i,n1,n2,n3,i1,i2,i3,nu1,nu2,nu3
  real(gp) :: x,y,z
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:,:), allocatable :: psir

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i

  if (at%geocode == 'F') then
     nl1=14
     nu1=15
     nl2=14
     nu2=15
     nl3=14
     nu3=15
  else if (at%geocode == 'S') then
     nl1=0
     nu1=0
     nl2=14
     nu2=15
     nl3=0
     nu3=0
  else if (at%geocode == 'P') then
     nl1=0
     nu1=0
     nl2=0
     nu2=0
     nl3=0
     nu3=0
  end if

  call initialize_work_arrays_sumrho(1,lr,.true.,w)
 
  psir = f_malloc((/ -nl1.to.2*n1+1+nu1, -nl2.to.2*n2+1+nu2, -nl3.to.2*n3+1+nu3 /),id='psir')
  !initialisation
  if (lr%geocode == 'F') then
     call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if
 
  call daub_to_isf(lr,w,psi,psir)
 
  open(unit=22,file=trim(orbname),status='unknown')

  do i3=-nl3,2*n3+1+nu3
     tmp = i3 + 2*lr%ns3
     z=hzh*real(tmp,gp)
     do i2=-nl2,2*n2+1+nu2
        tmp = i2 + 2*lr%ns2
        y=hyh*real(tmp,gp)
        do i1=-nl1,2*n1+1+nu1
           tmp = i1 + 2*lr%ns1
           x=hxh*real(tmp,gp)
           if (y == rxyz(2,1) .and. z ==  rxyz(3,1)) then
!              print *,'Localization:',x,y,z,rxyz(1,1),rxyz(2,1),rxyz(3,1)
              write(22,'(1x,f9.5,1pe18.10)')x,psir(i1,i2,i3)**nexpo
           end if
        end do
     end do
  end do
  close(22)

  call f_free(psir)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf_sandbox
