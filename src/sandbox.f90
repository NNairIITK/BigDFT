!!****p* BigDFT/sandbox
!!
!! COPYRIGHT
!!    Copyright (C) 2011 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
program sandbox
  use BigDFT_API
  use Poisson_Solver
  implicit none
  character(len=*), parameter :: subname='sandbox'
  logical :: dokernel=.false.
  logical :: endloop
  integer :: n1i,n2i,n3i,iproc,nproc,ifine,i_stat,i_all,nelec
  integer :: n3d,n3p,n3pi,i3xcsh,i3s,n1,n2,n3,ndegree_ip
  integer :: idsx_actual,ndiis_sd_sw,idsx_actual_before,iter,iorb,jorb
  real(gp) :: hxh,hyh,hzh
  real(gp) :: tt,gnrm,gnrm_zero,epot_sum,eexctX,ekin_sum,eproj_sum,alpha
  real(gp) :: energy,energy_min,energy_old,energybs,evsum,scprsum
  type(atoms_data) :: atoms
  type(input_variables) :: in
  type(orbitals_data) :: orbs
  type(locreg_descriptors) :: Glr
  type(nonlocal_psp_descriptors) :: nlpspd
  type(communications_arrays) :: comms
  type(GPU_pointers) :: GPU
  type(diis_objects) :: diis
  character(len=4) :: itername
  real(gp), dimension(3) :: shift
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(gp), dimension(:,:), allocatable :: radii_cf
  real(wp), dimension(:), pointer :: hpsi,psit,psi,psidst,hpsidst,proj
  real(dp), dimension(:), pointer :: pkernel,pot_ion
  real(gp), dimension(:,:), pointer :: rxyz
  ! arrays for DIIS convergence accelerator
  real(wp), dimension(:,:,:), pointer :: ads
  !##########################################
  ! My variables
  !##########################################
  integer :: ilr
  integer, parameter :: nlr=2,alr=1,blr=2
  integer :: ldim   ! dimension of lpsi
  integer :: fmin   ! min(1,nseg_f)
  integer :: isovrlp
  real(wp) :: scpr
  integer,dimension(:),allocatable :: projflg
  integer,dimension(3,nlr) :: outofzone
  real(gp), dimension(3,2) :: cxyz
!  real(gp), dimension(nlr),parameter :: locrad=(/16.0, 16.0 /)
  real(wp),dimension(nlr),parameter :: locrad=(/ 2.565213453, 2.565213453  /) !it has dimension nlr-1, because the last element of Llr is overlap
  real(wp), dimension(:),allocatable :: lpsi    ! local projection of |Psi>
  real(wp), dimension(:),allocatable :: lhpsi   ! local projection of H|Psi>
  real(wp), dimension(:),allocatable :: lppsi   ! local projection of H|Psi>
  real(wp), dimension(:),allocatable :: ppsi   ! local projection of H|Psi>
  real(wp), dimension(:), pointer :: Lproj
  type(locreg_descriptors), dimension(nlr) :: Llr
  type(locreg_descriptors), allocatable :: Olr(:)
  real(wp), dimension(:), pointer :: potential
  real :: sum_pot !debug
  integer :: ii    ! debug  
  type(nonlocal_psp_descriptors) :: Lnlpspd
  type(atoms_data) :: Latoms

  !for the moment no need to have parallelism
  iproc=0
  nproc=1

  !initalise the variables for the calculation

  call read_input_variables(iproc,'posinp', &
       & "input.dft", "input.kpt","input.mix", "input.geopt", "input.perf", in, atoms, rxyz)

  if (iproc == 0) then
     call print_general_parameters(in,atoms)
  end if
       
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call system_properties(iproc,nproc,in,atoms,orbs,radii_cf,nelec)

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,atoms,rxyz,radii_cf,in%crmult,in%frmult,in%hx,in%hy,in%hz,&
       Glr,shift)

  !variables substitution for the PSolver part
  hxh=0.5d0*in%hx
  hyh=0.5d0*in%hy
  hzh=0.5d0*in%hz
  n1i=Glr%d%n1i
  n2i=Glr%d%n2i
  n3i=Glr%d%n3i

  n1=Glr%d%n1
  n2=Glr%d%n2
  n3=Glr%d%n3

  call timing(iproc,'CrtDescriptors','ON')
  call createWavefunctionsDescriptors(iproc,in%hx,in%hy,in%hz,&
       atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr)
  call timing(iproc,'CrtDescriptors','OF')

  ! Calculate all projectors, or allocate array for on-the-fly calculation
  call timing(iproc,'CrtProjectors ','ON')
  call createProjectorsArrays(iproc,Glr,rxyz,atoms,orbs,&
       radii_cf,in%frmult,in%frmult,in%hx,in%hy,in%hz,nlpspd,proj)
  call timing(iproc,'CrtProjectors ','OF')

  !allocate communications arrays
  call orbitals_communicators(iproc,nproc,Glr,orbs,comms)  

  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ngatherarr,'ngatherarr',subname)

  call createDensPotDescriptors(iproc,nproc,atoms%geocode,'D',&
       Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  !commented out, to be used in the future
  if (dokernel) then
     ndegree_ip=16 !default value 
     call createKernel(iproc,nproc,atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,&
          pkernel)
  else
     nullify(pkernel)
  end if

  !allocate ionic potential
  if (n3pi > 0) then
     allocate(pot_ion(n1i*n2i*n3pi+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  else
     allocate(pot_ion(1+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  end if

  allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)
  if (nproc == 1) then
     nullify(psit)
  else
     stop 'for the moment only sequential runs'
  end if

  ! allocate arrays necessary for DIIS convergence acceleration
  !the allocation with npsidim is not necessary here since DIIS arrays
  !are always calculated in the transposed form
  if (in%idsx > 0) then
     call allocate_diis_objects(in%idsx,sum(comms%ncntt(0:nproc-1)),orbs%nkptsp,diis,subname)  
  endif

  sum_pot = 0.0d0
  do ii=1,n1i*n2i*n3pi
     sum_pot = sum_pot + pot_ion(ii)
  end do
  print *,'Sum of the potential :',sum_pot

  !create input guess wavefunction
  call psi_from_gaussians(iproc,nproc,atoms,orbs,Glr,rxyz,in%hx,in%hy,in%hz,in%nspin,psi)

  !calculate scalar product between wavefunctions
  !this scheme works only in sequential
  write(*,'(A26)') 'Overlap matrix with ddot:'
  do iorb=1,orbs%norb
     do jorb=iorb,orbs%norb
        write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
             dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
             psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
     end do
  end do
 
!#######################################################################################################################
! Test the overlap calculation of wavefunctions defined in different localisation regions
!########################################################################################################################
 write(*,'(A32)') 'Entering the localisation code.'

  cxyz(:,1) = rxyz(:,1)
  cxyz(:,2) = rxyz(:,size(rxyz,2))

! Write some physical information on the Glr
  write(*,'(a24,3i4)')'Global region n1,n2,n3:',Glr%d%n1,Glr%d%n2,Glr%d%n3
  write(*,'(a27,f6.2,f6.2,f6.2)')'Global dimension (x,y,z):',Glr%d%n1*in%hx,Glr%d%n2*in%hy,Glr%d%n3*in%hz
  write(*,'(a17,f12.2)')'Global volume: ',Glr%d%n1*in%hx*Glr%d%n2*in%hy*Glr%d%n3*in%hz
  print *,'Global statistics:',Glr%wfd%nseg_c,Glr%wfd%nseg_f,Glr%wfd%nvctr_c,Glr%wfd%nvctr_f
   
! First, determine the localisation regions
! NOTE: we only pass nlr-1 as the dimension of Llr, while it is nrl 
  call determine_locreg2(nlr,rxyz,locrad,in%hx,in%hy,in%hz,Glr,Llr,outofzone)

  print *,'Outside determine_locreg2'
! Plot the localization regions
!  call draw_locregs(2,in%hx,in%hy,in%hz,Llr)

! Second, calculate the number of overlap regions using periodicity
  call get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr,outofzone)

  write(*,*)'Outside get_number_of_overlap_region:',isovrlp

  if(isovrlp > 0) then

     ! allocate the overlap region descriptors (Olr)
     allocate(Olr(isovrlp),stat=i_stat)
!     call memocc(i_stat,Olr,'Olr',subname)
   
     ! Third, construct the overlap region descriptors
     call get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr,outofzone)
  
!    Write some physical information on the overlap
     do ilr=1,isovrlp
        write(*,'(a25,3i4)')'Overlap region',ilr
        write(*,*)' n1,n2,n3:',Olr(ilr)%d%n1,Olr(ilr)%d%n2,Olr(ilr)%d%n3
        write(*,'(a27,f6.2,f6.2,f6.2)')'Overlap dimension (x,y,z):',Olr(ilr)%d%n1*in%hx,Olr(ilr)%d%n2*in%hy,Olr(ilr)%d%n3*in%hz
        write(*,'(a17,f12.2)')'Overlap volume: ',Olr(ilr)%d%n1*in%hx*Olr(ilr)%d%n2*in%hy*Olr(ilr)%d%n3*in%hz
        write(*,*)'Overlap starting point:',Olr(ilr)%ns1,Olr(ilr)%ns2,Olr(ilr)%ns3
     end do

!    Plot an overlap region
     call draw_locregs(1,in%hx,in%hy,in%hz,Olr(1))

! Third, transform the wavefunction to overlap regions
     do ilr=1,isovrlp
       ldim = (Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*orbs%norb*orbs%nspinor

       ! Allocate the local wavefunction (in one overlap region)
       allocate(lpsi(ldim+ndebug), stat=i_stat)
       call memocc(i_stat,lpsi,'lpsi',subname)

       ! Project the wavefunction inside the overlap region
       call psi_to_locreg(Glr,ilr,ldim,Olr,lpsi,nlr,orbs,psi)

! Calculate overlap in localization region using ddot
     write(*,'(A26)') 'Overlap matrix with ddot:'
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
                dot(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f,lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
                lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
        end do
     end do

! Calculate overlap in localization region using wpdot
     write(*,'(A27)')'Overlap matrix with wpdot:'
     fmin = min(Olr(ilr)%wfd%nseg_f,1)  ! checks if there is some fine_grid in the region, if not, do not shift keyv
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           call wpdot(Olr(ilr)%wfd%nvctr_c,Olr(ilr)%wfd%nvctr_f,Olr(ilr)%wfd%nseg_c,Olr(ilr)%wfd%nseg_f,&
&                  Olr(ilr)%wfd%keyv(1),Olr(ilr)%wfd%keyv(Olr(ilr)%wfd%nseg_c+fmin),Olr(ilr)%wfd%keyg(1,1),&
&                  Olr(ilr)%wfd%keyg(1,Olr(ilr)%wfd%nseg_c+fmin),lpsi(1+(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)),&
&                  lpsi(Olr(ilr)%wfd%nvctr_c+fmin+(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)),&
&                  Olr(ilr)%wfd%nvctr_c,Olr(ilr)%wfd%nvctr_f,Olr(ilr)%wfd%nseg_c,Olr(ilr)%wfd%nseg_f,&
&                  Olr(ilr)%wfd%keyv(1),Olr(ilr)%wfd%keyv(Olr(ilr)%wfd%nseg_c+fmin),Olr(ilr)%wfd%keyg(1,1),&
&                  Olr(ilr)%wfd%keyg(1,Olr(ilr)%wfd%nseg_c+fmin),lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),&
&                  lpsi(Olr(ilr)%wfd%nvctr_c+(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+fmin),scpr)
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,scpr
        end do
     end do

! Plot the orbitals that are inside the overlap region
! plot first orbital
!     call plot_wf_sandbox('orb1_ovrlp',1,atoms,Olr(ilr),hxh,hyh,hzh,rxyz,lpsi,'')
! plot second orbital
!     call plot_wf_sandbox('orb2_ovrlp',1,atoms,Olr(ilr),hxh,hyh,hzh,rxyz,&
!&                     lpsi(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f+1),'')

! Plot first orbital
     call plot_wf_sandbox('iter0-1',1,atoms,Glr,hxh,hyh,hzh,rxyz,psi,'')
! Plot second orbital
     call plot_wf_sandbox('iter0-2',1,atoms,Glr,hxh,hyh,hzh,rxyz,psi(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f+1),'')

     allocate(orbs%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
     call memocc(i_stat,orbs%eval,'eval',subname)

! Fake Eigenvalues for orbitals
     orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0


! #################
!  Kinetic part
! #################

  !allocate the potential in the full box
     call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3pi,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,in%nspin,&
          orbs%norb,orbs%norbp,ngatherarr,pot_ion,potential)

  !put in hpsi the kinetic operator (pot_ion = 0d.0)
     call HamiltonianApplication(iproc,nproc,atoms,orbs,in%hx,in%hy,in%hz,rxyz,&
         nlpspd,proj,Glr,ngatherarr,potential,psi,hpsi,ekin_sum,epot_sum,eexctX,&
         eproj_sum,in%nspin,GPU)

  !calculate scalar product between wavefunctions
  !this scheme works only in sequential
     write(*,*)' '
     write(*,'(A34)') 'Kinetic overlap matrix with ddot:'
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
                dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
                hpsi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
        end do
     end do

! Now transform hpsi to overlap region
     allocate(lhpsi(ldim+ndebug), stat=i_stat)
     call memocc(i_stat,lhpsi,'lhpsi',subname)

     call psi_to_locreg(Glr,ilr,ldim,Olr,lhpsi,nlr,orbs,hpsi)

     write(*,'(A51)') 'Kinetic overlap matrix with ddot and localization:'
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
               dot(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f,lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
               lhpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
        end do
     end do

! ##################
! LOCAL: Non-local Projectors part
! ##################

! allocate  projflg
  allocate(projflg(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,projflg,'projflg',subname)

! Need to restrict non-local projectors to overlap region
! First make the descriptors
  call nlpspd_to_locreg(in,iproc,Glr,Olr(ilr),rxyz,atoms,orbs,&
&       radii_cf,in%frmult,in%frmult,in%hx,in%hy,in%hz,nlpspd,Lnlpspd,projflg)

  print *,'outside of nlpspd_to_locreg'
  print *,'nvctrp',Lnlpspd%nvctr_p
  print *,'nseg_p',Lnlpspd%nseg_p

! Allocate Lproj and Lppsi
  allocate(Lproj(Lnlpspd%nprojel+ndebug),stat=i_stat)
  call memocc(i_stat,Lproj,'Lproj',subname)
  allocate(lppsi(ldim+ndebug), stat=i_stat)
  call memocc(i_stat,lppsi,'lppsi',subname)
  allocate(ppsi((Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f)*orbs%norb*orbs%nspinor+ndebug), stat=i_stat)
  call memocc(i_stat,ppsi,'ppsi',subname)

! initialise lppsi and ppsi
  call razero((Olr(ilr)%wfd%nvctr_c + 7*Olr(ilr)%wfd%nvctr_f)*orbs%norb*orbs%nspinor,lppsi)
  call razero((Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f)*orbs%norb*orbs%nspinor,ppsi)

! Apply without the localization (for comparison)
  eproj_sum = 0.0_gp
  call applyprojectorsonthefly(iproc,orbs,atoms,Glr,&
          rxyz,in%hx,in%hy,in%hz,Glr%wfd,nlpspd,proj,psi,ppsi,eproj_sum)

! Calculate dotprod: <Psi_a|p><p|Psi_b>
   write(*,'(A37)') 'NL-Operator overlap matrix with ddot:'
   do iorb=1,orbs%norb
      do jorb=iorb,orbs%norb
         write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
             dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
             ppsi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
      end do
   end do

  print *,'Lnlpspd:'
! Fill and Apply the projectors on the wavefunctions
  call apply_local_projectors(atoms,in,Olr(ilr),Lnlpspd,Lproj,orbs,projflg,lpsi,rxyz,lppsi)

! Calculate dotprod: <Psi_a|p><p|Psi_b>
 write(*,'(A54)') 'NL-Operator overlap matrix with ddot and localization:'
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
               dot(Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f,lpsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
               lppsi((Olr(ilr)%wfd%nvctr_c+7*Olr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
        end do
     end do

! Deallocations
     deallocate(lpsi,stat=i_stat)
     call memocc(i_stat,i_all,'lpsi',subname)

     deallocate(lhpsi,stat=i_stat)
     call memocc(i_stat,i_all,'lhpsi',subname)
   end do ! ilr
  end if   !isovrlp

  stop  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------------------------------------------
!#######################################################################################################################
! End of Localization part
!######################################################################################################################



  !othogonalise them
  !transpose the psi wavefunction
  call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,&
       psi,work=hpsi)
  call orthogonalize(iproc,nproc,orbs,comms,Glr%wfd,psi,in)
  !untranspose psi
  call untranspose_v(iproc,nproc,orbs,Glr%wfd,comms,psi,work=hpsi)


  alpha=2.d0
  energy=1.d10
  gnrm=1.d10
  gnrm_zero=0.0_gp
  ekin_sum=0.d0 
  epot_sum=0.d0 
  eproj_sum=0.d0
  !minimum value of the energy during the minimisation procedure
  energy_min=1.d10
  !local variable for the diis history
  idsx_actual=in%idsx
  !number of switching betweed DIIS and SD during self-consistent loop
  ndiis_sd_sw=0
  !previous value of idsx_actual to control if switching has appeared
  idsx_actual_before=idsx_actual

  wfn_loop: do iter=1,in%itermax

     if (iproc == 0 .and. verbose > 0) then 
        write( *,'(1x,a,i0)') &
             & repeat('~',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
     endif
     !control whether the minimisation iterations ended
     endloop= gnrm <= in%gnrm_cv .or. iter == in%itermax
     
     !control how many times the DIIS has switched into SD
     if (idsx_actual /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1

     !terminate SCF loop if forced to switch more than once from DIIS to SD
     endloop=endloop .or. ndiis_sd_sw > 2

     call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3pi,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,in%nspin,&
          orbs%norb,orbs%norbp,ngatherarr,pot_ion,potential)

     call HamiltonianApplication(iproc,nproc,atoms,orbs,in%hx,in%hy,in%hz,rxyz,&
          nlpspd,proj,Glr,ngatherarr,potential,psi,hpsi,ekin_sum,epot_sum,eexctX,&
          eproj_sum,in%nspin,GPU)

     energybs=ekin_sum+epot_sum+eproj_sum
     energy_old=energy
     energy=energybs-eexctX

     !check for convergence or whether max. numb. of iterations exceeded
     if (endloop) then 
        if (iproc == 0) then 
           if (verbose > 1) write( *,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write( *,'(1x,a)') &
                '------------------------------------------- End of Virtual Wavefunction Optimisation'
           write( *,'(1x,a,3(1x,1pe18.11))') &
                'final  ekin,  epot,  eproj ',ekin_sum,epot_sum,eproj_sum
           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') &
                'FINAL iter,total "energy",gnrm',iter,energy,gnrm
           !write(61,*)hx,hy,hz,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
           if (energy > energy_min) write( *,'(1x,a,1pe9.2)')&
                'WARNING: Found an "energy" value lower than the FINAL "energy", delta:',energy-energy_min
        end if
        exit wfn_loop 
     endif

     !control the previous value of idsx_actual
     idsx_actual_before=idsx_actual

     call hpsitopsi(iproc,nproc,orbs,in%hx,in%hy,in%hz,Glr,comms,in%ncong,&
          iter,diis,in%idsx,gnrm,gnrm_zero,scprsum,psi,psit,hpsi,in%nspin,GPU,in)

     write(itername,'(i4.4)')iter
     call plot_wf_sandbox('iter'//itername,1,atoms,Glr,hxh,hyh,hzh,rxyz,psi,'')

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
        end if
        write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total "energy",gnrm',iter,energy,gnrm
     endif

  end do wfn_loop
  if (iter == in%itermax .and. iproc == 0 ) &
       write( *,'(1x,a)')'No convergence within the allowed number of minimization steps'

  !this deallocates also hpsivirt and psitvirt
  call last_orthon(iproc,nproc,orbs,Glr%wfd,in%nspin,&
       comms,psi,hpsi,psit,evsum)
  
  if (in%idsx > 0) then
     call deallocate_diis_objects(diis,subname)
  end if

  if (nproc > 1) then
     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  end if

  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

  i_all=-product(shape(proj))*kind(proj)
  deallocate(proj,stat=i_stat)
  call memocc(i_stat,i_all,'proj',subname)

  i_all=-product(shape(orbs%eval))*kind(orbs%eval)
  deallocate(orbs%eval,stat=i_stat)
  call memocc(i_stat,i_all,'eval',subname)


  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
  deallocate(nscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'nscatterarr',subname)
  i_all=-product(shape(ngatherarr))*kind(ngatherarr)
  deallocate(ngatherarr,stat=i_stat)
  call memocc(i_stat,i_all,'ngatherarr',subname)

  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)


  call deallocate_lr(Glr,subname)
  call deallocate_comms(comms,subname)
  call deallocate_orbs(orbs,subname)
  call deallocate_atoms_scf(atoms,subname) 
  call deallocate_proj_descr(nlpspd,subname)

  if (dokernel) then
     i_all=-product(shape(pkernel))*kind(pkernel)
     deallocate(pkernel,stat=i_stat)
     call memocc(i_stat,i_all,'kernel',subname)
  end if

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)
  
  call deallocate_atoms(atoms,subname) 
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)
  call free_input_variables(in)

  !finalize memory counting
  call memocc(0,0,'count','stop')


end program sandbox
!!***

!!****f* BigDFT/psi_from_gaussians
!! FUNCTION
!!
!! SOURCE
!!
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

  allocate(gaucoeffs(G%ncoeff,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)

  !in view of complete gaussian calculation
  allocate(ovrlp(G%ncoeff,G%ncoeff),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)


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
     allocate(work(nwork+ndebug),stat=i_stat)
     call memocc(i_stat,work,'work',subname)
     allocate(ev(G%ncoeff+ndebug),stat=i_stat)
     call memocc(i_stat,ev,'ev',subname)

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
     call razero(G%ncoeff*orbs%norbp*orbs%nspinor,gaucoeffs)
     if (orbs%norb > G%ncoeff) stop 'wrong gaussian basis'
     jorb=mod(orbs%isorb,orbs%norb)
     do iorb=1,orbs%norbp
        jorb=jorb+1
        if (jorb == orbs%norb+1) jorb=1 !for k-points calculation
        call dcopy(G%ncoeff,ovrlp(1,jorb),1,gaucoeffs(1,orbs%nspinor*(iorb-1)+1),orbs%nspinor)
     end do

     i_all=-product(shape(work))*kind(work)
     deallocate(work,stat=i_stat)
     call memocc(i_stat,i_all,'work',subname)
     i_all=-product(shape(ev))*kind(ev)
     deallocate(ev,stat=i_stat)
     call memocc(i_stat,i_all,'ev',subname)

     !call MPI_BARRIER(MPI_COMM_WORLD,info)
     !stop

  end if

  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)


  call gaussians_to_wavelets_new(iproc,nproc,lr,orbs,hx,hy,hz,G,&
       gaucoeffs,psi)
  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)

  !deallocate gaussian array
  i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
  deallocate(gaucoeffs,stat=i_stat)
  call memocc(i_stat,i_all,'gaucoeffs',subname)
  i_all=-product(shape(gbd_occ))*kind(gbd_occ)
  deallocate(gbd_occ,stat=i_stat)
  call memocc(i_stat,i_all,'gbd_occ',subname)

  
END SUBROUTINE psi_from_gaussians
!!***


!!****f* BigDFT/plot_wf_sandbox
!! FUNCTION
!!
!! SOURCE
!!
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

  call initialize_work_arrays_sumrho(lr,w)
 
  allocate(psir(-nl1:2*n1+1+nu1,-nl2:2*n2+1+nu2,-nl3:2*n3+1+nu3+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
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

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf_sandbox
!%***

!#############################################################################################################################################
!!****f* BigDFT/psi_to_locreg
!#############################################################################################################################################
!! FUNCTION: Tranform wavefunction between Global region and localisation region
!!
!! WARNING: 
!!         Only coded for sequential, not parallel cases !! For parallel should change increment and loc_psi dimensions
!! SOURCE:
!!
subroutine psi_to_locreg(Glr,ilr,ldim,Olr,lpsi,nlr,orbs,psi)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: nlr                  ! number of localization regions
  integer :: ilr           ! index of the localization region we are considering
  integer :: ldim          ! dimension of lpsi 
  type(orbitals_data),intent(in) :: orbs      ! orbital descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Olr  ! Localization grid descriptors 
  real(wp),dimension(orbs%npsidim),intent(in) :: psi       !Wavefunction (compressed format)
  real(wp),dimension(ldim),intent(inout) :: lpsi !Wavefunction in localization region
  !#############################################
  !local variables
  !############################################
  integer :: igrid,isegloc,isegG,ix,iorbs
  integer :: lmin,lmax,Gmin,Gmax
  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
  integer :: length      ! Length of the overlap between Lseg and Gseg
  integer :: lincrement  ! Increment for writing orbitals in loc_psi
  integer :: Gincrement  ! Increment for reading orbitals in psi
  integer :: nseg        ! total number of segments in Llr
  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
  character(len=*), parameter :: subname='psi_to_locreg'
  integer :: i_stat,i_all
  integer :: start,Gstart
  integer :: lfinc,Gfinc

! Define integers
  nseg = Olr(ilr)%wfd%nseg_c + Olr(ilr)%wfd%nseg_f
  lincrement = Olr(ilr)%wfd%nvctr_c + 7*Olr(ilr)%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0

! Initialize loc_psi
  call razero(lincrement*orbs%norb*orbs%nspinor,lpsi)
 
! Get the keymask: shift for every segment of Llr (with respect to Glr)
  allocate(keymask(2,nseg),stat=i_stat)
  call memocc(i_stat,keymask,'keymask',subname)

  call shift_locreg_indexes(Glr,Olr(ilr),keymask,nseg)

!####################################################
! Do coarse region
!####################################################
  do isegloc = 1,Olr(ilr)%wfd%nseg_c
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
 
! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
     do isegG = 1,Glr%wfd%nseg_c
        Gmin = Glr%wfd%keyg(1,isegG)
        Gmax = Glr%wfd%keyg(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle
        
        ! Define the offset between the two segments
        offset = lmin - Gmin
        if(offset < 0) then
           offset = 0
        end if
    
        ! Define the length of the two segments
        length = min(lmax,Gmax)-max(lmin,Gmin)
 
        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           ! loop over the orbitals
           do iorbs=1,orbs%norb*orbs%nspinor
              lpsi(icheck+lincrement*(iorbs-1))=psi(Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1))
           end do
        end do
     end do
  end do

! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Olr(ilr)%wfd%nvctr_c) then
    write(*,*)'There is an error in psi_to_locreg: number of coarse points used',icheck
    write(*,*)'is not equal to the number of coarse points in the region',Olr(ilr)%wfd%nvctr_c
  end if

!##############################################################
! Now do fine region
!##############################################################

  icheck = 0
  start = Olr(ilr)%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c
  lfinc  = Olr(ilr)%wfd%nvctr_f
  Gfinc = Glr%wfd%nvctr_f

  do isegloc = Olr(ilr)%wfd%nseg_c+1,nseg
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
 
! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
     do isegG = Glr%wfd%nseg_c+1,Glr%wfd%nseg_c+Glr%wfd%nseg_f

        Gmin = Glr%wfd%keyg(1,isegG)
        Gmax = Glr%wfd%keyg(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle

        offset = lmin - Gmin
        if(offset < 0) offset = 0

        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           do igrid=0,6
              do iorbs=1,orbs%norb*orbs%nspinor
                 lpsi(start+icheck+lincrement*(iorbs-1)+igrid*lfinc)=&
&                psi(Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc)
              end do
           end do
        end do
     end do
  end do
  
 ! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Olr(ilr)%wfd%nvctr_f) then
    write(*,*)'There is an error in psi_to_locreg: number of fine points used',icheck
    write(*,*)'is not equal to the number of fine points in the region',Olr(ilr)%wfd%nvctr_f
  end if

  deallocate(keymask,stat=i_stat)
  call memocc(i_stat,i_all,'keymask',subname)

END SUBROUTINE psi_to_locreg
!%***

!##############################################################################################################################################
!!****f* BigDFT/get_overlap_region
!##############################################################################################################################################
!! FUNCTION Given two localization regions, A and B, this routine returns a localization region corresponding to the intersection of A & B. 
!!
!! SOURCE
!!
subroutine get_overlap_region_free(alr,blr,Glr,isovrlp,Llr,nlr,Olr)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(out) :: isovrlp             ! True if there is an overlap
  type(locreg_descriptors),intent(out) :: Olr ! Overlap localization regions 
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  !#############################################
  !local variables
  !############################################
  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
  character(len=*), parameter :: subname='get_overlap_region'

! Set the bounds of region A
  axmin = Llr(alr)%ns1
  aymin = Llr(alr)%ns2
  azmin = Llr(alr)%ns3
  axmax = Llr(alr)%ns1 + Llr(alr)%d%n1
  aymax = Llr(alr)%ns2 + Llr(alr)%d%n2
  azmax = Llr(alr)%ns3 + Llr(alr)%d%n3

! Set the bounds of region B
  bxmin = Llr(blr)%ns1
  bymin = Llr(blr)%ns2
  bzmin = Llr(blr)%ns3
  bxmax = Llr(blr)%ns1 + Llr(blr)%d%n1
  bymax = Llr(blr)%ns2 + Llr(blr)%d%n2
  bzmax = Llr(blr)%ns3 + Llr(blr)%d%n3

! Set initial value of isovrlp
  isovrlp = 0  

! Determine if there is an overlap
! To do this, we compare axis by axis if there is an overlap.
! The cubes overlap if they overlap on all axis.
  if(((axmin .le. bxmax).and.(bxmin .le. axmax)) .and. &
&    ((aymin .le. bymax).and.(bymin .le. aymax)) .and. &
&    ((azmin .le. bzmax).and.(bzmin .le. azmax))) then
     isovrlp = 1
  end if

! Now construct the Overlap localization region descriptor
! only if there is an overlap. The following only works
! when the previous test is successful. Note also that
! isx, isy and isz are necessarily in the Glr by construction
! of the Llrs, so don't need to test them.
  if(isovrlp > 0) then
 
!   Determine the limits of the overlap region
    isx = max(axmin,bxmin)
    isy = max(aymin,bymin)
    isz = max(azmin,bymin)

    iex = min(axmax,bxmax)
    iey = min(aymax,bymax)
    iez = min(azmax,bzmax)

!   Checks to assign the geometric code of the overlap region (TO DO?)
!   This could change the values of the bounds, so do it here
!   for now, is sandbox, deal only with free boundary
    Olr%geocode = 'F'  

!   Values for the starting point of the cube
    Olr%ns1 = isx
    Olr%ns2 = isy
    Olr%ns3 = isz

!   Dimensions of the overlap region
    Olr%d%n1 = iex - isx 
    Olr%d%n2 = iey - isy 
    Olr%d%n3 = iez - isz 
    
!   Dimensions of the fine grid inside the overlap region
    if (isx < iex) then
       Olr%d%nfl1=max(isx,Glr%d%nfl1)-isx
       Olr%d%nfu1=min(iex,Glr%d%nfu1)-isx
    else
       write(*,*)'Yet to be implemented (little effort?)'
       stop
    end if

    if (isy < iey) then
       Olr%d%nfl2=max(isy,Glr%d%nfl2)-isy
       Olr%d%nfu2=min(iey,Glr%d%nfu2)-isy
    else
       write(*,*)'Yet to be implemented (little effort?)'
       stop
    end if

    if (isz < iez) then
       Olr%d%nfl3=max(isz,Glr%d%nfl3)-isz
       Olr%d%nfu3=min(iez,Glr%d%nfu3)-isz
    else
       write(*,*)'Yet to be implemented (little effort?)'
       stop
    end if

!   Dimensions of the interpolating scaling function grid 
!   (geocode already taken into acount because it is simple)
    select case(Olr%geocode)
    case('F')
       Olr%d%n1i=2*Olr%d%n1+31
       Olr%d%n2i=2*Olr%d%n2+31
       Olr%d%n3i=2*Olr%d%n3+31
    case('S')
       Olr%d%n1i=2*Olr%d%n1+2
       Olr%d%n2i=2*Olr%d%n2+31
       Olr%d%n3i=2*Olr%d%n3+2
    case('P')
       Olr%d%n1i=2*Olr%d%n1+2
       Olr%d%n2i=2*Olr%d%n2+2
       Olr%d%n3i=2*Olr%d%n3+2
    end select
 
!   Now define the wavefunction descriptors inside the overlap region
!   First calculate the number of points and segments for the region
!   Coarse part:
    call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
         Olr%wfd%nseg_c,Olr%wfd%nvctr_c)
!   Fine part:
    call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr%wfd%nseg_f,Olr%wfd%nvctr_f)

!   Now allocate the wavefunction descriptors (keyg,keyv) following the needs
    call allocate_wfd(Olr%wfd,subname)

!   At last, fill the wavefunction descriptors
!   Coarse part
     call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
          Olr%wfd%nseg_c,Olr%wfd%nvctr_c,&
          Olr%wfd%keyg(1,1),Olr%wfd%keyv(1))
!   Fine part
     call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr%wfd%nseg_f,Olr%wfd%nvctr_f,&
          Olr%wfd%keyg(1,Olr%wfd%nseg_c+min(1,Olr%wfd%nseg_f)),&
          Olr%wfd%keyv(Olr%wfd%nseg_c+min(1,Olr%wfd%nseg_f)))

!    If the localisation region is isolated build also the bounds
     if (Olr%geocode=='F') then
        call locreg_bounds(Olr%d%n1,Olr%d%n2,Olr%d%n3,&
             Olr%d%nfl1,Olr%d%nfu1,Olr%d%nfl2,Olr%d%nfu2,&
             Olr%d%nfl3,Olr%d%nfu3,Olr%wfd,Olr%bounds)
     
     end if

! If there is no overlap, write it and end routine
  else
     write(*,*)'There is no overlap between regions:',alr, blr
  end if

END SUBROUTINE get_overlap_region_free
!%***


!#############################################################################################################################################
!!****f* BigDFT/shift_locreg_indexes
!#############################################################################################################################################
!! FUNCTION:
!!        Find the shift necessary for the indexes of every segment of Blr
!!        to make them compatible with the indexes of Alr. These shifts are
!!        returned in the array keymask(nseg), where nseg should be the number
!!        of segments in Blr.
!! WARNING: 
!!         This routine supposes that the region Blr is contained in the region Alr.
!!         This should always be the case, if we concentrate on the overlap between two regions.
!! SOURCE:
!!
subroutine shift_locreg_indexes(Alr,Blr,keymask,nseg)

  use module_base
  use module_types
 
 implicit none

!############################
! Arguments
!############################
 type(locreg_descriptors),intent(in) :: Alr,Blr   ! The two localization regions
 integer,intent(in) :: nseg
 integer,intent(out) :: keymask(2,nseg)
!#############################
! Local variable
!#############################
 integer :: iseg      !integer for the loop
 integer :: Bindex    !starting index of segments in Blr
 integer :: x,y,z     !coordinates of start of segments in Blr 
 integer :: shift(3)  !shift between the beginning of the segment in Blr and the origin of Alr
 integer ::  tmp

!Big loop on all segments
 do iseg=1,nseg

!##########################################
! For the Starting index
    Bindex = Blr%wfd%keyg(1,iseg)
    tmp = Bindex -1
    z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
    tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
    y   = tmp / (Blr%d%n1+1)
    x   = tmp - y * (Blr%d%n1+1)
 
! Shift between the beginning of the segment and the start of the Alr region
    shift(1) = x + Blr%ns1 - Alr%ns1
    shift(2) = y + Blr%ns2 - Alr%ns2
    shift(3) = z + Blr%ns3 - Alr%ns3

! Write the shift in index form
    keymask(1,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1

!######################################
! For the ending index

    Bindex = Blr%wfd%keyg(2,iseg)
    tmp = Bindex -1
    z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
    tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
    y   = tmp / (Blr%d%n1+1)
    x   = tmp - y * (Blr%d%n1+1)

! Shift between the beginning of the segment and the start of the Alr region
    shift(1) = x + Blr%ns1 - Alr%ns1
    shift(2) = y + Blr%ns2 - Alr%ns2
    shift(3) = z + Blr%ns3 - Alr%ns3

! Write the shift in index form
    keymask(2,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1
 end do

END SUBROUTINE shift_locreg_indexes
!%***

!#############################################################################################################################################
!!****f* BigDFT/nlpspd_to_locreg
!#############################################################################################################################################
!! FUNCTION: Transform projectors between Global region and localisation region
!!
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine nlpspd_to_locreg(input_parameters,iproc,Glr,Llr,rxyz,atoms,orbs,&
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,Lnlpspd,projflg)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  type(input_variables),intent(in) :: input_parameters
  integer,intent(in) :: iproc
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(atoms_data),intent(in) :: atoms        ! atom descriptors
  type(orbitals_data),intent(in) :: orbs      ! orbital descriptors
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz  ! grid descriptions
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  type(nonlocal_psp_descriptors),intent(out) :: Lnlpspd  ! local descriptors for the projectors 
  !########################################
  !Subroutine Array Arguments
  !########################################
  integer,dimension(atoms%nat),intent(out) :: projflg
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz !atomic positions
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf  ! radii of the different atom types
  !#############################################
  !local variables
  !############################################
   integer :: iatom,igrid
   integer :: ii,jj,i1,i2,i3      !integers for loops
   integer :: mproj,natp
   integer :: mseg_c,mvctr_c,mseg_f,mvctr_f
   integer :: mseg !total number of segments
   integer :: iat  ! index of the atoms
   integer :: nprojelat ! total number of elements
   integer :: isx,isy,isz,iex,iey,iez
   integer :: iseg,jseg,Gseg,Gvctr
   integer :: nl1,nl2,nl3,nu1,nu2,nu3 ! bounds of projectors around atom iatom
   real(gp) :: hhx,hhy,hhz   !reals to shorten name of variables
   integer,dimension(1:2,1:2,1:3) :: bounds
   logical,dimension(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3) :: logrid
   character(len=*),parameter :: subname='nlpspd_to_locreg'

! Rename some variables
  hhx = input_parameters%hx
  hhy = input_parameters%hy
  hhz = input_parameters%hz

!Determine the number of projectors with components in locreg
! and also which atoms have such projectors and number of atoms
   call  number_of_projectors_in_locreg(atoms,cpmult,fpmult,Glr,hx,hy,hz,Llr,nlpspd,&
&        mproj,projflg,natp,radii_cf,rxyz)

   Lnlpspd%nproj = mproj

!TESTS
  print *,'Llr check:',Llr%ns1,Llr%ns2,Llr%ns3,Llr%d%n1,Llr%d%n2,Llr%d%n3
  print *,'Number of projectors:', mproj
  print *,'Projflg', projflg
!ENDTESTS

!Allocate the arrays of Lnlpspd, except keyg_p and keyv_p
 call allocate_Lnlpspd(natp,Lnlpspd,subname)

  iat = 0
  nprojelat = 0
  Lnlpspd%nprojel = 0
  mseg = 0
  do iatom = 1,atoms%nat
     if(projflg(iatom) == 0) cycle 
     iat = iat + 1  

!    Determine the bounds of the projectors
     call projector_box_in_locreg(iatom,Glr,Llr,nlpspd,bounds)

     do ii = 1,2
        do jj = 1,3
           Lnlpspd%nboxp_c(ii,jj,iat) = bounds(1,ii,jj)
           Lnlpspd%nboxp_f(ii,jj,iat) = bounds(2,ii,jj)
        end do
     end do

!    Now we can determine the number of segments and elements of coarse grid
     call fill_logrid(atoms%geocode,Glr%d%n1,Glr%d%n2,Glr%d%n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,1,&
&                     atoms%iatype(iatom),rxyz(1,iatom),radii_cf(atoms%iatype(iatom),3),cpmult,hhx,hhy,hhz,logrid)

     call number_of_projector_elements_in_locreg(iatom,1,atoms,Glr,Llr,logrid,nlpspd,mproj,mseg_c,mvctr_c)
 
     Lnlpspd%nseg_p(2*iat-1) = mseg_c
     Lnlpspd%nvctr_p(2*iat-1) = mvctr_c 

     print *,'radii:',radii_cf(atoms%iatype(iatom),3),radii_cf(1,3)*cpmult
     do i1=40,41
     do i2=40,41
     do i3=40,53
        print *,'for point:',i1,i2,i3,'logrid:',logrid(i1,i2,i3)
     end do
     end do
     end do

! Do the same for fine grid
     call fill_logrid(atoms%geocode,Glr%d%n1,Glr%d%n2,Glr%d%n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,1,&
&                     atoms%iatype(iatom),rxyz(1,iatom),radii_cf(atoms%iatype(iatom),2),fpmult,hhx,hhy,hhz,logrid)

     call number_of_projector_elements_in_locreg(iatom,2,atoms,Glr,Llr,logrid,nlpspd,mproj,mseg_f,mvctr_f)

     Lnlpspd%nseg_p(2*iat) = mseg_f
     Lnlpspd%nvctr_p(2*iat) = mvctr_f

     print *,'radii:',radii_cf(atoms%iatype(iatom),2),radii_cf(1,2)
     do i1=40,41
     do i2=40,41
     do i3=40,53
        print *,'for point:',i1,i2,i3,'logrid:',logrid(i1,i2,i3)
     end do
     end do
     end do

!    Should not be useful, because if projflg is > 0 there should be some elements
     if(mvctr_c == 0 .and. mvctr_f == 0) then
        projflg(iatom) = 0 
     end if    

     print *,'Second projflg',projflg

     nprojelat = mvctr_c*projflg(iatom) + 7*mvctr_f*projflg(iatom)
     Lnlpspd%nprojel = max(Lnlpspd%nprojel,nprojelat)
     mseg = mseg + mseg_c + mseg_f
  end do

!   Now allocate keyg_p,keyv_p following the needs
    call allocate_projd(mseg,Lnlpspd,subname)

! Renaming some variables to simply calling of routines
   !starting point of locreg
   isx = Llr%ns1
   isy = Llr%ns2
   isz = Llr%ns3
   !ending point of locreg
   iex = Llr%ns1 + Llr%d%n1
   iey = Llr%ns2 + Llr%d%n2
   iez = Llr%ns3 + Llr%d%n3
  
!  At last, fill the projector descriptors (keyg_p,keyv_p)
   iseg = 1
   iat = 0
   do iatom= 1,atoms%nat
      if(projflg(iatom) == 0) cycle
      iat = iat + 1
      print *,'nseg_c',Lnlpspd%nseg_p
      print *,'entering segkeys_loc, for iatom:',iatom,'with nseg_c',Lnlpspd%nseg_p(2*iat-1)
!     number of segments for coarse
      jseg = nlpspd%nseg_p(2*iatom-2)+1 ! index where to start in keyg for global region (nlpspd)
      Gseg = nlpspd%nseg_p(2*iatom-1)-nlpspd%nseg_p(2*iatom-2) ! number of segments for global region
      Gvctr = nlpspd%nvctr_p(2*iatom-1)-nlpspd%nvctr_p(2*iatom-2)!number of elements for global region

!     Coarse part      
      call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
           Gseg,Gvctr,nlpspd%keyg_p(1,jseg),nlpspd%keyv_p(jseg),&
           Lnlpspd%nseg_p(2*iat-1),Lnlpspd%nvctr_p(2*iat-1),&
           Lnlpspd%keyg_p(1,iseg),Lnlpspd%keyv_p(iseg))

      print *,'after first segkeys_loc(coarse)'      

      iseg = iseg + Lnlpspd%nseg_p(2*iat-1)      
      if(Lnlpspd%nseg_p(2*iat) > 0) then  !only do fine grid if present
!     Number of segments for fine
         jseg = nlpspd%nseg_p(2*iatom-1)+1 ! index where to start in keyg for global region (nlpspd)
         Gseg = nlpspd%nseg_p(2*iatom)-nlpspd%nseg_p(2*iatom-1) ! number of segments for global region
         Gvctr = nlpspd%nvctr_p(2*iatom)-nlpspd%nvctr_p(2*iatom-1)!number of elements for global region

!     Fine part
         call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
              Gseg,Gvctr,nlpspd%keyg_p(1,jseg),nlpspd%keyv_p(jseg),&
              Lnlpspd%nseg_p(2*iat),Lnlpspd%nvctr_p(2*iat),&
              Lnlpspd%keyg_p(1,iseg),Lnlpspd%keyv_p(iseg))

         iseg = iseg + Lnlpspd%nseg_p(2*iat)
      end if 
   end do

END SUBROUTINE nlpspd_to_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/number_of_projectors_in_locreg
!#############################################################################################################################################
!! FUNCTION: Calculates the number of projectors with components in the locreg
!!           It also returns a vector, projflg, which identifies the atoms with projectors inside the region
!!                      projflg = 0, no projectors in locreg
!!                      projflg = nproj, nproj projectors from atom iatom in locreg
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine number_of_projectors_in_locreg(atoms,cpmult,fpmult,Glr,hx,hy,hz,Llr,nlpspd,&
&          mproj,projflg,natp,radii_cf,rxyz)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  real(gp),intent(in) :: cpmult,fpmult,hx,hy,hz  ! grid descriptions
  type(atoms_data),intent(in) :: atoms        ! atoms descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  integer,intent(out) :: mproj  ! number of projectors
  integer,intent(out) :: natp   ! number of atoms having projectors in region
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  integer,dimension(atoms%nat),intent(out) :: projflg ! flag which is equal to the number of projectors with components inside locreg for each atom
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz !atomic positions
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf  ! radii of the different atom types
  !#######################################
  ! Local Variables
  !#######################################
  integer :: iatom,ii,izone                 ! integer for loop
  integer :: bound(1:2,1:3)           ! bound of locreg
  integer :: nproj                    ! temporary number of projectors
  integer :: i_stat                   ! allocation error
  logical :: intersect                ! logical for intersect of projector with locreg
  character(len=*), parameter :: subname='number_of_projectors_in_locreg'

  projflg = 0
  mproj = 0
  natp = 0
  do iatom=1,atoms%nat
!       check if projector of atom iatom overlap the locreg (coarse grid)        
        call check_projector_intersect_with_locreg(atoms,cpmult,hx,hy,hz,Llr,radii_cf,rxyz,intersect)

        if(intersect) then
           call numb_proj(atoms%iatype(iatom),atoms%ntypes,atoms%psppar,atoms%npspcode,nproj)
           mproj = mproj + nproj
           if(nproj > 0) then
              projflg(iatom) = nproj
              natp = natp + 1
           end if
        end if        

! Only have to do it if the atom is not yet selected
     if (projflg(iatom) .eq. 0) then
!         check if projector of atom iatom overlap the locreg (fine grid)
          call check_projector_intersect_with_locreg(Llr,intersect)

          if(intersect) then        
             call numb_proj(atoms%iatype(iatom),atoms%ntypes,atoms%psppar,atoms%npspcode,nproj)
             mproj = mproj + nproj
             if(nproj > 0) projflg(iatom) = nproj
          end if
     end if        
  end do

END SUBROUTINE number_of_projectors_in_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/fracture_projector
!#############################################################################################################################################
!! FUNCTION: Returns the limits of the various folded projector zones.
!!           
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine fracture_projector(Glr,natom,nboxp,nzones,pbox)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: natom   !number of atoms
  integer,intent(in) :: nzones  !number of zones
  type(locreg_descriptors),intent(in) :: Glr ! global region descriptor
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  integer,dimension(2,3),intent(in) :: nboxp           ! Limits of the projector box
  integer,dimension(2,3,nzones),intent(out) :: pbox    ! Limits of the projector box
  !#######################################
  ! Local Variables
  !#######################################
  integer :: izone             !integer for loops
  integer :: box(2,3,2)        !two possible starting and ending points in X,Y,Z
  
! Find limits in X
  box(1,1,2) = -99
  box(2,1,2) = -99
  if(nboxp(1,1) < 0) then
    box(1,1,1) = Glr%ns1
    box(2,1,1) = nboxp(2,1)
    box(1,1,2) = modulo(nboxp(1,1),Glr%d%n1+1)
    box(2,1,2) = Glr%ns1 + Glr%d%n1    
  else if (nboxp(2,1) > Glr%d%n1) then
    box(1,1,1) = Glr%ns1
    box(2,1,1) = modulo(nboxp(2,1),Glr%d%n1+1)
    box(1,1,2) = nboxp(1,1)
    box(2,1,2) = Glr%ns1 + Glr%d%n1
  else
    box(1,1,1) = nboxp(1,1)
    box(2,1,1) = nboxp(2,1)
  end if

!Find limits in Y
  box(1,2,2) = -99
  box(2,2,2) = -99
  if(nboxp(1,2) < 0) then
    box(1,2,1) = Glr%ns2
    box(2,2,1) = nboxp(2,1)
    box(1,2,2) = modulo(nboxp(1,2),Glr%d%n2+1)
    box(2,2,2) = Glr%ns2 + Glr%d%n2    
  else if (nboxp(2,2) > Glr%d%n2) then
    box(1,2,1) = Glr%ns2
    box(2,2,1) = modulo(nboxp(2,2),Glr%d%n2+1)
    box(1,2,2) = nboxp(1,2)
    box(2,2,2) = Glr%ns2 + Glr%d%n2
  else
    box(1,2,1) = nboxp(1,2)
    box(2,2,1) = nboxp(2,2)
  end if

! Find limits in Z
  box(1,3,2) = -99
  box(2,3,2) = -99
  if(nboxp(1,3) < 0) then
    box(1,3,1) = Glr%ns3
    box(2,3,1) = nboxp(2,3)
    box(1,3,2) = modulo(nboxp(1,3),Glr%d%n3+1)
    box(2,3,2) = Glr%ns3 + Glr%d%n3    
  else if (nboxp(2,3) > Glr%d%n3) then
    box(1,3,1) = Glr%ns3
    box(2,3,1) = modulo(nboxp(2,3),Glr%d%n3+1)
    box(1,3,2) = nboxp(1,3)
    box(2,3,2) = Glr%ns3 + Glr%d%n3
  else
    box(1,3,1) = nboxp(1,3)
    box(2,3,2) = nboxp(2,3)
  end if

  !Limits of first zone (trivial)
  pbox(1,1,1) = box(1,1,1)
  pbox(2,1,1) = box(2,1,1)
  pbox(1,2,1) = box(1,2,1)
  pbox(2,2,1) = box(2,2,1)
  pbox(1,3,1) = box(1,3,1)
  pbox(2,3,1) = box(2,3,1)
  izone = 1

  !Limits of second zone (Translate X if you must)
  if (box(1,1,2) > 0) then
    izone = izone + 1
    pbox(1,1,izone) = box(1,1,2)
    pbox(2,1,izone) = box(2,1,2)
    pbox(1,2,izone) = box(1,2,1)
    pbox(2,2,izone) = box(2,2,1)
    pbox(1,3,izone) = box(1,3,1)
    pbox(2,3,izone) = box(2,3,1)
  end if

  !Limits of third zone (Translate Y if you must)
  if (box(1,2,2) > 0) then
    izone = izone + 1
    pbox(1,1,izone) = box(1,1,1)
    pbox(2,1,izone) = box(2,1,1)
    pbox(1,2,izone) = box(1,2,2)
    pbox(2,2,izone) = box(2,2,2)
    pbox(1,3,izone) = box(1,3,1)
    pbox(2,3,izone) = box(2,3,1)
  end if
  
  !Limits of fourth zone (Translate Z if you must)
  if (box(1,3,2) > 0) then
    izone = izone + 1
    pbox(1,1,izone) = box(1,1,1)
    pbox(2,1,izone) = box(2,1,1)
    pbox(1,2,izone) = box(1,2,1)
    pbox(2,2,izone) = box(2,2,1)
    pbox(1,3,izone) = box(1,3,2)
    pbox(2,3,izone) = box(2,3,2)
  end if

  !Limits of fifth zone (Translate X + Y if you must)
  if (box(1,1,2) > 0 .and. box(1,2,2) > 0) then
    izone = izone + 1
    pbox(1,1,izone) = box(1,1,2)
    pbox(2,1,izone) = box(2,1,2)
    pbox(1,2,izone) = box(1,2,2)
    pbox(2,2,izone) = box(2,2,2)
    pbox(1,3,izone) = box(1,3,1)
    pbox(2,3,izone) = box(2,3,1)
  end if

  !Limits of sixth zone (Translate X + Z if you must)
  if (box(1,1,2) > 0 .and. box(1,3,2) > 0) then
    izone = izone + 1
    pbox(1,1,izone) = box(1,1,2)
    pbox(2,1,izone) = box(2,1,2)
    pbox(1,2,izone) = box(1,2,1)
    pbox(2,2,izone) = box(2,2,1)
    pbox(1,3,izone) = box(1,3,2)
    pbox(2,3,izone) = box(2,3,2)
  end if 
  
  !Limits of seventh zone (Translate Y + Z if you must)
  if (box(1,2,2) > 0 .and. box(1,3,2) > 0) then
    izone = izone + 1
    pbox(1,1,izone) = box(1,1,1)
    pbox(2,1,izone) = box(2,1,1)
    pbox(1,2,izone) = box(1,2,2)
    pbox(2,2,izone) = box(2,2,2)
    pbox(1,3,izone) = box(1,3,2)
    pbox(2,3,izone) = box(2,3,2)
  end if

  !Limits of eigth zone (Translate X + Y + Z if you must)
  if (box(1,1,2) > 0 .and. box(1,2,2) > 0 .and. box(1,3,2) > 0) then
    izone = izone + 1
    pbox(1,1,izone) = box(1,1,2)
    pbox(2,1,izone) = box(2,1,2)
    pbox(1,2,izone) = box(1,2,2)
    pbox(2,2,izone) = box(2,2,2)
    pbox(1,3,izone) = box(1,3,2)
    pbox(2,3,izone) = box(2,3,2)
  end if

END SUBROUTINE fracture_projector
!%***

!#############################################################################################################################################
!!****f* BigDFT/number_of_projector_elements_in_locreg
!#############################################################################################################################################
!! FUNCTION: Calculates the number of segments (mseg) and elements (mvctr) of projectors in locreg
!!           
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine number_of_projector_elements_in_locreg(iatom,igrid,atoms,Glr,Llr,logrid,nlpspd,mproj,mseg,mvctr)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: iatom  ! current atom
  integer,intent(in) :: igrid  ! treat coarse (1) or fine (2) grid
  type(atoms_data),intent(in) :: atoms        ! atoms descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  integer,intent(in) :: mproj  ! number of projectors
  integer,intent(out):: mseg   ! number of segments
  integer,intent(out):: mvctr  ! number of elements
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  logical, dimension(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3), intent(in) :: logrid
  !#######################################
  ! Local Variables
  !#######################################
  integer :: i1,i2,i3  ! integers for loops
  integer :: nl1,nl2,nl3,nu1,nu2,nu3   ! rename the bounds of projectors coarse grid
  integer :: nend,nsrt,mvctri,nsrti,nendi
  integer,dimension(1:2,1:3) :: bound  ! rename the bounds of locreg
  logical :: plogrid   ! logical to check start of new segment


! Set boundaries of projectors (coarse)
   if(igrid == 1) then
      nl1 = nlpspd%nboxp_c(1,1,iatom)
      nl2 = nlpspd%nboxp_c(1,2,iatom)
      nl3 = nlpspd%nboxp_c(1,3,iatom)

      nu1 = nlpspd%nboxp_c(2,1,iatom)
      nu2 = nlpspd%nboxp_c(2,2,iatom)
      nu3 = nlpspd%nboxp_c(2,3,iatom)

   end if

! Set boundaries of projectors (fine)
   if(igrid == 2) then
      nl1 = nlpspd%nboxp_f(1,1,iatom)
      nl2 = nlpspd%nboxp_f(1,2,iatom)
      nl3 = nlpspd%nboxp_f(1,3,iatom)

      nu1 = nlpspd%nboxp_f(2,1,iatom)
      nu2 = nlpspd%nboxp_f(2,2,iatom)
      nu3 = nlpspd%nboxp_f(2,3,iatom)
   end if

! bounds of the localization region
  !lower bound
  bound(1,1) = Llr%ns1 - Glr%ns1
  bound(1,2) = Llr%ns2 - Glr%ns2
  bound(1,3) = Llr%ns3 - Glr%ns3

  !upper bound  (WILL NOT WORK FOR PERIODICITY)
  bound(2,1) = Llr%ns1 + Llr%d%n1 - Glr%ns1
  bound(2,2) = Llr%ns2 + Llr%d%n2 - Glr%ns2
  bound(2,3) = Llr%ns3 + Llr%d%n3 - Glr%ns3 

  print *, 'atom:',iatom
  print *, 'limits of proj:',nl1,nl2,nl3,nu1,nu2,nu3
  print *,'limits of overlap:',bound(1,:),bound(2,:)

!  do i3=1,Glr%d%n3
!    do i2=1,Glr%d%n2
!     do i1=1,Glr%d%n1
!        if(logrid(i1,i2,i3)) print *, 'logrid true:',i1,i2,i3
!     end do
!    end do
!  end do

  if (igrid == 1) then
!Initialize counters
     mvctr=0
     nsrt=0
     nend=0
     mvctri=0
     nsrti=0
     nendi=0

! Do coarse grid
     do i3=nl3,nu3
        if(i3 > bound(2,3) .or. i3 < bound(1,3)) cycle
        do i2=nl2,nu2
           if(i2 > bound(2,2) .or. i2 < bound(1,2)) cycle
           plogrid=.false.
           do i1=nl1,nu1
              if(i1 > bound(2,1) .or. i1 < bound(1,1))cycle

              if(logrid(i1,i2,i3)) then
                 mvctri=mvctri+1
                 if (.not. plogrid) then
                    nsrti=nsrti+1
                 endif
              else
                 if(plogrid) then
                    nendi=nendi+1
                 endif
              endif
              plogrid=logrid(i1,i2,i3)
           enddo
           if (i2 .le. bound(2,2) .and. i2 .ge. bound(1,2) .and. &
&              i3 .le. bound(2,3) .and. i3 .ge. bound(1,3) .and. &
               plogrid .eqv. .true.) then
              nendi=nendi+1
           endif
        enddo
     enddo

     mvctr=mvctr+mvctri
     nsrt=nsrt+nsrti
     nend=nend+nendi

     if (nend /= nsrt) then
        write(*,*)' ERROR in number_of_projector_elements_in_locreg : nend <> nsrt',nend,nsrt
        stop
     endif
     mseg=nend
  end if

  if(igrid == 2) then

     !Initialize counters
     mvctr=0
     nsrt=0
     nend=0
     mvctri=0
     nsrti=0
     nendi=0

! Do fine grid
     do i3=nl3,nu3
        if(i3 > bound(2,3) .or. i3 < bound(1,3)) cycle
        do i2=nl2,nu2
           if(i2 > bound(2,2) .or. i2 < bound(1,2)) cycle
           plogrid=.false.
           do i1=nl1,nu1
              if(i1 > bound(2,1) .or. i1 < bound(1,1))cycle

              if(logrid(i1,i2,i3)) then
                 mvctri=mvctri+1
                 if (.not. plogrid) then
                    nsrti=nsrti+1
                 endif
              else
                 if(plogrid) then
                    nendi=nendi+1
                 endif
              endif
              plogrid=logrid(i1,i2,i3)
           enddo
           if (i2 .le. bound(2,2) .and. i2 .ge. bound(1,2) .and. &
&              i3 .le. bound(2,3) .and. i3 .ge. bound(1,3) .and. &
               plogrid .eqv. .true.) then
              nendi=nendi+1
           endif
        enddo
     enddo

     mvctr=mvctr+mvctri
     nsrt=nsrt+nsrti
     nend=nend+nendi

     if (nend /= nsrt) then
        write(*,*)' ERROR in number_of_projector_elements_in_locreg (fine) : nend <> nsrt',nend,nsrt
        stop
     endif
     mseg=nend
  end if

END SUBROUTINE number_of_projector_elements_in_locreg
!%***


!#############################################################################################################################################
!!****f* BigDFT/ projector_box_in_locreg
!#############################################################################################################################################
!! FUNCTION: Calculates the bounds of the box of the projector in locreg
!!           bounds(1,:,:) for coarse grid
!!           bounds(2,:,:) for fine grid
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine projector_box_in_locreg(iatom,Glr,Llr,nlpspd,bounds)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: iatom  ! current atom
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  integer,dimension(1:2,1:2,1:3),intent(out) :: bounds
  !#######################################
  ! Local Variables
  !#######################################
  integer :: ii
  integer,dimension(1:2,1:3) :: Cnl,Fnl,Lnl
  
! Set boundaries of projectors (coarse)
  !lower bounds
  Cnl(1,1) = nlpspd%nboxp_c(1,1,iatom)
  Cnl(1,2) = nlpspd%nboxp_c(1,2,iatom)
  Cnl(1,3) = nlpspd%nboxp_c(1,3,iatom)
  !upper bounds
  Cnl(2,1) = nlpspd%nboxp_c(2,1,iatom)
  Cnl(2,2) = nlpspd%nboxp_c(2,2,iatom)
  Cnl(2,3) = nlpspd%nboxp_c(2,3,iatom)

! Set boundaries of projectors (fine)
  !lower bounds
  Fnl(1,1) = nlpspd%nboxp_f(1,1,iatom)
  Fnl(1,2) = nlpspd%nboxp_f(1,2,iatom)
  Fnl(1,3) = nlpspd%nboxp_f(1,3,iatom)
  !upper bounds
  Fnl(2,1) = nlpspd%nboxp_f(2,1,iatom)
  Fnl(2,2) = nlpspd%nboxp_f(2,2,iatom)
  Fnl(2,3) = nlpspd%nboxp_f(2,3,iatom)

! bounds of the localization region
  !lower bounds
  Lnl(1,1) = Llr%ns1 - Glr%ns1
  Lnl(1,2) = Llr%ns2 - Glr%ns2
  Lnl(1,3) = Llr%ns3 - Glr%ns3
  !upper bounds
  Lnl(2,1) = Llr%ns1 + Llr%d%n1 - Glr%ns1
  Lnl(2,2) = Llr%ns2 + Llr%d%n2 - Glr%ns2
  Lnl(2,3) = Llr%ns3 + Llr%d%n3 - Glr%ns3

! Calculate the bounds for the grids
  do ii=1,3
     !lower bounds
     bounds(1,1,ii) = max(Lnl(1,ii),Cnl(1,ii))
     bounds(2,1,ii) = max(Lnl(1,ii),Fnl(1,ii))
     ! upper bounds
     bounds(1,2,ii) = min(Lnl(2,ii),Cnl(2,ii))
     bounds(2,2,ii) = min(Lnl(2,ii),Fnl(2,ii))
  end do

END SUBROUTINE projector_box_in_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/allocate_Lnlpspd
!#############################################################################################################################################
!! FUNCTION:  Allocates most of the arrays in Lnlpspd 
!!
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine allocate_Lnlpspd(natom,Lnlpspd,subname)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: natom
  type(nonlocal_psp_descriptors),intent(inout) :: Lnlpspd  ! Local descriptors for the projectors
  character(len=*), intent(in) :: subname
  !#######################################
  ! Local Variables 
  !#######################################
  integer :: i_stat

  allocate(Lnlpspd%nvctr_p(2*natom+ndebug),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nvctr_p,'nvctr_p',subname)
  allocate(Lnlpspd%nseg_p(2*natom+ndebug),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nseg_p,'nseg_p',subname)
  allocate(Lnlpspd%nboxp_c(2,3,2*natom),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nboxp_c,'nbox_c',subname)
  allocate(Lnlpspd%nboxp_f(2,3,2*natom),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nboxp_f,'nbox_f',subname)

END SUBROUTINE allocate_Lnlpspd
!%***


!#############################################################################################################################################
!!****f* BigDFT/allocate_projd
!#############################################################################################################################################
!! FUNCTION: allocates the keyg_p and keyv_p descriptors for the projectors
!!          
!!           
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine allocate_projd(mseg,Lnlpspd,subname)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: mseg
  type(nonlocal_psp_descriptors),intent(inout) :: Lnlpspd  ! Local descriptors for the projectors
  character(len=*), intent(in) :: subname
  !#######################################
  ! Local Variables 
  !#######################################
  integer :: i_stat

  allocate(Lnlpspd%keyg_p(2,mseg),stat=i_stat)
  call memocc(i_stat,Lnlpspd%keyg_p,'keyg_p',subname)
  allocate(Lnlpspd%keyv_p(mseg),stat=i_stat)
  call memocc(i_stat,Lnlpspd%keyv_p,'keyv_p',subname)

END SUBROUTINE allocate_projd
!%***

!#############################################################################################################################################
!!****f* BigDFT/apply_local_projectors
!#############################################################################################################################################
!! FUNCTION: Fills the projector pointer and applies the projectors to the wavefunctions
!!           
!!           
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine apply_local_projectors(atoms,in,Llr,Lnlpspd,Lproj,orbs,projflg,psi,rxyz,hpsi)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  type(atoms_data),intent(in) :: atoms
  type(input_variables),intent(in) :: in
  type(locreg_descriptors),intent(in) :: Llr
  type(nonlocal_psp_descriptors),intent(in) :: Lnlpspd  ! Local descriptors for the projectors
  type(orbitals_data),intent(in) :: orbs
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  integer,dimension(atoms%nat),intent(in) :: projflg
  real(wp),dimension(Lnlpspd%nprojel),intent(out):: Lproj  !local projectors
  real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),intent(in) :: psi  !local wavefunction
  real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),intent(out):: hpsi ! local |p><p|Psi>
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  !#######################################
  ! Local Variables 
  !#######################################
  integer :: ikpt,istart_c,ncplx,jseg_c,iproj,iat,ityp,l,i,nwarnings
  integer :: isorb,ieorb,nspinor,iorb,istart_o,ispinor
  integer :: nels,ipsi,ii,iatom
  real(gp) :: kx,ky,kz,eproj_spinor
  real(wp),dimension(orbs%norbp,(Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f),orbs%nspinor) :: psi_tmp
  real(wp),dimension(orbs%norbp,(Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f),orbs%nspinor) :: hpsi_tmp

!  First reshape the wavefunctions: psi_tmp(norb,nels,nspinor)
   nels = Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f
   print *,'nels:',nels, Llr%wfd%nvctr_c,Llr%wfd%nvctr_f
   psi_tmp = reshape(psi, (/ orbs%norbp, nels, orbs%nspinor /),order=(/ 2, 3, 1 /))
   hpsi_tmp = reshape(hpsi,(/ orbs%norbp, nels, orbs%nspinor /),order=(/ 2, 3, 1 /))


     ikpt=orbs%iokpt(1)
     loop_kpt: do
      
        !features of the k-point ikpt
        kx=orbs%kpts(1,ikpt)
        ky=orbs%kpts(2,ikpt)
        kz=orbs%kpts(3,ikpt)

        !evaluate the complexity of the k-point
        if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
           ncplx=1
        else
           ncplx=2
        end if

        jseg_c = 1
        iproj = 0
        iatom = 0
        do iat = 1,atoms%nat
           if(projflg(iat) == 0) cycle
           iatom = iatom +1
           istart_c = 1
           ityp=atoms%iatype(iat)

           do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
              do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
                 if (atoms%psppar(l,i,ityp) /= 0.0_gp) then

!                   Second fill the projectors
!                   NOTE : idir was set to 0 because we don't care for derivatives
                    call local_projector(atoms%geocode,atoms%atomnames(ityp),iat,0,l,i,&
                         atoms%psppar(l,0,ityp),rxyz(1,iat),Llr,&
                         in%hx,in%hy,in%hz,kx,ky,kz,ncplx,Lnlpspd%nvctr_p(2*iat-1),&
                         Lnlpspd%nvctr_p(2*iat),Lnlpspd%nseg_p(2*iat-1),Lnlpspd%nseg_p(2*iat),&
                         Lnlpspd%keyv_p(jseg_c),Lnlpspd%keyg_p(1,jseg_c),Lproj(istart_c),nwarnings)

                    iproj=iproj+2*l-1
                    istart_c=istart_c+(Lnlpspd%nvctr_p(2*iatom-1)+7*Lnlpspd%nvctr_p(2*iatom))*(2*l-1)*ncplx
                    !print *,'iproc,istart_c,nlpspd%nprojel',istart_c,Lnlpspd%nprojel,ncplx,nlpspd%nprojel
                    if (istart_c > Lnlpspd%nprojel+1) stop 'istart_c > nprojel+1'
                    if (iproj > Lnlpspd%nproj) stop 'iproj > nproj'
                 endif
              enddo
           enddo


!          Apply them on the wavefunctions in the overlap region
!          hpsi contains the new wavefunctions
           call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor) 

           do iorb=isorb,ieorb 
              istart_o=1
              do ispinor=1,nspinor,ncplx
                 if (ispinor >= 2) istart_o=1

                 !GTH and HGH pseudopotentials
                 do l=1,4
                    do i=1,3
                       if (atoms%psppar(l,i,ityp) /= 0.0_gp) then
                          call applyprojector(ncplx,l,i,atoms%psppar(0,0,ityp),atoms%npspcode(ityp),&
                               Llr%wfd%nvctr_c,Llr%wfd%nvctr_f,Llr%wfd%nseg_c,&
                               Llr%wfd%nseg_f,Llr%wfd%keyv,Llr%wfd%keyg,&
                               Lnlpspd%nvctr_p(2*iatom-1),Lnlpspd%nvctr_p(2*iatom),Lnlpspd%nseg_p(2*iatom-1),&
                               Lnlpspd%nseg_p(2*iatom),Lnlpspd%keyv_p(jseg_c),Lnlpspd%keyg_p(1,jseg_c),&
                               Lproj(istart_o),psi_tmp(iorb,:,ispinor),hpsi_tmp(iorb,:,ispinor),eproj_spinor)
                           
                           istart_o=istart_o+(Lnlpspd%nvctr_p(2*iatom-1)+7*Lnlpspd%nvctr_p(2*iatom))*(2*l-1)*ncplx
                       end if
                    enddo
                 enddo
              end do
           end do
           jseg_c = jseg_c + Lnlpspd%nseg_p(2*iatom - 1)+ Lnlpspd%nseg_p(2*iatom) 
        end do  !on iat
        
        ipsi = 0
        do iorb = isorb,ieorb
           do ispinor=1,nspinor,ncplx
              do ii = 1,Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f
                 hpsi(ipsi+ii) = hpsi_tmp(iorb,ii,ispinor)
              end do
              ipsi = ipsi + Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f
           end do
        end do
        
        if (iproj /= Lnlpspd%nproj) stop 'incorrect number of projectors created'
        if (ieorb == orbs%norbp) exit loop_kpt
        ikpt=ikpt+1
     end do loop_kpt

END SUBROUTINE apply_local_projectors
!%***

!> BigDFT/projector
!!
!!
subroutine local_projector(geocode,atomname,iat,idir,l,i,gau_a,rxyz,Llr,&
     hx,hy,hz,kx,ky,kz,ncplx,&
     mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj,nwarnings)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=20), intent(in) :: atomname
  type(locreg_descriptors),intent(in) :: Llr
  integer, intent(in) :: iat,idir,l,i,mbvctr_c,mbvctr_f,mseg_c,mseg_f,ncplx
  real(gp), intent(in) :: hx,hy,hz,gau_a,kx,ky,kz
  !integer, dimension(2,3), intent(in) :: nboxp_c,nboxp_f
  integer, dimension(mseg_c+mseg_f), intent(in) :: keyv_p
  integer, dimension(2,mseg_c+mseg_f), intent(in) :: keyg_p
  real(gp), dimension(3), intent(in) :: rxyz
  integer, intent(inout) :: nwarnings
  real(wp), dimension((mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx), intent(out) :: proj
  !local variables
  integer, parameter :: nterm_max=20 !if GTH nterm_max=4
  integer :: m,iterm
  !integer :: nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f
  integer :: istart_c,nterm
  real(gp) :: fpi,factor,rx,ry,rz
  real(dp) :: scpr
  integer, dimension(3) :: nterm_arr
  integer, dimension(nterm_max) :: lx,ly,lz
  integer, dimension(3,nterm_max,3) :: lxyz_arr
  real(gp), dimension(nterm_max) :: factors
  real(gp), dimension(nterm_max,3) :: fac_arr

  !this value can also be inserted as a parameter
  fpi=(4.0_gp*atan(1.0_gp))**(-.75_gp)

  rx=rxyz(1)
  ry=rxyz(2)
  rz=rxyz(3)
  
  print *, 'atom:',iat
  print *,'dim of proj',(mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx,mbvctr_c,mbvctr_f
  istart_c=1
  !start of the projectors expansion routine
  factor=sqrt(2.0_gp)*fpi/(sqrt(gau_a)**(2*(l-1)+4*i-1))
  do m=1,2*l-1

     if (idir==0) then !normal projector calculation case
        call calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,factors)

        factors(1:nterm)=factor*factors(1:nterm)
     else !calculation of projector derivative
        call calc_coeff_derproj(l,i,m,nterm_max,gau_a,nterm_arr,lxyz_arr,fac_arr)

        nterm=nterm_arr(idir)
        do iterm=1,nterm
           factors(iterm)=factor*fac_arr(iterm,idir)
           lx(iterm)=lxyz_arr(1,iterm,idir)
           ly(iterm)=lxyz_arr(2,iterm,idir)
           lz(iterm)=lxyz_arr(3,iterm,idir)
        end do
     end if
     
     call crtproj(geocode,nterm,Llr,hx,hy,hz,kx,ky,kz,ncplx,&
          gau_a,factors,rx,ry,rz,lx,ly,lz,&
          mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj(istart_c))

     ! testing
     if (idir == 0) then
        !here the norm should be done with the complex components
        call wnrm_wrap(ncplx,mbvctr_c,mbvctr_f,proj(istart_c),scpr)
        if (abs(1.d0-scpr) > 1.d-2) then
           if (abs(1.d0-scpr) > 1.d-1) then
              !if (iproc == 0) then
                write(*,'(1x,a,i4,a,a6,a,i1,a,i1,a,f6.3)')&
                      'The norm of the nonlocal PSP for atom n=',iat,&
                      ' (',trim(atomname),&
                      ') labeled by l=',l,' m=',m,' is ',scpr
                 write(*,'(1x,a)')&
                      'while it is supposed to be about 1.0.'
              !end if
           else
              nwarnings=nwarnings+1
           end if
        end if
     end if
     !end testing
     istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
  enddo
END SUBROUTINE local_projector

!determine a set of localisation regions from the centers and the radii.
!cut in cubes the global reference system
subroutine determine_locreg2(nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,outofzone)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  character(len=*), parameter :: subname='determine_locreg'
  logical :: perx,pery,perz
  integer :: ilr,isx,isy,isz,iex,iey,iez
  integer :: ln1,ln2,ln3
  integer :: ii !tests
  real(gp) :: rx,ry,rz,cutoff  

  write(*,*)'Inside determine_locreg2:'

  !initialize out of zone
  outofzone (:,:) = 0     
  
  !determine the limits of the different localisation regions
  do ilr=1,nlr

     rx=cxyz(1,ilr)
     ry=cxyz(2,ilr)
     rz=cxyz(3,ilr)

     cutoff=locrad(ilr)

     isx=floor((rx-cutoff)/hx)
     isy=floor((ry-cutoff)/hy)
     isz=floor((rz-cutoff)/hz)

     iex=ceiling((rx+cutoff)/hx)
     iey=ceiling((ry+cutoff)/hy)
     iez=ceiling((rz+cutoff)/hz)

     ln1 = iex-isx
     ln2 = iey-isy
     ln3 = iez-isz

     ! First check if localization region fits inside box
     if (iex - isx >= Glr%d%n1 - 14) then
        write(*,*)'Width of direction x :',(iex - isx)*hx,' of localization region:',ilr
        write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n1*hx,'.'
        write(*,*)'Increasing the simulation box is recommended. The code will use the box limits'
        write(*,*)'of the simulation box.'
     end if
     if (iey - isy >= Glr%d%n2 - 14) then
        write(*,*)'Width of direction z :',(iey - isy)*hy,' of localization region:',ilr
        write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n2*hy,'.'
        write(*,*)'Increasing the simulation box is recommended. The code will use the width'
        write(*,*)'of the simulation box.'
     end if
     if (iez - isz >= Glr%d%n3 - 14) then
        write(*,*)'Width of direction z :',(iez - isz)*hz,' of localization region:',ilr
        write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n3*hz,'.'
        write(*,*)'Increasing the simulation box is recommended. The code will use the width'
        write(*,*)'of the simulation box.'
     end if 

     ! Localization regions should always have free boundary conditions
     Llr(ilr)%geocode='F'

     !assign the starting/ending points and outofzone for the different
     ! geometries
     select case(Glr%geocode)
     case('F')
        isx=max(isx,Glr%ns1)
        isy=max(isy,Glr%ns2)
        isz=max(isz,Glr%ns3)

        iex=min(iex,Glr%ns1+Glr%d%n1)
        iey=min(iey,Glr%ns2+Glr%d%n2)
        iez=min(iez,Glr%ns3+Glr%d%n3)

     case('S')
        ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1,ilr)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        isy=max(isy,Glr%ns2)
        iey=min(iey,Glr%ns2 + Glr%d%n2)
        outofzone(2,ilr) = 0

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3,ilr)=modulo(iez,Glr%d%n3+1)
           end if 
        end if

     case('P')
         ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1,ilr)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        if (iey - isy >= Glr%d%n2) then       
           isy=Glr%ns2
           iey=Glr%ns2 + Glr%d%n2
         else
           isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
           iey= ln2 + isy
           if (iey > Glr%ns2+Glr%d%n2) then
              outofzone(2,ilr)=modulo(iey,Glr%d%n2+1)
           end if           
        end if

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3,ilr)=modulo(iez,Glr%d%n3+1)
           end if 
        end if
     end select

     !values for the starting point of the cube
     Llr(ilr)%ns1=isx
     Llr(ilr)%ns2=isy
     Llr(ilr)%ns3=isz
     !dimensions of the localisation region
     Llr(ilr)%d%n1=iex-isx
     Llr(ilr)%d%n2=iey-isy
     Llr(ilr)%d%n3=iez-isz

     !dimensions of the fine grid inside the localisation region
     Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
     Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
     Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
     
     !NOTE: This will not work with symmetries (must change it)
     Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
     Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
     Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz

     !dimensions of the interpolating scaling functions grid (reduce to +2?, check with Luigi)
     Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
     Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
     Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31

!DEBUG
     write(*,*)'Description of zone:',ilr
     write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
     write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
     write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
     write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
     write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
     write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
     write(*,*)'outofzone',ilr,':',outofzone(:,ilr)
!DEBUG

    ! construct the wavefunction descriptors (wfd)
     call determine_wfd_periodicity(ilr,nlr,Glr,Llr,outofzone)
     
!     print *,'Before locreg_bounds'
!     print *,'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!     print *,'nl,nu:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3
!     print *,'wfd(nseg):',Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nseg_f
!     print *,'wfd(nvctr):',Llr(ilr)%wfd%nvctr_c,Llr(ilr)%wfd%nvctr_f

     ! Sould check if nfu works properly... also relative to locreg!!
     !if the localisation region is isolated build also the bounds
     if (Llr(ilr)%geocode=='F') then
        call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
             Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
             Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
     end if
     print *,'Outside locreg_bounds'
  end do !on ilr

  !after all localisation regions are determined draw them
  !call draw_locregs(nlr,hx,hy,hz,Llr)

END SUBROUTINE determine_locreg2


!#############################################################################################################################################
!!****f* BigDFT/overlap_region
!#############################################################################################################################################
!! FUNCTION: Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg taking into account the pediodicity
!!          
!! WARNING: We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
!!         
!! SOURCE:
!!
subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr,outofzone)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  !#############################################
  !local variables
  !############################################
  integer :: ii
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period
  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements
  character(len=*), parameter :: subname='determine_wfd_periodicity'

   !starting point of locreg (always inside locreg)
   isdir(1) = Llr(ilr)%ns1
   isdir(2) = Llr(ilr)%ns2
   isdir(3) = Llr(ilr)%ns3
   !ending point of locreg (can be outside the simulation box)
   iedir(1) = Llr(ilr)%ns1 + Llr(ilr)%d%n1
   iedir(2) = Llr(ilr)%ns2 + Llr(ilr)%d%n2
   iedir(3) = Llr(ilr)%ns3 + Llr(ilr)%d%n3
   ! starting and ending point of fine grid in Global region
   Gifs(1) = Glr%d%nfl1 + Glr%ns1
   Gifs(2) = Glr%d%nfl2 + Glr%ns2
   Gifs(3) = Glr%d%nfl3 + Glr%ns3
   Gife(1) = Glr%d%nfu1 + Glr%ns1
   Gife(2) = Glr%d%nfu2 + Glr%ns2
   Gife(3) = Glr%d%nfu3 + Glr%ns3
   ! periodicity
   period(1) = Glr%d%n1
   period(2) = Glr%d%n2
   period(3) = Glr%d%n3

   ! Determine starting point of the fine grid in locreg
   do ii=1,3
      if (outofzone(ii,ilr) > 0) then
         ! When periodicity, we must check for 2 different situations:
         ! (1) : starting of locreg before or in fine grid zone
         if (isdir(ii) < Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
         ! (2) : starting point after fine grid
         if (isdir(ii) > Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii)+period(ii))-isdir(ii)
      else
          Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
      end if 
   end do

   ! Determine ending point of the fine grid in locreg
   do ii=1,3
      if(outofzone(ii,ilr) > 0) then
         !When periodicity, we must check for three different situations:
         ! (1) : ending of locreg before fine grid zone
         if(iedir(ii) < (Gifs(ii) + period(ii))) Life(ii) = Gife(ii)-isdir(ii)
         ! (2) : ending of locreg in fine grid zone
         if(iedir(ii) > (Gifs(ii) + period(ii)) .and. iedir(ii) < (Gife(ii) + period(ii))) then
           Life(ii) = iedir(ii)-isdir(ii)
         end if
         ! (3) : ending of locreg after ending of fine grid zone
         if(iedir(ii) > (Gife(ii)+period(ii))) Life(ii) = Gife(ii) + period(ii)-isdir(ii)
      else
         Life(ii) = min(iedir(ii),Gife(ii))-isdir(ii)
      end if
   end do

   ! Assign values to Llr
   Llr(ilr)%d%nfl1 = Lifs(1)
   Llr(ilr)%d%nfl2 = Lifs(2)
   Llr(ilr)%d%nfl3 = Lifs(3)
   Llr(ilr)%d%nfu1 = Life(1)
   Llr(ilr)%d%nfu2 = Life(2)
   Llr(ilr)%d%nfu3 = Life(3)

   ! define the wavefunction descriptors inside the localisation region
   !coarse part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_c,Glr%wfd%nvctr_c,&
          Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),nseg_c,nvctr_c,outofzone(:,ilr))
   !fine part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),nseg_f,nvctr_f,outofzone(:,ilr))

   ! Assign the values to Llr
   Llr(ilr)%wfd%nseg_c = nseg_c
   Llr(ilr)%wfd%nseg_f = nseg_f
   Llr(ilr)%wfd%nvctr_c= nvctr_c
   Llr(ilr)%wfd%nvctr_f= nvctr_f

   !allocate the wavefunction descriptors following the needs
   call allocate_wfd(Llr(ilr)%wfd,subname)

   !Now, fill the descriptors:
   !coarse part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),&
        Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
        Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
        Llr(ilr)%wfd%keyg(1,1),Llr(ilr)%wfd%keyv(1),outofzone(:,ilr))

   !fine part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
        Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
        Llr(ilr)%wfd%keyg(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),outofzone(:,ilr))

END SUBROUTINE determine_wfd_periodicity


!#############################################################################################################################################
!!****f* BigDFT/num_segkeys_periodic
!#############################################################################################################################################
!! FUNCTION: Calculates the number of segments and elements in localisation region
!!          
!! WARNING:   
!!         
!! SOURCE:
!!
subroutine num_segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,outofzone)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, intent(out) :: nseg_loc,nvctr_loc
  integer, dimension(3),intent(in) :: outofzone
  !local variables
  logical :: lseg,go1,go2,go3
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check

  nvctr_loc=0
  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0

  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     lseg=.false.
     ! overlap conditions if zone completely inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec)
     ! overlap conditions if zone as components in other periodic cells
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        nvctr_check=nvctr_check+1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)

        if (go1 .and. go2 .and. go3 ) then
           nvctr_loc=nvctr_loc+1
           if (.not. lseg) then
              nsrt=nsrt+1
           end if
           lseg=.true.
        else
           if (lseg) then
              nend=nend+1
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
     end if
  end do
  nseg_loc=nend

  !check
  if (nend /= nsrt) then
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif

  if (nvctr_check /= nvctr) then
     write(*,'(1x,a,2(i6))')&
          'ERROR: incorrect number of coarse points examined for reducing the localisation region',&
          nvctr_check,nvctr
     stop
  end if

END SUBROUTINE num_segkeys_periodic

subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keyg_loc,keyv_loc,outofzone)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(3), intent(in) :: outofzone
  integer, dimension(nseg_loc), intent(out) :: keyv_loc
  integer, dimension(2,nseg_loc), intent(out) :: keyg_loc
  !local variables
  logical :: go1,go2,go3,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l
  integer :: ngridp

  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0
  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     lseg=.false.

     ! intersection condition if zone inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec) 
     ! intersection condition if zone has components outside simulation box (periodic)
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)
        if (go1 .and. go2 .and. go3) then
          !index of the compressed function
          i1l=i-i1sc
          if(outofzone(1) > 0 .and. i <= outofzone(1))i1l = i - i1sc + n1 + 1  
          i2l=i2-i2sc
          if(outofzone(2) > 0 .and. i2 <= outofzone(2))i2l = i2 - i2sc + n2 + 1
          i3l=i3-i3sc
          if(outofzone(3) > 0 .and. i3 <= outofzone(3))i3l = i3 - i3sc + n3 + 1
          ngridp=i3l*((n1l+1)*(n2l+1)) + i2l*(n1l+1) + i1l+1

          nvctr_check=nvctr_check+1
          if (.not. lseg) then
!             print *,'         check:',i,i2,i3,i1l,i2l,i3l,ngridp
             nsrt=nsrt+1
             keyg_loc(1,nsrt)=ngridp
             keyv_loc(nsrt)=nvctr_check
          end if
          lseg=.true.
        else 
           if (lseg) then
!              print *,'in        else:',i,i2,i3,i1l,i2l,i3l,ngridp
              nend=nend+1
              keyg_loc(2,nend)=ngridp
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
!        print *,'in second else:',i,i2,i3,i1l,i2l,i3l,ngridp
        nend=nend+1
        keyg_loc(2,nend)=ngridp
     end if
  end do

  !check
  if (nvctr_check /= nvctr_loc .or. nend /= nsrt .or. nend /= nseg_loc) then
     print *,'global region statistics:',nseg,nvctr
     write(*,*)&
          'ERROR: problem in segkeys_periodic  ',&
          'nvctr_check:',nvctr_check,'nvctr_loc:',nvctr_loc,&
          'nend:',nend,'nsrt:',nsrt,'nseg_loc:',nseg_loc
     stop
  end if

END SUBROUTINE segkeys_periodic

!#############################################################################################################################################
!!****f* BigDFT/overlap_region
!#############################################################################################################################################
!! FUNCTION: Determines the number of intersection regions between locregs, taking into account the periodicity of the system.
!!          
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr,outofzone)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(out) :: isovrlp             ! Integer giving the number of overlaps (max 8 with periodicity)
  !########################################
  !Subroutine Array Arguments
  !########################################
  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  !#############################################
  !local variables
  !############################################
  integer :: ii,azones,bzones,i_stat
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  character(len=*), parameter :: subname='get_number_of_overlap_region'
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(outofzone(ii,alr) > 0) azones = azones * 2
     if(outofzone(ii,blr) > 0) bzones = bzones * 2
  end do

write(*,*)'azones,bzones',azones,bzones
write(*,*)'outofzone',alr,':',outofzone(:,alr)
write(*,*)'outofzone',blr,':',outofzone(:,blr)

!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),outofzone(:,alr),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),outofzone(:,blr),bstart,bend)

! Now check the number of overlapping zones
 isovrlp = 0
 do izones=1,azones
   do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
      if(go1 .and. go2 .and. go3) then
         isovrlp = isovrlp + 1
      end if
   end do
 end do
  
END SUBROUTINE get_number_of_overlap_region

!#############################################################################################################################################
!!****f* BigDFT/fracture_periodic_zone
!#############################################################################################################################################
!! FUNCTION: Divides the locreg into zones contained inside the simulation box, by applying the primitive vectors
!!           It returns: astart(3,nzones) which is the starting points of the different zones (max. 8)
!!                       aend(3,nzones) which is the ending points of the different zones (max. 8)
!!          
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine fracture_periodic_zone(nzones,Glr,Llr,outofzone,astart,aend)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: nzones
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
  integer,dimension(3),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  integer,dimension(3,nzones),intent(out) :: astart !
  integer,dimension(3,nzones),intent(out) :: aend !
  !#############################################
  !local variables
  !############################################
  integer :: ii,index,jj
  integer,dimension(3) :: alrs,alre,Gend,Gstart,period
  character(len=*), parameter :: subname='fracture_periodic_zone'
  
! Start and end of Global region
  Gstart(1) = Glr%ns1 
  Gstart(2) = Glr%ns2
  Gstart(3) = Glr%ns3  
  Gend(1) = Glr%ns1 + Glr%d%n1
  Gend(2) = Glr%ns2 + Glr%d%n2
  Gend(3) = Glr%ns3 + Glr%d%n3

! Periodicity of the system
  period(1) = Glr%d%n1 + 1
  period(2) = Glr%d%n2 + 1
  period(3) = Glr%d%n3 + 1

! Start and end of local region
  alrs(1) = Llr%ns1
  alrs(2) = Llr%ns2
  alrs(3) = Llr%ns3
  alre(1) = Llr%ns1 + Llr%d%n1
  alre(2) = Llr%ns2 + Llr%d%n2
  alre(3) = Llr%ns3 + Llr%d%n3

!assign the first zone (necessarily without shift) and initialize the rest
  do ii=1,3
     astart(ii,:) = alrs(ii)
     aend(ii,:) = min(Gend(ii),alre(ii))
  end do

!assign the other zones
  index = 2
  do ii=1,3
     if(outofzone(ii) > 0) then    !Translation: X,Y,Z
        astart(ii,index) =  Gstart(ii)
        aend(ii,index) = modulo(alre(ii),period(ii))
        index = index + 1
     end if 
     do jj=ii+1,3
        if(outofzone(ii) > 0 .and. outofzone(jj) > 0) then  !Translation: X+Y,X+Z,Y+Z
           astart(ii,index) = Gstart(ii)
           astart(jj,index) = Gstart(jj)
           aend(ii,index) = modulo(alre(ii),period(ii))
           aend(jj,index) = modulo(alre(jj),period(jj))
           index = index + 1
        end if
     end do
  end do

  if(outofzone(1) > 0 .and. outofzone(2) > 0 .and. outofzone(3) > 0 ) then ! Translation: X+Y+Z
     astart(1,index) = Gstart(1)
     astart(2,index) = Gstart(2)
     astart(3,index) = Gstart(3)
     aend(1,index) = modulo(alre(1),period(1))
     aend(2,index) = modulo(alre(2),period(2))
     aend(3,index) = modulo(alre(3),period(3))
  end if

END SUBROUTINE fracture_periodic_zone


!##############################################################################################################################################
!!****f* BigDFT/get_overlap_region
!##############################################################################################################################################
!! FUNCTION Given two localization regions, A and B, this routine returns a localization region corresponding to the intersection of A & B. 
!!
!! SOURCE
!!
subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr,outofzone)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(in) :: isovrlp              ! Number of overlap regions
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  integer,dimension(3,nlr),intent(in) :: outofzone 
  !#############################################
  !local variables
  !############################################
  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
  character(len=*), parameter :: subname='get_overlap_region_periodic'
  !# NEW
  integer :: ii,azones,bzones,i_stat,index
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(outofzone(ii,alr) > 0) azones = azones * 2
     if(outofzone(ii,blr) > 0) bzones = bzones * 2
  end do

!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),outofzone(:,alr),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),outofzone(:,blr),bstart,bend)

! Now check the number of overlapping zones
  index = 0
  do izones=1,azones
    do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
      if(go1 .and. go2 .and. go3) then
        index = index + 1

! Now construct the Overlap localization region descriptor
! only if there is an overlap. The following only works
! when the previous test is successful. Note also that
! isx, isy and isz are necessarily in the Glr by construction
! of the Llrs, so don't need to test them.
         
        ! Determine the limits of the overlap region
        isx = max(astart(1,izones),bstart(1,jzones))
        isy = max(astart(2,izones),bstart(2,jzones))
        isz = max(astart(3,izones),bstart(3,jzones))

        iex = min(aend(1,izones),bend(1,jzones))
        iey = min(aend(2,izones),bend(2,jzones))
        iez = min(aend(3,izones),bend(3,jzones))

!       Checks to assign the geometric code of the overlap region (TO DO,could be interesting for Pascal?)
!       This could change the values of the bounds, so do it here
!       for now, in sandbox,put free boundary to all zones
        Olr(index)%geocode = 'F'  

!       Values for the starting point of the cube
        Olr(index)%ns1 = isx
        Olr(index)%ns2 = isy
        Olr(index)%ns3 = isz

!       Dimensions of the overlap region
        Olr(index)%d%n1 = iex - isx 
        Olr(index)%d%n2 = iey - isy 
        Olr(index)%d%n3 = iez - isz 
    
!       Dimensions of the fine grid inside the overlap region
        if (isx < iex) then
           Olr(index)%d%nfl1=max(isx,Glr%d%nfl1)-isx
           Olr(index)%d%nfu1=min(iex,Glr%d%nfu1)-isx
        else
           write(*,*)'Yet to be implemented (little effort?)'
           stop
        end if

        if (isy < iey) then
           Olr(index)%d%nfl2=max(isy,Glr%d%nfl2)-isy
           Olr(index)%d%nfu2=min(iey,Glr%d%nfu2)-isy
        else
           write(*,*)'Yet to be implemented (little effort?)'
           stop
        end if

        if (isz < iez) then
           Olr(index)%d%nfl3=max(isz,Glr%d%nfl3)-isz
           Olr(index)%d%nfu3=min(iez,Glr%d%nfu3)-isz
        else
           write(*,*)'Yet to be implemented (little effort?)'
           stop
        end if

!       Dimensions of the interpolating scaling function grid 
!       (geocode already taken into acount because it is simple)
        select case(Olr(index)%geocode)
        case('F')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+31
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+31
        case('S')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        case('P')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+2
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        end select
 
!       Now define the wavefunction descriptors inside the overlap region
!       First calculate the number of points and segments for the region
!       Coarse part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!       Fine part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
          Olr(index)%wfd%keyg(1,1),Olr(index)%wfd%keyv(1))
!       Fine part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
          Olr(index)%wfd%keyg(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
          Olr(index)%wfd%keyv(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))

!       If the localisation region is isolated build also the bounds
        if (Olr(index)%geocode=='F') then
           call locreg_bounds(Olr(index)%d%n1,Olr(index)%d%n2,Olr(index)%d%n3,&
             Olr(index)%d%nfl1,Olr(index)%d%nfu1,Olr(index)%d%nfl2,Olr(index)%d%nfu2,&
             Olr(index)%d%nfl3,Olr(index)%d%nfu3,Olr(index)%wfd,Olr(index)%bounds)
     
        end if
     end if ! go1 .and. go2 .and. go3
   end do !jzones
 end do !izones

! Check on the number of zones
  if (index /= isovrlp) then
      write(*,*)&
          'ERROR: problem in get_overlap_region_periodic ',&
          'index:',index,'not equal to isovrlp:',isovrlp,&
          'The number of overlap descriptors constructed does not',&
          'correspond to the number of overlap regions.'
     stop
  end if

END SUBROUTINE get_overlap_region_periodic
!%***
