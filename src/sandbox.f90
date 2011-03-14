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
  integer, parameter :: nlr=3,alr=1,blr=2
  integer :: ldim   ! dimension of lpsi
  integer :: fmin   ! min(1,nseg_f)
  logical :: isovrlp
  real(wp) :: scpr
  real(wp),dimension(nlr-1),parameter :: locrad=(/ 22.0, 22.0 /) !it has dimension nlr-1, because the last element of Llr is overlap
  real(wp), dimension(:),allocatable :: lpsi    ! local projection of |Psi>
  real(wp), dimension(:),allocatable :: lhpsi   ! local projection of H|Psi>
  type(locreg_descriptors), dimension(nlr+1) :: Llr
  real(wp), dimension(:), pointer :: potential
  real :: sum_pot  ! debug
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
  call createProjectorsArrays(iproc,n1,n2,n3,rxyz,atoms,orbs,&
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

! First, determine the localisation regions
! NOTE: we only pass nlr as the dimension of Llr (when it really is nlr+1)
  call determine_locreg(nlr-1,rxyz,locrad,in%hx,in%hy,in%hz,Glr,Llr)

! Plot the localization regions
!  call draw_locregs(2,in%hx,in%hy,in%hz,Llr)

! Second, calculate the overlap region
  call get_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr,Llr(nlr))

! Write some physical information on the overlap
  write(*,'(a24,3i4)')'Global region n1,n2,n3:',Glr%d%n1,Glr%d%n2,Glr%d%n3
  write(*,'(a27,f6.2,f6.2,f6.2)')'Global dimension (x,y,z):',Glr%d%n1*in%hx,Glr%d%n2*in%hy,Glr%d%n3*in%hz
  write(*,'(a17,f12.2)')'Global volume: ',Glr%d%n1*in%hx*Glr%d%n2*in%hy*Glr%d%n3*in%hz
  write(*,'(a25,3i4)')'Overlap region n1,n2,n3:',Llr(nlr)%d%n1,Llr(nlr)%d%n2,Llr(nlr)%d%n3
  write(*,'(a27,f6.2,f6.2,f6.2)')'Overlap dimension (x,y,z):',Llr(nlr)%d%n1*in%hx,Llr(nlr)%d%n2*in%hy,Llr(nlr)%d%n3*in%hz
  write(*,'(a17,f12.2)')'Overlap volume: ',Llr(nlr)%d%n1*in%hx*Llr(nlr)%d%n2*in%hy*Llr(nlr)%d%n3*in%hz

! Plot the overlap region
  if(isovrlp) then
     call draw_locregs(1,in%hx,in%hy,in%hz,Llr(nlr))

! Third, transform the wavefunction to this overlap region
     ilr = nlr
     ldim = (Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*orbs%norb*orbs%nspinor
     allocate(lpsi(ldim), stat=i_stat)
     call memocc(i_stat,lpsi,'lpsi',subname)

     call psi_to_locreg(Glr,ilr,ldim,Llr,lpsi,nlr,orbs,psi)

! Calculate overlap in localization region using ddot
     write(*,'(A26)') 'Overlap matrix with ddot:'
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
                dot(Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f,lpsi((Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
                lpsi((Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
        end do
     end do

! Calculate overlap in localization region using wpdot
     write(*,'(A27)')'Overlap matrix with wpdot:'
     fmin = min(Llr(ilr)%wfd%nseg_f,1)  ! checks if there is some fine_grid in the region, if not, do not shift keyv
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           call wpdot(Llr(ilr)%wfd%nvctr_c,Llr(ilr)%wfd%nvctr_f,Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nseg_f,&
&                  Llr(ilr)%wfd%keyv(1),Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+fmin),Llr(ilr)%wfd%keyg(1,1),&
&                  Llr(ilr)%wfd%keyg(1,Llr(ilr)%wfd%nseg_c+fmin),lpsi(1+(Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(iorb-1)),&
&                  lpsi(Llr(ilr)%wfd%nvctr_c+fmin+(Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(iorb-1)),&
&                  Llr(ilr)%wfd%nvctr_c,Llr(ilr)%wfd%nvctr_f,Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nseg_f,&
&                  Llr(ilr)%wfd%keyv(1),Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+fmin),Llr(ilr)%wfd%keyg(1,1),&
&                  Llr(ilr)%wfd%keyg(1,Llr(ilr)%wfd%nseg_c+fmin),lpsi((Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(jorb-1)+1),&
&                  lpsi(Llr(ilr)%wfd%nvctr_c+(Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(jorb-1)+fmin),scpr)
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,scpr
        end do
     end do
           
! Plot the orbitals that are inside the overlap region
! plot first orbital
     call plot_wf_sandbox('orb1_ovrlp',1,atoms,Llr(ilr),hxh,hyh,hzh,rxyz,lpsi,'')
! plot second orbital
     call plot_wf_sandbox('orb2_ovrlp',1,atoms,Llr(ilr),hxh,hyh,hzh,rxyz,&
&                     lpsi(Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f+1),'')
 
! Fourth, calculate scalar product of the overlap
!  call scalar_product_in_overlap_region()
  
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
     allocate(lhpsi(ldim), stat=i_stat)
     call memocc(i_stat,lpsi,'lpsi',subname)

     call psi_to_locreg(Glr,ilr,ldim,Llr,lhpsi,nlr,orbs,hpsi)

     write(*,'(A51)') 'Kinetic overlap matrix with ddot and localization:'
     do iorb=1,orbs%norb
        do jorb=iorb,orbs%norb
           write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
               dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,lpsi((Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(iorb-1)+1),1,&
               lhpsi((Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*(jorb-1)+1),1)
        end do
     end do

! ##################
! Non-local Projectors part
! ##################

! Need to project non-local projectors to overlap region
!   First make the descriptors
     call nlpspd_to_locreg(in,iproc,Glr,Llr(ilr),rxyz,atoms,orbs,&
&       radii_cf,in%frmult,in%frmult,in%hx,in%hy,in%hz,nlpspd,Lnlpspd)
    
!   Second fill the projectors

! Apply them on the wavefunctions in the overlap region

! Calculate dotprod: <Psi_a|p><p|Psi_b>



! Deallocations
     deallocate(lpsi,stat=i_stat)
     call memocc(i_stat,i_all,'lpsi',subname)

     deallocate(lhpsi,stat=i_stat)
     call memocc(i_stat,i_all,'lhpsi',subname)
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
subroutine psi_to_locreg(Glr,ilr,ldim,Llr,lpsi,nlr,orbs,psi)

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
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
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
  nseg = Llr(ilr)%wfd%nseg_c + Llr(ilr)%wfd%nseg_f
  lincrement = Llr(ilr)%wfd%nvctr_c + 7*Llr(ilr)%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0

! Initialize loc_psi
  call razero(lincrement*orbs%norb*orbs%nspinor,lpsi)
 
! Get the keymask: shift for every segment of Llr (with respect to Glr)
  allocate(keymask(2,nseg),stat=i_stat)
  call memocc(i_stat,keymask,'keymask',subname)

  call shift_locreg_indexes(Glr,Llr(ilr),keymask,nseg)

!####################################################
! Do coarse region
!####################################################
  do isegloc = 1,Llr(ilr)%wfd%nseg_c
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
  if(icheck .ne. Llr(ilr)%wfd%nvctr_c) then
    write(*,*)'There is an error in psi_to_locreg: number of coarse points used',icheck
    write(*,*)'is not equal to the number of coarse points in the region',Llr(ilr)%wfd%nvctr_c
  end if

!##############################################################
! Now do fine region
!##############################################################

  icheck = 0
  start = Llr(ilr)%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c
  lfinc  = Llr(ilr)%wfd%nvctr_f
  Gfinc = Glr%wfd%nvctr_f

  do isegloc = Llr(ilr)%wfd%nseg_c+1,nseg
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
  if(icheck .ne. Llr(ilr)%wfd%nvctr_f) then
    write(*,*)'There is an error in psi_to_locreg: number of fine points used',icheck
    write(*,*)'is not equal to the number of fine points in the region',Llr(ilr)%wfd%nvctr_f
  end if

  deallocate(keymask,stat=i_stat)
  call memocc(i_stat,i_all,'keymask',subname)

END SUBROUTINE psi_to_locreg
!%***

!##############################################################################################################################################
!!****f* BigDFT/get_overlap_region
!##############################################################################################################################################
!! FUNCTION
!!
!! SOURCE
!!
subroutine get_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr,Olr)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  logical, intent(out) :: isovrlp             ! True if there is an overlap
  type(locreg_descriptors),intent(out) :: Olr ! Overlap localization region 
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
  isovrlp = .false.  

! Determine if there is an overlap
! To do this, we compare axis by axis if there is an overlap.
! The cubes overlap if they overlap on all axis.
  if(((axmin .le. bxmax).and.(bxmin .le. axmax)) .and. &
&    ((aymin .le. bymax).and.(bymin .le. aymax)) .and. &
&    ((azmin .le. bzmax).and.(bzmin .le. azmax))) then
     isovrlp = .true.
  end if

! Now construct the Overlap localization region descriptor
! only if there is an overlap. The following only works
! when the previous test is successful. Note also that
! isx, isy and isz are necessarily in the Glr by construction
! of the Llrs, so don't need to test them.
  if(isovrlp) then
 
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

END SUBROUTINE get_overlap_region
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
 integer :: shift(3)  !shift between the beginning of the segment in Blr and the origine of Alr
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
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,Lnlpspd)

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
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz !atomic positions
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf  ! radii of the different atom types
  !#############################################
  !local variables
  !############################################
   integer :: iatom,igrid
   integer :: mproj
   integer :: mseg_c,mvctr_c,mseg_f,mvctr_f
   integer :: iat  ! index of the atoms
   real(gp) :: hhx,hhy,hhz   !reals to shorten name of variables
   integer,dimension(atoms%nat) :: projflg
   logical,dimension(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3) :: logrid

! Rename some variables
  hhx = input_parameters%hx
  hhy = input_parameters%hy
  hhz = input_parameters%hz
  
!Determine the number of projectors with components in locreg
! and also which atoms have such projectors
   call  number_of_projectors_in_locreg(atoms,Glr,Llr,nlpspd,mproj,projflg)

   Lnlpspd%nproj = mproj

!TESTS
  print *,'Number of projectors:', mproj
  print *,'Projflg', projflg
!ENDTESTS

  iat = 0
  do iatom = 1,atoms%nat
     if(projflg(iatom) == 0) cycle 
     iat = iat + 1    
!    Now we can determine the number of segments and elements of coarse grid
     call fill_logrid(atoms%geocode,Glr%d%n1,Glr%d%n2,Glr%d%n3,nlpspd%nboxp_c(1,1,iatom),&
&                     nlpspd%nboxp_c(2,1,iatom),nlpspd%nboxp_c(1,2,iatom),nlpspd%nboxp_c(2,2,iatom),&
&                     nlpspd%nboxp_c(1,3,iatom),nlpspd%nboxp_c(2,3,iatom),0,1,atoms%ntypes,&
&                     atoms%iatype(iatom),rxyz(1,iatom),radii_cf(1,3),cpmult,hhx,hhy,hhz,logrid)

     call number_of_projector_elements_in_locreg(iatom,1,atoms,Glr,Llr,logrid,nlpspd,mproj,projflg,mseg_c,mvctr_c)

!     Lnlpspd%nseg_p(2*iat-1) = mseg_c
!     Lnlpspd%nvctr_p(2*iat-1) = mvctr_c 

!TEST     
     print *,'Number of segments,coarse:',nlpspd%nseg_p(2*iatom-1),'for atom:',iatom
     print *,'Number of elements,coarse:',nlpspd%nvctr_p(2*iatom-1),'for atom:',iatom
     print *,'Number of segments,fine:',nlpspd%nseg_p(2*iatom),'for atom:',iatom
     print *,'Number of elements,fine:',nlpspd%nvctr_p(2*iatom),'for atom:',iatom
     print *,'Calculated now:'
     print *,'Number of segments,coarse:',mseg_c,'for atom:',iatom
     print *,'Number of elements,coarse:',mvctr_c,'for atom:',iatom
!END TEST

! Do the same for fine grid
     call fill_logrid(atoms%geocode,Glr%d%n1,Glr%d%n2,Glr%d%n3,nlpspd%nboxp_f(1,1,iatom),&
&                     nlpspd%nboxp_f(2,1,iatom),nlpspd%nboxp_f(1,2,iatom),nlpspd%nboxp_f(2,2,iatom),&
&                     nlpspd%nboxp_f(1,3,iatom),nlpspd%nboxp_f(2,3,iatom),0,1,atoms%ntypes,&
&                     atoms%iatype(iatom),rxyz(1,iatom),radii_cf(1,2),fpmult,hhx,hhy,hhz,logrid)

     call number_of_projector_elements_in_locreg(iatom,2,atoms,Glr,Llr,logrid,nlpspd,mproj,projflg,mseg_f,mvctr_f)

!     Lnlpspd%nseg_p(2*iat) = mseg_f
!     Lnlpspd%nvctr_p(2*iat) = mvctr_f

!TEST     
     print *,'Number of segments,fine:',mseg_f,'for atom:',iatom
     print *,'Number of elements,fine:',mvctr_f,'for atom:',iatom
     print *,'-------------------------'
!END TEST   

! Allocate the arrays

  end do

END SUBROUTINE nlpspd_to_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/number_of_projectors_in_locreg
!#############################################################################################################################################
!! FUNCTION: Calculates the number of projectors with components in the locreg
!!           It also returns a vector, projflg, which identifies the atoms with projectors inside the region
!!                      projflg = 0, no projectors in locreg
!!                      projflg = 1, projectros in locreg
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine number_of_projectors_in_locreg(atoms,Glr,Llr,nlpspd,mproj,projflg)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  type(atoms_data),intent(in) :: atoms        ! atoms descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  integer,intent(out) :: mproj  ! number of projectors
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  integer,dimension(atoms%nat),intent(out) :: projflg ! flag which is 1 if atom as projector components inside locreg
  !#######################################
  ! Local Variables
  !#######################################
  integer :: iatom            ! integer for loop
  integer :: bound(1:2,1:3)   ! bound of locreg
  integer :: nproj            ! temporary number of projectors

! bounds of the localization region
  !lower bound
  bound(1,1) = Llr%ns1 - Glr%ns1
  bound(1,2) = Llr%ns2 - Glr%ns2
  bound(1,3) = Llr%ns3 - Glr%ns3

  !upper bound
  bound(2,1) = Llr%ns1 + Llr%d%n1 - Glr%ns1
  bound(2,2) = Llr%ns2 + Llr%d%n2 - Glr%ns2
  bound(2,3) = Llr%ns3 + Llr%d%n3 - Glr%ns3 

  projflg = 0
  mproj = 0
  do iatom=1,atoms%nat

! check if projector of atom iatom overlap the locreg (coarse grid)
     if(nlpspd%nboxp_c(1,1,iatom) < bound(2,1) .and. nlpspd%nboxp_c(2,1,iatom) > bound(1,1) .and. &
&       nlpspd%nboxp_c(1,2,iatom) < bound(2,2) .and. nlpspd%nboxp_c(2,2,iatom) > bound(1,2) .and. &
&       nlpspd%nboxp_c(1,3,iatom) < bound(2,3) .and. nlpspd%nboxp_c(2,3,iatom) > bound(1,3)) then
        
        call numb_proj(atoms%iatype(iatom),atoms%ntypes,atoms%psppar,atoms%npspcode,nproj)
        mproj = mproj + nproj
        if(nproj > 0) projflg(iatom) = 1
!        nbox_c(,iatom)
     end if

! check if projector of atom iatom overlap the locreg (fine grid)
! Only have to do it if the atom is not yet selected
     if (projflg(iatom) .eq. 0) then
        if(nlpspd%nboxp_f(1,1,iatom) < bound(2,1) .and. nlpspd%nboxp_f(2,1,iatom) > bound(1,1) .and. &
&          nlpspd%nboxp_f(1,2,iatom) < bound(2,2) .and. nlpspd%nboxp_f(2,2,iatom) > bound(1,2) .and.&
&          nlpspd%nboxp_f(1,3,iatom) < bound(2,3) .and. nlpspd%nboxp_f(2,3,iatom) > bound(1,3)) then
          
           call numb_proj(atoms%iatype(iatom),atoms%ntypes,atoms%psppar,atoms%npspcode,nproj)
           mproj = mproj + nproj
           if(nproj > 0) projflg(iatom) = 1
        end if
     end if        
  end do

END SUBROUTINE number_of_projectors_in_locreg
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
subroutine number_of_projector_elements_in_locreg(iatom,igrid,atoms,Glr,Llr,logrid,nlpspd,mproj,projflg,mseg,mvctr)

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
  integer,dimension(atoms%nat),intent(in) :: projflg ! flag which is 1 if atom as projector components inside locreg
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

  !upper bound
  bound(2,1) = Llr%ns1 + Llr%d%n1 - Glr%ns1
  bound(2,2) = Llr%ns2 + Llr%d%n2 - Glr%ns2
  bound(2,3) = Llr%ns3 + Llr%d%n3 - Glr%ns3 


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
