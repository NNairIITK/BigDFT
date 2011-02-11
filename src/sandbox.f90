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
  real(wp), dimension(:),allocatable :: apsi_c,bpsi_c
  real(wp), dimension(:,:),allocatable :: apsi_f,bpsi_f
  real(gp) :: scpr
  real(gp), dimension(:,:),allocatable :: keyag_c,keybg_c,keyag_f,keybg_f
  real(gp), dimension(:),allocatable :: keyav_c,keybv_c,keyav_f,keybv_f
  ! arrays for DIIS convergence accelerator
  real(wp), dimension(:,:,:), pointer :: ads

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

  pot_ion=0.0d0

  !create input guess wavefunction
  call psi_from_gaussians(iproc,nproc,atoms,orbs,Glr,rxyz,in%hx,in%hy,in%hz,in%nspin,psi)

  !calculate scalar product between wavefunctions
  !this scheme works only in sequential
  write(*,'(A25)') 'Overlap matrix with ddot:'
  do iorb=1,orbs%norb
     do jorb=iorb,orbs%norb
        write(*,'(2(i4),1x,1pe19.12)')iorb,jorb,&
             dot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(iorb-1)+1),1,&
             psi((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*(jorb-1)+1),1)
     end do
  end do

!#######################################################################################################################
!  write(*,'(A25)') 'Overlap matrix using wpdot:'

  allocate(apsi_c(Glr%wfd%nvctr_c))
  allocate(bpsi_c(Glr%wfd%nvctr_c))
  allocate(apsi_f(7,Glr%wfd%nvctr_f))
  allocate(bpsi_f(7,Glr%wfd%nvctr_f))
  allocate(keyag_c(2,Glr%wfd%nvctr_c))
  allocate(keybg_c(2,Glr%wfd%nvctr_c)) 
  allocate(keyag_f(2,Glr%wfd%nvctr_f)) 
  allocate(keybg_f(2,Glr%wfd%nvctr_f))
  allocate(keyav_c(Glr%wfd%nvctr_c)) 
  allocate(keybv_c(Glr%wfd%nvctr_c)) 
  allocate(keyav_f(Glr%wfd%nvctr_f)) 
  allocate(keybv_f(Glr%wfd%nvctr_f)) 

! First format the two orbitals 
  apsi_c(:) = psi(1:Glr%wfd%nvctr_c)
  bpsi_c(:) = psi(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f+1:2*Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
 
  keyag_c(:,:) = Glr%wfd%keyg(:,1:Glr%wfd%nseg_c)    ! dim are (2,nseg_c)
  keybg_c(:,:) = Glr%wfd%keyg(:,Glr%wfd%nseg_c+Glr%wfd%nseg_f+1:2*Glr%wfd%nseg_c+Glr%wfd%nseg_f)
  keyag_f(:,:) = Glr%wfd%keyg(:,Glr%wfd%nseg_c+1:Glr%wfd%nseg_c+Glr%wfd%nseg_f)
  keybg_f(:,:) = Glr%wfd%keyg(:,2*Glr%wfd%nseg_c+Glr%wfd%nseg_f+1:2*Glr%wfd%nseg_c+2*Glr%wfd%nseg_f)
  keyav_c(:) = Glr%wfd%keyv(1:Glr%wfd%nseg_c)
  keybv_c(:) = Glr%wfd%keyv(Glr%wfd%nseg_c+Glr%wfd%nseg_f+1:2*Glr%wfd%nseg_c+Glr%wfd%nseg_f)
  keyav_f(:) = Glr%wfd%keyv(Glr%wfd%nseg_c+1:Glr%wfd%nseg_c+Glr%wfd%nseg_f)
  keybv_f(:) = Glr%wfd%keyv(2*Glr%wfd%nseg_c+Glr%wfd%nseg_f+1:2*Glr%wfd%nseg_c+2*Glr%wfd%nseg_f)

  do ifine=1,7
     apsi_f(ifine,:) = psi(Glr%wfd%nvctr_c+(ifine-1)*Glr%wfd%nvctr_f+1:Glr%wfd%nvctr_c+ifine*Glr%wfd%nvctr_f)
     bpsi_f(ifine,:) = psi(2*Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f+(ifine-1)*Glr%wfd%nvctr_f+1:&
     &                     Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f+(ifine)*Glr%wfd%nvctr_f)
  end do
!  write(*,'(a,1x,1pe19.12)')'Overlap with dot',dot(Glr%wfd%nvctr_f,apsi_c(:),1,bpsi_c(:),1)
!  & dot(Glr%wfd%nvctr_f,apsi_f(1,:),1,bpsi_f(1,:),1) + dot(Glr%wfd%nvctr_f,apsi_f(2,:),1,bpsi_f(2,:),1) + &
!  & dot(Glr%wfd%nvctr_f,apsi_f(3,:),1,bpsi_f(3,:),1) + dot(Glr%wfd%nvctr_f,apsi_f(4,:),1,bpsi_f(4,:),1) + &
!  & dot(Glr%wfd%nvctr_f,apsi_f(5,:),1,bpsi_f(5,:),1) + dot(Glr%wfd%nvctr_f,apsi_f(6,:),1,bpsi_f(6,:),1) !+ &
!  & dot(Glr%wfd%nvctr_f,apsi_f(7,:),1,bpsi_f(7,:),1)

! Then call the subroutine calculating the overlap (only calculate the overlap and not the 4 permutations)
  call wpdot( Glr%wfd%nvctr_c,Glr%wfd%nvctr_f,Glr%wfd%nseg_c,Glr%wfd%nseg_f,keyav_c,keyav_f,&
& keyag_c,keyag_f,apsi_c,apsi_f,Glr%wfd%nvctr_c,Glr%wfd%nvctr_f,Glr%wfd%nseg_c,Glr%wfd%nseg_f,&
& keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
  stop
  write(*,'(2(i4),1x,1pe19.12)') 1,2,scpr

  deallocate(apsi_c)
  deallocate(bpsi_c)
  deallocate(apsi_f)
  deallocate(bpsi_f)

!#######################################################################################################################

!!$  psi=0.0d0
!!$  ttsum=0.0d0
!!$  do i=1,orbs%npsidim
!!$     do j=0,iproc-1
!!$        call random_number(ttr)
!!$     end do
!!$     call random_number(ttr)
!!$     psi(i)=real(ttr,wp)*0.01_wp
!!$     ttsum=ttsum+psi(i)
!!$     do j=iproc+1,nproc
!!$        call random_number(ttr)
!!$     end do
!!$  end do
!!$
!!$  psi=1.d0


  print *,'norbs',orbs%norb

  !plot first orbital
  call plot_wf_sandbox('iter0-1',1,atoms,Glr,hxh,hyh,hzh,rxyz,psi,'')
  !plot second orbital
  call plot_wf_sandbox('iter0-2',1,atoms,Glr,hxh,hyh,hzh,rxyz,psi(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f+1),'')

  

  stop
  !othogonalise them
  !transpose the psi wavefunction
  call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,&
       psi,work=hpsi)
  call orthogonalize(iproc,nproc,orbs,comms,Glr%wfd,psi,in)
  !untranspose psi
  call untranspose_v(iproc,nproc,orbs,Glr%wfd,comms,psi,work=hpsi)

  allocate(orbs%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%eval,'eval',subname)

  orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0

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

     call HamiltonianApplication(iproc,nproc,atoms,orbs,in%hx,in%hy,in%hz,rxyz,&
          nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
          pot_ion,psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,in%nspin,GPU)

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
  integer :: i_stat,i_all
  integer :: nl1,nl2,nl3,n1i,n2i,n3i,n1,n2,n3,i1,i2,i3,nu1,nu2,nu3,iat
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
     z=hzh*real(i3,gp)
     do i2=-nl2,2*n2+1+nu2
        y=hyh*real(i2,gp)
        do i1=-nl1,2*n1+1+nu1
           x=hxh*real(i1,gp)
           if (i3 == n3+1 .and. i2 ==  n2+1) then
           !print *,'value of the center',rxyz(1,iat)
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
