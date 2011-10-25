!> @file
!!  Application of the Hamiltonian + orthonormalize constraints
!! @author
!!    Copyright (C) 2007-2011 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Application of the Hamiltonian
subroutine HamiltonianApplication2(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
     proj,Lzd,ngatherarr,pot,psi,hpsi,&
     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel,orbsocc,psirocc)
  use module_base
  use module_types
  use module_interfaces, except_this_one => HamiltonianApplication2
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: Lzd
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(Lzd%Lnprojel), intent(in) :: proj
  real(wp), dimension(Lzd%Lpsidimtot), intent(in) :: psi
  real(wp), dimension(:), pointer :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
  real(wp), target, dimension(Lzd%Lpsidimtot), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  real(dp), dimension(*), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  !local variables
  real(gp), dimension(2,orbs%norbp) :: ekin
  real(gp), dimension(2,orbs%norbp) :: epot
  real(wp), dimension(:), pointer :: hpsi2
  character(len=*), parameter :: subname='HamiltonianApplication2'
  logical :: exctX,op2p
  integer :: i_all,i_stat,ierr,iorb,n3p,ispot,istart_c,iat
  integer :: istart_ck,isorb,ieorb,ikpt,ispsi_k,nspinor,ispsi
!OCL  integer, dimension(3) :: periodic
!OCL  real(wp) :: maxdiff
!OCL  real(gp) :: eproj,ek_fake,ep_fake
  real(gp), dimension(3,2) :: wrkallred
!OCL  real(wp), dimension(:), allocatable :: hpsi_OCL
 integer :: size_pot,size_potxc
 real(wp), dimension(:), allocatable :: potxc

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  !check if the potential has been associated
  if (.not. associated(pot)) then
     if (iproc ==0) then
        write(*,*)' ERROR, HamiltonianApplication2, potential not associated!'
        stop
     end if
  end if


!##################################################################################################
! Exact exchange potential calculation
! For Linear scaling, keep cubic code (using whole simulation box and total orbitals) for now
! Defined new subroutine cubic_exact_exchange to ease futur change
!##################################################################################################
  !initialise exact exchange energy 
  op2p=(eexctX == -99.0_gp)
  eexctX=0.0_gp

  exctX = libxc_functionals_exctXfac() /= 0.0_gp

  size_potxc = 0
  if (exctX) then
     size_potxc = max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%norbp,ngatherarr(0,1)*orbs%norb),1)
     allocate(potxc(size_potxc+ndebug),stat=i_stat)
     call memocc(i_stat,potxc,'potxc',subname)
     if (present(pkernel) .and. present(orbsocc) .and. present(psirocc)) then
        call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
             ngatherarr,psi,potxc,eexctX,pkernel,orbsocc,psirocc)
     else if(present(pkernel)) then
         call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
             ngatherarr,psi,potxc,eexctX,pkernel=pkernel)
     else
         call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
             ngatherarr,psi,potxc,eexctX)
     end if
  else
     allocate(potxc(1+ndebug),stat=i_stat)
     call memocc(i_stat,potxc,'potxc',subname)
  end if

!################################################################################################
! Application of the local potential
!###############################################################################################

  call timing(iproc,'ApplyLocPotKin','ON')

  !apply the local hamiltonian for each of the orbitals
  !given to each processor
  !pot=0.d0
  !psi=1.d0
  
!switch between GPU/CPU treatment
!  do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!       call random_number(psi(i))
!  end do
  if(OCLconv .and. ASYNCconv) then
    allocate(hpsi2((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),stat=i_stat)
    call memocc(i_stat,hpsi2,'hpsi2',subname)
    hpsi(:)=0.0
  else
    hpsi2 => hpsi
  end if

  if (GPUconv) then
     if(Lzd%linear) stop 'HamiltonianApplication: Linear scaling not implemented with GPU convolutions yet'
     call local_hamiltonian_GPU(iproc,orbs,Lzd%Glr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum,GPU)
  else if (OCLconv) then
     if(Lzd%linear) stop 'HamiltonianApplication: Linear scaling not implemented with OCL yet'
     call local_hamiltonian_OCL(iproc,orbs,Lzd%Glr,hx,hy,hz,nspin,pot,psi,hpsi2,ekin_sum,epot_sum,GPU,ekin,epot)
  else
     call local_hamiltonian2(iproc,exctX,orbs,Lzd,hx,hy,hz,nspin,pot,size_potxc,potxc,psi,hpsi,ekin_sum,epot_sum)
  end if
  !test part to check the results wrt OCL convolutions
!!$  if (OCLconv) then
!!$     allocate(hpsi_OCL((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
!!$     call memocc(i_stat,hpsi_OCL,'hpsi_OCL',subname)
!!$     print *,'fulllocam',GPU%full_locham
!!$     call local_hamiltonian_OCL(iproc,orbs,at%geocode,lr,hx,hy,hz,nspin,pot,psi,hpsi,ek_fake,ep_fake,GPU)
!!$     maxdiff=0.0_wp
!!$     do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!!$        maxdiff=max(maxdiff,abs(hpsi(i)-hpsi_OCL(i)))
!!$     end do
!!$     print *,'maxdiff',maxdiff
!!$     print *,'ekin_diff',abs(ek_fake-ekin_sum)
!!$     print *,'epot_diff',abs(ep_fake-epot_sum)
!!$     i_all=-product(shape(hpsi_OCL))*kind(hpsi_OCL)
!!$     deallocate(hpsi_OCL,stat=i_stat)
!!$     call memocc(i_stat,i_all,'hpsi_OCL',subname)
!!$  end if
  call timing(iproc,'ApplyLocPotKin','OF')

!#############################################################################################################################
! Applying the NonLocal part of the  pseudopotentials
!#############################################################################################################################

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')

  !here the localisation region should be changed, temporary only for cubic approach
  eproj_sum=0.0_gp
  !apply the projectors following the strategy (On-the-fly calculation or not)
  if (DistProjApply .and. .not.Lzd%linear) then
     call applyprojectorsonthefly(iproc,orbs,at,Lzd%Glr,&
          rxyz,hx,hy,hz,Lzd%Glr%wfd,Lzd%Gnlpspd,proj,psi,hpsi,eproj_sum)
  else if(orbs%norbp > 0 .and. .not.Lzd%linear) then
     !apply the projectors  k-point of the processor
     !starting k-point
     ikpt=orbs%iokpt(1)
     istart_ck=1
     ispsi_k=1
     loop_kpt: do
        call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

        ! loop over all my orbitals
        ispsi=ispsi_k
        do iorb=isorb,ieorb
           istart_c=istart_ck
           do iat=1,at%nat
              call apply_atproj_iorb(iat,iorb,istart_c,at,orbs,Lzd%Glr%wfd,Lzd%Gnlpspd,&
                   proj,psi(ispsi),hpsi(ispsi),eproj_sum)
           end do
           ispsi=ispsi+(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*nspinor
        end do
        istart_ck=istart_c
        if (ieorb == orbs%norbp) exit loop_kpt
        ikpt=ikpt+1
        ispsi_k=ispsi
     end do loop_kpt
     if (istart_ck-1 /= Lzd%Gnlpspd%nprojel) stop 'incorrect once-and-for-all psp application'
     if (ispsi-1 /= (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'

! Linear part of the NLPSP
  else if(orbs%norbp > 0 .and. Lzd%linear) then
     call ApplyProjectorsLinear(iproc,hx,hy,hz,at,Lzd,orbs,rxyz,psi,hpsi,eproj_sum)
  end if

  if(OCLconv .and. ASYNCconv) then
    call finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU,ekin,epot)
    call daxpy(size(hpsi), 1.0_wp, hpsi2(1), 1, hpsi(1),1)
    i_all=-product(shape(hpsi2))*kind(hpsi2)
    deallocate(hpsi2,stat=i_stat)
    call memocc(i_stat,i_all,'hpsi2',subname)
  endif

  call timing(iproc,'ApplyProj     ','OF')

  !energies reduction
  if (nproc > 1) then
     wrkallred(1,2)=ekin_sum
     wrkallred(2,2)=epot_sum
     wrkallred(3,2)=eproj_sum
     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
          mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
     ekin_sum=wrkallred(1,1)
     epot_sum=wrkallred(2,1)
     eproj_sum=wrkallred(3,1)
  endif



  !up to this point, the value of the potential energy is 
  !only taking into account the local potential part
  !whereas it should consider also the value coming from the 
  !exact exchange operator (twice the exact exchange energy)

  if (exctX) epot_sum=epot_sum+2.0_gp*eexctX

END SUBROUTINE HamiltonianApplication2





subroutine HamiltonianApplication3(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
     proj,Lzd,ngatherarr,Lpot,psi,hpsi,&
     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,withConfinement,energyReductionFlag, &
     pkernel,orbsocc,psirocc,lin, confinementCenter)
  use module_base
  use module_types
  use module_interfaces, except_this_one => HamiltonianApplication3
  use module_xc
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: Lzd
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(Lzd%Lnprojel), intent(in) :: proj
  real(wp), dimension(Lzd%Lpsidimtot), intent(in) :: psi
  real(wp), dimension(lzd%ndimpotisf) :: Lpot
  real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
  real(wp), target, dimension(Lzd%Lpsidimtot), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  logical,intent(in):: withConfinement
  logical,intent(in):: energyReductionFlag
  real(dp), dimension(*), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  type(linearParameters),intent(in),optional:: lin
  integer,dimension(orbs%norbp),intent(in),optional:: confinementCenter
  !local variables
  real(gp), dimension(2,orbs%norbp) :: ekin
  real(gp), dimension(2,orbs%norbp) :: epot
  real(wp), dimension(:), pointer :: hpsi2
  character(len=*), parameter :: subname='HamiltonianApplication3'
  logical :: exctX,op2p
  integer :: i_all,i_stat,ierr,iorb,n3p,ispot,istart_c,iat, ist
  integer :: istart_ck,isorb,ieorb,ikpt,ispsi_k,nspinor,ispsi
!OCL  integer, dimension(3) :: periodic
!OCL  real(wp) :: maxdiff
!OCL  real(gp) :: eproj,ek_fake,ep_fake
  real(gp), dimension(3,2) :: wrkallred
!OCL  real(wp), dimension(:), allocatable :: hpsi_OCL
 integer :: size_pot,size_potxc

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  !!!check if the potential has been associated
  !!if (.not. associated(Lpot)) then
  !!   if (iproc ==0) then
  !!      write(*,*)' ERROR, HamiltonianApplication2, potential not associated!'
  !!      stop
  !!   end if
  !!end if


!##################################################################################################
! Exact exchange potential calculation
! For Linear scaling, keep cubic code (using whole simulation box and total orbitals) for now
! Defined new subroutine cubic_exact_exchange to ease futur change
!##################################################################################################
  !initialise exact exchange energy 
  op2p=(eexctX ==  UNINITIALIZED(1.0_gp))
  eexctX=0.0_gp

  exctX = xc_exctXfac() /= 0.0_gp

  size_potxc = 0
  if (exctX .and. .not.present(lin)) then
     size_potxc = max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%norbp,ngatherarr(0,1)*orbs%norb),1)
     ist=lzd%ndimpotisf-size_potxc+1
     if (present(pkernel) .and. present(orbsocc) .and. present(psirocc)) then
        call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
             ngatherarr,psi,Lpot(ist),eexctX,pkernel,orbsocc,psirocc)
     else if(present(pkernel)) then

         call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
             ngatherarr,psi,Lpot(ist),eexctX,pkernel=pkernel)
     else
         call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
             ngatherarr,psi,Lpot(ist),eexctX)
     end if
  end if
!################################################################################################
! Application of the local potential
!###############################################################################################

  call timing(iproc,'ApplyLocPotKin','ON')

  !apply the local hamiltonian for each of the orbitals
  !given to each processor
  !pot=0.d0
  !psi=1.d0
  
!switch between GPU/CPU treatment
!  do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!       call random_number(psi(i))
!  end do
  if(OCLconv .and. ASYNCconv) then
    allocate(hpsi2((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),stat=i_stat)
    call memocc(i_stat,hpsi2,'hpsi2',subname)
    hpsi(:)=0.0
  else
    hpsi2 => hpsi
  end if

  if (GPUconv) then
     if(Lzd%linear) stop 'HamiltonianApplication: Linear scaling not implemented with GPU convolutions yet'
     call local_hamiltonian_GPU(iproc,orbs,Lzd%Glr,hx,hy,hz,nspin,Lpot,psi,hpsi,ekin_sum,epot_sum,GPU)
  else if (OCLconv) then
     if(Lzd%linear) stop 'HamiltonianApplication: Linear scaling not implemented with OCL yet'
     call local_hamiltonian_OCL(iproc,orbs,Lzd%Glr,hx,hy,hz,nspin,Lpot,psi,hpsi2,ekin_sum,epot_sum,GPU,ekin,epot)
  else
     !call local_hamiltonian2(iproc,exctX,orbs,Lzd,hx,hy,hz,nspin,pot,size_potxc,potxc,psi,hpsi,ekin_sum,epot_sum)
     if(present(lin)) then
         if(.not.present(confinementCenter)) then
             write(*,'(a,i0,a)') "ERROR on processes ",iproc,": 'confinementCenter' must be present when passing 'lin'!"
             stop
         end if
         call local_hamiltonian3(iproc,exctx,orbs,lzd,hx,hy,hz,&
              nspin,Lpot,psi,hpsi,ekin_sum,epot_sum,&
              withconfinement, at, rxyz, ist, lin, confinementCenter)
     else
         call local_hamiltonian3(iproc,exctx,orbs,lzd,hx,hy,hz,&
              nspin,Lpot,psi,hpsi,ekin_sum,epot_sum,&
              withconfinement, at, rxyz, ist)
     end if
  end if

  !test part to check the results wrt OCL convolutions
!!$  if (OCLconv) then
!!$     allocate(hpsi_OCL((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
!!$     call memocc(i_stat,hpsi_OCL,'hpsi_OCL',subname)
!!$     print *,'fulllocam',GPU%full_locham
!!$     call local_hamiltonian_OCL(iproc,orbs,at%geocode,lr,hx,hy,hz,nspin,pot,psi,hpsi,ek_fake,ep_fake,GPU)
!!$     maxdiff=0.0_wp
!!$     do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!!$        maxdiff=max(maxdiff,abs(hpsi(i)-hpsi_OCL(i)))
!!$     end do
!!$     print *,'maxdiff',maxdiff
!!$     print *,'ekin_diff',abs(ek_fake-ekin_sum)
!!$     print *,'epot_diff',abs(ep_fake-epot_sum)
!!$     i_all=-product(shape(hpsi_OCL))*kind(hpsi_OCL)
!!$     deallocate(hpsi_OCL,stat=i_stat)
!!$     call memocc(i_stat,i_all,'hpsi_OCL',subname)
!!$  end if
  call timing(iproc,'ApplyLocPotKin','OF')

!#############################################################################################################################
! Applying the NonLocal part of the  pseudopotentials
!#############################################################################################################################

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')

  !here the localisation region should be changed, temporary only for cubic approach
  eproj_sum=0.0_gp
  !apply the projectors following the strategy (On-the-fly calculation or not)
  if (DistProjApply .and. .not.Lzd%linear) then
     call applyprojectorsonthefly(iproc,orbs,at,Lzd%Glr,&
          rxyz,hx,hy,hz,Lzd%Glr%wfd,Lzd%Gnlpspd,proj,psi,hpsi,eproj_sum)
  else if(orbs%norbp > 0 .and. .not.Lzd%linear) then
     !apply the projectors  k-point of the processor
     !starting k-point
     ikpt=orbs%iokpt(1)
     istart_ck=1
     ispsi_k=1
     loop_kpt: do
        call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

        ! loop over all my orbitals
        ispsi=ispsi_k
        do iorb=isorb,ieorb
           istart_c=istart_ck
           do iat=1,at%nat
              call apply_atproj_iorb(iat,iorb,istart_c,at,orbs,Lzd%Glr%wfd,Lzd%Gnlpspd,&
                   proj,psi(ispsi),hpsi(ispsi),eproj_sum)
           end do
           ispsi=ispsi+(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*nspinor
        end do
        istart_ck=istart_c
        if (ieorb == orbs%norbp) exit loop_kpt
        ikpt=ikpt+1
        ispsi_k=ispsi
     end do loop_kpt
     if (istart_ck-1 /= Lzd%Gnlpspd%nprojel) stop 'incorrect once-and-for-all psp application'
     if (ispsi-1 /= (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'

! Linear part of the NLPSP
  else if(orbs%norbp > 0 .and. Lzd%linear) then
     call ApplyProjectorsLinear(iproc,hx,hy,hz,at,Lzd,orbs,rxyz,psi,hpsi,eproj_sum)
  end if

  if(OCLconv .and. ASYNCconv) then
    call finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU,ekin,epot)
    call daxpy(size(hpsi), 1.0_wp, hpsi2(1), 1, hpsi(1),1)
    i_all=-product(shape(hpsi2))*kind(hpsi2)
    deallocate(hpsi2,stat=i_stat)
    call memocc(i_stat,i_all,'hpsi2',subname)
  endif

  call timing(iproc,'ApplyProj     ','OF')


  if(energyReductionFlag) then
      !energies reduction
      if (nproc > 1) then
         wrkallred(1,2)=ekin_sum
         wrkallred(2,2)=epot_sum
         wrkallred(3,2)=eproj_sum
         call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
              mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
         ekin_sum=wrkallred(1,1)
         epot_sum=wrkallred(2,1)
         eproj_sum=wrkallred(3,1)
      endif
  else
      ! Do not sum up the energies, but set them to default values
      ekin_sum=0.d0
      epot_sum=0.d0
      eproj_sum=0.d0
  end if



  !up to this point, the value of the potential energy is 
  !only taking into account the local potential part
  !whereas it should consider also the value coming from the 
  !exact exchange operator (twice the exact exchange energy)

  if (exctX) epot_sum=epot_sum+2.0_gp*eexctX

END SUBROUTINE HamiltonianApplication3

!!  Routine to calculate the action of the hamiltonian
!! @author
!!   Copyright (C) 2005-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!>   Calculate the action of the local hamiltonian on the orbitals
subroutine local_hamiltonian2(iproc,exctX,orbs,Lzd,hx,hy,hz,&
     nspin,pot,size_potxc,potxc,psi,hpsi,ekin_sum,epot_sum)
  use module_base
  use module_types
  use module_interfaces, except_this_one => local_hamiltonian2
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nspin
  integer,intent(in) :: size_potxc
  real(gp), intent(in) :: hx,hy,hz
  logical, intent(in) :: exctX
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  real(wp), dimension(Lzd%Lpsidimtot), intent(in) :: psi
  real(wp), dimension(size_potxc),intent(in) :: potxc
  real(wp), dimension(*),target :: pot
  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension(Lzd%Lpsidimtot), intent(out) :: hpsi
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian2'
  integer :: i_all,i_stat,iorb,npot,nsoffset,oidx,dimwf
  integer :: ilr,size_pot,size_Lpot,ispot,iels
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,etest
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable,target :: psir
  real(wp), dimension(:),allocatable,target :: Lpot
  real(wp), dimension(:), pointer :: potu

  exctXcoeff=libxc_functionals_exctXfac()

  !components of the potential
  npot=orbs%nspinor
  if (orbs%nspinor == 2) npot=1

  ekin_sum=0.0_gp
  epot_sum=0.0_gp
  etest=0.0_gp

  ! if doing a cubic calculation, add the exact-exchange to the potential
  if(exctX .and. size_potxc > 0 .and. .not. Lzd%linear) then
     ispot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin+1
     ! fill in the exact-exchange part
     call dcopy(size_potxc,potxc(1),1,pot(ispot),1)
     !do iels = 1, size_potxc
     !   if(iproc==0) write(*,'(a,3i9)') 'size(pot), ispot+iels-1, iels', &
     !       Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin + size_potxc, ispot+iels-1, iels
     !   pot(ispot+iels-1) = potxc(iels)
     !end do
  end if
  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
     
     if(Lzd%linear) then
        !Determine size of local potential
        if (exctX) then
           size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin + size_potxc
           size_Lpot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin + size_potxc
        else
           size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin
           size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin
        end if

        !Allocate the local potential
        allocate(Lpot(size_Lpot+ndebug), stat=i_stat)
        call memocc(i_stat,Lpot,'Lpot',subname)
        call razero(size_Lpot,Lpot)
   
        ! Cut the potential into locreg pieces
        call global_to_local(Lzd%Glr,Lzd%Llr(ilr),nspin,size_pot,size_Lpot,pot,Lpot)
   
        ! Set some quantities: ispot=shift for potential, dimwf=dimension of wavefunction
        ispot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin+1
        dimwf = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
   
        ! fill in the exact-exchange part
        if (size_potxc > 0) then
           do iels = 1, size_potxc
              Lpot(ispot+iels-1) = potxc(iels)
           end do
        end if
     end if

     !initialise the work arrays
     call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,wrk_lh)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     if(orbs%spinsgn(iorb+orbs%isorb)>0.0_gp .or. nspin == 1 .or. nspin == 4 ) then
        nsoffset=1
     else
        nsoffset=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i+1
     end if

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     call daub_to_isf_locham(orbs%nspinor,Lzd%Llr(ilr),wrk_lh,psi(1+oidx),psir)

     !apply the potential to the psir wavefunction and calculate potential energy
     select case(Lzd%Llr(ilr)%geocode)
     case('F')
           if(.not. Lzd%linear) then
              call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
                pot(nsoffset),epot,Lzd%Llr(ilr)%bounds%ibyyzz_r)
           else
              call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
                Lpot(nsoffset),epot,Lzd%Llr(ilr)%bounds%ibyyzz_r) 
           end if
     case('P')
        !here the hybrid BC act the same way
        if(.not. Lzd%linear) then
           call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,0,0,0,0,orbs%nspinor,npot,psir,&
             pot(nsoffset),epot)
        else
           call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,0,0,0,0,orbs%nspinor,npot,psir,&
             Lpot(nsoffset),epot)
        end if
     case('S')
        if(.not. Lzd%linear) then
           call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,0,1,0,0,orbs%nspinor,npot,psir,&
             pot(nsoffset),epot)
        else
           call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,0,1,0,0,orbs%nspinor,npot,psir,&
             Lpot(nsoffset),epot)
        end if
     end select

     !k-point values, if present
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))

     if (exctXcoeff /= 0.0_gp) then
        ispot=1+Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*(nspin+iorb-1)
        !add to the psir function the part of the potential coming from the exact exchange
        if(.not. Lzd%linear) then
        call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
        else
        call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,exctXcoeff,Lpot(ispot),1,psir(1,1),1)
        end if
     end if

     !apply the kinetic term, sum with the potential and transform back to Daubechies basis
     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,Lzd%Llr(ilr),wrk_lh,&
          psir,hpsi(1+oidx),ekin)

     ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot

     if(Lzd%linear) then
        i_all=-product(shape(Lpot))*kind(Lpot)
        deallocate(Lpot,stat=i_stat)
        call memocc(i_stat,i_all,'Lpot',subname)
     end if

     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     call deallocate_work_arrays_locham(Lzd%Llr(ilr),wrk_lh)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
  enddo

END SUBROUTINE local_hamiltonian2

subroutine local_hamiltonian3(iproc,exctX,orbs,Lzd,hx,hy,hz,&
     nspin,Lpot,psi,hpsi,ekin_sum,epot_sum,&
     withConfinement, at, rxyz, istexct, lin, confinementCenter)
  use module_base
  use module_types
  use module_interfaces, except_this_one => local_hamiltonian3
  use module_xc
  implicit none
  integer, intent(in) :: iproc,nspin, istexct
  real(gp), intent(in) :: hx,hy,hz
  logical, intent(in) :: exctX
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  real(wp), dimension(Lzd%Lpsidimtot), intent(in) :: psi
  real(wp), dimension(lzd%ndimpotisf),target :: Lpot
  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension(Lzd%Lpsidimtot), intent(out) :: hpsi
  logical,intent(in):: withConfinement
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(linearParameters),intent(in),optional:: lin
  integer,dimension(orbs%norbp),intent(in),optional:: confinementCenter
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian3'
  integer :: i_all,i_stat,iorb,npot,oidx, icenter
  integer :: ilr,size_pot,size_Lpot,ispot, ilrold, nsoffset
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,etest, hxh, hyh, hzh
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable,target :: psir

  if(withConfinement .and. .not.present(lin)) then
      stop "'lin' must be present when 'withConfinement' is true!"
  end if

  exctXcoeff=xc_exctXfac()

  !components of the potential
  npot=orbs%nspinor
  if (orbs%nspinor == 2) npot=1

  ekin_sum=0.0_gp
  epot_sum=0.0_gp
  etest=0.0_gp

  oidx = 0
  ilrold=-1
  nsoffset=1
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
     if(.not.lzd%doHamAppl(ilr)) cycle
  
     !initialise the work arrays
     call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,wrk_lh)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     call daub_to_isf_locham(orbs%nspinor,Lzd%Llr(ilr),wrk_lh,psi(1+oidx),psir)
     !apply the potential to the psir wavefunction and calculate potential energy
     select case(Lzd%Llr(ilr)%geocode)
     case('F')
        if(withConfinement) then
            hxh=.5d0*hx
            hyh=.5d0*hy
            hzh=.5d0*hz
            icenter=confinementCenter(iorb)
            !!call apply_potentialConfinement2(lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
            !!     Lpot(orbs%ispot(iorb)),epot, rxyz(1,ilr), hxh, hyh, hzh,&
            !!     lin%potentialprefac(at%iatype(ilr)),lin%confpotorder, &
            !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3, &
            !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
            call apply_potentialConfinement2(lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
                 Lpot(orbs%ispot(iorb)),epot, rxyz(1,icenter), hxh, hyh, hzh,&
                 lin%potentialprefac(at%iatype(icenter)),lin%confpotorder, &
                 lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3, &
                 lzd%llr(ilr)%bounds%ibyyzz_r) !optional

        else
            call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
                 Lpot(orbs%ispot(iorb)),epot,Lzd%Llr(ilr)%bounds%ibyyzz_r)
        end if
     case('P')
        !here the hybrid BC act the same way
        if(withConfinement) stop 'withConfinement and periodic not possible at the same time!'
        call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,0,0,0,0,orbs%nspinor,npot,psir,&
             Lpot(orbs%ispot(iorb)),epot)
     case('S')
        if(withConfinement) stop 'withConfinement and surface not possible at the same time!'
        call apply_potential(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,0,1,0,0,orbs%nspinor,npot,psir,&
             Lpot(orbs%ispot(iorb)),epot)
     end select
     !k-point values, if present
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))

     if (exctXcoeff /= 0.0_gp) then
        ispot = istexct + Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*(iorb-1)
        !add to the psir function the part of the potential coming from the exact exchange
        call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,exctXcoeff,Lpot(ispot),1,psir(1,1),1)
     end if

     !apply the kinetic term, sum with the potential and transform back to Daubechies basis
     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,Lzd%Llr(ilr),wrk_lh,&
          psir,hpsi(1+oidx),ekin)

     ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot


     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     call deallocate_work_arrays_locham(Lzd%Llr(ilr),wrk_lh)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

     ilrold=ilr
  enddo


END SUBROUTINE local_hamiltonian3

subroutine full_local_potential2(iproc,nproc,ndimpot,ndimgrid,ndimrhopot,nspin,orbs,lzd,ngatherarr,potential,Lpot,flag,comgp)
  use module_base
  use module_types
  use module_xc
  implicit none
  integer, intent(in) :: iproc,nproc,ndimpot,ndimgrid, flag,nspin,ndimrhopot
  type(orbitals_data),intent(inout):: orbs
  type(local_zone_descriptors),intent(inout):: lzd
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(wp), dimension(ndimrhopot), intent(in), target :: potential
  real(wp), dimension(:), pointer, intent(out) :: Lpot
  type(p2pCommsGatherPot),intent(inout), optional:: comgp
  !local variables
  character(len=*), parameter :: subname='full_local_potential'
  logical :: exctX, newvalue
  integer :: npot,ispot,ispotential,ispin,ierr,i_stat, iilr, ilr, ii, iorb, iorb2, nilr, i_all
  integer:: istl, ist, size_Lpot, i3s, i3e
  integer,dimension(:,:),allocatable:: ilrtable
  real(8),dimension(:),pointer:: pot

  call timing(iproc,'Rho_commun    ','ON')

  exctX = xc_exctXfac() /= 0.0_gp

!#########################################################################################################################
! Build the potential on the whole simulation box
! NOTE: in the linear scaling case this should be done for a given localisation
!       region this routine should then be modified or integrated in HamiltonianApplication
! WARNING : orbs%nspin and nspin are not the same !! Check if orbs%nspin should be replaced everywhere
!###########################################################################################################################
  if (flag<2) then
      if (exctX) then
         npot=ndimgrid*nspin+max(max(ndimgrid*orbs%norbp,ngatherarr(0,1)*orbs%norb),1) !part which refers to exact exchange
      else
         npot=ndimgrid*nspin
      end if
      if (nproc > 1) then
         allocate(pot(npot+ndebug),stat=i_stat)
         call memocc(i_stat,pot,'pot',subname)
!         call razero(npot,pot)
         ispot=1
         ispotential=1
         do ispin=1,nspin
            call MPI_ALLGATHERV(potential(ispotential),ndimpot,&
                 mpidtypw,pot(ispot),ngatherarr(0,1),&
                 ngatherarr(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
            ispot=ispot+ndimgrid
            ispotential=ispotential+max(1,ndimpot)
         end do
      else
         if (exctX) then
            allocate(pot(npot+ndebug),stat=i_stat)
            call memocc(i_stat,pot,'pot',subname)
            call dcopy(ndimgrid*nspin,potential,1,pot,1)
         else
            pot => potential
         end if
         ispot=ndimgrid*nspin+1
      end if
  else
      call gatherPotential(iproc, nproc, comgp)
  end if

!#####################################################################################################################
! Determine the dimension of the potential array and orbs%ispot
!#####################################################################################################################
  if(associated(orbs%ispot)) then
     nullify(orbs%ispot)
!     i_all=-product(shape(orbs%ispot))*kind(orbs%ispot)
!     deallocate(orbs%ispot,stat=i_stat)
!     call memocc(i_stat,i_all,'orbs%ispot',subname)
  end if
  allocate(orbs%ispot(orbs%norbp),stat=i_stat)
  call memocc(i_stat,orbs%ispot,'orbs%ispot',subname)

  if(Lzd%nlr > 1) then
     print *,'Using linear scaling'
     allocate(ilrtable(orbs%norbp,2),stat=i_stat)
     call memocc(i_stat,ilrtable,'ilrtable',subname)
     !call to_zero(orbs%norbp*2,ilrtable(1,1))
     ilrtable=0
     ii=0
     do iorb=1,orbs%norbp
        newvalue=.true.
        !localization region to which the orbital belongs
        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
        !spin state of the orbital
        if (orbs%spinsgn(orbs%isorb+iorb) > 0.0_gp) then
           ispin = 1       
        else
           ispin=2
        end if
        !check if the orbitals already visited have the same conditions
        loop_iorb2: do iorb2=1,orbs%norbp
           if(ilrtable(iorb2,1) == ilr .and. ilrtable(iorb2,2)==ispin) then
              newvalue=.false.
              exit loop_iorb2
           end if
        end do loop_iorb2
        if (newvalue) then
          ii = ii + 1
          ilrtable(ii,1)=ilr
          ilrtable(ii,2)=ispin    !SOMETHING IS NOT WORKING IN THE CONCEPT HERE... ispin is not a property of the locregs, but of the orbitals
        end if
     end do
     !number of inequivalent potential regions
     nilr = ii
  else 
     allocate(ilrtable(1,2),stat=i_stat)
     call memocc(i_stat,ilrtable,'ilrtable',subname)
     nilr = 1
     ilrtable=1
  end if 

  lzd%ndimpotisf=0
  do iilr=1,nilr
     ilr=ilrtable(iilr,1)
     do iorb=1,orbs%norbp
       !put the starting point
       if (orbs%inWhichLocreg(iorb+orbs%isorb) == ilr) then
            !assignment of ispot array to the value of the starting address of inequivalent
            orbs%ispot(iorb)=lzd%ndimpotisf + 1
            if(orbs%spinsgn(orbs%isorb+iorb) <= 0.0_gp) then
               orbs%ispot(iorb)=lzd%ndimpotisf + 1 + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
            end if
       end if
     end do
       lzd%ndimpotisf = lzd%ndimpotisf + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i*nspin
  end do 

  if (exctX) then
     lzd%ndimpotisf = lzd%ndimpotisf + max(max(ndimgrid*orbs%norbp,ngatherarr(0,1)*orbs%norb),1) !part which refers to exact exchange
  end if

!#################################################################################################################################################
! Depending on the scheme, cut out the local pieces of the potential
!#################################################################################################################################################
  if(flag==0) then
!       allocate(Lpot(lzd%ndimpotisf+ndebug),stat=i_stat)
!       call dcopy(lzd%ndimpotisf,pot,1,Lpot,1) 
       Lpot=>pot
  else if(flag<2 .and. flag>0) then
     allocate(Lpot(lzd%ndimpotisf+ndebug),stat=i_stat)
     call memocc(i_stat,Lpot,'Lpot',subname)
     ! Cut potential
     istl=1
     do iorb=1,nilr
        ilr = ilrtable(iorb,1)
        
        ! Cut the potential into locreg pieces
        call global_to_local(Lzd%Glr,Lzd%Llr(ilr),orbs%nspin,npot,lzd%ndimpotisf,pot,Lpot(istl))
        istl = istl + Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin
     end do
  else
     allocate(Lpot(lzd%ndimpotisf+ndebug),stat=i_stat)
     call memocc(i_stat,Lpot,'Lpot',subname)
     ist=1
     do iorb=1,nilr
         ilr = ilrtable(iorb,1)
         !determine the dimension of the potential array (copied from full_local_potential)
         if (exctX) then
            stop 'exctX not yet implemented!'
         else
            size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
         end if

         ! Extract the part of the potential which is needed for the current localization region.
         i3s=lzd%Llr(ilr)%nsi3-comgp%ise3(1,iproc)+2 ! starting index of localized  potential with respect to total potential in comgp%recvBuf
         i3e=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i-comgp%ise3(1,iproc)+1 ! ending index of localized potential with respect to total potential in comgp%recvBuf
         if(i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i) then
             write(*,'(a,i0,3x,i0)') 'ERROR: i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i',i3e-i3s+1, Lzd%Llr(ilr)%d%n3i
             stop
         end if

         call global_to_local_parallel(lzd%Glr, lzd%Llr(ilr), orbs%nspin, comgp%nrecvBuf, size_Lpot,&
              comgp%recvBuf, Lpot(ist), i3s, i3e)

         ist = ist + size_lpot
     end do
  end if

  i_all=-product(shape(ilrtable))*kind(ilrtable)
  deallocate(ilrtable,stat=i_stat)
  call memocc(i_stat,i_all,'ilrtable',subname)

  ! Deallocate pot.
  if (flag<2 .and. flag>0) then
      if (nproc > 1) then
         i_all=-product(shape(pot))*kind(pot)
         deallocate(pot,stat=i_stat)
         call memocc(i_stat,i_all,'pot',subname)
      else
         if (exctX) then
            i_all=-product(shape(pot))*kind(pot)
            deallocate(pot,stat=i_stat)
            call memocc(i_stat,i_all,'pot',subname)
         else
            nullify(pot)
         end if
      end if
  end if
  call timing(iproc,'Rho_commun    ','OF')

END SUBROUTINE full_local_potential2

