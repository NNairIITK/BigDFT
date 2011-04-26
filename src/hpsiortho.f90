!> @file
!!  Application of the Hamiltonian + orthonormalize constraints
!! @author
!!    Copyright (C) 2007-2011 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Application of the Hamiltonian
subroutine HamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
     nlpspd,proj,lr,ngatherarr,pot,psi,hpsi,&
     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel,orbsocc,psirocc)
  use module_base
  use module_types
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: lr 
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(in) :: psi
  real(wp), dimension(:), pointer :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
  real(wp), target, dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  real(dp), dimension(*), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  !local variables
  real(gp), dimension(2,orbs%norbp) :: ekin
  real(gp), dimension(2,orbs%norbp) :: epot
  real(wp), dimension(:), pointer :: hpsi2
  character(len=*), parameter :: subname='HamiltonianApplication'
  logical :: exctX,op2p
  integer :: i_all,i_stat,ierr,iorb,n3p,ispot,istart_c,iat
  integer :: istart_ck,isorb,ieorb,ikpt,ispsi_k,nspinor,ispsi
!OCL  integer, dimension(3) :: periodic
!OCL  real(wp) :: maxdiff
!OCL  real(gp) :: eproj,ek_fake,ep_fake
  real(gp), dimension(3,2) :: wrkallred
!OCL  real(wp), dimension(:), allocatable :: hpsi_OCL

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  !check if the potential has been associated
  if (.not. associated(pot)) then
     if (iproc ==0) then
        write(*,*)' ERROR, HamiltonianApplication, potential not associated!'
        stop
     end if
  end if

  !initialise exact exchange energy 
  op2p=(eexctX == -99.0_gp)
  eexctX=0.0_gp

  exctX = libxc_functionals_exctXfac() /= 0.0_gp

  ispot=lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin+1

  !fill the rest of the potential with the exact-exchange terms
  if (present(pkernel) .and. exctX) then
     n3p=ngatherarr(iproc,1)/(lr%d%n1i*lr%d%n2i)
     !exact exchange for virtual orbitals (needs psirocc)

     !here we have to add the round part
     if (present(psirocc) .and. present(orbsocc)) then
        call exact_exchange_potential_virt(iproc,nproc,at%geocode,nspin,&
             lr,orbsocc,orbs,ngatherarr(0,1),n3p,&
             0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psirocc,psi,pot(ispot))
        eexctX = 0._gp
     else
!!$        call exact_exchange_potential_round(iproc,nproc,at%geocode,nspin,lr,orbs,&
!!$             0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi,pot(ispot),eexctX)

        !here the condition for the scheme should be chosen
        if (.not. op2p) then
           call exact_exchange_potential(iproc,nproc,at%geocode,nspin,&
                lr,orbs,ngatherarr(0,1),n3p,&
                0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi,pot(ispot),eexctX)
        else
           !the psi should be transformed in real space
           call exact_exchange_potential_round(iproc,nproc,at%geocode,nspin,lr,orbs,&
                0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi,pot(ispot),eexctX)

        end if
     end if
  else
     eexctX = 0._gp
     !print *,'iproc,eexctX',iproc,eexctX
  end if

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
    allocate(hpsi2((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),stat=i_stat)
    call memocc(i_stat,hpsi2,'hpsi2',subname)
    hpsi(:)=0.0
  else
    hpsi2 => hpsi
  end if
  if (GPUconv) then
     call local_hamiltonian_GPU(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum,GPU)
  else if (OCLconv) then
     call local_hamiltonian_OCL(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi2,ekin_sum,epot_sum,GPU,ekin,epot)
  else
     call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum)
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

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')

  !here the localisation region should be changed, temporary only for cubic approach
  eproj_sum=0.0_gp
  !apply the projectors following the strategy (On-the-fly calculation or not)
  if (DistProjApply) then
     call applyprojectorsonthefly(iproc,orbs,at,lr%d%n1,lr%d%n2,lr%d%n3,&
          rxyz,hx,hy,hz,lr%wfd,nlpspd,proj,psi,hpsi,eproj_sum)
  else if(orbs%norbp > 0) then
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
              call apply_atproj_iorb(iat,iorb,istart_c,at,orbs,lr%wfd,nlpspd,&
                   proj,psi(ispsi),hpsi(ispsi),eproj_sum)           
           end do
           ispsi=ispsi+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*nspinor
        end do
        istart_ck=istart_c
        if (ieorb == orbs%norbp) exit loop_kpt
        ikpt=ikpt+1
        ispsi_k=ispsi
     end do loop_kpt

     if (istart_ck-1 /= nlpspd%nprojel) stop 'incorrect once-and-for-all psp application'
     if (ispsi-1 /= (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'
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

END SUBROUTINE HamiltonianApplication


!> Build the potential in the whole box
subroutine full_local_potential(iproc,nproc,ndimpot,ndimgrid,nspin,norb,norbp,ngatherarr,potential,pot)
  use module_base
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nproc,nspin,ndimpot,norb,norbp,ndimgrid
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(wp), dimension(max(ndimpot,1)*nspin), intent(in), target :: potential
  real(wp), dimension(:), pointer :: pot
  !local variables
  character(len=*), parameter :: subname='full_local_potential'
  logical :: exctX
  integer :: npot,ispot,ispotential,ispin,ierr,i_stat

  call timing(iproc,'Rho_commun    ','ON')
  
  exctX = libxc_functionals_exctXfac() /= 0.0_gp
  !determine the dimension of the potential array
  if (exctX) then
     npot=ndimgrid*nspin+&
          max(max(ndimgrid*norbp,ngatherarr(0,1)*norb),1) !part which refers to exact exchange
  else
     npot=ndimgrid*nspin
  end if

  !build the potential on the whole simulation box
  !in the linear scaling case this should be done for a given localisation region
  !this routine should then be modified or integrated in HamiltonianApplication
  if (nproc > 1) then
     allocate(pot(npot+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)
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

  call timing(iproc,'Rho_commun    ','OF') 

END SUBROUTINE full_local_potential


subroutine free_full_potential(nproc,pot,subname)
  use module_base
  use libxc_functionals
  implicit none
  character(len=*), intent(in) :: subname
  integer, intent(in) :: nproc
  real(wp), dimension(:), pointer :: pot
  !local variables
  logical :: exctX
  integer :: i_all,i_stat

  exctX = libxc_functionals_exctXfac() /= 0.0_gp
  if (nproc > 1 .or. exctX) then
     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot',subname)
  else
     nullify(pot)
  end if

END SUBROUTINE free_full_potential

!> Extract the energy (the quantity which has to be minimised by the wavefunction)
!! and calculate the corresponding gradient.
!! The energy can be the actual Kohn-Sham energy or the trace of the hamiltonian, 
!! depending of the functional we want to calculate. The gradient wrt the wavefucntion
!! Is the put in hpsi accordingly to the functional
subroutine calculate_energy_and_gradient(iter,iproc,nproc,orbs,comms,GPU,lr,hx,hy,hz,ncong,iscf,&
     ekin,epot,eproj,ehart,exc,evxc,eexctX,eion,edisp,psi,psit,hpsi,gnrm,gnrm_zero,energy)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,iscf,iter
  real(gp), intent(in) :: hx,hy,hz,ekin,epot,eproj,ehart,exc,evxc,eexctX,eion,edisp
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(locreg_descriptors), intent(in) :: lr
  type(GPU_pointers), intent(in) :: GPU
  real(gp), intent(out) :: gnrm,gnrm_zero,energy
  real(wp), dimension(:), pointer :: psi,psit,hpsi
  !local variables
  character(len=*), parameter :: subname='calculate_energy_and_gradient' 
  logical :: lcs
  integer :: ierr,ikpt,iorb,i_all,i_stat,k
  real(gp) :: energybs,trH,rzeroorbs,tt,energyKS
  real(wp), dimension(:,:,:), allocatable :: mom_vec

  !band structure energy calculated with occupation numbers
  energybs=ekin+epot+eproj !the potential energy contains also exctX
  !this is the Kohn-Sham energy
  energyKS=energybs-ehart+exc-evxc-eexctX+eion+edisp

  !calculate orbital poloarisation directions
  if(orbs%nspinor==4) then
     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
          lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,psi,mom_vec)
  end if


  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'done,  orthoconstraint...'
  end if

  !transpose the hpsi wavefunction
  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,hpsi,work=psi)

  if (nproc == 1) then
     !associate psit pointer for orthoconstraint and transpose it (for the non-collinear case)
     psit => psi
     call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psit)
  end if

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  !takes also into account parallel k-points distribution
  !here the orthogonality with respect to other occupied functions should be 
  !passed as an optional argument
  call orthoconstraint(iproc,nproc,orbs,comms,lr%wfd,psit,hpsi,trH)

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,hpsi,work=psi)

  !after having calcutated the trace of the hamiltonian, the functional have to be defined
  !new value without the trace, to be added in hpsitopsi
  if (iscf >1) then
     energy=trH
  else
     energy=trH-ehart+exc-evxc-eexctX+eion+edisp
  end if

  !check that the trace of the hamiltonian is compatible with the 
  !band structure energy 
  !this can be done only if the occupation numbers are all equal
  tt=(energybs-trH)/trH
  if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
       (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
     !write this warning only if the system is closed shell
     call check_closed_shell(orbs,lcs)
     if (lcs) then
        write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energybs,trH
     end if
  endif

  call timing(iproc,'Precondition  ','ON')
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'done,  preconditioning...'
  end if

  !Preconditions all orbitals belonging to iproc
  !and calculate the partial norm of the residue
  !switch between CPU and GPU treatment
  if (GPUconv) then
     call preconditionall_GPU(iproc,nproc,orbs,lr,hx,hy,hz,ncong,&
          hpsi,gnrm,gnrm_zero,GPU)
  else if (OCLconv) then
     call preconditionall_OCL(iproc,nproc,orbs,lr,hx,hy,hz,ncong,&
          hpsi,gnrm,gnrm_zero,GPU)
  else
     call preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
  end if

  !sum over all the partial residues
  if (nproc > 1) then
     call mpiallred(gnrm,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     call mpiallred(gnrm_zero,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif

  !count the number of orbitals which have zero occupation number
  !weight this with the corresponding k point weight
  rzeroorbs=0.0_gp
  do ikpt=1,orbs%nkpts
     do iorb=1,orbs%norb
        if (orbs%occup(iorb+(ikpt-1)*orbs%norb) == 0.0_gp) then
           rzeroorbs=rzeroorbs+orbs%kwgts(ikpt)
        end if
     end do
  end do
  !commented out, the kwgts sum already to one
  !if (orbs%nkpts > 1) nzeroorbs=nint(real(nzeroorbs,gp)/real(orbs%nkpts,gp))

  gnrm=sqrt(gnrm/(real(orbs%norb,gp)-rzeroorbs))

  if (rzeroorbs /= 0.0_gp) then
     gnrm_zero=sqrt(gnrm_zero/rzeroorbs)
  else
     gnrm_zero=0.0_gp
  end if

  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)')&
          'done.'
  end if
  call timing(iproc,'Precondition  ','OF')

  if (orbs%nspinor == 4) then
     !only the root process has the correct array
     if(iproc==0 .and. verbose > 0) then
        write(*,'(1x,a)')&
             'Magnetic polarization per orbital'
        write(*,'(1x,a)')&
             '  iorb    m_x       m_y       m_z'
        do iorb=1,orbs%norb
           write(*,'(1x,i5,3f10.5)') &
                iorb,(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)
        end do
     end if
     i_all=-product(shape(mom_vec))*kind(mom_vec)
     deallocate(mom_vec,stat=i_stat)
     call memocc(i_stat,i_all,'mom_vec',subname)
  end if


  !write the energy information
  if (iproc == 0) then
     if (verbose > 0 .and. iscf<1) then
        write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
             ekin,epot,eproj
        write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,exc,evxc
     end if
     if (iscf > 1) then
        if (gnrm_zero == 0.0_gp) then
           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter, tr(H),gnrm',iter,trH,gnrm
        else
           write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter, tr(H),gnrm,gnrm_zero',iter,trH,gnrm,gnrm_zero
        end if
     else
        if (gnrm_zero == 0.0_gp) then
           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energyKS,gnrm
        else
           write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter,total energy,gnrm,gnrm_zero',iter,energyKS,gnrm,gnrm_zero
        end if
     end if
  endif

end subroutine calculate_energy_and_gradient

!>   Operations after h|psi> 
!!   (transposition, orthonormalisation, inverse transposition)
subroutine hpsitopsi(iproc,nproc,orbs,lr,comms,iter,diis,idsx,psi,psit,hpsi,nspin,input)
  use module_base
  use module_types
  use module_interfaces, except_this_one_A => hpsitopsi
  implicit none
  integer, intent(in) :: iproc,nproc,idsx,iter,nspin
  type(locreg_descriptors), intent(in) :: lr
  type(communications_arrays), intent(in) :: comms
  type(orbitals_data), intent(in) :: orbs
  type(input_variables), intent(in) :: input
  type(diis_objects), intent(inout) :: diis
  real(wp), dimension(:), pointer :: psi,psit,hpsi
  !local variables
  character(len=*), parameter :: subname='hpsitopsi'
  integer :: ierr,iorb,k,i_stat,i_all,nzeroorbs

  !adjust the save variables for DIIS/SD switch
  if (iter == 1) then
     diis%ids=0
     diis%mids=1
     diis%idiistol=0
  end if
  !update variables at each iteration step
  if (idsx > 0) then
     diis%mids=mod(diis%ids,idsx)+1
     diis%ids=diis%ids+1
  end if

  diis%energy_min=min(diis%energy_min,diis%energy)

  !transpose the hpsi wavefunction
  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,&
       hpsi,work=psi)

  !apply the minimization method (DIIS or steepest descent)
  call timing(iproc,'Diis          ','ON')

  call psimix(iproc,nproc,orbs,comms,diis,hpsi,psit)

  call timing(iproc,'Diis          ','OF')

  if (iproc == 0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'Orthogonalization...'
  end if

  call orthogonalize(iproc,nproc,orbs,comms,lr%wfd,psit,input)

  !       call checkortho_p(iproc,nproc,norb,nvctrp,psit)
  
  call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,&
       psit,work=hpsi,outadd=psi(1))

  if (nproc == 1) then
     nullify(psit)
  end if
  
  if (iproc == 0 .and. verbose > 1) then
     write(*,'(1x,a)')&
          'done.'
  end if

!!$  if(orbs%nspinor==4) then
!!$     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
!!$     call memocc(i_stat,mom_vec,'mom_vec',subname)
!!$
!!$     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
!!$          lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,psi,mom_vec)
!!$     !only the root process has the correct array
!!$     if(iproc==0 .and. verbose > 0) then
!!$        write(*,'(1x,a)')&
!!$             'Magnetic polarization per orbital'
!!$        write(*,'(1x,a)')&
!!$             '  iorb    m_x       m_y       m_z'
!!$        do iorb=1,orbs%norb
!!$           write(*,'(1x,i5,3f10.5)') &
!!$                iorb,(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)
!!$        end do
!!$     end if
!!$
!!$     i_all=-product(shape(mom_vec))*kind(mom_vec)
!!$     deallocate(mom_vec,stat=i_stat)
!!$     call memocc(i_stat,i_all,'mom_vec',subname)
!!$  end if

  call diis_or_sd(iproc,idsx,orbs%nkptsp,diis)

  !previous value already fulfilled
  diis%energy_old=diis%energy

END SUBROUTINE hpsitopsi

!>Choose among the wavefunctions a subset of them
!! Rebuild orbital descriptors for the new space and allocate the psi_as wavefunction
!! By hypothesis the work array is big enough to contain both wavefunctions
subroutine select_active_space(iproc,nproc,orbs,comms,mask_array,Glr,orbs_as,comms_as,psi,psi_as)
  use module_base
  use module_types
  use module_interfaces, except_this_one => select_active_space
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  logical, dimension(orbs%norb*orbs%nkpts), intent(in) :: mask_array
  real(wp), dimension(orbs%npsidim), intent(in) :: psi
  type(orbitals_data), intent(out) :: orbs_as
  type(communications_arrays), intent(out) :: comms_as
  real(wp), dimension(:), pointer :: psi_as
  !local variables
  character(len=*), parameter :: subname='select_active_space'
  integer :: iorb,ikpt,norb_as,norbu_as,norbd_as,icnt,ikptp,ispsi,ispsi_as
  integer :: i_all,i_stat,nvctrp
    
  !count the number of orbitals of the active space
  norbu_as=-1
  norbd_as=-1
  do ikpt=1,orbs%nkpts
     icnt=0
     do iorb=1,orbs%norbu
        if (mask_array(iorb+(ikpt-1)*orbs%norb)) icnt=icnt+1
     end do
     if (norbu_as /= icnt .and. norbu_as /= -1) then
        write(*,*)'ERROR(select_active_space): the mask array should define always the same norbu'
        stop
     end if
     norbu_as=icnt
     icnt=0
     do iorb=orbs%norbu+1,orbs%norbu+orbs%norbd
        if (mask_array(iorb+(ikpt-1)*orbs%norb)) icnt=icnt+1
     end do
     if (norbd_as /= icnt .and. norbd_as /= -1) then
        write(*,*)'ERROR(select_active_space): the mask array should define always the same norbd'
        stop
     end if
     norbd_as=icnt
  end do

  !allocate the descriptors of the active space
  call orbitals_descriptors(iproc,nproc,norbu_as+norbd_as,norbu_as,norbd_as, &
       & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbs_as)
  !allocate communications arrays for virtual orbitals
  call orbitals_communicators(iproc,nproc,Glr,orbs_as,comms_as)  
  !allocate array of the eigenvalues
  allocate(orbs_as%eval(orbs_as%norb*orbs_as%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs_as%eval,'orbs_as%eval',subname)

  !fill the orbitals array with the values and the wavefunction in transposed form
  icnt=0
  do iorb=1,orbs%nkpts*orbs%norb
     if (mask_array(iorb)) then
        icnt=icnt+1
        orbs_as%eval(icnt)=orbs%eval(iorb)
     end if
  end do
  if (icnt/=orbs_as%norb*orbs_as%nkpts) stop 'ERROR(select_active_space): icnt/=orbs_as%norb*orbs_as%nkpts'

  allocate(psi_as(orbs_as%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi_as,'psi_as',subname)

  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp
     nvctrp=comms%nvctr_par(iproc,ikptp) 
     !this should be identical in both the distributions
     if (nvctrp /= comms_as%nvctr_par(iproc,ikptp)) then
        write(*,*)'ERROR(select_active_space): the component distrbution is not identical'
        stop
     end if
     
     !put all the orbitals which match the active space
     ispsi=1
     ispsi_as=1
     do iorb=1,orbs%norb
        if (mask_array(iorb+(ikpt-1)*orbs%norb)) then
           call dcopy(nvctrp,psi(ispsi),1,psi_as(ispsi_as),1)
           ispsi_as=ispsi_as+nvctrp*orbs_as%nspinor
        end if
        ispsi=ispsi+nvctrp*orbs%nspinor
     end do
  end do


end subroutine select_active_space



!>   First orthonormalisation
subroutine first_orthon(iproc,nproc,orbs,wfd,comms,psi,hpsi,psit,input)
  use module_base
  use module_types
  use module_interfaces, except_this_one_B => first_orthon
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(communications_arrays), intent(in) :: comms
  type(input_variables):: input
  real(wp), dimension(:) , pointer :: psi,hpsi,psit
  !local variables
  character(len=*), parameter :: subname='first_orthon'
  integer :: i_stat

!!!  if(nspin==4) then
!!!     nspinor=4
!!!  else
!!!     nspinor=1
!!!  end if
  
  if (nproc > 1) then
     !allocate hpsi array (used also as transposed)
     !allocated in the transposed way such as 
     !it can also be used as the transposed hpsi
     allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hpsi,'hpsi',subname)
     !allocate transposed principal wavefunction
     allocate(psit(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psit,'psit',subname)
  else
     psit => psi
  end if

  !to be substituted, must pass the wavefunction descriptors to the routine
  call transpose_v(iproc,nproc,orbs,wfd,comms,psi,&
       work=hpsi,outadd=psit(1))

  call orthogonalize(iproc,nproc,orbs,comms,wfd,psit,input)

  !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

  call untranspose_v(iproc,nproc,orbs,wfd,comms,psit,&
       work=hpsi,outadd=psi(1))

  if (nproc == 1) then
     nullify(psit)
     !allocate hpsi array
     allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hpsi,'hpsi',subname)
  end if

END SUBROUTINE first_orthon


!>   Transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,orbs,wfd,nspin,comms,psi,hpsi,psit,evsum, opt_keeppsit)
  use module_base
  use module_types
  use module_interfaces, except_this_one_C => last_orthon
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  integer, intent(in) :: iproc,nproc,nspin
  real(wp), intent(out) :: evsum
  real(wp), dimension(:) , pointer :: psi,hpsi,psit
  logical, optional :: opt_keeppsit
  !local variables

  logical :: keeppsit

  character(len=*), parameter :: subname='last_orthon'
  logical :: dowrite !write the screen output
  integer :: i_all,i_stat,iorb,jorb,md,ikpt,isorb
  real(gp) :: mpol
  real(wp), dimension(:,:,:), allocatable :: mom_vec

  
  if (present(opt_keeppsit)) then
     keeppsit=opt_keeppsit
  else
     keeppsit=.false.
  end if

  call transpose_v(iproc,nproc,orbs,wfd,comms,&
       hpsi,work=psi)
  if (nproc==1) then
     psit => psi
     call transpose_v(iproc,nproc,orbs,wfd,comms,psit)
  end if

  call subspace_diagonalisation(iproc,nproc,orbs,comms,psit,hpsi,evsum)

  call untranspose_v(iproc,nproc,orbs,wfd,comms,&
       psit,work=hpsi,outadd=psi(1))

  if(.not.  keeppsit) then
     if (nproc > 1  ) then
        i_all=-product(shape(psit))*kind(psit)
        deallocate(psit,stat=i_stat)
        call memocc(i_stat,i_all,'psit',subname)
     else
        nullify(psit)
     end if
     
     i_all=-product(shape(hpsi))*kind(hpsi)
     deallocate(hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi',subname)

  endif
  !for a non-collinear treatment,
  !we add the calculation of the moments for printing their value
  !close to the corresponding eigenvector
  if(orbs%nspinor==4) then
     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,wfd%nvctr_c+7*wfd%nvctr_f,&
          orbs%nspinor,psi,mom_vec)
  end if

  ! Send all eigenvalues to all procs.
  call broadcast_kpt_objects(nproc, orbs%nkpts, orbs%norb, &
       & orbs%eval(1), orbs%ikptproc)

  !print the found eigenvalues
  if (iproc == 0) then
     write(*,'(1x,a)')&
          '--------------------------------------- Kohn-Sham Eigenvalues and Occupation Numbers'
     ! Calculate and print the magnetisation
     if (nspin == 2) then
        mpol = 0._gp
        do ikpt=1,orbs%nkpts
           isorb = (ikpt - 1) * orbs%norb
           do iorb = 1, orbs%norbu
              mpol = mpol + orbs%occup(isorb + iorb) * orbs%kwgts(ikpt)
           end do
           do iorb = orbs%norbu + 1, orbs%norb, 1
              mpol = mpol - orbs%occup(isorb + iorb) * orbs%kwgts(ikpt)
           end do
        end do
        write(*,"(1x,A,f9.6)") "Total magnetisation: ", mpol
     end if
     if (orbs%nspinor ==4) then
        write(*,'(1x,a)')&
             '           Eigenvalue                                      m_x       m_y       m_z'
     end if
     do ikpt=1,orbs%nkpts
        if (orbs%nkpts > 1 .and. orbs%nspinor >= 2) write(*,"(1x,A,I4.4,A,3F12.6)") &
             & "Kpt #", ikpt, " BZ coord. = ", orbs%kpts(:, ikpt)
        isorb = (ikpt - 1) * orbs%norb
        if (nspin==1.or.orbs%nspinor==4) then
           do iorb=1,orbs%norb
              dowrite =(iorb <= 5 .or. iorb >= orbs%norb-5) .or. verbose > 0
              if (orbs%nspinor ==4) then
                 if (dowrite) & 
                      write(*,'(1x,a,i4,a,1x,1pe21.14,1x,0pf6.4,16x,(1x,3(0pf10.5)))') &
                      'e(',iorb,')=',orbs%eval(isorb + iorb),orbs%occup(isorb+iorb),(mom_vec(md,iorb,1)/mom_vec(1,iorb,1),md=2,4)
              else
                 if (dowrite) & 
                      write(*,'(1x,a,i4,a,1x,1pe21.14,1x,0pf6.4)') 'e(',iorb,')=',orbs%eval(isorb + iorb),orbs%occup(isorb+iorb)
              end if
           end do
        else
           do iorb=1,min(orbs%norbu,orbs%norbd)
              jorb=orbs%norbu+iorb
              dowrite =(iorb <= 5 .or. iorb >= min(orbs%norbu,orbs%norbd)-5)  .or. verbose > 0
              if (dowrite) & 
                   write(*,'(1x,a,i4,a,1x,1pe21.14,1x,0pf6.4,6x,0pf6.4,1x,a,i4,a,1x,1pe21.14)') &
                   'e(',iorb,',u)=',orbs%eval(isorb + iorb),orbs%occup(isorb + iorb),&
                   orbs%occup(isorb + jorb),'e(',iorb,',d)=',orbs%eval(isorb + jorb)
           end do
           if (orbs%norbu > orbs%norbd) then
              do iorb=orbs%norbd+1,orbs%norbu
                 dowrite =(iorb <= 5 .or. iorb >= orbs%norbu-5) .or. verbose > 0
                 if (dowrite) & 
                      write(*,'(1x,a,i4,a,1x,1pe21.14,1x,0pf6.4)') 'e(',iorb,',u)=',orbs%eval(isorb + iorb),orbs%occup(isorb+iorb)
              end do
           else if (orbs%norbd > orbs%norbu) then
              do iorb=2*orbs%norbu+1,orbs%norbu+orbs%norbd
                 dowrite =(iorb <= 5 .or. iorb >= orbs%norbd-5) .or. verbose > 0
                 if (dowrite) & 
                      write(*,'(46x,0pf6.4,1x,a,i4,a,1x,1pe21.14)') orbs%occup(isorb + iorb),&
                      'e(',iorb-orbs%norbu,',d)=',orbs%eval(isorb + iorb)
              end do
           end if
        end if
     end do
  end if

  if (orbs%nspinor ==4) then
     i_all=-product(shape(mom_vec))*kind(mom_vec)
     deallocate(mom_vec,stat=i_stat)
     call memocc(i_stat,i_all,'mom_vec',subname)
  end if


END SUBROUTINE last_orthon


!> Finds the fermi level ef for an error function distribution with a width wf
!! eval are the Kohn Sham eigenvalues and melec is the total number of electrons
subroutine evaltoocc(iproc,nproc,filewrite,wf,orbs)
 use module_base
 use module_types
 implicit none
 logical, intent(in) :: filewrite
 integer, intent(in) :: iproc, nproc
 real(gp), intent(in) :: wf
 type(orbitals_data), intent(inout) :: orbs
 !local variables
 integer :: ikpt,iorb,melec,ii
 real(gp) :: charge, chargef
 real(gp) :: ef,pi,electrons,dlectrons,factor,arg,argu,argd,corr,cutoffu,cutoffd,diff,full,res,resu,resd
 parameter(pi=3.1415926535897932d0)
 !write(*,*)  'ENTER Fermilevel',orbs%norbu,orbs%norbd
 
 if (orbs%norbd==0) then 
    full=2.d0   ! maximum occupation for closed shell  orbital
 else
    full=1.d0   ! maximum occupation for spin polarized orbital
 endif
 
 if (orbs%nkpts.ne.1 .and. filewrite) stop 'Fermilevel: CANNOT write input.occ with more than one k-point'
 charge=0.0_gp
 do ikpt=1,orbs%nkpts
    !number of zero orbitals for the given k-point
    !overall charge of the system
    do iorb=1,orbs%norb
       charge=charge+orbs%occup(iorb+(ikpt-1)*orbs%norb) * orbs%kwgts(ikpt)
    end do
 end do
 melec=nint(charge)
 !if (iproc == 0) write(*,*) 'charge',charge,melec

 ! Send all eigenvalues to all procs.
 call broadcast_kpt_objects(nproc, orbs%nkpts, orbs%norb, &
      & orbs%eval(1), orbs%ikptproc)

 if (wf > 0.0_gp) then
    ii=0
    if (orbs%efermi == UNINITIALISED) then
       ! Take initial value at gamma point.
       do iorb = 1, orbs%norbu
          if (orbs%occup(iorb) < 1.0_gp) then
             orbs%efermi = orbs%eval(iorb)
             exit
          end if
       end do
    end if
    ef=orbs%efermi
    factor=1.d0/(sqrt(pi)*wf)
    !print *,0,ef
    loop_fermi: do
       ii=ii+1
       if (ii > 10000) stop 'error Fermilevel'
       electrons=0.d0
       dlectrons=0.d0
       do ikpt=1,orbs%nkpts
          do iorb=1,orbs%norbd+orbs%norbu
             arg=(orbs%eval((ikpt-1)*orbs%norb+iorb)-ef)/wf
             if (occopt == SMEARING_DIST_ERF) then
                ! next 2 line error function distribution
                call derf_ab(res,arg)
                electrons=electrons+.5d0*(1.d0-res) * orbs%kwgts(ikpt)
                dlectrons=dlectrons-exp(-arg**2) * orbs%kwgts(ikpt)
             else if (occopt == SMEARING_DIST_FERMI) then
                !print *,iorb,ef,orbs%eval((ikpt-1)*orbs%norb+iorb),arg,electrons
             else if (occopt == SMEARING_DIST_FERMI) then
                !! next 2 line Fermi function distribution
                electrons=electrons+1.d0/(1.d0+exp(arg)) * orbs%kwgts(ikpt)
                dlectrons=dlectrons-1.d0/(2.d0+exp(arg)+exp(-arg)) * orbs%kwgts(ikpt)
             end if
          enddo
       enddo

       if (occopt == SMEARING_DIST_ERF) then
          ! next  line error function distribution
          dlectrons=dlectrons*factor
       else if (occopt == SMEARING_DIST_FERMI) then
          ! next  line Fermi function distribution
          dlectrons=dlectrons/wf
       end if
       
       diff=real(melec,gp)/full-electrons
       if (abs(diff) < 1.d-12) exit loop_fermi
       corr=diff/dlectrons
       !if (iproc==0) write(*,*) ii,electrons,ef,dlectrons,melec,corr
       if (corr > 1.d0*wf) corr=1.d0*wf
       if (corr < -1.d0*wf) corr=-1.d0*wf
       if (abs(dlectrons) < 1.d-18  .and. electrons > real(melec,gp)/full) corr=3.d0*wf
       if (abs(dlectrons) < 1.d-18  .and. electrons < real(melec,gp)/full) corr=-3.d0*wf
       ef=ef-corr
    end do loop_fermi
    
    do ikpt=1,orbs%nkpts
       argu=(orbs%eval((ikpt-1)*orbs%norb+orbs%norbu)-ef)/wf
       argd=(orbs%eval((ikpt-1)*orbs%norb+orbs%norbu+orbs%norbd)-ef)/wf
       if (occopt == SMEARING_DIST_ERF) then
          !error function
          call derf_ab(resu,argu)
          call derf_ab(resd,argd)
          cutoffu=.5d0*(1.d0-resu)
          cutoffd=.5d0*(1.d0-resd)
       else if (occopt == SMEARING_DIST_FERMI) then
          !Fermi function
          cutoffu=1.d0/(1.d0+exp(argu))
          cutoffd=1.d0/(1.d0+exp(argd))
       end if
    enddo
    if (iproc==0) write(*,'(1x,a,1pe21.14,2(1x,e8.1))') 'Fermi level, Fermi distribution cut off at:  ',ef,cutoffu,cutoffd
    orbs%efermi=ef
    
    !update the occupation number
    do ikpt=1,orbs%nkpts
       do iorb=1,orbs%norbu + orbs%norbd
          arg=(orbs%eval((ikpt-1)*orbs%norb+iorb)-ef)/wf
          if (occopt == SMEARING_DIST_ERF) then
             !error function
             call derf_ab(res,arg)
             orbs%occup((ikpt-1)*orbs%norb+iorb)=full*.5d0*(1.d0-res)
             !print *,'iorb,arg,res,full*.5d0*(1.d0-res)',iorb,arg,res,full*.5d0*(1.d0-res)
          else if (occopt == SMEARING_DIST_FERMI) then
             !Fermi function
             orbs%occup((ikpt-1)*orbs%norb+iorb)=full*1.d0/(1.d0+exp(arg))
          end if
       end do
    end do
    ! Sanity check on sum of occup.
    chargef=0.0_gp
    do ikpt=1,orbs%nkpts
       do iorb=1,orbs%norb
          chargef=chargef+orbs%kwgts(ikpt) * orbs%occup(iorb+(ikpt-1)*orbs%norb)
       end do
    end do
    if (abs(charge - chargef) > 1e-6)  stop 'error occupation update'
 else if(full==1.0_gp) then
    call eFermi_nosmearing(iproc,orbs)
 end if
 
 !write on file the results if needed
 if (filewrite) then
    open(unit=11,file='input.occ',status='unknown')
    write(11,*)orbs%norbu,orbs%norbd
    do iorb=1,orbs%norb
       write(11,'(i5,e19.12)')iorb,orbs%occup((ikpt-1)*orbs%norb+iorb)
       !    write(11,'(i5,e19.12)')iorb,full/(1.d0+exp(arg))  !,orbs%eval((ikpt-1)*orbs%norb+iorb)
    end do
    close(unit=11)
 end if

END SUBROUTINE evaltoocc


subroutine eFermi_nosmearing(iproc,orbs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  type(orbitals_data), intent(inout) :: orbs
  !local variables
  integer :: iu,id,n,nzeroorbs,ikpt,iorb
  real(gp) :: charge
  real(wp) :: eF

  iu=0
  id=0
  eF = 0._wp
  do ikpt=1,orbs%nkpts
     !number of zero orbitals for the given k-point
     nzeroorbs=0
     !overall charge of the system
     charge=0.0_gp
     do iorb=1,orbs%norb
        if (orbs%occup(iorb+(ikpt-1)*orbs%norb) == 0.0_gp) then
           nzeroorbs=nzeroorbs+1
        else
           charge=charge+orbs%occup(iorb+(ikpt-1)*orbs%norb)
        end if
     end do
     if (nzeroorbs /= 0 .and. orbs%norbd .gt.0) then
        do iorb=1,orbs%norbu-1
           if (orbs%eval((ikpt-1)*orbs%norb+iorb) > orbs%eval((ikpt-1)*orbs%norb+iorb+1)) &
                write(*,*) 'wrong ordering of up EVs',iorb,iorb+1
        end do
        do iorb=1,orbs%norbd-1
           if (orbs%eval((ikpt-1)*orbs%norb+iorb+orbs%norbu) > orbs%eval((ikpt-1)*orbs%norb+iorb+1+orbs%norbu))&
                write(*,*) 'wrong ordering of dw EVs',iorb+orbs%norbu,iorb+1+orbs%norbu
        enddo

        iu=0
        id=0
        n=0
        do while (real(n,gp) < charge)
           if (orbs%eval((ikpt-1)*orbs%norb+iu+1) <= orbs%eval((ikpt-1)*orbs%norb+id+1+orbs%norbu)) then
              iu=iu+1
              eF=orbs%eval((ikpt-1)*orbs%norb+iu+1)
           else
              id=id+1
              eF=orbs%eval((ikpt-1)*orbs%norb+id+1+orbs%norbu)
           endif
           n=n+1
        enddo
        if (iproc==0) write(*,'(1x,a,1pe21.14,a,i4)') 'Suggested Homo energy level',eF,', Spin polarization',iu-id
        !write(*,*) 'up,down, up-down',iu,id,iu-id
     end if
  end do
  orbs%efermi=eF
  !assign the values for the occupation numbers
  do iorb=1,iu
     orbs%occup(iorb)=1.0_gp
  end do
  do iorb=iu+1,orbs%norbu
     orbs%occup(iorb)=0.0_gp
  end do
  do iorb=1,id
     orbs%occup(iorb+orbs%norbu)=1.0_gp
  end do
  do iorb=id+1,orbs%norbd
     orbs%occup(iorb+orbs%norbu)=0.0_gp
  end do

END SUBROUTINE eFermi_nosmearing


!>   Calculate magnetic moments
subroutine calc_moments(iproc,nproc,norb,norb_par,nvctr,nspinor,psi,mom_vec)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctr,nspinor
  integer, dimension(0:nproc-1), intent(in) :: norb_par
  real(wp), dimension(nvctr,norb*nspinor), intent(in) :: psi
  real(wp), dimension(4,norb,min(nproc,2)), intent(out) :: mom_vec
  !local variables
  character(len=*), parameter :: subname='calc_moments'
  integer :: i_all,i_stat,ierr,iorb,jproc
  integer :: ndim,oidx
  integer, dimension(:), allocatable :: norb_displ
  real(wp) :: m00,m11,m13,m24,m14,m23
  !real(wp), dimension(:,:,:), allocatable :: mom_vec

  ndim=2
  if (nproc==1) ndim=1

  if(nspinor==4) then
     
     call razero(4*norb*ndim,mom_vec)
     
     do iorb=1,norb_par(iproc)
        oidx=(iorb-1)*nspinor+1
        m00=dot(2*nvctr,psi(1,oidx),1,psi(1,oidx),1)
        m11=dot(2*nvctr,psi(1,oidx+2),1,psi(1,oidx+2),1)
        m13=dot(nvctr,psi(1,oidx),1,psi(1,oidx+2),1)
        m24=dot(nvctr,psi(1,oidx+1),1,psi(1,oidx+3),1)
        !        m12=dot(nvctr,psi(1,oidx),1,psi(1,oidx+1),1)
        !        m34=dot(nvctr,psi(1,oidx+2),1,psi(1,oidx+3),1)
        m14=dot(nvctr,psi(1,oidx),1,psi(1,oidx+3),1)
        m23=dot(nvctr,psi(1,oidx+1),1,psi(1,oidx+2),1)

        mom_vec(1,iorb,ndim)=(m00+m11) !rho
        mom_vec(2,iorb,ndim)=2.0d0*(m13+m24)       !m_x
        !        mom_vec(3,iorb,ndim)=2.0d0*(m12-m34)       !m_y
        mom_vec(3,iorb,ndim)=2.0d0*(m14-m23)       !m_y
        mom_vec(4,iorb,ndim)=(m00-m11)             !m_z
     end do
     
     if(nproc>1) then
        allocate(norb_displ(0:nproc-1+ndebug),stat=i_stat)
        call memocc(i_stat,norb_displ,'norb_displ',subname)

        norb_displ(0)=0
        do jproc=1,nproc-1
           norb_displ(jproc)=norb_displ(jproc-1)+norb_par(jproc-1)
        end do
        
        call MPI_GATHERV(mom_vec(1,1,2),4*norb_par(iproc),mpidtypw,&
             mom_vec(1,1,1),4*norb_par,4*norb_displ,mpidtypw,&
             0,MPI_COMM_WORLD,ierr)

        i_all=-product(shape(norb_displ))*kind(norb_displ)
        deallocate(norb_displ,stat=i_stat)
        call memocc(i_stat,i_all,'norb_displ',subname)
     end if
    
  end if

END SUBROUTINE calc_moments


subroutine check_communications(iproc,nproc,orbs,lr,comms)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  type(communications_arrays), intent(in) :: comms
  !local variables
  character(len=*), parameter :: subname='check_communications'
  integer :: i,ispinor,iorb,indspin,indorb,jproc,i_stat,i_all,iscomp,idsx,index,ikptsp
  integer :: ikpt,ispsi,nspinor,nvctrp
  real(wp) :: psival,maxdiff,ierr
  real(wp), dimension(:), allocatable :: psi
  real(wp), dimension(:), pointer :: pwork
  real(wp) :: epsilon
  character(len = 25) :: filename
  logical :: abort

  !allocate the "wavefunction" amd fill it, and also the workspace
  allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(pwork(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,pwork,'pwork',subname)

  do iorb=1,orbs%norbp
     ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
     !valkpt=real(512*ikpt,wp)
     !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
     indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
     do ispinor=1,orbs%nspinor
        indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
        do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !vali=real(i,wp)/512.0_wp  ! *1.d-5
           call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
           psi(i+indspin+indorb)=psival!(valorb+vali)*(-1)**(ispinor-1)
        end do
     end do
  end do
 
  !transpose the hpsi wavefunction
  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psi,work=pwork)

  !check the results of the transposed wavefunction
  maxdiff=0.0_wp
  ispsi=0
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
     !valkpt=real(512*ikpt,wp)
     !calculate the starting point for the component distribution
     iscomp=0
     do jproc=0,iproc-1
        iscomp=iscomp+comms%nvctr_par(jproc,ikptsp)
     end do
     nvctrp=comms%nvctr_par(iproc,ikptsp)
     nspinor=orbs%nspinor

     do iorb=1,orbs%norb
        !valorb=real(iorb,wp)+valkpt
        indorb=(iorb-1)*nvctrp*nspinor
        do idsx=1,(nspinor-1)/2+1
           do i=1,nvctrp
              !vali=real(i+iscomp,wp)/512.d0  ! *1.d-5
              do ispinor=1,((2+nspinor)/4+1)
                 !psival=(-1)**(ispinor-1)*(valorb+vali)
                 call test_value(ikpt,iorb,ispinor,i+iscomp,psival)
                 !this is just to force the IEEE representation of psival
                 !              if (psival .lt. 0.d0) then  
                 !              write(321,*) psival,psival**2
                 !              endif
                 index=ispinor+(i-1)*((2+nspinor)/4+1)+&
                      (idsx-1)*((2+nspinor)/4+1)*nvctrp+indorb+ispsi
                 maxdiff=max(abs(psi(index)-psival),maxdiff)
              end do
           end do
        end do
     end do
     ispsi=ispsi+nvctrp*orbs%norb*nspinor
  end do

  abort = .false.
  if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
     write(*,*)'ERROR: process',iproc,'does not transpose wavefunctions correctly!'
     write(*,*)'       found an error of',maxdiff,'cannot continue.'
     write(*,*)'       data are written in the file transerror.log, exiting...'

     write(filename, "(A,I0,A)") 'transerror', iproc, '.log'
     open(unit=22,file=trim(filename),status='unknown')
     ispsi=0
     do ikptsp=1,orbs%nkptsp
        ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
        !valkpt=real(512*ikpt,wp)
        !calculate the starting point for the component distribution
        iscomp=0
        do jproc=0,iproc-1
           iscomp=iscomp+comms%nvctr_par(jproc,ikptsp)
        end do
        nvctrp=comms%nvctr_par(iproc,ikptsp)
        nspinor=orbs%nspinor

        do iorb=1,orbs%norb
           !valorb=real(iorb,wp)+valkpt
           indorb=(iorb-1)*nvctrp*nspinor
           do idsx=1,(nspinor-1)/2+1
              do i=1,nvctrp
                 !vali=real(i+iscomp,wp)/512.d0  !*1.d-5
                 do ispinor=1,((2+nspinor)/4+1)
                    !psival=(-1)**(ispinor-1)*(valorb+vali)
                    call test_value(ikpt,iorb,ispinor,i+iscomp,psival)
                    index=ispinor+(i-1)*((2+nspinor)/4+1)+&
                         (idsx-1)*((2+nspinor)/4+1)*nvctrp+indorb+ispsi
                    maxdiff=abs(psi(index)-psival)
                    if (maxdiff > 0.d0) then
                       write(22,'(i3,i6,2i4,3(1x,1pe13.6))')ispinor,i+iscomp,iorb,ikpt,psival,&
                            psi(index),maxdiff
                    end if
                 end do
              end do
           end do
        end do
        ispsi=ispsi+nvctrp*orbs%norb*nspinor
     end do
     close(unit=22)
     abort = .true.
     write(filename, "(A,I0,A)") 'distscheme', iproc, '.log'
     open(unit=22,file=trim(filename),status='unknown')
     call print_distribution_schemes(22,nproc,orbs%nkpts,orbs%norb_par,comms%nvctr_par)
     close(unit=22)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  if (abort) then
     if (iproc == 0) call print_distribution_schemes(6,nproc,orbs%nkpts,orbs%norb_par,comms%nvctr_par)
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
  end if

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,&
       psi,work=pwork)

  maxdiff=0.0_wp
  do iorb=1,orbs%norbp
     ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
     !valkpt=real(512*ikpt,wp)
     !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
     indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
     do ispinor=1,orbs%nspinor
        indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
        do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !vali=real(i,wp)/512.d0  !*1.d-5
           !psival=(valorb+vali)*(-1)**(ispinor-1)
           call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
           maxdiff=max(abs(psi(i+indspin+indorb)-psival),maxdiff)
        end do
     end do
  end do

  abort = .false.
  if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
     write(*,*)'ERROR: process',iproc,'does not untranspose wavefunctions correctly!'
     write(*,*)'       found an error of',maxdiff,'cannot continue.'
     write(*,*)'       data are written in the file transerror.log, exiting...'

     write(filename, "(A,I0,A)") 'transerror', iproc, '.log'
     open(unit=22,file=trim(filename),status='unknown')
     maxdiff=0.0_wp
     do iorb=1,orbs%norbp
        ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
        !valkpt=real(512*ikpt,wp)
        !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
        indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
        do ispinor=1,orbs%nspinor
           indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
           do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
              !vali=real(i,wp)/512.d0  !*1.d-5
              !psival=(valorb+vali)*(-1)**(ispinor-1)
              call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
              maxdiff=abs(psi(i+indspin+indorb)-psival)
              if (maxdiff > 0.d0) then
                 write(22,'(i3,i6,2i4,3(1x,1pe13.6))')ispinor,i,iorb,orbs%isorb,psival,&
                      psi(ispinor+(i-1)*orbs%nspinor+indorb),maxdiff
              end if
           end do
        end do
     end do
     close(unit=22)
     abort = .true.
  end if

  if (abort) call MPI_ABORT(MPI_COMM_WORLD,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(pwork))*kind(pwork)
  deallocate(pwork,stat=i_stat)
  call memocc(i_stat,i_all,'pwork',subname)



END SUBROUTINE check_communications


!> define a value for the wavefunction which is dependent of the indices
subroutine test_value(ikpt,iorb,ispinor,icomp,val)
  use module_base
  implicit none
  integer, intent(in) :: ikpt,icomp,iorb,ispinor
  real(wp), intent(out) :: val
  !local variables
  real(wp) :: valkpt,valorb,vali

! recognizable pattern, for debugging
! valkpt=real(10000*(ikpt-1),wp)!real(512*ikpt,wp)
! valorb=real(iorb,wp)+valkpt
! vali=real(icomp,wp)*1.e-5_wp  !real(icomp,wp)/512.0_wp  ! *1.d-5
!
! val=(valorb+vali)*(-1)**(ispinor-1)

  valkpt=real(512*ikpt,wp)
  valorb=real(iorb,wp)+valkpt
  vali=real(icomp,wp)/512.0_wp  ! *1.d-5

  val=(valorb+vali)*(-1)**(ispinor-1)

  
END SUBROUTINE test_value
  

subroutine broadcast_kpt_objects(nproc, nkpts, ndata, data, ikptproc)
  use module_base
  implicit none
  integer, intent(in) :: nproc, nkpts, ndata
  integer, intent(in) :: ikptproc(nkpts)
  real(gp), intent(inout) :: data(ndata, nkpts)

  integer :: ikpt, ierr

  if (nproc > 1) then
     do ikpt = 1, nkpts
        call MPI_BCAST(data(1,ikpt), ndata, MPI_DOUBLE_PRECISION, &
             & ikptproc(ikpt), MPI_COMM_WORLD, ierr)
     end do
  end if
END SUBROUTINE broadcast_kpt_objects
