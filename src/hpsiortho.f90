!!****f* BigDFT/HamiltonianApplication
!! FUNCTION
!!  Application of the Hamiltonian
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine HamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
     cpmult,fpmult,radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,psi,hpsi,&
     ekin_sum,epot_sum,eproj_sum,nspin,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,ndimpot,nspin
  real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: lr 
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf  
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
  real(wp), dimension(max(ndimpot,1),nspin), intent(in), target :: potential
  real(gp), intent(out) :: ekin_sum,epot_sum,eproj_sum
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  !local variables
  character(len=*), parameter :: subname='HamiltonianApplication'
  integer :: i_all,i_stat,ierr,iorb,ispin
  real(gp) :: eproj
  real(gp), dimension(3,2) :: wrkallred
  real(wp), dimension(:,:), pointer :: pot
  integer,parameter::lupfil=14

  !stream ptr array
!  real(kind=8), dimension(orbs%norbp) :: tab_stream_ptr
!  real(kind=8) :: stream_ptr_first_trsf

  call timing(iproc,'ApplyLocPotKin','ON')

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  !build the potential on the whole simulation box
  !in the linear scaling case this should be done for a given localisation region
  !cannot be deplaced due to the fact that n1i is not calculated
  if (nproc > 1) then
     allocate(pot(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)

     do ispin=1,nspin
        call MPI_ALLGATHERV(potential(1,ispin),ndimpot,&
             mpidtypw,pot(1,ispin),ngatherarr(0,1),&
             ngatherarr(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
     end do
  else
     pot => potential
  end if

  !apply the local hamiltonian for each of the orbitals
  !given to each processor

  !switch between GPU/CPU treatment
  if (GPUconv) then
     call local_hamiltonian_GPU(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum,GPU)
  else
     call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum)
  end if
  
  if (nproc > 1) then
     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot',subname)
  else
     nullify(pot)
  end if

  call timing(iproc,'ApplyLocPotKin','OF')

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')

  !here the localisation region should be changed, temporary only for cubic approach
  eproj_sum=0.0_gp
  !apply the projectors following the strategy (On-the-fly calculation or not)
  if (DistProjApply) then
     call applyprojectorsonthefly(iproc,orbs,at,lr%d%n1,lr%d%n2,lr%d%n3,&
          rxyz,hx,hy,hz,cpmult,fpmult,radii_cf,lr%wfd,nlpspd,proj,psi,hpsi,eproj_sum)
  else
     !one should add a flag here which states that it works only for global region
     ! loop over all my orbitals
     !should be changed in view of spin-orbit coupling
     do iorb=1,orbs%norbp*orbs%nspinor
        call applyprojectorsone(at%ntypes,at%nat,at%iatype,at%psppar,at%npspcode, &
             nlpspd%nprojel,nlpspd%nproj,nlpspd%nseg_p,nlpspd%keyg_p,nlpspd%keyv_p,&
             nlpspd%nvctr_p,&
             proj,lr%wfd%nseg_c,lr%wfd%nseg_f,lr%wfd%keyg,lr%wfd%keyv,&
             lr%wfd%nvctr_c,lr%wfd%nvctr_f, & 
             psi(1,iorb),hpsi(1,iorb),eproj)
        eproj_sum=eproj_sum+&
             orbs%kwgts(orbs%iokpt((iorb-1)/orbs%nspinor+1))*&
             orbs%occup((iorb+orbs%isorb-1)/orbs%nspinor+1)*eproj
     enddo
  end if

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

END SUBROUTINE HamiltonianApplication
!!***


!!****f* BigDFT/hpsitopsi
!! FUNCTION
!!   Operations after h|psi> 
!!   (transposition, orthonormalisation, inverse transposition)
!! SOURCE
!!
subroutine hpsitopsi(iproc,nproc,orbs,hx,hy,hz,lr,comms,&
     ncong,iter,idsx,idsx_actual,ads,energy,energy_old,energy_min,&
     alpha,gnrm,scprsum,psi,psit,hpsi,psidst,hpsidst,nspin,GPU)
  use module_base
  use module_types
  use module_interfaces, except_this_one_A => hpsitopsi
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,idsx,iter,nspin
  real(gp), intent(in) :: hx,hy,hz,energy,energy_old
  type(locreg_descriptors), intent(in) :: lr
  type(communications_arrays), intent(in) :: comms
  type(orbitals_data), intent(in) :: orbs
  integer, intent(inout) :: idsx_actual
  real(wp), intent(inout) :: alpha
  real(dp), intent(inout) :: gnrm,scprsum
  real(gp), intent(inout) :: energy_min
  real(wp), dimension(:), pointer :: psi,psit,hpsi,psidst,hpsidst
  real(wp), dimension(:,:,:), pointer :: ads
  type(GPU_pointers), intent(inout) :: GPU
  !local variables
  character(len=*), parameter :: subname='hpsitopsi'
  logical, save :: switchSD
  integer, save :: idiistol,mids,ids
  integer :: ierr,iorb,k,i_stat,i_all
  real(dp) :: tt,scprpart
  real(wp), dimension(:,:,:), allocatable :: mom_vec

  !stream ptr array
 ! real(kind=8), dimension(orbs%norbp) :: tab_stream_ptr

  !adjust the save variables for DIIS/SD switch
  if (iter == 1) then
     !logical control variable for switch DIIS-SD
     switchSD=.false.
     ids=0
     mids=1
     idiistol=0
  end if
  !update variables at each iteration step
  if (idsx > 0) then
     mids=mod(ids,idsx)+1
     ids=ids+1
  end if

  energy_min=min(energy_min,energy)

  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'done, orthoconstraint...'
  end if

  !transpose the hpsi wavefunction
  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,hpsi,work=psi)

  if (nproc == 1) then
     !associate psit pointer for orthoconstraint and transpose it (for the non-collinear case)
     psit => psi
     call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psit)
  end if

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  call orthoconstraint(iproc,nproc,orbs,comms,lr%wfd,psit,hpsi,scprsum)

!!$  call orthoconstraint_p(iproc,nproc,orbs%norbu,orbs%occup,comms%nvctr_par(iproc,1),psit,hpsi,&
!!$       scprsum,orbs%nspinor)
!!$  scprpart=0.0_dp
!!$  if(orbs%norbd > 0) then
!!$     scprpart=scprsum 
!!$     call orthoconstraint_p(iproc,nproc,orbs%norbd,orbs%occup(orbs%norbu+1),comms%nvctr_par(iproc,1),&
!!$          psit(1+comms%nvctr_par(iproc,1)*orbs%norbu),hpsi(1+comms%nvctr_par(iproc,1)*orbs%norbu),&
!!$          scprsum,orbs%nspinor)
!!$     scprsum=scprsum+scprpart
!!$  end if

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,hpsi,work=psi)

  call timing(iproc,'Precondition  ','ON')
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'done, preconditioning...'
  end if

  !Preconditions all orbitals belonging to iproc
  !and calculate the partial norm of the residue
  !switch between CPU and GPU treatment
  if (GPUconv) then
     call preconditionall_GPU(iproc,nproc,orbs,lr,hx,hy,hz,ncong,&
          hpsi,gnrm,GPU)
  else
     call preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm)
  end if

  !sum over all the partial residues
  if (nproc > 1) then
     tt=gnrm
     call MPI_ALLREDUCE(tt,gnrm,1,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif
  gnrm=sqrt(gnrm/real(orbs%norb,dp))

  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)')&
          'done.'
  end if
  call timing(iproc,'Precondition  ','OF')

  !apply the minimization method (DIIS or steepest descent)
  if (idsx_actual > 0) then
     !transpose the hpsi wavefunction into the diis array
     call transpose_v(iproc,nproc,orbs,lr%wfd,comms,&
          hpsi,work=psi,&
          outadd=hpsidst(1+comms%nvctr_par(iproc,1)*orbs%nspinor*orbs%norb*(mids-1)))

     call timing(iproc,'Diis          ','ON')

     !psidst=psit
     call dcopy(comms%nvctr_par(iproc,1)*orbs%norb*orbs%nspinor,&
          psit(1),1,&
          psidst(1+comms%nvctr_par(iproc,1)*orbs%nspinor*orbs%norb*(mids-1)),1) 

     if (nproc == 1) then
        !hpsidst=hpsi
        call dcopy(comms%nvctr_par(iproc,1)*orbs%norb*orbs%nspinor,&
             hpsi(1),1,&
             hpsidst(1+comms%nvctr_par(iproc,1)*orbs%nspinor*orbs%norb*(mids-1)),1) 
     endif

     call diisstp(orbs%norb,nproc,iproc,orbs%nspinor,  &
          ads,ids,mids,idsx_actual,comms%nvctr_par(iproc,1),&
          psit,psidst,hpsidst)
  else
     ! update all wavefunctions with the preconditioned gradient
     if (energy > energy_old) then
        alpha=max(.125_wp,.5_wp*alpha)
        if (alpha == .125_wp) write(*,*) 'Convergence problem or limit'
     else
        alpha=min(1.05_wp*alpha,1._wp)
     endif
     if (iproc == 0 .and. verbose > 0) write(*,'(1x,a,1pe11.3)') 'alpha=',alpha

     !transpose the hpsi wavefunction
     call transpose_v(iproc,nproc,orbs,lr%wfd,comms,&
          hpsi,work=psi)

     call timing(iproc,'Diis          ','ON')

!!$     do iorb=1,orbs%norb*orbs%nspinor
!!$        call axpy(comms%nvctr_par(iproc),&
!!$             -alpha,hpsi(1+comms%nvctr_par(iproc)*(iorb-1)),1,&
!!$             psit(1+comms%nvctr_par(iproc)*(iorb-1)),1)
!!$     enddo

     call axpy(comms%nvctr_par(iproc,1)*orbs%norb*orbs%nspinor,&
          -alpha,hpsi(1),1,psit(1),1)

  endif

  call timing(iproc,'Diis          ','OF')

  if (iproc == 0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'Orthogonalization...'
  end if

  call orthogonalize(iproc,nproc,orbs,comms,lr%wfd,psit)

!!$  call orthon_p(iproc,nproc,orbs%norbu,comms%nvctr_par(iproc,1),lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,&
!!$       psit,orbs%nspinor)
!!$  if(orbs%norbd > 0) then
!!$     call orthon_p(iproc,nproc,orbs%norbd,comms%nvctr_par(iproc,1),lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,&
!!$          psit(1+comms%nvctr_par(iproc,1)*orbs%norbu),orbs%nspinor)
!!$  end if
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

  if(orbs%nspinor==4) then
     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
          lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,psi,mom_vec)
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

  !here we should add the DIIS/SD switch
  !add switch between DIIS and SD
  !this section should be inserted into hpsitopsi
  if (energy == energy_min .and. .not. switchSD) idiistol=0
  if (energy > energy_min .and. idsx >0 .and. .not. switchSD) then
     idiistol=idiistol+1
  end if
  if (idiistol > idsx .and. .not. switchSD) then
     !the energy has not decreasing for too much steps, switching to SD for next steps
     if (iproc ==0) write(*,'(1x,a,1pe9.2,a)')&
          'WARNING: The energy value is growing (delta=',energy-energy_min,') switch to SD'
     switchSD=.true.
     i_all=-product(shape(psidst))*kind(psidst)
     deallocate(psidst,stat=i_stat)
     call memocc(i_stat,i_all,'psidst',subname)
     i_all=-product(shape(hpsidst))*kind(hpsidst)
     deallocate(hpsidst,stat=i_stat)
     call memocc(i_stat,i_all,'hpsidst',subname)
     i_all=-product(shape(ads))*kind(ads)
     deallocate(ads,stat=i_stat)
     call memocc(i_stat,i_all,'ads',subname)
     idsx_actual=0
     idiistol=0
  end if

  if ((energy == energy_min) .and. switchSD) then
     idiistol=idiistol+1
  end if
  if (idiistol > idsx .and. switchSD) then
     !restore the original DIIS
     if (iproc ==0) write(*,'(1x,a,1pe9.2)')&
          'WARNING: The energy value is now decreasing again, coming back to DIIS'
     switchSD=.false.
     idsx_actual=idsx
     ids=0
     idiistol=0

     allocate(psidst(comms%nvctr_par(iproc,1)*orbs%nspinor*orbs%norb*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,psidst,'psidst',subname)
     allocate(hpsidst(comms%nvctr_par(iproc,1)*orbs%nspinor*orbs%norb*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,hpsidst,'hpsidst',subname)
     allocate(ads(idsx+1,idsx+1,3+ndebug),stat=i_stat)
     call memocc(i_stat,ads,'ads',subname)
     call razero(3*(idsx+1)**2,ads)
  end if

END SUBROUTINE hpsitopsi
!!***

!!****f* BigDFT/first_orthon
!! FUNCTION
!!   First orthonormalisation
!! SOURCE
!!
subroutine first_orthon(iproc,nproc,orbs,wfd,comms,psi,hpsi,psit)
  use module_base
  use module_types
  use module_interfaces, except_this_one_B => first_orthon
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(communications_arrays), intent(in) :: comms
  real(wp), dimension(:) , pointer :: psi,hpsi,psit
  !local variables
  character(len=*), parameter :: subname='first_orthon'
  integer :: i_stat

!!$  if(nspin==4) then
!!$     nspinor=4
!!$  else
!!$     nspinor=1
!!$  end if
  
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

  call orthogonalize(iproc,nproc,orbs,comms,wfd,psit)

!!$  call orthon_p(iproc,nproc,orbs%norbu,comms%nvctr_par(iproc,1),wfd%nvctr_c+7*wfd%nvctr_f,&
!!$       psit,orbs%nspinor) 
!!$  if(orbs%norbd > 0) then
!!$     call orthon_p(iproc,nproc,orbs%norbd,comms%nvctr_par(iproc,1),wfd%nvctr_c+7*wfd%nvctr_f,&
!!$          psit(1+comms%nvctr_par(iproc,1)*orbs%norbu),orbs%nspinor) 
!!$  end if
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
!!***

!!****f* BigDFT/last_orthon
!! FUNCTION
!!   Transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
!! SOURCE
!!
subroutine last_orthon(iproc,nproc,orbs,wfd,nspin,comms,psi,hpsi,psit,evsum)
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
  !local variables
  character(len=*), parameter :: subname='last_orthon'
  logical :: dowrite !write the screen output
  integer :: i_all,i_stat,iorb,jorb,md
  real(wp) :: evpart
  real(wp), dimension(:,:,:), allocatable :: mom_vec

  call transpose_v(iproc,nproc,orbs,wfd,comms,&
       hpsi,work=psi)
  if (nproc==1) then
     psit => psi
     call transpose_v(iproc,nproc,orbs,wfd,comms,psit)
  end if

!!$  if(nspin==1.or.nspinor==4) then
!!$     call KStrans_p(iproc,nproc,norb,nvctrp,occup,hpsi,psit,evsum,eval,nspinor)
!!$  else
  call KStrans_p(iproc,nproc,orbs%norbu,comms%nvctr_par(iproc,1),orbs%occup,hpsi,psit,&
       evsum,orbs%eval,orbs%nspinor)
  evpart=evsum
  if(orbs%norbd > 0) then
     call KStrans_p(iproc,nproc,orbs%norbd,comms%nvctr_par(iproc,1),orbs%occup(orbs%norbu+1),&
          hpsi(1+comms%nvctr_par(iproc,1)*orbs%norbu),psit(1+comms%nvctr_par(iproc,1)*orbs%norbu),&
          evsum,orbs%eval(orbs%norbu+1),orbs%nspinor)
     evsum=evsum+evpart
  end if
!!$  end if

  call untranspose_v(iproc,nproc,orbs,wfd,comms,&
       psit,work=hpsi,outadd=psi(1))

  if (nproc > 1) then
     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  else
     nullify(psit)
  end if

  !for a non-collinear treatment,
  !we add the calculation of the moments for printing their value
  !close to the corresponding eigenvector
  if(orbs%nspinor==4) then
     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,wfd%nvctr_c+7*wfd%nvctr_f,&
          orbs%nspinor,psi,mom_vec)
  end if

  !print the found eigenvalues
  if (iproc == 0) then
     write(*,'(1x,a)')&
          '-------------------------------------------------------------- Kohn-Sham Eigenvalues'
     if (orbs%nspinor ==4) then
        write(*,'(1x,a)')&
             '           Eigenvalue                                      m_x       m_y       m_z'
     end if
     if (nspin==1.or.orbs%nspinor==4) then
        do iorb=1,orbs%norb
           dowrite =(iorb <= 5 .or. iorb >= orbs%norb-5) .or. verbose > 0
           if (orbs%nspinor ==4) then
              if (dowrite) & 
                   write(*,'(1x,a,i4,a,1x,1pe21.14,20x,(1x,3(0pf10.5)))') &
                   'eval(',iorb,')=',orbs%eval(iorb),(mom_vec(md,iorb,1)/mom_vec(1,iorb,1),md=2,4)
           else
              if (dowrite) & 
                   write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,')=',orbs%eval(iorb)
           end if
        end do
     else
        do iorb=1,min(orbs%norbu,orbs%norbd)
           jorb=orbs%norbu+iorb
           dowrite =(iorb <= 5 .or. iorb >= min(orbs%norbu,orbs%norbd)-5)  .or. verbose > 0
           if (dowrite) & 
                write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
                'eval(',iorb,',u)=',orbs%eval(iorb),'eval(',iorb,',d)=',orbs%eval(jorb)
        end do
        if (orbs%norbu > orbs%norbd) then
           do iorb=orbs%norbd+1,orbs%norbu
              dowrite =(iorb <= 5 .or. iorb >= orbs%norbu-5) .or. verbose > 0
              if (dowrite) & 
              write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,',u)=',orbs%eval(iorb)
           end do
        else if (orbs%norbd > orbs%norbu) then
           do iorb=2*orbs%norbu+1,orbs%norbu+orbs%norbd
              dowrite =(iorb <= 5 .or. iorb >= orbs%norbd-5) .or. verbose > 0
              if (dowrite) & 
                   write(*,'(50x,a,i4,a,1x,1pe21.14)') 'eval(',iorb-orbs%norbu,',d)=',orbs%eval(iorb)
           end do
        end if
     end if
  end if

  if (orbs%nspinor ==4) then
     i_all=-product(shape(mom_vec))*kind(mom_vec)
     deallocate(mom_vec,stat=i_stat)
     call memocc(i_stat,i_all,'mom_vec',subname)
  end if

  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)

end subroutine last_orthon
!!***


!!****f* BigDFT/calc_moments
!! FUNCTION
!!   Calculate magnetic moments
!! SOURCE
!!
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
!!***

!experimental routine for correcting the potential from a vacancy
subroutine correct_hartree_potential(at,iproc,nproc,n1i,n2i,n3i,n3p,n3pi,n3d,&
     i3s,i3xcsh,hxh,hyh,hzh,pkernel,ngatherarr,&
     rhoref,pkernel_ref,pot_ion,rhopot,ixc,nspin,ehart,eexcu,vexcu,PSquiet,correct_offset)
  use module_base
  use module_types
  use Poisson_Solver
  implicit none
  character(len=3), intent(in) :: PSquiet
  logical, intent(in) :: correct_offset
  integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,n3p,n3pi,n3d,nspin,ixc,i3xcsh,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(dp), dimension(n1i,n2i,max(n3d,1),nspin), intent(inout) :: rhoref
  real(dp), dimension(n1i,n2i,max(n3pi,1)), intent(inout) :: pot_ion
  real(dp), dimension(n1i,n2i,max(n3d,1),nspin), intent(inout) :: rhopot
  real(gp), intent(out) :: ehart,eexcu,vexcu
  real(dp), dimension(:), pointer :: pkernel_ref,pkernel
  !local variables
  character(len=*), parameter :: subname='correct_hartree_potential'
  integer :: i_all,i_stat,ierr,i1,i2,i3,ispin
  real(gp) :: ehart_fake,eexcu_fake,vexcu_fake,ehartA,ehartB,tt,offset
  real(dp), dimension(:,:,:,:), allocatable :: potref,drho,vxc

  allocate(potref(n1i,n2i,max(n3d,1),nspin+ndebug),stat=i_stat)
  call memocc(i_stat,potref,'potref',subname)

  allocate(drho(n1i,n2i,n3i,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,drho,'drho',subname)

  allocate(vxc(n1i,n2i,max(n3p,1),nspin+ndebug),stat=i_stat)
  call memocc(i_stat,vxc,'vxc',subname)

  !Delta rho = rho - rho_ref
  do ispin=1,nspin
     do i3=1,n3d
        do i2=1,n2i
           do i1=1,n1i
              drho(i1,i2,i3s-1+i3,ispin)=&
                   rhopot(i1,i2,i3,ispin)-rhoref(i1,i2,i3,ispin)
           end do
        end do
     end do
  end do

!!$  call plot_density(at%geocode,'deltarho.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
!!$       1,at%alat1,at%alat2,at%alat3,ngatherarr,drho(1,1,i3s+i3xcsh,1))

  !calculate the offset
  tt=0.d0
  do i3=1,n3p
     do i2=1,n2i
        do i1=1,n1i
           tt=tt+drho(i1,i2,i3+i3s-1+i3xcsh,1)
        enddo
     enddo
  enddo
  !tt=tt*hxh*hyh*hzh
  if (nproc > 1) then
     call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     offset=tt
  end if

  if (iproc==0) print *,'charge deltarho',offset*hxh*hyh*hzh

  !rho_tot -> VH_tot & VXC_tot 
  call PSolver(at%geocode,'D',iproc,nproc,n1i,n2i,n3i,&
       ixc,hxh,hyh,hzh,&
       rhopot,pkernel,vxc,ehart,eexcu,vexcu,0.d0,.false.,nspin,&
       quiet=PSquiet)

  !calculate the reference Hartree potential
  !here the offset should be specified for charged systems
  call dcopy(n1i*n2i*n3d,rhoref,1,potref,1) 
  !Rho_ref -> VH_ref
  call PSolver(at%geocode,'D',iproc,nproc,n1i,n2i,n3i,&
       0,hxh,hyh,hzh,&
       potref,pkernel,pot_ion,ehart_fake,eexcu_fake,vexcu_fake,0.d0,.false.,1,&
       quiet=PSquiet)
  
  !save the total density in rhopot
  do ispin=1,nspin
     do i3=1,n3p
        do i2=1,n2i
           do i1=1,n1i
              rhopot(i1,i2,i3,ispin)=rhoref(i1,i2,i3,ispin)+&
                   drho(i1,i2,i3s+i3xcsh-1+i3,ispin)
           end do
        end do
     end do
  end do
  !rhopot=rhoref+drho

  !gather the total density difference in drho
  if (nproc > 1) then
     do ispin=1,nspin
        call MPI_ALLGATHERV(MPI_IN_PLACE,n1i*n2i*n3p,&
             mpidtypd,drho(1,1,1,ispin),ngatherarr(0,1),&
             ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
     end do
  end if

  !calculate the vacancy hartree potential
  !Delta rho -> VH_drho(I)
  !use global distribution scheme for deltarho
  call PSolver('F','G',iproc,nproc,n1i,n2i,n3i,&
       0,hxh,hyh,hzh,&
       drho,pkernel_ref,pot_ion,ehart_fake,eexcu_fake,vexcu_fake,0.d0,.false.,1,&
       quiet=PSquiet)

!!$  call plot_density(at%geocode,'VHdeltarho.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
!!$       1,at%alat1,at%alat2,at%alat3,ngatherarr,drho(1,1,1+i3xcsh,1))


  !sum the complete hartree potential
  call axpy(n1i*n2i*n3p,1.0_dp,potref(1,1,1+i3xcsh,1),1,drho(1,1,i3s+i3xcsh,1),1)

  if (correct_offset) then
     !calculate the offset
     tt=0.d0
     do i3=1,n3p
        do i2=1,n2i
           do i1=1,n1i
              tt=tt+drho(i1,i2,i3+i3s-1+i3xcsh,1)
           end do
        end do
     end do
     !tt=tt*hxh*hyh*hzh
     if (nproc > 1) then
        call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
             MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        offset=tt
     end if

     if (iproc==0) print *,'offset to subtract',offset

     !now the hartree potential has zero integral
     drho=drho-offset/real(n1i*n2i*n3i,dp)

     !calculate the offset
     tt=0.d0
     do i3=1,n3p
        do i2=1,n2i
           do i1=1,n1i
              tt=tt+drho(i1,i2,i3+i3s-1+i3xcsh,1)
           enddo
        enddo
     enddo
     if (nproc > 1) then
        call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
             MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        offset=tt
     end if

     if (iproc==0) print *,'offset (verify)',offset*hxh*hyh*hzh
  end if

  !calculate total hartree energy
  ehartA=0.0_dp
  ehartA=0.5_dp*hxh*hyh*hzh*&
       dot(n1i*n2i*n3p,rhopot(1,1,1+i3xcsh,1),1,drho(1,1,i3s+i3xcsh,1),1)

  !reduce the result
  if (nproc > 1) then
     call MPI_ALLREDUCE(ehartA,ehartB,1,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     ehartB=ehartA
  end if
  ehart=ehartB

  !sum the different potentials (valid only for non-spin polarised systems)
  call axpy(n1i*n2i*n3p,1.0_dp,pot_ion(1,1,1),1,drho(1,1,i3s+i3xcsh,1),1)

!!$  call plot_density(at%geocode,'VHpVion.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
!!$       1,at%alat1,at%alat2,at%alat3,ngatherarr,drho(1,1,i3s+i3xcsh,1))
!!$
!!$  !calculate the offset
!!$  tt=0.d0
!!$  do i3=1,n3p
!!$     do i2=1,n2i
!!$        do i1=1,n1i
!!$           tt=tt+drho(i1,i2,i3+i3xcsh,1)
!!$        enddo
!!$     enddo
!!$  enddo
!!$  !tt=tt*hxh*hyh*hzh
!!$  if (nproc > 1) then
!!$     call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
!!$          MPI_SUM,MPI_COMM_WORLD,ierr)
!!$  else
!!$     offset=tt
!!$  end if
!!$
!!$  if (iproc==0) print *,'offset of Vh+Vion',offset

  call axpy(n1i*n2i*n3p*nspin,1.0_dp,vxc(1,1,1,1),1,drho(1,1,i3s+i3xcsh,1),1)

  !final result in rhopot
  call dcopy(n1i*n2i*n3d*nspin,drho(1,1,i3s,1),1,rhopot,1) 


!!$  call plot_density(at%geocode,'Vtot.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
!!$       1,at%alat1,at%alat2,at%alat3,ngatherarr,rhopot(1,1,1+i3xcsh,1))

  !calculate the offset
  tt=0.d0
  do i3=1,n3p
     do i2=1,n2i
        do i1=1,n1i
           tt=tt+rhopot(i1,i2,i3+i3xcsh,1)
        enddo
     enddo
  enddo
  !tt=tt*hxh*hyh*hzh
  if (nproc > 1) then
     call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     offset=tt
  end if

  if (iproc==0) print *,'offset of the total potential',offset

  i_all=-product(shape(potref))*kind(potref)
  deallocate(potref,stat=i_stat)
  call memocc(i_stat,i_all,'potref',subname)

  i_all=-product(shape(vxc))*kind(vxc)
  deallocate(vxc,stat=i_stat)
  call memocc(i_stat,i_all,'vxc',subname)

  i_all=-product(shape(drho))*kind(drho)
  deallocate(drho,stat=i_stat)
  call memocc(i_stat,i_all,'drho',subname)

end subroutine correct_hartree_potential


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
  integer :: i,ispinor,iorb,indspin,indorb,jproc,i_stat,i_all,iscomp,idsx,index
  real(wp) :: vali,valorb,psival,maxdiff,ierr
  real(wp), dimension(:), allocatable :: psi
  real(wp), dimension(:), pointer :: pwork
  real(wp) :: epsilon

  !allocate the "wavefunction" amd fill it, and also the workspace
  allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(pwork(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,pwork,'pwork',subname)

  do iorb=1,orbs%norbp
     valorb=real(orbs%isorb+iorb,wp)
     indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
     do ispinor=1,orbs%nspinor
        indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
        do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           vali=real(i,wp)*1.d-5
           psi(i+indspin+indorb)=(valorb+vali)*(-1)**(ispinor-1)
        end do
     end do
  end do

  !transpose the hpsi wavefunction
  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psi,work=pwork)

  !check the results of the transposed wavefunction
  maxdiff=0.0_wp
  !calculate the starting point for the component distribution
  iscomp=0
  if (orbs%nkptsp /=0) then
     do jproc=0,iproc-1
        iscomp=iscomp+comms%nvctr_par(jproc,1)
     end do
     do iorb=1,orbs%norb
        valorb=real(iorb,wp)
        indorb=(iorb-1)*(comms%nvctr_par(iproc,1))*orbs%nspinor
        do idsx=1,(orbs%nspinor-1)/2+1
           do i=1,comms%nvctr_par(iproc,1)
              vali=real(i+iscomp,wp)*1.d-5
              do ispinor=1,((2+orbs%nspinor)/4+1)
                 psival=(-1)**(ispinor-1)*(valorb+vali)
                 if (psival .lt. 0.d0) then  !this is just to force the IEEE representation of psival
                    write(321,*) psival,psival**2
                 endif
                 index=ispinor+(i-1)*((2+orbs%nspinor)/4+1)+&
                      (idsx-1)*((2+orbs%nspinor)/4+1)*comms%nvctr_par(iproc,1)+indorb
                 maxdiff=max(abs(psi(index)-psival),maxdiff)
              end do
           end do
        end do
     end do
  end if
  if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
     write(*,*)'ERROR: process',iproc,'does not transpose wavefunctions correctly!'
     write(*,*)'       found an error of',maxdiff,'cannot continue.'
     write(*,*)'       data are written in the file transerror.log, exiting...'

     open(unit=22,file='transerror.log',status='unknown')
     do iorb=1,orbs%norb
        valorb=real(iorb,wp)
        indorb=(iorb-1)*(comms%nvctr_par(iproc,1))*orbs%nspinor
        do idsx=1,(orbs%nspinor-1)/2+1
           do i=1,comms%nvctr_par(iproc,1)
              vali=real(i+iscomp,wp)*1.d-5
              do ispinor=1,((2+orbs%nspinor)/4+1)
                 psival=(-1)**(ispinor-1)*(valorb+vali)
                 index=ispinor+(i-1)*((2+orbs%nspinor)/4+1)+&
                      (idsx-1)*((2+orbs%nspinor)/4+1)*comms%nvctr_par(iproc,1)+indorb
                 maxdiff=abs(psi(index)-psival)
                 write(22,'(i3,i6,i5,3(1x,1pe13.6))')ispinor,i+iscomp,iorb,psival,&
                      psi(index),maxdiff
              end do
           end do
        end do
     end do
     close(unit=22)

     call MPI_ABORT(MPI_COMM_WORLD,ierr)

  end if

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,&
       psi,work=pwork)

  maxdiff=0.0_wp
  do iorb=1,orbs%norbp
     valorb=real(orbs%isorb+iorb,wp)
     indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
     do ispinor=1,orbs%nspinor
        indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
        do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           vali=real(i,wp)*1.d-5
           psival=(valorb+vali)*(-1)**(ispinor-1)
           maxdiff=max(abs(psi(i+indspin+indorb)-psival),maxdiff)
        end do
     end do
  end do

  if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
     write(*,*)'ERROR: process',iproc,'does not untranspose wavefunctions correctly!'
     write(*,*)'       found an error of',maxdiff,'cannot continue.'
     write(*,*)'       data are written in the file transerror.log, exiting...'

     open(unit=22,file='transerror.log',status='unknown')
     maxdiff=0.0_wp
     do iorb=1,orbs%norbp
        valorb=real(orbs%isorb+iorb,wp)
        indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
        do ispinor=1,orbs%nspinor
           indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
           do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
              vali=real(i,wp)*1.d-5
              psival=(valorb+vali)*(-1)**(ispinor-1)
              maxdiff=abs(psi(i+indspin+indorb)-psival)
              write(22,'(i3,i6,i5,3(1x,1pe13.6))')ispinor,i,iorb+orbs%isorb,psival,&
                   psi(ispinor+(i-1)*orbs%nspinor+indorb),maxdiff
           end do
        end do
     end do
     close(unit=22)

     call MPI_ABORT(MPI_COMM_WORLD,ierr)

  end if

  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(pwork))*kind(pwork)
  deallocate(pwork,stat=i_stat)
  call memocc(i_stat,i_all,'pwork',subname)



end subroutine check_communications
