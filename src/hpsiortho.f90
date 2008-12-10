subroutine HamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
     cpmult,fpmult,radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,psi,hpsi,&
     ekin_sum,epot_sum,eproj_sum,nspin)
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
  !local variables
  character(len=*), parameter :: subname='HamiltonianApplication'
  integer :: i_all,i_stat,ierr,n1i,n2i,n3i,iorb,ispin
  real(gp) :: eproj
  real(gp), dimension(3,2) :: wrkallred
  real(wp), dimension(:,:), pointer :: pot
  integer,parameter::lupfil=14

  call timing(iproc,'ApplyLocPotKin','ON')

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0) then
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
  call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum)
  
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
     call applyprojectorsonthefly(iproc,orbs%nspinor,orbs%norbp,&
          orbs%occup(min(orbs%iorbs+1,orbs%norb),at,lr%d%n1,lr%d%n2,lr%d%n3,&
          rxyz,hx,hy,hz,cpmult,fpmult,radii_cf,lr%wfd,nlpspd,proj,psi,hpsi,eproj_sum)
  else
     !one should add a flag here which states that it works only for global reion
     ! loop over all my orbitals
     do iorb=1,orbs%norbp*orbs%nspinor
        call applyprojectorsone(at%ntypes,at%nat,at%iatype,at%psppar,at%npspcode, &
             nlpspd%nprojel,nlpspd%nproj,nlpspd%nseg_p,nlpspd%keyg_p,nlpspd%keyv_p,&
             nlpspd%nvctr_p,&
             proj,lr%wfd%nseg_c,lr%wfd%nseg_f,lr%wfd%keyg,lr%wfd%keyv,&
             lr%wfd%nvctr_c,lr%wfd%nvctr_f, & 
             psi(1,iorb),hpsi(1,iorb),eproj)
        eproj_sum=eproj_sum+occup((iorb+orbs%iorbs-1)/orbs%nspinor+1)*eproj
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

end subroutine HamiltonianApplication


subroutine hpsitopsi(iproc,nproc,orbs,hx,hy,hz,nvctrp,lr,comms,&
     ncong,iter,idsx,idsx_actual,ads,energy,energy_old,energy_min,&
     alpha,gnrm,scprsum,psi,psit,hpsi,psidst,hpsidst,nspin)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,idsx,iter,nvctrp,nspin
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
  !local variables
  character(len=*), parameter :: subname='hpsitopsi'
  logical, save :: switchSD
  integer, save :: idiistol,mids,ids
  integer :: ierr,ind,i1,i2,iorb,i,k,i_stat,i_all,oidx,sidx
  real(wp) :: cprecr
  real(dp) :: tt,scpr,scprpart
  real(wp), dimension(:,:,:), allocatable :: mom_vec
  real(kind=4), dimension(:), allocatable :: psitcuda,hpsitcuda

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

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, orthoconstraint...'
  end if

  !transpose the hpsi wavefunction
  call transpose_v(iproc,nproc,orbs%norbp,orbs%nspinor,lr%wfd,nvctrp,comms,hpsi,work=psi)

  if (nproc == 1) then
     !associate psit pointer for orthoconstraint and transpose it (for the non-collinear case)
     psit => psi
     call transpose_v(iproc,nproc,orbs%norbp,orbs%nspinor,lr%wfd,nvctrp,comms,psit)
  end if

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  ! insert branching for CUDA section(experimental)
  ! once the mixed precision version is ready such part can be eliminated
  if (GPUblas) then
     allocate(psitcuda(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psitcuda,'psitcuda',subname)
     allocate(hpsitcuda(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hpsitcuda,'hpsitcuda',subname)

     do i=1,orbs%npsidim
        psitcuda(i)=real(psit(i),kind=4)
        hpsitcuda(i)=real(hpsi(i),kind=4)
     end do

     if(nspin==1.or.nspinor==4) then
        call orthoconstraint_cuda(iproc,nproc,orbs%norb,orbs%occup,nvctrp,&
             psitcuda,hpsitcuda,scprsum,orbs%nspinor)
     else
        call orthoconstraint_cuda(iproc,nproc,orbs%norbu,orbs%occup,nvctrp,&
             psitcuda,hpsitcuda,scprsum,orbs%nspinor)
        scprpart=0.0d0
        if(orbs%norbd > 0) then
           scprpart=scprsum 
           call orthoconstraint_cuda(iproc,nproc,orbs%norbd,orbs%occup(orbs%norbu+1),&
                nvctrp,psitcuda(1+nvctrp*orbs%norbu),hpsitcuda(1+nvctrp*orbs%norbu),&
                scprsum,orbs%nspinor)
        end if
        scprsum=scprsum+scprpart
     end if

     do i=1,orbs%npsidim
        psit(i)=real(psitcuda(i),wp)
        hpsi(i)=real(hpsitcuda(i),wp)
     end do

     i_all=-product(shape(psitcuda))*kind(psitcuda)
     deallocate(psitcuda,stat=i_stat)
     call memocc(i_stat,i_all,'psitcuda',subname)
     i_all=-product(shape(hpsitcuda))*kind(hpsitcuda)
     deallocate(hpsitcuda,stat=i_stat)
     call memocc(i_stat,i_all,'hpsitcuda',subname)
  else
!!$     if(nspin==1 .or. orbs%nspinor==4) then
!!$        call orthoconstraint_p(iproc,nproc,orbs%norb,orbs%occup,nvctrp,psit,hpsi,&
!!$             scprsum,orbs%nspinor)
!!$     else
     call orthoconstraint_p(iproc,nproc,orbs%norbu,orbs%occup,nvctrp,psit,hpsi,&
          scprsum,orbs%nspinor)
     scprpart=0.0_dp
     if(orbs%norbd > 0) then
        scprpart=scprsum 
        call orthoconstraint_p(iproc,nproc,orbs%norbd,orbs%occup(orbs%norbu+1),nvctrp,&
             psit(1+nvctrp*orbs%norbu),hpsi(1+nvctrp*orbs%norbu),scprsum,orbs%nspinor)
!!$        end if
        scprsum=scprsum+scprpart
     end if
  end if

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,orbs%norbp,orbs%nspinor,lr%wfd,nvctrp,comms,&
       hpsi,work=psi)

  call timing(iproc,'Precondition  ','ON')
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, preconditioning...'
  end if


  !Preconditions all orbitals belonging to iproc
  !and calculate the partial norm of the residue
  call preconditionall(iproc,nproc,orbs%norbp,lr,hx,hy,hz,ncong,orbs%nspinor,&
       orbs%eval(min(orbs%isorb+1,orbs%norb)),hpsi,gnrm,lr%hybrid_on)

  !sum over all the partial residues
  if (nproc > 1) then
     tt=gnrm
     call MPI_ALLREDUCE(tt,gnrm,1,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif
  gnrm=sqrt(gnrm/real(orbs%norb,dp))

  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if
  call timing(iproc,'Precondition  ','OF')

  !apply the minimization method (DIIS or steepest descent)
  if (idsx_actual > 0) then
     !transpose the hpsi wavefunction into the diis array
     call transpose_v(iproc,nproc,orbs%norbp,orbs%nspinor,lr%wfd,nvctrp,comms,&
          hpsi,work=psi,outadd=hpsidst(1+nvctrp*orbs%nspinor*orbs%norbp*nproc*(mids-1)))

     call timing(iproc,'Diis          ','ON')
     if (nproc > 1) then
        do i=1,nvctrp*orbs%norb*orbs%nspinor
           psidst(i+nvctrp*orbs%nspinor*orbs%norbp*nproc*(mids-1))= psit(i)
        enddo
     else
        do i=1,nvctrp*orbs%norb*orbs%nspinor
           psidst(i+nvctrp*orbs%nspinor*orbs%norbp*nproc*(mids-1))= psit(i)
           hpsidst(i+nvctrp*orbs%nspinor*orbs%norbp*nproc*(mids-1))=hpsi(i)
        enddo
     endif

     call diisstp(orbs%norb,orbs%norbp,nproc,iproc,orbs%nspinor,  &
          ads,ids,mids,idsx_actual,nvctrp,psit,psidst,hpsidst)
  else
     ! update all wavefunctions with the preconditioned gradient
     if (energy > energy_old) then
        alpha=max(.125_wp,.5_wp*alpha)
        if (alpha == .125_wp) write(*,*) 'Convergence problem or limit'
     else
        alpha=min(1.05_wp*alpha,1._wp)
     endif
     if (iproc == 0) write(*,'(1x,a,1pe11.3)') 'alpha=',alpha

     !transpose the hpsi wavefunction
     call transpose_v(iproc,nproc,orbs%norbp,orbs%nspinor,lr%wfd,nvctrp,comms,&
          hpsi,work=psi)

     call timing(iproc,'Diis          ','ON')
     do iorb=1,orbs%norb*orbs%nspinor
        call axpy(nvctrp,-alpha,hpsi(1+nvctrp*(iorb-1)),1,psit(1+nvctrp*(iorb-1)),1)
     enddo
  endif

  call timing(iproc,'Diis          ','OF')

  if (iproc == 0) then
     write(*,'(1x,a)',advance='no')&
          'Orthogonalization...'
  end if

!!$  if(nspin==1 .or. nspinor==4) then
!!$     call orthon_p(iproc,nproc,norb,nvctrp,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,psit,nspinor)
!!$  else
  call orthon_p(iproc,nproc,orbs%norbu,nvctrp,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,&
       psit,orbs%nspinor)
  if(norbd > 0) then
     call orthon_p(iproc,nproc,orbs%norbd,nvctrp,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,&
          psit(1+nvctrp*norbu),orbs%nspinor)
  end if
!!$  end if
    !       call checkortho_p(iproc,nproc,norb,nvctrp,psit)
  
  call untranspose_v(iproc,nproc,orbs%norbp,orbs%nspinor,lr%wfd,nvctrp,comms,&
       psit,work=hpsi,outadd=psi(1))
  if (nproc == 1) then
     nullify(psit)
  end if
  
  if (iproc == 0) then
     write(*,'(1x,a)')&
          'done.'
  end if

  if(orbs%nspinor==4) then
     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
          lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,psi,mom_vec)
     !only the root process has the correct array
     if(iproc==0) then
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

     allocate(psidst(nvctrp*orbs%nspinor*orbs%norbp*nproc*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,psidst,'psidst',subname)
     allocate(hpsidst(nvctrp*orbs%nspinor*orbs%norbp*nproc*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,hpsidst,'hpsidst',subname)
     allocate(ads(idsx+1,idsx+1,3+ndebug),stat=i_stat)
     call memocc(i_stat,ads,'ads',subname)
     call razero(3*(idsx+1)**2,ads)
  end if


end subroutine hpsitopsi

!calculate the address to start from for calculating the 
!norm of the residue if hpsi is allocated in the transposed way
!it can be eliminated when including all this procedure in a subroutine
!in other terms, it takes the i1,i2 component of an array psi(nvctr,norbp) 
!from an array of the form psi(nvctrp,norb)
!for this routine norbp is not needed
subroutine trans_address(nvctrp,nvctr,i,iorb,i1,i2)
  implicit none
  integer, intent(in) :: nvctrp,nvctr,i,iorb
  integer, intent(out) :: i1,i2
  !local variables
  integer :: ind
!  if (nproc > 1) then
     ind=i+nvctr*(iorb-1)
     i1=mod(ind-1,nvctrp)+1
     i2=(ind-i1)/nvctrp+1
!!$  else
!!$     i1=1
!!$     i2=iorb
!!$  end if
end subroutine trans_address

subroutine first_orthon(iproc,nproc,orbs,wfd,nvctrp,comms,psi,hpsi,psit)
  use module_base
  use module_types
  use module_interfaces, except_this_one => first_orthon
  implicit none
  integer, intent(in) :: iproc,nproc,nvctrp
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(communications_arrays), intent(in) :: comms
  real(wp), dimension(:) , pointer :: psi,hpsi,psit
  !local variables
  character(len=*), parameter :: subname='first_orthon'
  integer :: i_all,i_stat,ierr,iorb

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
  call transpose_v(iproc,nproc,orbs%norbp,orbs%nspinor,wfd,nvctrp,comms,psi,&
       work=hpsi,outadd=psit(1))

!!$  if(nspin==1 .or. nspinor==4) then
!!$     call orthon_p(iproc,nproc,norb,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,nspinor) 
!!$  else
  call orthon_p(iproc,nproc,orbs%norbu,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,orbs%nspinor) 
  if(norbd > 0) then
     call orthon_p(iproc,nproc,orbs%norbd,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,&
          psit(1+nvctrp*orbs%norbu),orbs%nspinor) 
  end if
!!$  end if
  !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

  call untranspose_v(iproc,nproc,orbs%norbp,orbs%nspinor,wfd,nvctrp,comms,psit,&
       work=hpsi,outadd=psi(1))

  if (nproc == 1) then
     nullify(psit)
     !allocate hpsi array
     allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hpsi,'hpsi',subname)
  end if

end subroutine first_orthon

! transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,orbs,wfd,nvctrp,&
     nspin,comms,psi,hpsi,psit,evsum)
  use module_base
  use module_types
  use module_interfaces, except_this_one => last_orthon
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  integer, intent(in) :: iproc,nproc,nvctrp,nspin
  real(wp), intent(out) :: evsum
  real(wp), dimension(:) , pointer :: psi,hpsi,psit
  !local variables
  character(len=*), parameter :: subname='last_orthon'
  integer :: i_all,i_stat,ierr,iorb,jorb,md
  real(wp) :: evpart
  real(wp), dimension(:,:,:), allocatable :: mom_vec

  call transpose_v(iproc,nproc,orbs%norbp,orbs%nspinor,wfd,nvctrp,comms,&
       hpsi,work=psi)
  if (nproc==1) then
     psit => psi
     call transpose_v(iproc,nproc,orbs%norbp,orbs%nspinor,wfd,nvctrp,comms,psit)
  end if

!!$  if(nspin==1.or.nspinor==4) then
!!$     call KStrans_p(iproc,nproc,norb,nvctrp,occup,hpsi,psit,evsum,eval,nspinor)
!!$  else
  call KStrans_p(iproc,nproc,orbs%norbu,nvctrp,orbs%occup,hpsi,psit,&
       evsum,orbs%eval,orbs%nspinor)
  evpart=evsum
  if(orbs%norbd > 0) then
     call KStrans_p(iproc,nproc,orbs%norbd,nvctrp,orbs%occup(orbs%norbu+1),&
          hpsi(1+nvctrp*orbs%norbu),psit(1+nvctrp*orbs%norbu),&
          evsum,orbs%eval(orbs%norbu+1),orbs%nspinor)
     evsum=evsum+evpart
  end if
!!$  end if

  call untranspose_v(iproc,nproc,orbs%norbp,orbs%nspinor,wfd,nvctrp,comms,&
       psit,work=hpsi,outadd=psi(1))

  if (nproc > 1) then
     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  else
     nullify(psit)
  end if

  !for a non-collinear treatment,
  !here we can add the calculation of the moments for printing their value
  !close to the corresponding eigenvector
  if(orbs%nspinor==4) then
     allocate(mom_vec(4,norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,wfd%nvctr_c+7*wfd%nvctr_f,&
          orbs%nspinor,psi,mom_vec)
  end if

  !print the found eigenvalues
  if (iproc == 0) then
     write(*,'(1x,a)')&
          '-------------------------------------------------------------- Kohn-Sham Eigenvalues'
     if (nspin==1.or.orbs%nspinor==4) then
        if (orbs%nspinor ==4) then
        write(*,'(1x,a)')&
             '           Eigenvalue                                      m_x       m_y       m_z'
           do iorb=1,orbs%norb
              write(*,'(1x,a,i4,a,1x,1pe21.14,20x,(1x,3(0pf10.5)))') &
                   'eval(',iorb,')=',orbs%eval(iorb),(mom_vec(md,iorb,1)/mom_vec(1,iorb,1),md=2,4)
           end do
        else
           do iorb=1,orbs%norb
              write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,')=',orbs%eval(iorb)
           end do
        end if
     else
        do iorb=1,min(orbs%norbu,orbs%norbd)
           jorb=orbs%norbu+iorb
           write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
                'eval(',iorb,',u)=',orbs%eval(iorb),'eval(',iorb,',d)=',orbs%eval(jorb)
        end do
        if (orbs%norbu > orbs%norbd) then
           do iorb=orbs%norbd+1,orbs%norbu
              write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,',u)=',orbs%eval(iorb)
           end do
        else if (orbs%norbd > orbs%norbu) then
           do iorb=2*orbs%norbu+1,orbs%norbu+orbs%norbd
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


subroutine calc_moments(iproc,nproc,norb,norb_par,nvctr,nspinor,psi,mom_vec)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctr,nspino
  real(wp), dimension(0:nproc-1), intent(in) :: norb_par
  real(wp), dimension(nvctr,norb*nspinor), intent(in) :: psi
  real(wp), dimension(4,norb,min(nproc,2)), intent(out) :: mom_vec
  !local variables
  character(len=*), parameter :: subname='calc_moments'
  integer :: i_all,i_stat,ierr,iorb,jproc
  integer :: ispin,md,ndim,oidx
  integer, dimension(:), allocatable :: norb_displ
  real(wp) :: m00,m11,m13,m24,m12,m34,m14,m23,dot
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
           norb_displ(jproc)=norb_displ(jproc-1)+norb_par(jproc)
        end do
        
        call MPI_GATHERV(mom_vec(1,1,2),4*norbp,mpidtypw,&
             mom_vec(1,1,1),4*norb_par,4*norb_displ,mpidtypw,&
             0,MPI_COMM_WORLD,ierr)

        i_all=-product(shape(norb_displ))*kind(norb_displ)
        deallocate(norb_displ,stat=i_stat)
        call memocc(i_stat,i_all,'norb_displ',subname)
     end if
     
     
    
  end if

end subroutine calc_moments
