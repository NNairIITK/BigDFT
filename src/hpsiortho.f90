subroutine HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
     psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     wfd,bounds,nlpspd,proj,ngatherarr,n3p,&
     potential,psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,spinar)

  use module_types

  implicit none
  include 'mpif.h'
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(convolutions_bounds), intent(in) :: bounds
  logical, intent(in) :: parallel
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nat,ntypes,n3p
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspin
  real(kind=8), intent(in) :: hgrid
  integer, dimension(ntypes), intent(in) :: npspcode
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(kind=8), dimension(norb), intent(in) :: occup,spinar
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
  real(kind=8), dimension(*), intent(in) :: potential
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(out) :: hpsi
  real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum


  !local variables
  integer :: i_all,i_stat,ierr,iorb
  real(kind=8) :: eproj
  real(kind=8), dimension(:), allocatable :: pot
  real(kind=8), dimension(:,:), allocatable :: wrkallred

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  call timing(iproc,'ApplyLocPotKin','ON')

  if (datacode=='D') then
     !allocate full potential
     allocate(pot((2*n1+31)*(2*n2+31)*(2*n3+31)*nspin),stat=i_stat)
     call memocc(i_stat,product(shape(pot))*kind(pot),'pot','hamiltonianapplication')

     call MPI_ALLGATHERV(potential,(2*n1+31)*(2*n2+31)*n3p,&
          MPI_DOUBLE_PRECISION,pot,ngatherarr(0,1),&
          ngatherarr(0,2),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     if(nspin==2) then
        call MPI_ALLGATHERV(potential((2*n1+31)*(2*n2+31)*n3p+1),(2*n1+31)*(2*n2+31)*n3p,&
             MPI_DOUBLE_PRECISION,pot((2*n1+31)*(2*n2+31)*(2*n3+31)+1),ngatherarr(0,1),&
             ngatherarr(0,2),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     end if
     call applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,0, &
          hgrid,occup,wfd%nseg_c,wfd%nseg_f,wfd%nvctr_c,wfd%nvctr_f,wfd%keyg,wfd%keyv,&
          bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c,&
          bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f, &
          psi,pot,hpsi,epot_sum,ekin_sum,nspin,spinar,&
          bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
          bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
          bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
          bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,bounds%ibyyzz_r)

     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot','hamiltonianapplication')

  else

     call applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,0, &
          hgrid,occup,wfd%nseg_c,wfd%nseg_f,wfd%nvctr_c,wfd%nvctr_f,wfd%keyg,wfd%keyv,&
          bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c,&
          bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f, &
          psi,potential,hpsi,epot_sum,ekin_sum,nspin,spinar,&
          bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
          bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
          bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
          bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,bounds%ibyyzz_r)
  end if

  call timing(iproc,'ApplyLocPotKin','OF')

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')

  eproj_sum=0.d0
  ! loop over all my orbitals
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     call applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
          nlpspd%nprojel,nlpspd%nproj,nlpspd%nseg_p,nlpspd%keyg_p,nlpspd%keyv_p,nlpspd%nvctr_p,&
          proj,  &
          wfd%nseg_c,wfd%nseg_f,wfd%keyg,wfd%keyv,wfd%nvctr_c,wfd%nvctr_f,  & 
          psi(1,iorb-iproc*norbp),hpsi(1,iorb-iproc*norbp),eproj)
     eproj_sum=eproj_sum+occup(iorb)*eproj
     !     write(*,*) 'iorb,eproj',iorb,eproj
  enddo

  call timing(iproc,'ApplyProj     ','OF')

  if (parallel) then
     allocate(wrkallred(3,2),stat=i_stat)
     call memocc(i_stat,product(shape(wrkallred))*kind(wrkallred),'wrkallred','hamiltonianapplication')

     wrkallred(1,2)=ekin_sum 
     wrkallred(2,2)=epot_sum 
     wrkallred(3,2)=eproj_sum
     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     ekin_sum=wrkallred(1,1)
     epot_sum=wrkallred(2,1)
     eproj_sum=wrkallred(3,1) 

     i_all=-product(shape(wrkallred))*kind(wrkallred)
     deallocate(wrkallred,stat=i_stat)
     call memocc(i_stat,i_all,'wrkallred','hamiltonianapplication')
  endif

end subroutine HamiltonianApplication


subroutine hpsitopsi(iter,parallel,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,kbounds,&
     eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
     psi,psit,hpsi,psidst,hpsidst,nspin,spinar)

  use module_types

  implicit none
  include 'mpif.h'
  real(kind=8), parameter :: eps_mach=1.d-12
  logical, intent(in) :: parallel
  integer, intent(in) :: iter,iproc,nproc,n1,n2,n3,norb,norbp,ncong,mids,idsx
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,nspin
  real(kind=8), intent(in) :: hgrid,energy,energy_old
  real(kind=8), dimension(norb), intent(in) :: occup,eval,spinar
  real(kind=8), intent(inout) :: alpha
  real(kind=8), intent(inout) :: gnrm,scprsum
  real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
  real(kind=8), dimension(:,:,:), pointer :: psidst,hpsidst,ads
  type(kinetic_bounds), intent(in) :: kbounds
  type(wavefunctions_descriptors), intent(in) :: wfd
  !local variables
  integer :: ierr,ind,i1,i2,iorb,k,norbu,norbd
  real(kind=8) :: tt,scpr,dnrm2,scprpart

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, orthoconstraint...'
  end if

  !Calculate no. up and dw orbitals for spin-polarized starting guess
  norbu=0
  norbd=0
  do iorb=1,norb
     if(spinar(iorb)>0.0d0) norbu=norbu+1
     if(spinar(iorb)<0.0d0) norbd=norbd+1
  end do
  !write(*,'(1x,a,3i4,30f6.2)')'Spins: ',norb,norbu,norbd,(spinar(iorb),iorb=1,norb)

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  if (parallel) then
     !transpose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !here hpsi is the transposed array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition


     if(nspin==1) then
        call  orthoconstraint_p(iproc,nproc,norb,occup,nvctrp,psit,hpsi,scprsum)
     else
        call  orthoconstraint_p(iproc,nproc,norbu,occup,nvctrp,psit,hpsi,scprsum)
        scprpart=0.0d0
        if(norbd>0) then
           scprpart=scprsum 
           call  orthoconstraint_p(iproc,nproc,norbd,occup,nvctrp,psit(1,norbu+1),hpsi(1,norbu+1),scprsum)
        end if
        scprsum=scprsum+scprpart
     end if

     !retranspose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !here hpsi is the direct array
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,psi,hpsi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     if(nspin==1) then
        call orthoconstraint(norb,occup,nvctrp,psi,hpsi,scprsum)
     else
        call orthoconstraint(norbu,occup,nvctrp,psi,hpsi,scprsum)
        scprpart=0.0d0
        if(norbd>0) then
           scprpart=scprsum 
           call orthoconstraint(norbd,occup,nvctrp,psi(1,norbu+1),hpsi(1,norbu+1),scprsum)
        end if
        scprsum=scprsum+scprpart
     end if
  endif


  ! norm of gradient
  gnrm=0.d0
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     !calculate the address to start from for calculating the 
     !norm of the residue if hpsi is allocated in the transposed way
     !it can be eliminated when including all this procedure in a subroutine
     if (parallel) then
        ind=1+(wfd%nvctr_c+7*wfd%nvctr_f)*(iorb-iproc*norbp-1)
        i1=mod(ind-1,nvctrp)+1
        i2=(ind-i1)/nvctrp+1
     else
        i1=1
        i2=iorb-iproc*norbp
     end if
     scpr=dnrm2(wfd%nvctr_c+7*wfd%nvctr_f,hpsi(i1,i2),1)
     !scpr=dnrm2(nvctr_c+7*nvctr_f,hpsi(1,iorb-iproc*norbp),1)
     !lines for writing the residue following the orbitals
     !if (iorb <=5) write(83,*)iter,iorb,scpr
     gnrm=gnrm+scpr**2
  enddo
  if (parallel) then
     tt=gnrm
     call MPI_ALLREDUCE(tt,gnrm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif
  gnrm=sqrt(gnrm/real(norb,kind=8))

  call timing(iproc,'Precondition  ','ON')
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, preconditioning...'
  end if

  !   write(*,'(10f10.6)') (eval(iorb),iorb=1,norb)
  ! Preconditions all orbitals belonging to iproc
  call preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid, &
       ncong,wfd%nseg_c,wfd%nseg_f,wfd%nvctr_c,wfd%nvctr_f,wfd%keyg,wfd%keyv,eval,&
       kbounds%ibyz_c,kbounds%ibxz_c,kbounds%ibxy_c,&
       kbounds%ibyz_f,kbounds%ibxz_f,kbounds%ibxy_f,hpsi)

  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if
  call timing(iproc,'Precondition  ','OF')

  !apply the minimization method (DIIS or steepest descent)
  if (idsx.gt.0) then
     if (parallel) then
        !transpose the hpsi wavefunction into the diis array

        !here psi is used as a work array
        call timing(iproc,'Un-TransSwitch','ON')
        call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi)
        call timing(iproc,'Un-TransSwitch','OF')
        call timing(iproc,'Un-TransComm  ','ON')
        call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             hpsidst(1,1,mids),nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-TransComm  ','OF')
        !end of transposition

        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           do k=1,nvctrp
              psidst(k,iorb,mids)= psit(k,iorb) 
           enddo
        enddo

        call diisstp(parallel,norb,norbp,nproc,iproc,  &
             ads,iter,mids,idsx,nvctrp,psit,psidst,hpsidst)
     else
        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           do k=1,nvctrp
              psidst(k,iorb,mids)= psi(k,iorb)
              hpsidst(k,iorb,mids)=hpsi(k,iorb)
           enddo
        enddo

        call diisstp(parallel,norb,norbp,nproc,iproc,  &
             ads,iter,mids,idsx,nvctrp,psi,psidst,hpsidst)

     endif
     !    write(*,*) 'psi update done',iproc

  else

     ! update all wavefunctions with the preconditioned gradient
     if (energy.gt.energy_old) then
        alpha=max(.125d0,.5d0*alpha)
        if (alpha.eq..125d0) write(*,*) 'Convergence problem or limit'
     else
        alpha=min(1.05d0*alpha,1.d0)
     endif
     if (iproc.eq.0) write(*,'(1x,a,1pe11.3)') 'alpha=',alpha

     if (parallel) then
        !transpose the hpsi wavefunction

        !here psi is used as a work array
        call timing(iproc,'Un-TransSwitch','ON')
        call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi)
        call timing(iproc,'Un-TransSwitch','OF')
        !here hpsi is the transposed array
        call timing(iproc,'Un-TransComm  ','ON')
        call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-TransComm  ','OF')
        !end of transposition

        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psit(1,iorb),1)
        enddo
     else
        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psi(1,iorb),1)
        enddo
     endif

  endif

  call timing(iproc,'Diis          ','OF')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Orthogonalization...'
  end if

  if (parallel) then
     if(nspin==1) then
        call orthon_p(iproc,nproc,norb,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit)
     else
        call orthon_p(iproc,nproc,norbu,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit)
        if(norbd>0) then
           call orthon_p(iproc,nproc,norbd,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit(1,norbu+1))
        end if
     end if
     !       call checkortho_p(iproc,nproc,norb,nvctrp,psit)

     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     if(nspin==1) then
        call orthon(norb,nvctrp,psi)
        !          call checkortho(norb,nvctrp,psi)
     else
        call orthon(norbu,nvctrp,psi)
        !          call checkortho(norbu,nvctrp,psi)
        if(norbd>0) then
           call orthon(norbd,nvctrp,psi(1,norbu+1))
        end if
        !          call checkortho(norbd,nvctrp,psi(1,norbu+1))
     end if
  endif

  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if

end subroutine hpsitopsi

subroutine first_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
     nspin,psi,hpsi,psit)
  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
  real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
  !local variables
  integer :: i_all,i_stat,ierr

  if (parallel) then
     !allocate hpsi array (used also as transposed)
     !allocated in the transposed way such as 
     !it can also be used as the transposed hpsi
     allocate(hpsi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','first_orthon')

     !transpose the psi wavefunction
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,hpsi)
     call timing(iproc,'Un-TransSwitch','OF')
     !allocate transposed principal wavefunction
     allocate(psit(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psit))*kind(psit),'psit','first_orthon')
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition

     if(nspin==1) then
        call orthon_p(iproc,nproc,norb,nvctrp,nvctr_c+7*nvctr_f,psit) 
     else
        call orthon_p(iproc,nproc,norbu,nvctrp,nvctr_c+7*nvctr_f,psit) 
        if(norbd>0) then
           call orthon_p(iproc,nproc,norbd,nvctrp,nvctr_c+7*nvctr_f,psit(1,norbu+1)) 
        end if
     end if
     !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

     !retranspose the psit wavefunction into psi
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     if(nspin==1) then
        call orthon(norb,nvctrp,psi)
        !call checkortho(norb,norbp,nvctrp,psi)
     else
        call orthon(norbu,nvctrp,psi)
        !call checkortho(norbu,nvctrp,psi)
        if(norbd>0) then
           call orthon(norbd,nvctrp,psi(1,norbu+1))
           !call checkortho(norbd,nvctrp,psi(1,norbu+1))
        end if
     end if
     !allocate hpsi array
     allocate(hpsi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','first_orthon')
  endif

end subroutine first_orthon

! transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
     nspin,psi,hpsi,psit,occup,evsum,eval)
  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), intent(out) :: evsum
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
  !local variables
  integer :: i_all,i_stat,ierr,iorb,jorb
  real(kind=8) :: evpart

  if (parallel) then

     !transpose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !here hpsi is the transposed array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition

     if(nspin==1) then
        call KStrans_p(iproc,nproc,norb,norbp*nproc,nvctrp,occup,hpsi,psit,evsum,eval)
     else
        call KStrans_p(iproc,nproc,norbu,norbu,nvctrp,occup,hpsi,psit,evsum,eval)
        evpart=evsum
        if(norbd>0) then
           call KStrans_p(iproc,nproc,norbd,norbd,nvctrp,occup,&
                hpsi(1,norbu+1),psit(1,norbu+1),evsum,eval(norbu+1))
           evsum=evsum+evpart
        end if
     end if

     !retranspose the psit wavefunction into psi
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition

     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit','last_orthon')

  else
     if(nspin==1) then
        call KStrans(norb,nvctrp,occup,hpsi,psi,evsum,eval)
     else
        call KStrans(norbu,nvctrp,occup,hpsi,psi,evsum,eval)
        evpart=evsum
        if(norbd>0) then
           call KStrans(norbd,nvctrp,occup,hpsi(1,norbu+1),psi(1,norbu+1),&
                evsum,eval(norbu+1))
           evsum=evsum+evpart
        end if
     end if
  endif

  !print the found eigenvalues
  if (iproc == 0) then
     write(*,'(1x,a)')&
          '-------------------------------------------------------------- Kohn-Sham Eigenvalues'
     if (nspin==1) then
        do iorb=1,norb!/2
!!$           jorb=norb/2+iorb
!!$           write(*,'(1x,a,i4,a,1x,1pe21.14,17x,a,i4,a,1x,1pe21.14)') &
!!$                'eval(',iorb,')=',eval(iorb),'eval(',jorb,')=',eval(jorb)
!!$        end do
!!$        if (2*norb/2 /= norb) then
!!$           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',norb/2+1,')=',eval(norb/2+1)
!!$        end if
           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
        end do
     else
        do iorb=1,min(norbu,norbd)
           jorb=norbu+iorb
           write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
                'eval(',iorb,',u)=',eval(iorb),'eval(',iorb,',d)=',eval(jorb)
        end do
        if (norbu > norbd) then
           do iorb=norbd+1,norbu
              write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,',u)=',eval(iorb)
           end do
        else if (norbd > norbu) then
           do iorb=2*norbu+1,norbu+norbd
              write(*,'(48x,a,i0,a,1x,1pe21.14)') 'eval(',iorb-norbu,',d)=',eval(iorb)
           end do
        end if
     end if
  end if

  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi','last_orthon')

end subroutine last_orthon
