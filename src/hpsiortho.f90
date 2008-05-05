subroutine HamiltonianApplication(geocode,iproc,nproc,at,hx,hy,hz,&
     norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,bounds,nlpspd,proj,&
     ngatherarr,ndimpot,potential,psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,nspinor,spinsgn)
  use module_base
  use module_types
  implicit none
  include 'mpif.h'
  type(atoms_data), intent(in) :: at
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(convolutions_bounds), intent(in) :: bounds
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,ndimpot
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspin,nspinor
  real(kind=8), intent(in) :: hx,hy,hz
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(kind=8), dimension(norb), intent(in) :: occup,spinsgn
  real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor*norbp), intent(in) :: psi
  real(kind=8), dimension(max(ndimpot,1),nspin), intent(in), target :: potential
  real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor*norbp), intent(out) :: hpsi
  !local variables
  character(len=*), parameter :: subname='HamiltonianApplication'
  integer :: i_all,i_stat,ierr,iorb,n1i,n2i,n3i
  integer :: nw1,nw2,nsoffset,oidx,ispin,md
  real(kind=8) :: ekin,epot,eproj
  real(kind=8) :: m00,m11,m12,m13,m34,m24,dnrm2
  real(kind=8), dimension(3,2) :: wrkallred
  real(kind=8), dimension(:,:), allocatable :: w1,w2,psir
  !for the periodic BC case, these arrays substitute 
  !psifscf,psifscfk,psig,ww respectively
  real(kind=8), dimension(:,:,:,:), allocatable ::x_c,y_c,x_f1,x_f2,x_f3
  real(kind=8), dimension(:,:,:,:,:), allocatable ::x_f,x_fc,y_f
  real(kind=8), dimension(:,:), pointer :: pot
  real(kind=8), dimension(:,:,:), allocatable :: mom_vec

  call timing(iproc,'ApplyLocPotKin','ON')

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  select case(geocode)
     case('F')
        n1i=2*n1+31
        n2i=2*n2+31
        n3i=2*n3+31

        !dimensions of work arrays
        ! shrink convention: nw1>nw2
        nw1=max((n3+1)*(2*n1+31)*(2*n2+31),&
             (n1+1)*(2*n2+31)*(2*n3+31),&
             2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
             2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))

        nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
             4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
             (n1+1)*(n2+1)*(2*n3+31),&
             (2*n1+31)*(n2+1)*(n3+1))

        !allocation of work arrays
        allocate(y_c(0:n1,0:n2,0:n3,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,y_c,'y_c',subname)
        allocate(y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,y_f,'y_f',subname)
        allocate(x_c(0:n1,0:n2,0:n3,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_c,'x_c',subname)
        allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_f,'x_f',subname)
        allocate(w1(nw1,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,w1,'w1',subname)
        allocate(w2(nw2,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,w2,'w2',subname)
        allocate(x_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_f1,'x_f1',subname)
        allocate(x_f2(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_f2,'x_f2',subname)
        allocate(x_f3(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_f3,'x_f3',subname)

        !initialisation of the work arrays
        call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*nspinor,x_f1)
        call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*nspinor,x_f2)
        call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*nspinor,x_f3)
        call razero((n1+1)*(n2+1)*(n3+1)*nspinor,x_c)
        call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*nspinor,x_f)
        call razero((n1+1)*(n2+1)*(n3+1)*nspinor,y_c)
        call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*nspinor,y_f)

!!$        call razero(nw1*nspinor,w1)
!!$        call razero(nw2*nspinor,w2)

     case('S')
        n1i=2*n1+2
        n2i=2*n2+31
        n3i=2*n3+2
     case('P')
        n1i=2*n1+2
        n2i=2*n2+2
        n3i=2*n3+2

        !allocation of work arrays
        allocate(x_c(n1i,n2i,n3i,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_c,'x_c',subname)
        allocate(y_c(n1i,n2i,n3i,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,y_c,'y_c',subname)
        allocate(x_f1(n1i,n2i,n3i,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_f1,'x_f1',subname)
        allocate(x_f2(n1i,n2i,n3i,nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,x_f2,'x_f2',subname)

  end select

  !then build the potential on the whole simulation box
  if (nproc > 1) then
     allocate(pot(n1i*n2i*n3i,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)

     do ispin=1,nspin
        call MPI_ALLGATHERV(potential(1,ispin),ndimpot,&
             MPI_DOUBLE_PRECISION,pot(1,ispin),ngatherarr(0,1),&
             ngatherarr(0,2),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     end do
  else
     pot => potential
  end if

  ! Wavefunction in real space
  allocate(psir(n1i*n2i*n3i,nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  call razero(n1i*n2i*n3i*nspinor,psir)

  ekin_sum=0.d0
  epot_sum=0.d0
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

     if(spinsgn(iorb)>0.0d0) then
        nsoffset=1
     else
        nsoffset=2
     end if
     if(nspinor==4) nsoffset=1

     oidx=(iorb-1)*nspinor+1-iproc*norbp*nspinor
     select case(geocode)
        case('F')
           call applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,0, &
                hx,wfd%nseg_c,wfd%nseg_f,wfd%nvctr_c,wfd%nvctr_f,wfd%keyg,wfd%keyv,&
                bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c,&
                bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f,y_c,y_f,psir, &
                psi(1,oidx),pot(1,nsoffset),hpsi(1,oidx),epot,ekin,&
                x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
                bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
                bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
                bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
                bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,nw1,nw2,bounds%ibyyzz_r,&
                nspinor)
        case('P')
           call applylocpotkinone_per(n1,n2,n3,hx,hy,hz,wfd%nseg_c,wfd%nseg_f,&
                wfd%nvctr_c,wfd%nvctr_f,wfd%keyg,wfd%keyv,& 
                psir,x_c,y_c,x_f1,x_f2,psi(1,oidx),pot(1,nsoffset),&
                hpsi(1,oidx),epot,ekin) 
        end select

     ekin_sum=ekin_sum+occup(iorb)*ekin
     epot_sum=epot_sum+occup(iorb)*epot

  enddo

  !deallocations of work arrays
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  i_all=-product(shape(x_c))*kind(x_c)
  deallocate(x_c,stat=i_stat)
  call memocc(i_stat,i_all,'x_c',subname)
  i_all=-product(shape(y_c))*kind(y_c)
  deallocate(y_c,stat=i_stat)
  call memocc(i_stat,i_all,'y_c',subname)
  i_all=-product(shape(x_f1))*kind(x_f1)
  deallocate(x_f1,stat=i_stat)
  call memocc(i_stat,i_all,'x_f1',subname)
  i_all=-product(shape(x_f2))*kind(x_f2)
  deallocate(x_f2,stat=i_stat)
  call memocc(i_stat,i_all,'x_f2',subname)

  if (geocode == 'F') then
     i_all=-product(shape(x_f3))*kind(x_f3)
     deallocate(x_f3,stat=i_stat)
     call memocc(i_stat,i_all,'x_f3',subname)
     i_all=-product(shape(y_f))*kind(y_f)
     deallocate(y_f,stat=i_stat)
     call memocc(i_stat,i_all,'y_f',subname)
     i_all=-product(shape(x_f))*kind(x_f)
     deallocate(x_f,stat=i_stat)
     call memocc(i_stat,i_all,'x_f',subname)
     i_all=-product(shape(w1))*kind(w1)
     deallocate(w1,stat=i_stat)
     call memocc(i_stat,i_all,'w1',subname)
     i_all=-product(shape(w2))*kind(w2)
     deallocate(w2,stat=i_stat)
     call memocc(i_stat,i_all,'w2',subname)
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

  eproj_sum=0.d0
  ! loop over all my orbitals
  do iorb=iproc*norbp*nspinor+1,min((iproc+1)*norbp,norb)*nspinor
     call applyprojectorsone(at%ntypes,at%nat,at%iatype,at%psppar,at%npspcode, &
          nlpspd%nprojel,nlpspd%nproj,nlpspd%nseg_p,nlpspd%keyg_p,nlpspd%keyv_p,nlpspd%nvctr_p,&
          proj,wfd%nseg_c,wfd%nseg_f,wfd%keyg,wfd%keyv,wfd%nvctr_c,wfd%nvctr_f,  & 
          psi(1,iorb-iproc*norbp*nspinor),hpsi(1,iorb-iproc*norbp*nspinor),eproj)
        eproj_sum=eproj_sum+occup((iorb-1)/nspinor+1)*eproj
  enddo

  call timing(iproc,'ApplyProj     ','OF')

  if (nproc > 1) then
     wrkallred(1,2)=ekin_sum 
     wrkallred(2,2)=epot_sum 
     wrkallred(3,2)=eproj_sum
     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     ekin_sum=wrkallred(1,1)
     epot_sum=wrkallred(2,1)
     eproj_sum=wrkallred(3,1) 
  endif

!!$  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
!!$     oidx=(iorb-1)*nspinor+1-iproc*norbp*nspinor
!!$     !loop over the spinorial components
!!$     do ispin=oidx,oidx+nspinor-1
!!$        m00=dnrm2(wfd%nvctr_c+7*wfd%nvctr_f,hpsi(1,ispin),1)
!!$        if (iproc == 0) then
!!$           write(*,'(1x,i0,1pe24.17)'),iorb,m00
!!$        end if
!!$     end do
!!$  enddo


end subroutine HamiltonianApplication


subroutine hpsitopsi(geocode,iter,iproc,nproc,norb,norbp,occup,hx,hy,hz,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,kbounds,&
     eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
     psi,psit,hpsi,psidst,hpsidst,nspin,nspinor,spinsgn)
  use module_base
  use module_types
  implicit none
  type(kinetic_bounds), intent(in) :: kbounds
  type(wavefunctions_descriptors), intent(in) :: wfd
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iter,iproc,nproc,n1,n2,n3,norb,norbp,ncong,mids,idsx
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,nspin,nspinor
  real(kind=8), intent(in) :: hx,hy,hz,energy,energy_old
  real(kind=8), dimension(norb), intent(in) :: occup,eval,spinsgn
  real(kind=8), intent(inout) :: alpha
  real(kind=8), intent(inout) :: gnrm,scprsum
  real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
  real(kind=8), dimension(:,:,:), pointer :: psidst,hpsidst,ads
  !local variables
  include 'mpif.h'
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: ierr,ind,i1,i2,iorb,k,norbu,norbd,i_stat,i_all,oidx,sidx
  real(kind=8) :: tt,scpr,dnrm2,scprpart,cprecr
  real(kind=8), dimension(:), allocatable :: zalpha

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, orthoconstraint...'
  end if

  !Calculate no. up and dw orbitals for spin-polarized starting guess
  norbu=0
  norbd=0
  do iorb=1,norb
     if(spinsgn(iorb)>0.0d0) norbu=norbu+1
     if(spinsgn(iorb)<0.0d0) norbd=norbd+1
  end do

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  if (nproc > 1) then
     !transpose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     !here hpsi is the transposed array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition
  else
     psit => psi
     if(nspinor==4) then
        call psitransspi(nvctrp,norb,psit,.true.)
        call psitransspi(nvctrp,norb,hpsi,.true.)
     end if
  end if


  if(nspin==1.or.nspinor==4) then
     call orthoconstraint_p(iproc,nproc,norb,occup,nvctrp,psit,hpsi,scprsum,nspinor)
  else
     call orthoconstraint_p(iproc,nproc,norbu,occup,nvctrp,psit,hpsi,scprsum,nspinor)
     scprpart=0.0d0
     if(norbd>0) then
        scprpart=scprsum 
        call orthoconstraint_p(iproc,nproc,norbd,occup(norbu+1),nvctrp,&
             psit(1,norbu+1),hpsi(1,norbu+1),scprsum,nspinor)
     end if
     scprsum=scprsum+scprpart
  end if
  if (nproc > 1) then
     !retranspose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(hpsi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,  &
          psi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !here hpsi is the direct array
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,psi,hpsi,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     if(nspinor==4) call psitransspi(nvctrp,norb,hpsi,.false.)
!    if(nspinor==4) call psitransspi(nvctrp,norb,psit.true.)
  end if

  ! norm of gradient
  gnrm=0.d0

  call timing(iproc,'Precondition  ','ON')
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, preconditioning...'
  end if
  ! Preconditions all orbitals belonging to iproc
  !and calculate the norm of the residue
  !loop over the orbitals
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     oidx=(iorb-1)*nspinor+1-iproc*norbp*nspinor
     !loop over the spinorial components
     do sidx=oidx,oidx+nspinor-1
        !starting address of the direct array in transposed form
        call trans_address(nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,1,sidx,i1,i2)

        scpr=dnrm2(wfd%nvctr_c+7*wfd%nvctr_f,hpsi(i1,i2),1)
        !lines for writing the residue following the orbitals
        !if (iorb <=5) write(83,*)iter,iorb,scpr
        gnrm=gnrm+scpr**2
        select case(geocode)
        case('F')
           !in this case the grid spacings are uniform
           cprecr=-eval(iorb)
           if(scpr /=0.d0) then
              call precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
                   wfd%nseg_c,wfd%nvctr_c,wfd%nseg_f,wfd%nvctr_f,wfd%keyg,wfd%keyv, &
                   ncong,cprecr,hx,kbounds%ibyz_c,kbounds%ibxz_c,kbounds%ibxy_c,&
                   kbounds%ibyz_f,kbounds%ibxz_f,kbounds%ibxy_f,hpsi(i1,i2))
           end if
        case('P')
           cprecr=0.5d0
           call prec_fft(n1,n2,n3, &
                wfd%nseg_c,wfd%nvctr_c,wfd%nseg_f,wfd%nvctr_f,wfd%keyg,wfd%keyv, &
                ncong,cprecr,hx,hy,hz,hpsi(i1,i2))
        end select
     end do
  enddo
  if (nproc > 1) then
     tt=gnrm
     call MPI_ALLREDUCE(tt,gnrm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif
  gnrm=sqrt(gnrm/real(norb,kind=8))

  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if
  call timing(iproc,'Precondition  ','OF')

  !apply the minimization method (DIIS or steepest descent)
  if (idsx.gt.0) then
     if (nproc > 1) then
        !transpose the hpsi wavefunction into the diis array
        !here psi is used as a work array
        call timing(iproc,'Un-TransSwitch','ON')
        call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi,nspinor)
        call timing(iproc,'Un-TransSwitch','OF')
        call timing(iproc,'Un-TransComm  ','ON')
        call MPI_ALLTOALL(psi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,  &
             hpsidst(1,1,mids),nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-TransComm  ','OF')
        !end of transposition
        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb*nspinor
           do k=1,nvctrp
              psidst(k,iorb,mids)= psit(k,iorb)
           enddo
        enddo
     else
        if(nspinor==4) call psitransspi(nvctrp,norb,hpsi,.true.)
        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb*nspinor
           do k=1,nvctrp
              psidst(k,iorb,mids)= psit(k,iorb)
              hpsidst(k,iorb,mids)=hpsi(k,iorb)
           enddo
        enddo
     endif

     call diisstp(norb,norbp,nproc,iproc, nspinor,  &
          ads,iter,mids,idsx,nvctrp,psit,psidst,hpsidst)

  else

     ! update all wavefunctions with the preconditioned gradient
     if (energy.gt.energy_old) then
        alpha=max(.125d0,.5d0*alpha)
        if (alpha.eq..125d0) write(*,*) 'Convergence problem or limit'
     else
        alpha=min(1.05d0*alpha,1.d0)
     endif
     if (iproc.eq.0) write(*,'(1x,a,1pe11.3)') 'alpha=',alpha

     if (nproc > 1) then
        !transpose the hpsi wavefunction
        !here psi is used as a work array
        call timing(iproc,'Un-TransSwitch','ON')
        call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi,nspinor)
        call timing(iproc,'Un-TransSwitch','OF')
        !here hpsi is the transposed array
        call timing(iproc,'Un-TransComm  ','ON')
        call MPI_ALLTOALL(psi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,  &
             hpsi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-TransComm  ','OF')
        !end of transposition
     else if(nspinor==4) then
        call psitransspi(nvctrp,norb,hpsi,.true.)
     endif

     call timing(iproc,'Diis          ','ON')
     do iorb=1,norb*nspinor
        call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psit(1,iorb),1)
     enddo

  endif

  call timing(iproc,'Diis          ','OF')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Orthogonalization...'
  end if

  if(nspin==1.or.nspinor==4) then
     call orthon_p(iproc,nproc,norb,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,nspinor)
  else
     call orthon_p(iproc,nproc,norbu,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,nspinor)
     if(norbd>0) then
        call orthon_p(iproc,nproc,norbd,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit(1,norbu+1),nspinor)
     end if
  end if
    !       call checkortho_p(iproc,nproc,norb,nvctrp,psit)

  if (nproc > 1) then
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     nullify(psit)
     if(nspinor==4) call psitransspi(nvctrp,norb,psi,.false.)

!     if(nspinor==4) call psitransspi(nvctrp,norb,hpsi,.false.) !(needed?)

  endif

  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
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



subroutine first_orthon(iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
     nspin,psi,hpsi,psit)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
  real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
  !local variables
  character(len=*), parameter :: subname='first_orthon'
  include 'mpif.h'
  integer :: i_all,i_stat,ierr,nspinor,iorb

  if(nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if
  
  if (nproc > 1) then
     !allocate hpsi array (used also as transposed)
     !allocated in the transposed way such as 
     !it can also be used as the transposed hpsi
     allocate(hpsi(nvctrp,nspinor*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,hpsi,'hpsi',subname)

!     write(*,'(a,i3,30f10.5)') 'SWI',iproc,(sum(psi(1:nvctrp,iorb)),iorb=1,norbp*nspinor)
     !transpose the psi wavefunction
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,hpsi,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     !allocate transposed principal wavefunction
     allocate(psit(nvctrp,nspinor*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psit,'psit',subname)
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(hpsi,nvctrp*nspinor*norbp,MPI_DOUBLE_PRECISION,  &
          psit,nvctrp*nspinor*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition
  else
     psit => psi
  end if

  if(nspin==1.or.nspinor==4) then
     call orthon_p(iproc,nproc,norb,nvctrp,nvctr_c+7*nvctr_f,psit,nspinor) 
  else
     call orthon_p(iproc,nproc,norbu,nvctrp,nvctr_c+7*nvctr_f,psit,nspinor) 
     if(norbd>0) then
        call orthon_p(iproc,nproc,norbd,nvctrp,nvctr_c+7*nvctr_f,psit(1,norbu+1),nspinor) 
     end if
  end if
  !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

  if (nproc > 1) then
     !retranspose the psit wavefunction into psi
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*nspinor*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*nspinor*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     nullify(psit)
     !allocate hpsi array
     allocate(hpsi(nvctr_c+7*nvctr_f,nspinor*norbp+ndebug),stat=i_stat)
     call memocc(i_stat,hpsi,'hpsi',subname)
  end if

!!$  else
!!$     if(nspin==1) then
!!$        call orthon(norb,nvctrp,psi)
!!$        !call checkortho(norb,norbp,nvctrp,psi)
!!$     else
!!$        call orthon(norbu,nvctrp,psi)
!!$        !call checkortho(norbu,nvctrp,psi)
!!$        if(norbd>0) then
!!$           call orthon(norbd,nvctrp,psi(1,norbu+1))
!!$           !call checkortho(norbd,nvctrp,psi(1,norbu+1))
!!$        end if
!!$     end if
!!$     !allocate hpsi array
!!$     allocate(hpsi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
!!$     call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','first_orthon')
!!$  endif

end subroutine first_orthon

! transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
     nspin,psi,hpsi,psit,occup,evsum,eval)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), intent(out) :: evsum
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
  !local variables
  character(len=*), parameter :: subname='last_orthon'
  include 'mpif.h'
  integer :: i_all,i_stat,ierr,iorb,jorb,nspinor
  real(kind=8) :: evpart

  if(nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if

  if (nproc > 1) then
     !transpose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     !here hpsi is the transposed array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition
  else
     psit => psi
     if(nspinor==4) then
        call psitransspi(nvctrp,norb,psit,.true.)
        call psitransspi(nvctrp,norb,hpsi,.true.)
     end if
  end if

  if(nspin==1.or.nspinor==4) then
     call KStrans_p(iproc,nproc,norb,nvctrp,occup,hpsi,psit,evsum,eval,nspinor)
  else
     call KStrans_p(iproc,nproc,norbu,nvctrp,occup,hpsi,psit,evsum,eval,nspinor)
     evpart=evsum
     if(norbd>0) then
        call KStrans_p(iproc,nproc,norbd,nvctrp,occup(norbu+1),&
             hpsi(1,norbu+1),psit(1,norbu+1),evsum,eval(norbu+1),nspinor)
        evsum=evsum+evpart
     end if
  end if

  if (nproc > 1) then
     !retranspose the psit wavefunction into psi
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp*nspinor,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition

     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  else
     if(nspinor==4) then
        call psitransspi(nvctrp,norb,psit,.false.)
     end if
     nullify(psit)
  end if

!!$  else
!!$     if(nspin==1) then
!!$        call KStrans(norb,nvctrp,occup,hpsi,psi,evsum,eval)
!!$     else
!!$        call KStrans(norbu,nvctrp,occup,hpsi,psi,evsum,eval)
!!$        evpart=evsum
!!$        if(norbd>0) then
!!$           call KStrans(norbd,nvctrp,occup(norbu+1),hpsi(1,norbu+1),psi(1,norbu+1),&
!!$                evsum,eval(norbu+1))
!!$           evsum=evsum+evpart
!!$        end if
!!$     end if
!!$  endif

  !print the found eigenvalues
  if (iproc == 0) then
     write(*,'(1x,a)')&
          '-------------------------------------------------------------- Kohn-Sham Eigenvalues'
     if (nspin==1.or.nspinor==4) then
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
  call memocc(i_stat,i_all,'hpsi',subname)

end subroutine last_orthon


subroutine calc_moments(iproc,nproc,norb,norbp,nvctr,nspinor,psi)
  use module_base
  implicit none
  include 'mpif.h'
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr,nspinor
  real(kind=8), dimension(nvctr,norbp*nproc*nspinor), intent(in) :: psi
  !local variables
  character(len=*), parameter :: subname='calc_moments'
  integer :: i_all,i_stat,ierr,iorb
  integer :: oidx,ispin,md,ndim
  real(kind=8) :: m00,m11,m13,m24,m12,m34,m14,m23
  real(kind=8), dimension(:,:,:), allocatable :: mom_vec
  real(kind=8) :: ddot

  ndim=2
  if (nproc==1) ndim=1

  if(nspinor==4) then
     
     allocate(mom_vec(4,norbp*nproc,ndim+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)
     call razero(4*norbp*nproc*ndim,mom_vec)
     
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        oidx=(iorb-1)*nspinor+1-iproc*norbp*nspinor
        m00=ddot(2*nvctr,psi(1,oidx),1,psi(1,oidx),1)
        m11=ddot(2*nvctr,psi(1,oidx+2),1,psi(1,oidx+2),1)
        m13=ddot(nvctr,psi(1,oidx),1,psi(1,oidx+2),1)
        m24=ddot(nvctr,psi(1,oidx+1),1,psi(1,oidx+3),1)
!        m12=ddot(nvctr,psi(1,oidx),1,psi(1,oidx+1),1)
!        m34=ddot(nvctr,psi(1,oidx+2),1,psi(1,oidx+3),1)
        m14=ddot(nvctr,psi(1,oidx),1,psi(1,oidx+3),1)
        m23=ddot(nvctr,psi(1,oidx+1),1,psi(1,oidx+2),1)
        mom_vec(1,iorb-iproc*norbp,1)=(m00+m11) !rho
        mom_vec(2,iorb-iproc*norbp,1)=2.0d0*(m13+m24)       !m_x
!        mom_vec(3,iorb-iproc*norbp,1)=2.0d0*(m12-m34)       !m_y
        mom_vec(3,iorb-iproc*norbp,1)=2.0d0*(m14-m23)       !m_y
        mom_vec(4,iorb-iproc*norbp,1)=(m00-m11) !m_z
     end do
     
     if(nproc>1) then
           call MPI_GATHER(mom_vec(1,1,1),4*norbp,MPI_DOUBLE_PRECISION,mom_vec(1,1,2),4*norbp, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        oidx=2
     else
        oidx=1
     end if
     
     if(iproc==0) then
        write(*,'(1x,a)')&
             'Magnetic polarization per orbital'
        write(*,'(1x,a)')&
             '  iorb    m_x       m_y       m_z'
        do iorb=1,norb
           write(*,'(1x,i5,3f10.5)') &
                iorb,(mom_vec(md,iorb,oidx)/mom_vec(1,iorb,oidx),md=2,4)
        end do
     end if
     
     i_all=-product(shape(mom_vec))*kind(mom_vec)
     deallocate(mom_vec,stat=i_stat)
     call memocc(i_stat,i_all,'mom_vec',subname)
     

     
  end if

end subroutine calc_moments
