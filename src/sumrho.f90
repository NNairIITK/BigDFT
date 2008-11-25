subroutine sumrho(iproc,nproc,norb,norbp,lr,ixc,hxh,hyh,hzh,occup,  & 
     psi,rho,nrho,nscatterarr,nspin,nspinor,spinsgn,hybrid_on)
  ! Calculates the charge density by summing the square of all orbitals
  ! Input: psi
  ! Output: rho
  use module_base!, only: gp,dp,wp,ndebug,memocc
  use module_types
  implicit none
  logical, intent(in) :: hybrid_on
  integer, intent(in) :: iproc,nproc,norb,norbp,nrho,nspin,nspinor,ixc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors) :: lr 
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,norbp*nspinor), intent(in) :: psi
  real(dp), dimension(max(nrho,1),nspin), intent(out), target :: rho
  real(gp), dimension(norb), intent(in) :: occup,spinsgn
  !local variables
  character(len=*), parameter :: subname='sumrho'
  logical :: rsflag
  integer :: nw1,nw2,nrhotot,n3d,n1i,n2i,n3i,nxc,nxf,itmred
  integer :: ind1,ind2,ind3,ind1s,ind2s,ind3s,oidx,sidx,nspinn
  integer :: i00,i0,i1,i2,i3,i3off,i3s,isjmp,i,ispin,iorb,jproc,i_all,i_stat,ierr,j3,j3p,j
  real(dp) :: charge,tt
  real(wp), dimension(:,:), allocatable :: tmred
  real(dp), dimension(:,:), pointer :: rho_p

  call timing(iproc,'Rho_comput    ','ON')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Calculation of charge density...'
  end if

!!$  select case(geocode)
!!$  case('F')
!!$     n1i=2*n1+31
!!$     n2i=2*n2+31
!!$     n3i=2*n3+31
!!$  case('S')
!!$     n1i=2*n1+2
!!$     n2i=2*n2+31
!!$     n3i=2*n3+2
!!$  case('P')
!!$     n1i=2*n1+2
!!$     n2i=2*n2+2
!!$     n3i=2*n3+2
!!$  end select
!!$
!!$  call create_Glr(geocode,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i,wfd,bounds,&
!!$       lr)

  !flag for toggling the REDUCE_SCATTER stategy
  rsflag=.not. (ixc >= 11 .and. ixc <=16)

  !calculate dimensions of the complete array to be allocated before the reduction procedure
  if (rsflag) then
     nrhotot=0
     do jproc=0,nproc-1
        nrhotot=nrhotot+nscatterarr(jproc,1)
     end do
  else
     nrhotot=lr%d%n3i
  end if
  nspinn=max(nspin,nspinor)
  if (nproc > 1) then
     allocate(rho_p(lr%d%n1i*lr%d%n2i*nrhotot,nspinn+ndebug),stat=i_stat)
     call memocc(i_stat,rho_p,'rho_p',subname)
  else
     rho_p => rho
  end if

  !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
  if(nspinor==4) then 
     call razero(lr%d%n1i*lr%d%n2i*nrhotot*nspinor,rho_p)
     call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot,rho_p,nproc)
  else
     call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,rho_p,nproc)
  end if

  !for each of the orbitals treated by the processor build the partial densities
  call local_partial_density(iproc,nproc,hybrid_on,rsflag,nscatterarr,&
     nrhotot,lr,hxh,hyh,hzh,nspin,nspinor,norbp,norb,occup,spinsgn,psi,rho_p)

  !the density must be communicated to meet the shape of the poisson solver
  if (nproc > 1) then
     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     if (rsflag) then
        do ispin=1,nspin
          call MPI_REDUCE_SCATTER(rho_p(1,ispin),rho(1,ispin),&
               lr%d%n1i*lr%d%n2i*nscatterarr(:,1),&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        end do
     else
         call MPI_ALLREDUCE(MPI_IN_PLACE,rho_p,lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin,&
              MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
         !stop 'rsflag active in sumrho.f90, check MPI2 implementation'
     end if
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
     if (.not. rsflag) then
        !treatment which includes the periodic GGA
        !the density should meet the poisson solver distribution
        i3s=nscatterarr(iproc,3)-nscatterarr(iproc,4)
        n3d=nscatterarr(iproc,1)
        do ispin=1,nspin
           do i3=1,n3d
              j3=i3+i3s
              j3p=modulo(j3-1,lr%d%n3i)+1
              do i2=1,lr%d%n2i
                 do i1=1,lr%d%n1i
                    i=i1+(i2-1)*lr%d%n1i+lr%d%n1i*lr%d%n2i*(i3-1)
                    j=i1+(i2-1)*lr%d%n1i+lr%d%n1i*lr%d%n2i*(j3p-1)
                    rho(i,ispin)=rho_p(j,ispin)
                 end do
              end do
           end do
        end do
     end if
  end if

  ! Check
  tt=0.d0
  i3off=lr%d%n1i*lr%d%n2i*nscatterarr(iproc,4)

  !allocation of the magnetic density orientation array
  if (nproc > 1) then
     itmred=2
  else
     itmred=1
  end if

  !use this check also for the magnetic density orientation
  allocate(tmred(nspin+1,itmred+ndebug),stat=i_stat)
  call memocc(i_stat,tmred,'tmred',subname)

!!$  if(nspinor==4) then
!!$     nspinn=1
!!$  else
!!$     nspinn=nspin
!!$  end if


  !print *,'nspin,iproc',nspin,nspinn,nspinor,iproc
  

  tmred(nspin+1,itmred)=0.0_dp
  do ispin=1,nspin!n
     tmred(ispin,itmred)=0.0_dp
     do i=1,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2)
!!$        tt=tt+rho(i+i3off,ispin)
        tmred(ispin,itmred)=tmred(ispin,itmred)+rho(i+i3off,ispin)
!!$        !temporary check for debugging purposes
!!$        if (rho(i+i3off,ispin) < 9.d-21) then
!!$           print *,iproc,'error in density construction',rho(i+i3off,ispin)
!!$        end if
     enddo
     tmred(nspin+1,itmred)=tmred(nspin+1,itmred)+tmred(ispin,itmred)
  end do

  if (nproc > 1) then
     i_all=-product(shape(rho_p))*kind(rho_p)
     deallocate(rho_p,stat=i_stat)
     call memocc(i_stat,i_all,'rho_p',subname)

     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
!!$     call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(tmred(1,2),tmred(1,1),nspin+1,mpidtypd,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
  else
     !useless, only for completeness
     nullify(rho_p)

     !charge=tt
  end if

  !write the results
  if (iproc == 0) then
     if(nspin==4) then
        charge=tmred(1,1)
        tt=sqrt(tmred(2,1)**2+tmred(3,1)**2+tmred(4,1)**2)
     else
        charge=0._dp
        do ispin=1,nspin
           charge=charge+tmred(ispin,1)
        end do
     end if
     write(*,'(1x,a,f21.12)')&
          'done. Total electronic charge=',real(charge,gp)*hxh*hyh*hzh
     if(nspin == 4 .and. tt > 0._dp)&
          write(*,'(a,5f10.4)')'  Magnetic density orientation:',&
          (tmred(ispin,1)/tmred(1,1),ispin=2,nspin)
  end if
  
  i_all=-product(shape(tmred))*kind(tmred)
  deallocate(tmred,stat=i_stat)
  call memocc(i_stat,i_all,'tmred',subname)

  call timing(iproc,'Rho_comput    ','OF')

end subroutine sumrho

!here starts the routine for building partial density inside the localisation region
!this routine should be treated as a building-block for the linear scaling code
subroutine local_partial_density(iproc,nproc,hybrid_on,rsflag,nscatterarr,&
     nrhotot,lr,hxh,hyh,hzh,nspin,nspinor,norbp,norb,occup,spinsgn,psi,rho_p)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  logical, intent(in) :: hybrid_on,rsflag
  integer, intent(in) :: iproc,nproc,nrhotot
  integer, intent(in) :: nspin,nspinor,norb,norbp
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(gp), dimension(norb), intent(in) :: occup,spinsgn
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,norbp*nspinor), intent(in) :: psi
  real(dp), dimension(lr%d%n1i,lr%d%n2i,nrhotot,max(nspin,nspinor)), intent(inout) :: rho_p
  !local variables
  character(len=*), parameter :: subname='local_partial_density'
  integer :: nw1,nw2,nxc,nxf,iorb
  integer :: n1,n2,n3,n1i,n2i,n3i,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer :: oidx,sidx,nspinn
  integer :: i_all,i_stat,ierr,j3,j3p,j,i
  real(gp) :: hfac
  real(wp), dimension(0:3) :: scal
  real(wp), dimension(:,:), allocatable :: psir
  real(wp), dimension(:), allocatable :: x_c_psifscf,x_f_psig,w1,w2

  do i=0,3
     scal(i)=1.0_wp
  enddo

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i
  nfl1=lr%d%nfl1
  nfl2=lr%d%nfl2
  nfl3=lr%d%nfl3
  nfu1=lr%d%nfu1
  nfu2=lr%d%nfu2
  nfu3=lr%d%nfu3
  
  select case(lr%geocode)
  case('F')
     !dimension of the work arrays
     ! shrink convention: nw1>nw2
     nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
          (n1+1)*(2*n2+31)*(2*n3+31),&
          2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
          2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
     nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
          4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
          (n1+1)*(n2+1)*(2*n3+31),&
          (2*n1+31)*(n2+1)*(n3+1))
     nxc=(n1+1)*(n2+1)*(n3+1)!(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
  case('S')
     !dimension of the work arrays
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+31)*(2*n3+2)
     nxf=1
  case('P')
     if (hybrid_on) then
        ! hybrid case:
        nxc=(n1+1)*(n2+1)*(n3+1)
        nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

        nw1=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*n1+2),(2*n1+2)*(n2+2)*(n3+2))
        nw1=max(nw1,2*(n3+1)*(n1+1)*(n2+1))	   ! for the comb_shrink_hyb_c
        nw1=max(nw1,4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1)) ! for the _f

        nw2=max(2*(nfu3-nfl3+1)*(2*n1+2)*(2*n2+2),(n3+1)*(2*n1+2)*(2*n2+2))
        nw2=max(nw2,4*(n2+1)*(n3+1)*(n1+1))	! for the comb_shrink_hyb_c   
        nw2=max(nw2,2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1)) ! for the _f
     else
        !dimension of the work arrays, fully periodic case
        nw1=1
        nw2=1
        nxc=(2*n1+2)*(2*n2+2)*(2*n3+2)
        nxf=1
     endif

  end select
  !work arrays
  allocate(x_c_psifscf(nxc+ndebug),stat=i_stat)
  call memocc(i_stat,x_c_psifscf,'x_c_psifscf',subname)
  allocate(x_f_psig(nxf+ndebug),stat=i_stat)
  call memocc(i_stat,x_f_psig,'x_f_psig',subname)
  allocate(w1(nw1+ndebug),stat=i_stat)
  call memocc(i_stat,w1,'w1',subname)
  allocate(w2(nw2+ndebug),stat=i_stat)
  call memocc(i_stat,w2,'w2',subname)


  ! Wavefunction in real space
  nspinn=max(nspin,nspinor)
  allocate(psir(n1i*n2i*n3i,nspinn+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call razero(nxc,x_c_psifscf)
     call razero(nxf,x_f_psig)
     call razero(n1i*n2i*n3i*nspinn,psir)
  end if

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     hfac=(occup(iorb)/(hxh*hyh*hzh))

     oidx=(iorb-iproc*norbp-1)*nspinor

     if (hfac /= 0.d0) then

        select case(lr%geocode)
        case('F')

           do sidx=1,nspinor
              call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
                   lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),  & 
                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                   lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1), &
                   scal,psi(1,oidx+sidx),psi(lr%wfd%nvctr_c+1,oidx+sidx),x_c_psifscf,x_f_psig)

              call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,&
                   x_c_psifscf,x_f_psig,  & 
                   psir(1,sidx),lr%bounds%kb%ibyz_c,lr%bounds%gb%ibzxx_c,lr%bounds%gb%ibxxyy_c,&
                   lr%bounds%gb%ibyz_ff,lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f,lr%bounds%ibyyzz_r)
           end do

           call partial_density(rsflag,nproc,n1i,n2i,n3i,nspinor,nspinn,nrhotot,&
                hfac,nscatterarr,spinsgn(iorb),psir,rho_p,lr%bounds%ibyyzz_r)

        case('P')

           do sidx=1,nspinor
              if (hybrid_on) then
                 ! hybrid case
                 call uncompress_per_f_short(n1,n2,n3,lr%wfd%nseg_c,&
                      lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),& 
                      lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                      lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1), &
                      psi(1,oidx+sidx),psi(lr%wfd%nvctr_c+1,oidx+sidx),x_c_psifscf,x_f_psig,&
                      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

                 call comb_grow_all_hybrid(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nw1,nw2&
                      ,w1,w2,x_c_psifscf,x_f_psig,psir(1,sidx),lr%bounds%gb)
              else
                 call uncompress_per(n1,n2,n3,lr%wfd%nseg_c,&
                      lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
                      lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                      lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1),&
                      psi(1,oidx+sidx),psi(lr%wfd%nvctr_c+1,oidx+sidx),x_c_psifscf,psir(1,sidx))

                 call convolut_magic_n_per_self(2*n1+1,2*n2+1,2*n3+1,&
                      x_c_psifscf,psir(1,sidx)) 
              endif

           end do

           call partial_density(rsflag,nproc,n1i,n2i,n3i,nspinor,nspinn,nrhotot,&
                hfac,nscatterarr,spinsgn(iorb),psir,rho_p)

        case('S')

           do sidx=1,nspinor
              call uncompress_slab(n1,n2,n3,lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                   lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
                   lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keyg(1,lr%wfd%nseg_c+1),&
                   lr%wfd%keyv(lr%wfd%nseg_c+1),   &
                   psi(1,oidx+sidx),psi(lr%wfd%nvctr_c+1,oidx+sidx),x_c_psifscf,psir(1,sidx))

              call convolut_magic_n_slab_self(2*n1+1,2*n2+15,2*n3+1,x_c_psifscf,psir(1,sidx)) 
           end do

           call partial_density(rsflag,nproc,n1i,n2i,n3i,nspinor,nspinn,nrhotot,&
                hfac,nscatterarr,spinsgn(iorb),psir,rho_p)
        end select
     end if
  enddo

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)
  i_all=-product(shape(x_c_psifscf))*kind(x_c_psifscf)
  deallocate(x_c_psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'x_c_psifscf',subname)
  i_all=-product(shape(x_f_psig))*kind(x_f_psig)
  deallocate(x_f_psig,stat=i_stat)
  call memocc(i_stat,i_all,'x_f_psig',subname)
  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1',subname)
  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2',subname)

end subroutine local_partial_density


subroutine partial_density(rsflag,nproc,n1i,n2i,n3i,nspinor,nspinn,nrhotot,&
     hfac,nscatterarr,spinsgn,psir,rho_p,&
     ibyyzz_r) !optional argument
  use module_base
  use module_types
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinor,nspinn
  real(gp), intent(in) :: hfac,spinsgn
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
  real(wp), dimension(n1i,n2i,n3i,nspinn), intent(in) :: psir
  real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
  integer, dimension(:,:,:), pointer, optional :: ibyyzz_r 
  !local variables
  integer :: i3s,jproc,i3off,n3d,isjmp,i1,i2,i3,i1s,i1e,j3
  real(gp) :: hfac2
  real(dp) :: psisq,p1,p2,p3,p4,r1,r2,r3,r4
  !sum different slices by taking into account the overlap
  i3s=0
  hfac2=2.0_gp*hfac

  !case without bounds
  i1s=1
  i1e=n1i

  loop_xc_overlap: do jproc=0,nproc-1
     !case for REDUCE_SCATTER approach, not used for GGA since it enlarges the 
     !communication buffer
     if (rsflag) then
        i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
        n3d=nscatterarr(jproc,1)
        if (n3d==0) exit loop_xc_overlap
     else
        i3off=0
        n3d=n3i
     end if
     !here the condition for the MPI_ALLREDUCE should be entered
     if(spinsgn > 0.0d0) then
        isjmp=1
     else
        isjmp=2
     end if
     do i3=i3off+1,i3off+n3d
        !this allows the presence of GGA with non-isolated BC. If i3 is between 1 and n3i
        !j3=i3. This is useful only when dealing with rsflags and GGA, so we can comment it out
        !j3=modulo(i3-1,n3i)+1 
        j3=i3
        i3s=i3s+1
        do i2=1,n2i
           !this if statement is inserted here for avoiding code duplication
           !it is to be seen whether the code results to be too much unoptimised
           if (present(ibyyzz_r)) then
              i1s=ibyyzz_r(1,i2-15,j3-15)+1
              i1e=ibyyzz_r(2,i2-15,j3-15)+1
           end if
           if (nspinor == 1) then
              do i1=i1s,i1e
                 !conversion between the different types
                 psisq=real(psir(i1,i2,j3,1),dp)
                 psisq=psisq*psisq
                 rho_p(i1,i2,i3s,isjmp)=rho_p(i1,i2,i3s,isjmp)+real(hfac,dp)*psisq
              end do
           else  !similar loop for nspinor=4
              do i1=i1s,i1e
                 !conversion between the different types
                 p1=real(psir(i1,i2,j3,1),dp)
                 p2=real(psir(i1,i2,j3,2),dp)
                 p3=real(psir(i1,i2,j3,3),dp)
                 p4=real(psir(i1,i2,j3,4),dp)

                 !density values
                 r1=p1*p1+p2*p2+p3*p3+p4*p4
                 r2=p1*p3+p2*p4
                 r3=p1*p4-p2*p3
                 r4=p1*p1+p2*p2-p3*p3-p4*p4

                 rho_p(i1,i2,i3s,1)=rho_p(i1,i2,i3s,1)+real(hfac,dp)*r1
                 rho_p(i1,i2,i3s,2)=rho_p(i1,i2,i3s,2)+real(hfac2,dp)*r2
                 rho_p(i1,i2,i3s,3)=rho_p(i1,i2,i3s,3)+real(hfac2,dp)*r3
                 rho_p(i1,i2,i3s,4)=rho_p(i1,i2,i3s,4)+real(hfac,dp)*r4
              end do
           end if
        end do
     end do
     if (.not. rsflag) exit loop_xc_overlap !the whole range is already done
  end do loop_xc_overlap

  if (i3s /= nrhotot) then
     write(*,'(1x,a,i0,1x,i0)')'ERROR: problem with rho_p: i3s,nrhotot,',i3s,nrhotot
     stop
  end if

end subroutine partial_density


