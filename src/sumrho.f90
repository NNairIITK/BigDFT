subroutine sumrho(geocode,iproc,nproc,norb,norbp,n1,n2,n3,hxh,hyh,hzh,occup,  & 
     wfd,psi,rho,nrho,nscatterarr,nspin,spinar,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
  ! Calculates the charge density by summing the square of all orbitals
  ! Input: psi
  ! Output: rho
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(in) :: bounds
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,norb,norbp,nrho,nspin
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(kind=8), intent(in) :: hxh,hyh,hzh
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8), dimension(norb), intent(in) :: occup,spinar
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
  real(kind=8), dimension(max(nrho,1),nspin), intent(out), target :: rho
  !local variables
  include 'mpif.h'
  integer :: nw1,nw2,nrhotot,n3d,n1i,n2i,n3i
  integer :: ind1,ind2,ind3,ind1s,ind2s,ind3s
  integer :: i00,i0,i1,i2,i3,i3off,i3s,isjmp,i,ispin,iorb,jproc,i_all,i_stat,ierr
  real(kind=8) :: hfac,hgridh,tt,charge
  real(kind=8), dimension(0:3) :: scal
  real(kind=8), dimension(:), allocatable :: psir
  real(kind=8), dimension(:,:,:), allocatable :: x_c!input 
  real(kind=8), dimension(:,:,:,:), allocatable :: x_f !added for newmethod=.false.
  real(kind=8), allocatable, dimension(:) :: w1,w2
  real(kind=8), dimension(:,:), pointer :: rho_p


  call timing(iproc,'Rho_comput    ','ON')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Calculation of charge density...'
  end if

  do i=0,3
     scal(i)=1.d0
  enddo

  select case(geocode)
     case('F')
        n1i=2*n1+31
        n2i=2*n2+31
        n3i=2*n3+31
     case('S')
        n1i=2*n1+2
        n2i=2*n2+31
        n3i=2*n3+2
     case('P')
        n1i=2*n1+2
        n2i=2*n2+2
        n3i=2*n3+2
  end select

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
  !work arrays
  allocate(x_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(x_c))*kind(x_c),'x_c','sumrho')
  allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)! work
  call memocc(i_stat,product(shape(x_f))*kind(x_f),'x_f','sumrho')
  allocate(w1(nw1),stat=i_stat)
  call memocc(i_stat,product(shape(w1))*kind(w1),'w1','sumrho')
  allocate(w2(nw2),stat=i_stat) ! work
  call memocc(i_stat,product(shape(w2))*kind(w2),'w2','sumrho')

  ! Wavefunction in real space
  allocate(psir(n1i*n2i*n3i),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','sumrho')

  !initialisation
  call razero((n1+1)*(n2+1)*(n3+1),x_c)
  call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f)

  call razero(n1i*n2i*n3i,psir)


  !calculate dimensions of the complete array to be allocated before the reduction procedure
  nrhotot=0
  do jproc=0,nproc-1
     nrhotot=nrhotot+nscatterarr(jproc,1)
  end do

  if (nproc > 1) then
     allocate(rho_p(n1i*n2i*nrhotot,nspin),stat=i_stat)
     call memocc(i_stat,product(shape(rho_p))*kind(rho_p),'rho_p','sumrho')
  else
     rho_p => rho
  end if

  !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
  call tenminustwenty(n1i*n2i*nrhotot*nspin,rho_p,nproc)

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     hfac=(occup(iorb)/hxh*hyh*hzh)

     if (hfac /= 0.d0) then

        call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
             wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
             scal,psi(1,iorb-iproc*norbp),psi(wfd%nvctr_c+1,iorb-iproc*norbp),x_c,x_f)


        call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
             psir,bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
             bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,bounds%ibyyzz_r)

        !sum different slices by taking into account the overlap
        i3s=0
        loop_xc_overlap: do jproc=0,nproc-1
           i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
           n3d=nscatterarr(jproc,1)
           if (n3d==0) exit loop_xc_overlap
           if(spinar(iorb)>0.0d0) then
              isjmp=1
           else
              isjmp=2
           end if
           do i3=i3off+1,i3off+n3d
              i3s=i3s+1
              ind3=(i3-1)*(2*n1+31)*(2*n2+31)
              ind3s=(i3s-1)*(2*n1+31)*(2*n2+31)
              do i2=1,2*n2+31
                 ind2=(i2-1)*(2*n1+31)+ind3
                 ind2s=(i2-1)*(2*n1+31)+ind3s
                 !                 do i1=1,2*n1+31
                 do i1=bounds%ibyyzz_r(1,i2-15,i3-15)+1,bounds%ibyyzz_r(2,i2-15,i3-15)+1
                    ind1=i1+ind2
                    ind1s=i1+ind2s
                    !do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
                    rho_p(ind1s,isjmp)=rho_p(ind1s,isjmp)+hfac*psir(ind1)**2
                 end do
              end do
           end do
        end do loop_xc_overlap

        if (i3s /= nrhotot) then
           print *,'problem with rhopot array in sumrho,i3s,nrhotot,',i3s,nrhotot
           stop
        end if
     end if

  enddo

  if (nproc > 1) then
     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     call MPI_REDUCE_SCATTER(rho_p,rho,n1i*n2i*nscatterarr(:,1),&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     if(nspin>1) then
        call MPI_REDUCE_SCATTER(rho_p(1,2),rho(1,2),&
             n1i*n2i*nscatterarr(:,1),&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
  end if

  ! Check
  tt=0.d0
  i3off=n1i*n2i*nscatterarr(iproc,4)
  do ispin=1,nspin
     do i=1,n1i*n2i*nscatterarr(iproc,2)
        tt=tt+rho(i+i3off,ispin)
!!$        !temporary check for debugging purposes
!!$        if (rho(i+i3off) < 9.d-21) then
!!$           print *,iproc,'error in density construction',rho(i+i3off)
!!$        end if
     enddo
  end do

  if (nproc > 1) then
     i_all=-product(shape(rho_p))*kind(rho_p)
     deallocate(rho_p,stat=i_stat)
     call memocc(i_stat,i_all,'rho_p','sumrho')

     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
  else
     !useless, only for completeness
     nullify(rho_p)

     charge=tt
  end if

  if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
       'done. Total electronic charge=',charge*hxh*hyh*hzh

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir','sumrho')
  i_all=-product(shape(x_c))*kind(x_c)
  deallocate(x_c,stat=i_stat)
  call memocc(i_stat,i_all,'x_c','sumrho')
  i_all=-product(shape(x_f))*kind(x_f)
  deallocate(x_f,stat=i_stat)
  call memocc(i_stat,i_all,'x_f','sumrho')
  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1','sumrho')
  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2','sumrho')

  call timing(iproc,'Rho_comput    ','OF')

END SUBROUTINE sumrho





