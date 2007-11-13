subroutine sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
     wfd,psi,rho,nrho,nscatterarr,nspin,spinar,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
  ! Calculates the charge density by summing the square of all orbitals
  ! Input: psi
  ! Output: rho

  use module_types

  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(in) :: bounds
  logical, intent(in) ::  parallel
  integer, intent(in) :: iproc,nproc,norb,norbp,nrho,nspin
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(kind=8), intent(in) :: hgrid
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8), dimension(norb), intent(in) :: occup,spinar
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
  real(kind=8), dimension(max(nrho,1),nspin), intent(out) :: rho
  !local variables
  logical, parameter :: newmethod=.true.
  integer :: nw1,nw2,nrhotot,n3d
  integer :: ind1,ind2,ind3,ind1s,ind2s,ind3s
  integer :: i00,i0,i1,i2,i3,i3off,i3s,isjmp,i,ispin,iorb,jproc,i_all,i_stat,ierr
  real(kind=8) :: hfac,hgridh,tt,charge
  real(kind=8), dimension(:), allocatable :: psir,rho_p
  real(kind=8), dimension(:,:,:), allocatable :: x_c!input 
  real(kind=8), dimension(:,:,:,:), allocatable :: x_f,x_fc !added for newmethod=.false.
  real(kind=8), allocatable, dimension(:) :: w1,w2,scal
  include 'mpif.h'

  call timing(iproc,'Rho_comput    ','ON')

  hgridh=hgrid*.5d0 

  ! shrink convention: nw1>nw2
  nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
       (n1+1)*(2*n2+31)*(2*n3+31),&
       2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))

  nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
       4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
       (n1+1)*(n2+1)*(2*n3+31),&
	   (2*n1+31)*(n2+1)*(n3+1))

  allocate(scal(0:3),stat=i_stat)
  call memocc(i_stat,product(shape(scal))*kind(scal),'scal','sumrho')

  do i=0,3
     scal(i)=1.d0
  enddo

  allocate(x_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(x_c))*kind(x_c),'x_c','sumrho')
  
  allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)! work
  call memocc(i_stat,product(shape(x_f))*kind(x_f),'x_f','sumrho')
  allocate(w1(nw1),stat=i_stat)
  call memocc(i_stat,product(shape(w1))*kind(w1),'w1','sumrho')
  allocate(w2(nw2),stat=i_stat) ! work
  call memocc(i_stat,product(shape(w2))*kind(w2),'w2','sumrho')

  call razero((n1+1)*(n2+1)*(n3+1),x_c)
  call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f)
  
!  call razero(nw1,w1)
!  call razero(nw2,w2)

  ! Wavefunction in real space
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','sumrho')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Calculation of charge density...'
  end if

  !switch between the old and the new method
  if (newmethod) then
	  
     call razero((2*n1+31)*(2*n2+31)*(2*n3+31),psir)
     if (parallel) then
        !calculate dimensions of the complete array to be allocated before the reduction procedure
        nrhotot=0
        do jproc=0,nproc-1
           nrhotot=nrhotot+nscatterarr(jproc,1)
        end do
        allocate(rho_p((2*n1+31)*(2*n2+31)*nrhotot*nspin),stat=i_stat)
        call memocc(i_stat,product(shape(rho_p))*kind(rho_p),'rho_p','sumrho')

        !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        call tenminustwenty((2*n1+31)*(2*n2+31)*nrhotot*nspin,rho_p,nproc)

        do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
           hfac=(occup(iorb)/hgridh**3)

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
                 isjmp=0
              else
                 isjmp=(2*n1+31)*(2*n2+31)*nrhotot
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
                       rho_p(ind1s+isjmp)=rho_p(ind1s+isjmp)+hfac*psir(ind1)**2
                    end do
                 end do
              end do
           end do loop_xc_overlap

           if (i3s /= nrhotot) then
              print *,'problem with rhopot array in sumrho,i3s,nrhotot,',i3s,nrhotot
              stop
           end if

        enddo

        call timing(iproc,'Rho_comput    ','OF')
        call timing(iproc,'Rho_commun    ','ON')
        call MPI_REDUCE_SCATTER(rho_p,rho,(2*n1+31)*(2*n2+31)*nscatterarr(:,1),&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(nspin>1) then
           call MPI_REDUCE_SCATTER(rho_p((2*n1+31)*(2*n2+31)*nrhotot+1),rho(1,2),&
                (2*n1+31)*(2*n2+31)*nscatterarr(:,1),&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        end if
        call timing(iproc,'Rho_commun    ','OF')
        call timing(iproc,'Rho_comput    ','ON')

        ! Check
        tt=0.d0
        i3off=(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)
        do ispin=1,nspin
           do i=1,(2*n1+31)*(2*n2+31)*nscatterarr(iproc,2)
              tt=tt+rho(i+i3off,ispin)
!!$        !temporary check for debugging purposes
!!$        if (rho(i+i3off) < 9.d-21) then
!!$           print *,iproc,'error in density construction',rho(i+i3off)
!!$        end if
           enddo
        end do
        call timing(iproc,'Rho_comput    ','OF')
        call timing(iproc,'Rho_commun    ','ON')
        call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
             'done. Total electronic charge=',charge*hgridh**3
        call timing(iproc,'Rho_commun    ','OF')
        call timing(iproc,'Rho_comput    ','ON')
        i_all=-product(shape(rho_p))*kind(rho_p)
        deallocate(rho_p,stat=i_stat)
        call memocc(i_stat,i_all,'rho_p','sumrho')

     else
        !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        call tenminustwenty((2*n1+31)*(2*n2+31)*(2*n3+31)*nspin,rho,nproc)
        !     call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

        do iorb=1,norb

           call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
                wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
                scal,psi(1,iorb),psi(wfd%nvctr_c+1,iorb),x_c,x_f)

           call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
                psir,bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
                bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,bounds%ibyyzz_r)


           if(spinar(iorb)>0.0d0) then
              ispin=1
           else
              ispin=2
           end if
           ! i=(i3-1)*(2*n2+31)*(2*n1+31)+(i2-1)*(2*n1+31)+i1
           ! i1=1,2*n1+31
           ! i2=1,2*n2+31
           ! i3=1,2*n3+31
           do i3=1,2*n3+31
              i00=(i3-1)*(2*n2+31)*(2*n1+31)+1
              do i2=1,2*n2+31
                 i0=i00+(i2-1)*(2*n1+31)
                 do i=i0+bounds%ibyyzz_r(1,i2-15,i3-15),i0+bounds%ibyyzz_r(2,i2-15,i3-15)
                    rho(i,ispin)=rho(i,ispin)+(occup(iorb)/hgridh**3)*psir(i)**2
                 enddo
              enddo
           enddo
        enddo
        ! Check
        tt=0.d0
        do ispin=1,nspin
           do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
              tt=tt+rho(i,ispin)
           enddo
        end do
        tt=tt*hgridh**3
        if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
             'done. Total electronic charge=',tt

     endif

  else
     !allocate one array more
     allocate(x_fc(0:n1,0:n2,0:n3,3),stat=i_stat)
     call memocc(i_stat,product(shape(x_fc))*kind(x_fc),'x_fc','sumrho')

     if (parallel) then
        !calculate dimensions of the complete array to be allocated before the reduction procedure
        nrhotot=0
        do jproc=0,nproc-1
           nrhotot=nrhotot+nscatterarr(jproc,1)
        end do
        allocate(rho_p((2*n1+31)*(2*n2+31)*nrhotot*nspin),stat=i_stat)
        call memocc(i_stat,product(shape(rho_p))*kind(rho_p),'rho_p','sumrho')

        !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        call tenminustwenty((2*n1+31)*(2*n2+31)*nrhotot*nspin,rho_p,nproc)
        !call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)

        do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
           hfac=(occup(iorb)/hgridh**3)

           call uncompress_forstandard_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
                wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
                scal,psi(1,iorb-iproc*norbp),psi(wfd%nvctr_c+1,iorb-iproc*norbp),x_c,x_fc,x_f)

           call comb_grow_all_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
                psir,bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
                bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f)

           !sum different slices by taking into account the overlap
           i3s=0
           loop_xc_overlap_prev: do jproc=0,nproc-1
              i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
              n3d=nscatterarr(jproc,1)
              if (n3d==0) exit loop_xc_overlap_prev
              if(spinar(iorb)>0.0d0) then
                 isjmp=0
              else
                 isjmp=(2*n1+31)*(2*n2+31)*nrhotot
              end if
              do i3=i3off+1,i3off+n3d
                 i3s=i3s+1
                 ind3=(i3-1)*(2*n1+31)*(2*n2+31)
                 ind3s=(i3s-1)*(2*n1+31)*(2*n2+31)
                 do i2=1,2*n2+31
                    ind2=(i2-1)*(2*n1+31)+ind3
                    ind2s=(i2-1)*(2*n1+31)+ind3s
                    do i1=1,2*n1+31
                       ind1=i1+ind2
                       ind1s=i1+ind2s
                       rho_p(ind1s+isjmp)=rho_p(ind1s+isjmp)+hfac*psir(ind1)**2
                    end do
                 end do
              end do
           end do loop_xc_overlap_prev

           if (i3s /= nrhotot) then
              print *,'problem with rhopot array in sumrho,i3s,nrhotot,',i3s,nrhotot
              stop
           end if

        enddo

        call timing(iproc,'Rho_comput    ','OF')
        call timing(iproc,'Rho_commun    ','ON')
        call MPI_REDUCE_SCATTER(rho_p,rho,(2*n1+31)*(2*n2+31)*nscatterarr(:,1),&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(nspin>1) then
           call MPI_REDUCE_SCATTER(rho_p((2*n1+31)*(2*n2+31)*nrhotot+1:),rho(1:,2:),&
                (2*n1+31)*(2*n2+31)*nscatterarr(:,1),&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        end if
        call timing(iproc,'Rho_commun    ','OF')
        call timing(iproc,'Rho_comput    ','ON')

        ! Check
        tt=0.d0
        i3off=(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)
        do ispin=1,nspin
           do i=1,(2*n1+31)*(2*n2+31)*nscatterarr(iproc,2)
              tt=tt+rho(i+i3off,ispin)
           enddo
        end do
        call timing(iproc,'Rho_comput    ','OF')
        call timing(iproc,'Rho_commun    ','ON')
        call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
             'done. Total electronic charge=',charge*hgridh**3
        call timing(iproc,'Rho_commun    ','OF')
        call timing(iproc,'Rho_comput    ','ON')
        i_all=-product(shape(rho_p))*kind(rho_p)
        deallocate(rho_p,stat=i_stat)
        call memocc(i_stat,i_all,'rho_p','sumrho')

     else
        !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        call tenminustwenty((2*n1+31)*(2*n2+31)*(2*n3+31)*nspin,rho,nproc)
        !call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

        do iorb=1,norb

           call uncompress_forstandard_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
                wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
                scal,psi(1,iorb-iproc*norbp),psi(wfd%nvctr_c+1,iorb-iproc*norbp),x_c,x_fc,x_f)

           call comb_grow_all_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
                psir,bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
                bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f)

           if(spinar(iorb)>0.0d0) then
              ispin=1
           else
              ispin=2
           end if
           do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
              rho(i,ispin)=rho(i,ispin)+(occup(iorb)/hgridh**3)*psir(i)**2
           enddo
        enddo
        ! Check
        tt=0.d0
        do ispin=1,nspin
           do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
              tt=tt+rho(i,ispin)
           enddo
        end do
        tt=tt*hgridh**3
        if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
             'done. Total electronic charge=',tt

     endif

     i_all=-product(shape(x_fc))*kind(x_fc)
     deallocate(x_fc,stat=i_stat)
     call memocc(i_stat,i_all,'x_fc','sumrho')

  end if

  i_all=-product(shape(scal))*kind(scal)
  deallocate(scal,stat=i_stat)
  call memocc(i_stat,i_all,'scal','sumrho')
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





