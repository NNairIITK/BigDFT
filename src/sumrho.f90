subroutine sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
     nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho,nrho,nscatterarr,nspin,spinar,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
  ! Calculates the charge density by summing the square of all orbitals
  ! Input: psi
  ! Output: rho
  implicit real(kind=8) (a-h,o-z)
  logical parallel,withmpi2
  dimension rho(nrho,nspin),occup(norb),spinar(norb)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f,norbp)
  dimension nscatterarr(0:nproc-1,4)!n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8), allocatable :: psir(:),rho_p(:)
  !***************Alexey**************************************************************************
  real(kind=8), allocatable, dimension(:,:,:) :: x_c!input 
  real*8,allocatable::x_f(:,:,:,:)
  real(kind=8), allocatable, dimension(:) :: w1,w2

  real(kind=8) :: scal(0:3),hfac
  ! for grow:
  integer ibyz_c(2,0:n2,0:n3)
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
  !	for real space:
!  integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
  integer,intent(in):: ibyyzz_r(2,2*n2+31,2*n3+31)
  !***********************************************************************************************
  include 'mpif.h'

  call timing(iproc,'Rho_comput    ','ON')

  hgridh=hgrid*.5d0 

 !***************Alexey**************************************************************************

  ! shrink convention: nw1>nw2
  nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
       (n1+1)*(2*n2+31)*(2*n3+31),&
       2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))

  nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
       4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
       (n1+1)*(n2+1)*(2*n3+31))

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
 !***********************************************************************************************
  ! Wavefunction in real space
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','sumrho')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Calculation of charge density...'
  end if

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
!***************Alexey**************************************************************************

        call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
             nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
             scal,psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),x_c,x_f)

        call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
             psir,ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

!**********************************************************************************************
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
				 do	i1=ibyyzz_r(1,i2,i3)+1,ibyyzz_r(2,i2,i3)+1
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
       call MPI_REDUCE_SCATTER(rho_p((2*n1+31)*(2*n2+31)*nrhotot+1:),rho(1:,2:),(2*n1+31)*(2*n2+31)*nscatterarr(:,1),&
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

!***************Alexey**************************************************************************

        call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
             nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
             scal,psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),x_c,x_f)

        call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
             psir,ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

!***********************************************************************************************
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
					do i=i0+ibyyzz_r(1,i2,i3),i0+ibyyzz_r(2,i2,i3)
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

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir','sumrho')


  !*********************Alexey*********************************************************************
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
  !**********************************************************************************************

  call timing(iproc,'Rho_comput    ','OF')

END SUBROUTINE sumrho





