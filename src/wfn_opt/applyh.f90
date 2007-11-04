subroutine applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
     hgrid,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,&
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
     psi,pot,hpsi,epot_sum,ekin_sum,nspin,spinar,&
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
  !  Applies the local potential and kinetic energy operator to all wavefunctions belonging to processor
  ! Input: pot,psi
  ! Output: hpsi,epot,ekin
  implicit real(kind=8) (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension occup(norb),pot((2*n1+31)*(2*n2+31)*(2*n3+31)*nspin),spinar(norb)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension  psi(nvctr_c+7*nvctr_f,norbp)
  dimension hpsi(nvctr_c+7*nvctr_f,norbp)
  logical :: parameter, newmethod=.false.
  real(kind=8), allocatable :: psir(:),y_c(:,:,:),y_f(:,:,:,:)

  !********************Alexey***************************************************************
  ! for shrink:
  integer ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  ! for grow:
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  ! for real space:
  integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

  real(kind=8), dimension(:,:,:), allocatable ::x_c !input 
  real(kind=8), dimension(:,:,:,:), allocatable::x_f,x_fc ! input
  real(kind=8), dimension(:,:,:), allocatable :: x_f1,x_f2,x_f3 ! input
  real(kind=8), dimension(:), allocatable :: w1,w2 

  ! shrink convention: nw1>nw2
  nw1=max((n3+1)*(2*n1+31)*(2*n2+31),&
       (n1+1)*(2*n2+31)*(2*n3+31),&
       2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))

  nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
       4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
       (n1+1)*(n2+1)*(2*n3+31))

  allocate(y_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(y_c))*kind(y_c),'y_c','applylocpotkinall')
  allocate(y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
  call memocc(i_stat,product(shape(y_f))*kind(y_f),'y_f','applylocpotkinall')
  ! Wavefunction in real space
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','applylocpotkinall')

  allocate(x_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(x_c))*kind(x_c),'x_c','applylocpotkinall')
  allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)! work
  call memocc(i_stat,product(shape(x_f))*kind(x_f),'x_f','applylocpotkinall')
  allocate(w1(nw1),stat=i_stat)
  call memocc(i_stat,product(shape(w1))*kind(w1),'w1','applylocpotkinall')
    
  allocate(w2(nw2),stat=i_stat) ! work
  call memocc(i_stat,product(shape(w2))*kind(w2),'w2','applylocpotkinall')


  if (newmethod) then
     allocate(x_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
     call memocc(i_stat,product(shape(x_f1))*kind(x_f1),'x_f1','applylocpotkinall')
     allocate(x_f2(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
     call memocc(i_stat,product(shape(x_f2))*kind(x_f2),'x_f2','applylocpotkinall')
     allocate(x_f3(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
     call memocc(i_stat,product(shape(x_f3))*kind(x_f3),'x_f3','applylocpotkinall')

     call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f1)
     call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f2)
     call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f3)

     call razero((n1+1)*(n2+1)*(n3+1),x_c)
     call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f)

     !to be initialised
     call razero((n1+1)*(n2+1)*(n3+1),y_c)
     call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),y_f)


     call razero((2*n1+31)*(2*n2+31)*(2*n3+31),psir)

     ekin_sum=0.d0
     epot_sum=0.d0
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        if(spinar(iorb)>0.0d0) then
           nsoffset=1
        else
           nsoffset=(2*n1+31)*(2*n2+31)*(2*n3+31)+1
        end if
        call applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
             hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
             ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,y_c,y_f,psir,  &
             psi(1,iorb-iproc*norbp),pot(nsoffset),hpsi(1,iorb-iproc*norbp),epot,ekin, & 
             x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
             ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
             ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,nw1,nw2,ibyyzz_r)
        ekin_sum=ekin_sum+occup(iorb)*ekin
        epot_sum=epot_sum+occup(iorb)*epot

     enddo

     i_all=-product(shape(x_f1))*kind(x_f1)
     deallocate(x_f1,stat=i_stat)
     call memocc(i_stat,i_all,'x_f1','applylocpotkinall')

     i_all=-product(shape(x_f2))*kind(x_f2)
     deallocate(x_f2,stat=i_stat)
     call memocc(i_stat,i_all,'x_f2','applylocpotkinall')

     i_all=-product(shape(x_f3))*kind(x_f3)
     deallocate(x_f3,stat=i_stat)
     call memocc(i_stat,i_all,'x_f3','applylocpotkinall')

  else

     allocate(x_fc(0:n1,0:n2,0:n3,3),stat=i_stat)
     call memocc(i_stat,product(shape(x_fc))*kind(x_fc),'x_fc','applylocpotkinall')

     ekin_sum=0.d0
     epot_sum=0.d0
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        if(spinar(iorb)>0.0d0) then
           nsoffset=1
        else
           nsoffset=(2*n1+31)*(2*n2+31)*(2*n3+31)+1
        end if
        call applylocpotkinone_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
             hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
             ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,y_c,y_f,psir,  &
             psi(1,iorb-iproc*norbp),pot(nsoffset),hpsi(1,iorb-iproc*norbp),epot,ekin, & 
             x_c,x_fc,x_f,w1,w2,&
             ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
             ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,nw1,nw2,ibyyzz_r)
        ekin_sum=ekin_sum+occup(iorb)*ekin
        epot_sum=epot_sum+occup(iorb)*epot

     enddo

     i_all=-product(shape(x_fc))*kind(x_fc)
     deallocate(x_fc,stat=i_stat)
     call memocc(i_stat,i_all,'x_fc','applylocpotkinall')

  end if
  i_all=-product(shape(y_c))*kind(y_c)
  deallocate(y_c,stat=i_stat)
  call memocc(i_stat,i_all,'y_c','applylocpotkinall')
  
  i_all=-product(shape(y_f))*kind(y_f)
  deallocate(y_f,stat=i_stat)
  call memocc(i_stat,i_all,'y_f','applylocpotkinall')
  
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir','applylocpotkinall')

  i_all=-product(shape(x_c))*kind(x_c)
  deallocate(x_c,stat=i_stat)
  call memocc(i_stat,i_all,'x_c','applylocpotkinall')
 
  i_all=-product(shape(x_f))*kind(x_f)
  deallocate(x_f,stat=i_stat)
  call memocc(i_stat,i_all,'x_f','applylocpotkinall')

  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1','applylocpotkinall')

  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2','applylocpotkinall')

END SUBROUTINE applylocpotkinall

subroutine applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
     hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
     y_c,y_f,psir,  &
     psi,pot,hpsi,epot,ekin,x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,nw1,nw2,ibyyzz_r)!
  !  Applies the local potential and kinetic energy operator to one wavefunction 
  ! Input: pot,psi
  ! Output: hpsi,epot,ekin
  implicit real(kind=8) (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension pot((2*n1+31)*(2*n2+31)*(2*n3+31))
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f),scal(0:3)
  dimension hpsi(nvctr_c+7*nvctr_f)
  dimension y_c(0:n1,0:n2,0:n3)
  dimension y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  dimension psir((2*n1+31)*(2*n2+31)*(2*n3+31))
  !********************Alexey***************************************************************
  ! for shrink:
  integer ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  ! for grow:
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  ! for real space:
  integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
  !*****************************************************************************************
  real(kind=8) x_c(0:n1,0:n2,0:n3),x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)! input
    real(kind=8)::x_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
    real(kind=8)::x_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3)
    real(kind=8)::x_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2)
  real(kind=8) w1(nw1),w2(nw2) ! work
  !***********************************************************************************************
  do i=0,3
     scal(i)=1.d0
  enddo

!!$  tt=0.d0
!!$  do i=1,nvctr_c+7*nvctr_f
!!$     tt=tt+psi(i)**2
!!$  end do
!!$  print *,'before applylocpotkinone',tt

  
  call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
       nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       scal,psi(1),psi(nvctr_c+1),x_c,x_f,x_f1,x_f2,x_f3)

  call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
       psir,ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

!!$  tt=0.d0
!!$  do i3=0,n3
!!$     do i2=0,n2
!!$        do i1=0,n1
!!$           tt=tt+x_c(i1,i2,i3)**2
!!$        end do
!!$     end do
!!$  end do
!!$  do i3=nfl3,nfu3
!!$     do i2=nfl2,nfu2
!!$        do i1=nfl1,nfu1
!!$           do i=1,7
!!$              tt=tt+x_f(i,i1,i2,i3)**2
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  print *,'grow',tt,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3



  if (nbuf.eq.0) then
     call realspace(ibyyzz_r,pot,psir,epot,n1,n2,n3)
  else
     call realspace_nbuf(ibyyzz_r,pot,psir,epot,n1,n2,n3,nbuf)
  endif

!!$  tt=0.d0
!!$  do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
!!$     tt=tt+psir(i)**2
!!$  end do
!!$  print *,'psir',tt

  call comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,psir,&
       ibxy_c,ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,y_c,y_f)!,ibyz_c,ibyz_f)

!!$  tt=0.d0
!!$  do i3=0,n3
!!$     do i2=0,n2
!!$        do i1=0,n1
!!$           tt=tt+y_c(i1,i2,i3)**2
!!$        end do
!!$     end do
!!$  end do
!!$  do i3=nfl3,nfu3
!!$     do i2=nfl2,nfu2
!!$        do i1=nfl1,nfu1
!!$           do i=1,7
!!$              tt=tt+y_f(i,i1,i2,i3)**2
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  print *,'shrink',tt

  call ConvolkineticT(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,ekin,x_f1,x_f2,x_f3)


!!$  tt=0.d0
!!$  do i3=0,n3
!!$     do i2=0,n2
!!$        do i1=0,n1
!!$           tt=tt+y_c(i1,i2,i3)**2
!!$        end do
!!$     end do
!!$  end do
!!$  do i3=nfl3,nfu3
!!$     do i2=nfl2,nfu2
!!$        do i1=nfl1,nfu1
!!$           do i=1,7
!!$              tt=tt+y_f(i,i1,i2,i3)**2
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  print *,'kinetic',tt

     
  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       scal,y_c,y_f,hpsi(1),hpsi(nvctr_c+1))


!!$  tt=0.d0
!!$  do i=1,nvctr_c+7*nvctr_f
!!$     tt=tt+hpsi(i)**2
!!$  end do
!!$  print *,'after applylocpotkinone',tt


END SUBROUTINE applylocpotkinone

subroutine realspace(ibyyzz_r,pot,psir,epot,n1,n2,n3)
  implicit none
  integer,intent(in)::n1,n2,n3
  integer,intent(in)::ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

  real(kind=8),intent(in)::pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
  real(kind=8),intent(inout)::psir(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

  real(kind=8),intent(out)::epot
  real(kind=8) tt
  integer i1,i2,i3

  epot=0.d0
  do i3=-14,2*n3+16
     do i2=-14,2*n2+16
        do i1=max(ibyyzz_r(1,i2,i3)-14,-14),min(ibyyzz_r(2,i2,i3)-14,2*n1+16)
           tt=pot(i1,i2,i3)*psir(i1,i2,i3)
           epot=epot+tt*psir(i1,i2,i3)
           psir(i1,i2,i3)=tt
        enddo
     enddo
  enddo

end subroutine realspace

subroutine realspace_nbuf(ibyyzz_r,pot,psir,epot,nb1,nb2,nb3,nbuf)
  implicit none
  integer,intent(in)::nb1,nb2,nb3,nbuf
  integer,intent(in)::ibyyzz_r(2,-14:2*nb2+16,-14:2*nb3+16)
  real(kind=8),intent(in)::pot(-14:2*nb1+16-4*nbuf,-14:2*nb2+16-4*nbuf,-14:2*nb3+16-4*nbuf)
  real(kind=8),intent(inout)::psir(-14:2*nb1+16,-14:2*nb2+16,-14:2*nb3+16)

  real(kind=8),intent(out)::epot
  real(kind=8) tt,dnrm2
  integer i1,i2,i3

  epot=0.d0
  do i3=-14,2*nb3+16
     if (i3.ge.-14+2*nbuf .and. i3.le.2*nb3+16-2*nbuf) then
        do i2=-14,2*nb2+16
           if (i2.ge.-14+2*nbuf .and. i2.le.2*nb2+16-2*nbuf) then
              do i1=-14+2*nbuf,ibyyzz_r(1,i2,i3)-14-1
                 psir(i1,i2,i3)=0.d0
              enddo
              do i1=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf),min(ibyyzz_r(2,i2,i3)-14,2*nb1+16-2*nbuf)
                 tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf)*psir(i1,i2,i3)
                 epot=epot+tt*psir(i1,i2,i3)
                 psir(i1,i2,i3)=tt
              enddo
              do i1=ibyyzz_r(2,i2,i3)-14+1,2*nb1+16-2*nbuf
                 psir(i1,i2,i3)=0.d0
              enddo
           else
              do i1=-14,2*nb1+16
                 psir(i1,i2,i3)=0.d0
              enddo
           endif
        enddo
     else
        do i2=-14,2*nb2+16
           do i1=-14,2*nb1+16
              psir(i1,i2,i3)=0.d0
           enddo
        enddo
     endif
  enddo

end subroutine realspace_nbuf


subroutine applyprojectorsall(iproc,ntypes,nat,iatype,psppar,npspcode,occup, &
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
     norb,norbp,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)
  ! Applies all the projectors onto a wavefunction
  ! Input: psi_c,psi_f
  ! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initilized to zero at the beginning)
  implicit real(kind=8) (a-h,o-z)
  dimension psppar(0:4,0:6,ntypes),iatype(nat),npspcode(ntypes)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f,norbp),hpsi(nvctr_c+7*nvctr_f,norbp)
  dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
  dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
  dimension proj(nprojel),occup(norb)

  eproj_sum=0.d0
  ! loop over all my orbitals
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     call applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
          nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
          nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,  & 
          psi(1,iorb-iproc*norbp),hpsi(1,iorb-iproc*norbp),eproj)
     eproj_sum=eproj_sum+occup(iorb)*eproj
     !     write(*,*) 'iorb,eproj',iorb,eproj
  enddo

END SUBROUTINE applyprojectorsall

subroutine applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
     nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj)
  ! Applies all the projectors onto a single wavefunction
  ! Input: psi_c,psi_f
  ! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initilized to zero at the beginning)
  implicit real(kind=8) (a-h,o-z)
  dimension psppar(0:4,0:6,ntypes),iatype(nat),npspcode(ntypes)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f),hpsi(nvctr_c+7*nvctr_f)
  dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
  dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
  dimension proj(nprojel)

  ! loop over all projectors
  iproj=0
  eproj=0.d0
  istart_c=1
  do iat=1,nat
     mbseg_c=nseg_p(2*iat-1)-nseg_p(2*iat-2)
     mbseg_f=nseg_p(2*iat  )-nseg_p(2*iat-1)
     jseg_c=nseg_p(2*iat-2)+1
     jseg_f=nseg_p(2*iat-1)+1
     mbvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
     mbvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)
     ityp=iatype(iat)
     !GTH and HGH pseudopotentials
     do l=1,4
        do i=1,3
           if (psppar(l,i,ityp).ne.0.d0) then
              do m=1,2*l-1
                 iproj=iproj+1
                 istart_f=istart_c+mbvctr_c
                 call wpdot(  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),scpr)

                 ! test (will sometimes give wrong result)
                 call wpdot(  &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),tcpr)
                 if (scpr.ne.tcpr) then
                    print *,'projectors: scpr.ne.tcpr'
                    print *,'l,i,m,h_i^l=',l,i,m,psppar(l,i,ityp)
                    print *,'scpr,tcpr',scpr,tcpr
                    stop 
                 end if
                 ! testend

                 scprp=scpr*psppar(l,i,ityp)
                 eproj=eproj+scprp*scpr

                 call waxpy(&
                      scprp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                 istart_c=istart_f+7*mbvctr_f
              enddo
              if (npspcode(ityp) == 3 .and. l/=4 .and. i/=3) then !HGH case, offdiagonal terms
                 loop_j: do j=i+1,3
                    if (psppar(l,j,ityp) .eq. 0.d0) exit loop_j
                    !calculate the coefficients for the off-diagonal terms
                    if (l==1) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(3.d0/5.d0)
                          if (j==3) offdiagcoeff=0.5d0*sqrt(5.d0/21.d0)
                       else
                          offdiagcoeff=-0.5d0*sqrt(100.d0/63.d0)
                       end if
                    else if (l==2) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(5.d0/7.d0)
                          if (j==3) offdiagcoeff=1.d0/6.d0*sqrt(35.d0/11.d0)
                       else
                          offdiagcoeff=-7.d0/3.d0*sqrt(1.d0/11.d0)
                       end if
                    else if (l==3) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(7.d0/9.d0)
                          if (j==3) offdiagcoeff=0.5d0*sqrt(63.d0/143.d0)
                       else
                          offdiagcoeff=-9.d0*sqrt(1.d0/143.d0)
                       end if
                    end if
                    istart_c_i=istart_c-(2*l-1)*(mbvctr_c+7*mbvctr_f)
                    istart_c_j=istart_c_i+(j-i)*(2*l-1)*(mbvctr_c+7*mbvctr_f)
                    do m=1,2*l-1
                       !starting addresses of the projectors
                       istart_f_j=istart_c_j+mbvctr_c
                       istart_f_i=istart_c_i+mbvctr_c
                       call wpdot(&
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                            mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_j),proj(istart_f_j),scpr_j)

                       call wpdot(&
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                            mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_i),proj(istart_f_i),scpr_i)


                       scprp_j=scpr_j*offdiagcoeff*psppar(l,j,ityp)
                       scprp_i=scpr_i*offdiagcoeff*psppar(l,j,ityp)
                       !scpr_i*h_ij*scpr_j+scpr_j*h_ij*scpr_i
                       eproj=eproj+2.d0*scprp_j*scpr_i

                       !|hpsi>=|hpsi>+h_ij (<p_i|psi>|p_j>+<p_j|psi>|p_i>)
                       call waxpy(&
                            scprp_j,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                            keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_i),proj(istart_f_i),  &
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                       call waxpy(&
                            scprp_i,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                            keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_j),proj(istart_f_j),  &
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                       istart_c_j=istart_f_j+7*mbvctr_f
                       istart_c_i=istart_f_i+7*mbvctr_f
                    enddo
                 end do loop_j
              else if (npspcode(ityp) == 10 .and. i/=3) then !HGH-K case, offdiagonal terms
                 loop_jK: do j=i+1,3
                    if (psppar(l,j,ityp) .eq. 0.d0) exit loop_jK
                    istart_c_i=istart_c-(2*l-1)*(mbvctr_c+7*mbvctr_f)
                    istart_c_j=istart_c_i+(j-i)*(2*l-1)*(mbvctr_c+7*mbvctr_f)
                    do m=1,2*l-1
                       !starting addresses of the projectors
                       istart_f_j=istart_c_j+mbvctr_c
                       istart_f_i=istart_c_i+mbvctr_c
                       call wpdot(&
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                            mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_j),proj(istart_f_j),scpr_j)

                       call wpdot(&
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                            mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_i),proj(istart_f_i),scpr_i)

                       !scpr_i*h_ij*scpr_j+scpr_j*h_ij*scpr_i (with symmetric h_ij)
                       eproj=eproj+2.d0*scpr_i*psppar(l,i+j+1,ityp)*scpr_j
                       scprp_j=scpr_j*psppar(l,i+j+1,ityp)
                       scprp_i=scpr_i*psppar(l,i+j+1,ityp)

                       !|hpsi>=|hpsi>+h_ij (<p_i|psi>|p_j>+<p_j|psi>|p_i>)
                       call waxpy(&
                            scprp_j,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                            keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_i),proj(istart_f_i),  &
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                       call waxpy(&
                            scprp_i,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                            keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_j),proj(istart_f_j),  &
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                       istart_c_j=istart_f_j+7*mbvctr_f
                       istart_c_i=istart_f_i+7*mbvctr_f
                    enddo
                 end do loop_jK
              end if
           end if
        enddo
     enddo
  enddo
  if (iproj.ne.nproj) stop '1:applyprojectorsone'
  if (istart_c-1.ne.nprojel) stop '2:applyprojectorsone'
  return
END SUBROUTINE applyprojectorsone

