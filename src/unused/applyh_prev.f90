subroutine applylocpotkinone_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
     hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
     y_c,y_f,psir,  &
     psi,pot,hpsi,epot,ekin,x_c,x_fc,x_f,w1,w2,&
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,nw1,nw2,ibyyzz_r)
  !  Applies the local potential and kinetic energy operator to one wavefunction 
  ! Input: pot,psi
  ! Output: hpsi,epot,ekin
  use module_base
  implicit real(kind=8) (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension pot((2*n1+31)*(2*n2+31)*(2*n3+31))
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f),scal(0:3)
  dimension hpsi(nvctr_c+7*nvctr_f)
  dimension y_c(0:n1,0:n2,0:n3)
  dimension y_f(7,0:n1,0:n2,0:n3)
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
  real(kind=8) x_c(0:n1,0:n2,0:n3),  x_fc(0:n1,0:n2,0:n3,3), x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)! input
  real(kind=8) w1(nw1),w2(nw2) ! work
  !***********************************************************************************************
  character(len=*), parameter :: subname='applylocpotkinone_prev'
  do i=0,3
     scal(i)=1.d0
  enddo

  call uncompress_forstandard_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
       nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       scal,psi(1),psi(nvctr_c+1),x_c,x_fc,x_f)

  call comb_grow_all_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
       psir,ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f)

  if (nbuf.eq.0) then
     call realspace_prev(ibyyzz_r,pot,psir,epot,n1,n2,n3)
  else
     call realspace_nbuf_prev(ibyyzz_r,pot,psir,epot,n1,n2,n3,nbuf)
  endif

  call comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,psir,&
       ibxy_c,ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,y_c,y_f)!,ibyz_c,ibyz_f)

  call ConvolkineticT_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_fc,x_f,y_c,y_f,ekin)

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       scal,y_c,y_f,hpsi(1),hpsi(nvctr_c+1))

END SUBROUTINE applylocpotkinone_prev

subroutine realspace_prev(ibyyzz_r,pot,psir,epot,n1,n2,n3)
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
        do i1=max(ibyyzz_r(1,i2,i3),-14),min(ibyyzz_r(2,i2,i3),2*n1+16)
           tt=pot(i1,i2,i3)*psir(i1,i2,i3)
           epot=epot+tt*psir(i1,i2,i3)
           psir(i1,i2,i3)=tt
        enddo
     enddo
  enddo

end subroutine realspace_prev

subroutine realspace_nbuf_prev(ibyyzz_r,pot,psir,epot,nb1,nb2,nb3,nbuf)
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
              do i1=-14+2*nbuf,ibyyzz_r(1,i2,i3)-1
                 psir(i1,i2,i3)=0.d0
              enddo
              do i1=max(ibyyzz_r(1,i2,i3),-14+2*nbuf),min(ibyyzz_r(2,i2,i3),2*nb1+16-2*nbuf)
                 tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf)*psir(i1,i2,i3)
                 epot=epot+tt*psir(i1,i2,i3)
                 psir(i1,i2,i3)=tt
              enddo
              do i1=ibyyzz_r(2,i2,i3)+1,2*nb1+16-2*nbuf
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

end subroutine realspace_nbuf_prev


SUBROUTINE CALC_GRAD_REZA_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsi_c,xpsi_f,ypsi_c,ypsi_f)
  ! ypsi = (1/2) \Nabla^2 xpsi
  use module_base
  implicit none
  integer :: keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
  real(kind=8) :: xpsi_c(nvctr_c),xpsi_f(7,nvctr_f),scal(0:3)
  real(kind=8) :: ypsi_c(nvctr_c),ypsi_f(7,nvctr_f)
  integer :: ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  integer :: ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  character(len=*), parameter :: subname='CALC_GRAD_REZA_prev'
  real(kind=8),allocatable :: xpsig_c(:,:,:),ypsig_c(:,:,:)
  real(kind=8),allocatable :: xpsig_f(:,:,:,:),ypsig_f(:,:,:,:)
  real(kind=8),allocatable :: xpsig_fc(:,:,:,:)
  real(kind=8) :: cprecr,hgrid
  integer :: i_all,i_stat
  integer :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,nseg_c,nseg_f,nvctr_c,nvctr_f

  allocate(xpsig_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_c,'xpsig_c',subname)
  allocate(xpsig_fc(0:n1,0:n2,0:n3,3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_fc,'xpsig_fc',subname)
  allocate(xpsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_f,'xpsig_f',subname)
  allocate(ypsig_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,ypsig_c,'ypsig_c',subname)
  allocate(ypsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,ypsig_f,'ypsig_f',subname)

  call uncompress_forstandard_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,xpsi_c,xpsi_f,xpsig_c,xpsig_fc,xpsig_f)

  call Convolkinetic_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
       cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,xpsig_fc,xpsig_f,ypsig_c,ypsig_f)

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)

  i_all=-product(shape(xpsig_c))*kind(xpsig_c)
  deallocate(xpsig_c,stat=i_stat)
  call memocc(i_stat,i_all,'xpsig_c',subname)
  i_all=-product(shape(ypsig_c))*kind(ypsig_c)
  deallocate(ypsig_c,stat=i_stat)
  call memocc(i_stat,i_all,'ypsig_c',subname)
  i_all=-product(shape(xpsig_f))*kind(xpsig_f)
  deallocate(xpsig_f,stat=i_stat)
  call memocc(i_stat,i_all,'xpsig_f',subname)
  i_all=-product(shape(ypsig_f))*kind(ypsig_f)
  deallocate(ypsig_f,stat=i_stat)
  call memocc(i_stat,i_all,'ypsig_f',subname)
  i_all=-product(shape(xpsig_fc))*kind(xpsig_fc)
  deallocate(xpsig_fc,stat=i_stat)
  call memocc(i_stat,i_all,'xpsig_fc',subname)

END SUBROUTINE CALC_GRAD_REZA_prev
