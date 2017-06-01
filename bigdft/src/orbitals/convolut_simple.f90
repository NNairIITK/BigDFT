!> @file
!!  Simple convolution routines
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> A analysis wavelet transformation where the size of the data is forced to shrink
!! The input array y is overwritten
subroutine analyse_shrink(n1,n2,n3,ww,y,x)
  implicit real(kind=8) (a-h,o-z)
  dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
  dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
  dimension x(0:n1,2,0:n2,2,0:n3,2)

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+16)*(2*n3+16)
  call  ana_rot_shrink(n1,nt,y,ww)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+16)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_shrink(n3,nt,y,x)

  return
END SUBROUTINE analyse_shrink


!> A synthesis wavelet transformation where the size of the data is allowed to grow
!! The input array x is not overwritten
subroutine synthese_grow(n1,n2,n3,ww,x,y)
  implicit real(kind=8) (a-h,o-z)
  dimension x(0:n1,2,0:n2,2,0:n3,2)
  dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
  dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

  ! i1,i2,i3 -> i2,i3,I1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,ww)
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_grow(n3,nt,ww,y)

END SUBROUTINE synthese_grow


subroutine ana_rot_shrink(n,ndat,x,y)
  implicit real(kind=8) (a-h,o-z)
  dimension x(-7:2*n+8,ndat),y(ndat,0:2*n+1)
  real(kind=8) ch(-7:8) ,cg(-7:8)
  !       Daubechy S16
  data ch  /  -0.0033824159510050025955D0, & 
       -0.00054213233180001068935D0, 0.031695087811525991431D0, & 
       0.0076074873249766081919D0, -0.14329423835127266284D0, & 
       -0.061273359067811077843D0, 0.48135965125905339159D0,  & 
       0.77718575169962802862D0,0.36444189483617893676D0, &
       -0.051945838107881800736D0,-0.027219029917103486322D0, &
       0.049137179673730286787D0,0.0038087520138944894631D0, &
       -0.014952258337062199118D0,-0.00030292051472413308126D0, &
       0.0018899503327676891843D0 /
  data cg  / -0.0018899503327676891843D0, &
       -0.00030292051472413308126D0, 0.014952258337062199118D0, &
       0.0038087520138944894631D0, -0.049137179673730286787D0, &
       -0.027219029917103486322D0, 0.051945838107881800736D0, &
       0.36444189483617893676D0, -0.77718575169962802862D0, &
       0.48135965125905339159D0, 0.061273359067811077843D0, &
       -0.14329423835127266284D0, -0.0076074873249766081919D0, &
       0.031695087811525991431D0, 0.00054213233180001068935D0, &
       -0.0033824159510050025955D0  /

  do j=1,ndat

     do i=0,n
        ci=0.d0
        di=0.d0
        do l=-7,8
           ci=ci+ch(l)*x(l+2*i,j)
           di=di+cg(l)*x(l+2*i,j)
        enddo
        y(j,i)=ci
        y(j,n+1+i)=di
     enddo

  enddo

  return
END SUBROUTINE ana_rot_shrink


subroutine syn_rot_grow(n,ndat,x,y)
  implicit real(kind=8) (a-h,o-z)
  dimension x(0:2*n+1,ndat),y(ndat,-7:2*n+8)
  real(kind=8) ch(-8:9) ,cg(-8:9)
  !       Daubechy S16
  data ch  /  0.d0 , -0.0033824159510050025955D0, & 
       -0.00054213233180001068935D0, 0.031695087811525991431D0, & 
       0.0076074873249766081919D0, -0.14329423835127266284D0, & 
       -0.061273359067811077843D0, 0.48135965125905339159D0,  & 
       0.77718575169962802862D0,0.36444189483617893676D0, &
       -0.051945838107881800736D0,-0.027219029917103486322D0, &
       0.049137179673730286787D0,0.0038087520138944894631D0, &
       -0.014952258337062199118D0,-0.00030292051472413308126D0, &
       0.0018899503327676891843D0 , 0.d0 /
  data cg  / 0.d0 , -0.0018899503327676891843D0, &
       -0.00030292051472413308126D0, 0.014952258337062199118D0, &
       0.0038087520138944894631D0, -0.049137179673730286787D0, &
       -0.027219029917103486322D0, 0.051945838107881800736D0, &
       0.36444189483617893676D0, -0.77718575169962802862D0, &
       0.48135965125905339159D0, 0.061273359067811077843D0, &
       -0.14329423835127266284D0, -0.0076074873249766081919D0, &
       0.031695087811525991431D0, 0.00054213233180001068935D0, &
       -0.0033824159510050025955D0 , 0.d0 /

  do j=1,ndat

     i=-4
     so=0.d0
     do l=max(i-n,-4),min(i,4)
        so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+1+i-l,j)
     enddo
     y(j,2*i+1)=so

     do i=-3,n+3
        se=0.d0
        so=0.d0
        do l=max(i-n,-4),min(i,4)
           se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+1+i-l,j)
           so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+1+i-l,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

     i=n+4
     se=0.d0
     do l=max(i-n,-4),min(i,4)
        se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+1+i-l,j)
     enddo
     y(j,2*i  )=se

  enddo

  return
END SUBROUTINE syn_rot_grow


!> Simple non-optimized version of the major convolution routines
subroutine convrot_grow(n1,ndat,x,y)
  implicit real(kind=8) (a-h,o-z)
  parameter(lowfil=-8,lupfil=7)
  dimension x(0:n1,ndat),y(ndat,-lupfil:n1-lowfil)
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(kind=8) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816D-7,&
       -0.1290557201342060969516786758559028D-4,&
       0.8762984476210559564689161894116397D-4,&
       -0.30158038132690463167163703826169879D-3,&
       0.174723713672993903449447812749852942D-2,&
       -0.942047030201080385922711540948195075D-2,&
       0.2373821463724942397566389712597274535D-1,&
       0.612625895831207982195380597D-1,&
       0.9940415697834003993178616713D0,&
       -0.604895289196983516002834636D-1, &
       -0.2103025160930381434955489412839065067D-1,&
       0.1337263414854794752733423467013220997D-1,&
       -0.344128144493493857280881509686821861D-2,&
       0.49443227688689919192282259476750972D-3,&
       -0.5185986881173432922848639136911487D-4,&
       2.72734492911979659657715313017228D-6 /



  do j=1,ndat
     do i=-lupfil,n1-lowfil

        tt=0.d0
        do l=max(-i,lowfil),min(lupfil,n1-i)
           tt=tt+x(i+l,j)*fil(l)
        enddo
 y(j,i)=tt

     enddo
  enddo

  return
END SUBROUTINE convrot_grow


!> Simple non-optimized version of the major convolution routines
subroutine convrot_shrink(n1,ndat,x,y)
  implicit real(kind=8) (a-h,o-z)
  parameter(lowfil=-7,lupfil=8)
  dimension x(lowfil:n1+lupfil,ndat),y(ndat,0:n1)
  ! the filtered output data structure has shrunk by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(kind=8) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228D-6,&
       -0.5185986881173432922848639136911487D-4,&
       0.49443227688689919192282259476750972D-3,&
       -0.344128144493493857280881509686821861D-2,&
       0.1337263414854794752733423467013220997D-1,&
       -0.2103025160930381434955489412839065067D-1,&
       -0.604895289196983516002834636D-1,&
       0.9940415697834003993178616713D0,&
       0.612625895831207982195380597D-1,&
       0.2373821463724942397566389712597274535D-1,&
       -0.942047030201080385922711540948195075D-2,&
       0.174723713672993903449447812749852942D-2,&
       -0.30158038132690463167163703826169879D-3,&
       0.8762984476210559564689161894116397D-4,&
       -0.1290557201342060969516786758559028D-4,&
       8.4334247333529341094733325815816D-7 /


  do j=1,ndat
     do i=0,n1

        tt=0.d0
        do l=lowfil,lupfil
           tt=tt+x(i+l,j)*fil(l)
        enddo
 y(j,i)=tt

     enddo
  enddo

  return
END SUBROUTINE convrot_shrink
