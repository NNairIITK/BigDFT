!> @file
!!  Kinetic convolution routines (unused)
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine convolut_kinetic(n1,n2,n3,hgrid,x,y)
  !use module_profile
  !   applies the kinetic energy operator onto x to get y
  implicit real(kind=8) (a-h,o-z)
  parameter(lowfil=-14,lupfil=14)
  dimension fil(lowfil:lupfil),x(0:n1,0:n2,0:n3),y(0:n1,0:n2,0:n3)

  scale=-.5d0/hgrid**2

  !! second derivative filters for Daubechies 8
  !        fil(0)=-342643.d0/82248.d0*scale
  !        fil(1)=2852128.d0/1079505.d0*scale
  !        fil(2)=-12053651.d0/17272080.d0*scale
  !        fil(3)=162976.d0/1079505.d0*scale
  !        fil(4)=-60871.d0/5757360.d0*scale
  !        fil(5)=-352.d0/215901.d0*scale
  !        fil(6)=55.d0/3454416.d0*scale
  
  ! Routine works only for n > 26 because it is unrolled according to 2*(filter length)
  if (n1 .lt. 27 .or. n2 .lt. 27 .or. n3 .lt. 27) then
    stop 'In convolut_kinetic optimized routines it is required that: n1 > 26 ,n2 > 26 ,n3 > 26. Please run convolut_simple routines.'
  endif

  ! second derivative filters for Daubechies 16
  fil(-14)=  -6.924474940639200152025730585882d-18*scale
  fil(-13)=   2.70800493626319438269856689037647576d-13*scale
  fil(-12)=  -5.813879830282540547959250667d-11*scale
  fil(-11)=  -1.05857055496741470373494132287d-8*scale
  fil(-10)=  -3.7230763047369275848791496973044d-7*scale
  fil(-9)=    2.0904234952920365957922889447361d-6*scale
  fil(-8)=   -0.2398228524507599670405555359023135d-4*scale
  fil(-7)=    0.45167920287502235349480037639758496d-3*scale
  fil(-6)=   -0.409765689342633823899327051188315485d-2*scale
  fil(-5)=    0.02207029188482255523789911295638968409d0*scale
  fil(-4)=   -0.0822663999742123340987663521d0*scale
  fil(-3)=    0.2371780582153805636239247476d0*scale
  fil(-2)=   -0.6156141465570069496314853949d0*scale
  fil(-1)=    2.2191465938911163898794546405d0*scale
  fil(0)=    -3.5536922899131901941296809374d0*scale
  fil(1)=     2.2191465938911163898794546405d0*scale
  fil(2)=    -0.6156141465570069496314853949d0*scale
  fil(3)=     0.2371780582153805636239247476d0*scale
  fil(4)=    -0.0822663999742123340987663521d0*scale
  fil(5)=     0.02207029188482255523789911295638968409d0*scale
  fil(6)=    -0.409765689342633823899327051188315485d-2*scale
  fil(7)=     0.45167920287502235349480037639758496d-3*scale
  fil(8)=    -0.2398228524507599670405555359023135d-4*scale
  fil(9)=     2.0904234952920365957922889447361d-6*scale
  fil(10)=   -3.7230763047369275848791496973044d-7*scale
  fil(11)=   -1.05857055496741470373494132287d-8*scale
  fil(12)=   -5.813879830282540547959250667d-11*scale
  fil(13)=    2.70800493626319438269856689037647576d-13*scale
  fil(14)=   -6.924474940639200152025730585882d-18*scale

  iend1 = n1 - lupfil
  iend2 = n2 - lupfil
  iend3 = n3 - lupfil

!$omp  parallel do default(private) shared(x,y,fil,n3,n2,n1,iend1,iend2,iend3)
  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2

        do i1=0,lupfil-4,4
           tt0=0.d0
           tt1=x(0,i2,i3)*fil(-i1-1)
           tt2=x(0,i2,i3)*fil(-i1-2)
           tt3=x(0,i2,i3)*fil(-i1-3)
           tt2=tt2+x(1,i2,i3)*fil(-i1-1)
           tt3=tt3+x(1,i2,i3)*fil(-i1-2)
           tt3=tt3+x(2,i2,i3)*fil(-i1-1)
           do l=-i1,lupfil
              tt0=tt0+x(i1+l,i2,i3)*fil(l)
              tt1=tt1+x(i1+1+l,i2,i3)*fil(l)
              tt2=tt2+x(i1+2+l,i2,i3)*fil(l)
              tt3=tt3+x(i1+3+l,i2,i3)*fil(l)
           enddo
           y(i1,i2,i3)=tt0
           y(i1+1,i2,i3)=tt1
           y(i1+2,i2,i3)=tt2
           y(i1+3,i2,i3)=tt3
        enddo
        do i1=i1,lupfil-1
           tt0=0.d0
           do l=-i1,lupfil
              tt0=tt0+x(i1+l,i2,i3)*fil(l)
           enddo
           y(i1,i2,i3)=tt0
        enddo

        do i1=i1,iend1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=lowfil,lupfil
              tt0=tt0+x(i1+l,i2,i3)*fil(l)
              tt1=tt1+x(i1+1+l,i2,i3)*fil(l)
              tt2=tt2+x(i1+2+l,i2,i3)*fil(l)
              tt3=tt3+x(i1+3+l,i2,i3)*fil(l)
           enddo
           y(i1,i2,i3)=tt0
           y(i1+1,i2,i3)=tt1
           y(i1+2,i2,i3)=tt2
           y(i1+3,i2,i3)=tt3
        enddo
        do i1=i1,iend1
           tt0=0.d0
           do l=lowfil,lupfil
              tt0=tt0+x(i1+l,i2,i3)*fil(l)
           enddo
           y(i1,i2,i3)=tt0
        enddo

        do i1=i1,n1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=lowfil,n1-i1-3
              tt0=tt0+x(i1+l,i2,i3)*fil(l)
              tt1=tt1+x(i1+1+l,i2,i3)*fil(l)
              tt2=tt2+x(i1+2+l,i2,i3)*fil(l)
              tt3=tt3+x(i1+3+l,i2,i3)*fil(l)
           enddo
           tt2=tt2+x(n1,i2,i3)*fil(n1-i1-2)
           tt1=tt1+x(n1,i2,i3)*fil(n1-i1-1)
           tt0=tt0+x(n1,i2,i3)*fil(n1-i1)
           tt1=tt1+x(n1-1,i2,i3)*fil(n1-i1-2)
           tt0=tt0+x(n1-1,i2,i3)*fil(n1-i1-1)
           tt0=tt0+x(n1-2,i2,i3)*fil(n1-i1-2)
           y(i1,i2,i3)=tt0
           y(i1+1,i2,i3)=tt1
           y(i1+2,i2,i3)=tt2
           y(i1+3,i2,i3)=tt3
        enddo
        do i1=i1,n1
           tt0=0.d0
           do l=lowfil,n1-i1
              tt0=tt0+x(i1+l,i2,i3)*fil(l)
           enddo
           y(i1,i2,i3)=tt0
        enddo

     enddo


     ! + (1/2) d^2/dy^2
     do i2=0,lupfil-1
        do i1=0,n1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=-i2,lupfil
              tt0=tt0+x(i1,i2+l,i3)*fil(l)
              tt1=tt1+x(i1+1,i2+l,i3)*fil(l)
              tt2=tt2+x(i1+2,i2+l,i3)*fil(l)
              tt3=tt3+x(i1+3,i2+l,i3)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
           y(i1+1,i2,i3)=y(i1+1,i2,i3)+tt1
           y(i1+2,i2,i3)=y(i1+2,i2,i3)+tt2
           y(i1+3,i2,i3)=y(i1+3,i2,i3)+tt3
        enddo
        do i1=i1,n1
           tt0=0.d0
           do l=-i2,lupfil
              tt0=tt0+x(i1,i2+l,i3)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
        enddo
     enddo

     do i2=i2,iend2
        do i1=0,n1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=lowfil,lupfil
              tt0=tt0+x(i1,i2+l,i3)*fil(l)
              tt1=tt1+x(i1+1,i2+l,i3)*fil(l)
              tt2=tt2+x(i1+2,i2+l,i3)*fil(l)
              tt3=tt3+x(i1+3,i2+l,i3)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
           y(i1+1,i2,i3)=y(i1+1,i2,i3)+tt1
           y(i1+2,i2,i3)=y(i1+2,i2,i3)+tt2
           y(i1+3,i2,i3)=y(i1+3,i2,i3)+tt3
        enddo
        do i1=i1,n1
           tt0=0.d0
           do l=lowfil,lupfil
              tt0=tt0+x(i1,i2+l,i3)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
        enddo
     enddo

     do i2=i2,n2
        do i1=0,n1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=lowfil,n2-i2
              tt0=tt0+x(i1,i2+l,i3)*fil(l)
              tt1=tt1+x(i1+1,i2+l,i3)*fil(l)
              tt2=tt2+x(i1+2,i2+l,i3)*fil(l)
              tt3=tt3+x(i1+3,i2+l,i3)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
           y(i1+1,i2,i3)=y(i1+1,i2,i3)+tt1
           y(i1+2,i2,i3)=y(i1+2,i2,i3)+tt2
           y(i1+3,i2,i3)=y(i1+3,i2,i3)+tt3
        enddo
        do i1=i1,n1
           tt0=0.d0
           do l=lowfil,n2-i2
              tt0=tt0+x(i1,i2+l,i3)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
        enddo
     enddo
  enddo
!$omp  end parallel do

!$omp  parallel do default(private) shared(x,y,fil,n3,n2,n1,iend1,iend2,iend3)
 ! (1/2) d^2/dz^2
  do i2=0,n2
     do i3=0,lupfil-1
        do i1=0,n1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=-i3,lupfil
              tt0=tt0+x(i1,i2,i3+l)*fil(l)
              tt1=tt1+x(i1+1,i2,i3+l)*fil(l)
              tt2=tt2+x(i1+2,i2,i3+l)*fil(l)
              tt3=tt3+x(i1+3,i2,i3+l)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
           y(i1+1,i2,i3)=y(i1+1,i2,i3)+tt1
           y(i1+2,i2,i3)=y(i1+2,i2,i3)+tt2
           y(i1+3,i2,i3)=y(i1+3,i2,i3)+tt3
        enddo
        do i1=i1,n1
           tt0=0.d0
           do l=-i3,lupfil
              tt0=tt0+x(i1,i2,i3+l)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
        enddo
     enddo

     do i3=i3,iend3
        do i1=0,n1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=lowfil,lupfil
              tt0=tt0+x(i1,i2,i3+l)*fil(l)
              tt1=tt1+x(i1+1,i2,i3+l)*fil(l)
              tt2=tt2+x(i1+2,i2,i3+l)*fil(l)
              tt3=tt3+x(i1+3,i2,i3+l)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
           y(i1+1,i2,i3)=y(i1+1,i2,i3)+tt1
           y(i1+2,i2,i3)=y(i1+2,i2,i3)+tt2
           y(i1+3,i2,i3)=y(i1+3,i2,i3)+tt3
        enddo
        do i1=i1,n1
           tt0=0.d0
           do l=lowfil,lupfil
              tt0=tt0+x(i1,i2,i3+l)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
        enddo
     enddo

     do i3=i3,n3
        do i1=0,n1-3,4
           tt0=0.d0
           tt1=0.d0
           tt2=0.d0
           tt3=0.d0
           do l=lowfil,n3-i3
              tt0=tt0+x(i1,i2,i3+l)*fil(l)
              tt1=tt1+x(i1+1,i2,i3+l)*fil(l)
              tt2=tt2+x(i1+2,i2,i3+l)*fil(l)
              tt3=tt3+x(i1+3,i2,i3+l)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
           y(i1+1,i2,i3)=y(i1+1,i2,i3)+tt1
           y(i1+2,i2,i3)=y(i1+2,i2,i3)+tt2
           y(i1+3,i2,i3)=y(i1+3,i2,i3)+tt3
        enddo
        do i1=i1,n1
           tt0=0.d0
           do l=lowfil,n3-i3
              tt0=tt0+x(i1,i2,i3+l)*fil(l)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt0
        enddo
     enddo
  enddo
!$omp  end parallel do

  return


END SUBROUTINE convolut_kinetic


subroutine convrot_grow(n1,ndat,x,y)
  implicit real(kind=8) (a-h,o-z)
  parameter(lowfil=-8,lupfil=7)
  dimension x(0:n1,ndat),y(ndat,-lupfil:n1-lowfil)    

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
       
  ! Routine works only for n1 > 13 because it is unrolled according to 2*(filter length)
  if (n1 .lt. 14) then
    stop 'In convrot_grow optimized routines it is required that n1 > 13. Please run convolut_simple routines.'
  endif

!$omp  parallel do default(private) shared(x,y,fil,ndat,n1)
  do j=1,ndat

     do i=-lupfil,lupfil-4,4
        tt0=0.d0
        tt1=x(0,j)*fil(-i-1)
        tt2=x(0,j)*fil(-i-2)
        tt3=x(0,j)*fil(-i-3)
        tt2=tt2+x(1,j)*fil(-i-1)
        tt3=tt3+x(1,j)*fil(-i-2)
        tt3=tt3+x(2,j)*fil(-i-1)
        do l=-i,lupfil
           tt0=tt0+x(i+l,j)*fil(l)
           tt1=tt1+x(i+1+l,j)*fil(l)
           tt2=tt2+x(i+2+l,j)*fil(l)
           tt3=tt3+x(i+3+l,j)*fil(l)
        enddo
        y(j,i)=tt0
        y(j,i+1)=tt1
        y(j,i+2)=tt2
        y(j,i+3)=tt3
     enddo
     do i=i,lupfil
        tt0=0.d0
        do l=-i,lupfil
           tt0=tt0+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt0
     enddo

     do i=i,n1-3-lupfil,4
        tt0=0.d0
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        do l=lowfil,lupfil
           tt0=tt0+x(i+l,j)*fil(l)
           tt1=tt1+x(i+1+l,j)*fil(l)
           tt2=tt2+x(i+2+l,j)*fil(l)
           tt3=tt3+x(i+3+l,j)*fil(l)
        enddo
        y(j,i)=tt0
        y(j,i+1)=tt1
        y(j,i+2)=tt2
        y(j,i+3)=tt3
     enddo
     do i=i,n1-lupfil
        tt0=0.d0
        do l=lowfil,lupfil
           tt0=tt0+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt0
     enddo

     do i=i,n1-lowfil-3,4
        tt0=0.d0
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        do l=lowfil,n1-i-3
           tt0=tt0+x(i+l,j)*fil(l)
           tt1=tt1+x(i+1+l,j)*fil(l)
           tt2=tt2+x(i+2+l,j)*fil(l)
           tt3=tt3+x(i+3+l,j)*fil(l)
        enddo
        tt2=tt2+x(n1,j)*fil(n1-i-2)
        tt1=tt1+x(n1,j)*fil(n1-i-1)
        tt0=tt0+x(n1,j)*fil(n1-i)
        tt1=tt1+x(n1-1,j)*fil(n1-i-2)
        tt0=tt0+x(n1-1,j)*fil(n1-i-1)
        tt0=tt0+x(n1-2,j)*fil(n1-i-2)
        y(j,i)=tt0
        y(j,i+1)=tt1
        y(j,i+2)=tt2
        y(j,i+3)=tt3
     enddo
     do i=i,n1-lowfil
        tt0=0.d0
        do l=lowfil,n1-i
           tt0=tt0+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt0
     enddo
  enddo
!$omp end parallel do
  return
END SUBROUTINE convrot_grow


subroutine convrot_shrink(n1,ndat,x,y)
  implicit real(kind=8) (a-h,o-z)
  parameter(lowfil=-8,lupfil=7)
  dimension x(lowfil:n1+lupfil,ndat),y(ndat,0:n1)

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
       
  ! Routine works only for n1 > 2 because it is unrolled by 4
  if (n1 .lt. 3) then
    stop 'In convrot_shrink optimized routines it is required that n1 > 2. Please run convolut_simple routines.'
  endif

!$omp  parallel do default(private) shared(x,y,fil,n1,ndat)
  do j=1,ndat
     do i=0,n1-3,4
        tt0=0.d0
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        do l=lowfil,lupfil
           tt0=tt0+x(i+l,j)*fil(l)
           tt1=tt1+x(i+1+l,j)*fil(l)
           tt2=tt2+x(i+2+l,j)*fil(l)
           tt3=tt3+x(i+3+l,j)*fil(l)
        enddo
        y(j,i)=tt0
        y(j,i+1)=tt1
        y(j,i+2)=tt2
        y(j,i+3)=tt3
     enddo
     do i=i,n1
        tt0=0.d0
        do l=lowfil,lupfil
           tt0=tt0+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt0
     enddo
  enddo
!$omp end parallel do
  
  return
END SUBROUTINE convrot_shrink


subroutine syn_rot_grow(n,nt,x,y)
  implicit real(kind=8) (a-h,o-z)
  !  Daubechies_16 symetric 
  !        real(kind=8) CH(-M:M)/0.d0,&
  !        -0.0033824159510050025955D0,-0.00054213233180001068935D0,&
  !        0.031695087811525991431D0,0.0076074873249766081919D0,&
  !        -0.14329423835127266284D0,-0.061273359067811077843D0,&
  !        0.48135965125905339159D0,0.77718575169962802862D0,0.36444189483617893676D0,&
  !        -0.051945838107881800736D0,-0.027219029917103486322D0,&
  !        0.049137179673730286787D0,0.0038087520138944894631D0,&
  !        -0.014952258337062199118D0,-0.00030292051472413308126D0,&
  !        0.0018899503327676891843D0&
  !       daubechy S6
  parameter(m=8,mm=m+2)
  
  dimension x(0:2*n+1,nt),y(nt,-m+1:2*n+m)
          
    !******** coefficints for wavelet transform *********************
    ! h coefficients
    real(kind=8) ch(-mm:mm)
    DATA ch / &
    0.D0, & 
    0.D0, &
    0.d0, &
    -0.0033824159510050025955D0, &
    -0.00054213233180001068935D0, &
    0.031695087811525991431D0, &
    0.0076074873249766081919D0, &
    -0.14329423835127266284D0, &
    -0.061273359067811077843D0, &
    0.48135965125905339159D0, &
    0.77718575169962802862D0, &
    0.36444189483617893676D0, &
    -0.051945838107881800736D0, &
    -0.027219029917103486322D0, &
    0.049137179673730286787D0, &
    0.0038087520138944894631D0, &
    -0.014952258337062199118D0, &
    -0.00030292051472413308126D0, &
    0.0018899503327676891843D0, &
    0.D0, &
    0.D0 /

! g coefficients
    real(kind=8) cg(-mm:mm)
    DATA cg / &
    0.D0, & 
    0.D0, &
    0.d0, &
    -0.0018899503327676891843D0, &
    -0.00030292051472413308126D0, &
    0.014952258337062199118D0, &
    0.0038087520138944894631D0, &
    -0.049137179673730286787D0, &
    -0.027219029917103486322D0, &
    0.051945838107881800736D0, &
    0.36444189483617893676D0, &
    -0.77718575169962802862D0, &
    0.48135965125905339159D0, &
    0.061273359067811077843D0, &
    -0.14329423835127266284D0, &
    -0.0076074873249766081919D0, &
    0.031695087811525991431D0, &
    0.00054213233180001068935D0, &
    -0.0033824159510050025955D0, &
    0.D0, &
    0.D0 /      
    
  ! Routine works only for n > 5 because it is unrolled according to 2*(filter length)
  if (n .lt. 6) then
    stop 'In syn_rot_grow optimized routines it is required that n > 5. Please run convolut_simple routines.'
  endif


!$omp  parallel do default(private) shared(x,y,cg,ch,n,nt)
  do l=1,nt

     do i=-M+1,M-5,4
        tt0 = 0.d0
        tt1 = 0.d0  
        tt2 = ch(i+2)*x(0,l) + cg(i+2)*x(n+1,l)
        tt3 = ch(i+3)*x(0,l) + cg(i+3)*x(n+1,l)
        do j=-M+1,i,2
           tt0 = tt0 + ch(j)*x((i-j)/2,l) + cg(j)*x(n+1+(i-j)/2,l)
           tt1 = tt1 + ch(j+1)*x((i-j)/2,l) + cg(j+1)*x(n+1+(i-j)/2,l)  
           tt2 = tt2 + ch(j)*x((i-j)/2+1,l) + cg(j)*x(n+1+(i-j)/2+1,l)
           tt3 = tt3 + ch(j+1)*x((i-j)/2+1,l) + cg(j+1)*x(n+1+(i-j)/2+1,l)
        enddo
        y(l,i) = tt0
        y(l,i+1) = tt1
        y(l,i+2) = tt2
        y(l,i+3) = tt3
     enddo
     do i=i,M-3,2
        tt0 = 0.d0
        tt1 = 0.d0          
        do j=-M+1,i,2
           tt0 = tt0 + ch(j)*x((i-j)/2,l) + cg(j)*x(n+1+(i-j)/2,l)
           tt1 = tt1 + ch(j+1)*x((i-j)/2,l) + cg(j+1)*x(n+1+(i-j)/2,l)             
        enddo
        y(l,i) = tt0
        y(l,i+1) = tt1        
     enddo

     do i=i,2*n-M-1,4
        tt0 = 0.d0
        tt1 = 0.d0
        tt2 = 0.d0
        tt3 = 0.d0
        do j=-M+1,M-1,2
           tt0 = tt0 + ch(j)*x((i-j)/2,l) + cg(j)*x(n+1+(i-j)/2,l)
           tt1 = tt1 + ch(j+1)*x((i-j)/2,l) + cg(j+1)*x(n+1+(i-j)/2,l)
           tt2 = tt2 + ch(j)*x((i-j)/2+1,l) + cg(j)*x(n+1+(i-j)/2+1,l)
           tt3 = tt3 + ch(j+1)*x((i-j)/2+1,l) + cg(j+1)*x(n+1+(i-j)/2+1,l)
        enddo
        y(l,i) = tt0
        y(l,i+1) = tt1
        y(l,i+2) = tt2
        y(l,i+3) = tt3
     enddo
     do i=i,2*n-M+1,2
        tt0 = 0.d0
        tt1 = 0.d0       
        do j=-M+1,M-1,2
           tt0 = tt0 + ch(j)*x((i-j)/2,l) + cg(j)*x(n+1+(i-j)/2,l)
           tt1 = tt1 + ch(j+1)*x((i-j)/2,l) + cg(j+1)*x(n+1+(i-j)/2,l)           
        enddo
        y(l,i) = tt0
        y(l,i+1) = tt1        
     enddo

     do i=i,2*n+6,4
        tt0 = ch(i-2*n)*x(n,l) + cg(i-2*n)*x(n+1+n,l)
        tt1 = ch(i-2*n+1)*x(n,l) + cg(i-2*n+1)*x(n+1+n,l)
        tt2 = 0.d0
        tt3 = 0.d0
        do j=i-2*n+2,M-1,2
           tt0 = tt0 + ch(j)*x((i-j)/2,l) + cg(j)*x(n+1+(i-j)/2,l)
           tt1 = tt1 + ch(j+1)*x((i-j)/2,l) + cg(j+1)*x(n+1+(i-j)/2,l)
           tt2 = tt2 + ch(j)*x((i-j)/2+1,l) + cg(j)*x(n+1+(i-j)/2+1,l)
           tt3 = tt3 + ch(j+1)*x((i-j)/2+1,l) + cg(j+1)*x(n+1+(i-j)/2+1,l)
        enddo
        y(l,i) = tt0
        y(l,i+1) = tt1
        y(l,i+2) = tt2
        y(l,i+3) = tt3
     enddo
     do i=i,2*n+8,2
        tt0 = 0.d0
        tt1 = 0.d0       
        do j=i-2*n,M-1,2
           tt0 = tt0 + ch(j)*x((i-j)/2,l) + cg(j)*x(n+1+(i-j)/2,l)
           tt1 = tt1 + ch(j+1)*x((i-j)/2,l) + cg(j+1)*x(n+1+(i-j)/2,l)           
        enddo
        y(l,i) = tt0
        y(l,i+1) = tt1        
     enddo

  enddo
!$omp end parallel do
           
  return
  
END SUBROUTINE syn_rot_grow



subroutine ana_rot_shrink(n,nt,x,y)
  implicit real(kind=8) (a-h,o-z)
  !  Daubechies_8 symetric 
  !        parameter(chm3=.230377813308896501d0, chm2=.714846570552915647d0, & 
  !                  chm1=.630880767929858908d0, ch00=-.279837694168598542d-1, & 
  !                  chp1=-.187034811719093084d0, chp2=.308413818355607636d-1, & 
  !                  chp3=.328830116668851997d-1, chp4=-.105974017850690321d-1)
  !       daubechy S6
  
  parameter(m=8,mm=m+2)

  dimension x(-m+1:2*n+m,nt),y(nt,0:2*n+1)
  
        
    !******** coefficients for wavelet transform *********************
    ! h coefficients
    real(kind=8) ch(-mm:mm)
    DATA ch / &
    0.D0, & 
    0.D0, &
    0.d0, &
    -0.0033824159510050025955D0, &
    -0.00054213233180001068935D0, &
    0.031695087811525991431D0, &
    0.0076074873249766081919D0, &
    -0.14329423835127266284D0, &
    -0.061273359067811077843D0, &
    0.48135965125905339159D0, &
    0.77718575169962802862D0, &
    0.36444189483617893676D0, &
    -0.051945838107881800736D0, &
    -0.027219029917103486322D0, &
    0.049137179673730286787D0, &
    0.0038087520138944894631D0, &
    -0.014952258337062199118D0, &
    -0.00030292051472413308126D0, &
    0.0018899503327676891843D0, &
    0.D0, &
    0.D0 /

! g coefficients
    real(kind=8) cg(-mm:mm)
    DATA cg / &
    0.D0, & 
    0.D0, &
    0.d0, &
    -0.0018899503327676891843D0, &
    -0.00030292051472413308126D0, &
    0.014952258337062199118D0, &
    0.0038087520138944894631D0, &
    -0.049137179673730286787D0, &
    -0.027219029917103486322D0, &
    0.051945838107881800736D0, &
    0.36444189483617893676D0, &
    -0.77718575169962802862D0, &
    0.48135965125905339159D0, &
    0.061273359067811077843D0, &
    -0.14329423835127266284D0, &
    -0.0076074873249766081919D0, &
    0.031695087811525991431D0, &
    0.00054213233180001068935D0, &
    -0.0033824159510050025955D0, &
    0.D0, &
    0.D0 /      
    
  ! Routine works only for n > 0 because it is unrolled by 2
  if (n .lt. 1) then
    stop 'In syn_rot_grow optimized routines it is required that n > 0. Please run convolut_simple routines.'
  endif

!$omp  parallel do default(private) shared(x,y,ch,cg,n,nt)
  do l=1,nt

     do i=0,n-1,2
        i2=2*i
        ci0=0.D0
        di0=0.D0
        ci1=0.D0
        di1=0.D0
        do j=-M+1,M
           ci0=ci0+ch(j)*x(j+i2,l)
           ci1=ci1+ch(j)*x(j+i2+2,l)
           di0=di0+cg(j)*x(j+i2,l)
           di1=di1+cg(j)*x(j+i2+2,l)
        enddo
        y(l,i)=ci0
        y(l,n+1+i)=di0
        y(l,i+1)=ci1
        y(l,n+2+i)=di1
     enddo
     do i=i,n
        i2=2*i
        ci0=0.D0
        di0=0.D0        
        do j=-M+1,M
           ci0=ci0+ch(j)*x(j+i2,l)
           di0=di0+cg(j)*x(j+i2,l)           
        enddo
        y(l,i)=ci0
        y(l,n+1+i)=di0       
     enddo

  enddo
!$omp end parallel do

  return
END SUBROUTINE ana_rot_shrink

