!> @file
!!  Common convolutions
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> A analysis wavelet transformation where the size of the data is forced to shrink
!! The input array y is overwritten
subroutine analyse_shrink(n1,n2,n3,ww,y,x)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8), intent(inout) :: ww
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(inout) :: y
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: x
  !local variables
  integer :: nt

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+16)*(2*n3+16)
  call  ana_rot_shrink(n1,nt,y,ww)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+16)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_shrink(n3,nt,y,x)

END SUBROUTINE analyse_shrink


!> A synthesis wavelet transformation where the size of the data is allowed to grow
!! The input array x is not overwritten
subroutine synthese_grow(n1,n2,n3,ww,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(in) :: x
  real(wp), dimension(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8), intent(inout) :: ww
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(inout) :: y
  !local variables
  integer :: nt

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
  use module_base
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(-7:2*n+8,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: ci,di
  real(wp), dimension(-7:8) :: ch,cg
  !       Daubechy S16
  data ch  /  -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp /
  data cg  / -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp  /

  do j=1,ndat
     do i=0,n
        ci=0.e0_wp
        di=0.e0_wp
        do l=-7,8
           ci=ci+ch(l)*x(l+2*i,j)
           di=di+cg(l)*x(l+2*i,j)
        enddo
        y(j,i)=ci
        y(j,n+1+i)=di
     enddo
  enddo

END SUBROUTINE ana_rot_shrink


subroutine syn_rot_grow(n,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-7:2*n+8), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: so,se
  real(wp), dimension(-8:9) :: ch,cg
  !       Daubechy S16
  data ch  /  0.e0_wp , -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp , 0.e0_wp /
  data cg  / 0.e0_wp , -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp , 0.e0_wp /

  do j=1,ndat

     i=-4
     so=0.e0_wp
     do l=max(i-n,-4),min(i,4)
        so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+1+i-l,j)
     enddo
     y(j,2*i+1)=so

     do i=-3,n+3
        se=0.e0_wp
        so=0.e0_wp
        do l=max(i-n,-4),min(i,4)
           se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+1+i-l,j)
           so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+1+i-l,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

     i=n+4
     se=0.e0_wp
     do l=max(i-n,-4),min(i,4)
        se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+1+i-l,j)
     enddo
     y(j,2*i  )=se

  enddo

END SUBROUTINE syn_rot_grow


!> Simple non-optimized version of the major convolution routines
subroutine convrot_grow(n1,ndat,x,y)
  use module_base
  implicit none
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-lupfil:n1-lowfil), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /

  do j=1,ndat
     do i=-lupfil,n1-lowfil
        tt=0.e0_wp
        do l=max(-i,lowfil),min(lupfil,n1-i)
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE convrot_grow


!> Simple non-optimized version of the major convolution routines
subroutine convrot_shrink(n1,ndat,x,y)
  use module_base
  implicit none
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(lowfil:n1+lupfil,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228e-6_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.2103025160930381434955489412839065067e-1_wp,&
       -0.604895289196983516002834636e-1_wp,&
       0.9940415697834003993178616713_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       8.4334247333529341094733325815816e-7_wp /

  do j=1,ndat
     do i=0,n1
        tt=0.e0_wp
        do l=lowfil,lupfil
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE convrot_shrink
