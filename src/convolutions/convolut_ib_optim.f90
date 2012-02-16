!> @file
!! Optimized convolution routines
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> One of the most CPU intensive routines
!!   y = (kinetic energy operator)x + (cprec*I)x 
subroutine Convolkinetic(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3)
  use module_base
  implicit none
!dee
  !integer :: iend_test,count_rate_test,count_max_test,istart_test
  
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(wp), intent(in) :: cprecr
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: x_f3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer(8) :: clock0,clock1,clock2,ncount_max,ncount_rate

!$ integer(8) :: ncount1,ncount2,ncount3,ncount4,ncount5,ncount6,ncount0 !< variables used by OpenMP

  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  !real(kind=8) :: tel
  real(wp), dimension(-3+lowfil:lupfil+3) :: a,b,c
  real(wp), dimension(lowfil:lupfil) :: e

  scale=-.5_wp/real(hgrid**2,wp)
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0)=   -3.5536922899131901941296809374_wp*scale
  a(1)=    2.2191465938911163898794546405_wp*scale
  a(2)=   -0.6156141465570069496314853949_wp*scale
  a(3)=    0.2371780582153805636239247476_wp*scale
  a(4)=   -0.0822663999742123340987663521_wp*scale
  a(5)=    0.02207029188482255523789911295638968409_wp*scale
  a(6)=   -0.409765689342633823899327051188315485e-2_wp*scale
  a(7)=    0.45167920287502235349480037639758496e-3_wp*scale
  a(8)=   -0.2398228524507599670405555359023135e-4_wp*scale
  a(9)=    2.0904234952920365957922889447361e-6_wp*scale
  a(10)=  -3.7230763047369275848791496973044e-7_wp*scale
  a(11)=  -1.05857055496741470373494132287e-8_wp*scale
  a(12)=  -5.813879830282540547959250667e-11_wp*scale
  a(13)=   2.70800493626319438269856689037647576e-13_wp*scale
  a(14)=  -6.924474940639200152025730585882e-18_wp*scale

  a(15)=0.0_wp
  a(16)=0.0_wp 
  a(17)=0.0_wp
  
  do i=1,14+3
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
  c(-17)=0.0_wp
  c(-16)=0.0_wp
  c(-15)=0.0_wp
  
  c(-14)=     -3.869102413147656535541850057188e-18_wp*scale
  c(-13)=      1.5130616560866154733900029272077362e-13_wp*scale
  c(-12)=     -3.2264702314010525539061647271983988409e-11_wp*scale
  c(-11)=     -5.96264938781402337319841002642e-9_wp*scale
  c(-10)=     -2.1656830629214041470164889350342e-7_wp*scale
  c(-9 )=      8.7969704055286288323596890609625e-7_wp*scale
  c(-8 )=     -0.00001133456724516819987751818232711775_wp*scale
  c(-7 )=      0.00021710795484646138591610188464622454_wp*scale
  c(-6 )=     -0.0021356291838797986414312219042358542_wp*scale
  c(-5 )=      0.00713761218453631422925717625758502986_wp*scale
  c(-4 )=     -0.0284696165863973422636410524436931061_wp*scale
  c(-3 )=      0.14327329352510759457155821037742893841_wp*scale
  c(-2 )=     -0.42498050943780130143385739554118569733_wp*scale
  c(-1 )=      0.65703074007121357894896358254040272157_wp*scale
  c( 0 )=     -0.42081655293724308770919536332797729898_wp*scale
  c( 1 )=     -0.21716117505137104371463587747283267899_wp*scale
  c( 2 )=      0.63457035267892488185929915286969303251_wp*scale
  c( 3 )=     -0.53298223962800395684936080758073568406_wp*scale
  c( 4 )=      0.23370490631751294307619384973520033236_wp*scale
  c( 5 )=     -0.05657736973328755112051544344507997075_wp*scale
  c( 6 )=      0.0080872029411844780634067667008050127_wp*scale
  c( 7 )=     -0.00093423623304808664741804536808932984_wp*scale
  c( 8 )=      0.00005075807947289728306309081261461095_wp*scale
  c( 9 )=     -4.62561497463184262755416490048242e-6_wp*scale
  c( 10)=      6.3919128513793415587294752371778e-7_wp*scale
  c( 11)=      1.87909235155149902916133888931e-8_wp*scale
  c( 12)=      1.04757345962781829480207861447155543883e-10_wp*scale
  c( 13)=     -4.84665690596158959648731537084025836e-13_wp*scale
  c( 14)=      1.2392629629188986192855777620877e-17_wp*scale

  c(15)=0.0_wp
  c(16)=0.0_wp
  c(17)=0.0_wp
  !  <psi|D^2|phi_i>
  do i=-14-3,14+3
     b(i)=c(-i)
  enddo
  !<psi|D^2|psi_i>
  e(0)=   -24.875846029392331358907766562_wp*scale
  e(1)=   -7.1440597663471719869313377994_wp*scale
  e(2)=   -0.04251705323669172315864542163525830944_wp*scale
  e(3)=   -0.26995931336279126953587091167128839196_wp*scale
  e(4)=    0.08207454169225172612513390763444496516_wp*scale
  e(5)=   -0.02207327034586634477996701627614752761_wp*scale
  e(6)=    0.00409765642831595181639002667514310145_wp*scale
  e(7)=   -0.00045167920287507774929432548999880117_wp*scale
  e(8)=    0.00002398228524507599670405555359023135_wp*scale
  e(9)=   -2.0904234952920365957922889447361e-6_wp*scale
  e(10)=   3.7230763047369275848791496973044e-7_wp*scale
  e(11)=   1.05857055496741470373494132287e-8_wp*scale
  e(12)=   5.8138798302825405479592506674648873655e-11_wp*scale
  e(13)=  -2.70800493626319438269856689037647576e-13_wp*scale
  e(14)=   6.924474940639200152025730585882e-18_wp*scale
  do i=1,14
    e(-i)=e(i)
  enddo


! uncomment by Huan --- begin
!
!  if (firstcall) then
!
!     ! (1/2) d^2/dx^2
!     mflop1=0
!     do i3=0,n3
!        do i2=0,n2
!           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
!              do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!                 mflop1=mflop1+2
!              enddo
!              mflop1=mflop1+2
!           enddo
!            do i1=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil),&
!                  min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)
!                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_f(2,i2,i3)-i1)
!                    mflop1=mflop1+2
!                enddo
!                mflop1=mflop1+1
!            enddo
!            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!               do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!                  mflop1=mflop1+2
!               enddo
!            enddo
!     enddo
!  enddo
!     ! + (1/2) d^2/dy^2
!    mflop2=0
!    do i3=0,n3
!        do i1=0,n1
!            do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
!                   do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                   mflop2=mflop2+2       
!                   enddo
!                mflop2=mflop2+1
!            enddo
!            do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!               do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                    mflop2=mflop2+2
!               enddo
!               mflop2=mflop2+1
!            enddo
!            do i2=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil),&
!                  min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)
!               do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
!                  mflop2=mflop2+2
!               enddo
!            enddo
!        enddo
!    enddo
!     ! + (1/2) d^2/dz^2
!
!    mflop3=0
!    do i2=0,n2
!        do i1=0,n1
!            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
!                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+2
!                   enddo
!                mflop3=mflop3+1
!            enddo
!            do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!               do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+2
!               enddo
!               mflop3=mflop3+1
!            enddo
!            do i3=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil),&
!                  min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)
!               do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
!                  mflop3=mflop3+2
!               enddo
!            enddo
!
!        enddo
!    enddo
!
!     ! wavelet part
!     ! (1/2) d^2/dx^2
!     nflop1=0
!     do i3=nfl3,nfu3
!        do i2=nfl2,nfu2
!           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
!                 nflop1=nflop1+26
!              enddo
!              nflop1=nflop1+17
!           enddo
!        enddo
!     enddo
!
!     ! + (1/2) d^2/dy^2
!     nflop2=0
!     do i3=nfl3,nfu3
!        do i1=nfl1,nfu1
!           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
!                 nflop2=nflop2+26
!              enddo
!              nflop2=nflop2+7
!           enddo
!        enddo
!     enddo
!
!     ! + (1/2) d^2/dz^2
!     nflop3=0
!     do i2=nfl2,nfu2
!        do i1=nfl1,nfu1
!           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
!                 nflop3=nflop3+26
!              enddo
!              nflop3=nflop3+7
!           enddo
!        enddo
!     enddo
!
!     firstcall=.false.
!  endif

! uncomment by Huan -- end

  !---------------------------------------------------------------------------

  ! Scaling function part

  ! write(*,*) 'ncount0-1',ncount0

  ! (1/2) d^2/dx^2

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)


!call system_clock(ncount0,ncount_rate,ncount_max)

!$omp parallel default(private) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!$omp shared(x_f1,x_f2,x_f3,a,b,c,e,ncount0)&
!$omp private(ncount1,ncount2,ncount3,ncount4,ncount5,ncount6)

  !$omp do schedule(static,1)
!  !$omp parallel do collapse(2)

  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*a(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*a(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*a(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*a(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=dyi0+cprecr*x_c(i1+0,i2,i3)
              y_c(i1+1,i2,i3)=dyi1+cprecr*x_c(i1+1,i2,i3)
              y_c(i1+2,i2,i3)=dyi2+cprecr*x_c(i1+2,i2,i3)
              y_c(i1+3,i2,i3)=dyi3+cprecr*x_c(i1+3,i2,i3)

           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0.0_wp 
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*a(t-i1)
           enddo
           y_c(i1,i2,i3)=dyi+cprecr*x_c(i1,i2,i3)
        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + x_f1(t,i2,i3)*b(t-i1-0)
                 dyi1=dyi1 + x_f1(t,i2,i3)*b(t-i1-1)
                 dyi2=dyi2 + x_f1(t,i2,i3)*b(t-i1-2)
                 dyi3=dyi3 + x_f1(t,i2,i3)*b(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=y_c(i1+0,i2,i3)+dyi0
              y_c(i1+1,i2,i3)=y_c(i1+1,i2,i3)+dyi1
              y_c(i1+2,i2,i3)=y_c(i1+2,i2,i3)+dyi2
              y_c(i1+3,i2,i3)=y_c(i1+3,i2,i3)+dyi3
           enddo
           istart=i1
        endif

        do i1=istart,iend
           dyi=0.0_wp
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + x_f1(t,i2,i3)*b(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*c(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*c(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*c(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*c(t-i1-3)
              enddo
              y_f(1,i1+0,i2,i3)=dyi0
              y_f(1,i1+1,i2,i3)=dyi1
              y_f(1,i1+2,i2,i3)=dyi2
              y_f(1,i1+3,i2,i3)=dyi3
           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0.0_wp 
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*c(t-i1)
           enddo
           y_f(1,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !$omp enddo !nowait
  
    !call system_clock(ncount1,ncount_rate,ncount_max)
    !tel=dble(ncount1-ncount0)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2

  !$omp do schedule(static,1)
!  !$omp parallel do collapse(2)
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*a(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*a(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*a(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*a(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           dyi=0.0_wp 
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*a(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo
        istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
        iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i2=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp

              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + x_f2(t,i1,i3)*b(t-i2-0)
                 dyi1=dyi1 + x_f2(t,i1,i3)*b(t-i2-1)
                 dyi2=dyi2 + x_f2(t,i1,i3)*b(t-i2-2)
                 dyi3=dyi3 + x_f2(t,i1,i3)*b(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
           enddo
           istart=i2
        endif

        do i2=istart,iend
           dyi=0.0_wp
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + x_f2(t,i1,i3)*b(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*c(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*c(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*c(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*c(t-i2-3)
              enddo
              y_f(2,i1,i2+0,i3)=dyi0
              y_f(2,i1,i2+1,i3)=dyi1
              y_f(2,i1,i2+2,i3)=dyi2
              y_f(2,i1,i2+3,i3)=dyi3
           enddo
           icur=i2
        else
           icur=ibxz_f(1,i1,i3)
        endif

        do i2=icur,ibxz_f(2,i1,i3)
           dyi=0.0_wp 
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*c(t-i2)
           enddo
           y_f(2,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !$omp enddo !nowait

    !call system_clock(ncount2,ncount_rate,ncount_max)
    !tel=dble(ncount2-ncount1)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

!   !$omp parallel do collapse(2)
  !$omp do schedule(static,1)
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*a(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*a(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*a(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*a(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi=0.0_wp
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*a(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo
        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + x_f3(t,i1,i2)*b(t-i3-0)
                 dyi1=dyi1 + x_f3(t,i1,i2)*b(t-i3-1)
                 dyi2=dyi2 + x_f3(t,i1,i2)*b(t-i3-2)
                 dyi3=dyi3 + x_f3(t,i1,i2)*b(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi=0.0_wp
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + x_f3(t,i1,i2)*b(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*c(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*c(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*c(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*c(t-i3-3)
              enddo
              y_f(4,i1,i2,i3+0)=dyi0
              y_f(4,i1,i2,i3+1)=dyi1
              y_f(4,i1,i2,i3+2)=dyi2
              y_f(4,i1,i2,i3+3)=dyi3
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi=0.0_wp 
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*c(t-i3)
           enddo
           y_f(4,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !$omp enddo !nowait
  
    !call system_clock(ncount3,ncount_rate,ncount_max)
    !tel=dble(ncount3-ncount2)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2


!   !omp parallel do collapse(3)
  !$omp do schedule(static,1)
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l) + x_f(5,i1+l,i2,i3)*b(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l) + x_f(3,i1+l,i2,i3)*b(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l) + x_f(7,i1+l,i2,i3)*b(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*c(l) + x_f(5,i1+l,i2,i3)*e(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*c(l) + x_f(3,i1+l,i2,i3)*e(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*c(l) + x_f(7,i1+l,i2,i3)*e(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112+cprecr*x_f(4,i1,i2,i3)
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121+cprecr*x_f(2,i1,i2,i3)
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211+cprecr*x_f(1,i1,i2,i3)
           y_f(6,i1,i2,i3)=t122+cprecr*x_f(6,i1,i2,i3)
           y_f(5,i1,i2,i3)=t212+cprecr*x_f(5,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221+cprecr*x_f(3,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222+cprecr*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo !nowait

    !call system_clock(ncount4,ncount_rate,ncount_max)
    !tel=dble(ncount4-ncount3)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
!  !$omp parallel do collapse(2)
  !$omp do schedule(static,1)
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l) + x_f(6,i1,i2+l,i3)*b(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l) + x_f(3,i1,i2+l,i3)*b(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*c(l) + x_f(6,i1,i2+l,i3)*e(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l) + x_f(7,i1,i2+l,i3)*b(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*c(l) + x_f(3,i1,i2+l,i3)*e(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*c(l) + x_f(7,i1,i2+l,i3)*e(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
        enddo
     enddo
  enddo
  !$omp enddo !nowait

    !call system_clock(ncount5,ncount_rate,ncount_max)
    !tel=dble(ncount5-ncount4)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !$omp do schedule(static,1)
!  !$omp parallel do collapse(3)
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l) + x_f(6,i1,i2,i3+l)*b(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l) + x_f(5,i1,i2,i3+l)*b(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*c(l) + x_f(6,i1,i2,i3+l)*e(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*c(l) + x_f(5,i1,i2,i3+l)*e(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l) + x_f(7,i1,i2,i3+l)*b(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*c(l) + x_f(7,i1,i2,i3+l)*e(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
        enddo
     enddo
  enddo
  !$omp enddo !nowait

  !call system_clock(ncount6,ncount_rate,ncount_max)
  !tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

  !tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'K:ALL PART',  &
  !tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

  !$omp end parallel

!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

END SUBROUTINE Convolkinetic

!!$!> One of the most CPU intensive routines
!!$!!   y = (kinetic energy operator)x + (cprec*I)x 
!!$subroutine Convolkinetic_SSE(n1,n2,n3, &
!!$     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
!!$     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)
!!$  use module_base
!!$  implicit none
!!$!dee
!!$  !integer :: iend_test,count_rate_test,count_max_test,istart_test
!!$  
!!$  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
!!$  real(wp), intent(in) :: cprecr
!!$  real(gp), intent(in) :: hgrid
!!$  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
!!$  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
!!$  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
!!$  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
!!$  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
!!$  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
!!$  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
!!$  !local variables
!!$  integer, parameter :: lowfil=-14,lupfil=14
!!$  !logical :: firstcall=.true. 
!!$  integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
!!$  !integer(8) :: clock0,clock1,clock2,ncount_max,ncount_rate
!!$
!!$  integer(8) :: ncount1,ncount2,ncount3,ncount4,ncount5,ncount6,ncount0,ncount_max,ncount_rate
!!$
!!$  integer :: i,t,i1,i2,i3
!!$  integer :: icur,istart,iend,l
!!$  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211,tel
!!$  !real(kind=8) :: tel
!!$  real(wp), dimension(-3+lowfil:lupfil+3) :: a,b,c
!!$  real(wp), dimension(lowfil:lupfil) :: e
!!$  !------------------ REZA ------------------------------
!!$  integer, save::ntimes=0
!!$  real(8)::time1,time2
!!$  real(8), save::time=0.d0,count=0.d0,time_f=0.d0
!!$  real(8), save::dtime1_p1=0.d0,dtime2_p1=0.d0,dtime3_p1=0.d0
!!$  real(8), save::dtime1_p2=0.d0,dtime2_p2=0.d0,dtime3_p2=0.d0
!!$  real(8), save::dtime1_p3=0.d0,dtime2_p3=0.d0,dtime3_p3=0.d0
!!$  real(8), save::dtime1_p4=0.d0,dtime2_p4=0.d0,dtime3_p4=0.d0
!!$  real(4), allocatable, dimension(:,:,:) :: x_c_s
!!$  real(4), allocatable, dimension(:,:,:,:) :: x_f_s
!!$  real(4), allocatable, dimension(:,:,:) :: y_c_s
!!$  real(4), allocatable, dimension(:,:,:,:) :: y_f_s
!!$  real(4), dimension(-3+lowfil:lupfil+3) :: a_s,b_s,c_s
!!$  real(4), dimension(lowfil:lupfil) :: e_s
!!$  real(4) :: cprecr_s
!!$  !------------------------------------------------------
!!$
!!$  scale=-.5_wp/real(hgrid**2,wp)
!!$  !---------------------------------------------------------------------------
!!$  ! second derivative filters for Daubechies 16
!!$  !  <phi|D^2|phi_i>
!!$  a(0)=   -3.5536922899131901941296809374_wp*scale
!!$  a(1)=    2.2191465938911163898794546405_wp*scale
!!$  a(2)=   -0.6156141465570069496314853949_wp*scale
!!$  a(3)=    0.2371780582153805636239247476_wp*scale
!!$  a(4)=   -0.0822663999742123340987663521_wp*scale
!!$  a(5)=    0.02207029188482255523789911295638968409_wp*scale
!!$  a(6)=   -0.409765689342633823899327051188315485e-2_wp*scale
!!$  a(7)=    0.45167920287502235349480037639758496e-3_wp*scale
!!$  a(8)=   -0.2398228524507599670405555359023135e-4_wp*scale
!!$  a(9)=    2.0904234952920365957922889447361e-6_wp*scale
!!$  a(10)=  -3.7230763047369275848791496973044e-7_wp*scale
!!$  a(11)=  -1.05857055496741470373494132287e-8_wp*scale
!!$  a(12)=  -5.813879830282540547959250667e-11_wp*scale
!!$  a(13)=   2.70800493626319438269856689037647576e-13_wp*scale
!!$  a(14)=  -6.924474940639200152025730585882e-18_wp*scale
!!$
!!$  a(15)=0.0_wp
!!$  a(16)=0.0_wp 
!!$  a(17)=0.0_wp
!!$  
!!$  do i=1,14+3
!!$     a(-i)=a(i)
!!$  enddo
!!$  !  <phi|D^2|psi_i>
!!$  c(-17)=0.0_wp
!!$  c(-16)=0.0_wp
!!$  c(-15)=0.0_wp
!!$  
!!$  c(-14)=     -3.869102413147656535541850057188e-18_wp*scale
!!$  c(-13)=      1.5130616560866154733900029272077362e-13_wp*scale
!!$  c(-12)=     -3.2264702314010525539061647271983988409e-11_wp*scale
!!$  c(-11)=     -5.96264938781402337319841002642e-9_wp*scale
!!$  c(-10)=     -2.1656830629214041470164889350342e-7_wp*scale
!!$  c(-9 )=      8.7969704055286288323596890609625e-7_wp*scale
!!$  c(-8 )=     -0.00001133456724516819987751818232711775_wp*scale
!!$  c(-7 )=      0.00021710795484646138591610188464622454_wp*scale
!!$  c(-6 )=     -0.0021356291838797986414312219042358542_wp*scale
!!$  c(-5 )=      0.00713761218453631422925717625758502986_wp*scale
!!$  c(-4 )=     -0.0284696165863973422636410524436931061_wp*scale
!!$  c(-3 )=      0.14327329352510759457155821037742893841_wp*scale
!!$  c(-2 )=     -0.42498050943780130143385739554118569733_wp*scale
!!$  c(-1 )=      0.65703074007121357894896358254040272157_wp*scale
!!$  c( 0 )=     -0.42081655293724308770919536332797729898_wp*scale
!!$  c( 1 )=     -0.21716117505137104371463587747283267899_wp*scale
!!$  c( 2 )=      0.63457035267892488185929915286969303251_wp*scale
!!$  c( 3 )=     -0.53298223962800395684936080758073568406_wp*scale
!!$  c( 4 )=      0.23370490631751294307619384973520033236_wp*scale
!!$  c( 5 )=     -0.05657736973328755112051544344507997075_wp*scale
!!$  c( 6 )=      0.0080872029411844780634067667008050127_wp*scale
!!$  c( 7 )=     -0.00093423623304808664741804536808932984_wp*scale
!!$  c( 8 )=      0.00005075807947289728306309081261461095_wp*scale
!!$  c( 9 )=     -4.62561497463184262755416490048242e-6_wp*scale
!!$  c( 10)=      6.3919128513793415587294752371778e-7_wp*scale
!!$  c( 11)=      1.87909235155149902916133888931e-8_wp*scale
!!$  c( 12)=      1.04757345962781829480207861447155543883e-10_wp*scale
!!$  c( 13)=     -4.84665690596158959648731537084025836e-13_wp*scale
!!$  c( 14)=      1.2392629629188986192855777620877e-17_wp*scale
!!$
!!$  c(15)=0.0_wp
!!$  c(16)=0.0_wp
!!$  c(17)=0.0_wp
!!$  !  <psi|D^2|phi_i>
!!$  do i=-14-3,14+3
!!$     b(i)=c(-i)
!!$  enddo
!!$  !<psi|D^2|psi_i>
!!$  e(0)=   -24.875846029392331358907766562_wp*scale
!!$  e(1)=   -7.1440597663471719869313377994_wp*scale
!!$  e(2)=   -0.04251705323669172315864542163525830944_wp*scale
!!$  e(3)=   -0.26995931336279126953587091167128839196_wp*scale
!!$  e(4)=    0.08207454169225172612513390763444496516_wp*scale
!!$  e(5)=   -0.02207327034586634477996701627614752761_wp*scale
!!$  e(6)=    0.00409765642831595181639002667514310145_wp*scale
!!$  e(7)=   -0.00045167920287507774929432548999880117_wp*scale
!!$  e(8)=    0.00002398228524507599670405555359023135_wp*scale
!!$  e(9)=   -2.0904234952920365957922889447361e-6_wp*scale
!!$  e(10)=   3.7230763047369275848791496973044e-7_wp*scale
!!$  e(11)=   1.05857055496741470373494132287e-8_wp*scale
!!$  e(12)=   5.8138798302825405479592506674648873655e-11_wp*scale
!!$  e(13)=  -2.70800493626319438269856689037647576e-13_wp*scale
!!$  e(14)=   6.924474940639200152025730585882e-18_wp*scale
!!$  do i=1,14
!!$    e(-i)=e(i)
!!$  enddo
!!$
!!$
!!$  ! Scaling function part
!!$
!!$!call system_clock(istart_test,count_rate_test,count_max_test)
!!$
!!$  !------------------ REZA ------------------------------
!!$    allocate(x_c_s(0:n1,0:n2,0:n3),x_f_s(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),y_c_s(0:n1,0:n2,0:n3))
!!$    allocate(y_f_s(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3))
!!$    x_c_s=real(x_c,4)
!!$    x_f_s=real(x_f,4)
!!$    cprecr_s=real(cprecr,4)
!!$    a_s=real(a,4)
!!$    b_s=real(b,4)
!!$    c_s=real(c,4)
!!$    e_s=real(e,4)
!!$    !---------------------------------------------------------------------------
!!$    !call cpu_time(time1)
!!$    !---------------------------------------------------------------------------
!!$    call system_clock(ncount0,ncount_rate,ncount_max)
!!$    call convolut_ib_sse_cxyz(n1,n2,n3,ibyz_c,ibxz_c,ibxy_c,cprecr_s,x_c_s,y_c_s,a_s(-14))
!!$    !---------------------------------------------------------------------------
!!$    call convolut_ib_sse_fx_wss(n1,n2,n3,ibyz_f,cprecr_s, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fx_sws_wws(n1,n2,n3,ibyz_f,cprecr_s, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fx_ssw_wsw(n1,n2,n3,ibyz_f,cprecr_s, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fx_sww_www(n1,n2,n3,ibyz_f,cprecr_s, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    !---------------------------------------------------------------------------
!!$    call convolut_ib_sse_fy_sws(n1,n2,n3,ibxz_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fy_wss_wws(n1,n2,n3,ibxz_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fy_ssw_sww(n1,n2,n3,ibxz_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fy_wsw_www(n1,n2,n3,ibxz_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    !---------------------------------------------------------------------------
!!$    call convolut_ib_sse_fz_ssw(n1,n2,n3,ibxy_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fz_wss_wsw(n1,n2,n3,ibxy_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fz_sws_sww(n1,n2,n3,ibxy_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    call convolut_ib_sse_fz_wws_www(n1,n2,n3,ibxy_f, &
!!$      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x_f_s,y_f_s,a_s(-14),b_s(-14),c_s(-14),e_s(-14))
!!$    !---------------------------------------------------------------------------
!!$    !call cpu_time(time2)
!!$    call system_clock(ncount1,ncount_rate,ncount_max)
!!$    tel=dble(ncount1-ncount0)/dble(ncount_rate)
!!$    write(99,*) 'Convolkinetic',tel,y_c_s,y_f_s
!!$    y_c=real(y_c_s,8)
!!$    y_f=real(y_f_s,8)
!!$    deallocate(x_c_s,x_f_s,y_c_s,y_f_s)
!!$  !------------------------------------------------------
!!$
!!$
!!$  END SUBROUTINE Convolkinetic_SSE


subroutine ConvolkineticT(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,ekinout,x_f1,x_f2,x_f3)
  !   y = y+(kinetic energy operator)x 
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: x_f3
  real(gp), intent(out) :: ekinout
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount_max,ncount_rate

!$ integer :: ncount1,ncount2,ncount3,ncount4,ncount5,ncount6,ncount0

  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l

  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211,ekin
  !real(kind=8) :: tel
  real(wp), dimension(-3+lowfil:lupfil+3) :: a,b,c
  real(wp), dimension(lowfil:lupfil) :: e
  real(wp)::ekinp

  scale=-.5_wp/real(hgrid**2,wp)
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0)=   -3.5536922899131901941296809374_wp*scale
  a(1)=    2.2191465938911163898794546405_wp*scale
  a(2)=   -0.6156141465570069496314853949_wp*scale
  a(3)=    0.2371780582153805636239247476_wp*scale
  a(4)=   -0.0822663999742123340987663521_wp*scale
  a(5)=    0.02207029188482255523789911295638968409_wp*scale
  a(6)=   -0.409765689342633823899327051188315485e-2_wp*scale
  a(7)=    0.45167920287502235349480037639758496e-3_wp*scale
  a(8)=   -0.2398228524507599670405555359023135e-4_wp*scale
  a(9)=    2.0904234952920365957922889447361e-6_wp*scale
  a(10)=  -3.7230763047369275848791496973044e-7_wp*scale
  a(11)=  -1.05857055496741470373494132287e-8_wp*scale
  a(12)=  -5.813879830282540547959250667e-11_wp*scale
  a(13)=   2.70800493626319438269856689037647576e-13_wp*scale
  a(14)=  -6.924474940639200152025730585882e-18_wp*scale

  a(15)=0._wp
  a(16)=0._wp 
  a(17)=0._wp
  
  do i=1,14+3
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
  c(-17)=0._wp
  c(-16)=0._wp
  c(-15)=0._wp
  
  c(-14)=     -3.869102413147656535541850057188e-18_wp*scale
  c(-13)=      1.5130616560866154733900029272077362e-13_wp*scale
  c(-12)=     -3.2264702314010525539061647271983988409e-11_wp*scale
  c(-11)=     -5.96264938781402337319841002642e-9_wp*scale
  c(-10)=     -2.1656830629214041470164889350342e-7_wp*scale
  c(-9 )=      8.7969704055286288323596890609625e-7_wp*scale
  c(-8 )=     -0.00001133456724516819987751818232711775_wp*scale
  c(-7 )=      0.00021710795484646138591610188464622454_wp*scale
  c(-6 )=     -0.0021356291838797986414312219042358542_wp*scale
  c(-5 )=      0.00713761218453631422925717625758502986_wp*scale
  c(-4 )=     -0.0284696165863973422636410524436931061_wp*scale
  c(-3 )=      0.14327329352510759457155821037742893841_wp*scale
  c(-2 )=     -0.42498050943780130143385739554118569733_wp*scale
  c(-1 )=      0.65703074007121357894896358254040272157_wp*scale
  c( 0 )=     -0.42081655293724308770919536332797729898_wp*scale
  c( 1 )=     -0.21716117505137104371463587747283267899_wp*scale
  c( 2 )=      0.63457035267892488185929915286969303251_wp*scale
  c( 3 )=     -0.53298223962800395684936080758073568406_wp*scale
  c( 4 )=      0.23370490631751294307619384973520033236_wp*scale
  c( 5 )=     -0.05657736973328755112051544344507997075_wp*scale
  c( 6 )=      0.0080872029411844780634067667008050127_wp*scale
  c( 7 )=     -0.00093423623304808664741804536808932984_wp*scale
  c( 8 )=      0.00005075807947289728306309081261461095_wp*scale
  c( 9 )=     -4.62561497463184262755416490048242e-6_wp*scale
  c( 10)=      6.3919128513793415587294752371778e-7_wp*scale
  c( 11)=      1.87909235155149902916133888931e-8_wp*scale
  c( 12)=      1.04757345962781829480207861447155543883e-10_wp*scale
  c( 13)=     -4.84665690596158959648731537084025836e-13_wp*scale
  c( 14)=      1.2392629629188986192855777620877e-17_wp*scale

  c(15)=0._wp
  c(16)=0._wp
  c(17)=0._wp
  !  <psi|D^2|phi_i>
  do i=-14-3,14+3
     b(i)=c(-i)
  enddo
  !<psi|D^2|psi_i>
  e(0)=   -24.875846029392331358907766562_wp*scale
  e(1)=   -7.1440597663471719869313377994_wp*scale
  e(2)=   -0.04251705323669172315864542163525830944_wp*scale
  e(3)=   -0.26995931336279126953587091167128839196_wp*scale
  e(4)=    0.08207454169225172612513390763444496516_wp*scale
  e(5)=   -0.02207327034586634477996701627614752761_wp*scale
  e(6)=    0.00409765642831595181639002667514310145_wp*scale
  e(7)=   -0.00045167920287507774929432548999880117_wp*scale
  e(8)=    0.00002398228524507599670405555359023135_wp*scale
  e(9)=   -2.0904234952920365957922889447361e-6_wp*scale
  e(10)=   3.7230763047369275848791496973044e-7_wp*scale
  e(11)=   1.05857055496741470373494132287e-8_wp*scale
  e(12)=   5.8138798302825405479592506674648873655e-11_wp*scale
  e(13)=  -2.70800493626319438269856689037647576e-13_wp*scale
  e(14)=   6.924474940639200152025730585882e-18_wp*scale
  do i=1,14
     e(-i)=e(i)
  enddo

!  if (firstcall) then
!
!     ! (1/2) d^2/dx^2
!     mflop1=0
!     do i3=0,n3
!        do i2=0,n2
!           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
!              do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!                 mflop1=mflop1+2
!              enddo
!              mflop1=mflop1+3
!           enddo
!            do i1=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil),&
!                  min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)
!                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_f(2,i2,i3)-i1)
!                    mflop1=mflop1+2
!                enddo
!                mflop1=mflop1+3
!            enddo
!        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!           do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!              mflop1=mflop1+2
!           enddo
!           mflop1=mflop1+3
!        enddo
!     enddo
!  enddo
!  
!     ! + (1/2) d^2/dy^2
!    mflop2=0
!    do i3=0,n3
!        do i1=0,n1
!            do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
!                   do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                   mflop2=mflop2+2       
!                   enddo
!                mflop2=mflop2+3
!            enddo
!            do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!               do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                    mflop2=mflop2+2
!               enddo
!               mflop2=mflop2+3
!            enddo
!            do i2=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil),&
!                  min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)
!               do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
!                  mflop2=mflop2+2
!               enddo
!               mflop2=mflop2+3
!            enddo
!        enddo
!    enddo
!     ! + (1/2) d^2/dz^2
!
!    mflop3=0
!    do i2=0,n2
!        do i1=0,n1
!            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
!                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+2
!                   enddo
!                mflop3=mflop3+3
!            enddo
!            do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!               do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+2
!               enddo
!               mflop3=mflop3+3
!            enddo
!            do i3=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil),&
!                  min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)
!               do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
!                  mflop3=mflop3+2
!               enddo
!               mflop3=mflop3+3
!            enddo
!
!        enddo
!    enddo
!  
!     ! wavelet part
!     ! (1/2) d^2/dx^2
!     nflop1=0
!     do i3=nfl3,nfu3
!        do i2=nfl2,nfu2
!           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
!                 nflop1=nflop1+26
!              enddo
!              nflop1=nflop1+21
!           enddo
!        enddo
!     enddo
!
!     ! + (1/2) d^2/dy^2
!     nflop2=0
!     do i3=nfl3,nfu3
!        do i1=nfl1,nfu1
!           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
!                 nflop2=nflop2+26
!              enddo
!              nflop2=nflop2+21
!           enddo
!        enddo
!     enddo
!
!     ! + (1/2) d^2/dz^2
!     nflop3=0
!     do i2=nfl2,nfu2
!        do i1=nfl1,nfu1
!           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
!                 nflop3=nflop3+26
!              enddo
!              nflop3=nflop3+21
!           enddo
!        enddo
!     enddo
!
!     firstcall=.false.
!  endif
!
!  !---------------------------------------------------------------------------
!
  ! Scaling function part

  !call system_clock(ncount0,ncount_rate,ncount_max)

!  ! (1/2) d^2/dx^2
!

  ekin=0._wp
!$omp parallel default(private) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(ekin,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!$omp shared(x_f1,x_f2,x_f3,a,b,c,e,ncount0)&
!$omp private(ncount1,ncount2,ncount3,ncount4,ncount5,ncount6)
  ekinp=0._wp

  !$omp do schedule(static,1)
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0._wp 
              dyi1=0._wp 
              dyi2=0._wp 
              dyi3=0._wp 
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*a(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*a(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*a(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*a(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=y_c(i1+0,i2,i3)+dyi0
              y_c(i1+1,i2,i3)=y_c(i1+1,i2,i3)+dyi1
              y_c(i1+2,i2,i3)=y_c(i1+2,i2,i3)+dyi2
              y_c(i1+3,i2,i3)=y_c(i1+3,i2,i3)+dyi3

              ekinp=ekinp+dyi0*x_c(i1+0,i2,i3)
              ekinp=ekinp+dyi1*x_c(i1+1,i2,i3)
              ekinp=ekinp+dyi2*x_c(i1+2,i2,i3)
              ekinp=ekinp+dyi3*x_c(i1+3,i2,i3)
           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0._wp 
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*a(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_c(i1,i2,i3)
        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0._wp
              dyi1=0._wp
              dyi2=0._wp
              dyi3=0._wp
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + x_f1(t,i2,i3)*b(t-i1-0)
                 dyi1=dyi1 + x_f1(t,i2,i3)*b(t-i1-1)
                 dyi2=dyi2 + x_f1(t,i2,i3)*b(t-i1-2)
                 dyi3=dyi3 + x_f1(t,i2,i3)*b(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=y_c(i1+0,i2,i3)+dyi0
              y_c(i1+1,i2,i3)=y_c(i1+1,i2,i3)+dyi1
              y_c(i1+2,i2,i3)=y_c(i1+2,i2,i3)+dyi2
              y_c(i1+3,i2,i3)=y_c(i1+3,i2,i3)+dyi3

              ekinp=ekinp+dyi0*x_c(i1+0,i2,i3)
              ekinp=ekinp+dyi1*x_c(i1+1,i2,i3)
              ekinp=ekinp+dyi2*x_c(i1+2,i2,i3)
              ekinp=ekinp+dyi3*x_c(i1+3,i2,i3)
           enddo
           istart=i1
        endif

        do i1=istart,iend
           dyi=0._wp
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + x_f1(t,i2,i3)*b(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_c(i1,i2,i3)
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0._wp 
              dyi1=0._wp 
              dyi2=0._wp 
              dyi3=0._wp 
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*c(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*c(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*c(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*c(t-i1-3)
              enddo
              y_f(1,i1+0,i2,i3)=y_f(1,i1+0,i2,i3)+dyi0
              y_f(1,i1+1,i2,i3)=y_f(1,i1+1,i2,i3)+dyi1
              y_f(1,i1+2,i2,i3)=y_f(1,i1+2,i2,i3)+dyi2
              y_f(1,i1+3,i2,i3)=y_f(1,i1+3,i2,i3)+dyi3

              ekinp=ekinp+dyi0*x_f(1,i1+0,i2,i3)
              ekinp=ekinp+dyi1*x_f(1,i1+1,i2,i3)
              ekinp=ekinp+dyi2*x_f(1,i1+2,i2,i3)
              ekinp=ekinp+dyi3*x_f(1,i1+3,i2,i3)
           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0._wp 
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*c(t-i1)
           enddo
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_f(1,i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo
  
    !call system_clock(ncount1,ncount_rate,ncount_max)
    !tel=dble(ncount1-ncount0)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:FIRST PART:x',tel,1.d-6*mflop1/tel
  !!
  !!  ! + (1/2) d^2/dy^2
  !!
  !$omp do schedule(static,1)
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0._wp 
              dyi1=0._wp 
              dyi2=0._wp 
              dyi3=0._wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*a(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*a(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*a(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*a(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3

              ekinp=ekinp+dyi0*x_c(i1,i2+0,i3)
              ekinp=ekinp+dyi1*x_c(i1,i2+1,i3)
              ekinp=ekinp+dyi2*x_c(i1,i2+2,i3)
              ekinp=ekinp+dyi3*x_c(i1,i2+3,i3)
           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           dyi=0._wp 
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*a(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_c(i1,i2,i3)
        enddo

        istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
        iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i2=istart,iend-4,4
              dyi0=0._wp
              dyi1=0._wp
              dyi2=0._wp
              dyi3=0._wp

              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + x_f2(t,i1,i3)*b(t-i2-0)
                 dyi1=dyi1 + x_f2(t,i1,i3)*b(t-i2-1)
                 dyi2=dyi2 + x_f2(t,i1,i3)*b(t-i2-2)
                 dyi3=dyi3 + x_f2(t,i1,i3)*b(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3

              ekinp=ekinp+dyi0*x_c(i1,i2+0,i3)
              ekinp=ekinp+dyi1*x_c(i1,i2+1,i3)
              ekinp=ekinp+dyi2*x_c(i1,i2+2,i3)
              ekinp=ekinp+dyi3*x_c(i1,i2+3,i3)
           enddo
           istart=i2
        endif

        do i2=istart,iend
           dyi=0._wp
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + x_f2(t,i1,i3)*b(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_c(i1,i2,i3)
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0._wp 
              dyi1=0._wp 
              dyi2=0._wp 
              dyi3=0._wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*c(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*c(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*c(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*c(t-i2-3)
              enddo
              y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
              y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+dyi1
              y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+dyi2
              y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+dyi3

              ekinp=ekinp+dyi0*x_f(2,i1,i2+0,i3)
              ekinp=ekinp+dyi1*x_f(2,i1,i2+1,i3)
              ekinp=ekinp+dyi2*x_f(2,i1,i2+2,i3)
              ekinp=ekinp+dyi3*x_f(2,i1,i2+3,i3)
           enddo
           icur=i2
        else
           icur=ibxz_f(1,i1,i3)
        endif

        do i2=icur,ibxz_f(2,i1,i3)
           dyi=0._wp 
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*c(t-i2)
           enddo
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_f(2,i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo
  !    
  !!
    !call system_clock(ncount2,ncount_rate,ncount_max)
    !tel=dble(ncount2-ncount1)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:FIRST PART:y',tel,1.d-6*mflop2/tel
  !!
  !!  ! + (1/2) d^2/dz^2
  !!
  !$omp do schedule(static,1)
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0._wp 
              dyi1=0._wp 
              dyi2=0._wp 
              dyi3=0._wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*a(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*a(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*a(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*a(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3

              ekinp=ekinp+dyi0*x_c(i1,i2,i3+0)
              ekinp=ekinp+dyi1*x_c(i1,i2,i3+1)
              ekinp=ekinp+dyi2*x_c(i1,i2,i3+2)
              ekinp=ekinp+dyi3*x_c(i1,i2,i3+3)
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi=0._wp
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*a(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_c(i1,i2,i3)
        enddo

        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0._wp
              dyi1=0._wp
              dyi2=0._wp
              dyi3=0._wp
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + x_f3(t,i1,i2)*b(t-i3-0)
                 dyi1=dyi1 + x_f3(t,i1,i2)*b(t-i3-1)
                 dyi2=dyi2 + x_f3(t,i1,i2)*b(t-i3-2)
                 dyi3=dyi3 + x_f3(t,i1,i2)*b(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3

              ekinp=ekinp+dyi0*x_c(i1,i2,i3+0)
              ekinp=ekinp+dyi1*x_c(i1,i2,i3+1)
              ekinp=ekinp+dyi2*x_c(i1,i2,i3+2)
              ekinp=ekinp+dyi3*x_c(i1,i2,i3+3)
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi=0._wp
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + x_f3(t,i1,i2)*b(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_c(i1,i2,i3)
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0._wp 
              dyi1=0._wp 
              dyi2=0._wp 
              dyi3=0._wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*c(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*c(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*c(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*c(t-i3-3)
              enddo
              y_f(4,i1,i2,i3+0)=y_f(4,i1,i2,i3+0)+dyi0
              y_f(4,i1,i2,i3+1)=y_f(4,i1,i2,i3+1)+dyi1
              y_f(4,i1,i2,i3+2)=y_f(4,i1,i2,i3+2)+dyi2
              y_f(4,i1,i2,i3+3)=y_f(4,i1,i2,i3+3)+dyi3

              ekinp=ekinp+dyi0*x_f(4,i1,i2,i3+0)
              ekinp=ekinp+dyi1*x_f(4,i1,i2,i3+1)
              ekinp=ekinp+dyi2*x_f(4,i1,i2,i3+2)
              ekinp=ekinp+dyi3*x_f(4,i1,i2,i3+3)
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi=0._wp 
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*c(t-i3)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+dyi
           ekinp=ekinp+dyi*x_f(4,i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo
   !call system_clock(ncount3,ncount_rate,ncount_max)
   !tel=dble(ncount3-ncount2)/dble(ncount_rate)
   !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
  !$omp do schedule(static,1)
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l) + x_f(5,i1+l,i2,i3)*b(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l) + x_f(3,i1+l,i2,i3)*b(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l) + x_f(7,i1+l,i2,i3)*b(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*c(l) + x_f(5,i1+l,i2,i3)*e(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*c(l) + x_f(3,i1+l,i2,i3)*e(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*c(l) + x_f(7,i1+l,i2,i3)*e(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekinp=ekinp+t112*x_f(4,i1,i2,i3)
           ekinp=ekinp+t121*x_f(2,i1,i2,i3)
           ekinp=ekinp+t211*x_f(1,i1,i2,i3)
           ekinp=ekinp+t122*x_f(6,i1,i2,i3)
           ekinp=ekinp+t212*x_f(5,i1,i2,i3)
           ekinp=ekinp+t221*x_f(3,i1,i2,i3)
           ekinp=ekinp+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo

    !call system_clock(ncount4,ncount_rate,ncount_max)
    !tel=dble(ncount4-ncount3)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  !nb=16
  !$omp do schedule(static,1)
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l) + x_f(6,i1,i2+l,i3)*b(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l) + x_f(3,i1,i2+l,i3)*b(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*c(l) + x_f(6,i1,i2+l,i3)*e(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l) + x_f(7,i1,i2+l,i3)*b(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*c(l) + x_f(3,i1,i2+l,i3)*e(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*c(l) + x_f(7,i1,i2+l,i3)*e(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekinp=ekinp+t112*x_f(4,i1,i2,i3)
           ekinp=ekinp+t121*x_f(2,i1,i2,i3)
           ekinp=ekinp+t211*x_f(1,i1,i2,i3)
           ekinp=ekinp+t122*x_f(6,i1,i2,i3)
           ekinp=ekinp+t212*x_f(5,i1,i2,i3)
           ekinp=ekinp+t221*x_f(3,i1,i2,i3)
           ekinp=ekinp+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo

    !call system_clock(ncount5,ncount_rate,ncount_max)
    !tel=dble(ncount5-ncount4)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !nb=16
  !$omp do schedule(static,1)
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l) + x_f(6,i1,i2,i3+l)*b(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l) + x_f(5,i1,i2,i3+l)*b(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*c(l) + x_f(6,i1,i2,i3+l)*e(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*c(l) + x_f(5,i1,i2,i3+l)*e(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l) + x_f(7,i1,i2,i3+l)*b(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*c(l) + x_f(7,i1,i2,i3+l)*e(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekinp=ekinp+t112*x_f(4,i1,i2,i3)
           ekinp=ekinp+t121*x_f(2,i1,i2,i3)
           ekinp=ekinp+t211*x_f(1,i1,i2,i3)
           ekinp=ekinp+t122*x_f(6,i1,i2,i3)
           ekinp=ekinp+t212*x_f(5,i1,i2,i3)
           ekinp=ekinp+t221*x_f(3,i1,i2,i3)
           ekinp=ekinp+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo
  
    !call system_clock(ncount6,ncount_rate,ncount_max)
    !tel=dble(ncount6-ncount5)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:SECND PART:z',tel,1.d-6*nflop3/tel

    !tel=dble(ncount6-ncount0)/dble(ncount_rate)
    !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:ALL   PART',  & 
    !tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

  !$omp critical
  ekin=ekin+ekinp
  !$omp end critical
  
  !$omp end parallel

  !for them oment put the conversion of the kinetic energy here to avoid 
  !intensive conversions in the loops
  ekinout=real(ekin,gp)


  !  call system_clock(ncount6,ncount_rate,ncount_max)
  !  tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:SECND PART:z',tel,1.d-6*nflop3/tel

  !  tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:ALL   PART',  & 
  !  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

END SUBROUTINE ConvolkineticT
