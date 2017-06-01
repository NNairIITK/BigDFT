!> @file
!!  Simple convolution routines
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module for convolution operators
module convSimpleBench
  implicit none
  integer :: conv_f_nflop1,conv_f_nflop2,conv_f_nflop3
end module convSimpleBench


!> y = (kinetic energy operator)x + (cprec*I)x 
subroutine Convolkinetic_f(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_f,ibxz_f,ibxy_f,x_f,y_f)
  use module_defs, only: wp
  !n(c) use convSimpleBench, nflop1 => conv_f_nflop1, nflop2 => conv_f_nflop2, nflop3 => conv_f_nflop3
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(wp), intent(in) :: cprecr
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_f
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !n(c) logical :: firstcall=.true. 
  integer :: i,i1,i2,i3,l !n(c) t,ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  real(wp) :: scale,t112, t121,t122,t212,t221,t222,t211
  !n(c) real(kind=8) :: tel
  real(wp), dimension(lowfil:lupfil) :: a,e

  scale=-.5_wp/real(hgrid**2,wp)
  !---------------------------------------------------------------------------
  !< second derivative filters for Daubechies 16
  !!  <phi|D^2|phi_i>
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
  do i=1,14
     a(-i)=a(i)
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

  !call system_clock(ncount3,ncount_rate,ncount_max)
  !tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.e-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l) 
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l) 
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l) 
              t212=t212 + x_f(5,i1+l,i2,i3)*e(l)
              t221=t221 + x_f(3,i1+l,i2,i3)*e(l)
              t222=t222 + x_f(7,i1+l,i2,i3)*e(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=t112+cprecr*x_f(4,i1,i2,i3)
           y_f(2,i1,i2,i3)=t121+cprecr*x_f(2,i1,i2,i3)
           y_f(1,i1,i2,i3)=t211+cprecr*x_f(1,i1,i2,i3)
           y_f(6,i1,i2,i3)=t122+cprecr*x_f(6,i1,i2,i3)
           y_f(5,i1,i2,i3)=t212+cprecr*x_f(5,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221+cprecr*x_f(3,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222+cprecr*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo

  !call system_clock(ncount4,ncount_rate,ncount_max)
  !tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.e-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l)  
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l)  
              t122=t122 + x_f(6,i1,i2+l,i3)*e(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l)
              t221=t221 + x_f(3,i1,i2+l,i3)*e(l)
              t222=t222 + x_f(7,i1,i2+l,i3)*e(l)
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

  !call system_clock(ncount5,ncount_rate,ncount_max)
  !tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.e-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l)
              t122=t122 + x_f(6,i1,i2,i3+l)*e(l)
              t212=t212 + x_f(5,i1,i2,i3+l)*e(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l) 
              t222=t222 + x_f(7,i1,i2,i3+l)*e(l)
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

  !call system_clock(ncount6,ncount_rate,ncount_max)
  !tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.e-6*nflop3/tel

  !tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
  !tel,1.e-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

END SUBROUTINE Convolkinetic_f



!<   y = (kinetic energy operator)x + (cprec*I)x 
subroutine Convolkinetic_sep(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f) !n(c) x_f1,x_f2,x_f3 (arg:l, l-1)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(wp), intent(in) :: cprecr
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
  !n(c) real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f1
  !n(c) real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: x_f2
  !n(c) real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: x_f3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !n(c) logical :: firstcall=.true. 
  !n(c) integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  integer :: i,t,i1,i2,i3,l !n(c) ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  real(wp) :: scale,dyi,t112, t121,t122,t212,t221,t222,t211
  !n(c) real(kind=8) :: tel
  real(wp), dimension(lowfil:lupfil) :: a,c,e !n(c) b,d

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
  do i=1,14
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
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
  !  <psi|D^2|phi_i>
  !n(c) do i=-14,14
  !n(c)    b(i)=c(-i)
  !n(c) enddo
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

  ! Scaling function part

  !call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

    do i3=0,n3
        do i2=0,n2
            do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
                dyi=0._wp 
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                      dyi=dyi + x_c(t,i2,i3)*a(t-i1)
                   enddo
                   y_c(i1,i2,i3)=dyi+cprecr*x_c(i1,i2,i3)
            enddo
     enddo
  enddo
  !call system_clock(ncount1,ncount_rate,ncount_max)
  !tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.e-6*mflop1/tel

 ! + (1/2) d^2/dy^2
   do i3=0,n3
       do i1=0,n1
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
               dyi=0._wp 
                  do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                     dyi=dyi + x_c(i1,t,i3)*a(t-i2)
                  enddo
                  y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           enddo
    enddo
 enddo

 !call system_clock(ncount2,ncount_rate,ncount_max)
 !tel=dble(ncount2-ncount1)/dble(ncount_rate)
 !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.e-6*mflop2/tel

 ! ! + (1/2) d^2/dz^2

    do i2=0,n2
        do i1=0,n1
            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                dyi=0._wp
                do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                    dyi=dyi + x_c(i1,i2,t)*a(t-i3)
                   enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
            enddo
     enddo
  enddo

  !call system_clock(ncount3,ncount_rate,ncount_max)
  !tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.e-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l) 
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l) 
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l) 
              t212=t212 + x_f(5,i1+l,i2,i3)*e(l)
              t221=t221 + x_f(3,i1+l,i2,i3)*e(l)
              t222=t222 + x_f(7,i1+l,i2,i3)*e(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=t112+cprecr*x_f(4,i1,i2,i3)
           y_f(2,i1,i2,i3)=t121+cprecr*x_f(2,i1,i2,i3)
           y_f(1,i1,i2,i3)=t211+cprecr*x_f(1,i1,i2,i3)
           y_f(6,i1,i2,i3)=t122+cprecr*x_f(6,i1,i2,i3)
           y_f(5,i1,i2,i3)=t212+cprecr*x_f(5,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221+cprecr*x_f(3,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222+cprecr*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo

  !call system_clock(ncount4,ncount_rate,ncount_max)
  !tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.e-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l)  
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l)  
              t122=t122 + x_f(6,i1,i2+l,i3)*e(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l)
              t221=t221 + x_f(3,i1,i2+l,i3)*e(l)
              t222=t222 + x_f(7,i1,i2+l,i3)*e(l)
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

  !call system_clock(ncount5,ncount_rate,ncount_max)
  !tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.e-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l)
              t122=t122 + x_f(6,i1,i2,i3+l)*e(l)
              t212=t212 + x_f(5,i1,i2,i3+l)*e(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l) 
              t222=t222 + x_f(7,i1,i2,i3+l)*e(l)
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

  !call system_clock(ncount6,ncount_rate,ncount_max)
  !tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.e-6*nflop3/tel

  !tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
  !tel,1.e-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

END SUBROUTINE Convolkinetic_sep

!>  y = (kinetic energy operator)x + (cprec*I)x 
subroutine Convolkinetic(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3)
  use module_defs, only: wp
  implicit none
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
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  logical :: firstcall=.true. 
  integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  integer :: i,t,i1,i2,i3,l !n(c) ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  real(wp) :: scale,dyi,t112, t121,t122,t212,t221,t222,t211
  !n(c) real(kind=8) :: tel
  real(wp), dimension(lowfil:lupfil) :: a,b,c,e !n(c) d

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
  do i=1,14
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
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
  !  <psi|D^2|phi_i>
  do i=-14,14
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


  if (firstcall) then

     ! (1/2) d^2/dx^2
     mflop1=0
     do i3=0,n3
        do i2=0,n2
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
              do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
                 mflop1=mflop1+2
              enddo
              mflop1=mflop1+2
           enddo
            do i1=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil),&
                  min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)
                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_f(2,i2,i3)-i1)
                    mflop1=mflop1+2
                enddo
                mflop1=mflop1+1
            enddo
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
               do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
                  mflop1=mflop1+2
               enddo
            enddo
     enddo
  enddo
     ! + (1/2) d^2/dy^2
    mflop2=0
    do i3=0,n3
        do i1=0,n1
            do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
                   do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
                   mflop2=mflop2+2       
                   enddo
                mflop2=mflop2+1
            enddo
            do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
               do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
                    mflop2=mflop2+2
               enddo
               mflop2=mflop2+1
            enddo
            do i2=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil),&
                  min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)
               do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
                  mflop2=mflop2+2
               enddo
            enddo
        enddo
    enddo
     ! + (1/2) d^2/dz^2

    mflop3=0
    do i2=0,n2
        do i1=0,n1
            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
                    mflop3=mflop3+2
                   enddo
                mflop3=mflop3+1
            enddo
            do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
               do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
                    mflop3=mflop3+2
               enddo
               mflop3=mflop3+1
            enddo
            do i3=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil),&
                  min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)
               do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
                  mflop3=mflop3+2
               enddo
            enddo

        enddo
    enddo

     ! wavelet part
     ! (1/2) d^2/dx^2
     nflop1=0
     do i3=nfl3,nfu3
        do i2=nfl2,nfu2
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                 nflop1=nflop1+26
              enddo
              nflop1=nflop1+17
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     nflop2=0
     do i3=nfl3,nfu3
        do i1=nfl1,nfu1
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                 nflop2=nflop2+26
              enddo
              nflop2=nflop2+7
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dz^2
     nflop3=0
     do i2=nfl2,nfu2
        do i1=nfl1,nfu1
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 nflop3=nflop3+26
              enddo
              nflop3=nflop3+7
           enddo
        enddo
     enddo

     firstcall=.false.
  endif


  !---------------------------------------------------------------------------

  ! Scaling function part

  !call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

    do i3=0,n3
        do i2=0,n2
            do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
                dyi=0._wp 
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                      dyi=dyi + x_c(t,i2,i3)*a(t-i1)
                   enddo
                   y_c(i1,i2,i3)=dyi+cprecr*x_c(i1,i2,i3)
            enddo
            do i1=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil),&
                  min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)
                   dyi=0._wp
                do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
                    dyi=dyi + x_f1(t,i2,i3)*b(t-i1)
                enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
            enddo
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
               t211=0._wp 
               do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                  t211=t211 + x_c(t,i2,i3)*c(t-i1)
               enddo
               y_f(1,i1,i2,i3)=t211
            enddo
     enddo
  enddo
  !call system_clock(ncount1,ncount_rate,ncount_max)
  !tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.e-6*mflop1/tel

 ! + (1/2) d^2/dy^2
   do i3=0,n3
       do i1=0,n1
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
               dyi=0._wp 
                  do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                     dyi=dyi + x_c(i1,t,i3)*a(t-i2)
                  enddo
                  y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           enddo
           do i2=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil),&
                 min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)
              dyi=0._wp
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                 dyi=dyi + x_f2(t,i1,i3)*b(t-i2)
              enddo
              y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           enddo
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              t121=0._wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                 t121=t121 + x_c(i1,t,i3)*c(t-i2)
              enddo
              y_f(2,i1,i2,i3)=t121
           enddo
    enddo
 enddo


 !call system_clock(ncount2,ncount_rate,ncount_max)
 !tel=dble(ncount2-ncount1)/dble(ncount_rate)
 !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.e-6*mflop2/tel

 ! ! + (1/2) d^2/dz^2

    do i2=0,n2
        do i1=0,n1
            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                dyi=0._wp
                do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                    dyi=dyi + x_c(i1,i2,t)*a(t-i3)
                   enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
            enddo
            do i3=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil),&
                  min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)
               dyi=0._wp
               do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
                  dyi=dyi + x_f3(t,i1,i2)*b(t-i3)
               enddo
               y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
            enddo
            do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
               t112=0._wp 
               do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                  t112=t112 + x_c(i1,i2,t)*c(t-i3)
               enddo
               y_f(4,i1,i2,i3)=t112
            enddo
     enddo
  enddo

  !call system_clock(ncount3,ncount_rate,ncount_max)
  !tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.e-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
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

  !call system_clock(ncount4,ncount_rate,ncount_max)
  !tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.e-6*nflop1/tel


  ! + (1/2) d^2/dy^2
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
        enddo
     enddo
  enddo

  !call system_clock(ncount5,ncount_rate,ncount_max)
  !tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.e-6*nflop2/tel

  ! + (1/2) d^2/dz^2
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

        enddo
     enddo
  enddo

  !call system_clock(ncount6,ncount_rate,ncount_max)
  !tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.e-6*nflop3/tel

  !tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
  !tel,1.e-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

END SUBROUTINE Convolkinetic


!   y = (kinetic energy operator)x + (cprec*I)x 
subroutine Convolkinetic_c(n1,n2,n3,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,x_c,y_c,fac)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), intent(in) :: cprecr
  real(gp), intent(in) :: hgrid,fac
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i,t,i1,i2,i3 !n(c) l
  real(wp) :: scale,dyi
  real(wp), dimension(lowfil:lupfil) :: a

  scale=-.5_wp*fac/real(hgrid**2,wp)
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
  do i=1,14
     a(-i)=a(i)
  enddo

  !---------------------------------------------------------------------------

  ! Scaling function part

  ! (1/2) d^2/dx^2

    do i3=0,n3
        do i2=0,n2
            do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
                dyi=0._wp 
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                      dyi=dyi + x_c(t,i2,i3)*a(t-i1)
                   enddo
                   y_c(i1,i2,i3)=dyi+cprecr*x_c(i1,i2,i3)
            enddo
     enddo
  enddo

 ! + (1/2) d^2/dy^2
   do i3=0,n3
       do i1=0,n1
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
               dyi=0._wp 
                  do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                     dyi=dyi + x_c(i1,t,i3)*a(t-i2)
                  enddo
                  y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           enddo
    enddo
 enddo

  ! + (1/2) d^2/dz^2

    do i2=0,n2
        do i1=0,n1
            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                dyi=0._wp
                do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                    dyi=dyi + x_c(i1,i2,t)*a(t-i3)
                   enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
            enddo
     enddo
  enddo

END SUBROUTINE Convolkinetic_c



!> idir gives the dimension along which the convolutions shall be performed. It
!! is a three digit number (one for each dimension), each digit being either 1 (do the convolution) or 0 (don't do).
!! Example: 100 -> do in x, don't do in y and z
!!          010 -> don't do in x and z, do in y
!!          111 -> do in all three dimensions
subroutine ConvolkineticT(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hx,hy,hz,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,ekinout,x_f1,x_f2,x_f3,idir)
  !   y = y+(kinetic energy operator)x 
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,idir
  real(gp), intent(in) :: hx,hy,hz
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: x_f3
  real(gp), intent(out) :: ekinout
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  logical :: firstcall=.true. 
  integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  integer :: i,t,i1,i2,i3,ii,jj,kk,num !n(c) nb,ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: l !n(c) icur,istart,iend,j
  real(wp) :: scale,dyi,t112,t121,t122,t212,t221,t222,t211,ekin !n(c) dyi0,dyi1,dyi2,dyi3 
  !n(c) real(kind=8) :: tel
  real(wp), dimension(lowfil:lupfil) :: a,ax,ay,az,b,bx,by,bz,c,cx,cy,cz,e,ex,ey,ez !n(c) d
  logical,dimension(3) :: dodim

  !scale=-.5_wp/hgrid**2
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0)=   -3.5536922899131901941296809374_wp!*scale
  a(1)=    2.2191465938911163898794546405_wp!*scale
  a(2)=   -0.6156141465570069496314853949_wp!*scale
  a(3)=    0.2371780582153805636239247476_wp!*scale
  a(4)=   -0.0822663999742123340987663521_wp!*scale
  a(5)=    0.02207029188482255523789911295638968409_wp!*scale
  a(6)=   -0.409765689342633823899327051188315485e-2_wp!*scale
  a(7)=    0.45167920287502235349480037639758496e-3_wp!*scale
  a(8)=   -0.2398228524507599670405555359023135e-4_wp!*scale
  a(9)=    2.0904234952920365957922889447361e-6_wp!*scale
  a(10)=  -3.7230763047369275848791496973044e-7_wp!*scale
  a(11)=  -1.05857055496741470373494132287e-8_wp!*scale
  a(12)=  -5.813879830282540547959250667e-11_wp!*scale
  a(13)=   2.70800493626319438269856689037647576e-13_wp!*scale
  a(14)=  -6.924474940639200152025730585882e-18_wp!*scale
  do i=1,14
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
  c(-14)=     -3.869102413147656535541850057188e-18_wp!*scale
  c(-13)=      1.5130616560866154733900029272077362e-13_wp!*scale
  c(-12)=     -3.2264702314010525539061647271983988409e-11_wp!*scale
  c(-11)=     -5.96264938781402337319841002642e-9_wp!*scale
  c(-10)=     -2.1656830629214041470164889350342e-7_wp!*scale
  c(-9 )=      8.7969704055286288323596890609625e-7_wp!*scale
  c(-8 )=     -0.00001133456724516819987751818232711775_wp!*scale
  c(-7 )=      0.00021710795484646138591610188464622454_wp!*scale
  c(-6 )=     -0.0021356291838797986414312219042358542_wp!*scale
  c(-5 )=      0.00713761218453631422925717625758502986_wp!*scale
  c(-4 )=     -0.0284696165863973422636410524436931061_wp!*scale
  c(-3 )=      0.14327329352510759457155821037742893841_wp!*scale
  c(-2 )=     -0.42498050943780130143385739554118569733_wp!*scale
  c(-1 )=      0.65703074007121357894896358254040272157_wp!*scale
  c( 0 )=     -0.42081655293724308770919536332797729898_wp!*scale
  c( 1 )=     -0.21716117505137104371463587747283267899_wp!*scale
  c( 2 )=      0.63457035267892488185929915286969303251_wp!*scale
  c( 3 )=     -0.53298223962800395684936080758073568406_wp!*scale
  c( 4 )=      0.23370490631751294307619384973520033236_wp!*scale
  c( 5 )=     -0.05657736973328755112051544344507997075_wp!*scale
  c( 6 )=      0.0080872029411844780634067667008050127_wp!*scale
  c( 7 )=     -0.00093423623304808664741804536808932984_wp!*scale
  c( 8 )=      0.00005075807947289728306309081261461095_wp!*scale
  c( 9 )=     -4.62561497463184262755416490048242e-6_wp!*scale
  c( 10)=      6.3919128513793415587294752371778e-7_wp!*scale
  c( 11)=      1.87909235155149902916133888931e-8_wp!*scale
  c( 12)=      1.04757345962781829480207861447155543883e-10_wp!*scale
  c( 13)=     -4.84665690596158959648731537084025836e-13_wp!*scale
  c( 14)=      1.2392629629188986192855777620877e-17_wp!*scale
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i)=c(-i)
  enddo
  !<psi|D^2|psi_i>
  e(0)=   -24.875846029392331358907766562_wp!*scale
  e(1)=   -7.1440597663471719869313377994_wp!*scale
  e(2)=   -0.04251705323669172315864542163525830944_wp!*scale
  e(3)=   -0.26995931336279126953587091167128839196_wp!*scale
  e(4)=    0.08207454169225172612513390763444496516_wp!*scale
  e(5)=   -0.02207327034586634477996701627614752761_wp!*scale
  e(6)=    0.00409765642831595181639002667514310145_wp!*scale
  e(7)=   -0.00045167920287507774929432548999880117_wp!*scale
  e(8)=    0.00002398228524507599670405555359023135_wp!*scale
  e(9)=   -2.0904234952920365957922889447361e-6_wp!*scale
  e(10)=   3.7230763047369275848791496973044e-7_wp!*scale
  e(11)=   1.05857055496741470373494132287e-8_wp!*scale
  e(12)=   5.8138798302825405479592506674648873655e-11_wp!*scale
  e(13)=  -2.70800493626319438269856689037647576e-13_wp!*scale
  e(14)=   6.924474940639200152025730585882e-18_wp!*scale
  do i=1,14
     e(-i)=e(i)
  enddo

  scale=-.5_wp/real(hx**2,wp)
  ax = scale*a(:)
  bx = scale*b(:)
  cx = scale*c(:)
  ex = scale*e(:)

  scale=-.5_wp/real(hy**2,wp)
  ay = scale*a(:)
  by = scale*b(:)
  cy = scale*c(:)
  ey = scale*e(:)

  scale=-.5_wp/real(hz**2,wp)
  az = scale*a(:)
  bz = scale*b(:)
  cz = scale*c(:)
  ez = scale*e(:)

  ! Determine in which direction the convolutions should be performed
  num = idir
  do i=1,3
      jj = 10**i
      kk = 10**(i-1)
      ii = mod(num,jj)
      if (ii==kk) then
          dodim(4-i) = .true.
      else if (ii==0) then
          dodim(4-i) = .false.
      else
          stop 'wrong value of idir'
      end if
      num = num - ii
  end do

  if (firstcall) then

  if (dodim(1)) then
     ! (1/2) d^2/dx^2
     mflop1=0
     do i3=0,n3
        do i2=0,n2
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
              do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
                 mflop1=mflop1+2
              enddo
              mflop1=mflop1+3
           enddo
            do i1=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil),&
                  min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)
                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_f(2,i2,i3)-i1)
                    mflop1=mflop1+2
                enddo
                mflop1=mflop1+3
            enddo
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
              mflop1=mflop1+2
           enddo
           mflop1=mflop1+3
        enddo
     enddo
  enddo
  end if
  
  if (dodim(2)) then
     ! + (1/2) d^2/dy^2
    mflop2=0
    do i3=0,n3
        do i1=0,n1
            do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
                   do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
                   mflop2=mflop2+2       
                   enddo
                mflop2=mflop2+3
            enddo
            do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
               do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
                    mflop2=mflop2+2
               enddo
               mflop2=mflop2+3
            enddo
            do i2=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil),&
                  min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)
               do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
                  mflop2=mflop2+2
               enddo
               mflop2=mflop2+3
            enddo
        enddo
    enddo
  end if
     ! + (1/2) d^2/dz^2

  if (dodim(3)) then
    mflop3=0
    do i2=0,n2
        do i1=0,n1
            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
                    mflop3=mflop3+2
                   enddo
                mflop3=mflop3+3
            enddo
            do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
               do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
                    mflop3=mflop3+2
               enddo
               mflop3=mflop3+3
            enddo
            do i3=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil),&
                  min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)
               do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
                  mflop3=mflop3+2
               enddo
               mflop3=mflop3+3
            enddo

        enddo
    enddo
  end if
  
     ! wavelet part
     ! (1/2) d^2/dx^2
   if (dodim(1)) then
     nflop1=0
     do i3=nfl3,nfu3
        do i2=nfl2,nfu2
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                 nflop1=nflop1+26
              enddo
              nflop1=nflop1+21
           enddo
        enddo
     enddo
   end if

     ! + (1/2) d^2/dy^2
  if (dodim(2)) then
     nflop2=0
     do i3=nfl3,nfu3
        do i1=nfl1,nfu1
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                 nflop2=nflop2+26
              enddo
              nflop2=nflop2+21
           enddo
        enddo
     enddo
   end if

     ! + (1/2) d^2/dz^2
  if (dodim(3)) then
     nflop3=0
     do i2=nfl2,nfu2
        do i1=nfl1,nfu1
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 nflop3=nflop3+26
              enddo
              nflop3=nflop3+21
           enddo
        enddo
     enddo
   end if

     firstcall=.false.
  endif

  !---------------------------------------------------------------------------

  ! Scaling function part

  !call system_clock(ncount0,ncount_rate,ncount_max)
  ekin=0._gp

!  ! (1/2) d^2/dx^2
!

if (dodim(1)) then
    do i3=0,n3
        do i2=0,n2
            do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
                dyi=0._wp 
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                      dyi=dyi + x_c(t,i2,i3)*ax(t-i1)
                   enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
                   ekin=ekin+dyi*x_c(i1,i2,i3)
            enddo
            
            do i1=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil),&
                  min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)
                   dyi=0._wp
                do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
                    dyi=dyi + x_f1(t,i2,i3)*bx(t-i1)
                enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
                   ekin=ekin+dyi*x_c(i1,i2,i3)
            enddo
            
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
               t211=0._wp 
               do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                  t211=t211 + x_c(t,i2,i3)*cx(t-i1)
               enddo
               y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
               ekin=ekin+t211*x_f(1,i1,i2,i3)
            enddo
     enddo
  enddo
end if
  !call system_clock(ncount1,ncount_rate,ncount_max)
  !tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:FIRST PART:x',tel,1.e-6*mflop1/tel
!!
!!  ! + (1/2) d^2/dy^2
!!
if (dodim(2)) then
    do i3=0,n3
        do i1=0,n1
            do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
                dyi=0._wp 
                   do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                      dyi=dyi + x_c(i1,t,i3)*ay(t-i2)
                   enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
                   ekin=ekin+dyi*x_c(i1,i2,i3)
            enddo
            
            do i2=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil),&
                  min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)
               dyi=0._wp
               do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                  dyi=dyi + x_f2(t,i1,i3)*by(t-i2)
               enddo
               y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
               ekin=ekin+dyi*x_c(i1,i2,i3)
            enddo
            
            do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
               t121=0._wp 
               do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                  t121=t121 + x_c(i1,t,i3)*cy(t-i2)
               enddo
               y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
               ekin=ekin+t121*x_f(2,i1,i2,i3)
            enddo
     enddo
  enddo
end if
!    
!!
  !call system_clock(ncount2,ncount_rate,ncount_max)
  !tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:FIRST PART:y',tel,1.e-6*mflop2/tel
!!
!!  ! + (1/2) d^2/dz^2
!!
if (dodim(3)) then
    do i2=0,n2
        do i1=0,n1
            do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                dyi=0._wp
                do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                    dyi=dyi + x_c(i1,i2,t)*az(t-i3)
                   enddo
                   y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
                   ekin=ekin+dyi*x_c(i1,i2,i3)
            enddo
            
            do i3=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil),&
                  min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)
               dyi=0._wp
               do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
                  dyi=dyi + x_f3(t,i1,i2)*bz(t-i3)
               enddo
               y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
               ekin=ekin+dyi*x_c(i1,i2,i3)
            enddo

            do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
               t112=0._wp 
               do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                  t112=t112 + x_c(i1,i2,t)*cz(t-i3)
               enddo
               y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
               ekin=ekin+t112*x_f(4,i1,i2,i3)
            enddo
     enddo
  enddo
end if
  !call system_clock(ncount3,ncount_rate,ncount_max)
  !tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:FIRST PART:z',tel,1.e-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
if (dodim(1)) then
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*ax(l) + x_f(5,i1+l,i2,i3)*bx(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*ax(l) + x_f(3,i1+l,i2,i3)*bx(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*ax(l) + x_f(7,i1+l,i2,i3)*bx(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*cx(l) + x_f(5,i1+l,i2,i3)*ex(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*cx(l) + x_f(3,i1+l,i2,i3)*ex(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*cx(l) + x_f(7,i1+l,i2,i3)*ex(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*ex(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin=ekin+t112*x_f(4,i1,i2,i3)
           ekin=ekin+t121*x_f(2,i1,i2,i3)
           ekin=ekin+t211*x_f(1,i1,i2,i3)
           ekin=ekin+t122*x_f(6,i1,i2,i3)
           ekin=ekin+t212*x_f(5,i1,i2,i3)
           ekin=ekin+t221*x_f(3,i1,i2,i3)
           ekin=ekin+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
end if

  !call system_clock(ncount4,ncount_rate,ncount_max)
  !tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:SECND PART:x',tel,1.e-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  !n(c) nb=16
if (dodim(2)) then
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*ay(l) + x_f(6,i1,i2+l,i3)*by(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*ay(l) + x_f(3,i1,i2+l,i3)*by(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*cy(l) + x_f(6,i1,i2+l,i3)*ey(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*ay(l) + x_f(7,i1,i2+l,i3)*by(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*cy(l) + x_f(3,i1,i2+l,i3)*ey(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*cy(l) + x_f(7,i1,i2+l,i3)*ey(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*ey(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin=ekin+t112*x_f(4,i1,i2,i3)
           ekin=ekin+t121*x_f(2,i1,i2,i3)
           ekin=ekin+t211*x_f(1,i1,i2,i3)
           ekin=ekin+t122*x_f(6,i1,i2,i3)
           ekin=ekin+t212*x_f(5,i1,i2,i3)
           ekin=ekin+t221*x_f(3,i1,i2,i3)
           ekin=ekin+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
end if

  !call system_clock(ncount5,ncount_rate,ncount_max)
  !tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:SECND PART:y',tel,1.e-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !n(c) nb=16
if (dodim(3)) then
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*az(l) + x_f(6,i1,i2,i3+l)*bz(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*az(l) + x_f(5,i1,i2,i3+l)*bz(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*cz(l) + x_f(6,i1,i2,i3+l)*ez(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*cz(l) + x_f(5,i1,i2,i3+l)*ez(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*az(l) + x_f(7,i1,i2,i3+l)*bz(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*cz(l) + x_f(7,i1,i2,i3+l)*ez(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*ez(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin=ekin+t112*x_f(4,i1,i2,i3)
           ekin=ekin+t121*x_f(2,i1,i2,i3)
           ekin=ekin+t211*x_f(1,i1,i2,i3)
           ekin=ekin+t122*x_f(6,i1,i2,i3)
           ekin=ekin+t212*x_f(5,i1,i2,i3)
           ekin=ekin+t221*x_f(3,i1,i2,i3)
           ekin=ekin+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
end if

  !call system_clock(ncount6,ncount_rate,ncount_max)
  !tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:SECND PART:z',tel,1.e-6*nflop3/tel
  !
  !tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !write(99,'(a40,1x,e10.3,1x,f6.1)') 'T:ALL   PART',  & 
  !tel,1.e-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

  ekinout=real(ekin,gp)

END SUBROUTINE ConvolkineticT
