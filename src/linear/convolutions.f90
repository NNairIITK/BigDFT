subroutine getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, eff, filterCode)
!
! Purpose:
! ========
!   Calculates the effective filter for the operator [kineticEnergy + (x-x0)^4].
!   
! Calling arguments:
! ==================
!   Input arguments:
!     hgrid  grid spacing
!     x0     the center of the parabolic potential (x-x0)^2
!   Output arguments:
!     aeff   the effective filter for <phi|Op|phi>
!     beff   the effective filter for <psi|Op|phi>
!     ceff   the effective filter for <phi|Op|psi>
!     eeff   the effective filter for <psi|Op|psi>
!
use filterModule
implicit none

! Calling arguments
integer, intent(in):: it
real(8),intent(in):: parabPrefac, hgrid, x0
real(8),dimension(lb:ub),intent(out):: eff
character(len=*):: filterCode

! Local variables
integer:: i
real(8):: fac, fac2, prefac1, prefac2a, hgrid2, hgrid3, x02, x03
real(8):: scale
!scale=1.d0
scale=1.d-1
!scale=0.d-1
!scale=5.d-2
prefac1=-.5d0/hgrid**2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac*scale
fac2=parabPrefac*hgrid*scale
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate
select case(trim(filterCode))
case('a')
    do i=lb,ub
        !eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
        eff(i)=prefac1*a(i) + fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**4
case('b')
    do i=lb,ub
        !eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
        eff(i)=prefac1*b(i) + fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
case('c')
    do i=lb,ub
        !eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
        eff(i)=prefac1*c(i) + fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
case('e')
    do i=lb,ub
        !eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
        eff(i)=prefac1*e(i) + fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**4
case default
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end select


end subroutine getEffectiveFilterQuartic



subroutine getEffectiveFilterSextic(it,parabPrefac,hgrid, x0, eff, filterCode)
!
! Purpose:
! ========
!   Calculates the effective filter for the operator [kineticEnergy + (x-x0)^4].
!   
! Calling arguments:
! ==================
!   Input arguments:
!     hgrid  grid spacing
!     x0     the center of the parabolic potential (x-x0)^2
!   Output arguments:
!     aeff   the effective filter for <phi|Op|phi>
!     beff   the effective filter for <psi|Op|phi>
!     ceff   the effective filter for <phi|Op|psi>
!     eeff   the effective filter for <psi|Op|psi>
!
!use filterModule2
use filterModule
implicit none

! Calling arguments
integer, intent(in):: it
real(8),intent(in):: parabPrefac, hgrid, x0
real(8),dimension(lb:ub),intent(out):: eff
character(len=*):: filterCode

! Local variables
integer:: i
real(8):: fac, fac2, prefac1, prefac2a, hgrid2, hgrid3, hgrid4, hgrid5, x02, x03, x04, x05
real(8):: scale

!scale=1.d0
scale=5.d-2

prefac1=-.5d0/hgrid**2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac*scale
fac2=parabPrefac*hgrid*scale
hgrid2=hgrid**2
hgrid3=hgrid**3
hgrid4=hgrid**4
hgrid5=hgrid**5
x02=x0**2
x03=x0**3
x04=x0**4
x05=x0**5

! Determine which filter we have to calculate
select case(trim(filterCode))
case('a')
    do i=lb,ub
        !eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
        eff(i)=prefac1*a(i) + fac2*( hgrid5*a6(i) + 6.d0*hgrid4*x0*a5(i) + 15.d0*hgrid3*x02*a4(i) &
               + 20.d0*hgrid2*x03*a3(i) + 15.d0*hgrid*x04*a2(i) + 6.d0*x05*a1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**6
case('b')
    do i=lb,ub
        !eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
        eff(i)=prefac1*b(i) + fac2*( hgrid5*b6(i) + 6.d0*hgrid4*x0*b5(i) + 15.d0*hgrid3*x02*b4(i) &
               + 20.d0*hgrid2*x03*b3(i) + 15.d0*hgrid*x04*b2(i) + 6.d0*x05*b1(i))
        !eff(i)=prefac1*b(i) + fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
case('c')
    do i=lb,ub
        !eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
        eff(i)=prefac1*c(i) + fac2*( hgrid5*c6(i) + 6.d0*hgrid4*x0*c5(i) + 15.d0*hgrid3*x02*c4(i) &
               + 20.d0*hgrid2*x03*c3(i) + 15.d0*hgrid*x04*c2(i) + 6.d0*x05*c1(i))
        !eff(i)=prefac1*c(i) + fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
case('e')
    do i=lb,ub
        !eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
        eff(i)=prefac1*e(i) + fac2*( hgrid5*e6(i) + 6.d0*hgrid4*x0*e5(i) + 15.d0*hgrid3*x02*e4(i) &
               + 20.d0*hgrid2*x03*e3(i) + 15.d0*hgrid*x04*e2(i) + 6.d0*x05*e1(i))
        !eff(i)=prefac1*e(i) + fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**6
case default
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end select


end subroutine getEffectiveFilterSextic




subroutine getFilterQuartic(it,parabPrefac,hgrid, x0, eff, filterCode)
!
! Purpose:
! ========
!   Calculates the effective filter for the operator (x-x0)^4
!   
! Calling arguments:
! ==================
!   Input arguments:
!     hgrid  grid spacing
!     x0     the center of the parabolic potential (x-x0)^2
!   Output arguments:
!     aeff   the effective filter for <phi|Op|phi>
!     beff   the effective filter for <psi|Op|phi>
!     ceff   the effective filter for <phi|Op|psi>
!     eeff   the effective filter for <psi|Op|psi>
!
use filterModule
implicit none

! Calling arguments
integer, intent(in):: it
real(8),intent(in):: parabPrefac, hgrid, x0
real(8),dimension(lb:ub),intent(out):: eff
character(len=*):: filterCode

! Local variables
integer:: i
real(8):: fac, fac2, prefac2a, hgrid2, hgrid3, x02, x03
real(8):: scale
scale=1.d0
!scale=1.d-1
!scale=0.d-1
!scale=5.d-2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac*scale
fac2=parabPrefac*hgrid*scale
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate
select case(trim(filterCode))
case('a')
    do i=lb,ub
        !eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
        eff(i) = fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**4
case('b')
    do i=lb,ub
        !eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
        eff(i) = fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
case('c')
    do i=lb,ub
        !eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
        eff(i) = fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
case('e')
    do i=lb,ub
        !eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
        eff(i) = fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**4
case default
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end select


end subroutine getFilterQuartic



subroutine getFilterQuadratic(it,parabPrefac,hgrid, x0, eff, filterCode)
!
! Purpose:
! ========
!   Calculates the effective filter for the operator (x-x0)^2
!   
! Calling arguments:
! ==================
!   Input arguments:
!     hgrid  grid spacing
!     x0     the center of the parabolic potential (x-x0)^2
!   Output arguments:
!     aeff   the effective filter for <phi|Op|phi>
!     beff   the effective filter for <psi|Op|phi>
!     ceff   the effective filter for <phi|Op|psi>
!     eeff   the effective filter for <psi|Op|psi>
!
use filterModule
implicit none

! Calling arguments
integer, intent(in):: it
real(8),intent(in):: parabPrefac, hgrid, x0
real(8),dimension(lb:ub),intent(out):: eff
character(len=*):: filterCode

! Local variables
integer:: i
real(8):: fac, fac2, prefac2a, hgrid2, hgrid3, x02, x03
real(8):: scale
scale=1.d0
!scale=1.d-1
!scale=0.d-1
!scale=5.d-2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac*scale
fac2=parabPrefac*hgrid*scale
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate
select case(trim(filterCode))
case('a')
    do i=lb,ub
        eff(i) = fac2*( hgrid*a2(i) + 2.d0*x0*a1(i) )
        !eff(i) = fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
    !eff(0)=eff(0)+fac*x0**4
case('b')
    do i=lb,ub
        eff(i) = fac2*( hgrid*b2(i) + 2.d0*x0*b1(i) )
        !eff(i) = fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
case('c')
    do i=lb,ub
        eff(i) = fac2*( hgrid*c2(i) + 2.d0*x0*c1(i) )
        !eff(i) = fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
case('e')
    do i=lb,ub
        eff(i) = fac2*( hgrid*e2(i) + 2.d0*x0*e1(i) )
        !eff(i) = fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
    !eff(0)=eff(0)+fac*x0**4
case default
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end select


end subroutine getFilterQuadratic








!>  Applies the following operation: 
!!  y = [kinetic energy operator) + (cprec*I) + ((r-r0)^4)]*x
subroutine ConvolkineticQuartic(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzConf, potentialPrefac, it)
  use module_base
  implicit none

  ! Calling arguments
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
real(8),dimension(3):: rxyzConf
real(8):: potentialPrefac
integer:: it
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  !real(kind=8) :: tel
  real(wp), dimension(-3+lowfil:lupfil+3) :: a, aeff0, aeff1, aeff2, aeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: b, beff0, beff1, beff2, beff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: c, ceff0, ceff1, ceff2, ceff3
  real(wp), dimension(lowfil:lupfil) :: e, eeff0, eeff1, eeff2, eeff3
real(8):: x0, y0, z0
real(8):: x1, y1, z1
real(8):: x2, y2, z2
real(8):: x3, y3, z3
integer:: ii


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
!
!
!  !---------------------------------------------------------------------------

  ! Scaling function part

!  call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)



aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0


!!$!$omp parallel default(private) &
!!$!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!!$!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!!$!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
  !!!$omp do  
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the x dimension
              x0=hgrid*(i1+0)-rxyzConf(1)
              x1=hgrid*(i1+1)-rxyzConf(1)
              x2=hgrid*(i1+2)-rxyzConf(1)
              x3=hgrid*(i1+3)-rxyzConf(1)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x1, aeff1(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x2, aeff2(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x3, aeff3(lowfil), 'a')
!write(*,'(a,2i5,3x,2i6)') 'HERE: i2, i3, max(ibyz_c(1,i2,i3),lowfil+i1), min(lupfil+i1+3,ibyz_c(2,i2,i3))', i2, i3, max(ibyz_c(1,i2,i3),lowfil+i1), min(lupfil+i1+3,ibyz_c(2,i2,i3))
!write(*,'(a,2i6)') 'HERE: ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)', ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
!write(*,'(a,3i8)') 'HERE: t, i1, t-i1', t, i1, t-i1
                 dyi0=dyi0 + x_c(t,i2,i3)*aeff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*aeff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*aeff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*aeff3(t-i1-3)
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
           ! Get the effective a-filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*aeff0(t-i1)
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
              ! Get the effective b-filters for the x dimension
              x0=hgrid*(i1+0)-rxyzConf(1)
              x1=hgrid*(i1+1)-rxyzConf(1)
              x2=hgrid*(i1+2)-rxyzConf(1)
              x3=hgrid*(i1+3)-rxyzConf(1)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x1, beff1(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x2, beff2(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x3, beff3(lowfil), 'b')
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + x_f1(t,i2,i3)*beff0(t-i1-0)
                 dyi1=dyi1 + x_f1(t,i2,i3)*beff1(t-i1-1)
                 dyi2=dyi2 + x_f1(t,i2,i3)*beff2(t-i1-2)
                 dyi3=dyi3 + x_f1(t,i2,i3)*beff3(t-i1-3)
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
           ! Get the effective b-filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + x_f1(t,i2,i3)*beff0(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective c-filters for the x dimension
              x0=hgrid*(i1+0)-rxyzConf(1)
              x1=hgrid*(i1+1)-rxyzConf(1)
              x2=hgrid*(i1+2)-rxyzConf(1)
              x3=hgrid*(i1+3)-rxyzConf(1)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x1, ceff1(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x2, ceff2(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x3, ceff3(lowfil), 'c')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*ceff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*ceff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*ceff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*ceff3(t-i1-3)
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
           ! Get the effective c-filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*ceff0(t-i1)
           enddo
           y_f(1,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo
  
  !  call system_clock(ncount1,ncount_rate,ncount_max)
  !  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the y dimension
              y0=hgrid*(i2+0)-rxyzConf(2)
              y1=hgrid*(i2+1)-rxyzConf(2)
              y2=hgrid*(i2+2)-rxyzConf(2)
              y3=hgrid*(i2+3)-rxyzConf(2)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y1, aeff1(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y2, aeff2(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y3, aeff3(lowfil), 'a')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*aeff0(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*aeff1(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*aeff2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*aeff3(t-i2-3)
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
           ! Get the effective a-filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*aeff0(t-i2)
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
              ! Get the effective b-filters for the y dimension
              y0=hgrid*(i2+0)-rxyzConf(2)
              y1=hgrid*(i2+1)-rxyzConf(2)
              y2=hgrid*(i2+2)-rxyzConf(2)
              y3=hgrid*(i2+3)-rxyzConf(2)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y1, beff1(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y2, beff2(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y3, beff3(lowfil), 'b')
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + x_f2(t,i1,i3)*beff0(t-i2-0)
                 dyi1=dyi1 + x_f2(t,i1,i3)*beff1(t-i2-1)
                 dyi2=dyi2 + x_f2(t,i1,i3)*beff2(t-i2-2)
                 dyi3=dyi3 + x_f2(t,i1,i3)*beff3(t-i2-3)
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
           ! Get the effective b-filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + x_f2(t,i1,i3)*beff0(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective c-filters for the y dimension
              y0=hgrid*(i2+0)-rxyzConf(2)
              y1=hgrid*(i2+1)-rxyzConf(2)
              y2=hgrid*(i2+2)-rxyzConf(2)
              y3=hgrid*(i2+3)-rxyzConf(2)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y1, ceff1(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y2, ceff2(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y3, ceff3(lowfil), 'c')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*ceff0(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*ceff1(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*ceff2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*ceff3(t-i2-3)
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
           ! Get the effective c-filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*ceff0(t-i2)
           enddo
           y_f(2,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo


  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the z dimension
              z0=hgrid*(i3+0)-rxyzConf(3)
              z1=hgrid*(i3+1)-rxyzConf(3)
              z2=hgrid*(i3+2)-rxyzConf(3)
              z3=hgrid*(i3+3)-rxyzConf(3)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z1, aeff1(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z2, aeff2(lowfil), 'a')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z3, aeff3(lowfil), 'a')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*aeff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*aeff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*aeff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*aeff3(t-i3-3)
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
           ! Get the effective a-filters for the y dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*aeff0(t-i3)
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
              ! Get the effective b-filters for the z dimension
              z0=hgrid*(i3+0)-rxyzConf(3)
              z1=hgrid*(i3+1)-rxyzConf(3)
              z2=hgrid*(i3+2)-rxyzConf(3)
              z3=hgrid*(i3+3)-rxyzConf(3)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z1, beff1(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z2, beff2(lowfil), 'b')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z3, beff3(lowfil), 'b')
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + x_f3(t,i1,i2)*beff0(t-i3-0)
                 dyi1=dyi1 + x_f3(t,i1,i2)*beff1(t-i3-1)
                 dyi2=dyi2 + x_f3(t,i1,i2)*beff2(t-i3-2)
                 dyi3=dyi3 + x_f3(t,i1,i2)*beff3(t-i3-3)
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
           ! Get the effective b-filters for the y dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + x_f3(t,i1,i2)*beff0(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective c-filters for the z dimension
              z0=hgrid*(i3+0)-rxyzConf(3)
              z1=hgrid*(i3+1)-rxyzConf(3)
              z2=hgrid*(i3+2)-rxyzConf(3)
              z3=hgrid*(i3+3)-rxyzConf(3)
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, ceff0(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z1, ceff1(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z2, ceff2(lowfil), 'c')
              call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z3, ceff3(lowfil), 'c')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*ceff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*ceff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*ceff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*ceff3(t-i3-3)
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
           ! Get the effective c-filters for the z dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid,  z0, ceff0(lowfil), 'c')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*ceff0(t-i3)
           enddo
           y_f(4,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo


  
  !  call system_clock(ncount3,ncount_rate,ncount_max)
  !  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2

  !!!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           ! Get the effective filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, eeff0(lowfil), 'e')
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*aeff0(l) + x_f(5,i1+l,i2,i3)*beff0(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*aeff0(l) + x_f(3,i1+l,i2,i3)*beff0(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*aeff0(l) + x_f(7,i1+l,i2,i3)*beff0(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*ceff0(l) + x_f(5,i1+l,i2,i3)*eeff0(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*ceff0(l) + x_f(3,i1+l,i2,i3)*eeff0(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*ceff0(l) + x_f(7,i1+l,i2,i3)*eeff0(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*eeff0(l)
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
  !!!$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           ! Get the effective filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, eeff0(lowfil), 'e')
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*aeff0(l) + x_f(6,i1,i2+l,i3)*beff0(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*aeff0(l) + x_f(3,i1,i2+l,i3)*beff0(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*ceff0(l) + x_f(6,i1,i2+l,i3)*eeff0(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*aeff0(l) + x_f(7,i1,i2+l,i3)*beff0(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*ceff0(l) + x_f(3,i1,i2+l,i3)*eeff0(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*ceff0(l) + x_f(7,i1,i2+l,i3)*eeff0(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*eeff0(l)
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
  !!!$omp enddo

  !  call system_clock(ncount5,ncount_rate,ncount_max)
  !  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           ! Get the effective filters for the z dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, ceff0(lowfil), 'c')
           call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, eeff0(lowfil), 'e')
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*aeff0(l) + x_f(6,i1,i2,i3+l)*beff0(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*aeff0(l) + x_f(5,i1,i2,i3+l)*beff0(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*ceff0(l) + x_f(6,i1,i2,i3+l)*eeff0(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*ceff0(l) + x_f(5,i1,i2,i3+l)*eeff0(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*aeff0(l) + x_f(7,i1,i2,i3+l)*beff0(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*ceff0(l) + x_f(7,i1,i2,i3+l)*eeff0(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*eeff0(l)
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
  !!!$omp enddo

!  !!!$omp end parallel
!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel




END SUBROUTINE ConvolkineticQuartic






!>  Applies the following operation: 
!!  y = [kinetic energy operator) + (cprec*I) + ((r-r0)^6)]*x
subroutine ConvolkineticSextic(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzConf, potentialPrefac, it)
  use module_base
  implicit none

  ! Calling arguments
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
real(8),dimension(3):: rxyzConf
real(8):: potentialPrefac
integer:: it
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  !real(kind=8) :: tel
  real(wp), dimension(-3+lowfil:lupfil+3) :: a, aeff0, aeff1, aeff2, aeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: b, beff0, beff1, beff2, beff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: c, ceff0, ceff1, ceff2, ceff3
  real(wp), dimension(lowfil:lupfil) :: e, eeff0, eeff1, eeff2, eeff3
real(8):: x0, y0, z0
real(8):: x1, y1, z1
real(8):: x2, y2, z2
real(8):: x3, y3, z3
integer:: ii


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
!
!
!  !---------------------------------------------------------------------------

  ! Scaling function part

!  call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)



aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0


!!!$omp parallel default(private) &
!!!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!!!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!!!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
  !!!$omp do  
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the x dimension
              x0=hgrid*(i1+0)-rxyzConf(1)
              x1=hgrid*(i1+1)-rxyzConf(1)
              x2=hgrid*(i1+2)-rxyzConf(1)
              x3=hgrid*(i1+3)-rxyzConf(1)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x1, aeff1(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x2, aeff2(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x3, aeff3(lowfil), 'a')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*aeff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*aeff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*aeff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*aeff3(t-i1-3)
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
           ! Get the effective a-filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*aeff0(t-i1)
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
              ! Get the effective b-filters for the x dimension
              x0=hgrid*(i1+0)-rxyzConf(1)
              x1=hgrid*(i1+1)-rxyzConf(1)
              x2=hgrid*(i1+2)-rxyzConf(1)
              x3=hgrid*(i1+3)-rxyzConf(1)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x1, beff1(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x2, beff2(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x3, beff3(lowfil), 'b')
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + x_f1(t,i2,i3)*beff0(t-i1-0)
                 dyi1=dyi1 + x_f1(t,i2,i3)*beff1(t-i1-1)
                 dyi2=dyi2 + x_f1(t,i2,i3)*beff2(t-i1-2)
                 dyi3=dyi3 + x_f1(t,i2,i3)*beff3(t-i1-3)
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
           ! Get the effective b-filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + x_f1(t,i2,i3)*beff0(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective c-filters for the x dimension
              x0=hgrid*(i1+0)-rxyzConf(1)
              x1=hgrid*(i1+1)-rxyzConf(1)
              x2=hgrid*(i1+2)-rxyzConf(1)
              x3=hgrid*(i1+3)-rxyzConf(1)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x1, ceff1(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x2, ceff2(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x3, ceff3(lowfil), 'c')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*ceff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*ceff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*ceff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*ceff3(t-i1-3)
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
           ! Get the effective c-filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*ceff0(t-i1)
           enddo
           y_f(1,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo
  
  !  call system_clock(ncount1,ncount_rate,ncount_max)
  !  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the y dimension
              y0=hgrid*(i2+0)-rxyzConf(2)
              y1=hgrid*(i2+1)-rxyzConf(2)
              y2=hgrid*(i2+2)-rxyzConf(2)
              y3=hgrid*(i2+3)-rxyzConf(2)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y1, aeff1(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y2, aeff2(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y3, aeff3(lowfil), 'a')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*aeff0(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*aeff1(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*aeff2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*aeff3(t-i2-3)
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
           ! Get the effective a-filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*aeff0(t-i2)
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
              ! Get the effective b-filters for the y dimension
              y0=hgrid*(i2+0)-rxyzConf(2)
              y1=hgrid*(i2+1)-rxyzConf(2)
              y2=hgrid*(i2+2)-rxyzConf(2)
              y3=hgrid*(i2+3)-rxyzConf(2)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y1, beff1(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y2, beff2(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y3, beff3(lowfil), 'b')
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + x_f2(t,i1,i3)*beff0(t-i2-0)
                 dyi1=dyi1 + x_f2(t,i1,i3)*beff1(t-i2-1)
                 dyi2=dyi2 + x_f2(t,i1,i3)*beff2(t-i2-2)
                 dyi3=dyi3 + x_f2(t,i1,i3)*beff3(t-i2-3)
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
           ! Get the effective b-filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + x_f2(t,i1,i3)*beff0(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective c-filters for the y dimension
              y0=hgrid*(i2+0)-rxyzConf(2)
              y1=hgrid*(i2+1)-rxyzConf(2)
              y2=hgrid*(i2+2)-rxyzConf(2)
              y3=hgrid*(i2+3)-rxyzConf(2)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y1, ceff1(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y2, ceff2(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y3, ceff3(lowfil), 'c')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*ceff0(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*ceff1(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*ceff2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*ceff3(t-i2-3)
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
           ! Get the effective c-filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*ceff0(t-i2)
           enddo
           y_f(2,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo


  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the z dimension
              z0=hgrid*(i3+0)-rxyzConf(3)
              z1=hgrid*(i3+1)-rxyzConf(3)
              z2=hgrid*(i3+2)-rxyzConf(3)
              z3=hgrid*(i3+3)-rxyzConf(3)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z1, aeff1(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z2, aeff2(lowfil), 'a')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z3, aeff3(lowfil), 'a')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*aeff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*aeff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*aeff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*aeff3(t-i3-3)
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
           ! Get the effective a-filters for the y dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*aeff0(t-i3)
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
              ! Get the effective b-filters for the z dimension
              z0=hgrid*(i3+0)-rxyzConf(3)
              z1=hgrid*(i3+1)-rxyzConf(3)
              z2=hgrid*(i3+2)-rxyzConf(3)
              z3=hgrid*(i3+3)-rxyzConf(3)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z1, beff1(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z2, beff2(lowfil), 'b')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z3, beff3(lowfil), 'b')
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + x_f3(t,i1,i2)*beff0(t-i3-0)
                 dyi1=dyi1 + x_f3(t,i1,i2)*beff1(t-i3-1)
                 dyi2=dyi2 + x_f3(t,i1,i2)*beff2(t-i3-2)
                 dyi3=dyi3 + x_f3(t,i1,i2)*beff3(t-i3-3)
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
           ! Get the effective b-filters for the y dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + x_f3(t,i1,i2)*beff0(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective c-filters for the z dimension
              z0=hgrid*(i3+0)-rxyzConf(3)
              z1=hgrid*(i3+1)-rxyzConf(3)
              z2=hgrid*(i3+2)-rxyzConf(3)
              z3=hgrid*(i3+3)-rxyzConf(3)
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, ceff0(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z1, ceff1(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z2, ceff2(lowfil), 'c')
              call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z3, ceff3(lowfil), 'c')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*ceff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*ceff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*ceff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*ceff3(t-i3-3)
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
           ! Get the effective c-filters for the z dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid,  z0, ceff0(lowfil), 'c')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*ceff0(t-i3)
           enddo
           y_f(4,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo


  
  !  call system_clock(ncount3,ncount_rate,ncount_max)
  !  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2

  !!!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           ! Get the effective filters for the x dimension
           x0=hgrid*(i1+0)-rxyzConf(1)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, x0, eeff0(lowfil), 'e')
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*aeff0(l) + x_f(5,i1+l,i2,i3)*beff0(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*aeff0(l) + x_f(3,i1+l,i2,i3)*beff0(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*aeff0(l) + x_f(7,i1+l,i2,i3)*beff0(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*ceff0(l) + x_f(5,i1+l,i2,i3)*eeff0(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*ceff0(l) + x_f(3,i1+l,i2,i3)*eeff0(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*ceff0(l) + x_f(7,i1+l,i2,i3)*eeff0(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*eeff0(l)
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
  !!!$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           ! Get the effective filters for the y dimension
           y0=hgrid*(i2+0)-rxyzConf(2)
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, y0, eeff0(lowfil), 'e')
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*aeff0(l) + x_f(6,i1,i2+l,i3)*beff0(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*aeff0(l) + x_f(3,i1,i2+l,i3)*beff0(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*ceff0(l) + x_f(6,i1,i2+l,i3)*eeff0(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*aeff0(l) + x_f(7,i1,i2+l,i3)*beff0(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*ceff0(l) + x_f(3,i1,i2+l,i3)*eeff0(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*ceff0(l) + x_f(7,i1,i2+l,i3)*eeff0(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*eeff0(l)
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
  !!!$omp enddo

  !  call system_clock(ncount5,ncount_rate,ncount_max)
  !  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           ! Get the effective filters for the z dimension
           z0=hgrid*(i3+0)-rxyzConf(3)
           !call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, ceff0(lowfil), 'c')
           call getEffectiveFilterSextic(it,potentialPrefac,hgrid, z0, eeff0(lowfil), 'e')
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*aeff0(l) + x_f(6,i1,i2,i3+l)*beff0(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*aeff0(l) + x_f(5,i1,i2,i3+l)*beff0(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*ceff0(l) + x_f(6,i1,i2,i3+l)*eeff0(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*ceff0(l) + x_f(5,i1,i2,i3+l)*eeff0(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*aeff0(l) + x_f(7,i1,i2,i3+l)*beff0(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*ceff0(l) + x_f(7,i1,i2,i3+l)*eeff0(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*eeff0(l)
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
  !!!$omp enddo

  !!!$omp end parallel
!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel




END SUBROUTINE ConvolkineticSextic





!>  Applies the following operation: 
!!  y = [((r-r0)^4)]*x
subroutine ConvolQuartic(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid, offsetx, offsety, offsetz, &
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzConf, potentialPrefac, it)
  use module_base
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, offsetx, offsety, offsetz
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
real(8),dimension(3):: rxyzConf
real(8):: potentialPrefac
integer:: it
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp):: tt112, tt121, tt122, tt212, tt221, tt222, tt211
  real(wp):: tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7
  !real(kind=8) :: tel
  real(wp), dimension(-3+lowfil:lupfil+3) :: a, aeff0, aeff1, aeff2, aeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: b, beff0, beff1, beff2, beff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: c, ceff0, ceff1, ceff2, ceff3
  real(wp), dimension(lowfil:lupfil) :: e, eeff0, eeff1, eeff2, eeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: aeff0_2, aeff1_2, aeff2_2, aeff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: beff0_2, beff1_2, beff2_2, beff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: ceff0_2, ceff1_2, ceff2_2, ceff3_2
  real(wp), dimension(lowfil:lupfil) :: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  real(wp),dimension(:,:,:),allocatable:: x2_c, y2_c, z2_c
  real(wp),dimension(:,:,:,:),allocatable:: x2_f, y2_f, z2_f
real(8):: x0, y0, z0
real(8):: x1, y1, z1
real(8):: x2, y2, z2
real(8):: x3, y3, z3
integer:: ii, istat, iall
character(len=*),parameter:: subname='ConvolQuartic'


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
!
!
!  !---------------------------------------------------------------------------

  ! Scaling function part

!  call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)





allocate(x2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, x2_c, 'x2_c', subname)
allocate(y2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, y2_c, 'y2_c', subname)
allocate(z2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, z2_c, 'z2_c', subname)
allocate(x2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, x2_f, 'x2_f', subname)
allocate(y2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, y2_f, 'y2_f', subname)
allocate(z2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, z2_f, 'z2_f', subname)

x2_c=0.d0
y2_c=0.d0
z2_c=0.d0
x2_f=0.d0
y2_f=0.d0
z2_f=0.d0


aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0

aeff0_2=0.d0 ; beff0_2=0.d0 ; ceff0_2=0.d0 ; eeff0_2=0.0
aeff1_2=0.d0 ; beff1_2=0.d0 ; ceff1_2=0.d0 ; eeff1_2=0.0
aeff2_2=0.d0 ; beff2_2=0.d0 ; ceff2_2=0.d0 ; eeff2_2=0.0
aeff3_2=0.d0 ; beff3_2=0.d0 ; ceff3_2=0.d0 ; eeff3_2=0.0


!!$!$omp parallel default(private) &
!!$!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!!$!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!!$!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
  !!!$omp do  
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective a-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, aeff3_2(lowfil), 'a')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*aeff0(t-i1-0)
                 write(3100,'(a,3i9,2es16.7,7i10)') 't, i2, i3, x_c(t,i2,i3), aeff0(t-i1-0), n1, n2, n3, ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)-4, max(ibyz_c(1,i2,i3),lowfil+i1), min(lupfil+i1+3,ibyz_c(2,i2,i3))', &
                     t, i2, i3, x_c(t,i2,i3), aeff0(t-i1-0), n1, n2, n3, ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)-4, max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi1=dyi1 + x_c(t,i2,i3)*aeff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*aeff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*aeff3(t-i1-3)
                 tt0=tt0 + x_c(t,i2,i3)*aeff0_2(t-i1-0)
                 tt1=tt1 + x_c(t,i2,i3)*aeff1_2(t-i1-1)
                 tt2=tt2 + x_c(t,i2,i3)*aeff2_2(t-i1-2)
                 tt3=tt3 + x_c(t,i2,i3)*aeff3_2(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=dyi0!+cprecr*x_c(i1+0,i2,i3)
              write(3000,'(a,3i9,es16.7)') 'i1, i2, i3, y_c(i1+0,i2,i3)', i1, i2, i3, y_c(i1+0,i2,i3)
              y_c(i1+1,i2,i3)=dyi1!+cprecr*x_c(i1+1,i2,i3)
              y_c(i1+2,i2,i3)=dyi2!+cprecr*x_c(i1+2,i2,i3)
              y_c(i1+3,i2,i3)=dyi3!+cprecr*x_c(i1+3,i2,i3)
              x2_c(i1+0,i2,i3)=tt0!+cprecr*x_c(i1+0,i2,i3)
              x2_c(i1+1,i2,i3)=tt1!+cprecr*x_c(i1+1,i2,i3)
              x2_c(i1+2,i2,i3)=tt2!+cprecr*x_c(i1+2,i2,i3)
              x2_c(i1+3,i2,i3)=tt3!+cprecr*x_c(i1+3,i2,i3)
           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0.0_wp 
           tt0=0.0_wp
           ! Get the effective a-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*aeff0(t-i1)
              tt0=tt0 + x_c(t,i2,i3)*aeff0_2(t-i1)
           enddo
           y_c(i1,i2,i3)=dyi!+cprecr*x_c(i1,i2,i3)
           x2_c(i1,i2,i3)=tt0!+cprecr*x_c(i1,i2,i3)
        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective b-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, beff3(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, beff3_2(lowfil), 'b')
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + x_f1(t,i2,i3)*beff0(t-i1-0)
                 dyi1=dyi1 + x_f1(t,i2,i3)*beff1(t-i1-1)
                 dyi2=dyi2 + x_f1(t,i2,i3)*beff2(t-i1-2)
                 dyi3=dyi3 + x_f1(t,i2,i3)*beff3(t-i1-3)
                 tt0=tt0 + x_f1(t,i2,i3)*beff0_2(t-i1-0)
                 tt1=tt1 + x_f1(t,i2,i3)*beff1_2(t-i1-1)
                 tt2=tt2 + x_f1(t,i2,i3)*beff2_2(t-i1-2)
                 tt3=tt3 + x_f1(t,i2,i3)*beff3_2(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=y_c(i1+0,i2,i3)+dyi0
              y_c(i1+1,i2,i3)=y_c(i1+1,i2,i3)+dyi1
              y_c(i1+2,i2,i3)=y_c(i1+2,i2,i3)+dyi2
              y_c(i1+3,i2,i3)=y_c(i1+3,i2,i3)+dyi3
              x2_c(i1+0,i2,i3)=x2_c(i1+0,i2,i3)+tt0
              x2_c(i1+1,i2,i3)=x2_c(i1+1,i2,i3)+tt1
              x2_c(i1+2,i2,i3)=x2_c(i1+2,i2,i3)+tt2
              x2_c(i1+3,i2,i3)=x2_c(i1+3,i2,i3)+tt3
           enddo
           istart=i1
        endif

        do i1=istart,iend
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective b-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + x_f1(t,i2,i3)*beff0(t-i1)
              tt0=tt0 + x_f1(t,i2,i3)*beff0_2(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           x2_c(i1,i2,i3)=x2_c(i1,i2,i3)+tt0
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective c-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, ceff3(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, ceff3_2(lowfil), 'c')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*ceff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*ceff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*ceff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*ceff3(t-i1-3)
                 tt0=tt0 + x_c(t,i2,i3)*ceff0_2(t-i1-0)
                 tt1=tt1 + x_c(t,i2,i3)*ceff1_2(t-i1-1)
                 tt2=tt2 + x_c(t,i2,i3)*ceff2_2(t-i1-2)
                 tt3=tt3 + x_c(t,i2,i3)*ceff3_2(t-i1-3)
              enddo
              y_f(1,i1+0,i2,i3)=dyi0
              y_f(1,i1+1,i2,i3)=dyi1
              y_f(1,i1+2,i2,i3)=dyi2
              y_f(1,i1+3,i2,i3)=dyi3
              x2_f(1,i1+0,i2,i3)=tt0
              x2_f(1,i1+1,i2,i3)=tt1
              x2_f(1,i1+2,i2,i3)=tt2
              x2_f(1,i1+3,i2,i3)=tt3
           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective c-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*ceff0(t-i1)
              tt0=tt0 + x_c(t,i2,i3)*ceff0_2(t-i1)
           enddo
           y_f(1,i1,i2,i3)=dyi
           x2_f(1,i1,i2,i3)=tt0
        enddo
     enddo
  enddo
  !!!$omp enddo
  
  !  call system_clock(ncount1,ncount_rate,ncount_max)
  !  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective a-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, aeff3_2(lowfil), 'a')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*aeff0(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*aeff1(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*aeff2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*aeff3(t-i2-3)
                 tt0=tt0 + x_c(i1,t,i3)*aeff0_2(t-i2-0)
                 tt1=tt1 + x_c(i1,t,i3)*aeff1_2(t-i2-1)
                 tt2=tt2 + x_c(i1,t,i3)*aeff2_2(t-i2-2)
                 tt3=tt3 + x_c(i1,t,i3)*aeff3_2(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
              y2_c(i1,i2+0,i3)=tt0
              y2_c(i1,i2+1,i3)=tt1
              y2_c(i1,i2+2,i3)=tt2
              y2_c(i1,i2+3,i3)=tt3
           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective a-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*aeff0(t-i2)
              tt0=tt0 + x_c(i1,t,i3)*aeff0_2(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           y2_c(i1,i2,i3)=tt0
        enddo
        istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
        iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i2=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective b-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, beff3(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, beff3_2(lowfil), 'b')
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + x_f2(t,i1,i3)*beff0(t-i2-0)
                 dyi1=dyi1 + x_f2(t,i1,i3)*beff1(t-i2-1)
                 dyi2=dyi2 + x_f2(t,i1,i3)*beff2(t-i2-2)
                 dyi3=dyi3 + x_f2(t,i1,i3)*beff3(t-i2-3)
                 tt0=tt0 + x_f2(t,i1,i3)*beff0_2(t-i2-0)
                 tt1=tt1 + x_f2(t,i1,i3)*beff1_2(t-i2-1)
                 tt2=tt2 + x_f2(t,i1,i3)*beff2_2(t-i2-2)
                 tt3=tt3 + x_f2(t,i1,i3)*beff3_2(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
              y2_c(i1,i2+0,i3)=y2_c(i1,i2+0,i3)+tt0
              y2_c(i1,i2+1,i3)=y2_c(i1,i2+1,i3)+tt1
              y2_c(i1,i2+2,i3)=y2_c(i1,i2+2,i3)+tt2
              y2_c(i1,i2+3,i3)=y2_c(i1,i2+3,i3)+tt3
           enddo
           istart=i2
        endif

        do i2=istart,iend
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective b-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + x_f2(t,i1,i3)*beff0(t-i2)
              tt0=tt0 + x_f2(t,i1,i3)*beff0_2(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           y2_c(i1,i2,i3)=y2_c(i1,i2,i3)+tt0
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective c-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, ceff3(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, ceff3_2(lowfil), 'c')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*ceff0(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*ceff1(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*ceff2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*ceff3(t-i2-3)
                 tt0=tt0 + x_c(i1,t,i3)*ceff0_2(t-i2-0)
                 tt1=tt1 + x_c(i1,t,i3)*ceff1_2(t-i2-1)
                 tt2=tt2 + x_c(i1,t,i3)*ceff2_2(t-i2-2)
                 tt3=tt3 + x_c(i1,t,i3)*ceff3_2(t-i2-3)
              enddo
              y_f(2,i1,i2+0,i3)=dyi0
              y_f(2,i1,i2+1,i3)=dyi1
              y_f(2,i1,i2+2,i3)=dyi2
              y_f(2,i1,i2+3,i3)=dyi3
              y2_f(2,i1,i2+0,i3)=tt0
              y2_f(2,i1,i2+1,i3)=tt1
              y2_f(2,i1,i2+2,i3)=tt2
              y2_f(2,i1,i2+3,i3)=tt3
           enddo
           icur=i2
        else
           icur=ibxz_f(1,i1,i3)
        endif

        do i2=icur,ibxz_f(2,i1,i3)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective c-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*ceff0(t-i2)
              tt0=tt0 + x_c(i1,t,i3)*ceff0_2(t-i2)
           enddo
           y_f(2,i1,i2,i3)=dyi
           y2_f(2,i1,i2,i3)=tt0
        enddo
     enddo
  enddo
  !!!$omp enddo


  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective a-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, aeff3_2(lowfil), 'a')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*aeff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*aeff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*aeff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*aeff3(t-i3-3)
                 tt0=tt0 + x_c(i1,i2,t)*aeff0_2(t-i3-0)
                 tt1=tt1 + x_c(i1,i2,t)*aeff1_2(t-i3-1)
                 tt2=tt2 + x_c(i1,i2,t)*aeff2_2(t-i3-2)
                 tt3=tt3 + x_c(i1,i2,t)*aeff3_2(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
              z2_c(i1,i2,i3+0)=tt0
              z2_c(i1,i2,i3+1)=tt1
              z2_c(i1,i2,i3+2)=tt2
              z2_c(i1,i2,i3+3)=tt3
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective a-filters for the y dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*aeff0(t-i3)
              tt0=tt0 + x_c(i1,i2,t)*aeff0_2(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           z2_c(i1,i2,i3)=tt0
        enddo
        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective b-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, beff3(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, beff3_2(lowfil), 'b')
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + x_f3(t,i1,i2)*beff0(t-i3-0)
                 dyi1=dyi1 + x_f3(t,i1,i2)*beff1(t-i3-1)
                 dyi2=dyi2 + x_f3(t,i1,i2)*beff2(t-i3-2)
                 dyi3=dyi3 + x_f3(t,i1,i2)*beff3(t-i3-3)
                 tt0=tt0 + x_f3(t,i1,i2)*beff0_2(t-i3-0)
                 tt1=tt1 + x_f3(t,i1,i2)*beff1_2(t-i3-1)
                 tt2=tt2 + x_f3(t,i1,i2)*beff2_2(t-i3-2)
                 tt3=tt3 + x_f3(t,i1,i2)*beff3_2(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
              z2_c(i1,i2,i3+0)=z2_c(i1,i2,i3+0)+tt0
              z2_c(i1,i2,i3+1)=z2_c(i1,i2,i3+1)+tt1
              z2_c(i1,i2,i3+2)=z2_c(i1,i2,i3+2)+tt2
              z2_c(i1,i2,i3+3)=z2_c(i1,i2,i3+3)+tt3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective b-filters for the y dimension
           !!z0=hgrid*(i3+0)-rxyzConf(3)
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + x_f3(t,i1,i2)*beff0(t-i3)
              tt0=tt0 + x_f3(t,i1,i2)*beff0_2(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           z2_c(i1,i2,i3)=z2_c(i1,i2,i3)+tt0
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective c-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, ceff3(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, ceff3_2(lowfil), 'c')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*ceff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*ceff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*ceff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*ceff3(t-i3-3)
                 tt0=tt0 + x_c(i1,i2,t)*ceff0_2(t-i3-0)
                 tt1=tt1 + x_c(i1,i2,t)*ceff1_2(t-i3-1)
                 tt2=tt2 + x_c(i1,i2,t)*ceff2_2(t-i3-2)
                 tt3=tt3 + x_c(i1,i2,t)*ceff3_2(t-i3-3)
              enddo
              y_f(4,i1,i2,i3+0)=dyi0
              y_f(4,i1,i2,i3+1)=dyi1
              y_f(4,i1,i2,i3+2)=dyi2
              y_f(4,i1,i2,i3+3)=dyi3
              z2_f(4,i1,i2,i3+0)=tt0
              z2_f(4,i1,i2,i3+1)=tt1
              z2_f(4,i1,i2,i3+2)=tt2
              z2_f(4,i1,i2,i3+3)=tt3
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective c-filters for the z dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid,  z0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid,  z0, ceff0_2(lowfil), 'c')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*ceff0(t-i3)
              tt0=tt0 + x_c(i1,i2,t)*ceff0_2(t-i3)
           enddo
           y_f(4,i1,i2,i3)=dyi
           z2_f(4,i1,i2,i3)=tt0
        enddo
     enddo
  enddo
  !!!$omp enddo


  
  !  call system_clock(ncount3,ncount_rate,ncount_max)
  !  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2

  !!!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
           ! Get the effective filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0_2(lowfil), 'e')
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*aeff0(l) + x_f(5,i1+l,i2,i3)*beff0(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*aeff0(l) + x_f(3,i1+l,i2,i3)*beff0(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*aeff0(l) + x_f(7,i1+l,i2,i3)*beff0(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*ceff0(l) + x_f(5,i1+l,i2,i3)*eeff0(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*ceff0(l) + x_f(3,i1+l,i2,i3)*eeff0(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*ceff0(l) + x_f(7,i1+l,i2,i3)*eeff0(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*eeff0(l)
              tt112=tt112 + x_f(4,i1+l,i2,i3)*aeff0_2(l) + x_f(5,i1+l,i2,i3)*beff0_2(l)
              tt121=tt121 + x_f(2,i1+l,i2,i3)*aeff0_2(l) + x_f(3,i1+l,i2,i3)*beff0_2(l)
              tt122=tt122 + x_f(6,i1+l,i2,i3)*aeff0_2(l) + x_f(7,i1+l,i2,i3)*beff0_2(l)
              tt212=tt212 + x_f(4,i1+l,i2,i3)*ceff0_2(l) + x_f(5,i1+l,i2,i3)*eeff0_2(l)
              tt221=tt221 + x_f(2,i1+l,i2,i3)*ceff0_2(l) + x_f(3,i1+l,i2,i3)*eeff0_2(l)
              tt222=tt222 + x_f(6,i1+l,i2,i3)*ceff0_2(l) + x_f(7,i1+l,i2,i3)*eeff0_2(l)
              tt211=tt211 + x_f(1,i1+l,i2,i3)*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112!+cprecr*x_f(4,i1,i2,i3)
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121!+cprecr*x_f(2,i1,i2,i3)
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
           y_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
           y_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
           x2_f(4,i1,i2,i3)=+t112!+cprecr*x_f(4,i1,i2,i3)
           x2_f(2,i1,i2,i3)=+t121!+cprecr*x_f(2,i1,i2,i3)
           x2_f(1,i1,i2,i3)=x2_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
           x2_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
           x2_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
           x2_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
           x2_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
  !!!$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
           ! Get the effective filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, eeff0_2(lowfil), 'e')
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*aeff0(l) + x_f(6,i1,i2+l,i3)*beff0(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*aeff0(l) + x_f(3,i1,i2+l,i3)*beff0(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*ceff0(l) + x_f(6,i1,i2+l,i3)*eeff0(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*aeff0(l) + x_f(7,i1,i2+l,i3)*beff0(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*ceff0(l) + x_f(3,i1,i2+l,i3)*eeff0(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*ceff0(l) + x_f(7,i1,i2+l,i3)*eeff0(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*eeff0(l)
              tt112=tt112 + x_f(4,i1,i2+l,i3)*aeff0_2(l) + x_f(6,i1,i2+l,i3)*beff0_2(l)
              tt211=tt211 + x_f(1,i1,i2+l,i3)*aeff0_2(l) + x_f(3,i1,i2+l,i3)*beff0_2(l)
              tt122=tt122 + x_f(4,i1,i2+l,i3)*ceff0_2(l) + x_f(6,i1,i2+l,i3)*eeff0_2(l)
              tt212=tt212 + x_f(5,i1,i2+l,i3)*aeff0_2(l) + x_f(7,i1,i2+l,i3)*beff0_2(l)
              tt221=tt221 + x_f(1,i1,i2+l,i3)*ceff0_2(l) + x_f(3,i1,i2+l,i3)*eeff0_2(l)
              tt222=tt222 + x_f(5,i1,i2+l,i3)*ceff0_2(l) + x_f(7,i1,i2+l,i3)*eeff0_2(l)
              tt121=tt121 + x_f(2,i1,i2+l,i3)*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           y2_f(4,i1,i2,i3)=tt112
           y2_f(2,i1,i2,i3)=y2_f(2,i1,i2,i3)+tt121
           y2_f(1,i1,i2,i3)=tt211
           y2_f(6,i1,i2,i3)=tt122
           y2_f(5,i1,i2,i3)=tt212
           y2_f(3,i1,i2,i3)=tt221
           y2_f(7,i1,i2,i3)=tt222
        enddo
     enddo
  enddo
  !!!$omp enddo

  !  call system_clock(ncount5,ncount_rate,ncount_max)
  !  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
           ! Get the effective filters for the z dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           !call getFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, eeff0_2(lowfil), 'e')
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*aeff0(l) + x_f(6,i1,i2,i3+l)*beff0(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*aeff0(l) + x_f(5,i1,i2,i3+l)*beff0(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*ceff0(l) + x_f(6,i1,i2,i3+l)*eeff0(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*ceff0(l) + x_f(5,i1,i2,i3+l)*eeff0(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*aeff0(l) + x_f(7,i1,i2,i3+l)*beff0(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*ceff0(l) + x_f(7,i1,i2,i3+l)*eeff0(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*eeff0(l)
              tt121=tt121 + x_f(2,i1,i2,i3+l)*aeff0_2(l) + x_f(6,i1,i2,i3+l)*beff0_2(l)
              tt211=tt211 + x_f(1,i1,i2,i3+l)*aeff0_2(l) + x_f(5,i1,i2,i3+l)*beff0_2(l)
              tt122=tt122 + x_f(2,i1,i2,i3+l)*ceff0_2(l) + x_f(6,i1,i2,i3+l)*eeff0_2(l)
              tt212=tt212 + x_f(1,i1,i2,i3+l)*ceff0_2(l) + x_f(5,i1,i2,i3+l)*eeff0_2(l)
              tt221=tt221 + x_f(3,i1,i2,i3+l)*aeff0_2(l) + x_f(7,i1,i2,i3+l)*beff0_2(l)
              tt222=tt222 + x_f(3,i1,i2,i3+l)*ceff0_2(l) + x_f(7,i1,i2,i3+l)*eeff0_2(l)
              tt112=tt112 + x_f(4,i1,i2,i3+l)*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           z2_f(4,i1,i2,i3)=z2_f(4,i1,i2,i3)+tt112
           z2_f(2,i1,i2,i3)=tt121
           z2_f(1,i1,i2,i3)=tt211
           z2_f(6,i1,i2,i3)=tt122
           z2_f(5,i1,i2,i3)=tt212
           z2_f(3,i1,i2,i3)=tt221
           z2_f(7,i1,i2,i3)=tt222

        enddo
     enddo
  enddo
  !!!$omp enddo



  ! Sum up all the contributions
  do i3=0,n3
      do i2=0,n2
          do i1=0,n1
              tt0 = x2_c(i1,i2,i3)*y2_c(i1,i2,i3) + x2_c(i1,i2,i3)*z2_c(i1,i2,i3) + y2_c(i1,i2,i3)*z2_c(i1,i2,i3)
              y_c(i1,i2,i3) = y_c(i1,i2,i3) + 2.d0*tt0
          end do
      end do
  end do
  do i3=nfl3,nfu3
      do i2=nfl2,nfu2
          do i1=nfl1,nfu1
              tt1 = x2_f(1,i1,i2,i3)*y2_f(1,i1,i2,i3) + x2_f(1,i1,i2,i3)*z2_f(1,i1,i2,i3) + y2_f(1,i1,i2,i3)*z2_f(1,i1,i2,i3)
              tt2 = x2_f(2,i1,i2,i3)*y2_f(2,i1,i2,i3) + x2_f(2,i1,i2,i3)*z2_f(2,i1,i2,i3) + y2_f(2,i1,i2,i3)*z2_f(2,i1,i2,i3)
              tt3 = x2_f(3,i1,i2,i3)*y2_f(3,i1,i2,i3) + x2_f(3,i1,i2,i3)*z2_f(3,i1,i2,i3) + y2_f(3,i1,i2,i3)*z2_f(3,i1,i2,i3)
              tt4 = x2_f(4,i1,i2,i3)*y2_f(4,i1,i2,i3) + x2_f(4,i1,i2,i3)*z2_f(4,i1,i2,i3) + y2_f(4,i1,i2,i3)*z2_f(4,i1,i2,i3)
              tt5 = x2_f(5,i1,i2,i3)*y2_f(5,i1,i2,i3) + x2_f(5,i1,i2,i3)*z2_f(5,i1,i2,i3) + y2_f(5,i1,i2,i3)*z2_f(5,i1,i2,i3)
              tt6 = x2_f(6,i1,i2,i3)*y2_f(6,i1,i2,i3) + x2_f(6,i1,i2,i3)*z2_f(6,i1,i2,i3) + y2_f(6,i1,i2,i3)*z2_f(6,i1,i2,i3)
              tt7 = x2_f(7,i1,i2,i3)*y2_f(7,i1,i2,i3) + x2_f(7,i1,i2,i3)*z2_f(7,i1,i2,i3) + y2_f(7,i1,i2,i3)*z2_f(7,i1,i2,i3)
              y_f(1,i1,i2,i3) = y_f(1,i1,i2,i3) + 2.d0*tt1
              y_f(2,i1,i2,i3) = y_f(2,i1,i2,i3) + 2.d0*tt2
              y_f(3,i1,i2,i3) = y_f(3,i1,i2,i3) + 2.d0*tt3
              y_f(4,i1,i2,i3) = y_f(4,i1,i2,i3) + 2.d0*tt4
              y_f(5,i1,i2,i3) = y_f(5,i1,i2,i3) + 2.d0*tt5
              y_f(6,i1,i2,i3) = y_f(6,i1,i2,i3) + 2.d0*tt6
              y_f(7,i1,i2,i3) = y_f(7,i1,i2,i3) + 2.d0*tt7
          end do
      end do
  end do
  !!! Sum up all the contributions
  !!do i3=0,n3
  !!    do i2=0,n2
  !!        do i1=0,n1
  !!            y_c(i1,i2,i3) = x2_c(i1,i2,i3) + y2_c(i1,i2,i3) + z2_c(i1,i2,i3)
  !!        end do
  !!    end do
  !!end do
  !!do i3=nfl3,nfu3
  !!    do i2=nfl2,nfu2
  !!        do i1=nfl1,nfu1
  !!            y_f(1,i1,i2,i3) = x2_f(1,i1,i2,i3) + y2_f(1,i1,i2,i3) + z2_f(1,i1,i2,i3)
  !!            y_f(2,i1,i2,i3) = x2_f(2,i1,i2,i3) + y2_f(2,i1,i2,i3) + z2_f(2,i1,i2,i3)
  !!            y_f(3,i1,i2,i3) = x2_f(3,i1,i2,i3) + y2_f(3,i1,i2,i3) + z2_f(3,i1,i2,i3)
  !!            y_f(4,i1,i2,i3) = x2_f(4,i1,i2,i3) + y2_f(4,i1,i2,i3) + z2_f(4,i1,i2,i3)
  !!            y_f(5,i1,i2,i3) = x2_f(5,i1,i2,i3) + y2_f(5,i1,i2,i3) + z2_f(5,i1,i2,i3)
  !!            y_f(6,i1,i2,i3) = x2_f(6,i1,i2,i3) + y2_f(6,i1,i2,i3) + z2_f(6,i1,i2,i3)
  !!            y_f(7,i1,i2,i3) = x2_f(7,i1,i2,i3) + y2_f(7,i1,i2,i3) + z2_f(7,i1,i2,i3)
  !!        end do
  !!    end do
  !!end do
  


  iall=-product(shape(x2_c))*kind(x2_c)
  deallocate(x2_c, stat=istat)
  call memocc(istat, iall, 'x2_c', subname)

  iall=-product(shape(y2_c))*kind(y2_c)
  deallocate(y2_c, stat=istat)
  call memocc(istat, iall, 'y2_c', subname)

  iall=-product(shape(z2_c))*kind(z2_c)
  deallocate(z2_c, stat=istat)
  call memocc(istat, iall, 'z2_c', subname)

  iall=-product(shape(x2_f))*kind(x2_f)
  deallocate(x2_f, stat=istat)
  call memocc(istat, iall, 'x2_f', subname)

  iall=-product(shape(y2_f))*kind(y2_f)
  deallocate(y2_f, stat=istat)
  call memocc(istat, iall, 'y2_f', subname)

  iall=-product(shape(z2_f))*kind(z2_f)
  deallocate(z2_f, stat=istat)
  call memocc(istat, iall, 'z2_f', subname)




!  !!!$omp end parallel
!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel




END SUBROUTINE ConvolQuartic





!>  Applies the following operation: 
!!  y = [((r-r0)^4)]*x
subroutine ConvolQuartic2(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid, offsetx, offsety, offsetz, &
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzConf, potentialPrefac, it)
  use module_base
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, offsetx, offsety, offsetz
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
real(8),dimension(3):: rxyzConf
real(8):: potentialPrefac
integer:: it
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp):: tt112, tt121, tt122, tt212, tt221, tt222, tt211
  real(wp):: tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7
  !real(kind=8) :: tel
  real(wp), dimension(-3+lowfil:lupfil+3) :: a, aeff0, aeff1, aeff2, aeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: b, beff0, beff1, beff2, beff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: c, ceff0, ceff1, ceff2, ceff3
  real(wp), dimension(lowfil:lupfil) :: e, eeff0, eeff1, eeff2, eeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: aeff0_2, aeff1_2, aeff2_2, aeff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: beff0_2, beff1_2, beff2_2, beff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: ceff0_2, ceff1_2, ceff2_2, ceff3_2
  real(wp), dimension(lowfil:lupfil) :: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  real(wp),dimension(:,:,:),allocatable:: x2_c, y2_c, z2_c, x2_f2, x2_f3, y2_f3
  real(wp),dimension(:,:,:,:),allocatable:: x2_f, y2_f, z2_f
real(8):: x0, y0, z0
real(8):: x1, y1, z1
real(8):: x2, y2, z2
real(8):: x3, y3, z3
integer:: ii, istat, iall
character(len=*),parameter:: subname='ConvolQuartic2'


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
!
!
!  !---------------------------------------------------------------------------

  ! Scaling function part

!  call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)





allocate(x2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, x2_c, 'x2_c', subname)
allocate(y2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, y2_c, 'y2_c', subname)
allocate(z2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, z2_c, 'z2_c', subname)
allocate(x2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, x2_f, 'x2_f', subname)
allocate(y2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, y2_f, 'y2_f', subname)
allocate(z2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, z2_f, 'z2_f', subname)

allocate(x2_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2_f2, 'x2_f2', subname)
allocate(x2_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2_f3, 'x2_f3', subname)
allocate(y2_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y2_f3, 'y2_f3', subname)



x2_c=0.d0
y2_c=0.d0
z2_c=0.d0
x2_f=0.d0
y2_f=0.d0
z2_f=0.d0
x2_f2=0.d0
x2_f3=0.d0
y2_f3=0.d0


aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0

aeff0_2=0.d0 ; beff0_2=0.d0 ; ceff0_2=0.d0 ; eeff0_2=0.0
aeff1_2=0.d0 ; beff1_2=0.d0 ; ceff1_2=0.d0 ; eeff1_2=0.0
aeff2_2=0.d0 ; beff2_2=0.d0 ; ceff2_2=0.d0 ; eeff2_2=0.0
aeff3_2=0.d0 ; beff3_2=0.d0 ; ceff3_2=0.d0 ; eeff3_2=0.0


!!$!$omp parallel default(private) &
!!$!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!!$!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!!$!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
  !!!$omp do  
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective a-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, aeff3_2(lowfil), 'a')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*aeff0(t-i1-0)
                 !!write(3100,'(a,3i9,2es16.7,7i10)') 't, i2, i3, x_c(t,i2,i3), aeff0(t-i1-0), n1, n2, n3, ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)-4, max(ibyz_c(1,i2,i3),lowfil+i1), min(lupfil+i1+3,ibyz_c(2,i2,i3))', &
                 !!    t, i2, i3, x_c(t,i2,i3), aeff0(t-i1-0), n1, n2, n3, ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)-4, max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi1=dyi1 + x_c(t,i2,i3)*aeff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*aeff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*aeff3(t-i1-3)
                 tt0=tt0 + x_c(t,i2,i3)*aeff0_2(t-i1-0)
                 tt1=tt1 + x_c(t,i2,i3)*aeff1_2(t-i1-1)
                 tt2=tt2 + x_c(t,i2,i3)*aeff2_2(t-i1-2)
                 tt3=tt3 + x_c(t,i2,i3)*aeff3_2(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=dyi0!+cprecr*x_c(i1+0,i2,i3)
              !!write(3000,'(a,3i9,es16.7)') 'i1, i2, i3, y_c(i1+0,i2,i3)', i1, i2, i3, y_c(i1+0,i2,i3)
              y_c(i1+1,i2,i3)=dyi1!+cprecr*x_c(i1+1,i2,i3)
              y_c(i1+2,i2,i3)=dyi2!+cprecr*x_c(i1+2,i2,i3)
              y_c(i1+3,i2,i3)=dyi3!+cprecr*x_c(i1+3,i2,i3)
              x2_c(i1+0,i2,i3)=tt0!+cprecr*x_c(i1+0,i2,i3)
              x2_c(i1+1,i2,i3)=tt1!+cprecr*x_c(i1+1,i2,i3)
              x2_c(i1+2,i2,i3)=tt2!+cprecr*x_c(i1+2,i2,i3)
              x2_c(i1+3,i2,i3)=tt3!+cprecr*x_c(i1+3,i2,i3)
           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0.0_wp 
           tt0=0.0_wp
           ! Get the effective a-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*aeff0(t-i1)
              tt0=tt0 + x_c(t,i2,i3)*aeff0_2(t-i1)
           enddo
           y_c(i1,i2,i3)=dyi!+cprecr*x_c(i1,i2,i3)
           x2_c(i1,i2,i3)=tt0!+cprecr*x_c(i1,i2,i3)
        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective b-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, beff3(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, beff3_2(lowfil), 'b')
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + x_f1(t,i2,i3)*beff0(t-i1-0)
                 dyi1=dyi1 + x_f1(t,i2,i3)*beff1(t-i1-1)
                 dyi2=dyi2 + x_f1(t,i2,i3)*beff2(t-i1-2)
                 dyi3=dyi3 + x_f1(t,i2,i3)*beff3(t-i1-3)
                 tt0=tt0 + x_f1(t,i2,i3)*beff0_2(t-i1-0)
                 tt1=tt1 + x_f1(t,i2,i3)*beff1_2(t-i1-1)
                 tt2=tt2 + x_f1(t,i2,i3)*beff2_2(t-i1-2)
                 tt3=tt3 + x_f1(t,i2,i3)*beff3_2(t-i1-3)
              enddo
              y_c(i1+0,i2,i3)=y_c(i1+0,i2,i3)+dyi0
              y_c(i1+1,i2,i3)=y_c(i1+1,i2,i3)+dyi1
              y_c(i1+2,i2,i3)=y_c(i1+2,i2,i3)+dyi2
              y_c(i1+3,i2,i3)=y_c(i1+3,i2,i3)+dyi3
              x2_c(i1+0,i2,i3)=x2_c(i1+0,i2,i3)+tt0
              x2_c(i1+1,i2,i3)=x2_c(i1+1,i2,i3)+tt1
              x2_c(i1+2,i2,i3)=x2_c(i1+2,i2,i3)+tt2
              x2_c(i1+3,i2,i3)=x2_c(i1+3,i2,i3)+tt3
           enddo
           istart=i1
        endif

        do i1=istart,iend
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective b-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + x_f1(t,i2,i3)*beff0(t-i1)
              tt0=tt0 + x_f1(t,i2,i3)*beff0_2(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           x2_c(i1,i2,i3)=x2_c(i1,i2,i3)+tt0
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective c-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, ceff3(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, ceff3_2(lowfil), 'c')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*ceff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*ceff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*ceff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*ceff3(t-i1-3)
                 tt0=tt0 + x_c(t,i2,i3)*ceff0_2(t-i1-0)
                 tt1=tt1 + x_c(t,i2,i3)*ceff1_2(t-i1-1)
                 tt2=tt2 + x_c(t,i2,i3)*ceff2_2(t-i1-2)
                 tt3=tt3 + x_c(t,i2,i3)*ceff3_2(t-i1-3)
              enddo
              y_f(1,i1+0,i2,i3)=dyi0
              y_f(1,i1+1,i2,i3)=dyi1
              y_f(1,i1+2,i2,i3)=dyi2
              y_f(1,i1+3,i2,i3)=dyi3
              x2_f(1,i1+0,i2,i3)=tt0
              x2_f(1,i1+1,i2,i3)=tt1
              x2_f(1,i1+2,i2,i3)=tt2
              x2_f(1,i1+3,i2,i3)=tt3
           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective c-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*ceff0(t-i1)
              tt0=tt0 + x_c(t,i2,i3)*ceff0_2(t-i1)
           enddo
           y_f(1,i1,i2,i3)=dyi
           x2_f(1,i1,i2,i3)=tt0
        enddo
     enddo
  enddo
  !!!$omp enddo





  ! wavelet part

  !!!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
           ! Get the effective filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0_2(lowfil), 'e')
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*aeff0(l) + x_f(5,i1+l,i2,i3)*beff0(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*aeff0(l) + x_f(3,i1+l,i2,i3)*beff0(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*aeff0(l) + x_f(7,i1+l,i2,i3)*beff0(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*ceff0(l) + x_f(5,i1+l,i2,i3)*eeff0(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*ceff0(l) + x_f(3,i1+l,i2,i3)*eeff0(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*ceff0(l) + x_f(7,i1+l,i2,i3)*eeff0(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*eeff0(l)
              tt112=tt112 + x_f(4,i1+l,i2,i3)*aeff0_2(l) + x_f(5,i1+l,i2,i3)*beff0_2(l)
              tt121=tt121 + x_f(2,i1+l,i2,i3)*aeff0_2(l) + x_f(3,i1+l,i2,i3)*beff0_2(l)
              tt122=tt122 + x_f(6,i1+l,i2,i3)*aeff0_2(l) + x_f(7,i1+l,i2,i3)*beff0_2(l)
              tt212=tt212 + x_f(4,i1+l,i2,i3)*ceff0_2(l) + x_f(5,i1+l,i2,i3)*eeff0_2(l)
              tt221=tt221 + x_f(2,i1+l,i2,i3)*ceff0_2(l) + x_f(3,i1+l,i2,i3)*eeff0_2(l)
              tt222=tt222 + x_f(6,i1+l,i2,i3)*ceff0_2(l) + x_f(7,i1+l,i2,i3)*eeff0_2(l)
              tt211=tt211 + x_f(1,i1+l,i2,i3)*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=t112!+cprecr*x_f(4,i1,i2,i3)
           y_f(2,i1,i2,i3)=t121!+cprecr*x_f(2,i1,i2,i3)
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
           y_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
           y_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
           x2_f(4,i1,i2,i3)=+t112!+cprecr*x_f(4,i1,i2,i3)
           x2_f(2,i1,i2,i3)=+t121!+cprecr*x_f(2,i1,i2,i3)
           x2_f(1,i1,i2,i3)=x2_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
           x2_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
           x2_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
           x2_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
           x2_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
  !!!$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel



  ! copy to auxiliary array for better memory access
  do i3=nfl3,nfu3
      do i2=nfl2,nfu2
          do i1=nfl1,nfu1
              !x2_f1(i1,i2,i3)=x2_f(1,i1,i2,i3)
              x2_f2(i2,i1,i3)=x2_f(2,i1,i2,i3)
              x2_f3(i3,i1,i2)=x2_f(4,i1,i2,i3)
          end do
      end do
  end do







  
  !  call system_clock(ncount1,ncount_rate,ncount_max)
  !  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective a-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, aeff3_2(lowfil), 'a')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*aeff0(t-i2-0) + 2.d0*x2_c(i1,t,i3)*aeff0_2(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*aeff1(t-i2-1) + 2.d0*x2_c(i1,t,i3)*aeff1_2(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*aeff2(t-i2-2) + 2.d0*x2_c(i1,t,i3)*aeff2_2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*aeff3(t-i2-3) + 2.d0*x2_c(i1,t,i3)*aeff3_2(t-i2-3)
                 tt0=tt0 + x_c(i1,t,i3)*aeff0_2(t-i2-0)
                 tt1=tt1 + x_c(i1,t,i3)*aeff1_2(t-i2-1)
                 tt2=tt2 + x_c(i1,t,i3)*aeff2_2(t-i2-2)
                 tt3=tt3 + x_c(i1,t,i3)*aeff3_2(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
              y2_c(i1,i2+0,i3)=tt0
              y2_c(i1,i2+1,i3)=tt1
              y2_c(i1,i2+2,i3)=tt2
              y2_c(i1,i2+3,i3)=tt3
           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective a-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*aeff0(t-i2) + 2.d0*x2_c(i1,t,i3)*aeff0_2(t-i2)
              tt0=tt0 + x_c(i1,t,i3)*aeff0_2(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           y2_c(i1,i2,i3)=tt0
        enddo
        istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
        iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i2=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective b-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, beff3(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, beff3_2(lowfil), 'b')
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + x_f2(t,i1,i3)*beff0(t-i2-0) + 2.d0*x2_f2(t,i1,i3)*beff0_2(t-i2-0)
                 dyi1=dyi1 + x_f2(t,i1,i3)*beff1(t-i2-1) + 2.d0*x2_f2(t,i1,i3)*beff1_2(t-i2-1)
                 dyi2=dyi2 + x_f2(t,i1,i3)*beff2(t-i2-2) + 2.d0*x2_f2(t,i1,i3)*beff2_2(t-i2-2)
                 dyi3=dyi3 + x_f2(t,i1,i3)*beff3(t-i2-3) + 2.d0*x2_f2(t,i1,i3)*beff3_2(t-i2-3)
                 tt0=tt0 + x_f2(t,i1,i3)*beff0_2(t-i2-0)
                 tt1=tt1 + x_f2(t,i1,i3)*beff1_2(t-i2-1)
                 tt2=tt2 + x_f2(t,i1,i3)*beff2_2(t-i2-2)
                 tt3=tt3 + x_f2(t,i1,i3)*beff3_2(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
              y2_c(i1,i2+0,i3)=y2_c(i1,i2+0,i3)+tt0
              y2_c(i1,i2+1,i3)=y2_c(i1,i2+1,i3)+tt1
              y2_c(i1,i2+2,i3)=y2_c(i1,i2+2,i3)+tt2
              y2_c(i1,i2+3,i3)=y2_c(i1,i2+3,i3)+tt3
           enddo
           istart=i2
        endif

        do i2=istart,iend
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective b-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + x_f2(t,i1,i3)*beff0(t-i2) + 2.d0*x2_f2(t,i1,i3)*beff0_2(t-i2)
              tt0=tt0 + x_f2(t,i1,i3)*beff0_2(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           y2_c(i1,i2,i3)=y2_c(i1,i2,i3)+tt0
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective c-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, ceff3(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, ceff3_2(lowfil), 'c')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*ceff0(t-i2-0) + 2.d0*x2_c(i1,t,i3)*ceff0_2(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*ceff1(t-i2-1) + 2.d0*x2_c(i1,t,i3)*ceff1_2(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*ceff2(t-i2-2) + 2.d0*x2_c(i1,t,i3)*ceff2_2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*ceff3(t-i2-3) + 2.d0*x2_c(i1,t,i3)*ceff3_2(t-i2-3)
                 tt0=tt0 + x_c(i1,t,i3)*ceff0_2(t-i2-0)
                 tt1=tt1 + x_c(i1,t,i3)*ceff1_2(t-i2-1)
                 tt2=tt2 + x_c(i1,t,i3)*ceff2_2(t-i2-2)
                 tt3=tt3 + x_c(i1,t,i3)*ceff3_2(t-i2-3)
              enddo
              y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
              y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+dyi1
              y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+dyi2
              y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+dyi3
              y2_f(2,i1,i2+0,i3)=tt0
              y2_f(2,i1,i2+1,i3)=tt1
              y2_f(2,i1,i2+2,i3)=tt2
              y2_f(2,i1,i2+3,i3)=tt3
           enddo
           icur=i2
        else
           icur=ibxz_f(1,i1,i3)
        endif

        do i2=icur,ibxz_f(2,i1,i3)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective c-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*ceff0(t-i2) + 2.d0*x2_c(i1,t,i3)*ceff0_2(t-i2)
              tt0=tt0 + x_c(i1,t,i3)*ceff0_2(t-i2)
           enddo
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+dyi
           y2_f(2,i1,i2,i3)=tt0
        enddo
     enddo
  enddo
  !!!$omp enddo


  ! wavelet part

  !!!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
           ! Get the effective filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, eeff0_2(lowfil), 'e')
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*aeff0(l) + 2.d0*x2_f(4,i1,i2+l,i3)*aeff0_2(l) + &
                          x_f(6,i1,i2+l,i3)*beff0(l) + 2.d0*x2_f(6,i1,i2+l,i3)*beff0_2(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*aeff0(l) + 2.d0*x2_f(1,i1,i2+l,i3)*aeff0_2(l) + &
                          x_f(3,i1,i2+l,i3)*beff0(l) + 2.d0*x2_f(3,i1,i2+l,i3)*beff0_2(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*ceff0(l) + 2.d0*x2_f(4,i1,i2+l,i3)*ceff0_2(l) + &
                          x_f(6,i1,i2+l,i3)*eeff0(l) + 2.d0*x2_f(6,i1,i2+l,i3)*eeff0_2(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*aeff0(l) + 2.d0*x2_f(5,i1,i2+l,i3)*aeff0_2(l) + &
                          x_f(7,i1,i2+l,i3)*beff0(l) + 2.d0*x2_f(7,i1,i2+l,i3)*beff0_2(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*ceff0(l) + 2.d0*x2_f(1,i1,i2+l,i3)*ceff0_2(l) + &
                          x_f(3,i1,i2+l,i3)*eeff0(l) + 2.d0*x2_f(3,i1,i2+l,i3)*eeff0_2(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*ceff0(l) + 2.d0*x2_f(5,i1,i2+l,i3)*ceff0_2(l) + &
                          x_f(7,i1,i2+l,i3)*eeff0(l) + 2.d0*x2_f(7,i1,i2+l,i3)*eeff0_2(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*eeff0(l) + 2.d0*x2_f(2,i1,i2+l,i3)*eeff0_2(l)
              tt112=tt112 + x_f(4,i1,i2+l,i3)*aeff0_2(l) + x_f(6,i1,i2+l,i3)*beff0_2(l)
              tt211=tt211 + x_f(1,i1,i2+l,i3)*aeff0_2(l) + x_f(3,i1,i2+l,i3)*beff0_2(l)
              tt122=tt122 + x_f(4,i1,i2+l,i3)*ceff0_2(l) + x_f(6,i1,i2+l,i3)*eeff0_2(l)
              tt212=tt212 + x_f(5,i1,i2+l,i3)*aeff0_2(l) + x_f(7,i1,i2+l,i3)*beff0_2(l)
              tt221=tt221 + x_f(1,i1,i2+l,i3)*ceff0_2(l) + x_f(3,i1,i2+l,i3)*eeff0_2(l)
              tt222=tt222 + x_f(5,i1,i2+l,i3)*ceff0_2(l) + x_f(7,i1,i2+l,i3)*eeff0_2(l)
              tt121=tt121 + x_f(2,i1,i2+l,i3)*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           y2_f(4,i1,i2,i3)=tt112
           y2_f(2,i1,i2,i3)=y2_f(2,i1,i2,i3)+tt121
           y2_f(1,i1,i2,i3)=tt211
           y2_f(6,i1,i2,i3)=tt122
           y2_f(5,i1,i2,i3)=tt212
           y2_f(3,i1,i2,i3)=tt221
           y2_f(7,i1,i2,i3)=tt222
        enddo
     enddo
  enddo
  !!!$omp enddo




  ! copy to auxiliary array for better memory access
  do i3=nfl3,nfu3
      do i2=nfl2,nfu2
          do i1=nfl1,nfu1
              !y2_f1(i1,i2,i3)=y2_f(1,i1,i2,i3)
              !y2_f2(i2,i1,i3)=y2_f(2,i1,i2,i3)
              y2_f3(i3,i1,i2)=y2_f(4,i1,i2,i3)
          end do
      end do
  end do




  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective a-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, aeff3_2(lowfil), 'a')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*aeff0(t-i3-0) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*aeff0_2(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*aeff1(t-i3-1) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*aeff1_2(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*aeff2(t-i3-2) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*aeff2_2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*aeff3(t-i3-3) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*aeff3_2(t-i3-3)
                 !!tt0=tt0 + x_c(i1,i2,t)*aeff0_2(t-i3-0)
                 !!tt1=tt1 + x_c(i1,i2,t)*aeff1_2(t-i3-1)
                 !!tt2=tt2 + x_c(i1,i2,t)*aeff2_2(t-i3-2)
                 !!tt3=tt3 + x_c(i1,i2,t)*aeff3_2(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
              !!z2_c(i1,i2,i3+0)=tt0
              !!z2_c(i1,i2,i3+1)=tt1
              !!z2_c(i1,i2,i3+2)=tt2
              !!z2_c(i1,i2,i3+3)=tt3
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective a-filters for the y dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*aeff0(t-i3) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*aeff0_2(t-i3)
              !!tt0=tt0 + x_c(i1,i2,t)*aeff0_2(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           !!z2_c(i1,i2,i3)=tt0
        enddo
        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective b-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, beff3(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, beff3_2(lowfil), 'b')
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + x_f3(t,i1,i2)*beff0(t-i3-0) + 2.d0*(x2_f3(t,i1,i2)+y2_f3(t,i1,i2))*beff0_2(t-i3-0)
                 dyi1=dyi1 + x_f3(t,i1,i2)*beff1(t-i3-1) + 2.d0*(x2_f3(t,i1,i2)+y2_f3(t,i1,i2))*beff1_2(t-i3-1)
                 dyi2=dyi2 + x_f3(t,i1,i2)*beff2(t-i3-2) + 2.d0*(x2_f3(t,i1,i2)+y2_f3(t,i1,i2))*beff2_2(t-i3-2)
                 dyi3=dyi3 + x_f3(t,i1,i2)*beff3(t-i3-3) + 2.d0*(x2_f3(t,i1,i2)+y2_f3(t,i1,i2))*beff3_2(t-i3-3)
                 !!tt0=tt0 + x_f3(t,i1,i2)*beff0_2(t-i3-0)
                 !!tt1=tt1 + x_f3(t,i1,i2)*beff1_2(t-i3-1)
                 !!tt2=tt2 + x_f3(t,i1,i2)*beff2_2(t-i3-2)
                 !!tt3=tt3 + x_f3(t,i1,i2)*beff3_2(t-i3-3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
              !!z2_c(i1,i2,i3+0)=z2_c(i1,i2,i3+0)+tt0
              !!z2_c(i1,i2,i3+1)=z2_c(i1,i2,i3+1)+tt1
              !!z2_c(i1,i2,i3+2)=z2_c(i1,i2,i3+2)+tt2
              !!z2_c(i1,i2,i3+3)=z2_c(i1,i2,i3+3)+tt3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi=0.0_wp
           tt0=0.0_wp
           ! Get the effective b-filters for the y dimension
           !!z0=hgrid*(i3+0)-rxyzConf(3)
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + x_f3(t,i1,i2)*beff0(t-i3) + 2.d0*(x2_f3(t,i1,i2)+y2_f3(t,i1,i2))*beff0_2(t-i3)
              !!tt0=tt0 + x_f3(t,i1,i2)*beff0_2(t-i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
           !!z2_c(i1,i2,i3)=z2_c(i1,i2,i3)+tt0
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp 
              tt1=0.0_wp 
              tt2=0.0_wp 
              tt3=0.0_wp 
              ! Get the effective c-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, ceff3(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, ceff3_2(lowfil), 'c')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*ceff0(t-i3-0) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*ceff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*ceff1(t-i3-1) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*ceff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*ceff2(t-i3-2) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*ceff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*ceff3(t-i3-3) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*ceff3(t-i3-3)
                 !!tt0=tt0 + x_c(i1,i2,t)*ceff0_2(t-i3-0)
                 !!tt1=tt1 + x_c(i1,i2,t)*ceff1_2(t-i3-1)
                 !!tt2=tt2 + x_c(i1,i2,t)*ceff2_2(t-i3-2)
                 !!tt3=tt3 + x_c(i1,i2,t)*ceff3_2(t-i3-3)
              enddo
              y_f(4,i1,i2,i3+0)=dyi0
              y_f(4,i1,i2,i3+1)=dyi1
              y_f(4,i1,i2,i3+2)=dyi2
              y_f(4,i1,i2,i3+3)=dyi3
              !!z2_f(4,i1,i2,i3+0)=tt0
              !!z2_f(4,i1,i2,i3+1)=tt1
              !!z2_f(4,i1,i2,i3+2)=tt2
              !!z2_f(4,i1,i2,i3+3)=tt3
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective c-filters for the z dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid,  z0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid,  z0, ceff0_2(lowfil), 'c')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + x_c(i1,i2,t)*ceff0(t-i3) + 2.d0*(x2_c(i1,i2,t)+y2_c(i1,i2,t))*ceff0(t-i3)
              !!tt0=tt0 + x_c(i1,i2,t)*ceff0_2(t-i3)
           enddo
           y_f(4,i1,i2,i3)=dyi
           !!z2_f(4,i1,i2,i3)=tt0
        enddo
     enddo
  enddo
  !!!$omp enddo


  ! wavelet part

  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
           ! Get the effective filters for the z dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           !call getFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, eeff0_2(lowfil), 'e')
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*aeff0(l) + 2.d0*(x2_f(2,i1,i2,i3+l)+y2_f(2,i1,i2,i3+l))*aeff0_2(l) + &
                          x_f(6,i1,i2,i3+l)*beff0(l) + 2.d0*(x2_f(6,i1,i2,i3+l)+x2_f(6,i1,i2,i3+l))*beff0_2(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*aeff0(l) + 2.d0*(x2_f(1,i1,i2,i3+l)+y2_f(1,i1,i2,i3+l))*aeff0_2(l) + &
                          x_f(5,i1,i2,i3+l)*beff0(l) + 2.d0*(x2_f(5,i1,i2,i3+l)+x2_f(5,i1,i2,i3+l))*beff0_2(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*ceff0(l) + 2.d0*(x2_f(2,i1,i2,i3+l)+y2_f(2,i1,i2,i3+l))*ceff0_2(l) + &
                          x_f(6,i1,i2,i3+l)*eeff0(l) + 2.d0*(x2_f(6,i1,i2,i3+l)+x2_f(6,i1,i2,i3+l))*eeff0_2(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*ceff0(l) + 2.d0*(x2_f(1,i1,i2,i3+l)+y2_f(1,i1,i2,i3+l))*ceff0_2(l) + &
                          x_f(5,i1,i2,i3+l)*eeff0(l) + 2.d0*(x2_f(5,i1,i2,i3+l)+x2_f(5,i1,i2,i3+l))*eeff0_2(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*aeff0(l) + 2.d0*(x2_f(3,i1,i2,i3+l)+y2_f(3,i1,i2,i3+l))*aeff0_2(l) + &
                          x_f(7,i1,i2,i3+l)*beff0(l) + 2.d0*(x2_f(7,i1,i2,i3+l)+x2_f(7,i1,i2,i3+l))*beff0_2(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*ceff0(l) + 2.d0*(x2_f(3,i1,i2,i3+l)+y2_f(3,i1,i2,i3+l))*ceff0_2(l) + &
                          x_f(7,i1,i2,i3+l)*eeff0(l) + 2.d0*(x2_f(7,i1,i2,i3+l)+x2_f(7,i1,i2,i3+l))*eeff0_2(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*eeff0(l) + 2.d0*(x2_f(4,i1,i2,i3+l)+y2_f(4,i1,i2,i3+l))*eeff0_2(l)
              !!tt121=tt121 + x_f(2,i1,i2,i3+l)*aeff0_2(l) + x_f(6,i1,i2,i3+l)*beff0_2(l)
              !!tt211=tt211 + x_f(1,i1,i2,i3+l)*aeff0_2(l) + x_f(5,i1,i2,i3+l)*beff0_2(l)
              !!tt122=tt122 + x_f(2,i1,i2,i3+l)*ceff0_2(l) + x_f(6,i1,i2,i3+l)*eeff0_2(l)
              !!tt212=tt212 + x_f(1,i1,i2,i3+l)*ceff0_2(l) + x_f(5,i1,i2,i3+l)*eeff0_2(l)
              !!tt221=tt221 + x_f(3,i1,i2,i3+l)*aeff0_2(l) + x_f(7,i1,i2,i3+l)*beff0_2(l)
              !!tt222=tt222 + x_f(3,i1,i2,i3+l)*ceff0_2(l) + x_f(7,i1,i2,i3+l)*eeff0_2(l)
              !!tt112=tt112 + x_f(4,i1,i2,i3+l)*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           !!z2_f(4,i1,i2,i3)=z2_f(4,i1,i2,i3)+tt112
           !!z2_f(2,i1,i2,i3)=tt121
           !!z2_f(1,i1,i2,i3)=tt211
           !!z2_f(6,i1,i2,i3)=tt122
           !!z2_f(5,i1,i2,i3)=tt212
           !!z2_f(3,i1,i2,i3)=tt221
           !!z2_f(7,i1,i2,i3)=tt222

        enddo
     enddo
  enddo
  !!!$omp enddo






  
!!!!!  !  call system_clock(ncount3,ncount_rate,ncount_max)
!!!!!  !  tel=dble(ncount3-ncount2)/dble(ncount_rate)
!!!!!  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel
!!!!!
!!!!!  ! wavelet part
!!!!!  ! (1/2) d^2/dx^2
!!!!!
!!!!!  !!!$omp do
!!!!!  do i3=nfl3,nfu3
!!!!!     do i2=nfl2,nfu2
!!!!!        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!!!!!           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
!!!!!           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
!!!!!           ! Get the effective filters for the x dimension
!!!!!           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, x0, eeff0(lowfil), 'e')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0_2(lowfil), 'e')
!!!!!           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
!!!!!              t112=t112 + x_f(4,i1+l,i2,i3)*aeff0(l) + x_f(5,i1+l,i2,i3)*beff0(l)
!!!!!              t121=t121 + x_f(2,i1+l,i2,i3)*aeff0(l) + x_f(3,i1+l,i2,i3)*beff0(l)
!!!!!              t122=t122 + x_f(6,i1+l,i2,i3)*aeff0(l) + x_f(7,i1+l,i2,i3)*beff0(l)
!!!!!              t212=t212 + x_f(4,i1+l,i2,i3)*ceff0(l) + x_f(5,i1+l,i2,i3)*eeff0(l)
!!!!!              t221=t221 + x_f(2,i1+l,i2,i3)*ceff0(l) + x_f(3,i1+l,i2,i3)*eeff0(l)
!!!!!              t222=t222 + x_f(6,i1+l,i2,i3)*ceff0(l) + x_f(7,i1+l,i2,i3)*eeff0(l)
!!!!!              t211=t211 + x_f(1,i1+l,i2,i3)*eeff0(l)
!!!!!              tt112=tt112 + x_f(4,i1+l,i2,i3)*aeff0_2(l) + x_f(5,i1+l,i2,i3)*beff0_2(l)
!!!!!              tt121=tt121 + x_f(2,i1+l,i2,i3)*aeff0_2(l) + x_f(3,i1+l,i2,i3)*beff0_2(l)
!!!!!              tt122=tt122 + x_f(6,i1+l,i2,i3)*aeff0_2(l) + x_f(7,i1+l,i2,i3)*beff0_2(l)
!!!!!              tt212=tt212 + x_f(4,i1+l,i2,i3)*ceff0_2(l) + x_f(5,i1+l,i2,i3)*eeff0_2(l)
!!!!!              tt221=tt221 + x_f(2,i1+l,i2,i3)*ceff0_2(l) + x_f(3,i1+l,i2,i3)*eeff0_2(l)
!!!!!              tt222=tt222 + x_f(6,i1+l,i2,i3)*ceff0_2(l) + x_f(7,i1+l,i2,i3)*eeff0_2(l)
!!!!!              tt211=tt211 + x_f(1,i1+l,i2,i3)*eeff0_2(l)
!!!!!           enddo
!!!!!           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112!+cprecr*x_f(4,i1,i2,i3)
!!!!!           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121!+cprecr*x_f(2,i1,i2,i3)
!!!!!           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
!!!!!           y_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
!!!!!           y_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
!!!!!           y_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
!!!!!           y_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
!!!!!           x2_f(4,i1,i2,i3)=+t112!+cprecr*x_f(4,i1,i2,i3)
!!!!!           x2_f(2,i1,i2,i3)=+t121!+cprecr*x_f(2,i1,i2,i3)
!!!!!           x2_f(1,i1,i2,i3)=x2_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
!!!!!           x2_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
!!!!!           x2_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
!!!!!           x2_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
!!!!!           x2_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
!!!!!        enddo
!!!!!     enddo
!!!!!  enddo
!!!!!  !!!$omp enddo
!!!!!
!!!!!  !  call system_clock(ncount4,ncount_rate,ncount_max)
!!!!!  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
!!!!!  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel
!!!!!
!!!!!
!!!!!  ! + (1/2) d^2/dy^2
!!!!!  !!!$omp do
!!!!!  do i3=nfl3,nfu3
!!!!!     do i1=nfl1,nfu1
!!!!!        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!!!!!           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
!!!!!           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
!!!!!           ! Get the effective filters for the y dimension
!!!!!           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, y0, eeff0(lowfil), 'e')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, eeff0_2(lowfil), 'e')
!!!!!           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
!!!!!              t112=t112 + x_f(4,i1,i2+l,i3)*aeff0(l) + x_f(6,i1,i2+l,i3)*beff0(l)
!!!!!              t211=t211 + x_f(1,i1,i2+l,i3)*aeff0(l) + x_f(3,i1,i2+l,i3)*beff0(l)
!!!!!              t122=t122 + x_f(4,i1,i2+l,i3)*ceff0(l) + x_f(6,i1,i2+l,i3)*eeff0(l)
!!!!!              t212=t212 + x_f(5,i1,i2+l,i3)*aeff0(l) + x_f(7,i1,i2+l,i3)*beff0(l)
!!!!!              t221=t221 + x_f(1,i1,i2+l,i3)*ceff0(l) + x_f(3,i1,i2+l,i3)*eeff0(l)
!!!!!              t222=t222 + x_f(5,i1,i2+l,i3)*ceff0(l) + x_f(7,i1,i2+l,i3)*eeff0(l)
!!!!!              t121=t121 + x_f(2,i1,i2+l,i3)*eeff0(l)
!!!!!              tt112=tt112 + x_f(4,i1,i2+l,i3)*aeff0_2(l) + x_f(6,i1,i2+l,i3)*beff0_2(l)
!!!!!              tt211=tt211 + x_f(1,i1,i2+l,i3)*aeff0_2(l) + x_f(3,i1,i2+l,i3)*beff0_2(l)
!!!!!              tt122=tt122 + x_f(4,i1,i2+l,i3)*ceff0_2(l) + x_f(6,i1,i2+l,i3)*eeff0_2(l)
!!!!!              tt212=tt212 + x_f(5,i1,i2+l,i3)*aeff0_2(l) + x_f(7,i1,i2+l,i3)*beff0_2(l)
!!!!!              tt221=tt221 + x_f(1,i1,i2+l,i3)*ceff0_2(l) + x_f(3,i1,i2+l,i3)*eeff0_2(l)
!!!!!              tt222=tt222 + x_f(5,i1,i2+l,i3)*ceff0_2(l) + x_f(7,i1,i2+l,i3)*eeff0_2(l)
!!!!!              tt121=tt121 + x_f(2,i1,i2+l,i3)*eeff0_2(l)
!!!!!           enddo
!!!!!           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
!!!!!           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
!!!!!           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
!!!!!           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
!!!!!           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
!!!!!           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
!!!!!           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
!!!!!           y2_f(4,i1,i2,i3)=tt112
!!!!!           y2_f(2,i1,i2,i3)=y2_f(2,i1,i2,i3)+tt121
!!!!!           y2_f(1,i1,i2,i3)=tt211
!!!!!           y2_f(6,i1,i2,i3)=tt122
!!!!!           y2_f(5,i1,i2,i3)=tt212
!!!!!           y2_f(3,i1,i2,i3)=tt221
!!!!!           y2_f(7,i1,i2,i3)=tt222
!!!!!        enddo
!!!!!     enddo
!!!!!  enddo
!!!!!  !!!$omp enddo
!!!!!
!!!!!  !  call system_clock(ncount5,ncount_rate,ncount_max)
!!!!!  !  tel=dble(ncount5-ncount4)/dble(ncount_rate)
!!!!!  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.d-6*nflop2/tel
!!!!!
!!!!!  ! + (1/2) d^2/dz^2
!!!!!  !!!$omp do
!!!!!  do i2=nfl2,nfu2
!!!!!     do i1=nfl1,nfu1
!!!!!        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!!!!!           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
!!!!!           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
!!!!!           ! Get the effective filters for the z dimension
!!!!!           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
!!!!!           !call getFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
!!!!!           call getFilterQuartic(it, potentialPrefac, hgrid, z0, eeff0(lowfil), 'e')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
!!!!!           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, eeff0_2(lowfil), 'e')
!!!!!           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
!!!!!              t121=t121 + x_f(2,i1,i2,i3+l)*aeff0(l) + x_f(6,i1,i2,i3+l)*beff0(l)
!!!!!              t211=t211 + x_f(1,i1,i2,i3+l)*aeff0(l) + x_f(5,i1,i2,i3+l)*beff0(l)
!!!!!              t122=t122 + x_f(2,i1,i2,i3+l)*ceff0(l) + x_f(6,i1,i2,i3+l)*eeff0(l)
!!!!!              t212=t212 + x_f(1,i1,i2,i3+l)*ceff0(l) + x_f(5,i1,i2,i3+l)*eeff0(l)
!!!!!              t221=t221 + x_f(3,i1,i2,i3+l)*aeff0(l) + x_f(7,i1,i2,i3+l)*beff0(l)
!!!!!              t222=t222 + x_f(3,i1,i2,i3+l)*ceff0(l) + x_f(7,i1,i2,i3+l)*eeff0(l)
!!!!!              t112=t112 + x_f(4,i1,i2,i3+l)*eeff0(l)
!!!!!              tt121=tt121 + x_f(2,i1,i2,i3+l)*aeff0_2(l) + x_f(6,i1,i2,i3+l)*beff0_2(l)
!!!!!              tt211=tt211 + x_f(1,i1,i2,i3+l)*aeff0_2(l) + x_f(5,i1,i2,i3+l)*beff0_2(l)
!!!!!              tt122=tt122 + x_f(2,i1,i2,i3+l)*ceff0_2(l) + x_f(6,i1,i2,i3+l)*eeff0_2(l)
!!!!!              tt212=tt212 + x_f(1,i1,i2,i3+l)*ceff0_2(l) + x_f(5,i1,i2,i3+l)*eeff0_2(l)
!!!!!              tt221=tt221 + x_f(3,i1,i2,i3+l)*aeff0_2(l) + x_f(7,i1,i2,i3+l)*beff0_2(l)
!!!!!              tt222=tt222 + x_f(3,i1,i2,i3+l)*ceff0_2(l) + x_f(7,i1,i2,i3+l)*eeff0_2(l)
!!!!!              tt112=tt112 + x_f(4,i1,i2,i3+l)*eeff0_2(l)
!!!!!           enddo
!!!!!           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
!!!!!           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
!!!!!           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
!!!!!           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
!!!!!           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
!!!!!           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
!!!!!           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
!!!!!           z2_f(4,i1,i2,i3)=z2_f(4,i1,i2,i3)+tt112
!!!!!           z2_f(2,i1,i2,i3)=tt121
!!!!!           z2_f(1,i1,i2,i3)=tt211
!!!!!           z2_f(6,i1,i2,i3)=tt122
!!!!!           z2_f(5,i1,i2,i3)=tt212
!!!!!           z2_f(3,i1,i2,i3)=tt221
!!!!!           z2_f(7,i1,i2,i3)=tt222
!!!!!
!!!!!        enddo
!!!!!     enddo
!!!!!  enddo
!!!!!  !!!$omp enddo
!!!!!
!!!!!
!!!!!
!!!!!  ! Sum up all the contributions
!!!!!  do i3=0,n3
!!!!!      do i2=0,n2
!!!!!          do i1=0,n1
!!!!!              tt0 = x2_c(i1,i2,i3)*y2_c(i1,i2,i3) + x2_c(i1,i2,i3)*z2_c(i1,i2,i3) + y2_c(i1,i2,i3)*z2_c(i1,i2,i3)
!!!!!              y_c(i1,i2,i3) = y_c(i1,i2,i3) + 2.d0*tt0
!!!!!          end do
!!!!!      end do
!!!!!  end do
!!!!!  do i3=nfl3,nfu3
!!!!!      do i2=nfl2,nfu2
!!!!!          do i1=nfl1,nfu1
!!!!!              tt1 = x2_f(1,i1,i2,i3)*y2_f(1,i1,i2,i3) + x2_f(1,i1,i2,i3)*z2_f(1,i1,i2,i3) + y2_f(1,i1,i2,i3)*z2_f(1,i1,i2,i3)
!!!!!              tt2 = x2_f(2,i1,i2,i3)*y2_f(2,i1,i2,i3) + x2_f(2,i1,i2,i3)*z2_f(2,i1,i2,i3) + y2_f(2,i1,i2,i3)*z2_f(2,i1,i2,i3)
!!!!!              tt3 = x2_f(3,i1,i2,i3)*y2_f(3,i1,i2,i3) + x2_f(3,i1,i2,i3)*z2_f(3,i1,i2,i3) + y2_f(3,i1,i2,i3)*z2_f(3,i1,i2,i3)
!!!!!              tt4 = x2_f(4,i1,i2,i3)*y2_f(4,i1,i2,i3) + x2_f(4,i1,i2,i3)*z2_f(4,i1,i2,i3) + y2_f(4,i1,i2,i3)*z2_f(4,i1,i2,i3)
!!!!!              tt5 = x2_f(5,i1,i2,i3)*y2_f(5,i1,i2,i3) + x2_f(5,i1,i2,i3)*z2_f(5,i1,i2,i3) + y2_f(5,i1,i2,i3)*z2_f(5,i1,i2,i3)
!!!!!              tt6 = x2_f(6,i1,i2,i3)*y2_f(6,i1,i2,i3) + x2_f(6,i1,i2,i3)*z2_f(6,i1,i2,i3) + y2_f(6,i1,i2,i3)*z2_f(6,i1,i2,i3)
!!!!!              tt7 = x2_f(7,i1,i2,i3)*y2_f(7,i1,i2,i3) + x2_f(7,i1,i2,i3)*z2_f(7,i1,i2,i3) + y2_f(7,i1,i2,i3)*z2_f(7,i1,i2,i3)
!!!!!              y_f(1,i1,i2,i3) = y_f(1,i1,i2,i3) + 2.d0*tt1
!!!!!              y_f(2,i1,i2,i3) = y_f(2,i1,i2,i3) + 2.d0*tt2
!!!!!              y_f(3,i1,i2,i3) = y_f(3,i1,i2,i3) + 2.d0*tt3
!!!!!              y_f(4,i1,i2,i3) = y_f(4,i1,i2,i3) + 2.d0*tt4
!!!!!              y_f(5,i1,i2,i3) = y_f(5,i1,i2,i3) + 2.d0*tt5
!!!!!              y_f(6,i1,i2,i3) = y_f(6,i1,i2,i3) + 2.d0*tt6
!!!!!              y_f(7,i1,i2,i3) = y_f(7,i1,i2,i3) + 2.d0*tt7
!!!!!          end do
!!!!!      end do
!!!!!  end do
!!!!!  !!! Sum up all the contributions
!!!!!  !!do i3=0,n3
!!!!!  !!    do i2=0,n2
!!!!!  !!        do i1=0,n1
!!!!!  !!            y_c(i1,i2,i3) = x2_c(i1,i2,i3) + y2_c(i1,i2,i3) + z2_c(i1,i2,i3)
!!!!!  !!        end do
!!!!!  !!    end do
!!!!!  !!end do
!!!!!  !!do i3=nfl3,nfu3
!!!!!  !!    do i2=nfl2,nfu2
!!!!!  !!        do i1=nfl1,nfu1
!!!!!  !!            y_f(1,i1,i2,i3) = x2_f(1,i1,i2,i3) + y2_f(1,i1,i2,i3) + z2_f(1,i1,i2,i3)
!!!!!  !!            y_f(2,i1,i2,i3) = x2_f(2,i1,i2,i3) + y2_f(2,i1,i2,i3) + z2_f(2,i1,i2,i3)
!!!!!  !!            y_f(3,i1,i2,i3) = x2_f(3,i1,i2,i3) + y2_f(3,i1,i2,i3) + z2_f(3,i1,i2,i3)
!!!!!  !!            y_f(4,i1,i2,i3) = x2_f(4,i1,i2,i3) + y2_f(4,i1,i2,i3) + z2_f(4,i1,i2,i3)
!!!!!  !!            y_f(5,i1,i2,i3) = x2_f(5,i1,i2,i3) + y2_f(5,i1,i2,i3) + z2_f(5,i1,i2,i3)
!!!!!  !!            y_f(6,i1,i2,i3) = x2_f(6,i1,i2,i3) + y2_f(6,i1,i2,i3) + z2_f(6,i1,i2,i3)
!!!!!  !!            y_f(7,i1,i2,i3) = x2_f(7,i1,i2,i3) + y2_f(7,i1,i2,i3) + z2_f(7,i1,i2,i3)
!!!!!  !!        end do
!!!!!  !!    end do
!!!!!  !!end do
  


  iall=-product(shape(x2_c))*kind(x2_c)
  deallocate(x2_c, stat=istat)
  call memocc(istat, iall, 'x2_c', subname)

  iall=-product(shape(y2_c))*kind(y2_c)
  deallocate(y2_c, stat=istat)
  call memocc(istat, iall, 'y2_c', subname)

  iall=-product(shape(z2_c))*kind(z2_c)
  deallocate(z2_c, stat=istat)
  call memocc(istat, iall, 'z2_c', subname)

  iall=-product(shape(x2_f))*kind(x2_f)
  deallocate(x2_f, stat=istat)
  call memocc(istat, iall, 'x2_f', subname)

  iall=-product(shape(y2_f))*kind(y2_f)
  deallocate(y2_f, stat=istat)
  call memocc(istat, iall, 'y2_f', subname)

  iall=-product(shape(z2_f))*kind(z2_f)
  deallocate(z2_f, stat=istat)
  call memocc(istat, iall, 'z2_f', subname)

  iall=-product(shape(x2_f2))*kind(x2_f2)
  deallocate(x2_f2, stat=istat)
  call memocc(istat, iall, 'x2_f2', subname)

  iall=-product(shape(x2_f3))*kind(x2_f3)
  deallocate(x2_f3, stat=istat)
  call memocc(istat, iall, 'x2_f3', subname)

  iall=-product(shape(y2_f3))*kind(y2_f3)
  deallocate(y2_f3, stat=istat)
  call memocc(istat, iall, 'y2_f3', subname)



!  !!!$omp end parallel
!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel




END SUBROUTINE ConvolQuartic2






!>  Applies the following operation: 
!!  y = [((r-r0)^4)]*x
subroutine ConvolQuartic3(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid, offsetx, offsety, offsetz, &
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzConf, potentialPrefac, it)
  use module_base
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, offsetx, offsety, offsetz
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
real(8),dimension(3):: rxyzConf
real(8):: potentialPrefac
integer:: it
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp):: tt112, tt121, tt122, tt212, tt221, tt222, tt211
  real(wp):: tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7
  !real(kind=8) :: tel
  real(wp), dimension(-3+lowfil:lupfil+3) :: a, aeff0, aeff1, aeff2, aeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: b, beff0, beff1, beff2, beff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: c, ceff0, ceff1, ceff2, ceff3
  real(wp), dimension(lowfil:lupfil) :: e, eeff0, eeff1, eeff2, eeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: aeff0_2, aeff1_2, aeff2_2, aeff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: beff0_2, beff1_2, beff2_2, beff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: ceff0_2, ceff1_2, ceff2_2, ceff3_2
  real(wp), dimension(lowfil:lupfil) :: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  real(wp),dimension(:,:,:),allocatable:: x2_c, y2_c, z2_c, x2_f2, x2_f3, y2_f3
  real(wp),dimension(:,:,:,:),allocatable:: x2_f, y2_f, z2_f
  real(wp),dimension(:,:,:),allocatable:: xa_c, xb_c, xc_c, xe_c, ya_c, yb_c, yc_c, ye_c
  real(wp),dimension(:,:,:,:),allocatable:: xa_f, xb_f, xc_f, xe_f, ya_f, yb_f, yc_f, ye_f
real(8):: x0, y0, z0
real(8):: x1, y1, z1
real(8):: x2, y2, z2
real(8):: x3, y3, z3
integer:: ii, istat, iall
character(len=*),parameter:: subname='ConvolQuartic3'
real(8):: tt00, tt01, tt02, tt03
real(8):: tt10, tt11, tt12, tt13
real(8):: tt20, tt21, tt22, tt23
real(8):: tt30, tt31, tt32, tt33
real(8):: tt40, tt41, tt42, tt43
real(8):: tt50, tt51, tt52, tt53
real(8):: tt60, tt61, tt62, tt63
real(8):: tt70, tt71, tt72, tt73


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
!
!
!  !---------------------------------------------------------------------------

  ! Scaling function part

!  call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)





allocate(x2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, x2_c, 'x2_c', subname)
allocate(y2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, y2_c, 'y2_c', subname)
allocate(z2_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, z2_c, 'z2_c', subname)
allocate(x2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, x2_f, 'x2_f', subname)
allocate(y2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, y2_f, 'y2_f', subname)
allocate(z2_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, z2_f, 'z2_f', subname)

allocate(x2_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2_f2, 'x2_f2', subname)
allocate(x2_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2_f3, 'x2_f3', subname)
allocate(y2_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y2_f3, 'y2_f3', subname)

allocate(xa_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, xa_c, 'xa_c', subname)
allocate(xb_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, xb_c, 'xb_c', subname)
allocate(xc_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, xc_c, 'xc_c', subname)
allocate(xe_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, xe_c, 'xe_c', subname)

allocate(xa_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, xa_f, 'xa_f', subname)
allocate(xb_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, xb_f, 'xb_f', subname)
allocate(xc_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, xc_f, 'xc_f', subname)
allocate(xe_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, xe_f, 'xe_f', subname)

allocate(ya_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, ya_c, 'ya_c', subname)
allocate(yb_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, yb_c, 'yb_c', subname)
allocate(yc_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, yc_c, 'yc_c', subname)
allocate(ye_c(0:n1,0:n2,0:n3), stat=istat)
call memocc(istat, ye_c, 'ye_c', subname)

allocate(ya_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, ya_f, 'ya_f', subname)
allocate(yb_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, yb_f, 'yb_f', subname)
allocate(yc_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, yc_f, 'yc_f', subname)
allocate(ye_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
call memocc(istat, ye_f, 'ye_f', subname)



! Build x^2 and y^2 convolutions and store them in xa_c, xb_c, xc_c, xe_c, ya_f, yb_f, yc_f, ye_f
call auxiliary_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, offsetx, offsety, offsetz, &
     hgrid, rxyzConf, potentialPrefac, ibyz_c, ibyz_f, ibxz_c, ibxz_f, ibxy_c, ibxy_f, &
     x_c, x_f, &
     xa_c, xb_c, xc_c, xe_c, ya_c, yb_c, yc_c, ye_c, xa_f, xb_f, xc_f, xe_f, ya_f, yb_f, yc_f, ye_f)
     


x2_c=0.d0
y2_c=0.d0
z2_c=0.d0
x2_f=0.d0
y2_f=0.d0
z2_f=0.d0
x2_f2=0.d0
x2_f3=0.d0
y2_f3=0.d0


aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0

aeff0_2=0.d0 ; beff0_2=0.d0 ; ceff0_2=0.d0 ; eeff0_2=0.0
aeff1_2=0.d0 ; beff1_2=0.d0 ; ceff1_2=0.d0 ; eeff1_2=0.0
aeff2_2=0.d0 ; beff2_2=0.d0 ; ceff2_2=0.d0 ; eeff2_2=0.0
aeff3_2=0.d0 ; beff3_2=0.d0 ; ceff3_2=0.d0 ; eeff3_2=0.0


!!$!$omp parallel default(private) &
!!$!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!!$!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!!$!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
  !!!$omp do  
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt0=0.0_wp
              tt1=0.0_wp
              tt2=0.0_wp
              tt3=0.0_wp
              ! Get the effective a-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, aeff3(lowfil), 'a')

!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, aeff1_2(lowfil), 'a')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, aeff2_2(lowfil), 'a')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, aeff3_2(lowfil), 'a')
!!
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, beff1_2(lowfil), 'b')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, beff2_2(lowfil), 'b')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, beff3_2(lowfil), 'b')
!!
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, ceff1_2(lowfil), 'c')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, ceff2_2(lowfil), 'c')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, ceff3_2(lowfil), 'c')
!!
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0_2(lowfil), 'e')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, eeff1_2(lowfil), 'e')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, eeff2_2(lowfil), 'e')
!!              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, eeff3_2(lowfil), 'e')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*aeff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*aeff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*aeff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*aeff3(t-i1-3)
!!
!!                 ! sss coefficients
!!                 tt0a0=tt0a0 + x_c(t,i2,i3)*aeff0_2(t-i1-0)
!!                 tt0a1=tt0a1 + x_c(t,i2,i3)*aeff1_2(t-i1-1)
!!                 tt0a2=tt0a2 + x_c(t,i2,i3)*aeff2_2(t-i1-2)
!!                 tt0a3=tt0a3 + x_c(t,i2,i3)*aeff3_2(t-i1-3)
!!
!!                 tt0b0=tt0b0 + x_c(t,i2,i3)*beff0_2(t-i1-0)
!!                 tt0b1=tt0b1 + x_c(t,i2,i3)*beff1_2(t-i1-1)
!!                 tt0b2=tt0b2 + x_c(t,i2,i3)*beff2_2(t-i1-2)
!!                 tt0b3=tt0b3 + x_c(t,i2,i3)*beff3_2(t-i1-3)
!!
!!                 tt0c0=tt0c0 + x_c(t,i2,i3)*ceff0_2(t-i1-0)
!!                 tt0c1=tt0c1 + x_c(t,i2,i3)*ceff1_2(t-i1-1)
!!                 tt0c2=tt0c2 + x_c(t,i2,i3)*ceff2_2(t-i1-2)
!!                 tt0c3=tt0c3 + x_c(t,i2,i3)*ceff3_2(t-i1-3)
!!
!!                 tt0e0=tt0e0 + x_c(t,i2,i3)*eeff0_2(t-i1-0)
!!                 tt0e1=tt0e1 + x_c(t,i2,i3)*eeff1_2(t-i1-1)
!!                 tt0e2=tt0e2 + x_c(t,i2,i3)*eeff2_2(t-i1-2)
!!                 tt0e3=tt0e3 + x_c(t,i2,i3)*eeff3_2(t-i1-3)
!!
              enddo
              y_c(i1+0,i2,i3)=dyi0!+cprecr*x_c(i1+0,i2,i3)
              y_c(i1+1,i2,i3)=dyi1!+cprecr*x_c(i1+1,i2,i3)
              y_c(i1+2,i2,i3)=dyi2!+cprecr*x_c(i1+2,i2,i3)
              y_c(i1+3,i2,i3)=dyi3!+cprecr*x_c(i1+3,i2,i3)
!!!
!!!              xa_c(i1+0,i2,i3)=tt0a0
!!!              xa_c(i1+1,i2,i3)=tt0a1
!!!              xa_c(i1+2,i2,i3)=tt0a2
!!!              xa_c(i1+3,i2,i3)=tt0a3
!!!
!!!              xb_c(i1+0,i2,i3)=tt0b0
!!!              xb_c(i1+1,i2,i3)=tt0b1
!!!              xb_c(i1+2,i2,i3)=tt0b2
!!!              xb_c(i1+3,i2,i3)=tt0b3
!!!
!!!              xc_c(i1+0,i2,i3)=tt0c0
!!!              xc_c(i1+1,i2,i3)=tt0c1
!!!              xc_c(i1+2,i2,i3)=tt0c2
!!!              xc_c(i1+3,i2,i3)=tt0c3
!!!
!!!              xe_c(i1+0,i2,i3)=tt0e0
!!!              xe_c(i1+1,i2,i3)=tt0e1
!!!              xe_c(i1+2,i2,i3)=tt0e2
!!!              xe_c(i1+3,i2,i3)=tt0e3
           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0.0_wp 
           !!tt0=0.0_wp
           ! Get the effective a-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0_2(lowfil), 'a')
!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0_2(lowfil), 'b')
!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
!!           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0_2(lowfil), 'e')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*aeff0(t-i1)
!!              ! sss coefficients
!!              tt0a0=tt0a0 + x_c(t,i2,i3)*aeff0_2(t-i1-0)
!!
!!              tt0b0=tt0b0 + x_c(t,i2,i3)*beff0_2(t-i1-0)
!!
!!              tt0c0=tt0c0 + x_c(t,i2,i3)*ceff0_2(t-i1-0)
!!
!!              tt0e0=tt0e0 + x_c(t,i2,i3)*eeff0_2(t-i1-0)
           enddo
           y_c(i1,i2,i3)=dyi!+cprecr*x_c(i1,i2,i3)
!!           xa_c(i1+0,i2,i3)=tt0a0
!!
!!           xb_c(i1+0,i2,i3)=tt0b0
!!
!!           xc_c(i1+0,i2,i3)=tt0c0
!!
!!           xe_c(i1+0,i2,i3)=tt0e0
        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              ! Get the effective b-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, beff3(lowfil), 'b')
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + x_f1(t,i2,i3)*beff0(t-i1-0)
                 dyi1=dyi1 + x_f1(t,i2,i3)*beff1(t-i1-1)
                 dyi2=dyi2 + x_f1(t,i2,i3)*beff2(t-i1-2)
                 dyi3=dyi3 + x_f1(t,i2,i3)*beff3(t-i1-3)
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
           tt0=0.0_wp
           ! Get the effective b-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + x_f1(t,i2,i3)*beff0(t-i1)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
!!              tt0=0.0_wp 
!!              tt1=0.0_wp 
!!              tt2=0.0_wp 
!!              tt3=0.0_wp 
              ! Get the effective c-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)
              call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, x3, ceff3(lowfil), 'c')
              !!call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
              !!call getFilterQuadratic(it, potentialPrefac, hgrid, x1, ceff1_2(lowfil), 'c')
              !!call getFilterQuadratic(it, potentialPrefac, hgrid, x2, ceff2_2(lowfil), 'c')
              !!call getFilterQuadratic(it, potentialPrefac, hgrid, x3, ceff3_2(lowfil), 'c')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + x_c(t,i2,i3)*ceff0(t-i1-0)
                 dyi1=dyi1 + x_c(t,i2,i3)*ceff1(t-i1-1)
                 dyi2=dyi2 + x_c(t,i2,i3)*ceff2(t-i1-2)
                 dyi3=dyi3 + x_c(t,i2,i3)*ceff3(t-i1-3)
                 !!tt0=tt0 + x_c(t,i2,i3)*ceff0_2(t-i1-0)
                 !!tt1=tt1 + x_c(t,i2,i3)*ceff1_2(t-i1-1)
                 !!tt2=tt2 + x_c(t,i2,i3)*ceff2_2(t-i1-2)
                 !!tt3=tt3 + x_c(t,i2,i3)*ceff3_2(t-i1-3)
              enddo
              y_f(1,i1+0,i2,i3)=dyi0
              y_f(1,i1+1,i2,i3)=dyi1
              y_f(1,i1+2,i2,i3)=dyi2
              y_f(1,i1+3,i2,i3)=dyi3
              !!x2_f(1,i1+0,i2,i3)=tt0
              !!x2_f(1,i1+1,i2,i3)=tt1
              !!x2_f(1,i1+2,i2,i3)=tt2
              !!x2_f(1,i1+3,i2,i3)=tt3
           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0.0_wp 
           tt0=0.0_wp 
           ! Get the effective c-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           !!call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0_2(lowfil), 'c')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + x_c(t,i2,i3)*ceff0(t-i1)
              !!tt0=tt0 + x_c(t,i2,i3)*ceff0_2(t-i1)
           enddo
           y_f(1,i1,i2,i3)=dyi
           !!x2_f(1,i1,i2,i3)=tt0
        enddo
     enddo
  enddo
  !!!$omp enddo





  ! wavelet part

  !!!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
           ! Get the effective filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, x0, eeff0(lowfil), 'e')
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*aeff0(l) + x_f(5,i1+l,i2,i3)*beff0(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*aeff0(l) + x_f(3,i1+l,i2,i3)*beff0(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*aeff0(l) + x_f(7,i1+l,i2,i3)*beff0(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*ceff0(l) + x_f(5,i1+l,i2,i3)*eeff0(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*ceff0(l) + x_f(3,i1+l,i2,i3)*eeff0(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*ceff0(l) + x_f(7,i1+l,i2,i3)*eeff0(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*eeff0(l)
           enddo
           y_f(4,i1,i2,i3)=t112!+cprecr*x_f(4,i1,i2,i3)
           y_f(2,i1,i2,i3)=t121!+cprecr*x_f(2,i1,i2,i3)
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
           y_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
           y_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
  !!!$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel



!!  ! copy to auxiliary array for better memory access
!!  do i3=nfl3,nfu3
!!      do i2=nfl2,nfu2
!!          do i1=nfl1,nfu1
!!              !x2_f1(i1,i2,i3)=x2_f(1,i1,i2,i3)
!!              x2_f2(i2,i1,i3)=x2_f(2,i1,i2,i3)
!!              x2_f3(i3,i1,i2)=x2_f(4,i1,i2,i3)
!!          end do
!!      end do
!!  end do







  
  !  call system_clock(ncount1,ncount_rate,ncount_max)
  !  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, aeff3_2(lowfil), 'a')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*aeff0(t-i2-0) + 2.d0*xa_c(i1,t,i3)*aeff0_2(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*aeff1(t-i2-1) + 2.d0*xa_c(i1,t,i3)*aeff1_2(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*aeff2(t-i2-2) + 2.d0*xa_c(i1,t,i3)*aeff2_2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*aeff3(t-i2-3) + 2.d0*xa_c(i1,t,i3)*aeff3_2(t-i2-3)
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
           ! Get the effective a-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + x_c(i1,t,i3)*aeff0(t-i2) + 2.d0*xa_c(i1,t,i3)*aeff0_2(t-i2)
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
              ! Get the effective b-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, beff3(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, aeff3_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, beff3_2(lowfil), 'b')
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + x_f2(t,i1,i3)*beff0(t-i2-0) + 2.d0*xb_f(1,i1,t,i3)*aeff0_2(t-i2-0) + 2.d0*(xa_f(2,i1,t,i3)+xb_f(3,i1,t,i3))*beff0_2(t-i2-0)
                 dyi1=dyi1 + x_f2(t,i1,i3)*beff1(t-i2-1) + 2.d0*xb_f(1,i1,t,i3)*aeff1_2(t-i2-1) + 2.d0*(xa_f(2,i1,t,i3)+xb_f(3,i1,t,i3))*beff1_2(t-i2-1)
                 dyi2=dyi2 + x_f2(t,i1,i3)*beff2(t-i2-2) + 2.d0*xb_f(1,i1,t,i3)*aeff2_2(t-i2-2) + 2.d0*(xa_f(2,i1,t,i3)+xb_f(3,i1,t,i3))*beff2_2(t-i2-2)
                 dyi3=dyi3 + x_f2(t,i1,i3)*beff3(t-i2-3) + 2.d0*xb_f(1,i1,t,i3)*aeff3_2(t-i2-3) + 2.d0*(xa_f(2,i1,t,i3)+xb_f(3,i1,t,i3))*beff3_2(t-i2-3)
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
           ! Get the effective b-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi0=dyi0 + x_f2(t,i1,i3)*beff0(t-i2-0) + 2.d0*xb_f(1,i1,t,i3)*aeff0_2(t-i2-0) + 2.d0*(xa_f(2,i1,t,i3)+xb_f(3,i1,t,i3))*beff0_2(t-i2-0)
           enddo
           y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt10=0.0_wp 
              tt11=0.0_wp 
              tt12=0.0_wp 
              tt13=0.0_wp 
              tt20=0.0_wp 
              tt21=0.0_wp 
              tt22=0.0_wp 
              tt23=0.0_wp 
              tt30=0.0_wp 
              tt31=0.0_wp 
              tt32=0.0_wp 
              tt33=0.0_wp 
              ! Get the effective c-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)
              call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, y3, ceff3(lowfil), 'c')

              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, aeff3_2(lowfil), 'a')

              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, ceff3_2(lowfil), 'c')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + x_c(i1,t,i3)*ceff0(t-i2-0)
                 dyi1=dyi1 + x_c(i1,t,i3)*ceff1(t-i2-1)
                 dyi2=dyi2 + x_c(i1,t,i3)*ceff2(t-i2-2)
                 dyi3=dyi3 + x_c(i1,t,i3)*ceff3(t-i2-3)

                 tt10=tt10 + xc_c(i1,t,i3)*aeff0_2(t-i2-0)
                 tt11=tt11 + xc_c(i1,t,i3)*aeff1_2(t-i2-1)
                 tt12=tt12 + xc_c(i1,t,i3)*aeff2_2(t-i2-2)
                 tt13=tt13 + xc_c(i1,t,i3)*aeff3_2(t-i2-3)

                 tt20=tt20 + xa_c(i1,t,i3)*ceff0_2(t-i2-0)
                 tt21=tt21 + xa_c(i1,t,i3)*ceff1_2(t-i2-1)
                 tt22=tt22 + xa_c(i1,t,i3)*ceff2_2(t-i2-2)
                 tt23=tt23 + xa_c(i1,t,i3)*ceff3_2(t-i2-3)

                 tt30=tt30 + xc_c(i1,t,i3)*ceff0_2(t-i2-0)
                 tt31=tt31 + xc_c(i1,t,i3)*ceff1_2(t-i2-1)
                 tt32=tt32 + xc_c(i1,t,i3)*ceff2_2(t-i2-2)
                 tt33=tt33 + xc_c(i1,t,i3)*ceff3_2(t-i2-3)
              enddo
              y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
              y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+dyi1
              y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+dyi2
              y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+dyi3

              y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
              y_f(1,i1,i2+1,i3)=y_f(1,i1,i2+1,i3)+tt11
              y_f(1,i1,i2+2,i3)=y_f(1,i1,i2+2,i3)+tt12
              y_f(1,i1,i2+3,i3)=y_f(1,i1,i2+3,i3)+tt13

              y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
              y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+tt21
              y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+tt22
              y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+tt23

              y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
              y_f(3,i1,i2+1,i3)=y_f(3,i1,i2+1,i3)+tt31
              y_f(3,i1,i2+2,i3)=y_f(3,i1,i2+2,i3)+tt32
              y_f(3,i1,i2+3,i3)=y_f(3,i1,i2+3,i3)+tt33
           enddo
           icur=i2
        else
           icur=ibxz_f(1,i1,i3)
        endif

        do i2=icur,ibxz_f(2,i1,i3)
           dyi0=0.0_wp 
           tt10=0.0_wp 
           tt20=0.0_wp 
           tt30=0.0_wp 
           ! Get the effective c-filters for the y dimension
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi0=dyi0 + x_c(i1,t,i3)*ceff0(t-i2-0)
              tt10=tt10 + xc_c(i1,t,i3)*aeff0_2(t-i2-0)
              tt20=tt20 + xa_c(i1,t,i3)*ceff0_2(t-i2-0)
              tt30=tt30 + xc_c(i1,t,i3)*ceff0_2(t-i2-0)
           enddo
           y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
           y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
           y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
           y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
        enddo
     enddo
  enddo
  !!!$omp enddo


  ! wavelet part

  !!!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           ! Get the effective filters for the y dimension
           tt10 = 0.d0
           tt20 = 0.d0
           tt30 = 0.d0
           tt40 = 0.d0
           tt50 = 0.d0
           tt60 = 0.d0
           tt70 = 0.d0
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, y0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, eeff0_2(lowfil), 'e')
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              tt10 = tt10 + x_f(1,i1,i2+l,i3)*aeff0(l) + x_f(3,i1,i2+l,i3)*beff0(l) + 2.d0*xe_f(1,i1,i2+l,i3)*aeff0_2(l) + 2.d0*(xc_f(2,i1,i2+l,i3)+xe_f(3,i1,i2+l,i3))*beff0_2(l)

              tt20 = tt20 + x_f(2,i1,i2+l,i3)*eeff0(l) +                              2.d0*xb_f(1,i1,i2+l,i3)*ceff0_2(l) + 2.d0*(xa_f(2,i1,i2+l,i3)+xb_f(3,i1,i2+l,i3))*eeff0_2(l)

              tt30 = tt30 + x_f(1,i1,i2+l,i3)*ceff0(l) + x_f(3,i1,i2+l,i3)*eeff0(l) + 2.d0*xe_f(1,i1,i2+l,i3)*ceff0_2(l) + 2.d0*(xc_f(2,i1,i2+l,i3)+xe_f(3,i1,i2+l,i3))*eeff0_2(l)

              tt40 = tt40 + x_f(4,i1,i2+l,i3)*aeff0(l) + x_f(6,i1,i2+l,i3)*beff0(l) + 2.d0*(xa_f(4,i1,i2+l,i3)+xb_f(5,i1,i2+l,i3))*aeff0_2(l) + 2.d0*(xa_f(6,i1,i2+l,i3)+xb_f(7,i1,i2+l,i3))*beff0_2(l)

              tt50 = tt50 + x_f(5,i1,i2+l,i3)*aeff0(l) + x_f(7,i1,i2+l,i3)*beff0(l) + 2.d0*(xc_f(4,i1,i2+l,i3)+xe_f(5,i1,i2+l,i3))*aeff0_2(l) + 2.d0*(xc_f(6,i1,i2+l,i3)+xe_f(7,i1,i2+l,i3))*beff0_2(l)
              
              tt60 = tt60 + x_f(4,i1,i2+l,i3)*ceff0(l) + x_f(6,i1,i2+l,i3)*eeff0(l) + 2.d0*(xa_f(4,i1,i2+l,i3)+xb_f(5,i1,i2+l,i3))*ceff0_2(l) + 2.d0*(xa_f(6,i1,i2+l,i3)+xb_f(7,i1,i2+l,i3))*eeff0_2(l)

              tt70 = tt70 + x_f(5,i1,i2+l,i3)*ceff0(l) + x_f(7,i1,i2+l,i3)*eeff0(l) + 2.d0*(xc_f(4,i1,i2+l,i3)+xe_f(5,i1,i2+l,i3))*ceff0_2(l) + 2.d0*(xc_f(6,i1,i2+l,i3)+xe_f(7,i1,i2+l,i3))*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70
        enddo
     enddo
  enddo
  !!!$omp enddo





  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              ! Get the effective a-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, aeff1(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, aeff2(lowfil), 'a')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, aeff3(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, aeff3_2(lowfil), 'a')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*aeff0(t-i3-0) + 2.d0*(xa_c(i1,i2,t)+ya_c(i1,i2,t))*aeff0_2(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*aeff1(t-i3-1) + 2.d0*(xa_c(i1,i2,t)+ya_c(i1,i2,t))*aeff0_2(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*aeff2(t-i3-2) + 2.d0*(xa_c(i1,i2,t)+ya_c(i1,i2,t))*aeff0_2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*aeff3(t-i3-3) + 2.d0*(xa_c(i1,i2,t)+ya_c(i1,i2,t))*aeff0_2(t-i3-3)
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
           ! Get the effective a-filters for the y dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + x_c(i1,i2,t)*aeff0(t-i3-0) + 2.d0*(xa_c(i1,i2,t)+ya_c(i1,i2,t))*aeff0_2(t-i3-0)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
        enddo
        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              ! Get the effective b-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, beff1(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, beff2(lowfil), 'b')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, beff3(lowfil), 'b')

              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, aeff3_2(lowfil), 'a')

              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, beff1_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, beff2_2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, beff3_2(lowfil), 'b')
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + x_f3(t,i1,i2)*beff0(t-i3-0) + 2.d0*(xb_f(1,i1,i2,t)+yb_f(2,i1,i2,t))*aeff0_2(t-i3-0) + 2.d0*(xa_f(4,i1,i2,t)+xb_f(5,i1,i2,t)+ya_f(4,i1,i2,t)+yb_f(6,i1,i2,t))*beff0_2(t-i3-0)
                 dyi1=dyi1 + x_f3(t,i1,i2)*beff1(t-i3-1) + 2.d0*(xb_f(1,i1,i2,t)+yb_f(2,i1,i2,t))*aeff0_2(t-i3-1) + 2.d0*(xa_f(4,i1,i2,t)+xb_f(5,i1,i2,t)+ya_f(4,i1,i2,t)+yb_f(6,i1,i2,t))*beff0_2(t-i3-1)
                 dyi2=dyi2 + x_f3(t,i1,i2)*beff2(t-i3-2) + 2.d0*(xb_f(1,i1,i2,t)+yb_f(2,i1,i2,t))*aeff0_2(t-i3-2) + 2.d0*(xa_f(4,i1,i2,t)+xb_f(5,i1,i2,t)+ya_f(4,i1,i2,t)+yb_f(6,i1,i2,t))*beff0_2(t-i3-2)
                 dyi3=dyi3 + x_f3(t,i1,i2)*beff3(t-i3-3) + 2.d0*(xb_f(1,i1,i2,t)+yb_f(2,i1,i2,t))*aeff0_2(t-i3-3) + 2.d0*(xa_f(4,i1,i2,t)+xb_f(5,i1,i2,t)+ya_f(4,i1,i2,t)+yb_f(6,i1,i2,t))*beff0_2(t-i3-3)
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
           tt0=0.0_wp
           ! Get the effective b-filters for the y dimension
           !!z0=hgrid*(i3+0)-rxyzConf(3)
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi0=dyi0 + x_f3(t,i1,i2)*beff0(t-i3-0) + 2.d0*(xb_f(1,i1,i2,t)+yb_f(2,i1,i2,t))*aeff0_2(t-i3-0) + 2.d0*(xa_f(4,i1,i2,t)+xb_f(5,i1,i2,t)+ya_f(4,i1,i2,t)+yb_f(6,i1,i2,t))*beff0_2(t-i3-0)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt10 = 0.d0
              tt11 = 0.d0
              tt12 = 0.d0
              tt13 = 0.d0
              tt40 = 0.d0
              tt41 = 0.d0
              tt42 = 0.d0
              tt43 = 0.d0
              tt50 = 0.d0
              tt51 = 0.d0
              tt52 = 0.d0
              tt53 = 0.d0
              tt20 = 0.d0
              tt21 = 0.d0
              tt22 = 0.d0
              tt23 = 0.d0
              tt60 = 0.d0
              tt61 = 0.d0
              tt62 = 0.d0
              tt63 = 0.d0
              ! Get the effective c-filters for the z dimension
              z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
              z1=hgrid*(i3+offsetz+1)-rxyzConf(3)
              z2=hgrid*(i3+offsetz+2)-rxyzConf(3)
              z3=hgrid*(i3+offsetz+3)-rxyzConf(3)
              call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z1, ceff1(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z2, ceff2(lowfil), 'c')
              call getFilterQuartic(it, potentialPrefac, hgrid, z3, ceff3(lowfil), 'c')

              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, aeff1_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, aeff2_2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, aeff3_2(lowfil), 'a')

              call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z1, ceff1_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z2, ceff2_2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, z3, ceff3_2(lowfil), 'c')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + x_c(i1,i2,t)*ceff0(t-i3-0)
                 dyi1=dyi1 + x_c(i1,i2,t)*ceff1(t-i3-1)
                 dyi2=dyi2 + x_c(i1,i2,t)*ceff2(t-i3-2)
                 dyi3=dyi3 + x_c(i1,i2,t)*ceff3(t-i3-3)

                 tt10 = tt10 + 2.d0*xc_c(i1,i2,t)*aeff0_2(t-i3-0)
                 tt11 = tt11 + 2.d0*xc_c(i1,i2,t)*aeff1_2(t-i3-1)
                 tt12 = tt12 + 2.d0*xc_c(i1,i2,t)*aeff2_2(t-i3-2)
                 tt13 = tt13 + 2.d0*xc_c(i1,i2,t)*aeff3_2(t-i3-3)

                 tt40 = tt40 + 2.d0*xa_c(i1,i2,t)*ceff0_2(t-i3-0)
                 tt41 = tt41 + 2.d0*xa_c(i1,i2,t)*ceff1_2(t-i3-1)
                 tt42 = tt42 + 2.d0*xa_c(i1,i2,t)*ceff2_2(t-i3-2)
                 tt43 = tt43 + 2.d0*xa_c(i1,i2,t)*ceff3_2(t-i3-3)
                 
                 tt50 = tt50 + 2.d0*xc_c(i1,i2,t)*ceff0_2(t-i3-0)
                 tt51 = tt51 + 2.d0*xc_c(i1,i2,t)*ceff1_2(t-i3-1)
                 tt52 = tt52 + 2.d0*xc_c(i1,i2,t)*ceff2_2(t-i3-2)
                 tt53 = tt53 + 2.d0*xc_c(i1,i2,t)*ceff3_2(t-i3-3)
                 
                 tt20 = tt20 + 2.d0*yc_c(i1,i2,t)*aeff0_2(t-i3-0)
                 tt21 = tt21 + 2.d0*yc_c(i1,i2,t)*aeff1_2(t-i3-1)
                 tt22 = tt22 + 2.d0*yc_c(i1,i2,t)*aeff2_2(t-i3-2)
                 tt23 = tt23 + 2.d0*yc_c(i1,i2,t)*aeff3_2(t-i3-3)

                 tt40 = tt40 + 2.d0*ya_c(i1,i2,t)*ceff0_2(t-i3-0)
                 tt41 = tt41 + 2.d0*ya_c(i1,i2,t)*ceff1_2(t-i3-1)
                 tt42 = tt42 + 2.d0*ya_c(i1,i2,t)*ceff2_2(t-i3-2)
                 tt43 = tt43 + 2.d0*ya_c(i1,i2,t)*ceff3_2(t-i3-3)

                 tt60 = tt60 + 2.d0*yc_c(i1,i2,t)*ceff0_2(t-i3-0)
                 tt61 = tt61 + 2.d0*yc_c(i1,i2,t)*ceff1_2(t-i3-1)
                 tt62 = tt62 + 2.d0*yc_c(i1,i2,t)*ceff2_2(t-i3-2)
                 tt63 = tt63 + 2.d0*yc_c(i1,i2,t)*ceff3_2(t-i3-3)

              enddo
              y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
              y_f(4,i1,i2,i3+1) = y_f(4,i1,i2,i3+1) + dyi1 + tt41
              y_f(4,i1,i2,i3+2) = y_f(4,i1,i2,i3+2) + dyi2 + tt42
              y_f(4,i1,i2,i3+3) = y_f(4,i1,i2,i3+3) + dyi3 + tt43

              y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
              y_f(1,i1,i2,i3+1) = y_f(1,i1,i2,i3+1) + tt11
              y_f(1,i1,i2,i3+2) = y_f(1,i1,i2,i3+2) + tt12
              y_f(1,i1,i2,i3+3) = y_f(1,i1,i2,i3+3) + tt13

              y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
              y_f(5,i1,i2,i3+1) = y_f(5,i1,i2,i3+1) + tt51
              y_f(5,i1,i2,i3+2) = y_f(5,i1,i2,i3+2) + tt52
              y_f(5,i1,i2,i3+3) = y_f(5,i1,i2,i3+3) + tt53

              y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
              y_f(2,i1,i2,i3+1) = y_f(2,i1,i2,i3+1) + tt21
              y_f(2,i1,i2,i3+2) = y_f(2,i1,i2,i3+2) + tt22
              y_f(2,i1,i2,i3+3) = y_f(2,i1,i2,i3+3) + tt23

              y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
              y_f(6,i1,i2,i3+1) = y_f(6,i1,i2,i3+1) + tt61
              y_f(6,i1,i2,i3+2) = y_f(6,i1,i2,i3+2) + tt62
              y_f(6,i1,i2,i3+3) = y_f(6,i1,i2,i3+3) + tt63
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi0=0.0_wp 
           tt10 = 0.d0
           tt40 = 0.d0
           tt50 = 0.d0
           tt20 = 0.d0
           tt40 = 0.d0
           tt60 = 0.d0
           ! Get the effective c-filters for the z dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + x_c(i1,i2,t)*ceff0(t-i3-0)
              tt10 = tt10 + 2.d0*x_c(i1,i2,t)*aeff0_2(t-i3-0)
              tt40 = tt40 + 2.d0*xa_c(i1,i2,t)*ceff0_2(t-i3-0)
              tt50 = tt50 + 2.d0*xc_c(i1,i2,t)*ceff0_2(t-i3-0)
              tt20 = tt20 + 2.d0*yc_c(i1,i2,t)*aeff0_2(t-i3-0)
              tt40 = tt40 + 2.d0*ya_c(i1,i2,t)*ceff0_2(t-i3-0)
              tt60 = tt60 + 2.d0*yc_c(i1,i2,t)*ceff0_2(t-i3-0)
           enddo
           y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
           y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
           y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
           y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
           y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
        enddo
     enddo
  enddo
  !!!$omp enddo


  ! wavelet part

  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           tt10 = 0.d0
           tt20 = 0.d0
           tt30 = 0.d0
           tt40 = 0.d0
           tt50 = 0.d0
           tt60 = 0.d0
           tt70 = 0.d0
           ! Get the effective filters for the z dimension
           z0=hgrid*(i3+offsetz+0)-rxyzConf(3)
           !call getFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, aeff0(lowfil), 'a')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, beff0(lowfil), 'b')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, ceff0(lowfil), 'c')
           call getFilterQuartic(it, potentialPrefac, hgrid, z0, eeff0(lowfil), 'e')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, aeff0_2(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, beff0_2(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, ceff0_2(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, z0, eeff0_2(lowfil), 'e')
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)

              tt10 = tt10 + x_f(1,i1,i2,i3+l)*aeff0(l) + x_f(5,i1,i2,i3+l)*beff0(l) + 2.d0*                    xe_f(1,i1,i2,i3+l) *aeff0_2(l) + 2.d0*(xc_f(4,i1,i2,i3+l)+xe_f(5,i1,i2,i3+l))*beff0_2(l) &
                                                                                    + 2.d0*(ya_f(1,i1,i2,i3+l)+yb_f(3,i1,i2,i3+l))*aeff0_2(l) + 2.d0*(ya_f(5,i1,i2,i3+l)+yb_f(7,i1,i2,i3+l))*beff0_2(l)

              tt20 = tt20 + x_f(2,i1,i2,i3+l)*aeff0(l) + x_f(6,i1,i2,i3+l)*beff0(l) + 2.d0*(xa_f(2,i1,i2,i3+l)+xb_f(3,i1,i2,i3+l))*aeff0_2(l) + 2.d0*(xa_f(6,i1,i2,i3+l)+xb_f(7,i1,i2,i3+l))*beff0_2(l) &
                                                                                    + 2.d0*                    ye_f(2,i1,i2,i3+l) *aeff0_2(l) + 2.d0*(yc_f(4,i1,i2,i3+l)+ye_f(7,i1,i2,i3+l))*beff0_2(l)

              tt30 = tt30 + x_f(3,i1,i2,i3+l)*aeff0(l) + x_f(7,i1,i2,i3+l)*beff0(l) + 2.d0*(xc_f(2,i1,i2,i3+l)+xe_f(3,i1,i2,i3+l))*aeff0_2(l) + 2.d0*(xc_f(6,i1,i2,i3+l)+xe_f(7,i1,i2,i3+l))*beff0_2(l) &
                                                                                    + 2.d0*(yc_f(1,i1,i2,i3+l)+ye_f(3,i1,i2,i3+l))*aeff0_2(l) + 2.d0*(yc_f(5,i1,i2,i3+l)+ye_f(7,i1,i2,i3+l))*beff0_2(l)

              tt40 = tt40 + x_f(4,i1,i2,i3+l)*eeff0(l)                              + 2.d0*                    xb_f(1,i1,i2,i3+l) *ceff0_2(l) + 2.d0*(xa_f(4,i1,i2,i3+l)+xb_f(5,i1,i2,i3+l))*eeff0_2(l) &
                                                                                    + 2.d0*                    yb_f(2,i1,i2,i3+l) *ceff0_2(l) + 2.d0*(ya_f(4,i1,i2,i3+l)+yb_f(6,i1,i2,i3+l))*eeff0_2(l)

              tt50 = tt50 + x_f(1,i1,i2,i3+l)*ceff0(l) + x_f(5,i1,i2,i3+l)*eeff0(l) + 2.d0*                    xe_f(1,i1,i2,i3+l) *ceff0_2(l) + 2.d0*(xc_f(4,i1,i2,i3+l)+xe_f(5,i1,i2,i3+l))*eeff0_2(l) &
                                                                                    + 2.d0*(ya_f(1,i1,i2,i3+l)+yb_f(3,i1,i2,i3+l))*ceff0_2(l) + 2.d0*(ya_f(5,i1,i2,i3+l)+yb_f(7,i1,i2,i3+l))*eeff0_2(l)

              tt60 = tt60 + x_f(2,i1,i2,i3+l)*ceff0(l) + x_f(6,i1,i2,i3+l)*eeff0(l) + 2.d0*(xa_f(2,i1,i2,i3+l)+xb_f(3,i1,i2,i3+l))*ceff0_2(l) + 2.d0*(xa_f(6,i1,i2,i3+l)+xb_f(7,i1,i2,i3+l))*eeff0_2(l) &
                                                                                    + 2.d0*                    ye_f(2,i1,i2,i3+l) *ceff0_2(l) + 2.d0*(yc_f(4,i1,i2,i3+l)+ye_f(6,i1,i2,i3+l))*eeff0_2(l)

              tt70 = tt70 + x_f(3,i1,i2,i3+l)*ceff0(l) + x_f(7,i1,i2,i3+l)*eeff0(l) + 2.d0*(xc_f(2,i1,i2,i3+l)+xe_f(3,i1,i2,i3+l))*ceff0_2(l) + 2.d0*(xc_f(6,i1,i2,i3+l)+xe_f(7,i1,i2,i3+l))*eeff0_2(l) &
                                                                                    + 2.d0*(yc_f(1,i1,i2,i3+l)+ye_f(3,i1,i2,i3+l))*ceff0_2(l) + 2.d0*(yc_f(5,i1,i2,i3+l)+ye_f(7,i1,i2,i3+l))*eeff0_2(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70
        enddo
     enddo
  enddo
  !!!$omp enddo







  iall=-product(shape(x2_c))*kind(x2_c)
  deallocate(x2_c, stat=istat)
  call memocc(istat, iall, 'x2_c', subname)

  iall=-product(shape(y2_c))*kind(y2_c)
  deallocate(y2_c, stat=istat)
  call memocc(istat, iall, 'y2_c', subname)

  iall=-product(shape(z2_c))*kind(z2_c)
  deallocate(z2_c, stat=istat)
  call memocc(istat, iall, 'z2_c', subname)

  iall=-product(shape(x2_f))*kind(x2_f)
  deallocate(x2_f, stat=istat)
  call memocc(istat, iall, 'x2_f', subname)

  iall=-product(shape(y2_f))*kind(y2_f)
  deallocate(y2_f, stat=istat)
  call memocc(istat, iall, 'y2_f', subname)

  iall=-product(shape(z2_f))*kind(z2_f)
  deallocate(z2_f, stat=istat)
  call memocc(istat, iall, 'z2_f', subname)

  iall=-product(shape(x2_f2))*kind(x2_f2)
  deallocate(x2_f2, stat=istat)
  call memocc(istat, iall, 'x2_f2', subname)

  iall=-product(shape(x2_f3))*kind(x2_f3)
  deallocate(x2_f3, stat=istat)
  call memocc(istat, iall, 'x2_f3', subname)

  iall=-product(shape(y2_f3))*kind(y2_f3)
  deallocate(y2_f3, stat=istat)
  call memocc(istat, iall, 'y2_f3', subname)


  iall=-product(shape(xa_c))*kind(xa_c)
  deallocate(xa_c, stat=istat)
  call memocc(istat, iall, 'xa_c', subname)

  iall=-product(shape(xb_c))*kind(xb_c)
  deallocate(xb_c, stat=istat)
  call memocc(istat, iall, 'xb_c', subname)

  iall=-product(shape(xc_c))*kind(xc_c)
  deallocate(xc_c, stat=istat)
  call memocc(istat, iall, 'xc_c', subname)

  iall=-product(shape(xe_c))*kind(xe_c)
  deallocate(xe_c, stat=istat)
  call memocc(istat, iall, 'xe_c', subname)

  iall=-product(shape(ya_c))*kind(ya_c)
  deallocate(ya_c, stat=istat)
  call memocc(istat, iall, 'ya_c', subname)

  iall=-product(shape(yb_c))*kind(yb_c)
  deallocate(yb_c, stat=istat)
  call memocc(istat, iall, 'yb_c', subname)

  iall=-product(shape(yc_c))*kind(yc_c)
  deallocate(yc_c, stat=istat)
  call memocc(istat, iall, 'yc_c', subname)

  iall=-product(shape(ye_c))*kind(ye_c)
  deallocate(ye_c, stat=istat)
  call memocc(istat, iall, 'ye_c', subname)


  iall=-product(shape(xa_f))*kind(xa_f)
  deallocate(xa_f, stat=istat)
  call memocc(istat, iall, 'xa_f', subname)

  iall=-product(shape(xb_f))*kind(xb_f)
  deallocate(xb_f, stat=istat)
  call memocc(istat, iall, 'xb_f', subname)

  iall=-product(shape(xc_f))*kind(xc_f)
  deallocate(xc_f, stat=istat)
  call memocc(istat, iall, 'xc_f', subname)

  iall=-product(shape(xe_f))*kind(xe_f)
  deallocate(xe_f, stat=istat)
  call memocc(istat, iall, 'xe_f', subname)

  iall=-product(shape(ya_f))*kind(ya_f)
  deallocate(ya_f, stat=istat)
  call memocc(istat, iall, 'ya_f', subname)

  iall=-product(shape(yb_f))*kind(yb_f)
  deallocate(yb_f, stat=istat)
  call memocc(istat, iall, 'yb_f', subname)

  iall=-product(shape(yc_f))*kind(yc_f)
  deallocate(yc_f, stat=istat)
  call memocc(istat, iall, 'yc_f', subname)

  iall=-product(shape(ye_f))*kind(ye_f)
  deallocate(ye_f, stat=istat)
  call memocc(istat, iall, 'ye_f', subname)


!  !!!$omp end parallel
!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel




END SUBROUTINE ConvolQuartic3












subroutine createDerivativeBasis(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     w_c, w_f, w_f1, w_f2, w_f3, x_c, x_f, y_c, y_f, z_c, z_f)
     
  use module_base
  use filterModule
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: w_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: w_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: w_f3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: z_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: z_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp), dimension(-3+lowfil:lupfil+3) :: ad1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: bd1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: cd1_ext
  !real(kind=8) :: tel
  !real(wp), dimension(-3+lowfil:lupfil+3) :: a, aeff0, aeff1, aeff2, aeff3
  !real(wp), dimension(-3+lowfil:lupfil+3) :: b, beff0, beff1, beff2, beff3
  !real(wp), dimension(-3+lowfil:lupfil+3) :: c, ceff0, ceff1, ceff2, ceff3
  !real(wp), dimension(lowfil:lupfil) :: e, eeff0, eeff1, eeff2, eeff3
!real(8):: x0, y0, z0
!real(8):: x1, y1, z1
!real(8):: x2, y2, z2
!real(8):: x3, y3, z3
!integer:: ii


!!!$$  scale=-.5_wp/real(hgrid**2,wp)
!!!$$
!!!$$  !---------------------------------------------------------------------------
!!!$$  ! second derivative filters for Daubechies 16
!!!$$  !  <phi|D^2|phi_i>
!!!$$  a(0)=   -3.5536922899131901941296809374_wp*scale
!!!$$  a(1)=    2.2191465938911163898794546405_wp*scale
!!!$$  a(2)=   -0.6156141465570069496314853949_wp*scale
!!!$$  a(3)=    0.2371780582153805636239247476_wp*scale
!!!$$  a(4)=   -0.0822663999742123340987663521_wp*scale
!!!$$  a(5)=    0.02207029188482255523789911295638968409_wp*scale
!!!$$  a(6)=   -0.409765689342633823899327051188315485e-2_wp*scale
!!!$$  a(7)=    0.45167920287502235349480037639758496e-3_wp*scale
!!!$$  a(8)=   -0.2398228524507599670405555359023135e-4_wp*scale
!!!$$  a(9)=    2.0904234952920365957922889447361e-6_wp*scale
!!!$$  a(10)=  -3.7230763047369275848791496973044e-7_wp*scale
!!!$$  a(11)=  -1.05857055496741470373494132287e-8_wp*scale
!!!$$  a(12)=  -5.813879830282540547959250667e-11_wp*scale
!!!$$  a(13)=   2.70800493626319438269856689037647576e-13_wp*scale
!!!$$  a(14)=  -6.924474940639200152025730585882e-18_wp*scale
!!!$$
!!!$$  a(15)=0.0_wp
!!!$$  a(16)=0.0_wp 
!!!$$  a(17)=0.0_wp
!!!$$  
!!!$$  do i=1,14+3
!!!$$     a(-i)=a(i)
!!!$$  enddo
!!!$$  !  <phi|D^2|psi_i>
!!!$$  c(-17)=0.0_wp
!!!$$  c(-16)=0.0_wp
!!!$$  c(-15)=0.0_wp
!!!$$  
!!!$$  c(-14)=     -3.869102413147656535541850057188e-18_wp*scale
!!!$$  c(-13)=      1.5130616560866154733900029272077362e-13_wp*scale
!!!$$  c(-12)=     -3.2264702314010525539061647271983988409e-11_wp*scale
!!!$$  c(-11)=     -5.96264938781402337319841002642e-9_wp*scale
!!!$$  c(-10)=     -2.1656830629214041470164889350342e-7_wp*scale
!!!$$  c(-9 )=      8.7969704055286288323596890609625e-7_wp*scale
!!!$$  c(-8 )=     -0.00001133456724516819987751818232711775_wp*scale
!!!$$  c(-7 )=      0.00021710795484646138591610188464622454_wp*scale
!!!$$  c(-6 )=     -0.0021356291838797986414312219042358542_wp*scale
!!!$$  c(-5 )=      0.00713761218453631422925717625758502986_wp*scale
!!!$$  c(-4 )=     -0.0284696165863973422636410524436931061_wp*scale
!!!$$  c(-3 )=      0.14327329352510759457155821037742893841_wp*scale
!!!$$  c(-2 )=     -0.42498050943780130143385739554118569733_wp*scale
!!!$$  c(-1 )=      0.65703074007121357894896358254040272157_wp*scale
!!!$$  c( 0 )=     -0.42081655293724308770919536332797729898_wp*scale
!!!$$  c( 1 )=     -0.21716117505137104371463587747283267899_wp*scale
!!!$$  c( 2 )=      0.63457035267892488185929915286969303251_wp*scale
!!!$$  c( 3 )=     -0.53298223962800395684936080758073568406_wp*scale
!!!$$  c( 4 )=      0.23370490631751294307619384973520033236_wp*scale
!!!$$  c( 5 )=     -0.05657736973328755112051544344507997075_wp*scale
!!!$$  c( 6 )=      0.0080872029411844780634067667008050127_wp*scale
!!!$$  c( 7 )=     -0.00093423623304808664741804536808932984_wp*scale
!!!$$  c( 8 )=      0.00005075807947289728306309081261461095_wp*scale
!!!$$  c( 9 )=     -4.62561497463184262755416490048242e-6_wp*scale
!!!$$  c( 10)=      6.3919128513793415587294752371778e-7_wp*scale
!!!$$  c( 11)=      1.87909235155149902916133888931e-8_wp*scale
!!!$$  c( 12)=      1.04757345962781829480207861447155543883e-10_wp*scale
!!!$$  c( 13)=     -4.84665690596158959648731537084025836e-13_wp*scale
!!!$$  c( 14)=      1.2392629629188986192855777620877e-17_wp*scale
!!!$$
!!!$$  c(15)=0.0_wp
!!!$$  c(16)=0.0_wp
!!!$$  c(17)=0.0_wp
!!!$$  !  <psi|D^2|phi_i>
!!!$$  do i=-14-3,14+3
!!!$$     b(i)=c(-i)
!!!$$  enddo
!!!$$  !<psi|D^2|psi_i>
!!!$$  e(0)=   -24.875846029392331358907766562_wp*scale
!!!$$  e(1)=   -7.1440597663471719869313377994_wp*scale
!!!$$  e(2)=   -0.04251705323669172315864542163525830944_wp*scale
!!!$$  e(3)=   -0.26995931336279126953587091167128839196_wp*scale
!!!$$  e(4)=    0.08207454169225172612513390763444496516_wp*scale
!!!$$  e(5)=   -0.02207327034586634477996701627614752761_wp*scale
!!!$$  e(6)=    0.00409765642831595181639002667514310145_wp*scale
!!!$$  e(7)=   -0.00045167920287507774929432548999880117_wp*scale
!!!$$  e(8)=    0.00002398228524507599670405555359023135_wp*scale
!!!$$  e(9)=   -2.0904234952920365957922889447361e-6_wp*scale
!!!$$  e(10)=   3.7230763047369275848791496973044e-7_wp*scale
!!!$$  e(11)=   1.05857055496741470373494132287e-8_wp*scale
!!!$$  e(12)=   5.8138798302825405479592506674648873655e-11_wp*scale
!!!$$  e(13)=  -2.70800493626319438269856689037647576e-13_wp*scale
!!!$$  e(14)=   6.924474940639200152025730585882e-18_wp*scale
!!!$$  do i=1,14
!!!$$     e(-i)=e(i)
!!!$$  enddo



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
!
!
!  !---------------------------------------------------------------------------

  ! Scaling function part

!  call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)



!aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
!aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
!aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
!aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0

! Copy the filters to the 'extended filters', i.e. add some zeros.
! This seems to be required since we use loop unrolling.
ad1_ext=0.d0
bd1_ext=0.d0
cd1_ext=0.d0
do i=lowfil,lupfil
    ad1_ext(i)=ad1(i)
    !write(*,'(a,i4,2es15.8)') 'i, ad1_ext(i), ad1(i)', i, ad1_ext(i), ad1(i)
    bd1_ext(i)=bd1(i)
    !write(*,'(a,i4,2es15.8)') 'i, bd1_ext(i), bd1(i)', i, bd1_ext(i), bd1(i)
    cd1_ext(i)=cd1(i)
    !write(*,'(a,i4,2es15.8)') 'i, cd1_ext(i), cd1(i)', i, cd1_ext(i), cd1(i)
end do


!!!$omp parallel default(private) &
!!!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!!!$omp shared(ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,w_c,w_f,y_c,y_f)& 
!!!$omp shared(w_f1,w_f2,w_f3,ad1_ext,bd1_ext,cd1_ext)
  !!!$omp do  
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              !! Get the effective a-filters for the x dimension
              !x0=hgrid*(i1+0)-rxyzConf(1)
              !x1=hgrid*(i1+1)-rxyzConf(1)
              !x2=hgrid*(i1+2)-rxyzConf(1)
              !x3=hgrid*(i1+3)-rxyzConf(1)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x1, aeff1(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x2, aeff2(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x3, aeff3(lowfil), 'a')
!write(*,'(a,2i5,3x,2i6)') 'i2, i3, max(ibyz_c(1,i2,i3),lowfil+i1), min(lupfil+i1+3,ibyz_c(2,i2,i3))', i2, i3, max(ibyz_c(1,i2,i3),lowfil+i1), min(lupfil+i1+3,ibyz_c(2,i2,i3))
!write(*,'(a,2i6)') 'ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)', ibyz_c(1,i2,i3), ibyz_c(2,i2,i3)
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
!write(*,*) 't, i1, t-i1', t, i1, t-i1
                 dyi0=dyi0 + w_c(t,i2,i3)*ad1_ext(t-i1-0)
                 dyi1=dyi1 + w_c(t,i2,i3)*ad1_ext(t-i1-1)
                 dyi2=dyi2 + w_c(t,i2,i3)*ad1_ext(t-i1-2)
                 dyi3=dyi3 + w_c(t,i2,i3)*ad1_ext(t-i1-3)
              enddo
              x_c(i1+0,i2,i3)=dyi0!+cprecr*w_c(i1+0,i2,i3)
              x_c(i1+1,i2,i3)=dyi1!+cprecr*w_c(i1+1,i2,i3)
              x_c(i1+2,i2,i3)=dyi2!+cprecr*w_c(i1+2,i2,i3)
              x_c(i1+3,i2,i3)=dyi3!+cprecr*w_c(i1+3,i2,i3)
           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0.0_wp 
           !! Get the effective a-filters for the x dimension
           !x0=hgrid*(i1+0)-rxyzConf(1)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + w_c(t,i2,i3)*ad1_ext(t-i1)
           enddo
           x_c(i1,i2,i3)=dyi!+cprecr*w_c(i1,i2,i3)
        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              !! Get the effective b-filters for the x dimension
              !x0=hgrid*(i1+0)-rxyzConf(1)
              !x1=hgrid*(i1+1)-rxyzConf(1)
              !x2=hgrid*(i1+2)-rxyzConf(1)
              !x3=hgrid*(i1+3)-rxyzConf(1)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x1, beff1(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x2, beff2(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x3, beff3(lowfil), 'b')
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + w_f1(t,i2,i3)*bd1_ext(t-i1-0)
                 dyi1=dyi1 + w_f1(t,i2,i3)*bd1_ext(t-i1-1)
                 dyi2=dyi2 + w_f1(t,i2,i3)*bd1_ext(t-i1-2)
                 dyi3=dyi3 + w_f1(t,i2,i3)*bd1_ext(t-i1-3)
              enddo
              x_c(i1+0,i2,i3)=x_c(i1+0,i2,i3)+dyi0
              x_c(i1+1,i2,i3)=x_c(i1+1,i2,i3)+dyi1
              x_c(i1+2,i2,i3)=x_c(i1+2,i2,i3)+dyi2
              x_c(i1+3,i2,i3)=x_c(i1+3,i2,i3)+dyi3
           enddo
           istart=i1
        endif

        do i1=istart,iend
           dyi=0.0_wp
           !! Get the effective b-filters for the x dimension
           !x0=hgrid*(i1+0)-rxyzConf(1)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + w_f1(t,i2,i3)*bd1_ext(t-i1)
           enddo
           x_c(i1,i2,i3)=x_c(i1,i2,i3)+dyi
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              !! Get the effective c-filters for the x dimension
              !x0=hgrid*(i1+0)-rxyzConf(1)
              !x1=hgrid*(i1+1)-rxyzConf(1)
              !x2=hgrid*(i1+2)-rxyzConf(1)
              !x3=hgrid*(i1+3)-rxyzConf(1)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x1, ceff1(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x2, ceff2(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x3, ceff3(lowfil), 'c')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + w_c(t,i2,i3)*cd1_ext(t-i1-0)
                 dyi1=dyi1 + w_c(t,i2,i3)*cd1_ext(t-i1-1)
                 dyi2=dyi2 + w_c(t,i2,i3)*cd1_ext(t-i1-2)
                 dyi3=dyi3 + w_c(t,i2,i3)*cd1_ext(t-i1-3)
              enddo
              x_f(1,i1+0,i2,i3)=dyi0
              x_f(1,i1+1,i2,i3)=dyi1
              x_f(1,i1+2,i2,i3)=dyi2
              x_f(1,i1+3,i2,i3)=dyi3
           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0.0_wp 
           !! Get the effective c-filters for the x dimension
           !x0=hgrid*(i1+0)-rxyzConf(1)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + w_c(t,i2,i3)*cd1_ext(t-i1)
           enddo
           x_f(1,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo
  
  !  call system_clock(ncount1,ncount_rate,ncount_max)
  !  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              !! Get the effective a-filters for the y dimension
              !y0=hgrid*(i2+0)-rxyzConf(2)
              !y1=hgrid*(i2+1)-rxyzConf(2)
              !y2=hgrid*(i2+2)-rxyzConf(2)
              !y3=hgrid*(i2+3)-rxyzConf(2)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y1, aeff1(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y2, aeff2(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y3, aeff3(lowfil), 'a')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + w_c(i1,t,i3)*ad1_ext(t-i2-0)
                 dyi1=dyi1 + w_c(i1,t,i3)*ad1_ext(t-i2-1)
                 dyi2=dyi2 + w_c(i1,t,i3)*ad1_ext(t-i2-2)
                 dyi3=dyi3 + w_c(i1,t,i3)*ad1_ext(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=dyi0
              y_c(i1,i2+1,i3)=dyi1
              y_c(i1,i2+2,i3)=dyi2
              y_c(i1,i2+3,i3)=dyi3
           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           dyi=0.0_wp 
           !! Get the effective a-filters for the y dimension
           !y0=hgrid*(i2+0)-rxyzConf(2)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + w_c(i1,t,i3)*ad1_ext(t-i2)
           enddo
           y_c(i1,i2,i3)=dyi
        enddo
        istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
        iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i2=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              !! Get the effective b-filters for the y dimension
              !y0=hgrid*(i2+0)-rxyzConf(2)
              !y1=hgrid*(i2+1)-rxyzConf(2)
              !y2=hgrid*(i2+2)-rxyzConf(2)
              !y3=hgrid*(i2+3)-rxyzConf(2)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y1, beff1(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y2, beff2(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y3, beff3(lowfil), 'b')
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + w_f2(t,i1,i3)*bd1_ext(t-i2-0)
                 dyi1=dyi1 + w_f2(t,i1,i3)*bd1_ext(t-i2-1)
                 dyi2=dyi2 + w_f2(t,i1,i3)*bd1_ext(t-i2-2)
                 dyi3=dyi3 + w_f2(t,i1,i3)*bd1_ext(t-i2-3)
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
           !! Get the effective b-filters for the y dimension
           !y0=hgrid*(i2+0)-rxyzConf(2)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + w_f2(t,i1,i3)*bd1_ext(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              !! Get the effective c-filters for the y dimension
              !y0=hgrid*(i2+0)-rxyzConf(2)
              !y1=hgrid*(i2+1)-rxyzConf(2)
              !y2=hgrid*(i2+2)-rxyzConf(2)
              !y3=hgrid*(i2+3)-rxyzConf(2)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y1, ceff1(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y2, ceff2(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y3, ceff3(lowfil), 'c')
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + w_c(i1,t,i3)*cd1_ext(t-i2-0)
                 dyi1=dyi1 + w_c(i1,t,i3)*cd1_ext(t-i2-1)
                 dyi2=dyi2 + w_c(i1,t,i3)*cd1_ext(t-i2-2)
                 dyi3=dyi3 + w_c(i1,t,i3)*cd1_ext(t-i2-3)
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
           !! Get the effective c-filters for the y dimension
           !y0=hgrid*(i2+0)-rxyzConf(2)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + w_c(i1,t,i3)*cd1_ext(t-i2)
           enddo
           y_f(2,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo


  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              !! Get the effective a-filters for the z dimension
              !z0=hgrid*(i3+0)-rxyzConf(3)
              !z1=hgrid*(i3+1)-rxyzConf(3)
              !z2=hgrid*(i3+2)-rxyzConf(3)
              !z3=hgrid*(i3+3)-rxyzConf(3)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z1, aeff1(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z2, aeff2(lowfil), 'a')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z3, aeff3(lowfil), 'a')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + w_c(i1,i2,t)*ad1_ext(t-i3-0)
                 dyi1=dyi1 + w_c(i1,i2,t)*ad1_ext(t-i3-1)
                 dyi2=dyi2 + w_c(i1,i2,t)*ad1_ext(t-i3-2)
                 dyi3=dyi3 + w_c(i1,i2,t)*ad1_ext(t-i3-3)
              enddo
              z_c(i1,i2,i3+0)=dyi0
              z_c(i1,i2,i3+1)=dyi1
              z_c(i1,i2,i3+2)=dyi2
              z_c(i1,i2,i3+3)=dyi3
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi=0.0_wp
           !! Get the effective a-filters for the y dimension
           !z0=hgrid*(i3+0)-rxyzConf(3)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + w_c(i1,i2,t)*ad1_ext(t-i3)
           enddo
           z_c(i1,i2,i3)=dyi
        enddo
        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              !! Get the effective b-filters for the z dimension
              !z0=hgrid*(i3+0)-rxyzConf(3)
              !z1=hgrid*(i3+1)-rxyzConf(3)
              !z2=hgrid*(i3+2)-rxyzConf(3)
              !z3=hgrid*(i3+3)-rxyzConf(3)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z1, beff1(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z2, beff2(lowfil), 'b')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z3, beff3(lowfil), 'b')
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + w_f3(t,i1,i2)*bd1_ext(t-i3-0)
                 dyi1=dyi1 + w_f3(t,i1,i2)*bd1_ext(t-i3-1)
                 dyi2=dyi2 + w_f3(t,i1,i2)*bd1_ext(t-i3-2)
                 dyi3=dyi3 + w_f3(t,i1,i2)*bd1_ext(t-i3-3)
              enddo
              z_c(i1,i2,i3+0)=z_c(i1,i2,i3+0)+dyi0
              z_c(i1,i2,i3+1)=z_c(i1,i2,i3+1)+dyi1
              z_c(i1,i2,i3+2)=z_c(i1,i2,i3+2)+dyi2
              z_c(i1,i2,i3+3)=z_c(i1,i2,i3+3)+dyi3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi=0.0_wp
           !! Get the effective b-filters for the y dimension
           !z0=hgrid*(i3+0)-rxyzConf(3)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + w_f3(t,i1,i2)*bd1_ext(t-i3)
           enddo
           z_c(i1,i2,i3)=z_c(i1,i2,i3)+dyi
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              !! Get the effective c-filters for the z dimension
              !z0=hgrid*(i3+0)-rxyzConf(3)
              !z1=hgrid*(i3+1)-rxyzConf(3)
              !z2=hgrid*(i3+2)-rxyzConf(3)
              !z3=hgrid*(i3+3)-rxyzConf(3)
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, ceff0(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z1, ceff1(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z2, ceff2(lowfil), 'c')
              !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z3, ceff3(lowfil), 'c')
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + w_c(i1,i2,t)*cd1_ext(t-i3-0)
                 dyi1=dyi1 + w_c(i1,i2,t)*cd1_ext(t-i3-1)
                 dyi2=dyi2 + w_c(i1,i2,t)*cd1_ext(t-i3-2)
                 dyi3=dyi3 + w_c(i1,i2,t)*cd1_ext(t-i3-3)
              enddo
              z_f(4,i1,i2,i3+0)=dyi0
              z_f(4,i1,i2,i3+1)=dyi1
              z_f(4,i1,i2,i3+2)=dyi2
              z_f(4,i1,i2,i3+3)=dyi3
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi=0.0_wp 
           !! Get the effective c-filters for the z dimension
           !z0=hgrid*(i3+0)-rxyzConf(3)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid,  z0, ceff0(lowfil), 'c')
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + w_c(i1,i2,t)*cd1_ext(t-i3)
           enddo
           z_f(4,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !!!$omp enddo


  
  !  call system_clock(ncount3,ncount_rate,ncount_max)
  !  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2

  !!!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           !! Get the effective filters for the x dimension
           !x0=hgrid*(i1+0)-rxyzConf(1)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, aeff0(lowfil), 'a')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, beff0(lowfil), 'b')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, ceff0(lowfil), 'c')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, x0, eeff0(lowfil), 'e')
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              !t112=t112 + w_f(4,i1+l,i2,i3)*ad1_ext(l) + w_f(5,i1+l,i2,i3)*bd1_ext(l)
              t121=t121 + w_f(2,i1+l,i2,i3)*ad1_ext(l) + w_f(3,i1+l,i2,i3)*bd1_ext(l)
              !t122=t122 + w_f(6,i1+l,i2,i3)*ad1_ext(l) + w_f(7,i1+l,i2,i3)*bd1_ext(l)
              !t212=t212 + w_f(4,i1+l,i2,i3)*cd1_ext(l) + w_f(5,i1+l,i2,i3)*ed1(l)
              t221=t221 + w_f(2,i1+l,i2,i3)*cd1_ext(l) + w_f(3,i1+l,i2,i3)*ed1(l)
              !t222=t222 + w_f(6,i1+l,i2,i3)*cd1_ext(l) + w_f(7,i1+l,i2,i3)*ed1(l)
              !t211=t211 + w_f(1,i1+l,i2,i3)*ed1(l)
           enddo
           x_f(4,i1,i2,i3)=t112
           x_f(2,i1,i2,i3)=t121
           x_f(1,i1,i2,i3)=x_f(1,i1,i2,i3)+t211
           x_f(6,i1,i2,i3)=t122
           x_f(5,i1,i2,i3)=t212
           x_f(3,i1,i2,i3)=t221
           x_f(7,i1,i2,i3)=t222
        enddo
     enddo
  enddo
  !!!$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  !!!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           !! Get the effective filters for the y dimension
           !y0=hgrid*(i2+0)-rxyzConf(2)
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, aeff0(lowfil), 'a')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, beff0(lowfil), 'b')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, ceff0(lowfil), 'c')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, y0, eeff0(lowfil), 'e')
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + w_f(4,i1,i2+l,i3)*ad1_ext(l) + w_f(6,i1,i2+l,i3)*bd1_ext(l)
              t211=t211 + w_f(1,i1,i2+l,i3)*ad1_ext(l) + w_f(3,i1,i2+l,i3)*bd1_ext(l)
              t122=t122 + w_f(4,i1,i2+l,i3)*cd1_ext(l) + w_f(6,i1,i2+l,i3)*ed1(l)
              t212=t212 + w_f(5,i1,i2+l,i3)*ad1_ext(l) + w_f(7,i1,i2+l,i3)*bd1_ext(l)
              t221=t221 + w_f(1,i1,i2+l,i3)*cd1_ext(l) + w_f(3,i1,i2+l,i3)*ed1(l)
              t222=t222 + w_f(5,i1,i2+l,i3)*cd1_ext(l) + w_f(7,i1,i2+l,i3)*ed1(l)
              t121=t121 + w_f(2,i1,i2+l,i3)*ed1(l)
           enddo
           y_f(4,i1,i2,i3)=t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=t211
           y_f(6,i1,i2,i3)=t122
           y_f(5,i1,i2,i3)=t212
           y_f(3,i1,i2,i3)=t221
           y_f(7,i1,i2,i3)=t222
        enddo
     enddo
  enddo
  !!!$omp enddo

  !  call system_clock(ncount5,ncount_rate,ncount_max)
  !  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           !! Get the effective filters for the z dimension
           !z0=hgrid*(i3+0)-rxyzConf(3)
           !!call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, aeff0(lowfil), 'a')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, beff0(lowfil), 'b')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, ceff0(lowfil), 'c')
           !call getEffectiveFilterQuartic(it,potentialPrefac,hgrid, z0, eeff0(lowfil), 'e')
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + w_f(2,i1,i2,i3+l)*ad1_ext(l) + w_f(6,i1,i2,i3+l)*bd1_ext(l)
              t211=t211 + w_f(1,i1,i2,i3+l)*ad1_ext(l) + w_f(5,i1,i2,i3+l)*bd1_ext(l)
              t122=t122 + w_f(2,i1,i2,i3+l)*cd1_ext(l) + w_f(6,i1,i2,i3+l)*ed1(l)
              t212=t212 + w_f(1,i1,i2,i3+l)*cd1_ext(l) + w_f(5,i1,i2,i3+l)*ed1(l)
              t221=t221 + w_f(3,i1,i2,i3+l)*ad1_ext(l) + w_f(7,i1,i2,i3+l)*bd1_ext(l)
              t222=t222 + w_f(3,i1,i2,i3+l)*cd1_ext(l) + w_f(7,i1,i2,i3+l)*ed1(l)
              t112=t112 + w_f(4,i1,i2,i3+l)*ed1(l)
           enddo
           z_f(4,i1,i2,i3)=z_f(4,i1,i2,i3)+t112
           z_f(2,i1,i2,i3)=t121
           z_f(1,i1,i2,i3)=t211
           z_f(6,i1,i2,i3)=t122
           z_f(5,i1,i2,i3)=t212
           z_f(3,i1,i2,i3)=t221
           z_f(7,i1,i2,i3)=t222

        enddo
     enddo
  enddo
  !!!$omp enddo

  !!!$omp end parallel
!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel



END SUBROUTINE createDerivativeBasis




subroutine auxiliary_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, offsetx, offsety, offsetz, &
           hgrid, rxyzConf, potentialPrefac, ibyz_c, ibyz_f, ibxz_c, ibxz_f, ibxy_c, ibxy_f, &
           x_c, x_f, &
           xa_c, xb_c, xc_c, xe_c, ya_c, yb_c, yc_c, ye_c, xa_f, xb_f, xc_f, xe_f, ya_f, yb_f, yc_f, ye_f)
use module_base
implicit none

! Calling arguments
integer,intent(in):: n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, offsetx, offsety, offsetz
real(gp),intent(in):: hgrid
integer,dimension(2,0:n2,0:n3),intent(in):: ibyz_c,ibyz_f
integer,dimension(2,0:n1,0:n3),intent(in):: ibxz_c,ibxz_f
integer,dimension(2,0:n1,0:n2),intent(in):: ibxy_c,ibxy_f
real(wp),dimension(0:n1,0:n2,0:n3), intent(in):: x_c
real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in):: x_f
real(wp),dimension(0:n1,0:n2,0:n3),intent(out):: xa_c, xb_c, xc_c, xe_c, ya_c, yb_c, yc_c, ye_c
real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(out):: xa_f, xb_f, xc_f, xe_f, ya_f, yb_f, yc_f, ye_f
real(8),dimension(3):: rxyzConf
real(8):: potentialPrefac

! Local variables
integer,parameter:: lowfil=-14,lupfil=14, it=1! fake
real(wp),dimension(-3+lowfil:lupfil+3):: aeff0, aeff1, aeff2, aeff3
real(wp),dimension(-3+lowfil:lupfil+3):: beff0, beff1, beff2, beff3
real(wp),dimension(-3+lowfil:lupfil+3):: ceff0, ceff1, ceff2, ceff3
real(wp),dimension(-3+lowfil:lupfil+3):: eeff0, eeff1, eeff2, eeff3
integer:: i1, i2, i3, icur, t, l
real(8):: x0, x1, x2, x3, y0, y1, y2, y3
real(8):: tt0a0, tt0a1, tt0a2, tt0a3
real(8):: tt0b0, tt0b1, tt0b2, tt0b3
real(8):: tt0c0, tt0c1, tt0c2, tt0c3
real(8):: tt0e0, tt0e1, tt0e2, tt0e3
real(8):: tt1a0, tt1a1, tt1a2, tt1a3
real(8):: tt1b0, tt1b1, tt1b2, tt1b3
real(8):: tt1c0, tt1c1, tt1c2, tt1c3
real(8):: tt1e0, tt1e1, tt1e2, tt1e3
real(8):: tt2a0, tt2a1, tt2a2, tt2a3
real(8):: tt2b0, tt2b1, tt2b2, tt2b3
real(8):: tt2c0, tt2c1, tt2c2, tt2c3
real(8):: tt2e0, tt2e1, tt2e2, tt2e3
real(8):: tt3a0, tt3a1, tt3a2, tt3a3
real(8):: tt3b0, tt3b1, tt3b2, tt3b3
real(8):: tt3c0, tt3c1, tt3c2, tt3c3
real(8):: tt3e0, tt3e1, tt3e2, tt3e3
real(8):: tt4a0, tt4a1, tt4a2, tt4a3
real(8):: tt4b0, tt4b1, tt4b2, tt4b3
real(8):: tt4c0, tt4c1, tt4c2, tt4c3
real(8):: tt4e0, tt4e1, tt4e2, tt4e3
real(8):: tt5a0, tt5a1, tt5a2, tt5a3
real(8):: tt5b0, tt5b1, tt5b2, tt5b3
real(8):: tt5c0, tt5c1, tt5c2, tt5c3
real(8):: tt5e0, tt5e1, tt5e2, tt5e3
real(8):: tt6a0, tt6a1, tt6a2, tt6a3
real(8):: tt6b0, tt6b1, tt6b2, tt6b3
real(8):: tt6c0, tt6c1, tt6c2, tt6c3
real(8):: tt6e0, tt6e1, tt6e2, tt6e3
real(8):: tt7a0, tt7a1, tt7a2, tt7a3
real(8):: tt7b0, tt7b1, tt7b2, tt7b3
real(8):: tt7c0, tt7c1, tt7c2, tt7c3
real(8):: tt7e0, tt7e1, tt7e2, tt7e3


aeff0=0.d0 ; aeff1=0.d0 ; aeff2=0.d0 ; aeff3=0.d0
beff0=0.d0 ; beff1=0.d0 ; beff2=0.d0 ; beff3=0.d0
ceff0=0.d0 ; ceff1=0.d0 ; ceff2=0.d0 ; ceff3=0.d0
eeff0=0.d0 ; eeff1=0.d0 ; eeff2=0.d0 ; eeff3=0.d0


  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              tt0a0=0.d0
              tt0a1=0.d0
              tt0a2=0.d0
              tt0a3=0.d0
              tt0b0=0.d0
              tt0b1=0.d0
              tt0b2=0.d0
              tt0b3=0.d0
              tt0c0=0.d0
              tt0c1=0.d0
              tt0c2=0.d0
              tt0c3=0.d0
              tt0e0=0.d0
              tt0e1=0.d0
              tt0e2=0.d0
              tt0e3=0.d0
              ! Get the effective a-filters for the x dimension
              x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
              x1=hgrid*(i1+offsetx+1)-rxyzConf(1)
              x2=hgrid*(i1+offsetx+2)-rxyzConf(1)
              x3=hgrid*(i1+offsetx+3)-rxyzConf(1)

              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, aeff1(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, aeff2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, aeff3(lowfil), 'a')

              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, beff1(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, beff2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, beff3(lowfil), 'b')

              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, ceff1(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, ceff2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, ceff3(lowfil), 'c')

              call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0(lowfil), 'e')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x1, eeff1(lowfil), 'e')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x2, eeff2(lowfil), 'e')
              call getFilterQuadratic(it, potentialPrefac, hgrid, x3, eeff3(lowfil), 'e')
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))

                 ! sss coefficients
                 tt0a0=tt0a0 + x_c(t,i2,i3)*aeff0(t-i1-0)
                 tt0a1=tt0a1 + x_c(t,i2,i3)*aeff1(t-i1-1)
                 tt0a2=tt0a2 + x_c(t,i2,i3)*aeff2(t-i1-2)
                 tt0a3=tt0a3 + x_c(t,i2,i3)*aeff3(t-i1-3)

                 tt0b0=tt0b0 + x_c(t,i2,i3)*beff0(t-i1-0)
                 tt0b1=tt0b1 + x_c(t,i2,i3)*beff1(t-i1-1)
                 tt0b2=tt0b2 + x_c(t,i2,i3)*beff2(t-i1-2)
                 tt0b3=tt0b3 + x_c(t,i2,i3)*beff3(t-i1-3)

                 tt0c0=tt0c0 + x_c(t,i2,i3)*ceff0(t-i1-0)
                 tt0c1=tt0c1 + x_c(t,i2,i3)*ceff1(t-i1-1)
                 tt0c2=tt0c2 + x_c(t,i2,i3)*ceff2(t-i1-2)
                 tt0c3=tt0c3 + x_c(t,i2,i3)*ceff3(t-i1-3)

                 tt0e0=tt0e0 + x_c(t,i2,i3)*eeff0(t-i1-0)
                 tt0e1=tt0e1 + x_c(t,i2,i3)*eeff1(t-i1-1)
                 tt0e2=tt0e2 + x_c(t,i2,i3)*eeff2(t-i1-2)
                 tt0e3=tt0e3 + x_c(t,i2,i3)*eeff3(t-i1-3)

              enddo

              xa_c(i1+0,i2,i3)=tt0a0
              xa_c(i1+1,i2,i3)=tt0a1
              xa_c(i1+2,i2,i3)=tt0a2
              xa_c(i1+3,i2,i3)=tt0a3

              xb_c(i1+0,i2,i3)=tt0b0
              xb_c(i1+1,i2,i3)=tt0b1
              xb_c(i1+2,i2,i3)=tt0b2
              xb_c(i1+3,i2,i3)=tt0b3

              xc_c(i1+0,i2,i3)=tt0c0
              xc_c(i1+1,i2,i3)=tt0c1
              xc_c(i1+2,i2,i3)=tt0c2
              xc_c(i1+3,i2,i3)=tt0c3

              xe_c(i1+0,i2,i3)=tt0e0
              xe_c(i1+1,i2,i3)=tt0e1
              xe_c(i1+2,i2,i3)=tt0e2
              xe_c(i1+3,i2,i3)=tt0e3
           end do
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        end if

        do i1=icur,ibyz_c(2,i2,i3)
           tt0a0=0.d0
           tt0b0=0.d0
           tt0c0=0.d0
           tt0e0=0.d0
           ! Get the effective a-filters for the x dimension
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0(lowfil), 'e')
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              ! sss coefficients
              tt0a0=tt0a0 + x_c(t,i2,i3)*aeff0(t-i1)
              tt0b0=tt0b0 + x_c(t,i2,i3)*beff0(t-i1)
              tt0c0=tt0c0 + x_c(t,i2,i3)*ceff0(t-i1)
              tt0e0=tt0e0 + x_c(t,i2,i3)*eeff0(t-i1)
           enddo
           xa_c(i1,i2,i3)=tt0a0
           xb_c(i1,i2,i3)=tt0b0
           xc_c(i1,i2,i3)=tt0c0
           xe_c(i1,i2,i3)=tt0e0
        end do
     end do
  end do
  


  ! fine part
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           tt1a0=0.d0
           tt1b0=0.d0
           tt1c0=0.d0
           tt1e0=0.d0
           tt2a0=0.d0
           tt2b0=0.d0
           tt2c0=0.d0
           tt2e0=0.d0
           tt3a0=0.d0
           tt3b0=0.d0
           tt3c0=0.d0
           tt3e0=0.d0
           tt4a0=0.d0
           tt4b0=0.d0
           tt4c0=0.d0
           tt4e0=0.d0
           tt5a0=0.d0
           tt5b0=0.d0
           tt5c0=0.d0
           tt5e0=0.d0
           tt6a0=0.d0
           tt6b0=0.d0
           tt6c0=0.d0
           tt6e0=0.d0
           tt7a0=0.d0
           tt7b0=0.d0
           tt7c0=0.d0
           tt7e0=0.d0
           x0=hgrid*(i1+offsetx+0)-rxyzConf(1)
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, x0, eeff0(lowfil), 'e')
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              ! dss coefficients
              tt1a0=tt1a0 + x_f(1,i1+l,i2,i3)*aeff0(l)
              tt1b0=tt1b0 + x_f(1,i1+l,i2,i3)*beff0(l)
              tt1c0=tt1c0 + x_f(1,i1+l,i2,i3)*ceff0(l)
              tt1e0=tt1e0 + x_f(1,i1+l,i2,i3)*eeff0(l)
              ! sds coefficients
              tt2a0=tt2a0 + x_f(2,i1+l,i2,i3)*aeff0(l)
              tt2b0=tt2b0 + x_f(2,i1+l,i2,i3)*beff0(l)
              tt2c0=tt2c0 + x_f(2,i1+l,i2,i3)*ceff0(l)
              tt2e0=tt2e0 + x_f(2,i1+l,i2,i3)*eeff3(l)
              ! dds coefficients
              tt3a0=tt3a0 + x_f(3,i1+l,i2,i3)*aeff0(l)
              tt3b0=tt3b0 + x_f(3,i1+l,i2,i3)*beff0(l)
              tt3c0=tt3c0 + x_f(3,i1+l,i2,i3)*ceff0(l)
              tt3e0=tt3e0 + x_f(3,i1+l,i2,i3)*eeff0(l)
              ! ssd coefficients
              tt4a0=tt4a0 + x_f(4,i1+l,i2,i3)*aeff0(l)
              tt4b0=tt4b0 + x_f(4,i1+l,i2,i3)*beff0(l)
              tt4c0=tt4c0 + x_f(4,i1+l,i2,i3)*ceff0(l)
              tt4e0=tt4e0 + x_f(4,i1+l,i2,i3)*eeff0(l)
              ! dsd coefficients
              tt5a0=tt5a0 + x_f(5,i1+l,i2,i3)*aeff0(l)
              tt5b0=tt5b0 + x_f(5,i1+l,i2,i3)*beff0(l)
              tt5c0=tt5c0 + x_f(5,i1+l,i2,i3)*ceff0(l)
              tt5e0=tt5e0 + x_f(5,i1+l,i2,i3)*eeff0(l)
              ! sdd coefficients
              tt6a0=tt6a0 + x_f(6,i1+l,i2,i3)*aeff0(l)
              tt6b0=tt6b0 + x_f(6,i1+l,i2,i3)*beff0(l)
              tt6c0=tt6c0 + x_f(6,i1+l,i2,i3)*ceff0(l)
              tt6e0=tt6e0 + x_f(6,i1+l,i2,i3)*eeff0(l)
              ! ddd coefficients
              tt7a0=tt7a0 + x_f(7,i1+l,i2,i3)*aeff0(l)
              tt7b0=tt7b0 + x_f(7,i1+l,i2,i3)*beff0(l)
              tt7c0=tt7c0 + x_f(7,i1+l,i2,i3)*ceff0(l)
              tt7e0=tt7e0 + x_f(7,i1+l,i2,i3)*eeff0(l)
           end do
           ! dss coefficients
           xa_f(1,i1,i2,i3)=tt1a0
           xb_f(1,i1,i2,i3)=tt1b0
           xc_f(1,i1,i2,i3)=tt1c0
           xe_f(1,i1,i2,i3)=tt1e0
           ! sds coefficients
           xa_f(2,i1,i2,i3)=tt2a0
           xb_f(2,i1,i2,i3)=tt2b0
           xc_f(2,i1,i2,i3)=tt2c0
           xe_f(2,i1,i2,i3)=tt2e0
           ! dds coefficients
           xa_f(3,i1,i2,i3)=tt3a0
           xb_f(3,i1,i2,i3)=tt3b0
           xc_f(3,i1,i2,i3)=tt3c0
           xe_f(3,i1,i2,i3)=tt3e0
           ! ssd coefficients
           xa_f(4,i1,i2,i3)=tt4a0
           xb_f(4,i1,i2,i3)=tt4b0
           xc_f(4,i1,i2,i3)=tt4c0
           xe_f(4,i1,i2,i3)=tt4e0
           ! dsd coefficients
           xa_f(5,i1,i2,i3)=tt5a0
           xb_f(5,i1,i2,i3)=tt5b0
           xc_f(5,i1,i2,i3)=tt5c0
           xe_f(5,i1,i2,i3)=tt5e0
           ! sdd coefficients
           xa_f(6,i1,i2,i3)=tt6a0
           xb_f(6,i1,i2,i3)=tt6b0
           xc_f(6,i1,i2,i3)=tt6c0
           xe_f(6,i1,i2,i3)=tt6e0
           ! sdd coefficients
           xa_f(7,i1,i2,i3)=tt7a0
           xb_f(7,i1,i2,i3)=tt7b0
           xc_f(7,i1,i2,i3)=tt7c0
           xe_f(7,i1,i2,i3)=tt7e0
        end do
     end do
  end do




  ! y direction
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              tt0a0=0.d0
              tt0a1=0.d0
              tt0a2=0.d0
              tt0a3=0.d0
              tt0b0=0.d0
              tt0b1=0.d0
              tt0b2=0.d0
              tt0b3=0.d0
              tt0c0=0.d0
              tt0c1=0.d0
              tt0c2=0.d0
              tt0c3=0.d0
              tt0e0=0.d0
              tt0e1=0.d0
              tt0e2=0.d0
              tt0e3=0.d0
              ! Get the effective a-filters for the y dimension
              y0=hgrid*(i2+offsety+0)-rxyzConf(2)
              y1=hgrid*(i2+offsety+1)-rxyzConf(2)
              y2=hgrid*(i2+offsety+2)-rxyzConf(2)
              y3=hgrid*(i2+offsety+3)-rxyzConf(2)

              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, aeff1(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, aeff2(lowfil), 'a')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, aeff3(lowfil), 'a')

              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, beff1(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, beff2(lowfil), 'b')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, beff3(lowfil), 'b')

              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, ceff1(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, ceff2(lowfil), 'c')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, ceff3(lowfil), 'c')

              call getFilterQuadratic(it, potentialPrefac, hgrid, y0, eeff0(lowfil), 'e')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y1, eeff1(lowfil), 'e')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y2, eeff2(lowfil), 'e')
              call getFilterQuadratic(it, potentialPrefac, hgrid, y3, eeff3(lowfil), 'e')

              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 ! sss coefficients
                 tt0a0=tt0a0 + x_c(i1,t,i3)*aeff0(t-i2-0)
                 tt0a1=tt0a1 + x_c(i1,t,i3)*aeff1(t-i2-1)
                 tt0a2=tt0a2 + x_c(i1,t,i3)*aeff2(t-i2-2)
                 tt0a3=tt0a3 + x_c(i1,t,i3)*aeff3(t-i2-3)

                 tt0b0=tt0b0 + x_c(i1,t,i3)*beff0(t-i2-0)
                 tt0b1=tt0b1 + x_c(i1,t,i3)*beff1(t-i2-1)
                 tt0b2=tt0b2 + x_c(i1,t,i3)*beff2(t-i2-2)
                 tt0b3=tt0b3 + x_c(i1,t,i3)*beff3(t-i2-3)

                 tt0c0=tt0c0 + x_c(i1,t,i3)*ceff0(t-i2-0)
                 tt0c1=tt0c1 + x_c(i1,t,i3)*ceff1(t-i2-1)
                 tt0c2=tt0c2 + x_c(i1,t,i3)*ceff2(t-i2-2)
                 tt0c3=tt0c3 + x_c(i1,t,i3)*ceff3(t-i2-3)

                 tt0e0=tt0e0 + x_c(i1,t,i3)*eeff0(t-i2-0)
                 tt0e1=tt0e1 + x_c(i1,t,i3)*eeff1(t-i2-1)
                 tt0e2=tt0e2 + x_c(i1,t,i3)*eeff2(t-i2-2)
                 tt0e3=tt0e3 + x_c(i1,t,i3)*eeff3(t-i2-3)

              enddo
              ya_c(i1,i2+0,i3)=tt0a0
              ya_c(i1,i2+1,i3)=tt0a1
              ya_c(i1,i2+2,i3)=tt0a2
              ya_c(i1,i2+3,i3)=tt0a3
                        
              yb_c(i1,i2+0,i3)=tt0b0
              yb_c(i1,i2+1,i3)=tt0b1
              yb_c(i1,i2+2,i3)=tt0b2
              yb_c(i1,i2+3,i3)=tt0b3
                        
              yc_c(i1,i2+0,i3)=tt0c0
              yc_c(i1,i2+1,i3)=tt0c1
              yc_c(i1,i2+2,i3)=tt0c2
              yc_c(i1,i2+3,i3)=tt0c3
                        
              ye_c(i1,i2+0,i3)=tt0e0
              ye_c(i1,i2+1,i3)=tt0e1
              ye_c(i1,i2+2,i3)=tt0e2
              ye_c(i1,i2+3,i3)=tt0e3
           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           ! Get the effective a-filters for the y dimension
           tt0a0=0.d0
           tt0b0=0.d0
           tt0c0=0.d0
           tt0e0=0.d0
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, eeff0(lowfil), 'e')
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              ! sss coefficients
              tt0a0=tt0a0 + x_c(i1,t,i3)*aeff0(t-i2-0)
              tt0b0=tt0b0 + x_c(i1,t,i3)*beff0(t-i2-0)
              tt0c0=tt0c0 + x_c(i1,t,i3)*ceff0(t-i2-0)
              tt0e0=tt0e0 + x_c(i1,t,i3)*eeff0(t-i2-0)
           end do
           ya_c(i1,i2+0,i3)=tt0a0
           yb_c(i1,i2+0,i3)=tt0b0
           yc_c(i1,i2+0,i3)=tt0c0
           ye_c(i1,i2+0,i3)=tt0e0
        end do
     end do
  end do



  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           ! Get the effective filters for the y dimension
           tt1a0=0.d0
           tt1b0=0.d0
           tt1c0=0.d0
           tt1e0=0.d0
           tt2a0=0.d0
           tt2b0=0.d0
           tt2c0=0.d0
           tt2e0=0.d0
           tt3a0=0.d0
           tt3b0=0.d0
           tt3c0=0.d0
           tt3e0=0.d0
           tt4a0=0.d0
           tt4b0=0.d0
           tt4c0=0.d0
           tt4e0=0.d0
           tt5a0=0.d0
           tt5b0=0.d0
           tt5c0=0.d0
           tt5e0=0.d0
           tt6a0=0.d0
           tt6b0=0.d0
           tt6c0=0.d0
           tt6e0=0.d0
           tt7a0=0.d0
           tt7b0=0.d0
           tt7c0=0.d0
           tt7e0=0.d0
           y0=hgrid*(i2+offsety+0)-rxyzConf(2)
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, aeff0(lowfil), 'a')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, beff0(lowfil), 'b')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, ceff0(lowfil), 'c')
           call getFilterQuadratic(it, potentialPrefac, hgrid, y0, eeff0(lowfil), 'e')
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              ! dss coefficients
              tt1a0=tt1a0 + x_f(1,i1,i2+l,i3)*aeff0(l)
              tt1b0=tt1b0 + x_f(1,i1,i2+l,i3)*beff0(l)
              tt1c0=tt1c0 + x_f(1,i1,i2+l,i3)*ceff0(l)
              tt1e0=tt1e0 + x_f(1,i1,i2+l,i3)*eeff0(l)
              ! sds coefficients
              tt2a0=tt2a0 + x_f(2,i1,i2+l,i3)*aeff0(l)
              tt2b0=tt2b0 + x_f(2,i1,i2+l,i3)*beff0(l)
              tt2c0=tt2c0 + x_f(2,i1,i2+l,i3)*ceff0(l)
              tt2e0=tt2e0 + x_f(2,i1,i2+l,i3)*eeff3(l)
              ! dds coefficients
              tt3a0=tt3a0 + x_f(3,i1,i2+l,i3)*aeff0(l)
              tt3b0=tt3b0 + x_f(3,i1,i2+l,i3)*beff0(l)
              tt3c0=tt3c0 + x_f(3,i1,i2+l,i3)*ceff0(l)
              tt3e0=tt3e0 + x_f(3,i1,i2+l,i3)*eeff0(l)
              ! ssd coefficients
              tt4a0=tt4a0 + x_f(4,i1,i2+l,i3)*aeff0(l)
              tt4b0=tt4b0 + x_f(4,i1,i2+l,i3)*beff0(l)
              tt4c0=tt4c0 + x_f(4,i1,i2+l,i3)*ceff0(l)
              tt4e0=tt4e0 + x_f(4,i1,i2+l,i3)*eeff0(l)
              ! dsd coefficients
              tt5a0=tt5a0 + x_f(5,i1,i2+l,i3)*aeff0(l)
              tt5b0=tt5b0 + x_f(5,i1,i2+l,i3)*beff0(l)
              tt5c0=tt5c0 + x_f(5,i1,i2+l,i3)*ceff0(l)
              tt5e0=tt5e0 + x_f(5,i1,i2+l,i3)*eeff0(l)
              ! sdd coefficients
              tt6a0=tt6a0 + x_f(6,i1,i2+l,i3)*aeff0(l)
              tt6b0=tt6b0 + x_f(6,i1,i2+l,i3)*beff0(l)
              tt6c0=tt6c0 + x_f(6,i1,i2+l,i3)*ceff0(l)
              tt6e0=tt6e0 + x_f(6,i1,i2+l,i3)*eeff0(l)
              ! ddd coefficients
              tt7a0=tt7a0 + x_f(7,i1,i2+l,i3)*aeff0(l)
              tt7b0=tt7b0 + x_f(7,i1,i2+l,i3)*beff0(l)
              tt7c0=tt7c0 + x_f(7,i1,i2+l,i3)*ceff0(l)
              tt7e0=tt7e0 + x_f(7,i1,i2+l,i3)*eeff0(l)
           enddo
           ! dss coefficients
           ya_f(1,i1,i2,i3)=tt1a0
           yb_f(1,i1,i2,i3)=tt1b0
           yc_f(1,i1,i2,i3)=tt1c0
           ye_f(1,i1,i2,i3)=tt1e0
           ! sds coefficients
           ya_f(2,i1,i2,i3)=tt2a0
           yb_f(2,i1,i2,i3)=tt2b0
           yc_f(2,i1,i2,i3)=tt2c0
           ye_f(2,i1,i2,i3)=tt2e0
           ! dds coefficients
           ya_f(3,i1,i2,i3)=tt3a0
           yb_f(3,i1,i2,i3)=tt3b0
           yc_f(3,i1,i2,i3)=tt3c0
           ye_f(3,i1,i2,i3)=tt3e0
           ! ssd coefficients
           ya_f(4,i1,i2,i3)=tt4a0
           yb_f(4,i1,i2,i3)=tt4b0
           yc_f(4,i1,i2,i3)=tt4c0
           ye_f(4,i1,i2,i3)=tt4e0
           ! dsd coefficients
           ya_f(5,i1,i2,i3)=tt5a0
           yb_f(5,i1,i2,i3)=tt5b0
           yc_f(5,i1,i2,i3)=tt5c0
           ye_f(5,i1,i2,i3)=tt5e0
           ! sdd coefficients
           ya_f(6,i1,i2,i3)=tt6a0
           yb_f(6,i1,i2,i3)=tt6b0
           yc_f(6,i1,i2,i3)=tt6c0
           ye_f(6,i1,i2,i3)=tt6e0
           ! sdd coefficients
           ya_f(7,i1,i2,i3)=tt7a0
           yb_f(7,i1,i2,i3)=tt7b0
           yc_f(7,i1,i2,i3)=tt7c0
           ye_f(7,i1,i2,i3)=tt7e0
        end do
     end do
  end do

end subroutine auxiliary_convolutions
