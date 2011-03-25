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


prefac1=-.5d0/hgrid**2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac
fac2=parabPrefac*hgrid
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


!$omp parallel default(private) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
  !$omp do  
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
  !$omp enddo
  
  !  call system_clock(ncount1,ncount_rate,ncount_max)
  !  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  !$omp do
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
  !$omp enddo


  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2

  !$omp do
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
  !$omp enddo


  
  !  call system_clock(ncount3,ncount_rate,ncount_max)
  !  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2

  !$omp do
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
  !$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  !$omp do
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
  !$omp enddo

  !  call system_clock(ncount5,ncount_rate,ncount_max)
  !  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  !$omp do
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
  !$omp enddo

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




END SUBROUTINE ConvolkineticQuartic
