  !   y = (kinetic energy operator)x + (cprec*I)x 
! One of the most CPU intensive routines
subroutine Convolkinetic(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3)
  use module_base
  implicit none
!dee
  integer :: iend_test,count_rate_test,count_max_test,istart_test

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
  integer(8) :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6,ncount0,ncnt1
  integer(8) :: clock0,clock1,clock2
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(kind=8) :: tel
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
call system_clock(istart_test,count_rate_test,count_max_test)


call system_clock(ncount0,ncount_rate,ncount_max)

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
call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

!  call system_clock(ncount6,ncount_rate,ncount_max)
!  tel=dble(ncount6-ncount5)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:z',tel,1.d-6*nflop3/tel

!  tel=dble(ncount6-ncount0)/dble(ncount_rate)
!  write(99,'(a40,1x,e10.3,1x,f6.1)') 'ALL   PART',  & 
!  tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

END SUBROUTINE Convolkinetic

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
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  logical :: firstcall=.true. 
  integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6,ncount0
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211,ekin
  real(kind=8) :: tel
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

  call system_clock(ncount0,ncount_rate,ncount_max)

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


! **********************************************************************************************
! **********************************************************************************************
! **********************************************************************************************
! **********************************************************************************************
! **********************************************************************************************
! **********************************************************************************************









!!module filterModule
!!!
!!! Purpose:
!!! ========
!!!   This module contains various filters for the Daubechies 16 wavelet family.
!!!   The meaning is the following:
!!!     a_i=<phi|Op|phi_i>
!!!     b_i=<psi|Op|phi_i>
!!!     c_i=<phi|Op|psi_i>
!!!     e_i=<psi|Op|psi_i>
!!!   There are three different operators:
!!!     Op = second derivative -> a, b, c, e
!!!     Op = x -> a1, b1, c1, e1
!!!     Op = x^2 -> a2, b2, c2, e2
!!!
!!implicit none
!!
!!! The lower and upper bound for the filters.
!!integer,parameter:: lb=-14, ub=14
!!
!!real(8),dimension(lb:ub),parameter:: a = [ -6.92447494063920016E-18, &
!!                                            2.70800493626319433E-13, &
!!                                           -5.81387983028254035E-11, &
!!                                           -1.05857055496741463E-08, &
!!                                           -3.72307630473692755E-07, &
!!                                            2.09042349529203676E-06, &
!!                                           -2.39822852450759980E-05, &
!!                                            4.51679202875022361E-04, &
!!                                           -4.09765689342633799E-03, &
!!                                            2.20702918848225543E-02, &
!!                                           -8.22663999742123397E-02, &
!!                                            2.37178058215380572E-01, &
!!                                           -6.15614146557006969E-01, &
!!                                            2.21914659389111657E+00, &
!!                                           -3.55369228991319019E+00, &
!!                                            2.21914659389111657E+00, &
!!                                           -6.15614146557006969E-01, &
!!                                            2.37178058215380572E-01, &
!!                                           -8.22663999742123397E-02, &
!!                                            2.20702918848225543E-02, &
!!                                           -4.09765689342633799E-03, &
!!                                            4.51679202875022361E-04, &
!!                                           -2.39822852450759980E-05, &
!!                                            2.09042349529203676E-06, &
!!                                           -3.72307630473692755E-07, &
!!                                           -1.05857055496741463E-08, &
!!                                           -5.81387983028254035E-11, &
!!                                            2.70800493626319433E-13, &
!!                                           -6.92447494063920016E-18  ]
!!
!!real(8),dimension(lb:ub),parameter:: b = [  1.23926296291889880E-17, &
!!                                           -4.84665690596158973E-13, &
!!                                            1.04757345962781863E-10, &
!!                                            1.87909235155149931E-08, &
!!                                            6.39191285137934231E-07, &
!!                                           -4.62561497463184309E-06, &
!!                                            5.07580794728972886E-05, &
!!                                           -9.34236233048086760E-04, &
!!                                            8.08720294118447924E-03, &
!!                                           -5.65773697332875719E-02, &
!!                                            2.33704906317512895E-01, &
!!                                           -5.32982239628004284E-01, &
!!                                            6.34570352678925254E-01, &
!!                                           -2.17161175051370531E-01, &
!!                                           -4.20816552937241239E-01, &
!!                                            6.57030740071213981E-01, &
!!                                           -4.24980509437801468E-01, &
!!                                            1.43273293525107509E-01, &
!!                                           -2.84696165863973540E-02, &
!!                                            7.13761218453631829E-03, &
!!                                           -2.13562918387979854E-03, &
!!                                            2.17107954846461359E-04, &
!!                                           -1.13345672451682042E-05, &
!!                                            8.79697040552863148E-07, &
!!                                           -2.16568306292140375E-07, &
!!                                           -5.96264938781402246E-09, &
!!                                           -3.22647023140105305E-11, &
!!                                            1.51306165608661538E-13, &
!!                                           -3.86910241314765619E-18  ]
!!
!!real(8),dimension(lb:ub),parameter:: c = [ -3.86910241314765619E-18, &
!!                                            1.51306165608661538E-13, &
!!                                           -3.22647023140105240E-11, &
!!                                           -5.96264938781402246E-09, &
!!                                           -2.16568306292140402E-07, &
!!                                            8.79697040552863466E-07, &
!!                                           -1.13345672451681957E-05, &
!!                                            2.17107954846461495E-04, &
!!                                           -2.13562918387979854E-03, &
!!                                            7.13761218453631743E-03, &
!!                                           -2.84696165863973470E-02, &
!!                                            1.43273293525107648E-01, &
!!                                           -4.24980509437801190E-01, &
!!                                            6.57030740071214536E-01, &
!!                                           -4.20816552937243626E-01, &
!!                                           -2.17161175051370559E-01, &
!!                                            6.34570352678925476E-01, &
!!                                           -5.32982239628004395E-01, &
!!                                            2.33704906317513006E-01, &
!!                                           -5.65773697332875719E-02, &
!!                                            8.08720294118448098E-03, &
!!                                           -9.34236233048086760E-04, &
!!                                            5.07580794728972886E-05, &
!!                                           -4.62561497463184224E-06, &
!!                                            6.39191285137934125E-07, &
!!                                            1.87909235155149897E-08, &
!!                                            1.04757345962781863E-10, &
!!                                           -4.84665690596158973E-13, &
!!                                            1.23926296291889880E-17  ]
!!
!!real(8),dimension(lb:ub),parameter:: e = [  6.92447494063920016E-18, &
!!                                           -2.70800493626319433E-13, &
!!                                            5.81387983028254035E-11, &
!!                                            1.05857055496741463E-08, &
!!                                            3.72307630473692703E-07, &
!!                                           -2.09042349529203803E-06, &
!!                                            2.39822852450760014E-05, &
!!                                           -4.51679202875077818E-04, &
!!                                            4.09765642831595028E-03, &
!!                                           -2.20732703458663519E-02, &
!!                                            8.20745416922517429E-02, &
!!                                           -2.69959313362790998E-01, &
!!                                           -4.25170532366915405E-02, &
!!                                           -7.14405976634717543E+00, &
!!                                           -2.48758460293923314E+01, &
!!                                           -7.14405976634717366E+00, &
!!                                           -4.25170532366917972E-02, &
!!                                           -2.69959313362791276E-01, &
!!                                            8.20745416922516735E-02, &
!!                                           -2.20732703458663554E-02, &
!!                                            4.09765642831595115E-03, &
!!                                           -4.51679202875077709E-04, &
!!                                            2.39822852450760014E-05, &
!!                                           -2.09042349529203592E-06, &
!!                                            3.72307630473692544E-07, &
!!                                            1.05857055496741480E-08, &
!!                                            5.81387983028254099E-11, &
!!                                           -2.70800493626319433E-13, &
!!                                            6.92447494063920016E-18  ]
!!
!!real(8),dimension(lb:ub),parameter:: a1 = [ -4.27782405260083166E-22, &
!!                                             1.33836789073155300E-16, &
!!                                            -2.23077383438816909E-13, &
!!                                            -4.18701072588700730E-11, &
!!                                             3.84061420317539110E-10, &
!!                                             6.90677058415406752E-08, &
!!                                            -1.50667718130096206E-06, &
!!                                             1.43009564059536508E-05, &
!!                                            -7.38484428527907765E-05, &
!!                                             1.17437400174429801E-04, &
!!                                             4.36133077583543763E-03, &
!!                                            -2.15712170526503010E-02, &
!!                                             2.26678480716213104E-02, &
!!                                             9.18613942126320498E-02, &
!!                                            -3.43336140248461028E-02, &
!!                                             9.18613942126320082E-02, &
!!                                             2.26678480716213208E-02, &
!!                                            -2.15712170526503114E-02, &
!!                                             4.36133077583543936E-03, &
!!                                             1.17437400174429598E-04, &
!!                                            -7.38484428527907901E-05, &
!!                                             1.43009564059536593E-05, &
!!                                            -1.50667718130096206E-06, &
!!                                             6.90677058415406752E-08, &
!!                                             3.84061420317539886E-10, &
!!                                            -4.18701072588700924E-11, &
!!                                            -2.23077383438816909E-13, &
!!                                             1.33836789073155398E-16, &
!!                                            -4.27782405260083542E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: b1 = [  7.65595796896980290E-22, &
!!                                            -2.39526824360595103E-16, &
!!                                             3.99587762949889356E-13, &
!!                                             7.43506270483382512E-11, &
!!                                            -7.94750699421454253E-10, &
!!                                            -1.22260456716396509E-07, &
!!                                             2.87022345454858366E-06, &
!!                                              -3.01357019626700847E-05, &
!!                                               1.87517856451687890E-04, &
!!                                              -5.06159008672606094E-04, &
!!                                              -7.78069990955502968E-03, &
!!                                               5.60208710002205595E-02, &
!!                                              -1.59595684370184177E-01, &
!!                                               2.45149541999230086E-01, &
!!                                              -2.14819128526062691E-01, &
!!                                               9.90323036048479100E-02, &
!!                                              -1.30155486015622029E-02, &
!!                                              -6.91564877941363773E-03, &
!!                                               2.28647193534682299E-03, &
!!                                               3.94202903958867657E-06, &
!!                                              -2.62421504146951582E-05, &
!!                                               6.51532505438812717E-06, &
!!                                            -7.86922170343662117E-07, &
!!                                             3.89648346594732955E-08, &
!!                                             1.80799139957949220E-10, &
!!                                            -2.35773086204665242E-11, &
!!                                            -1.24537049640247525E-13, &
!!                                             7.47819629673571264E-17, &
!!                                            -2.39026633886714879E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: c1 = [ -2.39026633886714738E-22, &
!!                                             7.47819629673570894E-17, &
!!                                            -1.24537049640247500E-13, &
!!                                            -2.35773086204665145E-11, &
!!                                             1.80799139957948548E-10, &
!!                                             3.89648346594732955E-08, &
!!                                            -7.86922170343662328E-07, &
!!                                               6.51532505439451803E-06, &
!!                                              -2.62421504147356803E-05, &
!!                                               3.94202903959028594E-06, &
!!                                               2.28647193533762939E-03, &
!!                                              -6.91564877941468897E-03, &
!!                                              -1.30155486015782716E-02, &
!!                                               9.90323036046980160E-02, &
!!                                              -2.14819128525844588E-01, &
!!                                               2.45149542000502374E-01, &
!!                                              -1.59595684369265078E-01, &
!!                                               5.60208710003036389E-02, &
!!                                              -7.78069990956280558E-03, &
!!                                              -5.06159008670565192E-04, &
!!                                               1.87517856454080128E-04, &
!!                                              -3.01357019626822616E-05, &
!!                                             2.87022345454858366E-06, &
!!                                            -1.22260456716396615E-07, &
!!                                            -7.94750699421456321E-10, &
!!                                             7.43506270483382900E-11, &
!!                                             3.99587762949889406E-13, &
!!                                            -2.39526824360595251E-16, &
!!                                             7.65595796896980854E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: e1 = [  4.27782405260083260E-22, &
!!                                            -1.33836789073155324E-16, &
!!                                             2.23077383438816859E-13, &
!!                                             4.18701072588700601E-11, &
!!                                            -3.84061420317538800E-10, &
!!                                            -6.90677058415406884E-08, &
!!                                             1.50667718130096164E-06, &
!!                                              -1.43700241117742952E-05, &
!!                                               7.38480585687695598E-05, &
!!                                              -1.17436974242560681E-04, &
!!                                              -4.36283745276811082E-03, &
!!                                               2.14973686097599537E-02, &
!!                                              -1.83065172959537246E-02, &
!!                                              -6.91935461408035579E-02, &
!!                                              -7.64932673832824859E-08, &
!!                                              -6.91935461404422220E-02, &
!!                                              -1.83065172956272636E-02, &
!!                                               2.14973686097812283E-02, &
!!                                              -4.36283745279333023E-03, &
!!                                              -1.17436974242791575E-04, &
!!                                               7.38480585682706913E-05, &
!!                                              -1.43700241117716610E-05, &
!!                                             1.50667718130096227E-06, &
!!                                            -6.90677058415407281E-08, &
!!                                            -3.84061420317539989E-10, &
!!                                             4.18701072588700859E-11, &
!!                                             2.23077383438816859E-13, &
!!                                            -1.33836789073155374E-16, &
!!                                             4.27782405260083542E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: a2 = [  5.57613113797079647E-21, &
!!                                            -1.61540205057621221E-15, &
!!                                             2.44941827237665777E-12, &
!!                                             4.24629970074974494E-10, &
!!                                            -3.50622375577813727E-09, &
!!                                            -5.49912783753825352E-07, &
!!                                             1.07758069134433754E-05, &
!!                                            -8.87797723407857120E-05, &
!!                                             3.82019672542810440E-04, &
!!                                            -4.44838660769164562E-04, &
!!                                            -1.42199024558067322E-02, &
!!                                             3.99360917508602142E-02, &
!!                                             3.36528867483139038E-02, &
!!                                            -2.19314932823181152E-01, &
!!                                             1.65585339069366455E-01, &
!!                                            -3.55921462178230286E-02, &
!!                                             1.24324277043342590E-01, &
!!                                            -8.94912108778953552E-02, &
!!                                             2.06707436591386795E-02, &
!!                                             7.29535357095301151E-04, &
!!                                            -5.04161638673394918E-04, &
!!                                             1.11433619167655706E-04, &
!!                                            -1.33310277306009084E-05, &
!!                                             6.93305878485261928E-07, &
!!                                             4.17500478633314742E-09, &
!!                                            -4.96512386760628033E-10, &
!!                                            -2.90443884221058823E-12, &
!!                                             1.86435445701578929E-15, &
!!                                            -6.40177613495305920E-21  ]
!!
!!real(8),dimension(lb:ub),parameter:: b2 = [ -9.97951874177877463E-21, &
!!                                             2.89107560457921374E-15, &
!!                                            -4.38790578274431452E-12, &
!!                                            -7.53543776369184845E-10, &
!!                                             7.36262990148872225E-09, &
!!                                             9.71517587727683347E-07, &
!!                                            -2.06572685311829861E-05, &
!!                                               2.22835909285148021E-04, &
!!                                              -1.23566563086296836E-03, &
!!                                               2.81931936447991438E-03, &
!!                                               3.56570289060501519E-02, &
!!                                              -1.99870728151716748E-01, &
!!                                               4.09070814741372302E-01, &
!!                                              -3.83487892570934019E-01, &
!!                                               1.22706165091227781E-01, &
!!                                               4.07715967198992724E-02, &
!!                                              -1.73643886666752928E-02, &
!!                                              -1.68596647369580570E-02, &
!!                                               7.63552892817706119E-03, &
!!                                               5.38796919741555509E-05, &
!!                                              -1.33743588175554292E-04, &
!!                                               4.45539025400963482E-05, &
!!                                            -6.89934138768139274E-06, &
!!                                             3.91496795038837148E-07, &
!!                                             1.93241838971089117E-09, &
!!                                            -2.79800766779399761E-10, &
!!                                            -1.62135705651220271E-12, &
!!                                             1.04171682753349037E-15, &
!!                                            -3.57704058429328226E-21  ]
!!
!!real(8),dimension(lb:ub),parameter:: c2 = [  3.11570515711405105E-21, &
!!                                            -9.02614207064171034E-16, &
!!                                             1.36753209649678765E-12, &
!!                                             2.38900025108186001E-10, &
!!                                            -1.68356442529926349E-09, &
!!                                            -3.09870229133180094E-07, &
!!                                             5.69141334142564788E-06, &
!!                                              -4.40576577819613139E-05, &
!!                                               1.54894634362815977E-04, &
!!                                               1.84038098379026499E-05, &
!!                                              -8.36977460652474532E-03, &
!!                                               1.77185791599680058E-02, &
!!                                               2.16822571380031852E-02, &
!!                                              -5.82607068849176207E-02, &
!!                                              -9.21129634345036613E-02, &
!!                                               3.51960733427046313E-01, &
!!                                              -3.88907607124389210E-01, &
!!                                               1.92275368853645628E-01, &
!!                                              -3.43692702766296224E-02, &
!!                                              -2.74842897640559984E-03, &
!!                                               1.20206071746381668E-03, &
!!                                              -2.30446494864332518E-04, &
!!                                             2.52663067549429785E-05, &
!!                                            -1.22917063388938916E-06, &
!!                                            -8.53238416007029772E-09, &
!!                                             8.82170026218953071E-10, &
!!                                             5.20220040517289676E-12, &
!!                                            -3.33662182061700503E-15, &
!!                                             1.14571635475684414E-20  ]
!!
!!real(8),dimension(lb:ub),parameter:: e2 = [ -5.57613109685169262E-21, &
!!                                             1.61540333385508532E-15, &
!!                                            -2.44980885342920949E-12, &
!!                                            -4.23976970418981470E-10, &
!!                                             3.62502960829246053E-09, &
!!                                             5.48205480904632642E-07, &
!!                                            -1.09614968447903639E-05, &
!!                                               1.08182315140121479E-04, &
!!                                              -5.33852062456117776E-04, &
!!                                               6.96358245852797650E-04, &
!!                                               2.17704474730355292E-02, &
!!                                              -9.27586382821854322E-02, &
!!                                               1.17357772751682965E-01, &
!!                                              -1.44609613838661572E-02, &
!!                                               2.28332297326772876E-01, &
!!                                              -1.52848053667848732E-01, &
!!                                               4.41317035678250452E-02, &
!!                                               3.62255733762816071E-02, &
!!                                              -1.31322521565223685E-02, &
!!                                              -4.78012836417692096E-04, &
!!                                               3.52336930331104711E-04, &
!!                                              -9.07878558381537606E-05, &
!!                                             1.31453380630463766E-05, &
!!                                            -6.95013224719331636E-07, &
!!                                            -4.05619883260047394E-09, &
!!                                             4.97165393380446650E-10, &
!!                                             2.90404828044886745E-12, &
!!                                            -1.86435317747675210E-15, &
!!                                             6.40177623714996028E-21  ]
!!
!!end module filterModule





!!module filterModule
!!!
!!! Purpose:
!!! ========
!!!   This module contains various filters for the Daubechies 16 wavelet family.
!!!   The meaning is the following:
!!!     a_i=<phi|Op|phi_i>
!!!     b_i=<psi|Op|phi_i>
!!!     c_i=<phi|Op|psi_i>
!!!     e_i=<psi|Op|psi_i>
!!!   There are three different operators:
!!!     Op = second derivative -> a, b, c, e
!!!     Op = x -> a1, b1, c1, e1
!!!     Op = x^2 -> a2, b2, c2, e2
!!!
!!implicit none
!!
!!! The lower and upper bound for the filters.
!!integer,parameter:: lb=-14, ub=14
!!
!!real(8),dimension(lb:ub),parameter:: a = [ -6.92447494063920016E-18, &
!!                                            2.70800493626319433E-13, &
!!                                           -5.81387983028254035E-11, &
!!                                           -1.05857055496741463E-08, &
!!                                           -3.72307630473692755E-07, &
!!                                            2.09042349529203676E-06, &
!!                                           -2.39822852450759980E-05, &
!!                                            4.51679202875022361E-04, &
!!                                           -4.09765689342633799E-03, &
!!                                            2.20702918848225543E-02, &
!!                                           -8.22663999742123397E-02, &
!!                                            2.37178058215380572E-01, &
!!                                           -6.15614146557006969E-01, &
!!                                            2.21914659389111657E+00, &
!!                                           -3.55369228991319019E+00, &
!!                                            2.21914659389111657E+00, &
!!                                           -6.15614146557006969E-01, &
!!                                            2.37178058215380572E-01, &
!!                                           -8.22663999742123397E-02, &
!!                                            2.20702918848225543E-02, &
!!                                           -4.09765689342633799E-03, &
!!                                            4.51679202875022361E-04, &
!!                                           -2.39822852450759980E-05, &
!!                                            2.09042349529203676E-06, &
!!                                           -3.72307630473692755E-07, &
!!                                           -1.05857055496741463E-08, &
!!                                           -5.81387983028254035E-11, &
!!                                            2.70800493626319433E-13, &
!!                                           -6.92447494063920016E-18  ]
!!
!!real(8),dimension(lb:ub),parameter:: b = [  1.23926296291889880E-17, &
!!                                           -4.84665690596158973E-13, &
!!                                            1.04757345962781863E-10, &
!!                                            1.87909235155149931E-08, &
!!                                            6.39191285137934231E-07, &
!!                                           -4.62561497463184309E-06, &
!!                                            5.07580794728972886E-05, &
!!                                           -9.34236233048086760E-04, &
!!                                            8.08720294118447924E-03, &
!!                                           -5.65773697332875719E-02, &
!!                                            2.33704906317512895E-01, &
!!                                           -5.32982239628004284E-01, &
!!                                            6.34570352678925254E-01, &
!!                                           -2.17161175051370531E-01, &
!!                                           -4.20816552937241239E-01, &
!!                                            6.57030740071213981E-01, &
!!                                           -4.24980509437801468E-01, &
!!                                            1.43273293525107509E-01, &
!!                                           -2.84696165863973540E-02, &
!!                                            7.13761218453631829E-03, &
!!                                           -2.13562918387979854E-03, &
!!                                            2.17107954846461359E-04, &
!!                                           -1.13345672451682042E-05, &
!!                                            8.79697040552863148E-07, &
!!                                           -2.16568306292140375E-07, &
!!                                           -5.96264938781402246E-09, &
!!                                           -3.22647023140105305E-11, &
!!                                            1.51306165608661538E-13, &
!!                                           -3.86910241314765619E-18  ]
!!
!!real(8),dimension(lb:ub),parameter:: c = [ -3.86910241314765619E-18, &
!!                                            1.51306165608661538E-13, &
!!                                           -3.22647023140105240E-11, &
!!                                           -5.96264938781402246E-09, &
!!                                           -2.16568306292140402E-07, &
!!                                            8.79697040552863466E-07, &
!!                                           -1.13345672451681957E-05, &
!!                                            2.17107954846461495E-04, &
!!                                           -2.13562918387979854E-03, &
!!                                            7.13761218453631743E-03, &
!!                                           -2.84696165863973470E-02, &
!!                                            1.43273293525107648E-01, &
!!                                           -4.24980509437801190E-01, &
!!                                            6.57030740071214536E-01, &
!!                                           -4.20816552937243626E-01, &
!!                                           -2.17161175051370559E-01, &
!!                                            6.34570352678925476E-01, &
!!                                           -5.32982239628004395E-01, &
!!                                            2.33704906317513006E-01, &
!!                                           -5.65773697332875719E-02, &
!!                                            8.08720294118448098E-03, &
!!                                           -9.34236233048086760E-04, &
!!                                            5.07580794728972886E-05, &
!!                                           -4.62561497463184224E-06, &
!!                                            6.39191285137934125E-07, &
!!                                            1.87909235155149897E-08, &
!!                                            1.04757345962781863E-10, &
!!                                           -4.84665690596158973E-13, &
!!                                            1.23926296291889880E-17  ]
!!
!!real(8),dimension(lb:ub),parameter:: e = [  6.92447494063920016E-18, &
!!                                           -2.70800493626319433E-13, &
!!                                            5.81387983028254035E-11, &
!!                                            1.05857055496741463E-08, &
!!                                            3.72307630473692703E-07, &
!!                                           -2.09042349529203803E-06, &
!!                                            2.39822852450760014E-05, &
!!                                           -4.51679202875077818E-04, &
!!                                            4.09765642831595028E-03, &
!!                                           -2.20732703458663519E-02, &
!!                                            8.20745416922517429E-02, &
!!                                           -2.69959313362790998E-01, &
!!                                           -4.25170532366915405E-02, &
!!                                           -7.14405976634717543E+00, &
!!                                           -2.48758460293923314E+01, &
!!                                           -7.14405976634717366E+00, &
!!                                           -4.25170532366917972E-02, &
!!                                           -2.69959313362791276E-01, &
!!                                            8.20745416922516735E-02, &
!!                                           -2.20732703458663554E-02, &
!!                                            4.09765642831595115E-03, &
!!                                           -4.51679202875077709E-04, &
!!                                            2.39822852450760014E-05, &
!!                                           -2.09042349529203592E-06, &
!!                                            3.72307630473692544E-07, &
!!                                            1.05857055496741480E-08, &
!!                                            5.81387983028254099E-11, &
!!                                           -2.70800493626319433E-13, &
!!                                            6.92447494063920016E-18  ]
!!
!!real(8),dimension(lb:ub),parameter:: a1 = [ -4.27782405260083166E-22, &
!!                                             1.33836789073155300E-16, &
!!                                            -2.23077383438816909E-13, &
!!                                            -4.18701072588700730E-11, &
!!                                             3.84061420317539110E-10, &
!!                                             6.90677058415406752E-08, &
!!                                            -1.50667718130096206E-06, &
!!                                             1.43009564059536508E-05, &
!!                                            -7.38484428527907765E-05, &
!!                                             1.17437400174429801E-04, &
!!                                             4.36133077583543763E-03, &
!!                                            -2.15712170526503010E-02, &
!!                                             2.26678480716213104E-02, &
!!                                             9.18613942126320498E-02, &
!!                                            -3.43336140248461028E-02, &
!!                                             9.18613942126320082E-02, &
!!                                             2.26678480716213208E-02, &
!!                                            -2.15712170526503114E-02, &
!!                                             4.36133077583543936E-03, &
!!                                             1.17437400174429598E-04, &
!!                                            -7.38484428527907901E-05, &
!!                                             1.43009564059536593E-05, &
!!                                            -1.50667718130096206E-06, &
!!                                             6.90677058415406752E-08, &
!!                                             3.84061420317539886E-10, &
!!                                            -4.18701072588700924E-11, &
!!                                            -2.23077383438816909E-13, &
!!                                             1.33836789073155398E-16, &
!!                                            -4.27782405260083542E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: b1 = [  7.65595796896980290E-22, &
!!                                            -2.39526824360595103E-16, &
!!                                             3.99587762949889356E-13, &
!!                                             7.43506270483382512E-11, &
!!                                            -7.94750699421454253E-10, &
!!                                            -1.22260456716396509E-07, &
!!                                             2.87022345454858366E-06, &
!!                                            -3.01746667973454809E-05, &
!!                                             1.87517675653607282E-04, &
!!                                            -5.06158985095035209E-04, &
!!                                            -7.78069990944343145E-03, &
!!                                             5.60208710004729618E-02, &
!!                                            -1.59595684374156721E-01, &
!!                                             2.45149541999515802E-01, &
!!                                            -2.14819128526085146E-01, &
!!                                             9.90323036047134064E-02, &
!!                                            -1.30155486016067679E-02, &
!!                                            -6.91564877941106340E-03, &
!!                                             2.28647193494457184E-03, &
!!                                             3.94195468896341747E-06, &
!!                                            -2.62413556639357196E-05, &
!!                                             6.63758551111563805E-06, &
!!                                            -7.86922170343662117E-07, &
!!                                             3.89648346594732955E-08, &
!!                                             1.80799139957949220E-10, &
!!                                            -2.35773086204665242E-11, &
!!                                            -1.24537049640247525E-13, &
!!                                             7.47819629673571264E-17, &
!!                                            -2.39026633886714879E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: c1 = [ -2.39026633886714738E-22, &
!!                                             7.47819629673570894E-17, &
!!                                            -1.24537049640247500E-13, &
!!                                            -2.35773086204665145E-11, &
!!                                             1.80799139957948548E-10, &
!!                                             3.89648346594732955E-08, &
!!                                            -7.86922170343662328E-07, &
!!                                             6.63758551111563889E-06, &
!!                                            -2.62413556639357196E-05, &
!!                                             3.94195468896342340E-06, &
!!                                             2.28647193494457314E-03, &
!!                                            -6.91564877941106253E-03, &
!!                                            -1.30155486016067592E-02, &
!!                                             9.90323036047133926E-02, &
!!                                            -2.14819128526085035E-01, &
!!                                             2.45149541999515941E-01, &
!!                                            -1.59595684374156666E-01, &
!!                                             5.60208710004729687E-02, &
!!                                            -7.78069990944342624E-03, &
!!                                            -5.06158985095035317E-04, &
!!                                             1.87517675653607309E-04, &
!!                                            -3.01746667973454876E-05, &
!!                                             2.87022345454858366E-06, &
!!                                            -1.22260456716396615E-07, &
!!                                            -7.94750699421456321E-10, &
!!                                             7.43506270483382900E-11, &
!!                                             3.99587762949889406E-13, &
!!                                            -2.39526824360595251E-16, &
!!                                             7.65595796896980854E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: e1 = [  4.27782405260083260E-22, &
!!                                            -1.33836789073155324E-16, &
!!                                             2.23077383438816859E-13, &
!!                                             4.18701072588700601E-11, &
!!                                            -3.84061420317538800E-10, &
!!                                            -6.90677058415406884E-08, &
!!                                             1.50667718130096164E-06, &
!!                                            -1.43009564059536542E-05, &
!!                                             7.38484426297133877E-05, &
!!                                            -1.17437016113009345E-04, &
!!                                            -4.36283745301674006E-03, &
!!                                             2.14973686097975070E-02, &
!!                                            -1.83065172957859110E-02, &
!!                                            -6.91935461410107810E-02, &
!!                                            -3.05473962097790874E-17, &
!!                                            -6.91935461410106561E-02, &
!!                                            -1.83065172957858832E-02, &
!!                                             2.14973686097975382E-02, &
!!                                            -4.36283745301673746E-03, &
!!                                            -1.17437016113009453E-04, &
!!                                             7.38484426297134012E-05, &
!!                                            -1.43009564059536576E-05, &
!!                                             1.50667718130096227E-06, &
!!                                            -6.90677058415407281E-08, &
!!                                            -3.84061420317539989E-10, &
!!                                             4.18701072588700859E-11, &
!!                                             2.23077383438816859E-13, &
!!                                            -1.33836789073155374E-16, &
!!                                             4.27782405260083542E-22  ]
!!
!!real(8),dimension(lb:ub),parameter:: a2 = [  5.57613113797079647E-21, &
!!                                            -1.61540205057621221E-15, &
!!                                             2.44941827237665777E-12, &
!!                                             4.24629970074974494E-10, &
!!                                            -3.50622375577813727E-09, &
!!                                            -5.49912783753825352E-07, &
!!                                             1.07758069134433754E-05, &
!!                                            -8.87797723407857120E-05, &
!!                                             3.82019672542810440E-04, &
!!                                            -4.44838660769164562E-04, &
!!                                            -1.42199024558067322E-02, &
!!                                             3.99360917508602142E-02, &
!!                                             3.36528867483139038E-02, &
!!                                            -2.19314932823181152E-01, &
!!                                             1.65585339069366455E-01, &
!!                                            -3.55921462178230286E-02, &
!!                                             1.24324277043342590E-01, &
!!                                            -8.94912108778953552E-02, &
!!                                             2.06707436591386795E-02, &
!!                                             7.29535357095301151E-04, &
!!                                            -5.04161638673394918E-04, &
!!                                             1.11433619167655706E-04, &
!!                                            -1.33310277306009084E-05, &
!!                                             6.93305878485261928E-07, &
!!                                             4.17500478633314742E-09, &
!!                                            -4.96512386760628033E-10, &
!!                                            -2.90443884221058823E-12, &
!!                                             1.86435445701578929E-15, &
!!                                            -6.40177613495305920E-21  ]
!!
!!real(8),dimension(lb:ub),parameter:: b2 = [ -1.07451145386757544E-20, &
!!                                             3.13060242893980840E-15, &
!!                                            -4.78749354569420448E-12, &
!!                                            -8.27894403417522942E-10, &
!!                                             8.15738060091017371E-09, &
!!                                             1.09377804444408009E-06, &
!!                                            -2.35274919857315769E-05, &
!!                                             2.22483372452623123E-04, &
!!                                            -1.23566735179344547E-03, &
!!                                             2.81931953209789204E-03, &
!!                                             3.56570277578323047E-02, &
!!                                            -1.99870719674157843E-01, &
!!                                             4.09070790464668044E-01, &
!!                                            -3.83487855183457793E-01, &
!!                                             1.22706132268420240E-01, &
!!                                             4.07716118896357879E-02, &
!!                                            -1.73643906901317681E-02, &
!!                                            -1.68596657779329674E-02, &
!!                                             7.63552927760463660E-03, &
!!                                             5.38805216953898019E-05, &
!!                                            -1.33751749897266124E-04, &
!!                                             4.34601255610577065E-05, &
!!                                            -6.11241921733773052E-06, &
!!                                             3.52531960379363515E-07, &
!!                                             1.75161924975294303E-09, &
!!                                            -2.56223458158933282E-10, &
!!                                            -1.49682000687195529E-12, &
!!                                             9.66934864566133037E-16, &
!!                                            -3.33801395040656728E-21  ]
!!
!!real(8),dimension(lb:ub),parameter:: c2 = [  3.11570515711405105E-21, &
!!                                            -9.02614207064171034E-16, &
!!                                             1.36753209649678765E-12, &
!!                                             2.38900025108186001E-10, &
!!                                            -1.68356442529926349E-09, &
!!                                            -3.09870229133180094E-07, &
!!                                             5.69141334142564788E-06, &
!!                                            -4.28284860914085295E-05, &
!!                                             1.54903162339990210E-04, &
!!                                             1.84029297941124655E-05, &
!!                                            -8.36977426829569639E-03, &
!!                                             1.77185781285476296E-02, &
!!                                             2.16822550726257987E-02, &
!!                                            -5.82606916064312488E-02, &
!!                                            -9.21129964344663960E-02, &
!!                                             3.51960771004857065E-01, &
!!                                            -3.88907631541772436E-01, &
!!                                             1.92275377392105401E-01, &
!!                                            -3.43692714446171138E-02, &
!!                                            -2.74842930175147335E-03, &
!!                                             1.20206243182077255E-03, &
!!                                            -2.30136629585558173E-04, &
!!                                             2.52663067549429785E-05, &
!!                                            -1.22917063388938916E-06, &
!!                                            -8.53238416007029772E-09, &
!!                                             8.82170026218953071E-10, &
!!                                             5.20220040517289676E-12, &
!!                                            -3.33662182061700503E-15, &
!!                                             1.14571635475684414E-20  ]
!!
!!real(8),dimension(lb:ub),parameter:: e2 = [ -6.00391350211177579E-21, &
!!                                             1.74924012292824062E-15, &
!!                                            -2.67288623686802693E-12, &
!!                                            -4.65847077677851440E-10, &
!!                                             4.00909102860999587E-09, &
!!                                             6.17273186746173278E-07, &
!!                                            -1.24681740260913262E-05, &
!!                                             1.08808258331504672E-04, &
!!                                            -5.33848378123623628E-04, &
!!                                             6.96357768261355128E-04, &
!!                                             2.17704468186258569E-02, &
!!                                            -9.27586350364214274E-02, &
!!                                             1.17357770131266595E-01, &
!!                                            -1.44609723395050082E-02, &
!!                                             2.28332297556597369E-01, &
!!                                            -1.52848064721562760E-01, &
!!                                             4.41317008735308547E-02, &
!!                                             3.62255766421629780E-02, &
!!                                            -1.31322528112106044E-02, &
!!                                            -4.78012391870289409E-04, &
!!                                             3.52332933393174181E-04, &
!!                                            -9.14051313805750625E-05, &
!!                                             1.16386608817454143E-05, &
!!                                            -6.25945518877790789E-07, &
!!                                            -3.67213741228293157E-09, &
!!                                             4.55295286121576732E-10, &
!!                                             2.68097089701005081E-12, &
!!                                            -1.73051638840359641E-15, &
!!                                             5.97399383188987711E-21  ]
!!
!!end module filterModule






module filterModule
!
! Purpose:
! ========
!   This module contains various filters for the Daubechies 16 wavelet family.
!   The meaning is the following:
!     a_i=<phi|Op|phi_i>
!     b_i=<psi|Op|phi_i>
!     c_i=<phi|Op|psi_i>
!     e_i=<psi|Op|psi_i>
!   There are five different operators:
!     Op = second derivative -> a, b, c, e
!     Op = x -> a1, b1, c1, e1
!     Op = x^2 -> a2, b2, c2, e2
!     Op = x^3 -> a3, b3, c3, e3
!     Op = x^4 -> a4, b4, c4, e4
!
implicit none

! The lower and upper bound for the filters.
integer,parameter:: lb=-14, ub=14

real(8),dimension(lb:ub),parameter:: a = [ -6.92447490505951028E-18, &
                                            2.70800498995346639E-13, &
                                           -5.81387993303650319E-11, &
                                           -1.05857056453828591E-08, &
                                           -3.72307624729728559E-07, &
                                            2.09042354981647804E-06, &
                                           -2.39822857110993937E-05, &
                                            4.51679195975884795E-04, &
                                           -4.09765681251883507E-03, &
                                            2.20702923834323883E-02, &
                                           -8.22663977742195129E-02, &
                                            2.37178057432174683E-01, &
                                           -6.15614175796508789E-01, &
                                            2.21914649009704590E+00, &
                                           -3.55369234085083008E+00, &
                                            2.21914649009704590E+00, &
                                           -6.15614175796508789E-01, &
                                            2.37178057432174683E-01, &
                                           -8.22663977742195129E-02, &
                                            2.20702923834323883E-02, &
                                           -4.09765681251883507E-03, &
                                            4.51679195975884795E-04, &
                                           -2.39822857110993937E-05, &
                                            2.09042354981647804E-06, &
                                           -3.72307624729728559E-07, &
                                           -1.05857056453828591E-08, &
                                           -5.81387993303650319E-11, &
                                            2.70800498995346639E-13, &
                                           -6.92447490505951028E-18  ]

real(8),dimension(lb:ub),parameter:: b = [  1.23926298748915144E-17, &
                                           -4.84665694980683225E-13, &
                                            1.04757348540220623E-10, &
                                            1.87909231522041159E-08, &
                                            6.39191314318015870E-07, &
                                           -4.62561549613067740E-06, &
                                            5.07580797039626318E-05, &
                                           -9.34236181144083530E-04, &
                                            8.08720235045631147E-03, &
                                           -5.65773658045679445E-02, &
                                            2.33704891436108536E-01, &
                                           -5.32982207596253477E-01, &
                                            6.34570315769714122E-01, &
                                           -2.17161162613989767E-01, &
                                           -4.20816528718038985E-01, &
                                            6.57030701554058183E-01, &
                                           -4.24980483446119828E-01, &
                                            1.43273284054659222E-01, &
                                           -2.84696145431842136E-02, &
                                            7.13761173990859910E-03, &
                                           -2.13562905732165133E-03, &
                                            2.17107937709352976E-04, &
                                           -1.13345668089136984E-05, &
                                            8.79697176379788658E-07, &
                                           -2.16568315068759503E-07, &
                                           -5.96264927673747846E-09, &
                                           -3.22647031147176541E-11, &
                                            1.51306166977329782E-13, &
                                           -3.86910248985843272E-18  ]

real(8),dimension(lb:ub),parameter:: c = [ -3.86910248985843272E-18, &
                                            1.51306166977329782E-13, &
                                           -3.22647031147176476E-11, &
                                           -5.96264927673747763E-09, &
                                           -2.16568315068759423E-07, &
                                            8.79697176379788129E-07, &
                                           -1.13345668089136950E-05, &
                                            2.17107937709352895E-04, &
                                           -2.13562905732165177E-03, &
                                            7.13761173990859737E-03, &
                                           -2.84696145431841963E-02, &
                                            1.43273284054659222E-01, &
                                           -4.24980483446119439E-01, &
                                            6.57030701554058183E-01, &
                                           -4.20816528718039373E-01, &
                                           -2.17161162613989239E-01, &
                                            6.34570315769713789E-01, &
                                           -5.32982207596253588E-01, &
                                            2.33704891436108619E-01, &
                                           -5.65773658045679237E-02, &
                                            8.08720235045631668E-03, &
                                           -9.34236181144083638E-04, &
                                            5.07580797039626386E-05, &
                                           -4.62561549613067824E-06, &
                                            6.39191314318015764E-07, &
                                            1.87909231522041125E-08, &
                                            1.04757348540220623E-10, &
                                           -4.84665694980683124E-13, &
                                            1.23926298748915144E-17  ]

real(8),dimension(lb:ub),parameter:: e = [  6.92447507792733360E-18, &
                                           -2.70800496076004972E-13, &
                                            5.81387997394137159E-11, &
                                            1.05857053487793949E-08, &
                                            3.72307646478137301E-07, &
                                           -2.09042376230958751E-06, &
                                            2.39822849010384181E-05, &
                                           -4.51679172658164211E-04, &
                                            4.09765614832326326E-03, &
                                           -2.20732688908272662E-02, &
                                            8.20745364532223937E-02, &
                                           -2.69959298589643515E-01, &
                                           -4.25170863143317437E-02, &
                                           -7.14405969106864269E+00, &
                                           -2.48758457220106521E+01, &
                                           -7.14405969106864003E+00, &
                                           -4.25170863143322295E-02, &
                                           -2.69959298589643404E-01, &
                                            8.20745364532224075E-02, &
                                           -2.20732688908272801E-02, &
                                            4.09765614832326239E-03, &
                                           -4.51679172658164536E-04, &
                                            2.39822849010383978E-05, &
                                           -2.09042376230958878E-06, &
                                            3.72307646478137407E-07, &
                                            1.05857053487793949E-08, &
                                            5.81387997394137095E-11, &
                                           -2.70800496076004972E-13, &
                                            6.92447507792733360E-18  ]

real(8),dimension(lb:ub),parameter:: a1 = [ -4.27782380967095731E-22, &
                                             1.33836790789675412E-16, &
                                            -2.23077388809486687E-13, &
                                            -4.18701080751038290E-11, &
                                             3.84061421554449112E-10, &
                                             6.90677026682351425E-08, &
                                            -1.50667722209618660E-06, &
                                             1.43009565363172442E-05, &
                                            -7.38484450266696513E-05, &
                                             1.17437397420872003E-04, &
                                             4.36133099719882011E-03, &
                                            -2.15712171047925949E-02, &
                                             2.26678475737571716E-02, &
                                             9.18613970279693604E-02, &
                                            -3.43336127698421478E-02, &
                                             9.18613970279693604E-02, &
                                             2.26678475737571716E-02, &
                                            -2.15712171047925949E-02, &
                                             4.36133099719882011E-03, &
                                             1.17437397420872003E-04, &
                                            -7.38484450266696513E-05, &
                                             1.43009565363172442E-05, &
                                            -1.50667722209618660E-06, &
                                             6.90677026682351425E-08, &
                                             3.84061421554449112E-10, &
                                            -4.18701080751038290E-11, &
                                            -2.23077388809486687E-13, &
                                             1.33836790789675412E-16, &
                                            -4.27782380967095731E-22  ]

real(8),dimension(lb:ub),parameter:: b1 = [  7.65595806716108416E-22, &
                                            -2.39526829029937496E-16, &
                                             3.99587744885038038E-13, &
                                             7.43506281339207426E-11, &
                                            -7.94750732307691960E-10, &
                                            -1.22260456570108996E-07, &
                                             2.87022347262966871E-06, &
                                            -3.01746671145041032E-05, &
                                             1.87517678388117519E-04, &
                                            -5.06159000963081235E-04, &
                                            -7.78069985307987440E-03, &
                                             5.60208708828413976E-02, &
                                            -1.59595684240634528E-01, &
                                             2.45149541954812561E-01, &
                                            -2.14819128613303795E-01, &
                                             9.90323037445533516E-02, &
                                            -1.30155486978821185E-02, &
                                            -6.91564874276195691E-03, &
                                             2.28647192633767843E-03, &
                                             3.94195663735356918E-06, &
                                            -2.62413561954799259E-05, &
                                             6.63758559364208097E-06, &
                                            -7.86922176158694601E-07, &
                                             3.89648346451391753E-08, &
                                             1.80799149218002757E-10, &
                                            -2.35773089429423012E-11, &
                                            -1.24537043995966782E-13, &
                                             7.47819644251633525E-17, &
                                            -2.39026636952344337E-22  ]

real(8),dimension(lb:ub),parameter:: c1 = [ -2.39026636952344337E-22, &
                                             7.47819644251633525E-17, &
                                            -1.24537043995966782E-13, &
                                            -2.35773089429423012E-11, &
                                             1.80799149218002732E-10, &
                                             3.89648346451391819E-08, &
                                            -7.86922176158694495E-07, &
                                             6.63758559364208181E-06, &
                                            -2.62413561954799192E-05, &
                                             3.94195663735365473E-06, &
                                             2.28647192633767843E-03, &
                                            -6.91564874276196211E-03, &
                                            -1.30155486978821098E-02, &
                                             9.90323037445533100E-02, &
                                            -2.14819128613303739E-01, &
                                             2.45149541954812505E-01, &
                                            -1.59595684240634500E-01, &
                                             5.60208708828414045E-02, &
                                            -7.78069985307987266E-03, &
                                            -5.06159000963081343E-04, &
                                             1.87517678388117519E-04, &
                                            -3.01746671145040863E-05, &
                                             2.87022347262966871E-06, &
                                            -1.22260456570108996E-07, &
                                            -7.94750732307691960E-10, &
                                             7.43506281339207555E-11, &
                                             3.99587744885038089E-13, &
                                            -2.39526829029937447E-16, &
                                             7.65595806716108416E-22  ]

real(8),dimension(lb:ub),parameter:: e1 = [  4.27782410746594903E-22, &
                                            -1.33836791682177479E-16, &
                                             2.23077373341138652E-13, &
                                             4.18701078507265836E-11, &
                                            -3.84061437780811025E-10, &
                                            -6.90677057885295152E-08, &
                                             1.50667719157836454E-06, &
                                            -1.43009565685797195E-05, &
                                             7.38484438588798688E-05, &
                                            -1.17437022042463091E-04, &
                                            -4.36283743232178042E-03, &
                                             2.14973685507261897E-02, &
                                            -1.83065170245152523E-02, &
                                            -6.91935470160085636E-02, &
                                            -1.11786979288064706E-09, &
                                            -6.91935470160085220E-02, &
                                            -1.83065170245152939E-02, &
                                             2.14973685507262036E-02, &
                                            -4.36283743232178042E-03, &
                                            -1.17437022042463267E-04, &
                                             7.38484438588798823E-05, &
                                            -1.43009565685797144E-05, &
                                             1.50667719157836497E-06, &
                                            -6.90677057885295285E-08, &
                                            -3.84061437780811128E-10, &
                                             4.18701078507265836E-11, &
                                             2.23077373341138652E-13, &
                                            -1.33836791682177479E-16, &
                                             4.27782410746594903E-22  ]

real(8),dimension(lb:ub),parameter:: a2 = [  5.57613113797079647E-21, &
                                            -1.61540205057621221E-15, &
                                             2.44941827237665777E-12, &
                                             4.24629970074974494E-10, &
                                            -3.50622375577813727E-09, &
                                            -5.49912783753825352E-07, &
                                             1.07758069134433754E-05, &
                                            -8.87797723407857120E-05, &
                                             3.82019672542810440E-04, &
                                            -4.44838660769164562E-04, &
                                            -1.42199024558067322E-02, &
                                             3.99360917508602142E-02, &
                                             3.36528867483139038E-02, &
                                            -2.19314932823181152E-01, &
                                             1.65585339069366455E-01, &
                                            -3.55921462178230286E-02, &
                                             1.24324277043342590E-01, &
                                            -8.94912108778953552E-02, &
                                             2.06707436591386795E-02, &
                                             7.29535357095301151E-04, &
                                            -5.04161638673394918E-04, &
                                             1.11433619167655706E-04, &
                                            -1.33310277306009084E-05, &
                                             6.93305878485261928E-07, &
                                             4.17500478633314742E-09, &
                                            -4.96512386760628033E-10, &
                                            -2.90443884221058823E-12, &
                                             1.86435445701578929E-15, &
                                            -6.40177613495305920E-21  ]

real(8),dimension(lb:ub),parameter:: b2 = [ -1.07451146172287584E-20, &
                                             3.13060246628943967E-15, &
                                            -4.78749340108647041E-12, &
                                            -8.27894411800580517E-10, &
                                             8.15738084296199794E-09, &
                                             1.09377804413390616E-06, &
                                            -2.35274921407857352E-05, &
                                             2.22483374802824970E-04, &
                                            -1.23566737058074792E-03, &
                                             2.81931962789967955E-03, &
                                             3.56570274706771578E-02, &
                                            -1.99870719176705625E-01, &
                                             4.09070789978144278E-01, &
                                            -3.83487854930014471E-01, &
                                             1.22706132164835877E-01, &
                                             4.07716120232467089E-02, &
                                            -1.73643908349863270E-02, &
                                            -1.68596657049950964E-02, &
                                             7.63552926096438926E-03, &
                                             5.38805252160584285E-05, &
                                            -1.33751751720543845E-04, &
                                             4.34601260064193970E-05, &
                                            -6.11241925587226456E-06, &
                                             3.52531960342273372E-07, &
                                             1.75161930985329458E-09, &
                                            -2.56223460315090442E-10, &
                                            -1.49681996733928573E-12, &
                                             9.66934874772265743E-16, &
                                            -3.33801397186598064E-21  ]

real(8),dimension(lb:ub),parameter:: c2 = [  3.11570517857345124E-21, &
                                            -9.02614217267228957E-16, &
                                             1.36753205696331878E-12, &
                                             2.38900027287721405E-10, &
                                            -1.68356448485162691E-09, &
                                            -3.09870229248918517E-07, &
                                             5.69141338437313633E-06, &
                                            -4.28284866174613975E-05, &
                                             1.54903165637988938E-04, &
                                             1.84029190592225666E-05, &
                                            -8.36977423915950293E-03, &
                                             1.77185780796068407E-02, &
                                             2.16822550214746182E-02, &
                                            -5.82606911845323996E-02, &
                                            -9.21129972751009579E-02, &
                                             3.51960771815418294E-01, &
                                            -3.88907631878252169E-01, &
                                             1.92275377336910330E-01, &
                                            -3.43692713216321652E-02, &
                                            -2.74842935720907584E-03, &
                                             1.20206244492311446E-03, &
                                            -2.30136631629999558E-04, &
                                             2.52663068944202855E-05, &
                                            -1.22917063308547018E-06, &
                                            -8.53238440365876710E-09, &
                                             8.82170034527123711E-10, &
                                             5.20220026056775328E-12, &
                                            -3.33662185797648247E-15, &
                                             1.14571636261214845E-20  ]

real(8),dimension(lb:ub),parameter:: e2 = [ -6.00391354600385614E-21, &
                                             1.74924014379757953E-15, &
                                            -2.67288615604453199E-12, &
                                            -4.65847082273536332E-10, &
                                             4.00909115868513981E-09, &
                                             6.17273186723220167E-07, &
                                            -1.24681741129386358E-05, &
                                             1.08808259548054658E-04, &
                                            -5.33848387091264277E-04, &
                                             6.96357808776253747E-04, &
                                             2.17704467254226479E-02, &
                                            -9.27586349639902968E-02, &
                                             1.17357770175068044E-01, &
                                            -1.44609721723709447E-02, &
                                             2.28332297558082514E-01, &
                                            -1.52848064755532920E-01, &
                                             4.41317011994162914E-02, &
                                             3.62255764829988031E-02, &
                                            -1.31322527750372502E-02, &
                                            -4.78012401544689199E-04, &
                                             3.52332937508960181E-04, &
                                            -9.14051322471805153E-05, &
                                             1.16386609504637455E-05, &
                                            -6.25945518694966774E-07, &
                                            -3.67213752393608416E-09, &
                                             4.55295290054135794E-10, &
                                             2.68097082627805999E-12, &
                                            -1.73051640666939930E-15, &
                                             5.97399387029546905E-21  ]

real(8),dimension(lb:ub),parameter:: a3 = [ -5.45235389905989279E-20, &
                                             1.46358789074752665E-14, &
                                            -2.01811241329341584E-11, &
                                            -3.24547011487652526E-09, &
                                             2.42857414178843101E-08, &
                                             3.29380782204680145E-06, &
                                            -5.91745010751765221E-05, &
                                             4.36957285273820162E-04, &
                                            -1.67352124117314816E-03, &
                                             1.98640045709908009E-03, &
                                             4.50586304068565369E-02, &
                                            -6.27721697092056274E-02, &
                                            -1.87771692872047424E-01, &
                                             3.47782939672470093E-01, &
                                            -7.22572430968284607E-02, &
                                            -3.45776788890361786E-02, &
                                             2.86159813404083252E-01, &
                                            -2.85770207643508911E-01, &
                                             8.37636739015579224E-02, &
                                             4.12162579596042633E-03, &
                                            -2.77279899455606937E-03, &
                                             6.74822716973721981E-04, &
                                            -8.98371581570245326E-05, &
                                             5.22961454407777637E-06, &
                                             3.43174519912281539E-08, &
                                            -4.43153069795698684E-09, &
                                            -2.83714943899449068E-11, &
                                             1.94904500389536314E-14, &
                                            -7.18620847350200121E-20  ]

real(8),dimension(lb:ub),parameter:: b3 = [  1.13123436610365501E-19, &
                                            -3.07100166700629332E-14, &
                                             4.30376242080603316E-11, &
                                             6.94162274675861823E-09, &
                                            -6.34036165941134751E-08, &
                                            -7.35379995218091993E-06, &
                                             1.47173214489243947E-04, &
                                            -1.28034683673851303E-03, &
                                             6.59497106305552763E-03, &
                                            -1.36921541340832773E-02, &
                                            -1.46897049595819368E-01, &
                                             6.37478980242084381E-01, &
                                            -9.30556215252271834E-01, &
                                             5.59458823389911597E-01, &
                                            -1.72942680845484698E-01, &
                                             1.62694799746086138E-01, &
                                            -1.03342326611766899E-01, &
                                            -2.08251898321085600E-02, &
                                             2.27089101547702371E-02, &
                                             8.34163098143174539E-04, &
                                            -5.63993287410050016E-04, &
                                             2.23449248084776051E-04, &
                                            -3.63125559271239363E-05, &
                                             2.39828253857973301E-06, &
                                             1.28421664625344459E-08, &
                                            -2.09729277533997946E-09, &
                                            -1.34982620213705749E-11, &
                                             9.38388852480617667E-15, &
                                            -3.49671886957245755E-20  ]

real(8),dimension(lb:ub),parameter:: c3 = [ -3.04654359861049750E-20, &
                                             8.17787588062863243E-15, &
                                            -1.12680455877343733E-11, &
                                            -1.82444868673102709E-09, &
                                             1.18723849727281057E-08, &
                                             1.85434546705993780E-06, &
                                            -3.15762398492331706E-05, &
                                             2.17290763713708703E-04, &
                                            -7.38492447647642034E-04, &
                                             3.46249835019798774E-04, &
                                             2.65636962851717381E-02, &
                                            -2.40461111279289159E-02, &
                                            -1.13057521285956350E-01, &
                                             1.75811609409945502E-01, &
                                            -1.49997829824147166E-01, &
                                             4.88522886262723888E-01, &
                                            -8.54944372023607846E-01, &
                                             5.97603435176744080E-01, &
                                            -1.38204695319924187E-01, &
                                            -1.31073095034589335E-02, &
                                             6.26732305736694539E-03, &
                                            -1.36644597735980050E-03, &
                                             1.69343102660099410E-04, &
                                            -9.28314434210990361E-06, &
                                            -6.93099235960034988E-08, &
                                             7.87787729799378268E-09, &
                                             5.08133777611095430E-11, &
                                            -3.48819103637274581E-14, &
                                             1.28610501418155154E-19  ]

real(8),dimension(lb:ub),parameter:: e3 = [  6.32085703717368229E-20, &
                                            -1.71593847930409688E-14, &
                                             2.40296297678319407E-11, &
                                             3.90306169567011407E-09, &
                                            -3.16398262659731909E-08, &
                                            -4.14682388687797282E-06, &
                                             7.86772014015926069E-05, &
                                            -6.42360316804161161E-04, &
                                             3.03629694962606785E-03, &
                                            -4.23534856061406027E-03, &
                                            -9.01621008293588544E-02, &
                                             3.08186486078433797E-01, &
                                            -3.02648997344974635E-01, &
                                            -1.63291833641515064E-03, &
                                             2.40612326154245812E-02, &
                                            -2.52596470969653597E-01, &
                                             1.81819415713290494E-01, &
                                             5.37877230562302591E-02, &
                                            -3.83329370344230869E-02, &
                                            -2.59775806368091333E-03, &
                                             1.40265791479232313E-03, &
                                            -4.59627480916031711E-04, &
                                             6.87230434889209732E-05, &
                                            -4.26390036678194334E-06, &
                                            -2.65855222668025542E-08, &
                                             3.72895715741983598E-09, &
                                             2.41751537799131384E-11, &
                                            -1.67942725094584685E-14, &
                                             6.25802566107744753E-20  ]

real(8),dimension(lb:ub),parameter:: a4 = [  4.73982297509915827E-19, &
                                            -1.17970778003122223E-13, &
                                             1.47867218469599493E-10, &
                                             2.21538325462233843E-08, &
                                            -1.50909258422871062E-07, &
                                            -1.75753048097249120E-05, &
                                             2.94731522444635630E-04, &
                                            -2.00139172375202179E-03, &
                                             7.15284887701272964E-03, &
                                            -9.75563097745180130E-03, &
                                            -1.39017373323440552E-01, &
                                             3.75488102436065674E-02, &
                                             6.29527032375335693E-01, &
                                            -6.22455954551696777E-01, &
                                             2.39198490977287292E-01, &
                                            -1.79768219590187073E-01, &
                                             6.60393893718719482E-01, &
                                            -8.88859748840332031E-01, &
                                             3.33310693502426147E-01, &
                                             2.19652801752090454E-02, &
                                            -1.43004665151238441E-02, &
                                             3.75307234935462475E-03, &
                                            -5.46617549844086170E-04, &
                                             3.51455855707172304E-05, &
                                             2.53031799957170733E-07, &
                                            -3.52819533588899503E-08, &
                                            -2.46440173823359032E-10, &
                                             1.81234964102133800E-13, &
                                            -7.17145273459679160E-19  ]

real(8),dimension(lb:ub),parameter:: b4 = [ -1.05879107559698723E-18, &
                                             2.67975668991362517E-13, &
                                            -3.44037604990392778E-10, &
                                            -5.19398195117999243E-08, &
                                             4.41615034358724110E-07, &
                                             4.40071230383561900E-05, &
                                            -8.30762190462945734E-04, &
                                             6.77564167923146339E-03, &
                                            -3.32405637216197786E-02, &
                                             6.54637771233801125E-02, &
                                             5.95097451648795794E-01, &
                                            -2.00500735905367389E+00, &
                                             2.09342664719743654E+00, &
                                            -8.21319709469622516E-01, &
                                             1.49823989402030955E-01, &
                                             1.42193968646241675E-01, &
                                            -2.58174510470954799E-01, &
                                            -6.34201291534508289E-03, &
                                             6.66580611749569979E-02, &
                                             6.79364923822743660E-03, &
                                            -2.24574587498619586E-03, &
                                             1.06412843446368694E-03, &
                                            -1.95665138116746587E-04, &
                                             1.45486779922684797E-05, &
                                             8.45902475713472292E-08, &
                                            -1.53245377007307477E-08, &
                                            -1.08248520437343017E-10, &
                                             8.10102415572423451E-14, &
                                            -3.25649385965431680E-19  ]

real(8),dimension(lb:ub),parameter:: c4 = [  2.64841177354852064E-19, &
                                            -6.59168226006457639E-14, &
                                             8.25659985103387812E-11, &
                                             1.24446685547390312E-08, &
                                            -7.49608949097183663E-08, &
                                            -9.88743988180539929E-06, &
                                             1.58701212643154980E-04, &
                                            -1.01979783801389810E-03, &
                                             3.34978597853061990E-03, &
                                            -3.11164570844212974E-03, &
                                            -8.21852160064285853E-02, &
                                             1.90046894767117460E-03, &
                                             3.03170077824383655E-01, &
                                            -1.71554362097906837E-01, &
                                            -1.19411738947435606E-01, &
                                             6.67866005617359537E-01, &
                                            -1.84671114749237475E+00, &
                                             1.83677986289889206E+00, &
                                            -5.52785700365707955E-01, &
                                            -6.09059165980459158E-02, &
                                             3.09751752376641980E-02, &
                                            -7.46637518078680385E-03, &
                                             1.02466322093907919E-03, &
                                            -6.24487003129975296E-05, &
                                            -5.05322701667791702E-07, &
                                             6.27526586982343614E-08, &
                                             4.41347806410255288E-10, &
                                            -3.24354727486844677E-13, &
                                             1.28346416737257762E-18  ]

real(8),dimension(lb:ub),parameter:: e4 = [ -5.91607470707858549E-19, &
                                             1.49732874622216787E-13, &
                                            -1.92100738129978936E-10, &
                                            -2.91841511215122400E-08, &
                                             2.23394316572768950E-07, &
                                             2.48001476795504797E-05, &
                                            -4.47594374637332710E-04, &
                                             3.47128223750824399E-03, &
                                            -1.60161771714505095E-02, &
                                             2.38282684039843214E-02, &
                                             3.61599000680735583E-01, &
                                            -9.82482562143620841E-01, &
                                             7.52753810635772846E-01, &
                                            -1.14204696883252702E-01, &
                                             2.51501726273655912E-01, &
                                            -4.84276384425144313E-01, &
                                             5.62339757401709828E-01, &
                                             2.85047897727699270E-02, &
                                            -1.07918110020030947E-01, &
                                            -1.51435426900326242E-02, &
                                             5.34875354394725826E-03, &
                                            -2.14609073676606246E-03, &
                                             3.67972099703273354E-04, &
                                            -2.58921736961675790E-05, &
                                            -1.72989790542752900E-07, &
                                             2.72620366731215101E-08, &
                                             1.93858659954877238E-10, &
                                            -1.44983357227351277E-13, &
                                             5.82809854008927447E-19  ]

end module filterModule





subroutine getEffectiveFilter(it,parabPrefac,hgrid, x0, eff, filterCode)
!
! Purpose:
! ========
!   Calculates the effective filter for the operator (x-x0)^2. Using the notation
!   in the header of the module 'filterModule', this effective filter is given by
!   aeff_i = a2_i -2*x0*a1_i + x0^2*delta(i) (delta being the Kronecker delta).
!   The other filters (beff, ceff, eeff) are given in the analogous way.
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
real(8):: fac, fac2, prefac1, prefac2


prefac1=-.5d0/hgrid**2
fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
!fac=parabPrefac
!fac2=parabPrefac*hgrid

select case(trim(filterCode))
case('a')
    do i=lb,ub
        eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
case('b')
    do i=lb,ub
        eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
    end do
case('c')
    do i=lb,ub
        eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
    end do
case('e')
    do i=lb,ub
        eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
case default
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end select




end subroutine getEffectiveFilter









  !   y = (kinetic energy operator)x + (cprec*I)x + ((x-x0)^2*I)*x
! One of the most CPU intensive routines
subroutine ConvolkineticParabola(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzParabola, parabPrefac, it)
  use module_base
  implicit none
!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

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
real(8),dimension(3):: rxyzParabola
real(8):: parabPrefac
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

!write(901,*) x_c
!write(902,*) x_f
!write(903,*) x_f1
!write(904,*) x_f2
!write(905,*) x_f3
!call mpi_barrier(mpi_comm_world, ii)
!stop

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
              x0=hgrid*(i1+0)-rxyzParabola(1)
              x1=hgrid*(i1+1)-rxyzParabola(1)
              x2=hgrid*(i1+2)-rxyzParabola(1)
              x3=hgrid*(i1+3)-rxyzParabola(1)
              call getEffectiveFilter(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, x1, aeff1(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, x2, aeff2(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, x3, aeff3(lowfil), 'a')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilter(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a')
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
              x0=hgrid*(i1+0)-rxyzParabola(1)
              x1=hgrid*(i1+1)-rxyzParabola(1)
              x2=hgrid*(i1+2)-rxyzParabola(1)
              x3=hgrid*(i1+3)-rxyzParabola(1)
              call getEffectiveFilter(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, x1, beff1(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, x2, beff2(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, x3, beff3(lowfil), 'b')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilter(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b')
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
              x0=hgrid*(i1+0)-rxyzParabola(1)
              x1=hgrid*(i1+1)-rxyzParabola(1)
              x2=hgrid*(i1+2)-rxyzParabola(1)
              x3=hgrid*(i1+3)-rxyzParabola(1)
              call getEffectiveFilter(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, x1, ceff1(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, x2, ceff2(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, x3, ceff3(lowfil), 'c')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilter(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c')
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
              y0=hgrid*(i2+0)-rxyzParabola(2)
              y1=hgrid*(i2+1)-rxyzParabola(2)
              y2=hgrid*(i2+2)-rxyzParabola(2)
              y3=hgrid*(i2+3)-rxyzParabola(2)
              call getEffectiveFilter(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, y1, aeff1(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, y2, aeff2(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, y3, aeff3(lowfil), 'a')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilter(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a')
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
              y0=hgrid*(i2+0)-rxyzParabola(2)
              y1=hgrid*(i2+1)-rxyzParabola(2)
              y2=hgrid*(i2+2)-rxyzParabola(2)
              y3=hgrid*(i2+3)-rxyzParabola(2)
              call getEffectiveFilter(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, y1, beff1(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, y2, beff2(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, y3, beff3(lowfil), 'b')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilter(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b')
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
              y0=hgrid*(i2+0)-rxyzParabola(2)
              y1=hgrid*(i2+1)-rxyzParabola(2)
              y2=hgrid*(i2+2)-rxyzParabola(2)
              y3=hgrid*(i2+3)-rxyzParabola(2)
              call getEffectiveFilter(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, y1, ceff1(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, y2, ceff2(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, y3, ceff3(lowfil), 'c')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilter(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c')
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
              z0=hgrid*(i3+0)-rxyzParabola(3)
              z1=hgrid*(i3+1)-rxyzParabola(3)
              z2=hgrid*(i3+2)-rxyzParabola(3)
              z3=hgrid*(i3+3)-rxyzParabola(3)
              call getEffectiveFilter(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, z1, aeff1(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, z2, aeff2(lowfil), 'a')
              call getEffectiveFilter(it,parabPrefac,hgrid, z3, aeff3(lowfil), 'a')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           call getEffectiveFilter(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a')
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
              z0=hgrid*(i3+0)-rxyzParabola(3)
              z1=hgrid*(i3+1)-rxyzParabola(3)
              z2=hgrid*(i3+2)-rxyzParabola(3)
              z3=hgrid*(i3+3)-rxyzParabola(3)
              call getEffectiveFilter(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, z1, beff1(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, z2, beff2(lowfil), 'b')
              call getEffectiveFilter(it,parabPrefac,hgrid, z3, beff3(lowfil), 'b')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           call getEffectiveFilter(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b')
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
              z0=hgrid*(i3+0)-rxyzParabola(3)
              z1=hgrid*(i3+1)-rxyzParabola(3)
              z2=hgrid*(i3+2)-rxyzParabola(3)
              z3=hgrid*(i3+3)-rxyzParabola(3)
              call getEffectiveFilter(it,parabPrefac,hgrid, z0, ceff0(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, z1, ceff1(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, z2, ceff2(lowfil), 'c')
              call getEffectiveFilter(it,parabPrefac,hgrid, z3, ceff3(lowfil), 'c')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           call getEffectiveFilter(it,parabPrefac,hgrid,  z0, ceff0(lowfil), 'c')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilter(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a')
           call getEffectiveFilter(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b')
           call getEffectiveFilter(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c')
           call getEffectiveFilter(it,parabPrefac,hgrid, x0, eeff0(lowfil), 'e')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilter(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a')
           call getEffectiveFilter(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b')
           call getEffectiveFilter(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c')
           call getEffectiveFilter(it,parabPrefac,hgrid, y0, eeff0(lowfil), 'e')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           !call getEffectiveFilter(it,parabPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           call getEffectiveFilter(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a')
           call getEffectiveFilter(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b')
           call getEffectiveFilter(it,parabPrefac,hgrid, z0, ceff0(lowfil), 'c')
           call getEffectiveFilter(it,parabPrefac,hgrid, z0, eeff0(lowfil), 'e')
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




END SUBROUTINE ConvolkineticParabola








subroutine getEffectiveFilter2(it,parabPrefac,hgrid, x0, eff, filterCode, shift)
!
! Purpose:
! ========
!   Calculates the effective filter for the operator (x-x0)^2. Using the notation
!   in the header of the module 'filterModule', this effective filter is given by
!   aeff_i = a2_i -2*x0*a1_i + x0^2*delta(i) (delta being the Kronecker delta).
!   The other filters (beff, ceff, eeff) are given in the analogous way.
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
real(8),intent(in):: parabPrefac, hgrid, x0, shift
real(8),dimension(lb:ub),intent(out):: eff
character(len=*):: filterCode

! Local variables
integer:: i
real(8):: fac, fac2, prefac1, prefac2

!shift=0.d-1

!fac=2.d0*parabPrefac*shift
fac=2.d0*parabPrefac*shift*hgrid
!write(*,*) 'fac', fac

select case(trim(filterCode))
case('a')
    do i=lb,ub
        eff(i)=fac*a1(i)
    end do
    eff(0)=eff(0)+fac/hgrid*x0
case('b')
    do i=lb,ub
        eff(i)=fac*b1(i)
    end do
case('c')
    do i=lb,ub
        eff(i)=fac*c1(i)
    end do
case('e')
    do i=lb,ub
        eff(i)=fac*e1(i)
    end do
    eff(0)=eff(0)+fac/hgrid*x0
case default
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end select




end subroutine getEffectiveFilter2




  !   y = (kinetic energy operator)x + (cprec*I)x + ((x-x0)^2*I)*x
! One of the most CPU intensive routines
subroutine ConvolLinear(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzParabola, parabPrefac, it, direction, shift)
  use module_base
  implicit none
!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

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
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
real(8),dimension(3):: rxyzParabola
real(8):: parabPrefac
integer:: it
character(len=1):: direction
real(8):: shift
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

!write(901,*) x_c
!write(902,*) x_f
!write(903,*) x_f1
!write(904,*) x_f2
!write(905,*) x_f3
!call mpi_barrier(mpi_comm_world, ii)
!stop


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

! this is maybe not needed
y_c=0.d0
y_f=0.d0

if(direction=='x') then
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
                  x0=hgrid*(i1+0)-rxyzParabola(1)
                  x1=hgrid*(i1+1)-rxyzParabola(1)
                  x2=hgrid*(i1+2)-rxyzParabola(1)
                  x3=hgrid*(i1+3)-rxyzParabola(1)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x1, aeff1(lowfil), 'a', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x2, aeff2(lowfil), 'a', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x3, aeff3(lowfil), 'a', shift)
                  do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                     dyi0=dyi0 + x_c(t,i2,i3)*aeff0(t-i1-0)
                     dyi1=dyi1 + x_c(t,i2,i3)*aeff1(t-i1-1)
                     dyi2=dyi2 + x_c(t,i2,i3)*aeff2(t-i1-2)
                     dyi3=dyi3 + x_c(t,i2,i3)*aeff3(t-i1-3)
                  enddo
                  y_c(i1+0,i2,i3)=dyi0!+cprecr*x_c(i1+0,i2,i3)
                  y_c(i1+1,i2,i3)=dyi1!+cprecr*x_c(i1+1,i2,i3)
                  y_c(i1+2,i2,i3)=dyi2!+cprecr*x_c(i1+2,i2,i3)
                  y_c(i1+3,i2,i3)=dyi3!+cprecr*x_c(i1+3,i2,i3)
               enddo
               icur=i1
            else
               icur=ibyz_c(1,i2,i3)
            endif
    
            do i1=icur,ibyz_c(2,i2,i3)
               dyi=0.0_wp 
               ! Get the effective a-filters for the x dimension
               x0=hgrid*(i1+0)-rxyzParabola(1)
               call getEffectiveFilter2(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a', shift)
               do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                  dyi=dyi + x_c(t,i2,i3)*aeff0(t-i1)
               enddo
               y_c(i1,i2,i3)=dyi!+cprecr*x_c(i1,i2,i3)
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
                  x0=hgrid*(i1+0)-rxyzParabola(1)
                  x1=hgrid*(i1+1)-rxyzParabola(1)
                  x2=hgrid*(i1+2)-rxyzParabola(1)
                  x3=hgrid*(i1+3)-rxyzParabola(1)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x1, beff1(lowfil), 'b', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x2, beff2(lowfil), 'b', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x3, beff3(lowfil), 'b', shift)
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
               x0=hgrid*(i1+0)-rxyzParabola(1)
               call getEffectiveFilter2(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b', shift)
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
                  x0=hgrid*(i1+0)-rxyzParabola(1)
                  x1=hgrid*(i1+1)-rxyzParabola(1)
                  x2=hgrid*(i1+2)-rxyzParabola(1)
                  x3=hgrid*(i1+3)-rxyzParabola(1)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x1, ceff1(lowfil), 'c', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x2, ceff2(lowfil), 'c', shift)
                  call getEffectiveFilter2(it,parabPrefac,hgrid, x3, ceff3(lowfil), 'c', shift)
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
               x0=hgrid*(i1+0)-rxyzParabola(1)
               call getEffectiveFilter2(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c', shift)
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

else if(direction=='y') then
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
                y0=hgrid*(i2+0)-rxyzParabola(2)
                y1=hgrid*(i2+1)-rxyzParabola(2)
                y2=hgrid*(i2+2)-rxyzParabola(2)
                y3=hgrid*(i2+3)-rxyzParabola(2)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y1, aeff1(lowfil), 'a', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y2, aeff2(lowfil), 'a', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y3, aeff3(lowfil), 'a', shift)
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                   dyi0=dyi0 + x_c(i1,t,i3)*aeff0(t-i2-0)
                   dyi1=dyi1 + x_c(i1,t,i3)*aeff1(t-i2-1)
                   dyi2=dyi2 + x_c(i1,t,i3)*aeff2(t-i2-2)
                   dyi3=dyi3 + x_c(i1,t,i3)*aeff3(t-i2-3)
                enddo
                !y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
                !y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
                !y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
                !y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
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
             ! Get the effective a-filters for the y dimension
             y0=hgrid*(i2+0)-rxyzParabola(2)
             call getEffectiveFilter2(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a', shift)
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
                y0=hgrid*(i2+0)-rxyzParabola(2)
                y1=hgrid*(i2+1)-rxyzParabola(2)
                y2=hgrid*(i2+2)-rxyzParabola(2)
                y3=hgrid*(i2+3)-rxyzParabola(2)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y1, beff1(lowfil), 'b', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y2, beff2(lowfil), 'b', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y3, beff3(lowfil), 'b', shift)
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
             y0=hgrid*(i2+0)-rxyzParabola(2)
             call getEffectiveFilter2(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b', shift)
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
                y0=hgrid*(i2+0)-rxyzParabola(2)
                y1=hgrid*(i2+1)-rxyzParabola(2)
                y2=hgrid*(i2+2)-rxyzParabola(2)
                y3=hgrid*(i2+3)-rxyzParabola(2)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y1, ceff1(lowfil), 'c', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y2, ceff2(lowfil), 'c', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, y3, ceff3(lowfil), 'c', shift)
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
             y0=hgrid*(i2+0)-rxyzParabola(2)
             call getEffectiveFilter2(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c', shift)
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

else if(direction=='z') then
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
                z0=hgrid*(i3+0)-rxyzParabola(3)
                z1=hgrid*(i3+1)-rxyzParabola(3)
                z2=hgrid*(i3+2)-rxyzParabola(3)
                z3=hgrid*(i3+3)-rxyzParabola(3)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z1, aeff1(lowfil), 'a', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z2, aeff2(lowfil), 'a', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z3, aeff3(lowfil), 'a', shift)
                do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                   dyi0=dyi0 + x_c(i1,i2,t)*aeff0(t-i3-0)
                   dyi1=dyi1 + x_c(i1,i2,t)*aeff1(t-i3-1)
                   dyi2=dyi2 + x_c(i1,i2,t)*aeff2(t-i3-2)
                   dyi3=dyi3 + x_c(i1,i2,t)*aeff3(t-i3-3)
                enddo
                !y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
                !y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
                !y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
                !y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
                y_c(i1,i2,i3+0)=dyi0
                y_c(i1,i2,i3+1)=dyi1
                y_c(i1,i2,i3+2)=dyi2
                y_c(i1,i2,i3+3)=dyi3
             enddo
             icur=i3
          else
             icur=ibxy_c(1,i1,i2)
          endif
  
          do i3=icur,ibxy_c(2,i1,i2)
             dyi=0.0_wp
             ! Get the effective a-filters for the y dimension
             z0=hgrid*(i3+0)-rxyzParabola(3)
             call getEffectiveFilter2(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a', shift)
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
                z0=hgrid*(i3+0)-rxyzParabola(3)
                z1=hgrid*(i3+1)-rxyzParabola(3)
                z2=hgrid*(i3+2)-rxyzParabola(3)
                z3=hgrid*(i3+3)-rxyzParabola(3)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z1, beff1(lowfil), 'b', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z2, beff2(lowfil), 'b', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z3, beff3(lowfil), 'b', shift)
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
             z0=hgrid*(i3+0)-rxyzParabola(3)
             call getEffectiveFilter2(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b', shift)
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
                z0=hgrid*(i3+0)-rxyzParabola(3)
                z1=hgrid*(i3+1)-rxyzParabola(3)
                z2=hgrid*(i3+2)-rxyzParabola(3)
                z3=hgrid*(i3+3)-rxyzParabola(3)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z0, ceff0(lowfil), 'c', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z1, ceff1(lowfil), 'c', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z2, ceff2(lowfil), 'c', shift)
                call getEffectiveFilter2(it,parabPrefac,hgrid, z3, ceff3(lowfil), 'c', shift)
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
             z0=hgrid*(i3+0)-rxyzParabola(3)
             call getEffectiveFilter2(it,parabPrefac,hgrid,  z0, ceff0(lowfil), 'c', shift)
             do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                dyi=dyi + x_c(i1,i2,t)*ceff0(t-i3)
             enddo
             y_f(4,i1,i2,i3)=dyi
          enddo
       enddo
    enddo
    !$omp enddo
else
    write(*,*) "ERROR: allowed values for 'direction' are x, y or z, whereas we found ", trim(direction)
    stop
end if

  
  !  call system_clock(ncount3,ncount_rate,ncount_max)
  !  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2

if(direction=='x') then
    !$omp do
    do i3=nfl3,nfu3
       do i2=nfl2,nfu2
          do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
             t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
             ! Get the effective filters for the x dimension
             x0=hgrid*(i1+0)-rxyzParabola(1)
             call getEffectiveFilter2(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a', shift)
             call getEffectiveFilter2(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b', shift)
             call getEffectiveFilter2(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c', shift)
             call getEffectiveFilter2(it,parabPrefac,hgrid, x0, eeff0(lowfil), 'e', shift)
             do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                t112=t112 + x_f(4,i1+l,i2,i3)*aeff0(l) + x_f(5,i1+l,i2,i3)*beff0(l)
                t121=t121 + x_f(2,i1+l,i2,i3)*aeff0(l) + x_f(3,i1+l,i2,i3)*beff0(l)
                t122=t122 + x_f(6,i1+l,i2,i3)*aeff0(l) + x_f(7,i1+l,i2,i3)*beff0(l)
                t212=t212 + x_f(4,i1+l,i2,i3)*ceff0(l) + x_f(5,i1+l,i2,i3)*eeff0(l)
                t221=t221 + x_f(2,i1+l,i2,i3)*ceff0(l) + x_f(3,i1+l,i2,i3)*eeff0(l)
                t222=t222 + x_f(6,i1+l,i2,i3)*ceff0(l) + x_f(7,i1+l,i2,i3)*eeff0(l)
                t211=t211 + x_f(1,i1+l,i2,i3)*eeff0(l)
             enddo
             y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112!+cprecr*x_f(4,i1,i2,i3)
             y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121!+cprecr*x_f(2,i1,i2,i3)
             y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211!+cprecr*x_f(1,i1,i2,i3)
             y_f(6,i1,i2,i3)=t122!+cprecr*x_f(6,i1,i2,i3)
             y_f(5,i1,i2,i3)=t212!+cprecr*x_f(5,i1,i2,i3)
             y_f(3,i1,i2,i3)=t221!+cprecr*x_f(3,i1,i2,i3)
             y_f(7,i1,i2,i3)=t222!+cprecr*x_f(7,i1,i2,i3)
          enddo
       enddo
    enddo
    !$omp enddo

  !  call system_clock(ncount4,ncount_rate,ncount_max)
  !  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !  write(99,'(a40,1x,e10.3,1x,f6.1)') 'SECND PART:x',tel,1.d-6*nflop1/tel


else if(direction=='y') then
      ! + (1/2) d^2/dy^2
      !$omp do
      do i3=nfl3,nfu3
         do i1=nfl1,nfu1
            do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
               t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
               ! Get the effective filters for the y dimension
               y0=hgrid*(i2+0)-rxyzParabola(2)
               call getEffectiveFilter2(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a', shift)
               call getEffectiveFilter2(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b', shift)
               call getEffectiveFilter2(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c', shift)
               call getEffectiveFilter2(it,parabPrefac,hgrid, y0, eeff0(lowfil), 'e', shift)
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

else if(direction=='z') then
    ! + (1/2) d^2/dz^2
    !$omp do
    do i2=nfl2,nfu2
       do i1=nfl1,nfu1
          do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
             t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
             ! Get the effective filters for the z dimension
             z0=hgrid*(i3+0)-rxyzParabola(3)
             !call getEffectiveFilter2(it,parabPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
             call getEffectiveFilter2(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a', shift)
             call getEffectiveFilter2(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b', shift)
             call getEffectiveFilter2(it,parabPrefac,hgrid, z0, ceff0(lowfil), 'c', shift)
             call getEffectiveFilter2(it,parabPrefac,hgrid, z0, eeff0(lowfil), 'e', shift)
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
end if

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




END SUBROUTINE ConvolLinear




subroutine getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, eff, filterCode)
!
! Purpose:
! ========
!   Calculates the effective filter for the operator (x-x0)^2. Using the notation
!   in the header of the module 'filterModule', this effective filter is given by
!   aeff_i = a2_i -2*x0*a1_i + x0^2*delta(i) (delta being the Kronecker delta).
!   The other filters (beff, ceff, eeff) are given in the analogous way.
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









  !   y = (kinetic energy operator)x + (cprec*I)x + ((x-x0)^2*I)*x
! One of the most CPU intensive routines
subroutine ConvolkineticQuartic(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f,x_f1,x_f2,x_f3, &
     rxyzParabola, parabPrefac, it)
  use module_base
  implicit none
!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

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
real(8),dimension(3):: rxyzParabola
real(8):: parabPrefac
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

!write(901,*) x_c
!write(902,*) x_f
!write(903,*) x_f1
!write(904,*) x_f2
!write(905,*) x_f3
!call mpi_barrier(mpi_comm_world, ii)
!stop

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
              x0=hgrid*(i1+0)-rxyzParabola(1)
              x1=hgrid*(i1+1)-rxyzParabola(1)
              x2=hgrid*(i1+2)-rxyzParabola(1)
              x3=hgrid*(i1+3)-rxyzParabola(1)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x1, aeff1(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x2, aeff2(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x3, aeff3(lowfil), 'a')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a')
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
              x0=hgrid*(i1+0)-rxyzParabola(1)
              x1=hgrid*(i1+1)-rxyzParabola(1)
              x2=hgrid*(i1+2)-rxyzParabola(1)
              x3=hgrid*(i1+3)-rxyzParabola(1)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x1, beff1(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x2, beff2(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x3, beff3(lowfil), 'b')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b')
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
              x0=hgrid*(i1+0)-rxyzParabola(1)
              x1=hgrid*(i1+1)-rxyzParabola(1)
              x2=hgrid*(i1+2)-rxyzParabola(1)
              x3=hgrid*(i1+3)-rxyzParabola(1)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x1, ceff1(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x2, ceff2(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x3, ceff3(lowfil), 'c')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c')
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
              y0=hgrid*(i2+0)-rxyzParabola(2)
              y1=hgrid*(i2+1)-rxyzParabola(2)
              y2=hgrid*(i2+2)-rxyzParabola(2)
              y3=hgrid*(i2+3)-rxyzParabola(2)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y1, aeff1(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y2, aeff2(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y3, aeff3(lowfil), 'a')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a')
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
              y0=hgrid*(i2+0)-rxyzParabola(2)
              y1=hgrid*(i2+1)-rxyzParabola(2)
              y2=hgrid*(i2+2)-rxyzParabola(2)
              y3=hgrid*(i2+3)-rxyzParabola(2)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y1, beff1(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y2, beff2(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y3, beff3(lowfil), 'b')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b')
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
              y0=hgrid*(i2+0)-rxyzParabola(2)
              y1=hgrid*(i2+1)-rxyzParabola(2)
              y2=hgrid*(i2+2)-rxyzParabola(2)
              y3=hgrid*(i2+3)-rxyzParabola(2)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y1, ceff1(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y2, ceff2(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y3, ceff3(lowfil), 'c')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c')
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
              z0=hgrid*(i3+0)-rxyzParabola(3)
              z1=hgrid*(i3+1)-rxyzParabola(3)
              z2=hgrid*(i3+2)-rxyzParabola(3)
              z3=hgrid*(i3+3)-rxyzParabola(3)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z1, aeff1(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z2, aeff2(lowfil), 'a')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z3, aeff3(lowfil), 'a')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a')
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
              z0=hgrid*(i3+0)-rxyzParabola(3)
              z1=hgrid*(i3+1)-rxyzParabola(3)
              z2=hgrid*(i3+2)-rxyzParabola(3)
              z3=hgrid*(i3+3)-rxyzParabola(3)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z1, beff1(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z2, beff2(lowfil), 'b')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z3, beff3(lowfil), 'b')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b')
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
              z0=hgrid*(i3+0)-rxyzParabola(3)
              z1=hgrid*(i3+1)-rxyzParabola(3)
              z2=hgrid*(i3+2)-rxyzParabola(3)
              z3=hgrid*(i3+3)-rxyzParabola(3)
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, ceff0(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z1, ceff1(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z2, ceff2(lowfil), 'c')
              call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z3, ceff3(lowfil), 'c')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid,  z0, ceff0(lowfil), 'c')
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
           x0=hgrid*(i1+0)-rxyzParabola(1)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, aeff0(lowfil), 'a')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, beff0(lowfil), 'b')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, ceff0(lowfil), 'c')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, x0, eeff0(lowfil), 'e')
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
           y0=hgrid*(i2+0)-rxyzParabola(2)
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, aeff0(lowfil), 'a')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, beff0(lowfil), 'b')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, ceff0(lowfil), 'c')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, y0, eeff0(lowfil), 'e')
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
           z0=hgrid*(i3+0)-rxyzParabola(3)
           !call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, aeff0(lowfil), beff0(lowfil), ceff0(lowfil), eeff0(lowfil))
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, aeff0(lowfil), 'a')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, beff0(lowfil), 'b')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, ceff0(lowfil), 'c')
           call getEffectiveFilterQuartic(it,parabPrefac,hgrid, z0, eeff0(lowfil), 'e')
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
