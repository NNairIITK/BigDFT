!> @file
!!  Old convolution for kinetic operators
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine Convolkinetic(n1,n2,n3,&
               nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,&
               nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
               cprecr,hgrid,logrid_c,logrid_f,x,y)
!   y = (kinetic energy operator)x + (cprec*I)x 
    implicit real(kind=8) (a-h,o-z)
    dimension x(nl1_c:nu1_c,2,nl2_c:nu2_c,2,nl3_c:nu3_c,2)
    dimension y(nl1_c:nu1_c,2,nl2_c:nu2_c,2,nl3_c:nu3_c,2)
    logical logrid_c(0:n1,0:n2,0:n3),logrid_f(0:n1,0:n2,0:n3)

    parameter(lowfil=-14,lupfil=14)
    dimension a(lowfil:lupfil),b(lowfil:lupfil),c(lowfil:lupfil),e(lowfil:lupfil)
    scale=-.5d0/hgrid**2
!---------------------------------------------------------------------------
! second derivative filters for Daubechies 16
!  <phi|D^2|phi_i>
    a(0)=   -3.5536922899131901941296809374d0*scale
    a(1)=    2.2191465938911163898794546405d0*scale
    a(2)=   -0.6156141465570069496314853949d0*scale
    a(3)=    0.2371780582153805636239247476d0*scale
    a(4)=   -0.0822663999742123340987663521d0*scale
    a(5)=    0.02207029188482255523789911295638968409d0*scale
    a(6)=   -0.409765689342633823899327051188315485d-2*scale
    a(7)=    0.45167920287502235349480037639758496d-3*scale
    a(8)=   -0.2398228524507599670405555359023135d-4*scale
    a(9)=    2.0904234952920365957922889447361d-6*scale
    a(10)=  -3.7230763047369275848791496973044d-7*scale
    a(11)=  -1.05857055496741470373494132287d-8*scale
    a(12)=  -5.813879830282540547959250667d-11*scale
    a(13)=   2.70800493626319438269856689037647576d-13*scale
    a(14)=  -6.924474940639200152025730585882d-18*scale
    do i=1,14
        a(-i)=a(i)
    enddo
!  <phi|D^2|psi_i>
    c(-14)=     -3.869102413147656535541850057188d-18*scale
    c(-13)=      1.5130616560866154733900029272077362d-13*scale
    c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
    c(-11)=     -5.96264938781402337319841002642d-9*scale
    c(-10)=     -2.1656830629214041470164889350342d-7*scale
    c(-9 )=      8.7969704055286288323596890609625d-7*scale
    c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
    c(-7 )=      0.00021710795484646138591610188464622454d0*scale
    c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
    c(-5 )=      0.00713761218453631422925717625758502986d0*scale
    c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
    c(-3 )=      0.14327329352510759457155821037742893841d0*scale
    c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
    c(-1 )=      0.65703074007121357894896358254040272157d0*scale
    c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
    c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
    c( 2 )=      0.63457035267892488185929915286969303251d0*scale
    c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
    c( 4 )=      0.23370490631751294307619384973520033236d0*scale
    c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
    c( 6 )=      0.0080872029411844780634067667008050127d0*scale
    c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
    c( 8 )=      0.00005075807947289728306309081261461095d0*scale
    c( 9 )=     -4.62561497463184262755416490048242d-6*scale
    c( 10)=      6.3919128513793415587294752371778d-7*scale
    c( 11)=      1.87909235155149902916133888931d-8*scale
    c( 12)=      1.04757345962781829480207861447155543883d-10*scale
    c( 13)=     -4.84665690596158959648731537084025836d-13*scale
    c( 14)=      1.2392629629188986192855777620877d-17*scale
!  <psi|D^2|phi_i>
    do i=-14,14
        b(i)=c(-i)
    enddo
    !<psi|D^2|psi_i>
    e(0)=   -24.875846029392331358907766562d0*scale
    e(1)=   -7.1440597663471719869313377994d0*scale
    e(2)=   -0.04251705323669172315864542163525830944d0*scale
    e(3)=   -0.26995931336279126953587091167128839196d0*scale
    e(4)=    0.08207454169225172612513390763444496516d0*scale
    e(5)=   -0.02207327034586634477996701627614752761d0*scale
    e(6)=    0.00409765642831595181639002667514310145d0*scale
    e(7)=   -0.00045167920287507774929432548999880117d0*scale
    e(8)=    0.00002398228524507599670405555359023135d0*scale
    e(9)=   -2.0904234952920365957922889447361d-6*scale
    e(10)=   3.7230763047369275848791496973044d-7*scale
    e(11)=   1.05857055496741470373494132287d-8*scale
    e(12)=   5.8138798302825405479592506674648873655d-11*scale
    e(13)=  -2.70800493626319438269856689037647576d-13*scale
    e(14)=   6.924474940639200152025730585882d-18*scale
    do i=1,14
        e(-i)=e(i)
    enddo
!---------------------------------------------------------------------------
                                                            
! Scaling function part

       call cpu_time(tcpu0) ; call system_clock(ncount0,ncount_rate,ncount_max)

    do i3=nl3_c,nu3_c
        do i2=nl2_c,nu2_c
            do i1=nl1_c,nu1_c
            if (logrid_c(i1,i2,i3)) then
! (1/2) d^2/dx^2
                t111=0.d0
                s111=0.d0
                do l=max(nl1_c-i1,lowfil),min(lupfil,nu1_c-i1)
                    t111=t111 + x(i1+l,1,i2,1,i3,1)*a(l)
                    s111=s111 + x(i1+l,2,i2,1,i3,1)*b(l)
                enddo
! + (1/2) d^2/dy^2
                do l=max(nl2_c-i2,lowfil),min(lupfil,nu2_c-i2)
                    t111=t111 + x(i1,1,i2+l,1,i3,1)*a(l)
                    s111=s111 + x(i1,1,i2+l,2,i3,1)*b(l)
                enddo
! + (1/2) d^2/dz^2
                do l=max(nl3_c-i3,lowfil),min(lupfil,nu3_c-i3)
                    t111=t111 + x(i1,1,i2,1,i3+l,1)*a(l)
                    s111=s111 + x(i1,1,i2,1,i3+l,2)*b(l)
                enddo
                y(i1,1,i2,1,i3,1)=t111+s111+cprecr*x(i1,1,i2,1,i3,1)
            endif
            enddo
        enddo
    enddo

       call cpu_time(tcpu1) ; call system_clock(ncount1,ncount_rate,ncount_max)
       tel=dble(ncount1-ncount0)/dble(ncount_rate)
       write(*,'(a40,i4,2(1x,e10.3))') 'FIRST PART',iproc,tcpu1-tcpu0,tel

    do i3=nl3_f,nu3_f
        do i2=nl2_f,nu2_f
            do i1=nl1_f,nu1_f
            if (logrid_f(i1,i2,i3)) then
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                s211=0.d0 ; s121=0.d0 ; s112=0.d0
 ! (1/2) d^2/dx^2
                do l=max(nl1_f-i1,lowfil),min(lupfil,nu1_f-i1)
                    t112=t112 + x(i1+l,1,i2,1,i3,2)*a(l) + x(i1+l,2,i2,1,i3,2)*b(l)
                    t121=t121 + x(i1+l,1,i2,2,i3,1)*a(l) + x(i1+l,2,i2,2,i3,1)*b(l)
                    t122=t122 + x(i1+l,1,i2,2,i3,2)*a(l) + x(i1+l,2,i2,2,i3,2)*b(l)
                    t212=t212 + x(i1+l,1,i2,1,i3,2)*c(l) + x(i1+l,2,i2,1,i3,2)*e(l)
                    t221=t221 + x(i1+l,1,i2,2,i3,1)*c(l) + x(i1+l,2,i2,2,i3,1)*e(l)
                    t222=t222 + x(i1+l,1,i2,2,i3,2)*c(l) + x(i1+l,2,i2,2,i3,2)*e(l)
                enddo
                do l=max(nl1_c-i1,lowfil),min(lupfil,nu1_c-i1)
                    t211=t211 + x(i1+l,1,i2,1,i3,1)*c(l)
                    s211=s211 + x(i1+l,2,i2,1,i3,1)*e(l)
                enddo
 ! + (1/2) d^2/dy^2
                do l=max(nl2_f-i2,lowfil),min(lupfil,nu2_f-i2)
                    t112=t112 + x(i1,1,i2+l,1,i3,2)*a(l) + x(i1,1,i2+l,2,i3,2)*b(l)
                    t211=t211 + x(i1,2,i2+l,1,i3,1)*a(l) + x(i1,2,i2+l,2,i3,1)*b(l)
                    t122=t122 + x(i1,1,i2+l,1,i3,2)*c(l) + x(i1,1,i2+l,2,i3,2)*e(l)
                    t212=t212 + x(i1,2,i2+l,1,i3,2)*a(l) + x(i1,2,i2+l,2,i3,2)*b(l)
                    t221=t221 + x(i1,2,i2+l,1,i3,1)*c(l) + x(i1,2,i2+l,2,i3,1)*e(l)
                    t222=t222 + x(i1,2,i2+l,1,i3,2)*c(l) + x(i1,2,i2+l,2,i3,2)*e(l)
                enddo
                do l=max(nl2_c-i2,lowfil),min(lupfil,nu2_c-i2)
                    t121=t121 + x(i1,1,i2+l,1,i3,1)*c(l)
                    s121=s121 + x(i1,1,i2+l,2,i3,1)*e(l)
                enddo
 ! + (1/2) d^2/dz^2
                do l=max(nl3_f-i3,lowfil),min(lupfil,nu3_f-i3)
                    t121=t121 + x(i1,1,i2,2,i3+l,1)*a(l) + x(i1,1,i2,2,i3+l,2)*b(l)
                    t211=t211 + x(i1,2,i2,1,i3+l,1)*a(l) + x(i1,2,i2,1,i3+l,2)*b(l)
                    t122=t122 + x(i1,1,i2,2,i3+l,1)*c(l) + x(i1,1,i2,2,i3+l,2)*e(l)
                    t212=t212 + x(i1,2,i2,1,i3+l,1)*c(l) + x(i1,2,i2,1,i3+l,2)*e(l)
                    t221=t221 + x(i1,2,i2,2,i3+l,1)*a(l) + x(i1,2,i2,2,i3+l,2)*b(l)
                    t222=t222 + x(i1,2,i2,2,i3+l,1)*c(l) + x(i1,2,i2,2,i3+l,2)*e(l)
                enddo
                do l=max(nl3_c-i3,lowfil),min(lupfil,nu3_c-i3)
                    t112=t112 + x(i1,1,i2,1,i3+l,1)*c(l)
                    s112=s112 + x(i1,1,i2,1,i3+l,2)*e(l)
                enddo


                y(i1,1,i2,1,i3,2)=t112+s112+cprecr*x(i1,1,i2,1,i3,2)
                y(i1,1,i2,2,i3,1)=t121+s121+cprecr*x(i1,1,i2,2,i3,1)
                y(i1,2,i2,1,i3,1)=t211+s211+cprecr*x(i1,2,i2,1,i3,1)
                y(i1,1,i2,2,i3,2)=t122+cprecr*x(i1,1,i2,2,i3,2)
                y(i1,2,i2,1,i3,2)=t212+cprecr*x(i1,2,i2,1,i3,2)
                y(i1,2,i2,2,i3,1)=t221+cprecr*x(i1,2,i2,2,i3,1)
                y(i1,2,i2,2,i3,2)=t222+cprecr*x(i1,2,i2,2,i3,2)
            endif
            enddo
        enddo
    enddo

       call cpu_time(tcpu2) ; call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(*,'(a40,i4,2(1x,e10.3))') 'SECND PART',iproc,tcpu2-tcpu1,tel

    return
end




         subroutine ConvolkineticP(n1,n2,n3, &
               nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c, &
               nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f, &
               hgrid,logrid_c,logrid_f,x,y,ekin)
!   y = y + (kinetic energy operator)x 
    implicit real(kind=8) (a-h,o-z)
    dimension x(nl1_c:nu1_c,2,nl2_c:nu2_c,2,nl3_c:nu3_c,2)
    dimension y(nl1_c:nu1_c,2,nl2_c:nu2_c,2,nl3_c:nu3_c,2)
    logical logrid_c(0:n1,0:n2,0:n3),logrid_f(0:n1,0:n2,0:n3)

    parameter(lowfil=-14,lupfil=14)
    dimension a(lowfil:lupfil),b(lowfil:lupfil),c(lowfil:lupfil),e(lowfil:lupfil)
    scale=-.5d0/hgrid**2
!---------------------------------------------------------------------------
! second derivative filters for Daubechies 16
!  <phi|D^2|phi_i>
    a(0)=   -3.5536922899131901941296809374d0*scale
    a(1)=    2.2191465938911163898794546405d0*scale
    a(2)=   -0.6156141465570069496314853949d0*scale
    a(3)=    0.2371780582153805636239247476d0*scale
    a(4)=   -0.0822663999742123340987663521d0*scale
    a(5)=    0.02207029188482255523789911295638968409d0*scale
    a(6)=   -0.409765689342633823899327051188315485d-2*scale
    a(7)=    0.45167920287502235349480037639758496d-3*scale
    a(8)=   -0.2398228524507599670405555359023135d-4*scale
    a(9)=    2.0904234952920365957922889447361d-6*scale
    a(10)=  -3.7230763047369275848791496973044d-7*scale
    a(11)=  -1.05857055496741470373494132287d-8*scale
    a(12)=  -5.813879830282540547959250667d-11*scale
    a(13)=   2.70800493626319438269856689037647576d-13*scale
    a(14)=  -6.924474940639200152025730585882d-18*scale
    do i=1,14
        a(-i)=a(i)
    enddo
!  <phi|D^2|psi_i>
    c(-14)=     -3.869102413147656535541850057188d-18*scale
    c(-13)=      1.5130616560866154733900029272077362d-13*scale
    c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
    c(-11)=     -5.96264938781402337319841002642d-9*scale
    c(-10)=     -2.1656830629214041470164889350342d-7*scale
    c(-9 )=      8.7969704055286288323596890609625d-7*scale
    c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
    c(-7 )=      0.00021710795484646138591610188464622454d0*scale
    c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
    c(-5 )=      0.00713761218453631422925717625758502986d0*scale
    c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
    c(-3 )=      0.14327329352510759457155821037742893841d0*scale
    c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
    c(-1 )=      0.65703074007121357894896358254040272157d0*scale
    c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
    c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
    c( 2 )=      0.63457035267892488185929915286969303251d0*scale
    c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
    c( 4 )=      0.23370490631751294307619384973520033236d0*scale
    c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
    c( 6 )=      0.0080872029411844780634067667008050127d0*scale
    c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
    c( 8 )=      0.00005075807947289728306309081261461095d0*scale
    c( 9 )=     -4.62561497463184262755416490048242d-6*scale
    c( 10)=      6.3919128513793415587294752371778d-7*scale
    c( 11)=      1.87909235155149902916133888931d-8*scale
    c( 12)=      1.04757345962781829480207861447155543883d-10*scale
    c( 13)=     -4.84665690596158959648731537084025836d-13*scale
    c( 14)=      1.2392629629188986192855777620877d-17*scale
!  <psi|D^2|phi_i>
    do i=-14,14
        b(i)=c(-i)
    enddo
    !<psi|D^2|psi_i>
    e(0)=   -24.875846029392331358907766562d0*scale
    e(1)=   -7.1440597663471719869313377994d0*scale
    e(2)=   -0.04251705323669172315864542163525830944d0*scale
    e(3)=   -0.26995931336279126953587091167128839196d0*scale
    e(4)=    0.08207454169225172612513390763444496516d0*scale
    e(5)=   -0.02207327034586634477996701627614752761d0*scale
    e(6)=    0.00409765642831595181639002667514310145d0*scale
    e(7)=   -0.00045167920287507774929432548999880117d0*scale
    e(8)=    0.00002398228524507599670405555359023135d0*scale
    e(9)=   -2.0904234952920365957922889447361d-6*scale
    e(10)=   3.7230763047369275848791496973044d-7*scale
    e(11)=   1.05857055496741470373494132287d-8*scale
    e(12)=   5.8138798302825405479592506674648873655d-11*scale
    e(13)=  -2.70800493626319438269856689037647576d-13*scale
    e(14)=   6.924474940639200152025730585882d-18*scale
    do i=1,14
        e(-i)=e(i)
    enddo
!---------------------------------------------------------------------------


    ekin=0.d0
                                                            
! Scaling function part

    do i3=nl3_c,nu3_c
        do i2=nl2_c,nu2_c
            do i1=nl1_c,nu1_c
            if (logrid_c(i1,i2,i3)) then
! (1/2) d^2/dx^2
                t111=0.d0
                s111=0.d0
                do l=max(nl1_c-i1,lowfil),min(lupfil,nu1_c-i1)
                    t111=t111 + x(i1+l,1,i2,1,i3,1)*a(l)
                    s111=s111 + x(i1+l,2,i2,1,i3,1)*b(l)
                enddo
! + (1/2) d^2/dy^2
                do l=max(nl2_c-i2,lowfil),min(lupfil,nu2_c-i2)
                    t111=t111 + x(i1,1,i2+l,1,i3,1)*a(l)
                    s111=s111 + x(i1,1,i2+l,2,i3,1)*b(l)
                enddo
! + (1/2) d^2/dz^2
                do l=max(nl3_c-i3,lowfil),min(lupfil,nu3_c-i3)
                    t111=t111 + x(i1,1,i2,1,i3+l,1)*a(l)
                    s111=s111 + x(i1,1,i2,1,i3+l,2)*b(l)
                enddo
                ts111=t111+s111
                ekin=ekin+ts111*x(i1,1,i2,1,i3,1)
                y(i1,1,i2,1,i3,1)=ts111+y(i1,1,i2,1,i3,1)
            endif
            enddo
        enddo
    enddo

    do i3=nl3_f,nu3_f
        do i2=nl2_f,nu2_f
            do i1=nl1_f,nu1_f
            if (logrid_f(i1,i2,i3)) then
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                s211=0.d0 ; s121=0.d0 ; s112=0.d0
 ! (1/2) d^2/dx^2
                do l=max(nl1_f-i1,lowfil),min(lupfil,nu1_f-i1)
                    t112=t112 + x(i1+l,1,i2,1,i3,2)*a(l) + x(i1+l,2,i2,1,i3,2)*b(l)
                    t121=t121 + x(i1+l,1,i2,2,i3,1)*a(l) + x(i1+l,2,i2,2,i3,1)*b(l)
                    t122=t122 + x(i1+l,1,i2,2,i3,2)*a(l) + x(i1+l,2,i2,2,i3,2)*b(l)
                    t212=t212 + x(i1+l,1,i2,1,i3,2)*c(l) + x(i1+l,2,i2,1,i3,2)*e(l)
                    t221=t221 + x(i1+l,1,i2,2,i3,1)*c(l) + x(i1+l,2,i2,2,i3,1)*e(l)
                    t222=t222 + x(i1+l,1,i2,2,i3,2)*c(l) + x(i1+l,2,i2,2,i3,2)*e(l)
                enddo
                do l=max(nl1_c-i1,lowfil),min(lupfil,nu1_c-i1)
                    t211=t211 + x(i1+l,1,i2,1,i3,1)*c(l)
                    s211=s211 + x(i1+l,2,i2,1,i3,1)*e(l)
                enddo
 ! + (1/2) d^2/dy^2
                do l=max(nl2_f-i2,lowfil),min(lupfil,nu2_f-i2)
                    t112=t112 + x(i1,1,i2+l,1,i3,2)*a(l) + x(i1,1,i2+l,2,i3,2)*b(l)
                    t211=t211 + x(i1,2,i2+l,1,i3,1)*a(l) + x(i1,2,i2+l,2,i3,1)*b(l)
                    t122=t122 + x(i1,1,i2+l,1,i3,2)*c(l) + x(i1,1,i2+l,2,i3,2)*e(l)
                    t212=t212 + x(i1,2,i2+l,1,i3,2)*a(l) + x(i1,2,i2+l,2,i3,2)*b(l)
                    t221=t221 + x(i1,2,i2+l,1,i3,1)*c(l) + x(i1,2,i2+l,2,i3,1)*e(l)
                    t222=t222 + x(i1,2,i2+l,1,i3,2)*c(l) + x(i1,2,i2+l,2,i3,2)*e(l)
                enddo
                do l=max(nl2_c-i2,lowfil),min(lupfil,nu2_c-i2)
                    t121=t121 + x(i1,1,i2+l,1,i3,1)*c(l)
                    s121=s121 + x(i1,1,i2+l,2,i3,1)*e(l)
                enddo
 ! + (1/2) d^2/dz^2
                do l=max(nl3_f-i3,lowfil),min(lupfil,nu3_f-i3)
                    t121=t121 + x(i1,1,i2,2,i3+l,1)*a(l) + x(i1,1,i2,2,i3+l,2)*b(l)
                    t211=t211 + x(i1,2,i2,1,i3+l,1)*a(l) + x(i1,2,i2,1,i3+l,2)*b(l)
                    t122=t122 + x(i1,1,i2,2,i3+l,1)*c(l) + x(i1,1,i2,2,i3+l,2)*e(l)
                    t212=t212 + x(i1,2,i2,1,i3+l,1)*c(l) + x(i1,2,i2,1,i3+l,2)*e(l)
                    t221=t221 + x(i1,2,i2,2,i3+l,1)*a(l) + x(i1,2,i2,2,i3+l,2)*b(l)
                    t222=t222 + x(i1,2,i2,2,i3+l,1)*c(l) + x(i1,2,i2,2,i3+l,2)*e(l)
                enddo
                do l=max(nl3_c-i3,lowfil),min(lupfil,nu3_c-i3)
                    t112=t112 + x(i1,1,i2,1,i3+l,1)*c(l)
                    s112=s112 + x(i1,1,i2,1,i3+l,2)*e(l)
                enddo


                ts112=t112+s112
                ekin=ekin+ts112*x(i1,1,i2,1,i3,2)
                y(i1,1,i2,1,i3,2)=ts112+y(i1,1,i2,1,i3,2)

                ts121=t121+s121
                ekin=ekin+ts121*x(i1,1,i2,2,i3,1)
                y(i1,1,i2,2,i3,1)=ts121+y(i1,1,i2,2,i3,1)

                ts211=t211+s211
                ekin=ekin+ts211*x(i1,2,i2,1,i3,1)
                y(i1,2,i2,1,i3,1)=ts211+y(i1,2,i2,1,i3,1)

                ekin=ekin+t122*x(i1,1,i2,2,i3,2)
                y(i1,1,i2,2,i3,2)=t122+y(i1,1,i2,2,i3,2)

                ekin=ekin+t212*x(i1,2,i2,1,i3,2)
                y(i1,2,i2,1,i3,2)=t212+y(i1,2,i2,1,i3,2)

                ekin=ekin+t221*x(i1,2,i2,2,i3,1)
                y(i1,2,i2,2,i3,1)=t221+y(i1,2,i2,2,i3,1)

                ekin=ekin+t222*x(i1,2,i2,2,i3,2)
                y(i1,2,i2,2,i3,2)=t222+y(i1,2,i2,2,i3,2)
            endif
            enddo
        enddo
    enddo

    return
end

