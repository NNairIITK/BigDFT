!> @file
!!  New convolution routines
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Applies the operator (KE+cprecr*I)*x=y
!! array x is input, array y is output
!! See also the optimized version (apply_hp_sd_optim)
subroutine apply_hp_sd(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,x,y,psig_in,psig_out,modul1,modul2,modul3,a,b,c,e) !n(c) hx,hy,hz (arg:11,12,13)
  use module_base
  implicit none
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: cprecr
  !n(c) real(gp), intent(in) :: hx,hy,hz
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer,intent(in)::modul1(lowfil:n1+lupfil)
  integer,intent(in)::modul2(lowfil:n2+lupfil)
  integer,intent(in)::modul3(lowfil:n3+lupfil)
  real(gp),intent(in)::a(lowfil:lupfil,3)
  real(gp),intent(in)::b(lowfil:lupfil,3)
  real(gp),intent(in)::c(lowfil:lupfil,3)
  real(gp),intent(in)::e(lowfil:lupfil,3)
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f)  

  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f)
  !local variables
  !n(c) real(gp), dimension(3) :: hgrid
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)) :: psig_in,psig_out

  call uncompress_sd(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       x(1),x(nvctr_c+1),psig_in)

  !n(c) hgrid(1)=hx
  !n(c) hgrid(2)=hy
  !n(c) hgrid(3)=hz
  call convolut_kinetic_per_sdc(n1,n2,n3,psig_in,psig_out,cprecr,modul1,modul2,modul3,a,b,c,e)

  call compress_sd(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
       psig_out,y(1),y(nvctr_c+1))
END SUBROUTINE apply_hp_sd



!>   See also the optimized version (apply_hp_scal_optim)
subroutine apply_hp_scal(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,x,y,psig_in,psig_out,modul1,modul2,modul3,a,b,c,e,scal) !n(c) hx,hy,hz (arg:11,12,13)
  use module_base
  implicit none
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp),intent(in)::scal(0:7)
  real(gp), intent(in) :: cprecr
  !n(c) real(gp), intent(in) :: hx,hy,hz
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer,intent(in) :: modul1(lowfil:n1+lupfil)
  integer,intent(in) :: modul2(lowfil:n2+lupfil)
  integer,intent(in) :: modul3(lowfil:n3+lupfil)
  real(gp),intent(in) :: a(lowfil:lupfil,3)
  real(gp),intent(in) :: b(lowfil:lupfil,3)
  real(gp),intent(in) :: c(lowfil:lupfil,3)
  real(gp),intent(in) :: e(lowfil:lupfil,3)
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f)  

  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f)
  !local variables
  !n(c) real(gp), dimension(3) :: hgrid
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)) :: psig_in,psig_out

  call uncompress_sd_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       x(1),x(nvctr_c+1),psig_in,scal)

  !n(c) hgrid(1)=hx
  !n(c) hgrid(2)=hy
  !n(c) hgrid(3)=hz
  call convolut_kinetic_per_sdc(n1,n2,n3,psig_in,psig_out,cprecr,modul1,modul2,modul3,a,b,c,e)

  call compress_sd_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
       psig_out,y(1),y(nvctr_c+1),scal)
END SUBROUTINE apply_hp_scal



!>   Applies the kinetic energy operator onto x to get y. Works for periodic BC
!!   See also the optimized version (convolut_kinteic_per_sdc_optim)
subroutine convolut_kinetic_per_sdc(n1,n2,n3,x,y,cprecr,modul1,modul2,modul3,a,b,c,e)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::cprecr
  real(wp), dimension(8,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(8,0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,l,j
  real(wp) :: tt111,tt112,tt121,tt122,tt211,tt212,tt221,tt222
  !real(wp), dimension(3) :: scale
  integer,intent(in)::modul1(lowfil:n1+lupfil)
  integer,intent(in)::modul2(lowfil:n2+lupfil)
  integer,intent(in)::modul3(lowfil:n3+lupfil)
  real(gp),intent(in)::a(lowfil:lupfil,3)
  real(gp),intent(in)::b(lowfil:lupfil,3)
  real(gp),intent(in)::c(lowfil:lupfil,3)
  real(gp),intent(in)::e(lowfil:lupfil,3)
  
  !real(gp) :: tel
  !integer :: ncount0,ncount1,ncount2,ncount_rate,ncount_max
  !integer :: mflop1,mflop3
  
  ! filter length:29
  ! 8: wavelets+scfunction
  ! 4: flops for one iteration
!  mflop1=(n1+1)*(n2+1)*(n3+1)*29*8*4*2 ! convolution in the x and y direction
!  mflop3=(n1+1)*(n2+1)*(n3+1)*29*8*4   ! convolution in the z       direction
!  call system_clock(ncount0,ncount_rate,ncount_max)
  
!$omp parallel default(private) &
!$omp shared(x,y,n1,n2,n3,cprecr,modul1,modul2,modul3,a,b,c,e)

!$omp do
  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt111=x(1,i1,i2,i3)*cprecr
           tt211=x(2,i1,i2,i3)*cprecr
           tt121=x(3,i1,i2,i3)*cprecr
           tt221=x(4,i1,i2,i3)*cprecr
           tt112=x(5,i1,i2,i3)*cprecr
           tt212=x(6,i1,i2,i3)*cprecr
           tt122=x(7,i1,i2,i3)*cprecr
           tt222=x(8,i1,i2,i3)*cprecr
           do l=lowfil,lupfil
              j=modul1(i1+l)
              tt111=tt111+x(1,j,i2,i3)*a(l,1)+x(2,j,i2,i3)*b(l,1)
              tt211=tt211+x(2,j,i2,i3)*e(l,1)+x(1,j,i2,i3)*c(l,1)
              tt121=tt121+x(3,j,i2,i3)*a(l,1)+x(4,j,i2,i3)*b(l,1)
              tt221=tt221+x(4,j,i2,i3)*e(l,1)+x(3,j,i2,i3)*c(l,1)
              tt112=tt112+x(5,j,i2,i3)*a(l,1)+x(6,j,i2,i3)*b(l,1)
              tt212=tt212+x(6,j,i2,i3)*e(l,1)+x(5,j,i2,i3)*c(l,1)
              tt122=tt122+x(7,j,i2,i3)*a(l,1)+x(8,j,i2,i3)*b(l,1)
              tt222=tt222+x(8,j,i2,i3)*e(l,1)+x(7,j,i2,i3)*c(l,1)
           enddo
           y(1,i1,i2,i3)=tt111
           y(2,i1,i2,i3)=tt211
           y(3,i1,i2,i3)=tt121
           y(4,i1,i2,i3)=tt221
           y(5,i1,i2,i3)=tt112
           y(6,i1,i2,i3)=tt212
           y(7,i1,i2,i3)=tt122
           y(8,i1,i2,i3)=tt222
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt111=0.e0_wp
           tt211=0.e0_wp
           tt121=0.e0_wp
           tt221=0.e0_wp
           tt112=0.e0_wp
           tt212=0.e0_wp
           tt122=0.e0_wp
           tt222=0.e0_wp
           do l=lowfil,lupfil
              j=modul2(i2+l)
              tt111=tt111+x(1,i1,j,i3)*a(l,2)+x(3,i1,j,i3)*b(l,2)
              tt211=tt211+x(2,i1,j,i3)*a(l,2)+x(4,i1,j,i3)*b(l,2)
              tt121=tt121+x(3,i1,j,i3)*e(l,2)+x(1,i1,j,i3)*c(l,2)
              tt221=tt221+x(4,i1,j,i3)*e(l,2)+x(2,i1,j,i3)*c(l,2)
              tt112=tt112+x(5,i1,j,i3)*a(l,2)+x(7,i1,j,i3)*b(l,2)
              tt212=tt212+x(6,i1,j,i3)*a(l,2)+x(8,i1,j,i3)*b(l,2)
              tt122=tt122+x(7,i1,j,i3)*e(l,2)+x(5,i1,j,i3)*c(l,2)
              tt222=tt222+x(8,i1,j,i3)*e(l,2)+x(6,i1,j,i3)*c(l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt111
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt211
           y(3,i1,i2,i3)=y(3,i1,i2,i3)+tt121
           y(4,i1,i2,i3)=y(4,i1,i2,i3)+tt221
           y(5,i1,i2,i3)=y(5,i1,i2,i3)+tt112
           y(6,i1,i2,i3)=y(6,i1,i2,i3)+tt212
           y(7,i1,i2,i3)=y(7,i1,i2,i3)+tt122
           y(8,i1,i2,i3)=y(8,i1,i2,i3)+tt222
        enddo
     enddo
     
  enddo

  !$omp end do  

!  call system_clock(ncount1,ncount_rate,ncount_max)
!  tel=dble(ncount1-ncount0)/dble(ncount_rate)
!  write(97,'(a40,1x,e10.3,1x,f6.1)') 'x,y:',tel,1.d-6*mflop1/tel
 
!$omp do   
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt111=0.e0_wp
           tt211=0.e0_wp
           tt121=0.e0_wp
           tt221=0.e0_wp
           tt112=0.e0_wp
           tt212=0.e0_wp
           tt122=0.e0_wp
           tt222=0.e0_wp
           do l=lowfil,lupfil
              j=modul3(i3+l)
              tt111=tt111+x(1,i1,i2,j)*a(l,3)+x(5,i1,i2,j)*b(l,3)
              tt211=tt211+x(2,i1,i2,j)*a(l,3)+x(6,i1,i2,j)*b(l,3)
              tt121=tt121+x(3,i1,i2,j)*a(l,3)+x(7,i1,i2,j)*b(l,3)
              tt221=tt221+x(4,i1,i2,j)*a(l,3)+x(8,i1,i2,j)*b(l,3)
              tt112=tt112+x(5,i1,i2,j)*e(l,3)+x(1,i1,i2,j)*c(l,3)
              tt212=tt212+x(6,i1,i2,j)*e(l,3)+x(2,i1,i2,j)*c(l,3)
              tt122=tt122+x(7,i1,i2,j)*e(l,3)+x(3,i1,i2,j)*c(l,3)
              tt222=tt222+x(8,i1,i2,j)*e(l,3)+x(4,i1,i2,j)*c(l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt111
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt211
           y(3,i1,i2,i3)=y(3,i1,i2,i3)+tt121
           y(4,i1,i2,i3)=y(4,i1,i2,i3)+tt221
           y(5,i1,i2,i3)=y(5,i1,i2,i3)+tt112
           y(6,i1,i2,i3)=y(6,i1,i2,i3)+tt212
           y(7,i1,i2,i3)=y(7,i1,i2,i3)+tt122
           y(8,i1,i2,i3)=y(8,i1,i2,i3)+tt222
        enddo
     enddo
  enddo

  !$omp end do
  !$omp end parallel  
 
!  call system_clock(ncount2,ncount_rate,ncount_max)
!  tel=dble(ncount2-ncount1)/dble(ncount_rate)
!  write(97,'(a40,1x,e10.3,1x,f6.1)') 'z:',tel,1.d-6*mflop3/tel

END SUBROUTINE convolut_kinetic_per_sdc



subroutine prepare_sdc(n1,n2,n3,modul1,modul2,modul3,a,b,c,e,hx,hy,hz)
  use module_base
  implicit none
  integer,intent(in) :: n1,n2,n3
  real(gp),intent(in) :: hx,hy,hz

  integer, parameter :: lowfil=-14,lupfil=14

  integer,intent(out) :: modul1(lowfil:n1+lupfil)
  integer,intent(out) :: modul2(lowfil:n2+lupfil)
  integer,intent(out) :: modul3(lowfil:n3+lupfil)
  real(gp),intent(out) :: a(lowfil:lupfil,3)
  real(gp),intent(out) :: b(lowfil:lupfil,3)
  real(gp),intent(out) :: c(lowfil:lupfil,3)
  real(gp),intent(out) :: e(lowfil:lupfil,3)

  real(gp) :: hgrid(3)
  integer :: i
  real(gp) :: scale(3)

  call fill_mod_arr(modul1,lowfil,n1+lupfil,n1+1)
  call fill_mod_arr(modul2,lowfil,n2+lupfil,n2+1)
  call fill_mod_arr(modul3,lowfil,n3+lupfil,n3+1)
  
  hgrid(1)=hx
  hgrid(2)=hy
  hgrid(3)=hz

  scale(:)=real(-.5_gp/hgrid(:)**2,wp)
  
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0,:)=   -3.5536922899131901941296809374_wp*scale(:)
  a(1,:)=    2.2191465938911163898794546405_wp*scale(:)
  a(2,:)=   -0.6156141465570069496314853949_wp*scale(:)
  a(3,:)=    0.2371780582153805636239247476_wp*scale(:)
  a(4,:)=   -0.0822663999742123340987663521_wp*scale(:)
  a(5,:)=    0.02207029188482255523789911295638968409_wp*scale(:)
  a(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  a(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  a(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  a(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  a(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  a(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  a(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  a(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  a(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     a(-i,:)=a(i,:)
  enddo
  !  <phi|D^2|psi_i>
  c(-14,:)=     -3.869102413147656535541850057188e-18_wp*scale(:)
  c(-13,:)=      1.5130616560866154733900029272077362e-13_wp*scale(:)
  c(-12,:)=     -3.2264702314010525539061647271983988409e-11_wp*scale(:)
  c(-11,:)=     -5.96264938781402337319841002642e-9_wp*scale(:)
  c(-10,:)=     -2.1656830629214041470164889350342e-7_wp*scale(:)
  c(-9 ,:)=      8.7969704055286288323596890609625e-7_wp*scale(:)
  c(-8 ,:)=     -0.00001133456724516819987751818232711775_wp*scale(:)
  c(-7 ,:)=      0.00021710795484646138591610188464622454_wp*scale(:)
  c(-6 ,:)=     -0.0021356291838797986414312219042358542_wp*scale(:)
  c(-5 ,:)=      0.00713761218453631422925717625758502986_wp*scale(:)
  c(-4 ,:)=     -0.0284696165863973422636410524436931061_wp*scale(:)
  c(-3 ,:)=      0.14327329352510759457155821037742893841_wp*scale(:)
  c(-2 ,:)=     -0.42498050943780130143385739554118569733_wp*scale(:)
  c(-1 ,:)=      0.65703074007121357894896358254040272157_wp*scale(:)
  c( 0 ,:)=     -0.42081655293724308770919536332797729898_wp*scale(:)
  c( 1 ,:)=     -0.21716117505137104371463587747283267899_wp*scale(:)
  c( 2 ,:)=      0.63457035267892488185929915286969303251_wp*scale(:)
  c( 3 ,:)=     -0.53298223962800395684936080758073568406_wp*scale(:)
  c( 4 ,:)=      0.23370490631751294307619384973520033236_wp*scale(:)
  c( 5 ,:)=     -0.05657736973328755112051544344507997075_wp*scale(:)
  c( 6 ,:)=      0.0080872029411844780634067667008050127_wp*scale(:)
  c( 7 ,:)=     -0.00093423623304808664741804536808932984_wp*scale(:)
  c( 8 ,:)=      0.00005075807947289728306309081261461095_wp*scale(:)
  c( 9 ,:)=     -4.62561497463184262755416490048242e-6_wp*scale(:)
  c( 10,:)=      6.3919128513793415587294752371778e-7_wp*scale(:)
  c( 11,:)=      1.87909235155149902916133888931e-8_wp*scale(:)
  c( 12,:)=      1.04757345962781829480207861447155543883e-10_wp*scale(:)
  c( 13,:)=     -4.84665690596158959648731537084025836e-13_wp*scale(:)
  c( 14,:)=      1.2392629629188986192855777620877e-17_wp*scale(:)
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i,:)=c(-i,:)
  enddo
  !<psi|D^2|psi_i>
  e(0,:)=   -24.875846029392331358907766562_wp*scale(:)
  e(1,:)=   -7.1440597663471719869313377994_wp*scale(:)
  e(2,:)=   -0.04251705323669172315864542163525830944_wp*scale(:)
  e(3,:)=   -0.26995931336279126953587091167128839196_wp*scale(:)
  e(4,:)=    0.08207454169225172612513390763444496516_wp*scale(:)
  e(5,:)=   -0.02207327034586634477996701627614752761_wp*scale(:)
  e(6,:)=    0.00409765642831595181639002667514310145_wp*scale(:)
  e(7,:)=   -0.00045167920287507774929432548999880117_wp*scale(:)
  e(8,:)=    0.00002398228524507599670405555359023135_wp*scale(:)
  e(9,:)=   -2.0904234952920365957922889447361e-6_wp*scale(:)
  e(10,:)=   3.7230763047369275848791496973044e-7_wp*scale(:)
  e(11,:)=   1.05857055496741470373494132287e-8_wp*scale(:)
  e(12,:)=   5.8138798302825405479592506674648873655e-11_wp*scale(:)
  e(13,:)=  -2.70800493626319438269856689037647576e-13_wp*scale(:)
  e(14,:)=   6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     e(-i,:)=e(i,:)
  enddo
END SUBROUTINE prepare_sdc


subroutine convolut_kinetic_hyb_T(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,x_c,x_f,y_c,y_f,kstrten,x_f1,x_f2,x_f3,ibyz,ibxz,ibxy)
  !   y = y+(kinetic energy operator)x 
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer,intent(in)::ibyz(2,0:n2,0:n3),ibxz(2,0:n1,0:n3),ibxy(2,0:n1,0:n2)
  real(gp), intent(in) :: hgrid(3)
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: x_f3
  !real(gp), intent(out) :: ekinout
  real(gp), dimension(6), intent(out) :: kstrten
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: l
  real(wp) :: dyi,t112,t121,t122,t212,t221,t222,t211
  real(wp) :: scale(3)
  real(gp) :: ekin1,ekin2,ekin3,ekin4,ekin5,ekin6,ekin7,ekin8,ekin9
  !real(kind=8) :: tel
  real(wp), dimension(lowfil:lupfil,3) :: a,b,c,e
  integer :: ii1,ii2,ii3 ! for hybrid convolutions
  !n(c) integer :: ic
  integer :: mod_arr1(lowfil:n1+lupfil)
  integer :: mod_arr2(lowfil:n2+lupfil)
  integer :: mod_arr3(lowfil:n3+lupfil)

  !ekinout=0._gp
  kstrten=0._gp
!$omp parallel default(private) shared(mod_arr1,mod_arr2,mod_arr3,n1,n2,n3,hgrid,x_c,x_f,x_f1,x_f2,x_f3) &
!$omp shared(kstrten,y_c,y_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,ibyz,ibxz,ibxy)
        !n(c) ic=1
   call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
   call fill_mod_arr(mod_arr2,lowfil,n2+lupfil,n2+1)
   call fill_mod_arr(mod_arr3,lowfil,n3+lupfil,n3+1)

  scale(:)=-.5_wp/hgrid(:)**2
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0,:)=   -3.5536922899131901941296809374_wp*scale(:)
  a(1,:)=    2.2191465938911163898794546405_wp*scale(:)
  a(2,:)=   -0.6156141465570069496314853949_wp*scale(:)
  a(3,:)=    0.2371780582153805636239247476_wp*scale(:)
  a(4,:)=   -0.0822663999742123340987663521_wp*scale(:)
  a(5,:)=    0.02207029188482255523789911295638968409_wp*scale(:)
  a(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  a(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  a(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  a(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  a(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  a(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  a(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  a(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  a(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     a(-i,:)=a(i,:)
  enddo
  !  <phi|D^2|psi_i>
  c(-14,:)=     -3.869102413147656535541850057188e-18_wp*scale(:)
  c(-13,:)=      1.5130616560866154733900029272077362e-13_wp*scale(:)
  c(-12,:)=     -3.2264702314010525539061647271983988409e-11_wp*scale(:)
  c(-11,:)=     -5.96264938781402337319841002642e-9_wp*scale(:)
  c(-10,:)=     -2.1656830629214041470164889350342e-7_wp*scale(:)
  c(-9 ,:)=      8.7969704055286288323596890609625e-7_wp*scale(:)
  c(-8 ,:)=     -0.00001133456724516819987751818232711775_wp*scale(:)
  c(-7 ,:)=      0.00021710795484646138591610188464622454_wp*scale(:)
  c(-6 ,:)=     -0.0021356291838797986414312219042358542_wp*scale(:)
  c(-5 ,:)=      0.00713761218453631422925717625758502986_wp*scale(:)
  c(-4 ,:)=     -0.0284696165863973422636410524436931061_wp*scale(:)
  c(-3 ,:)=      0.14327329352510759457155821037742893841_wp*scale(:)
  c(-2 ,:)=     -0.42498050943780130143385739554118569733_wp*scale(:)
  c(-1 ,:)=      0.65703074007121357894896358254040272157_wp*scale(:)
  c( 0 ,:)=     -0.42081655293724308770919536332797729898_wp*scale(:)
  c( 1 ,:)=     -0.21716117505137104371463587747283267899_wp*scale(:)
  c( 2 ,:)=      0.63457035267892488185929915286969303251_wp*scale(:)
  c( 3 ,:)=     -0.53298223962800395684936080758073568406_wp*scale(:)
  c( 4 ,:)=      0.23370490631751294307619384973520033236_wp*scale(:)
  c( 5 ,:)=     -0.05657736973328755112051544344507997075_wp*scale(:)
  c( 6 ,:)=      0.0080872029411844780634067667008050127_wp*scale(:)
  c( 7 ,:)=     -0.00093423623304808664741804536808932984_wp*scale(:)
  c( 8 ,:)=      0.00005075807947289728306309081261461095_wp*scale(:)
  c( 9 ,:)=     -4.62561497463184262755416490048242e-6_wp*scale(:)
  c( 10,:)=      6.3919128513793415587294752371778e-7_wp*scale(:)
  c( 11,:)=      1.87909235155149902916133888931e-8_wp*scale(:)
  c( 12,:)=      1.04757345962781829480207861447155543883e-10_wp*scale(:)
  c( 13,:)=     -4.84665690596158959648731537084025836e-13_wp*scale(:)
  c( 14,:)=      1.2392629629188986192855777620877e-17_wp*scale(:)
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i,:)=c(-i,:)
  enddo
  !<psi|D^2|psi_i>
  e(0,:)=   -24.875846029392331358907766562_wp*scale(:)
  e(1,:)=   -7.1440597663471719869313377994_wp*scale(:)
  e(2,:)=   -0.04251705323669172315864542163525830944_wp*scale(:)
  e(3,:)=   -0.26995931336279126953587091167128839196_wp*scale(:)
  e(4,:)=    0.08207454169225172612513390763444496516_wp*scale(:)
  e(5,:)=   -0.02207327034586634477996701627614752761_wp*scale(:)
  e(6,:)=    0.00409765642831595181639002667514310145_wp*scale(:)
  e(7,:)=   -0.00045167920287507774929432548999880117_wp*scale(:)
  e(8,:)=    0.00002398228524507599670405555359023135_wp*scale(:)
  e(9,:)=   -2.0904234952920365957922889447361e-6_wp*scale(:)
  e(10,:)=   3.7230763047369275848791496973044e-7_wp*scale(:)
  e(11,:)=   1.05857055496741470373494132287e-8_wp*scale(:)
  e(12,:)=   5.8138798302825405479592506674648873655e-11_wp*scale(:)
  e(13,:)=  -2.70800493626319438269856689037647576e-13_wp*scale(:)
  e(14,:)=   6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     e(-i,:)=e(i,:)
  enddo

  ekin1=0._gp
  ekin2=0._gp
  ekin3=0._gp
  !---------------------------------------------------------------------------

  ! Scaling function part: 
  ! periodic convolutions for scaling functions, 
  ! free convolutions for wavelets


  call conv_kin_x1(x_c,y_c,n1,n2,n3,ekin1,mod_arr1,a)   
  call conv_kin_y1(x_c,y_c,n1,n2,n3,ekin2,mod_arr2,a)   
  call conv_kin_z1(x_c,y_c,n1,n2,n3,ekin3,mod_arr3,a)   
  
  
  ekin4=0._gp
  ekin5=0._gp
  ekin6=0._gp
  ekin7=0._gp
  ekin8=0._gp
  ekin9=0._gp
!!  ! (1/2) d^2/dx^2
!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        ! wavelet-scaling function:
        ! scaling function everywhere, also outside of the box
        ! wavelet restricted
            do i1=ibyz(1,i2,i3)-lupfil,ibyz(2,i2,i3)-lowfil
            ii1=mod_arr1(i1)
                dyi=0._wp
            do l=max(lowfil,ibyz(1,i2,i3)-i1),min(lupfil,ibyz(2,i2,i3)-i1)
                    dyi=dyi + x_f1(i1+l,i2,i3)*b(l,1)
                enddo
                y_c(ii1,i2,i3)=y_c(ii1,i2,i3)+dyi
                ekin4=ekin4+dyi*x_c(ii1,i2,i3)
            enddo
            
         ! scaling function-wavelet
         ! wavelet restricted, scaling function everywhere
            do i1=ibyz(1,i2,i3),ibyz(2,i2,i3)
               t211=0._wp 
            do l=lowfil,lupfil
                t=mod_arr1(i1+l)
                   t211=t211 + x_c(t,i2,i3)*c(l,1)
               enddo
               y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
               ekin4=ekin4+t211*x_f(1,i1,i2,i3)
            enddo
     enddo
  enddo
!$omp enddo  
  
!!
!!  ! + (1/2) d^2/dy^2
!!
!$omp do  
   do i3=nfl3,nfu3
      do i1=nfl1,nfu1
            do i2=ibxz(1,i1,i3)-lupfil,ibxz(2,i1,i3)-lowfil
            ii2=mod_arr2(i2)
                dyi=0._wp
            do l=max(lowfil,ibxz(1,i1,i3)-i2),min(lupfil,ibxz(2,i1,i3)-i2)
                   dyi=dyi + x_f2(i2+l,i1,i3)*b(l,2)
                enddo
                y_c(i1,ii2,i3)=y_c(i1,ii2,i3)+dyi
                ekin5=ekin5+dyi*x_c(i1,ii2,i3)
            enddo
            
            do i2=ibxz(1,i1,i3),ibxz(2,i1,i3)
                t121=0._wp 
            do l=lowfil,lupfil
               t=mod_arr2(i2+l)   
                    t121=t121 + x_c(i1,t,i3)*c(l,2)
                enddo
                y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
                ekin5=ekin5+t121*x_f(2,i1,i2,i3)
            enddo
      enddo
   enddo
!$omp enddo  
  
!!
!!  ! + (1/2) d^2/dz^2
!!
!$omp do  
   do i2=nfl2,nfu2
      do i1=nfl1,nfu1
            do i3=ibxy(1,i1,i2)-lupfil,ibxy(2,i1,i2)-lowfil
                dyi=0._wp
             ii3=mod_arr3(i3)
            do l=max(lowfil,ibxy(1,i1,i2)-i3),min(lupfil,ibxy(2,i1,i2)-i3)
                   dyi=dyi + x_f3(i3+l,i1,i2)*b(l,3)
               enddo
                y_c(i1,i2,ii3)=y_c(i1,i2,ii3)+dyi
                ekin6=ekin6+dyi*x_c(i1,i2,ii3)
            enddo

            do i3=ibxy(1,i1,i2),ibxy(2,i1,i2)
                t112=0._wp 
            do l=lowfil,lupfil
               t=mod_arr3(i3+l)   
                   t112=t112 + x_c(i1,i2,t)*c(l,3)
               enddo
                y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
                ekin6=ekin6+t112*x_f(4,i1,i2,i3)
            enddo
      enddo
   enddo
!$omp enddo   
  
  
  ! wavelet part
  ! completely similar to the free case
  ! (1/2) d^2/dx^2
!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz(1,i2,i3),ibyz(2,i2,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(ibyz(1,i2,i3)-i1,lowfil),min(lupfil,ibyz(2,i2,i3)-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l,1) + x_f(5,i1+l,i2,i3)*b(l,1)
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l,1) + x_f(3,i1+l,i2,i3)*b(l,1)
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l,1) + x_f(7,i1+l,i2,i3)*b(l,1)
              t212=t212 + x_f(4,i1+l,i2,i3)*c(l,1) + x_f(5,i1+l,i2,i3)*e(l,1)
              t221=t221 + x_f(2,i1+l,i2,i3)*c(l,1) + x_f(3,i1+l,i2,i3)*e(l,1)
              t222=t222 + x_f(6,i1+l,i2,i3)*c(l,1) + x_f(7,i1+l,i2,i3)*e(l,1)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l,1)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin7=ekin7+t112*x_f(4,i1,i2,i3)
           ekin7=ekin7+t121*x_f(2,i1,i2,i3)
           ekin7=ekin7+t211*x_f(1,i1,i2,i3)
           ekin7=ekin7+t122*x_f(6,i1,i2,i3)
           ekin7=ekin7+t212*x_f(5,i1,i2,i3)
           ekin7=ekin7+t221*x_f(3,i1,i2,i3)
           ekin7=ekin7+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
!$omp enddo
  ! + (1/2) d^2/dy^2
!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz(1,i1,i3),ibxz(2,i1,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(ibxz(1,i1,i3)-i2,lowfil),min(lupfil,ibxz(2,i1,i3)-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l,2) + x_f(6,i1,i2+l,i3)*b(l,2)
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l,2) + x_f(3,i1,i2+l,i3)*b(l,2)
              t122=t122 + x_f(4,i1,i2+l,i3)*c(l,2) + x_f(6,i1,i2+l,i3)*e(l,2)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l,2) + x_f(7,i1,i2+l,i3)*b(l,2)
              t221=t221 + x_f(1,i1,i2+l,i3)*c(l,2) + x_f(3,i1,i2+l,i3)*e(l,2)
              t222=t222 + x_f(5,i1,i2+l,i3)*c(l,2) + x_f(7,i1,i2+l,i3)*e(l,2)
              t121=t121 + x_f(2,i1,i2+l,i3)*e(l,2)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin8=ekin8+t112*x_f(4,i1,i2,i3)
           ekin8=ekin8+t121*x_f(2,i1,i2,i3)
           ekin8=ekin8+t211*x_f(1,i1,i2,i3)
           ekin8=ekin8+t122*x_f(6,i1,i2,i3)
           ekin8=ekin8+t212*x_f(5,i1,i2,i3)
           ekin8=ekin8+t221*x_f(3,i1,i2,i3)
           ekin8=ekin8+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
!$omp enddo
  ! + (1/2) d^2/dz^2
!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy(1,i1,i2),ibxy(2,i1,i2)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(ibxy(1,i1,i2)-i3,lowfil),min(lupfil,ibxy(2,i1,i2)-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l,3) + x_f(6,i1,i2,i3+l)*b(l,3)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l,3) + x_f(5,i1,i2,i3+l)*b(l,3)
              t122=t122 + x_f(2,i1,i2,i3+l)*c(l,3) + x_f(6,i1,i2,i3+l)*e(l,3)
              t212=t212 + x_f(1,i1,i2,i3+l)*c(l,3) + x_f(5,i1,i2,i3+l)*e(l,3)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l,3) + x_f(7,i1,i2,i3+l)*b(l,3)
              t222=t222 + x_f(3,i1,i2,i3+l)*c(l,3) + x_f(7,i1,i2,i3+l)*e(l,3)
              t112=t112 + x_f(4,i1,i2,i3+l)*e(l,3)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin9=ekin9+t112*x_f(4,i1,i2,i3)
           ekin9=ekin9+t121*x_f(2,i1,i2,i3)
           ekin9=ekin9+t211*x_f(1,i1,i2,i3)
           ekin9=ekin9+t122*x_f(6,i1,i2,i3)
           ekin9=ekin9+t212*x_f(5,i1,i2,i3)
           ekin9=ekin9+t221*x_f(3,i1,i2,i3)
           ekin9=ekin9+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo
!$omp enddo
!$omp critical
  !accumulate for multithread
  kstrten(1)=kstrten(1)+ekin1+ekin4+ekin7
  kstrten(2)=kstrten(2)+ekin2+ekin5+ekin8
  kstrten(3)=kstrten(3)+ekin3+ekin6+ekin9
  !ekinout=ekinout+ekin1+ekin2+ekin3+ekin4+ekin5+ekin6+ekin7+ekin8+ekin9
!$omp end critical
!$omp end parallel
END SUBROUTINE convolut_kinetic_hyb_T

  subroutine conv_kin_x1(x,y,n1,n2,n3,ekin,mod_arr1,a)
  use module_base
    implicit none
    integer, parameter :: lowfil=-14,lupfil=14
    integer, intent(in) :: n1,n2,n3
    integer, intent(in) :: mod_arr1(lowfil:n1+lupfil)
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: a
    integer ::ndat
    real(wp),intent(in):: x(0:n1,(n2+1)*(n3+1))
    real(wp),intent(inout)::y(0:n1,(n2+1)*(n3+1))
    real(gp),intent(inout)::ekin
    real(wp) tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
    integer ::i,l,i1,j
    ndat=(n2+1)*(n3+1) 
!$omp do
    do i=0,ndat/12-1
       do i1=0,n1
          tt1=0.e0_wp
          tt2=0.e0_wp
          tt3=0.e0_wp
          tt4=0.e0_wp
          tt5=0.e0_wp
          tt6=0.e0_wp
          tt7=0.e0_wp
          tt8=0.e0_wp
          tt9 =0.e0_wp
          tt10=0.e0_wp
          tt11=0.e0_wp
          tt12=0.e0_wp

          do l=lowfil,lupfil
             j=mod_arr1(i1+l)

             tt1=tt1+x(j,i*12+1)*a(l,1)
             tt2=tt2+x(j,i*12+2)*a(l,1)
             tt3=tt3+x(j,i*12+3)*a(l,1)
             tt4=tt4+x(j,i*12+4)*a(l,1)
             tt5=tt5+x(j,i*12+5)*a(l,1)
             tt6=tt6+x(j,i*12+6)*a(l,1)
             tt7=tt7+x(j,i*12+7)*a(l,1)
             tt8=tt8+x(j,i*12+8)*a(l,1)
             tt9 =tt9 +x(j,i*12+9 )*a(l,1)
             tt10=tt10+x(j,i*12+10)*a(l,1)
             tt11=tt11+x(j,i*12+11)*a(l,1)
             tt12=tt12+x(j,i*12+12)*a(l,1)
          enddo
          y(i1,i*12+1)=y(i1,i*12+1)+tt1;    ekin=ekin+tt1*x(i1,i*12+1)
          y(i1,i*12+2)=y(i1,i*12+2)+tt2;    ekin=ekin+tt2*x(i1,i*12+2)
          y(i1,i*12+3)=y(i1,i*12+3)+tt3;    ekin=ekin+tt3*x(i1,i*12+3)
          y(i1,i*12+4)=y(i1,i*12+4)+tt4;    ekin=ekin+tt4*x(i1,i*12+4)
          y(i1,i*12+5)=y(i1,i*12+5)+tt5;    ekin=ekin+tt5*x(i1,i*12+5)
          y(i1,i*12+6)=y(i1,i*12+6)+tt6;    ekin=ekin+tt6*x(i1,i*12+6)
          y(i1,i*12+7)=y(i1,i*12+7)+tt7;    ekin=ekin+tt7*x(i1,i*12+7)
          y(i1,i*12+8)=y(i1,i*12+8)+tt8;    ekin=ekin+tt8*x(i1,i*12+8)
          y(i1,i*12+9 )=y(i1,i*12+9 )+tt9 ;    ekin=ekin+tt9 *x(i1,i*12+9 )
          y(i1,i*12+10)=y(i1,i*12+10)+tt10;    ekin=ekin+tt10*x(i1,i*12+10)
          y(i1,i*12+11)=y(i1,i*12+11)+tt11;    ekin=ekin+tt11*x(i1,i*12+11)
          y(i1,i*12+12)=y(i1,i*12+12)+tt12;    ekin=ekin+tt12*x(i1,i*12+12)
       enddo
    enddo
!$omp enddo
!$omp do
    do i=(ndat/12)*12+1,ndat
       do i1=0,n1
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr1(i1+l)
             tt=tt+x(j   ,i)*a(l,1)
          enddo
          y(i1,i)=y(i1,i)+tt ; ekin=ekin+tt*x(i1,i)
       enddo
    enddo
!$omp enddo
  END SUBROUTINE conv_kin_x1
  
  subroutine conv_kin_y1(x,y,n1,n2,n3,ekin,mod_arr2,a)
  use module_base
    implicit none
    integer, parameter :: lowfil=-14,lupfil=14
    integer, intent(in) :: n1,n2,n3
    integer, intent(in) :: mod_arr2(lowfil:n1+lupfil)
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: a
    real(wp),intent(in):: x(0:n1,0:n2,0:n3)
    real(wp),intent(inout)::y(0:n1,0:n2,0:n3)
    real(gp),intent(inout)::ekin
    real(wp) tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt
    integer :: i3,i1,i2,l,j
!$omp do
    do i3=0,n3/8-1
       do i1=0,n1
          do i2=0,n2
             tt0=0.e0_wp
             tt1=0.e0_wp
             tt2=0.e0_wp
             tt3=0.e0_wp
             tt4=0.e0_wp
             tt5=0.e0_wp
             tt6=0.e0_wp
             tt7=0.e0_wp

             do l=lowfil,lupfil
                j=mod_arr2(i2+l)

                tt0=tt0+x(i1,j,i3*8+0)*a(l,2)
                tt1=tt1+x(i1,j,i3*8+1)*a(l,2)
                tt2=tt2+x(i1,j,i3*8+2)*a(l,2)
                tt3=tt3+x(i1,j,i3*8+3)*a(l,2)
                tt4=tt4+x(i1,j,i3*8+4)*a(l,2)
                tt5=tt5+x(i1,j,i3*8+5)*a(l,2)
                tt6=tt6+x(i1,j,i3*8+6)*a(l,2)
                tt7=tt7+x(i1,j,i3*8+7)*a(l,2)
             enddo
             y(i1,i2,i3*8+0)=y(i1,i2,i3*8+0)+tt0;    ekin=ekin+tt0*x(i1,i2,i3*8+0)
             y(i1,i2,i3*8+1)=y(i1,i2,i3*8+1)+tt1;    ekin=ekin+tt1*x(i1,i2,i3*8+1)
             y(i1,i2,i3*8+2)=y(i1,i2,i3*8+2)+tt2;    ekin=ekin+tt2*x(i1,i2,i3*8+2)
             y(i1,i2,i3*8+3)=y(i1,i2,i3*8+3)+tt3;    ekin=ekin+tt3*x(i1,i2,i3*8+3)
             y(i1,i2,i3*8+4)=y(i1,i2,i3*8+4)+tt4;    ekin=ekin+tt4*x(i1,i2,i3*8+4)
             y(i1,i2,i3*8+5)=y(i1,i2,i3*8+5)+tt5;    ekin=ekin+tt5*x(i1,i2,i3*8+5)
             y(i1,i2,i3*8+6)=y(i1,i2,i3*8+6)+tt6;    ekin=ekin+tt6*x(i1,i2,i3*8+6)
             y(i1,i2,i3*8+7)=y(i1,i2,i3*8+7)+tt7;    ekin=ekin+tt7*x(i1,i2,i3*8+7)
          enddo
       enddo
    enddo
!$omp enddo
!$omp do
    do i3=(n3/8)*8,n3
       do i1=0,n1
          do i2=0,n2
             tt=0.e0_wp
             do l=lowfil,lupfil
                j=mod_arr2(i2+l)
                tt=tt+x(i1,j   ,i3)*a(l,2)
             enddo
             y(i1,i2,i3)=y(i1,i2,i3)+tt;   ekin=ekin+tt*x(i1,i2,i3)
          enddo
       enddo
    enddo
!$omp enddo
  END SUBROUTINE conv_kin_y1


  subroutine conv_kin_z1(x,y,n1,n2,n3,ekin,mod_arr3,a)
  use module_base
    implicit none
    integer, parameter :: lowfil=-14,lupfil=14
    integer, intent(in) :: n1,n2,n3
    integer, intent(in) :: mod_arr3(lowfil:n1+lupfil)
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: a
    integer ::ndat
    real(wp),intent(in):: x((n1+1)*(n2+1),0:n1)
    real(wp),intent(inout)::y((n1+1)*(n2+1),0:n1)
    real(gp),intent(inout)::ekin
    real(wp) tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
    integer :: i,i3,l,j
    ndat=(n1+1)*(n2+1)
!$omp do
    do i=0,ndat/12-1
       do i3=0,n3
          tt1=0.e0_wp
          tt2=0.e0_wp
          tt3=0.e0_wp
          tt4=0.e0_wp
          tt5=0.e0_wp
          tt6=0.e0_wp
          tt7=0.e0_wp
          tt8=0.e0_wp
          tt9 =0.e0_wp
          tt10=0.e0_wp
          tt11=0.e0_wp
          tt12=0.e0_wp

          do l=lowfil,lupfil
             j=mod_arr3(i3+l)

             tt1=tt1+x(i*12+1,j)*a(l,3)
             tt2=tt2+x(i*12+2,j)*a(l,3)
             tt3=tt3+x(i*12+3,j)*a(l,3)
             tt4=tt4+x(i*12+4,j)*a(l,3)
             tt5=tt5+x(i*12+5,j)*a(l,3)
             tt6=tt6+x(i*12+6,j)*a(l,3)
             tt7=tt7+x(i*12+7,j)*a(l,3)
             tt8=tt8+x(i*12+8,j)*a(l,3)
             tt9 =tt9 +x(i*12+9 ,j)*a(l,3)
             tt10=tt10+x(i*12+10,j)*a(l,3)
             tt11=tt11+x(i*12+11,j)*a(l,3)
             tt12=tt12+x(i*12+12,j)*a(l,3)
          enddo

          y(i*12+1,i3)=y(i*12+1,i3)+tt1;    ekin=ekin+tt1*x(i*12+1,i3)
          y(i*12+2,i3)=y(i*12+2,i3)+tt2;    ekin=ekin+tt2*x(i*12+2,i3)
          y(i*12+3,i3)=y(i*12+3,i3)+tt3;    ekin=ekin+tt3*x(i*12+3,i3)
          y(i*12+4,i3)=y(i*12+4,i3)+tt4;    ekin=ekin+tt4*x(i*12+4,i3)
          y(i*12+5,i3)=y(i*12+5,i3)+tt5;    ekin=ekin+tt5*x(i*12+5,i3)
          y(i*12+6,i3)=y(i*12+6,i3)+tt6;    ekin=ekin+tt6*x(i*12+6,i3)
          y(i*12+7,i3)=y(i*12+7,i3)+tt7;    ekin=ekin+tt7*x(i*12+7,i3)
          y(i*12+8,i3)=y(i*12+8,i3)+tt8;    ekin=ekin+tt8*x(i*12+8,i3)
          y(i*12+9 ,i3)=y(i*12+9 ,i3)+tt9 ;    ekin=ekin+tt9*x(i*12+9 ,i3)
          y(i*12+10,i3)=y(i*12+10,i3)+tt10;    ekin=ekin+tt10*x(i*12+10,i3)
          y(i*12+11,i3)=y(i*12+11,i3)+tt11;    ekin=ekin+tt11*x(i*12+11,i3)
          y(i*12+12,i3)=y(i*12+12,i3)+tt12;    ekin=ekin+tt12*x(i*12+12,i3)
       enddo
    enddo
!$omp enddo
!$omp do
    do i=(ndat/12)*12+1,ndat
       do i3=0,n3
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr3(i3+l)
             tt=tt+x(i,j)*a(l,3)
          enddo
          y(i,i3)=y(i,i3)+tt; ekin=ekin+tt*x(i,i3)
       enddo
    enddo
!$omp enddo
  END SUBROUTINE conv_kin_z1


subroutine convolut_kinetic_hyb_c(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,x_c,x_f,y_c,y_f,cprecr,x_f1,x_f2,x_f3,ibyz,ibxz,ibxy)
  !   y = (kinetic energy operator+C)x 
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer,intent(in)::ibyz(2,0:n2,0:n3),ibxz(2,0:n1,0:n3),ibxy(2,0:n1,0:n2)
  real(gp), intent(in) :: hgrid(3)
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: x_f3
  real(gp), intent(in) :: cprecr
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  integer :: i,t,i1,i2,i3
  !integer :: ncount_max,ncount_rate,ncount1,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: l
  real(wp) :: dyi,t112,t121,t122,t212,t221,t222,t211
  real(wp)::scale(3)
  !real(kind=8) :: tel
  real(wp), dimension(lowfil:lupfil,3) :: a,b,c,e
  integer ii1,ii2,ii3 ! for hybrid convolutions
  !n(c) integer::ic
  integer :: mod_arr1(lowfil:n1+lupfil)
  integer :: mod_arr2(lowfil:n2+lupfil)
  integer :: mod_arr3(lowfil:n3+lupfil)
!$omp parallel default(private) shared(mod_arr1,mod_arr2,mod_arr3,n1,n2,n3,hgrid,x_c,x_f,x_f1,x_f2,x_f3) &
!$omp shared(cprecr,y_c,y_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,ibyz,ibxz,ibxy)
        !n(c) ic=1    

  call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
  call fill_mod_arr(mod_arr2,lowfil,n2+lupfil,n2+1)
  call fill_mod_arr(mod_arr3,lowfil,n3+lupfil,n3+1)

  scale(:)=-.5_wp/hgrid(:)**2
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0,:)=   -3.5536922899131901941296809374_wp*scale(:)
  a(1,:)=    2.2191465938911163898794546405_wp*scale(:)
  a(2,:)=   -0.6156141465570069496314853949_wp*scale(:)
  a(3,:)=    0.2371780582153805636239247476_wp*scale(:)
  a(4,:)=   -0.0822663999742123340987663521_wp*scale(:)
  a(5,:)=    0.02207029188482255523789911295638968409_wp*scale(:)
  a(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  a(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  a(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  a(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  a(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  a(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  a(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  a(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  a(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     a(-i,:)=a(i,:)
  enddo
  !  <phi|D^2|psi_i>
  c(-14,:)=     -3.869102413147656535541850057188e-18_wp*scale(:)
  c(-13,:)=      1.5130616560866154733900029272077362e-13_wp*scale(:)
  c(-12,:)=     -3.2264702314010525539061647271983988409e-11_wp*scale(:)
  c(-11,:)=     -5.96264938781402337319841002642e-9_wp*scale(:)
  c(-10,:)=     -2.1656830629214041470164889350342e-7_wp*scale(:)
  c(-9 ,:)=      8.7969704055286288323596890609625e-7_wp*scale(:)
  c(-8 ,:)=     -0.00001133456724516819987751818232711775_wp*scale(:)
  c(-7 ,:)=      0.00021710795484646138591610188464622454_wp*scale(:)
  c(-6 ,:)=     -0.0021356291838797986414312219042358542_wp*scale(:)
  c(-5 ,:)=      0.00713761218453631422925717625758502986_wp*scale(:)
  c(-4 ,:)=     -0.0284696165863973422636410524436931061_wp*scale(:)
  c(-3 ,:)=      0.14327329352510759457155821037742893841_wp*scale(:)
  c(-2 ,:)=     -0.42498050943780130143385739554118569733_wp*scale(:)
  c(-1 ,:)=      0.65703074007121357894896358254040272157_wp*scale(:)
  c( 0 ,:)=     -0.42081655293724308770919536332797729898_wp*scale(:)
  c( 1 ,:)=     -0.21716117505137104371463587747283267899_wp*scale(:)
  c( 2 ,:)=      0.63457035267892488185929915286969303251_wp*scale(:)
  c( 3 ,:)=     -0.53298223962800395684936080758073568406_wp*scale(:)
  c( 4 ,:)=      0.23370490631751294307619384973520033236_wp*scale(:)
  c( 5 ,:)=     -0.05657736973328755112051544344507997075_wp*scale(:)
  c( 6 ,:)=      0.0080872029411844780634067667008050127_wp*scale(:)
  c( 7 ,:)=     -0.00093423623304808664741804536808932984_wp*scale(:)
  c( 8 ,:)=      0.00005075807947289728306309081261461095_wp*scale(:)
  c( 9 ,:)=     -4.62561497463184262755416490048242e-6_wp*scale(:)
  c( 10,:)=      6.3919128513793415587294752371778e-7_wp*scale(:)
  c( 11,:)=      1.87909235155149902916133888931e-8_wp*scale(:)
  c( 12,:)=      1.04757345962781829480207861447155543883e-10_wp*scale(:)
  c( 13,:)=     -4.84665690596158959648731537084025836e-13_wp*scale(:)
  c( 14,:)=      1.2392629629188986192855777620877e-17_wp*scale(:)
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i,:)=c(-i,:)
  enddo
  !<psi|D^2|psi_i>
  e(0,:)=   -24.875846029392331358907766562_wp*scale(:)
  e(1,:)=   -7.1440597663471719869313377994_wp*scale(:)
  e(2,:)=   -0.04251705323669172315864542163525830944_wp*scale(:)
  e(3,:)=   -0.26995931336279126953587091167128839196_wp*scale(:)
  e(4,:)=    0.08207454169225172612513390763444496516_wp*scale(:)
  e(5,:)=   -0.02207327034586634477996701627614752761_wp*scale(:)
  e(6,:)=    0.00409765642831595181639002667514310145_wp*scale(:)
  e(7,:)=   -0.00045167920287507774929432548999880117_wp*scale(:)
  e(8,:)=    0.00002398228524507599670405555359023135_wp*scale(:)
  e(9,:)=   -2.0904234952920365957922889447361e-6_wp*scale(:)
  e(10,:)=   3.7230763047369275848791496973044e-7_wp*scale(:)
  e(11,:)=   1.05857055496741470373494132287e-8_wp*scale(:)
  e(12,:)=   5.8138798302825405479592506674648873655e-11_wp*scale(:)
  e(13,:)=  -2.70800493626319438269856689037647576e-13_wp*scale(:)
  e(14,:)=   6.924474940639200152025730585882e-18_wp*scale(:)
  do i=1,14
     e(-i,:)=e(i,:)
  enddo

  !---------------------------------------------------------------------------

  ! Scaling function part: 
  ! periodic convolutions for scaling functions, 
  ! free convolutions for wavelets


  call conv_kin_x2(x_c,y_c,n1,n2,n3,mod_arr1,a,cprecr)   
  call conv_kin_y2(x_c,y_c,n1,n2,n3,mod_arr2,a)   
  call conv_kin_z2(x_c,y_c,n1,n2,n3,mod_arr3,a)  
  
!  ! (1/2) d^2/dx^2
!  scaling function-scaling function 
!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        ! wavelet-scaling function:
        ! scaling function everywhere, also outside of the box
        ! wavelet restricted
        do i1=ibyz(1,i2,i3)-lupfil,ibyz(2,i2,i3)-lowfil
           ii1=mod_arr1(i1)
           dyi=0._wp
           do l=max(lowfil,ibyz(1,i2,i3)-i1),min(lupfil,ibyz(2,i2,i3)-i1)
              dyi=dyi + x_f1(i1+l,i2,i3)*b(l,1)
           enddo
           y_c(ii1,i2,i3)=y_c(ii1,i2,i3)+dyi
        enddo
        
        ! scaling function-wavelet
        ! wavelet restricted, scaling function everywhere
        do i1=ibyz(1,i2,i3),ibyz(2,i2,i3)
           t211=x_f(1,i1,i2,i3)*cprecr
           do l=lowfil,lupfil
              t=mod_arr1(i1+l)
              t211=t211 + x_c(t,i2,i3)*c(l,1)
           enddo
           y_f(1,i1,i2,i3)=t211
        enddo
     enddo
  enddo
!$omp enddo  
!!
!!  ! + (1/2) d^2/dy^2
!!
!$omp do  
   do i3=nfl3,nfu3
      do i1=nfl1,nfu1
         do i2=ibxz(1,i1,i3)-lupfil,ibxz(2,i1,i3)-lowfil
            ii2=mod_arr2(i2)
            dyi=0._wp
            do l=max(lowfil,ibxz(1,i1,i3)-i2),min(lupfil,ibxz(2,i1,i3)-i2)
               dyi=dyi + x_f2(i2+l,i1,i3)*b(l,2)
            enddo
            y_c(i1,ii2,i3)=y_c(i1,ii2,i3)+dyi
         enddo
         
         do i2=ibxz(1,i1,i3),ibxz(2,i1,i3)
            t121=x_f(2,i1,i2,i3)*cprecr 
            do l=lowfil,lupfil
               t=mod_arr2(i2+l)   
               t121=t121 + x_c(i1,t,i3)*c(l,2)
            enddo
            y_f(2,i1,i2,i3)=t121
         enddo
      enddo
   enddo
!$omp enddo  
  
!!
!!  ! + (1/2) d^2/dz^2
!!
!$omp do 
   do i2=nfl2,nfu2
      do i1=nfl1,nfu1
         do i3=ibxy(1,i1,i2)-lupfil,ibxy(2,i1,i2)-lowfil
            dyi=0._wp
            ii3=mod_arr3(i3)
            do l=max(lowfil,ibxy(1,i1,i2)-i3),min(lupfil,ibxy(2,i1,i2)-i3)
               dyi=dyi + x_f3(i3+l,i1,i2)*b(l,3)
            enddo
            y_c(i1,i2,ii3)=y_c(i1,i2,ii3)+dyi
         enddo

         do i3=ibxy(1,i1,i2),ibxy(2,i1,i2)
            t112=x_f(4,i1,i2,i3)*cprecr 
            do l=lowfil,lupfil
               t=mod_arr3(i3+l)   
               t112=t112 + x_c(i1,i2,t)*c(l,3)
            enddo
            y_f(4,i1,i2,i3)=t112
         enddo
      enddo
   enddo
!$omp enddo   
  
!  write(20,*) ibyz
!  stop
  ! wavelet part
  ! completely similar to the free case
  ! (1/2) d^2/dx^2
!$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz(1,i2,i3),ibyz(2,i2,i3)
           t112=0._wp
           t121=0._wp
           t122=x_f(6,i1,i2,i3)*cprecr
           t212=x_f(5,i1,i2,i3)*cprecr
           t221=x_f(3,i1,i2,i3)*cprecr
           t222=x_f(7,i1,i2,i3)*cprecr
           t211=0._wp 
           do l=max(ibyz(1,i2,i3)-i1,lowfil),min(lupfil,ibyz(2,i2,i3)-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l,1) + x_f(5,i1+l,i2,i3)*b(l,1)
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l,1) + x_f(3,i1+l,i2,i3)*b(l,1)
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l,1) + x_f(7,i1+l,i2,i3)*b(l,1)
              t212=t212 + x_f(4,i1+l,i2,i3)*c(l,1) + x_f(5,i1+l,i2,i3)*e(l,1)
              t221=t221 + x_f(2,i1+l,i2,i3)*c(l,1) + x_f(3,i1+l,i2,i3)*e(l,1)
              t222=t222 + x_f(6,i1+l,i2,i3)*c(l,1) + x_f(7,i1+l,i2,i3)*e(l,1)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l,1)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=t122
           y_f(5,i1,i2,i3)=t212
           y_f(3,i1,i2,i3)=t221
           y_f(7,i1,i2,i3)=t222
        enddo
     enddo
  enddo
!$omp enddo
  ! + (1/2) d^2/dy^2
!$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz(1,i1,i3),ibxz(2,i1,i3)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(ibxz(1,i1,i3)-i2,lowfil),min(lupfil,ibxz(2,i1,i3)-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l,2) + x_f(6,i1,i2+l,i3)*b(l,2)
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l,2) + x_f(3,i1,i2+l,i3)*b(l,2)
              t122=t122 + x_f(4,i1,i2+l,i3)*c(l,2) + x_f(6,i1,i2+l,i3)*e(l,2)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l,2) + x_f(7,i1,i2+l,i3)*b(l,2)
              t221=t221 + x_f(1,i1,i2+l,i3)*c(l,2) + x_f(3,i1,i2+l,i3)*e(l,2)
              t222=t222 + x_f(5,i1,i2+l,i3)*c(l,2) + x_f(7,i1,i2+l,i3)*e(l,2)
              t121=t121 + x_f(2,i1,i2+l,i3)*e(l,2)
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
  ! + (1/2) d^2/dz^2
!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy(1,i1,i2),ibxy(2,i1,i2)
           t112=0._wp;t121=0._wp;t122=0._wp;t212=0._wp;t221=0._wp;t222=0._wp;t211=0._wp 
           do l=max(ibxy(1,i1,i2)-i3,lowfil),min(lupfil,ibxy(2,i1,i2)-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l,3) + x_f(6,i1,i2,i3+l)*b(l,3)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l,3) + x_f(5,i1,i2,i3+l)*b(l,3)
              t122=t122 + x_f(2,i1,i2,i3+l)*c(l,3) + x_f(6,i1,i2,i3+l)*e(l,3)
              t212=t212 + x_f(1,i1,i2,i3+l)*c(l,3) + x_f(5,i1,i2,i3+l)*e(l,3)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l,3) + x_f(7,i1,i2,i3+l)*b(l,3)
              t222=t222 + x_f(3,i1,i2,i3+l)*c(l,3) + x_f(7,i1,i2,i3+l)*e(l,3)
              t112=t112 + x_f(4,i1,i2,i3+l)*e(l,3)
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
END SUBROUTINE convolut_kinetic_hyb_c

  subroutine conv_kin_y2(x,y,n1,n2,n3,mod_arr2,a)
  use module_base
    implicit none
    integer, parameter :: lowfil=-14,lupfil=14
    integer, intent(in) :: n1,n2,n3
    integer, intent(in) :: mod_arr2(lowfil:n1+lupfil)
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: a
    real(wp),intent(in):: x(0:n1,0:n2,0:n3)
    real(wp),intent(inout)::y(0:n1,0:n2,0:n3)
    real(wp) tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt
    integer :: i3,i1,i2,l,j
!$omp do
    do i3=0,n3/8-1
       do i1=0,n1
          do i2=0,n2
             tt0=0.e0_wp
             tt1=0.e0_wp
             tt2=0.e0_wp
             tt3=0.e0_wp
             tt4=0.e0_wp
             tt5=0.e0_wp
             tt6=0.e0_wp
             tt7=0.e0_wp
             do l=lowfil,lupfil
                j=mod_arr2(i2+l)

                tt0=tt0+x(i1,j,i3*8+0)*a(l,2)
                tt1=tt1+x(i1,j,i3*8+1)*a(l,2)
                tt2=tt2+x(i1,j,i3*8+2)*a(l,2)
                tt3=tt3+x(i1,j,i3*8+3)*a(l,2)
                tt4=tt4+x(i1,j,i3*8+4)*a(l,2)
                tt5=tt5+x(i1,j,i3*8+5)*a(l,2)
                tt6=tt6+x(i1,j,i3*8+6)*a(l,2)
                tt7=tt7+x(i1,j,i3*8+7)*a(l,2)
             enddo
             y(i1,i2,i3*8+0)=y(i1,i2,i3*8+0)+tt0;        
             y(i1,i2,i3*8+1)=y(i1,i2,i3*8+1)+tt1;        
             y(i1,i2,i3*8+2)=y(i1,i2,i3*8+2)+tt2;        
             y(i1,i2,i3*8+3)=y(i1,i2,i3*8+3)+tt3;        
             y(i1,i2,i3*8+4)=y(i1,i2,i3*8+4)+tt4;        
             y(i1,i2,i3*8+5)=y(i1,i2,i3*8+5)+tt5;        
             y(i1,i2,i3*8+6)=y(i1,i2,i3*8+6)+tt6;        
             y(i1,i2,i3*8+7)=y(i1,i2,i3*8+7)+tt7;        
          enddo
       enddo
    enddo
!$omp enddo
!$omp do
    do i3=(n3/8)*8,n3
       do i1=0,n1
          do i2=0,n2
             tt=0.e0_wp
             do l=lowfil,lupfil
                j=mod_arr2(i2+l)
                tt=tt+x(i1,j   ,i3)*a(l,2)
             enddo
             y(i1,i2,i3)=y(i1,i2,i3)+tt;        
          enddo
       enddo
    enddo
!$omp enddo
  END SUBROUTINE conv_kin_y2


  subroutine conv_kin_x2(x,y,n1,n2,n3,mod_arr1,a,cprecr)
  use module_base
    implicit none
    integer, parameter :: lowfil=-14,lupfil=14   
    integer, intent(in) :: n1,n2,n3
    integer, intent(in) :: mod_arr1(lowfil:n1+lupfil)
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: a
    integer ::ndat
    real(gp), intent(in) :: cprecr
    real(wp),intent(in):: x(0:n1,(n2+1)*(n3+1))
    real(wp),intent(inout)::y(0:n1,(n2+1)*(n3+1))
    real(wp) tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
    integer ::i,l,i1,j
    ndat=(n2+1)*(n3+1)
!$omp do
    do i=0,ndat/12-1
       do i1=0,n1
          tt1 =x(i1,i*12+1)*cprecr
          tt2 =x(i1,i*12+2)*cprecr
          tt3 =x(i1,i*12+3)*cprecr
          tt4 =x(i1,i*12+4)*cprecr
          tt5 =x(i1,i*12+5)*cprecr
          tt6 =x(i1,i*12+6)*cprecr
          tt7 =x(i1,i*12+7)*cprecr
          tt8 =x(i1,i*12+8)*cprecr
          tt9 =x(i1,i*12+9 )*cprecr
          tt10=x(i1,i*12+10)*cprecr
          tt11=x(i1,i*12+11)*cprecr
          tt12=x(i1,i*12+12)*cprecr

          do l=lowfil,lupfil
             j=mod_arr1(i1+l)

             tt1=tt1+x(j,i*12+1)*a(l,1)
             tt2=tt2+x(j,i*12+2)*a(l,1)
             tt3=tt3+x(j,i*12+3)*a(l,1)
             tt4=tt4+x(j,i*12+4)*a(l,1)
             tt5=tt5+x(j,i*12+5)*a(l,1)
             tt6=tt6+x(j,i*12+6)*a(l,1)
             tt7=tt7+x(j,i*12+7)*a(l,1)
             tt8=tt8+x(j,i*12+8)*a(l,1)
             tt9 =tt9 +x(j,i*12+9 )*a(l,1)
             tt10=tt10+x(j,i*12+10)*a(l,1)
             tt11=tt11+x(j,i*12+11)*a(l,1)
             tt12=tt12+x(j,i*12+12)*a(l,1)
          enddo
          y(i1,i*12+1)=tt1
          y(i1,i*12+2)=tt2
          y(i1,i*12+3)=tt3
          y(i1,i*12+4)=tt4
          y(i1,i*12+5)=tt5
          y(i1,i*12+6)=tt6
          y(i1,i*12+7)=tt7
          y(i1,i*12+8)=tt8
          y(i1,i*12+9 )=tt9 
          y(i1,i*12+10)=tt10
          y(i1,i*12+11)=tt11
          y(i1,i*12+12)=tt12
       enddo
    enddo
!$omp enddo
!$omp do
    do i=(ndat/12)*12+1,ndat
       do i1=0,n1
          tt=x(i1,i)*cprecr
          do l=lowfil,lupfil
             j=mod_arr1(i1+l)
             tt=tt+x(j   ,i)*a(l,1)
          enddo
          y(i1,i)=tt
       enddo
    enddo
!$omp enddo
  END SUBROUTINE conv_kin_x2

  subroutine conv_kin_z2(x,y,n1,n2,n3,mod_arr3,a)
  use module_base
    implicit none
    integer, parameter :: lowfil=-14,lupfil=14
    integer, intent(in) :: n1,n2,n3
    integer, intent(in) :: mod_arr3(lowfil:n1+lupfil)
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: a
    integer ::ndat
    real(wp),intent(in):: x((n1+1)*(n2+1),0:n1)
    real(wp),intent(inout)::y((n1+1)*(n2+1),0:n1)
    real(wp) tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
    integer :: i,i3,l,j
    ndat=(n1+1)*(n2+1)
!$omp do
    do i=0,ndat/12-1
       do i3=0,n3
          tt1=0.e0_wp
          tt2=0.e0_wp
          tt3=0.e0_wp
          tt4=0.e0_wp
          tt5=0.e0_wp
          tt6=0.e0_wp
          tt7=0.e0_wp
          tt8=0.e0_wp
          tt9 =0.e0_wp
          tt10=0.e0_wp
          tt11=0.e0_wp
          tt12=0.e0_wp

          do l=lowfil,lupfil
             j=mod_arr3(i3+l)

             tt1=tt1+x(i*12+1,j)*a(l,3)
             tt2=tt2+x(i*12+2,j)*a(l,3)
             tt3=tt3+x(i*12+3,j)*a(l,3)
             tt4=tt4+x(i*12+4,j)*a(l,3)
             tt5=tt5+x(i*12+5,j)*a(l,3)
             tt6=tt6+x(i*12+6,j)*a(l,3)
             tt7=tt7+x(i*12+7,j)*a(l,3)
             tt8=tt8+x(i*12+8,j)*a(l,3)
             tt9 =tt9 +x(i*12+9 ,j)*a(l,3)
             tt10=tt10+x(i*12+10,j)*a(l,3)
             tt11=tt11+x(i*12+11,j)*a(l,3)
             tt12=tt12+x(i*12+12,j)*a(l,3)
          enddo

          y(i*12+1,i3)=y(i*12+1,i3)+tt1;       
          y(i*12+2,i3)=y(i*12+2,i3)+tt2;       
          y(i*12+3,i3)=y(i*12+3,i3)+tt3;       
          y(i*12+4,i3)=y(i*12+4,i3)+tt4;       
          y(i*12+5,i3)=y(i*12+5,i3)+tt5;       
          y(i*12+6,i3)=y(i*12+6,i3)+tt6;       
          y(i*12+7,i3)=y(i*12+7,i3)+tt7;       
          y(i*12+8,i3)=y(i*12+8,i3)+tt8;       
          y(i*12+9 ,i3)=y(i*12+9 ,i3)+tt9 ;    
          y(i*12+10,i3)=y(i*12+10,i3)+tt10;    
          y(i*12+11,i3)=y(i*12+11,i3)+tt11;    
          y(i*12+12,i3)=y(i*12+12,i3)+tt12;    
       enddo
    enddo
!$omp enddo
!$omp do
    do i=(ndat/12)*12+1,ndat
       do i3=0,n3
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr3(i3+l)
             tt=tt+x(i,j)*a(l,3)
          enddo
          y(i,i3)=y(i,i3)+tt; 
       enddo
    enddo
!$omp enddo
END SUBROUTINE conv_kin_z2
