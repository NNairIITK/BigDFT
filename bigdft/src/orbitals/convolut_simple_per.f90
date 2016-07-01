!> @file
!! Simple routines of convolutions
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine ana_rot_per(n,ndat,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
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
           k=modulo(l+2*i,2*n+2)
            ci=ci+ch(l)*x(k    ,j)
            di=di+cg(l)*x(k    ,j)
        enddo
        y(j,i)=ci
        y(j,n+1+i)=di
     enddo

  enddo
END SUBROUTINE ana_rot_per


subroutine syn_rot_per(n,ndat,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
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

     do i=0,n
        se=0.e0_wp
        so=0.e0_wp
        do l=-4,4
           k=modulo(i-l,n+1)
           se=se+ch(2*l  )*x(  k,j)+cg(2*l  )*x(n+1+k  ,j)
           so=so+ch(2*l+1)*x(  k,j)+cg(2*l+1)*x(n+1+k  ,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

  enddo

  !first changement
  do i=0,n-1
     do l=-4,3
        so=so+ch(2*l+1)*x(i-1-l,j)+cg(2*l+1)*x(n+i-1-l  ,j)
        se=se+ch(2*(l+1))*x(i-1-l,j)+cg(2*(l+1))*x(n+i-1-l,j)
     end do
     y(j,2*(i-1)+1)=so
     y(j,2*i  )=se
  end do
 

  !second changement
  i=0
     do l=-4,3
        so=so+ch(2*l+1)*x(i-1-l,j)+cg(2*l+1)*x(n+i-1-l  ,j)
        se=se+ch(2*(l+1))*x(i-1-l,j)+cg(2*(l+1))*x(n+i-1-l,j)
     end do
  y(j,2*n-1)=so
  y(j,0  )=se

  do i=1,n-1
     do l=-4,3
        so=so+ch(2*l+1)*x(i-1-l,j)+cg(2*l+1)*x(n+i-1-l  ,j)
        se=se+ch(2*(l+1))*x(i-1-l,j)+cg(2*(l+1))*x(n+i-1-l,j)
     end do
     y(j,2*(i-1)+1)=so
     y(j,2*i  )=se
  end do

  !third changement
  i=0
     do l=-4,3
        so=so+ch(2*l+1)*x(i-1-l,j)+cg(2*l+1)*x(n+i-1-l  ,j)
        se=se+ch(2*(l+1))*x(i-1-l,j)+cg(2*(l+1))*x(n+i-1-l,j)
     end do
  y(j,2*n-1)=so
  y(j,0  )=se

  do i=0,n-2
     do l=-4,3
        so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+i-l  ,j)
        se=se+ch(2*(l+1))*x(i-l,j)+cg(2*(l+1))*x(n+i-l,j)
     end do
     y(j,2*i+1)=so
     y(j,2*i+2)=se
  end do

  !filter order
  do l=-3,4
     so=so+ch(-2*l+1)*x(i-1+l,j)+cg(-2*l+1)*x(n+i-1+l  ,j)
     se=se+ch(2*(-l+1))*x(i-1+l,j)+cg(2*(-l+1))*x(n+i-1+l,j)
  end do


  !fourth changement
  i=0
     do l=-3,4
        so=so+ch(-2*l+1)*x(i-1+l,j)+cg(-2*l+1)*x(n+i-1+l  ,j)
        se=se+ch(2*(-l+1))*x(i-1+l,j)+cg(2*(-l+1))*x(n+i-1+l,j)
     end do
  y(j,2*n-1)=so
  y(j,0  )=se

  do i=0,n-2
     do l=-3,4
        so=so+ch(-2*l+1)*x(i+l,j)+cg(-2*l+1)*x(n+i+l  ,j)
        se=se+ch(2*(-l+1))*x(i+l,j)+cg(2*(-l+1))*x(n+i+l,j)
     end do
     y(j,2*i+1)=so
     y(j,2*i+2)=se
  end do

  !fifth changement
  i=0
     do l=-3,4
        k=modulo(i-1+l,n)
        so=so+ch(-2*l+1)*x(k,j)+cg(-2*l+1)*x(n+k  ,j)
        se=se+ch(2*(-l+1))*x(k,j)+cg(2*(-l+1))*x(n+k,j)
     end do
  y(j,2*n-1)=so
  y(j,0  )=se

  do i=0,n-2
     do l=-3,4
        k=modulo(i+l,n)
        so=so+ch(-2*l+1)*x(k,j)+cg(-2*l+1)*x(n+k ,j)
        se=se+ch(2*(-l+1))*x(k,j)+cg(2*(-l+1))*x(n+k,j)
     end do
     y(j,2*i+1)=so
     y(j,2*i+2)=se
  end do


END SUBROUTINE syn_rot_per

subroutine convrot_n_per(n1,ndat,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l
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
     do i=0,n1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1+1)   
           tt=tt+x(  k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE convrot_n_per


subroutine convrot_t_per(n1,ndat,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-7,lupfil=8
  integer :: i,j,k,l
  real(wp) :: tt
  ! the filtered output data structure has shrunk by the filter length

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
           k=modulo(i+l,n1+1)
           tt=tt+x(k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE convrot_t_per


subroutine convolut_kinetic_per_c(n1,n2,n3,hgrid,x,y,c)
!   applies the kinetic energy operator onto x to get y. Works for periodic BC
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  
  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt=x(i1,i2,i3)*c
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt=tt+x(j   ,i2,i3)*fil(l,1)
           enddo
           y(i1,i2,i3)=tt
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i2+l,n2+1)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
     
  enddo

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt=tt+x(i1,i2,   j)*fil(l,3)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
  enddo
  
END SUBROUTINE convolut_kinetic_per_c


subroutine convolut_kinetic_per_T(n1,n2,n3,hgrid,x,y,ekin)
!   applies the kinetic energy operator onto x to get y. Works for periodic BC
  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  real(wp),intent(out)::ekin
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  ekin=0.0_wp

  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt=0.0_wp
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt=tt+x(j   ,i2,i3)*fil(l,1)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
           ekin=ekin+x(i1,i2,i3)*tt
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt=0.0_wp
           do l=lowfil,lupfil
              j=modulo(i2+l,n2+1)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
           ekin=ekin+x(i1,i2,i3)*tt
        enddo
     enddo
     
  enddo

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt=0.0_wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt=tt+x(i1,i2,   j)*fil(l,3)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
           ekin=ekin+x(i1,i2,i3)*tt
        enddo
     enddo
  enddo
  
END SUBROUTINE convolut_kinetic_per_T
