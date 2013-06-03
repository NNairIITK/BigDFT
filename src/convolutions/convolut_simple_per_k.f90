!> @file
!!  Simple routine of convolutions
!! @author
!!    Copyright (C) 2005-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>   Applies the modified kinetic energy operator onto x to get y. 
!!   Works for periodic BC.
!!   Modified kinetic energy operator:
!!   A=-1/2 exp(Ikr) Delta exp(-Ikr)+C
!!   where k=(k1,k2,k3); r=(x,y,z)
!!   We apply the Galerkin matrix of this operator
!!   in the scaling function basis phi(x/h1-i1)phi(y/h2-i2)phi(z/h3-i3)
!!   The only difference from the case k=0 is the first derivative operator
!!   2I( kx d/dx + ...+ kz d/dz)
!!   for which we need the 1-dimensional Galerkin matrix
!!   @f$ 1/h\int phi(x/h) d/dx phi(x/h-i) dx @f$ 
!!   that is stored in the arrays fil(2,:,1)..fil(2,:,3)
!!   multiplied by the factor scale1.
!!   The second derivative operator is stored in the arrays
!!   fil(1,:,1)..fil(1,..3) in the usual way.
!!   The whole array fil actually stores the full complex Galerkin matrix 
!!   of the operator A.
!!
!!   One can check (I did that) that the Galerkin matrix (almost) annihilates the array
!!   @f$ c_j=\int exp(Ikx) phi(x/h-j)dx @f$ that is the wavelet projection
!!   of the plane wave exp(Ikx), if the plane wave is periodic in the simulation
!!   box. 
!!   The annihilation accuracy improves with decreasing grid constant h, as expected,
!!   but the numerical noise spoils the convergence after some point.
subroutine convolut_kinetic_per_c_k(n1,n2,n3,hgrid,x,y,c_in,k1,k2,k3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c_in,k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,l,j1,j2
  real(wp) :: tt1,tt2,tt3,tt4,c
  real(wp) :: tt1_2,tt2_2,tt3_2,tt4_2
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,0:lupfil,3) :: fil  

  scale (:)=real(-.5_gp/hgrid(:)**2,wp)
  scale1(1)=real(k1/hgrid(1),wp)
  scale1(2)=real(k2/hgrid(2),wp)
  scale1(3)=real(k3/hgrid(3),wp)
  c=c_in+.5_wp*(k1*k1+k2*k2+k3*k3)

  ! second derivative filters for Daubechies 16
  fil(1,0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(1,2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(1,3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(1,4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(1,5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(1,6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(1,7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(1,8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(1,9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(1,10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(1,11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(1,12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(1,13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(1,14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  fil(2,0,:)= 0._wp
  fil(2,1,:)= 0.8834460460908270942785856e0_wp*scale1(:)
  fil(2,2,:)= -0.3032593514765938346887962e0_wp*scale1(:)
  fil(2,3,:)= 0.1063640682894442760934532e0_wp*scale1(:)
  fil(2,4,:)= -0.03129014783948023634381564e0_wp*scale1(:)
  fil(2,5,:)= 0.006958379116450707495020408e0_wp*scale1(:)
  fil(2,6,:)= -0.001031530213375445369097965e0_wp*scale1(:)
  fil(2,7,:)= 0.00007667706908380351933901775e0_wp*scale1(:)
  fil(2,8,:)= 2.451992111053665419191564e-7_wp*scale1(:)
  fil(2,9,:)= 3.993810456408053712133667e-8_wp*scale1(:)
  fil(2,10,:)=-7.207948238588481597101904e-8_wp*scale1(:)
  fil(2,11,:)=-9.697184925637300947553069e-10_wp*scale1(:)
  fil(2,12,:)=-7.252206916665149851135592e-13_wp*scale1(:)
  fil(2,13,:)=1.240078536096648534547439e-14_wp*scale1(:)
  fil(2,14,:)=-1.585464751677102510097179e-19_wp*scale1(:)


!$omp parallel default (private) shared(n1,n2,n3,x,y,c,fil)
!$omp do

  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,lupfil-1
           tt1=x(1,i1,i2,i3)*c
           tt2=x(2,i1,i2,i3)*c
           tt1=tt1+x(1,i1,i2,i3)*fil(1,0,1)
           tt2=tt2+x(2,i1,i2,i3)*fil(1,0,1)
           tt3=0
           tt4=0
           do l=1,14
              j1=i1+l
              j2=modulo(i1-l,n1+1)
              tt1=tt1+(x(1,j1,i2,i3)+x(1,j2,i2,i3))*fil(1,l,1)
              tt3=tt3+(x(2,j1,i2,i3)-x(2,j2,i2,i3))*fil(2,l,1)
              tt2=tt2+(x(2,j1,i2,i3)+x(2,j2,i2,i3))*fil(1,l,1)
              tt4=tt4+(x(1,j1,i2,i3)-x(1,j2,i2,i3))*fil(2,l,1)
           enddo
           y(1,i1,i2,i3)=tt1-tt3
           y(2,i1,i2,i3)=tt2+tt4
        enddo
        do i1=lupfil,n1-lupfil
           tt1=x(1,i1,i2,i3)*c
           tt2=x(2,i1,i2,i3)*c
           tt1=tt1+x(1,i1,i2,i3)*fil(1,0,1)
           tt2=tt2+x(2,i1,i2,i3)*fil(1,0,1)
           tt3=0
           tt4=0
           do l=1,14
              j1=i1+l
              j2=i1-l
              tt1=tt1+(x(1,j1,i2,i3)+x(1,j2,i2,i3))*fil(1,l,1)
              tt3=tt3+(x(2,j1,i2,i3)-x(2,j2,i2,i3))*fil(2,l,1)
              tt2=tt2+(x(2,j1,i2,i3)+x(2,j2,i2,i3))*fil(1,l,1)
              tt4=tt4+(x(1,j1,i2,i3)-x(1,j2,i2,i3))*fil(2,l,1)
           enddo
           y(1,i1,i2,i3)=tt1-tt3
           y(2,i1,i2,i3)=tt2+tt4
        enddo
        do i1=n1-lupfil+1,n1
           tt1=x(1,i1,i2,i3)*c
           tt2=x(2,i1,i2,i3)*c
           tt1=tt1+x(1,i1,i2,i3)*fil(1,0,1)
           tt2=tt2+x(2,i1,i2,i3)*fil(1,0,1)
           tt3=0
           tt4=0
           do l=1,14
              j1=modulo(i1+l,n1+1)
              j2=i1-l
              tt1=tt1+(x(1,j1,i2,i3)+x(1,j2,i2,i3))*fil(1,l,1)
              tt3=tt3+(x(2,j1,i2,i3)-x(2,j2,i2,i3))*fil(2,l,1)
              tt2=tt2+(x(2,j1,i2,i3)+x(2,j2,i2,i3))*fil(1,l,1)
              tt4=tt4+(x(1,j1,i2,i3)-x(1,j2,i2,i3))*fil(2,l,1)
           enddo
           y(1,i1,i2,i3)=tt1-tt3
           y(2,i1,i2,i3)=tt2+tt4
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1-1,2
        do i2=0,lupfil-1
           tt1=x(1,i1,i2,i3)*fil(1,0,2)
           tt1_2=x(1,i1+1,i2,i3)*fil(1,0,2)
           tt2=x(2,i1,i2,i3)*fil(1,0,2)
           tt2_2=x(2,i1+1,i2,i3)*fil(1,0,2)
           tt3=0
           tt3_2=0
           tt4=0
           tt4_2=0
           do l=1,lupfil
              j1=i2+l
              j2=modulo(i2-l,n2+1)
              tt1=tt1+(x(1,i1,j1,i3)+x(1,i1,j2,i3))*fil(1,l,2)
              tt1_2=tt1_2+(x(1,i1+1,j1,i3)+x(1,i1+1,j2,i3))*fil(1,l,2)
              tt3=tt3+(x(2,i1,j1,i3)-x(2,i1,j2,i3))*fil(2,l,2)
              tt3_2=tt3_2+(x(2,i1+1,j1,i3)-x(2,i1+1,j2,i3))*fil(2,l,2)
              tt2=tt2+(x(2,i1,j1,i3)+x(2,i1,j2,i3))*fil(1,l,2)
              tt2_2=tt2_2+(x(2,i1+1,j1,i3)+x(2,i1+1,j2,i3))*fil(1,l,2)
              tt4=tt4+(x(1,i1,j1,i3)-x(1,i1,j2,i3))*fil(2,l,2)
              tt4_2=tt4_2+(x(1,i1+1,j1,i3)-x(1,i1+1,j2,i3))*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(1,i1+1,i2,i3)=y(1,i1+1,i2,i3)+tt1_2-tt3_2
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
           y(2,i1+1,i2,i3)=y(2,i1+1,i2,i3)+tt2_2+tt4_2
        enddo
        do i2=lupfil,n2-lupfil
           tt1=x(1,i1,i2,i3)*fil(1,0,2)
           tt1_2=x(1,i1+1,i2,i3)*fil(1,0,2)
           tt2=x(2,i1,i2,i3)*fil(1,0,2)
           tt2_2=x(2,i1+1,i2,i3)*fil(1,0,2)
           tt3=0
           tt3_2=0
           tt4=0
           tt4_2=0
           do l=1,lupfil
              j1=i2+l
              j2=i2-l
              tt1=tt1+(x(1,i1,j1,i3)+x(1,i1,j2,i3))*fil(1,l,2)
              tt1_2=tt1_2+(x(1,i1+1,j1,i3)+x(1,i1+1,j2,i3))*fil(1,l,2)
              tt3=tt3+(x(2,i1,j1,i3)-x(2,i1,j2,i3))*fil(2,l,2)
              tt3_2=tt3_2+(x(2,i1+1,j1,i3)-x(2,i1+1,j2,i3))*fil(2,l,2)
              tt2=tt2+(x(2,i1,j1,i3)+x(2,i1,j2,i3))*fil(1,l,2)
              tt2_2=tt2_2+(x(2,i1+1,j1,i3)+x(2,i1+1,j2,i3))*fil(1,l,2)
              tt4=tt4+(x(1,i1,j1,i3)-x(1,i1,j2,i3))*fil(2,l,2)
              tt4_2=tt4_2+(x(1,i1+1,j1,i3)-x(1,i1+1,j2,i3))*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(1,i1+1,i2,i3)=y(1,i1+1,i2,i3)+tt1_2-tt3_2
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
           y(2,i1+1,i2,i3)=y(2,i1+1,i2,i3)+tt2_2+tt4_2
        enddo
        do i2=n2-lupfil+1,n2
           tt1=x(1,i1,i2,i3)*fil(1,0,2)
           tt1_2=x(1,i1+1,i2,i3)*fil(1,0,2)
           tt2=x(2,i1,i2,i3)*fil(1,0,2)
           tt2_2=x(2,i1+1,i2,i3)*fil(1,0,2)
           tt3=0
           tt3_2=0
           tt4=0
           tt4_2=0
           do l=1,lupfil
              j1=modulo(i2+l,n2+1)
              j2=i2-l
              tt1=tt1+(x(1,i1,j1,i3)+x(1,i1,j2,i3))*fil(1,l,2)
              tt1_2=tt1_2+(x(1,i1+1,j1,i3)+x(1,i1+1,j2,i3))*fil(1,l,2)
              tt3=tt3+(x(2,i1,j1,i3)-x(2,i1,j2,i3))*fil(2,l,2)
              tt3_2=tt3_2+(x(2,i1+1,j1,i3)-x(2,i1+1,j2,i3))*fil(2,l,2)
              tt2=tt2+(x(2,i1,j1,i3)+x(2,i1,j2,i3))*fil(1,l,2)
              tt2_2=tt2_2+(x(2,i1+1,j1,i3)+x(2,i1+1,j2,i3))*fil(1,l,2)
              tt4=tt4+(x(1,i1,j1,i3)-x(1,i1,j2,i3))*fil(2,l,2)
              tt4_2=tt4_2+(x(1,i1+1,j1,i3)-x(1,i1+1,j2,i3))*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(1,i1+1,i2,i3)=y(1,i1+1,i2,i3)+tt1_2-tt3_2
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
           y(2,i1+1,i2,i3)=y(2,i1+1,i2,i3)+tt2_2+tt4_2
        enddo

     enddo
     do i1=((n1+1)/2)*2,n1,1
        do i2=0,lupfil-1
           tt1=x(1,i1,i2,i3)*fil(1,0,2)
           tt2=x(2,i1,i2,i3)*fil(1,0,2)
           tt3=0
           tt4=0
           do l=1,lupfil
              j1=i2+l
              j2=modulo(i2-l,n2+1)
              tt1=tt1+(x(1,i1,j1,i3)+x(1,i1,j2,i3))*fil(1,l,2)
              tt3=tt3+(x(2,i1,j1,i3)-x(2,i1,j2,i3))*fil(2,l,2)
              tt2=tt2+(x(2,i1,j1,i3)+x(2,i1,j2,i3))*fil(1,l,2)
              tt4=tt4+(x(1,i1,j1,i3)-x(1,i1,j2,i3))*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
        enddo
        do i2=lupfil,n2-lupfil
           tt1=x(1,i1,i2,i3)*fil(1,0,2)
           tt2=x(2,i1,i2,i3)*fil(1,0,2)
           tt3=0
           tt4=0
           do l=1,lupfil
              j1=i2+l
              j2=i2-l
              tt1=tt1+(x(1,i1,j1,i3)+x(1,i1,j2,i3))*fil(1,l,2)
              tt3=tt3+(x(2,i1,j1,i3)-x(2,i1,j2,i3))*fil(2,l,2)
              tt2=tt2+(x(2,i1,j1,i3)+x(2,i1,j2,i3))*fil(1,l,2)
              tt4=tt4+(x(1,i1,j1,i3)-x(1,i1,j2,i3))*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
        enddo
        do i2=n2-lupfil+1,n2
           tt1=x(1,i1,i2,i3)*fil(1,0,2)
           tt2=x(2,i1,i2,i3)*fil(1,0,2)
           tt3=0
           tt4=0
           do l=1,lupfil
              j1=modulo(i2+l,n2+1)
              j2=i2-l
              tt1=tt1+(x(1,i1,j1,i3)+x(1,i1,j2,i3))*fil(1,l,2)
              tt3=tt3+(x(2,i1,j1,i3)-x(2,i1,j2,i3))*fil(2,l,2)
              tt2=tt2+(x(2,i1,j1,i3)+x(2,i1,j2,i3))*fil(1,l,2)
              tt4=tt4+(x(1,i1,j1,i3)-x(1,i1,j2,i3))*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
        enddo

     enddo
  enddo
!$omp enddo
!$omp do

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1-1,2
        do i3=0,lupfil-1
           tt1=x(1,i1,i2,i3)*fil(1,0,3)
           tt1_2=x(1,i1+1,i2,i3)*fil(1,0,3)
           tt2=x(2,i1,i2,i3)*fil(1,0,3)
           tt2_2=x(2,i1+1,i2,i3)*fil(1,0,3)
           tt3=0
           tt3_2=0
           tt4=0
           tt4_2=0
           do l=1,lupfil
              j2=modulo(i3-l,n3+1)
              j1=i3+l
              tt1=tt1+(x(1,i1,i2,j1)+x(1,i1,i2,j2))*fil(1,l,3)
              tt1_2=tt1_2+(x(1,i1+1,i2,j1)+x(1,i1+1,i2,j2))*fil(1,l,3)
              tt3=tt3+(x(2,i1,i2,j1)-x(2,i1,i2,j2))*fil(2,l,3)
              tt3_2=tt3_2+(x(2,i1+1,i2,j1)-x(2,i1+1,i2,j2))*fil(2,l,3)
              tt2=tt2+(x(2,i1,i2,j1)+x(2,i1,i2,j2))*fil(1,l,3)
              tt2_2=tt2_2+(x(2,i1+1,i2,j1)+x(2,i1+1,i2,j2))*fil(1,l,3)
              tt4=tt4+(x(1,i1,i2,j1)-x(1,i1,i2,j2))*fil(2,l,3)
              tt4_2=tt4_2+(x(1,i1+1,i2,j1)-x(1,i1+1,i2,j2))*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(1,i1+1,i2,i3)=y(1,i1+1,i2,i3)+tt1_2-tt3_2
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
           y(2,i1+1,i2,i3)=y(2,i1+1,i2,i3)+tt2_2+tt4_2
        enddo
        do i3=lupfil,n3-lupfil
           tt1=x(1,i1,i2,i3)*fil(1,0,3)
           tt1_2=x(1,i1+1,i2,i3)*fil(1,0,3)
           tt2=x(2,i1,i2,i3)*fil(1,0,3)
           tt2_2=x(2,i1+1,i2,i3)*fil(1,0,3)
           tt3=0
           tt3_2=0
           tt4=0
           tt4_2=0
           do l=1,lupfil
              j2=i3-l
              j1=i3+l
              tt1=tt1+(x(1,i1,i2,j1)+x(1,i1,i2,j2))*fil(1,l,3)
              tt1_2=tt1_2+(x(1,i1+1,i2,j1)+x(1,i1+1,i2,j2))*fil(1,l,3)
              tt3=tt3+(x(2,i1,i2,j1)-x(2,i1,i2,j2))*fil(2,l,3)
              tt3_2=tt3_2+(x(2,i1+1,i2,j1)-x(2,i1+1,i2,j2))*fil(2,l,3)
              tt2=tt2+(x(2,i1,i2,j1)+x(2,i1,i2,j2))*fil(1,l,3)
              tt2_2=tt2_2+(x(2,i1+1,i2,j1)+x(2,i1+1,i2,j2))*fil(1,l,3)
              tt4=tt4+(x(1,i1,i2,j1)-x(1,i1,i2,j2))*fil(2,l,3)
              tt4_2=tt4_2+(x(1,i1+1,i2,j1)-x(1,i1+1,i2,j2))*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(1,i1+1,i2,i3)=y(1,i1+1,i2,i3)+tt1_2-tt3_2
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
           y(2,i1+1,i2,i3)=y(2,i1+1,i2,i3)+tt2_2+tt4_2
        enddo
        do i3=n3-lupfil+1,n3
           tt1=x(1,i1,i2,i3)*fil(1,0,3)
           tt1_2=x(1,i1+1,i2,i3)*fil(1,0,3)
           tt2=x(2,i1,i2,i3)*fil(1,0,3)
           tt2_2=x(2,i1+1,i2,i3)*fil(1,0,3)
           tt3=0
           tt3_2=0
           tt4=0
           tt4_2=0
           do l=1,lupfil
              j1=modulo(i3+l,n3+1)
              j2=i3-l
              tt1=tt1+(x(1,i1,i2,j1)+x(1,i1,i2,j2))*fil(1,l,3)
              tt1_2=tt1_2+(x(1,i1+1,i2,j1)+x(1,i1+1,i2,j2))*fil(1,l,3)
              tt3=tt3+(x(2,i1,i2,j1)-x(2,i1,i2,j2))*fil(2,l,3)
              tt3_2=tt3_2+(x(2,i1+1,i2,j1)-x(2,i1+1,i2,j2))*fil(2,l,3)
              tt2=tt2+(x(2,i1,i2,j1)+x(2,i1,i2,j2))*fil(1,l,3)
              tt2_2=tt2_2+(x(2,i1+1,i2,j1)+x(2,i1+1,i2,j2))*fil(1,l,3)
              tt4=tt4+(x(1,i1,i2,j1)-x(1,i1,i2,j2))*fil(2,l,3)
              tt4_2=tt4_2+(x(1,i1+1,i2,j1)-x(1,i1+1,i2,j2))*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(1,i1+1,i2,i3)=y(1,i1+1,i2,i3)+tt1_2-tt3_2
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
           y(2,i1+1,i2,i3)=y(2,i1+1,i2,i3)+tt2_2+tt4_2
        enddo
     enddo
     do i1=((n1+1)/2)*2,n1,1
        do i3=0,lupfil-1
           tt1=x(1,i1,i2,i3)*fil(1,0,3)
           tt2=x(2,i1,i2,i3)*fil(1,0,3)
           tt3=0
           tt4=0
           do l=1,lupfil
              j2=modulo(i3-l,n3+1)
              j1=i3+l
              tt1=tt1+(x(1,i1,i2,j1)+x(1,i1,i2,j2))*fil(1,l,3)
              tt3=tt3+(x(2,i1,i2,j1)-x(2,i1,i2,j2))*fil(2,l,3)
              tt2=tt2+(x(2,i1,i2,j1)+x(2,i1,i2,j2))*fil(1,l,3)
              tt4=tt4+(x(1,i1,i2,j1)-x(1,i1,i2,j2))*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
        enddo
        do i3=lupfil,n3-lupfil
           tt1=x(1,i1,i2,i3)*fil(1,0,3)
           tt2=x(2,i1,i2,i3)*fil(1,0,3)
           tt3=0
           tt4=0
           do l=1,lupfil
              j2=i3-l
              j1=i3+l
              tt1=tt1+(x(1,i1,i2,j1)+x(1,i1,i2,j2))*fil(1,l,3)
              tt3=tt3+(x(2,i1,i2,j1)-x(2,i1,i2,j2))*fil(2,l,3)
              tt2=tt2+(x(2,i1,i2,j1)+x(2,i1,i2,j2))*fil(1,l,3)
              tt4=tt4+(x(1,i1,i2,j1)-x(1,i1,i2,j2))*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
        enddo
        do i3=n3-lupfil+1,n3
           tt1=x(1,i1,i2,i3)*fil(1,0,3)
           tt2=x(2,i1,i2,i3)*fil(1,0,3)
           tt3=0
           tt4=0
           do l=1,lupfil
              j1=modulo(i3+l,n3+1)
              j2=i3-l
              tt1=tt1+(x(1,i1,i2,j1)+x(1,i1,i2,j2))*fil(1,l,3)
              tt3=tt3+(x(2,i1,i2,j1)-x(2,i1,i2,j2))*fil(2,l,3)
              tt2=tt2+(x(2,i1,i2,j1)+x(2,i1,i2,j2))*fil(1,l,3)
              tt4=tt4+(x(1,i1,i2,j1)-x(1,i1,i2,j2))*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1-tt3
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2+tt4
        enddo
     enddo
  enddo
!$omp enddo
!$omp end parallel  
END SUBROUTINE convolut_kinetic_per_c_k

subroutine convolut_kinetic_per_c_k_notranspose(n1,n2,n3,hgrid,x,y,c_in,k1,k2,k3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c_in,k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3,2), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3,2), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,l,j1,j2
  real(wp) :: tt1,tt2,tt3,tt4,c
  real(wp) :: tt1_2,tt2_2,tt3_2,tt4_2
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,0:lupfil,3) :: fil  

  scale (:)=real(-.5_gp/hgrid(:)**2,wp)
  scale1(1)=real(k1/hgrid(1),wp)
  scale1(2)=real(k2/hgrid(2),wp)
  scale1(3)=real(k3/hgrid(3),wp)
  c=c_in+.5_wp*(k1*k1+k2*k2+k3*k3)

  ! second derivative filters for Daubechies 16
  fil(1,0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(1,2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(1,3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(1,4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(1,5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(1,6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(1,7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(1,8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(1,9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(1,10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(1,11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(1,12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(1,13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(1,14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  fil(2,0,:)= 0._wp
  fil(2,1,:)= 0.8834460460908270942785856e0_wp*scale1(:)
  fil(2,2,:)= -0.3032593514765938346887962e0_wp*scale1(:)
  fil(2,3,:)= 0.1063640682894442760934532e0_wp*scale1(:)
  fil(2,4,:)= -0.03129014783948023634381564e0_wp*scale1(:)
  fil(2,5,:)= 0.006958379116450707495020408e0_wp*scale1(:)
  fil(2,6,:)= -0.001031530213375445369097965e0_wp*scale1(:)
  fil(2,7,:)= 0.00007667706908380351933901775e0_wp*scale1(:)
  fil(2,8,:)= 2.451992111053665419191564e-7_wp*scale1(:)
  fil(2,9,:)= 3.993810456408053712133667e-8_wp*scale1(:)
  fil(2,10,:)=-7.207948238588481597101904e-8_wp*scale1(:)
  fil(2,11,:)=-9.697184925637300947553069e-10_wp*scale1(:)
  fil(2,12,:)=-7.252206916665149851135592e-13_wp*scale1(:)
  fil(2,13,:)=1.240078536096648534547439e-14_wp*scale1(:)
  fil(2,14,:)=-1.585464751677102510097179e-19_wp*scale1(:)


!$omp parallel default (private) shared(n1,n2,n3,x,y,c,fil)
!$omp do

  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,lupfil-1
           tt1=x(i1,i2,i3,1)*c
           tt2=x(i1,i2,i3,2)*c
           tt1=tt1+x(i1,i2,i3,1)*fil(1,0,1)
           tt2=tt2+x(i1,i2,i3,2)*fil(1,0,1)
           tt3=0
           tt4=0
           do l=1,14
              j1=i1+l
              j2=modulo(i1-l,n1+1)
              tt1=tt1+(x(j1,i2,i3,1)+x(j2,i2,i3,1))*fil(1,l,1)
              tt3=tt3+(x(j1,i2,i3,2)-x(j2,i2,i3,2))*fil(2,l,1)
              tt2=tt2+(x(j1,i2,i3,2)+x(j2,i2,i3,2))*fil(1,l,1)
              tt4=tt4+(x(j1,i2,i3,1)-x(j2,i2,i3,1))*fil(2,l,1)
           enddo
           y(i1,i2,i3,1)=tt1-tt3
           y(i1,i2,i3,2)=tt2+tt4
        enddo
        do i1=lupfil,n1-lupfil
           tt1=x(i1,i2,i3,1)*c
           tt2=x(i1,i2,i3,2)*c
           tt1=tt1+x(i1,i2,i3,1)*fil(1,0,1)
           tt2=tt2+x(i1,i2,i3,2)*fil(1,0,1)
           tt3=0
           tt4=0
           do l=1,14
              j1=i1+l
              j2=i1-l
              tt1=tt1+(x(j1,i2,i3,1)+x(j2,i2,i3,1))*fil(1,l,1)
              tt3=tt3+(x(j1,i2,i3,2)-x(j2,i2,i3,2))*fil(2,l,1)
              tt2=tt2+(x(j1,i2,i3,2)+x(j2,i2,i3,2))*fil(1,l,1)
              tt4=tt4+(x(j1,i2,i3,1)-x(j2,i2,i3,1))*fil(2,l,1)
           enddo
           y(i1,i2,i3,1)=tt1-tt3
           y(i1,i2,i3,2)=tt2+tt4
        enddo
        do i1=n1-lupfil+1,n1
           tt1=x(i1,i2,i3,1)*c
           tt2=x(i1,i2,i3,2)*c
           tt1=tt1+x(i1,i2,i3,1)*fil(1,0,1)
           tt2=tt2+x(i1,i2,i3,2)*fil(1,0,1)
           tt3=0.0_wp
           tt4=0.0_wp
           do l=1,14
              j1=modulo(i1+l,n1+1)
              j2=i1-l
              tt1=tt1+(x(j1,i2,i3,1)+x(j2,i2,i3,1))*fil(1,l,1)
              tt3=tt3+(x(j1,i2,i3,2)-x(j2,i2,i3,2))*fil(2,l,1)
              tt2=tt2+(x(j1,i2,i3,2)+x(j2,i2,i3,2))*fil(1,l,1)
              tt4=tt4+(x(j1,i2,i3,1)-x(j2,i2,i3,1))*fil(2,l,1)
           enddo
           y(i1,i2,i3,1)=tt1-tt3
           y(i1,i2,i3,2)=tt2+tt4
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1-1,2
        do i2=0,lupfil-1
           tt1  =x(i1,i2,i3  ,1)*fil(1,0,2)
           tt1_2=x(i1+1,i2,i3,1)*fil(1,0,2)
           tt2  =x(i1,i2,i3  ,2)*fil(1,0,2)
           tt2_2=x(i1+1,i2,i3,2)*fil(1,0,2)
           tt3  =0.0_wp
           tt3_2=0.0_wp
           tt4  =0.0_wp
           tt4_2=0.0_wp
           do l=1,lupfil
              j1=i2+l
              j2=modulo(i2-l,n2+1)
              tt1  =tt1  +(x(i1,j1,i3  ,1)+x(i1,j2,i3  ,1))*fil(1,l,2)
              tt1_2=tt1_2+(x(i1+1,j1,i3,1)+x(i1+1,j2,i3,1))*fil(1,l,2)
              tt3  =tt3  +(x(i1,j1,i3  ,2)-x(i1,j2,i3  ,2))*fil(2,l,2)
              tt3_2=tt3_2+(x(i1+1,j1,i3,2)-x(i1+1,j2,i3,2))*fil(2,l,2)
              tt2  =tt2  +(x(i1,j1,i3  ,2)+x(i1,j2,i3  ,2))*fil(1,l,2)
              tt2_2=tt2_2+(x(i1+1,j1,i3,2)+x(i1+1,j2,i3,2))*fil(1,l,2)
              tt4  =tt4  +(x(i1,j1,i3  ,1)-x(i1,j2,i3  ,1))*fil(2,l,2)
              tt4_2=tt4_2+(x(i1+1,j1,i3,1)-x(i1+1,j2,i3,1))*fil(2,l,2)
           enddo
           y(i1  ,i2,i3,1)=y(i1  ,i2,i3,1)+tt1-tt3
           y(i1+1,i2,i3,1)=y(i1+1,i2,i3,1)+tt1_2-tt3_2
           y(i1  ,i2,i3,2)=y(i1  ,i2,i3,2)+tt2+tt4
           y(i1+1,i2,i3,2)=y(i1+1,i2,i3,2)+tt2_2+tt4_2
        enddo
        do i2=lupfil,n2-lupfil
           tt1  =x(i1  ,i2,i3,1)*fil(1,0,2)
           tt1_2=x(i1+1,i2,i3,1)*fil(1,0,2)
           tt2  =x(i1  ,i2,i3,2)*fil(1,0,2)
           tt2_2=x(i1+1,i2,i3,2)*fil(1,0,2)
           tt3=0.0_wp
           tt3_2=0_wp
           tt4=0_wp
           tt4_2=0_wp
           do l=1,lupfil
              j1=i2+l
              j2=i2-l
              tt1  =tt1  +(x(i1  ,j1,i3,1)+x(i1  ,j2,i3,1))*fil(1,l,2)
              tt1_2=tt1_2+(x(i1+1,j1,i3,1)+x(i1+1,j2,i3,1))*fil(1,l,2)
              tt3  =tt3  +(x(i1  ,j1,i3,2)-x(i1  ,j2,i3,2))*fil(2,l,2)
              tt3_2=tt3_2+(x(i1+1,j1,i3,2)-x(i1+1,j2,i3,2))*fil(2,l,2)
              tt2  =tt2  +(x(i1  ,j1,i3,2)+x(i1  ,j2,i3,2))*fil(1,l,2)
              tt2_2=tt2_2+(x(i1+1,j1,i3,2)+x(i1+1,j2,i3,2))*fil(1,l,2)
              tt4  =tt4  +(x(i1  ,j1,i3,1)-x(i1  ,j2,i3,1))*fil(2,l,2)
              tt4_2=tt4_2+(x(i1+1,j1,i3,1)-x(i1+1,j2,i3,1))*fil(2,l,2)
           enddo
           y(i1  ,i2,i3,1)=y(i1  ,i2,i3,1)+tt1-tt3
           y(i1+1,i2,i3,1)=y(i1+1,i2,i3,1)+tt1_2-tt3_2
           y(i1  ,i2,i3,2)=y(i1  ,i2,i3,2)+tt2+tt4
           y(i1+1,i2,i3,2)=y(i1+1,i2,i3,2)+tt2_2+tt4_2
        enddo
        do i2=n2-lupfil+1,n2
           tt1  =x(i1  ,i2,i3,1)*fil(1,0,2)
           tt1_2=x(i1+1,i2,i3,1)*fil(1,0,2)
           tt2  =x(i1  ,i2,i3,2)*fil(1,0,2)
           tt2_2=x(i1+1,i2,i3,2)*fil(1,0,2)
           tt3  =0_wp
           tt3_2=0_wp
           tt4  =0_wp
           tt4_2=0_wp
           do l=1,lupfil
              j1=modulo(i2+l,n2+1)
              j2=i2-l
              tt1  =tt1  +(x(i1  ,j1,i3,1)+x(i1  ,j2,i3,1))*fil(1,l,2)
              tt1_2=tt1_2+(x(i1+1,j1,i3,1)+x(i1+1,j2,i3,1))*fil(1,l,2)
              tt3  =tt3  +(x(i1  ,j1,i3,2)-x(i1  ,j2,i3,2))*fil(2,l,2)
              tt3_2=tt3_2+(x(i1+1,j1,i3,2)-x(i1+1,j2,i3,2))*fil(2,l,2)
              tt2  =tt2  +(x(i1  ,j1,i3,2)+x(i1  ,j2,i3,2))*fil(1,l,2)
              tt2_2=tt2_2+(x(i1+1,j1,i3,2)+x(i1+1,j2,i3,2))*fil(1,l,2)
              tt4  =tt4  +(x(i1  ,j1,i3,1)-x(i1  ,j2,i3,1))*fil(2,l,2)
              tt4_2=tt4_2+(x(i1+1,j1,i3,1)-x(i1+1,j2,i3,1))*fil(2,l,2)
           enddo
           y(i1  ,i2,i3,1)=y(i1  ,i2,i3,1)+tt1-tt3
           y(i1+1,i2,i3,1)=y(i1+1,i2,i3,1)+tt1_2-tt3_2
           y(i1  ,i2,i3,2)=y(i1  ,i2,i3,2)+tt2+tt4
           y(i1+1,i2,i3,2)=y(i1+1,i2,i3,2)+tt2_2+tt4_2
        enddo

     enddo
     do i1=((n1+1)/2)*2,n1,1
        do i2=0,lupfil-1
           tt1=x(i1,i2,i3,1)*fil(1,0,2)
           tt2=x(i1,i2,i3,2)*fil(1,0,2)
           tt3=0_wp
           tt4=0_wp
           do l=1,lupfil
              j1=i2+l
              j2=modulo(i2-l,n2+1)
              tt1=tt1+(x(i1,j1,i3,1)+x(i1,j2,i3,1))*fil(1,l,2)
              tt3=tt3+(x(i1,j1,i3,2)-x(i1,j2,i3,2))*fil(2,l,2)
              tt2=tt2+(x(i1,j1,i3,2)+x(i1,j2,i3,2))*fil(1,l,2)
              tt4=tt4+(x(i1,j1,i3,1)-x(i1,j2,i3,1))*fil(2,l,2)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1-tt3
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2+tt4
        enddo
        do i2=lupfil,n2-lupfil
           tt1=x(i1,i2,i3,1)*fil(1,0,2)
           tt2=x(i1,i2,i3,2)*fil(1,0,2)
           tt3=0_wp
           tt4=0_wp
           do l=1,lupfil
              j1=i2+l
              j2=i2-l
              tt1=tt1+(x(i1,j1,i3,1)+x(i1,j2,i3,1))*fil(1,l,2)
              tt3=tt3+(x(i1,j1,i3,2)-x(i1,j2,i3,2))*fil(2,l,2)
              tt2=tt2+(x(i1,j1,i3,2)+x(i1,j2,i3,2))*fil(1,l,2)
              tt4=tt4+(x(i1,j1,i3,1)-x(i1,j2,i3,1))*fil(2,l,2)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1-tt3
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2+tt4
        enddo
        do i2=n2-lupfil+1,n2
           tt1=x(i1,i2,i3,1)*fil(1,0,2)
           tt2=x(i1,i2,i3,2)*fil(1,0,2)
           tt3=0
           tt4=0
           do l=1,lupfil
              j1=modulo(i2+l,n2+1)
              j2=i2-l
              tt1=tt1+(x(i1,j1,i3,1)+x(i1,j2,i3,1))*fil(1,l,2)
              tt3=tt3+(x(i1,j1,i3,2)-x(i1,j2,i3,2))*fil(2,l,2)
              tt2=tt2+(x(i1,j1,i3,2)+x(i1,j2,i3,2))*fil(1,l,2)
              tt4=tt4+(x(i1,j1,i3,1)-x(i1,j2,i3,1))*fil(2,l,2)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1-tt3
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2+tt4
        enddo

     enddo
  enddo
!$omp enddo
!$omp do

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1-1,2
        do i3=0,lupfil-1
           tt1  =x(i1  ,i2,i3,1)*fil(1,0,3)
           tt1_2=x(i1+1,i2,i3,1)*fil(1,0,3)
           tt2  =x(i1  ,i2,i3,2)*fil(1,0,3)
           tt2_2=x(i1+1,i2,i3,2)*fil(1,0,3)
           tt3=0_wp
           tt3_2=0_wp
           tt4=0_wp
           tt4_2=0_wp
           do l=1,lupfil
              j2=modulo(i3-l,n3+1)
              j1=i3+l
              tt1  =tt1  +(x(i1  ,i2,j1,1)+x(i1  ,i2,j2,1))*fil(1,l,3)
              tt1_2=tt1_2+(x(i1+1,i2,j1,1)+x(i1+1,i2,j2,1))*fil(1,l,3)
              tt3  =tt3  +(x(i1  ,i2,j1,2)-x(i1  ,i2,j2,2))*fil(2,l,3)
              tt3_2=tt3_2+(x(i1+1,i2,j1,2)-x(i1+1,i2,j2,2))*fil(2,l,3)
              tt2  =tt2  +(x(i1  ,i2,j1,2)+x(i1  ,i2,j2,2))*fil(1,l,3)
              tt2_2=tt2_2+(x(i1+1,i2,j1,2)+x(i1+1,i2,j2,2))*fil(1,l,3)
              tt4  =tt4  +(x(i1  ,i2,j1,1)-x(i1  ,i2,j2,1))*fil(2,l,3)
              tt4_2=tt4_2+(x(i1+1,i2,j1,1)-x(i1+1,i2,j2,1))*fil(2,l,3)
           enddo
           y(i1  ,i2,i3,1)=y(i1  ,i2,i3,1)+tt1-tt3
           y(i1+1,i2,i3,1)=y(i1+1,i2,i3,1)+tt1_2-tt3_2
           y(i1  ,i2,i3,2)=y(i1  ,i2,i3,2)+tt2+tt4
           y(i1+1,i2,i3,2)=y(i1+1,i2,i3,2)+tt2_2+tt4_2
        enddo
        do i3=lupfil,n3-lupfil
           tt1  =x(i1  ,i2,i3,1)*fil(1,0,3)
           tt1_2=x(i1+1,i2,i3,1)*fil(1,0,3)
           tt2  =x(i1  ,i2,i3,2)*fil(1,0,3)
           tt2_2=x(i1+1,i2,i3,2)*fil(1,0,3)
           tt3=0_wp
           tt3_2=0_wp
           tt4=0_wp
           tt4_2=0_wp
           do l=1,lupfil
              j2=i3-l
              j1=i3+l
              tt1  =tt1  +(x(i1  ,i2,j1,1)+x(i1  ,i2,j2,1))*fil(1,l,3)
              tt1_2=tt1_2+(x(i1+1,i2,j1,1)+x(i1+1,i2,j2,1))*fil(1,l,3)
              tt3  =tt3  +(x(i1  ,i2,j1,2)-x(i1  ,i2,j2,2))*fil(2,l,3)
              tt3_2=tt3_2+(x(i1+1,i2,j1,2)-x(i1+1,i2,j2,2))*fil(2,l,3)
              tt2  =tt2  +(x(i1  ,i2,j1,2)+x(i1  ,i2,j2,2))*fil(1,l,3)
              tt2_2=tt2_2+(x(i1+1,i2,j1,2)+x(i1+1,i2,j2,2))*fil(1,l,3)
              tt4  =tt4  +(x(i1  ,i2,j1,1)-x(i1  ,i2,j2,1))*fil(2,l,3)
              tt4_2=tt4_2+(x(i1+1,i2,j1,1)-x(i1+1,i2,j2,1))*fil(2,l,3)
           enddo
           y(i1  ,i2,i3,1)=y(i1  ,i2,i3,1)+tt1-tt3
           y(i1+1,i2,i3,1)=y(i1+1,i2,i3,1)+tt1_2-tt3_2
           y(i1  ,i2,i3,2)=y(i1  ,i2,i3,2)+tt2+tt4
           y(i1+1,i2,i3,2)=y(i1+1,i2,i3,2)+tt2_2+tt4_2
        enddo
        do i3=n3-lupfil+1,n3
           tt1  =x(i1  ,i2,i3,1)*fil(1,0,3)
           tt1_2=x(i1+1,i2,i3,1)*fil(1,0,3)
           tt2  =x(i1  ,i2,i3,2)*fil(1,0,3)
           tt2_2=x(i1+1,i2,i3,2)*fil(1,0,3)
           tt3=0_wp
           tt3_2=0_wp
           tt4=0_wp
           tt4_2=0_wp
           do l=1,lupfil
              j1=modulo(i3+l,n3+1)
              j2=i3-l
              tt1  =tt1  +(x(i1  ,i2,j1,1)+x(i1  ,i2,j2,1))*fil(1,l,3)
              tt1_2=tt1_2+(x(i1+1,i2,j1,1)+x(i1+1,i2,j2,1))*fil(1,l,3)
              tt3  =tt3  +(x(i1  ,i2,j1,2)-x(i1  ,i2,j2,2))*fil(2,l,3)
              tt3_2=tt3_2+(x(i1+1,i2,j1,2)-x(i1+1,i2,j2,2))*fil(2,l,3)
              tt2  =tt2  +(x(i1  ,i2,j1,2)+x(i1  ,i2,j2,2))*fil(1,l,3)
              tt2_2=tt2_2+(x(i1+1,i2,j1,2)+x(i1+1,i2,j2,2))*fil(1,l,3)
              tt4  =tt4  +(x(i1  ,i2,j1,1)-x(i1  ,i2,j2,1))*fil(2,l,3)
              tt4_2=tt4_2+(x(i1+1,i2,j1,1)-x(i1+1,i2,j2,1))*fil(2,l,3)
           enddo
           y(i1  ,i2,i3,1)=y(i1  ,i2,i3,1)+tt1-tt3
           y(i1+1,i2,i3,1)=y(i1+1,i2,i3,1)+tt1_2-tt3_2
           y(i1  ,i2,i3,2)=y(i1  ,i2,i3,2)+tt2+tt4
           y(i1+1,i2,i3,2)=y(i1+1,i2,i3,2)+tt2_2+tt4_2
        enddo
     enddo
     do i1=((n1+1)/2)*2,n1,1
        do i3=0,lupfil-1
           tt1=x(i1,i2,i3,1)*fil(1,0,3)
           tt2=x(i1,i2,i3,2)*fil(1,0,3)
           tt3=0_wp
           tt4=0_wp
           do l=1,lupfil
              j2=modulo(i3-l,n3+1)
              j1=i3+l
              tt1=tt1+(x(i1,i2,j1,1)+x(i1,i2,j2,1))*fil(1,l,3)
              tt3=tt3+(x(i1,i2,j1,2)-x(i1,i2,j2,2))*fil(2,l,3)
              tt2=tt2+(x(i1,i2,j1,2)+x(i1,i2,j2,2))*fil(1,l,3)
              tt4=tt4+(x(i1,i2,j1,1)-x(i1,i2,j2,1))*fil(2,l,3)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1-tt3
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2+tt4
        enddo
        do i3=lupfil,n3-lupfil
           tt1=x(i1,i2,i3,1)*fil(1,0,3)
           tt2=x(i1,i2,i3,2)*fil(1,0,3)
           tt3=0_wp
           tt4=0_wp
           do l=1,lupfil
              j2=i3-l
              j1=i3+l
              tt1=tt1+(x(i1,i2,j1,1)+x(i1,i2,j2,1))*fil(1,l,3)
              tt3=tt3+(x(i1,i2,j1,2)-x(i1,i2,j2,2))*fil(2,l,3)
              tt2=tt2+(x(i1,i2,j1,2)+x(i1,i2,j2,2))*fil(1,l,3)
              tt4=tt4+(x(i1,i2,j1,1)-x(i1,i2,j2,1))*fil(2,l,3)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1-tt3
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2+tt4
        enddo
        do i3=n3-lupfil+1,n3
           tt1=x(i1,i2,i3,1)*fil(1,0,3)
           tt2=x(i1,i2,i3,2)*fil(1,0,3)
           tt3=0_wp
           tt4=0_wp
           do l=1,lupfil
              j1=modulo(i3+l,n3+1)
              j2=i3-l
              tt1=tt1+(x(i1,i2,j1,1)+x(i1,i2,j2,1))*fil(1,l,3)
              tt3=tt3+(x(i1,i2,j1,2)-x(i1,i2,j2,2))*fil(2,l,3)
              tt2=tt2+(x(i1,i2,j1,2)+x(i1,i2,j2,2))*fil(1,l,3)
              tt4=tt4+(x(i1,i2,j1,1)-x(i1,i2,j2,1))*fil(2,l,3)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1-tt3
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2+tt4
        enddo
     enddo
  enddo
!$omp enddo
!$omp end parallel  
END SUBROUTINE convolut_kinetic_per_c_k_notranspose




!>   Applies the modified kinetic energy operator onto x to get y. 
!!   Computes the kinetic energy too.
!!   Works for periodic BC.
!!   Modified kinetic energy operator:
!!   A=-1/2 exp(Ikr) Delta exp(-Ikr)
!!   where k=(k1,k2,k3); r=(x,y,z)
!!
!! 
subroutine convolut_kinetic_per_T_k(n1,n2,n3,hgrid,x,y,kstrten,k1,k2,k3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(inout) :: y
  real(wp), dimension(6), intent(out):: kstrten
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt1,tt2,cx,cy,cz,kstrt1,kstrt2,kstrt3
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,lowfil:lupfil,3) :: fil  

  scale (:)=real(-.5_gp/hgrid(:)**2,wp)
  scale1(1)=real(k1/hgrid(1),wp)
  scale1(2)=real(k2/hgrid(2),wp)
  scale1(3)=real(k3/hgrid(3),wp)
  cx=0.5_wp*k1*k1
  cy=0.5_wp*k2*k2
  cz=0.5_wp*k3*k3

  ! second derivative filters for Daubechies 16
  fil(1,0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(1,2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(1,3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(1,4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(1,5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(1,6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(1,7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(1,8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(1,9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(1,10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(1,11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(1,12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(1,13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(1,14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(1,-i,:)=fil(1,i,:)
  enddo

  !first derivative filters for Daubechies 16
  fil(2,0,:)= 0._wp
  fil(2,1,:)= 0.8834460460908270942785856e0_wp*scale1(:)
  fil(2,2,:)= -0.3032593514765938346887962e0_wp*scale1(:)
  fil(2,3,:)= 0.1063640682894442760934532e0_wp*scale1(:)
  fil(2,4,:)= -0.03129014783948023634381564e0_wp*scale1(:)
  fil(2,5,:)= 0.006958379116450707495020408e0_wp*scale1(:)
  fil(2,6,:)= -0.001031530213375445369097965e0_wp*scale1(:)
  fil(2,7,:)= 0.00007667706908380351933901775e0_wp*scale1(:)
  fil(2,8,:)= 2.451992111053665419191564e-7_wp*scale1(:)
  fil(2,9,:)= 3.993810456408053712133667e-8_wp*scale1(:)
  fil(2,10,:)=-7.207948238588481597101904e-8_wp*scale1(:)
  fil(2,11,:)=-9.697184925637300947553069e-10_wp*scale1(:)
  fil(2,12,:)=-7.252206916665149851135592e-13_wp*scale1(:)
  fil(2,13,:)=1.240078536096648534547439e-14_wp*scale1(:)
  fil(2,14,:)=-1.585464751677102510097179e-19_wp*scale1(:)

  do i=1,14
     fil(2,-i,:)=-fil(2,i,:)
  enddo

  !sequence for kinetic stress tensors
  ! 11 22 33 12 13 23
  kstrten(1:6)=0.0_wp

  !ener=0._wp
!$omp parallel default (private) shared(x,y,kstrten,fil,cx,cy,cz,n1,n2,n3)
  kstrt1=0.0_wp
  kstrt2=0.0_wp
  kstrt3=0.0_wp

!$omp do 
  do i3=0,n3
     ! (1/2) |d/dx+ik_x)|^2
     do i2=0,n2
        do i1=0,n1
!           tt1=x(1,i1,i2,i3)*c
!           tt2=x(2,i1,i2,i3)*c
           tt1=x(1,i1,i2,i3)*cx
           tt2=x(2,i1,i2,i3)*cx
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt1=tt1+x(1,j,i2,i3)*fil(1,l,1)-x(2,j,i2,i3)*fil(2,l,1)
              tt2=tt2+x(2,j,i2,i3)*fil(1,l,1)+x(1,j,i2,i3)*fil(2,l,1)
           enddo
           !add the sum with the previous existing object
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2

           !ener=ener+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
           kstrt1=kstrt1+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
!           tt1=0._wp
!           tt2=0._wp
           tt1=x(1,i1,i2,i3)*cy
           tt2=x(2,i1,i2,i3)*cy
           do l=lowfil,lupfil
              j=modulo(i2+l,n2+1)
              tt1=tt1+x(1,i1,j,i3)*fil(1,l,2)-x(2,i1,j,i3)*fil(2,l,2)
              tt2=tt2+x(2,i1,j,i3)*fil(1,l,2)+x(1,i1,j,i3)*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2

           kstrt2=kstrt2+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
     
  enddo
!$omp enddo

!$omp do 
  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt1=x(1,i1,i2,i3)*cz
           tt2=x(2,i1,i2,i3)*cz
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt1=tt1+x(1,i1,i2,j)*fil(1,l,3)-x(2,i1,i2,j)*fil(2,l,3)
              tt2=tt2+x(2,i1,i2,j)*fil(1,l,3)+x(1,i1,i2,j)*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2

           kstrt3=kstrt3+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
  enddo
!$omp enddo

!$omp critical
  kstrten(1)=kstrten(1)+kstrt1
  kstrten(2)=kstrten(2)+kstrt2
  kstrten(3)=kstrten(3)+kstrt3
!$omp end critical

!$omp end parallel  
END SUBROUTINE convolut_kinetic_per_T_k

subroutine convolut_kinetic_per_T_k_notranspose(n1,n2,n3,hgrid,x,y,kstrten,k1,k2,k3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3,2), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3,2), intent(inout) :: y
  real(wp), dimension(6), intent(out):: kstrten
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt1,tt2,cx,cy,cz,kstrt1,kstrt2,kstrt3
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,lowfil:lupfil,3) :: fil  

  scale (:)=real(-.5_gp/hgrid(:)**2,wp)
  scale1(1)=real(k1/hgrid(1),wp)
  scale1(2)=real(k2/hgrid(2),wp)
  scale1(3)=real(k3/hgrid(3),wp)
  cx=0.5_wp*k1*k1
  cy=0.5_wp*k2*k2
  cz=0.5_wp*k3*k3

  ! second derivative filters for Daubechies 16
  fil(1,0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(1,2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(1,3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(1,4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(1,5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(1,6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(1,7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(1,8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(1,9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(1,10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(1,11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(1,12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(1,13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(1,14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(1,-i,:)=fil(1,i,:)
  enddo

  !first derivative filters for Daubechies 16
  fil(2,0,:)= 0._wp
  fil(2,1,:)= 0.8834460460908270942785856e0_wp*scale1(:)
  fil(2,2,:)= -0.3032593514765938346887962e0_wp*scale1(:)
  fil(2,3,:)= 0.1063640682894442760934532e0_wp*scale1(:)
  fil(2,4,:)= -0.03129014783948023634381564e0_wp*scale1(:)
  fil(2,5,:)= 0.006958379116450707495020408e0_wp*scale1(:)
  fil(2,6,:)= -0.001031530213375445369097965e0_wp*scale1(:)
  fil(2,7,:)= 0.00007667706908380351933901775e0_wp*scale1(:)
  fil(2,8,:)= 2.451992111053665419191564e-7_wp*scale1(:)
  fil(2,9,:)= 3.993810456408053712133667e-8_wp*scale1(:)
  fil(2,10,:)=-7.207948238588481597101904e-8_wp*scale1(:)
  fil(2,11,:)=-9.697184925637300947553069e-10_wp*scale1(:)
  fil(2,12,:)=-7.252206916665149851135592e-13_wp*scale1(:)
  fil(2,13,:)=1.240078536096648534547439e-14_wp*scale1(:)
  fil(2,14,:)=-1.585464751677102510097179e-19_wp*scale1(:)

  do i=1,14
     fil(2,-i,:)=-fil(2,i,:)
  enddo

  !sequence for kinetic stress tensors
  ! 11 22 33 12 13 23
  kstrten(1:6)=0.0_wp

  !ener=0._wp
!$omp parallel default (private) shared(x,y,kstrten,fil,cx,cy,cz,n1,n2,n3)
  kstrt1=0.0_wp
  kstrt2=0.0_wp
  kstrt3=0.0_wp

!$omp do 
  do i3=0,n3
     ! (1/2) |d/dx+ik_x)|^2
     do i2=0,n2
        do i1=0,n1
           tt1=x(i1,i2,i3,1)*cx
           tt2=x(i1,i2,i3,2)*cx
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt1=tt1+x(j,i2,i3,1)*fil(1,l,1)-x(j,i2,i3,2)*fil(2,l,1)
              tt2=tt2+x(j,i2,i3,2)*fil(1,l,1)+x(j,i2,i3,1)*fil(2,l,1)
           enddo
           !add the sum with the previous existing object
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2

           !ener=ener+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
           kstrt1=kstrt1+tt1*x(i1,i2,i3,1)+tt2*x(i1,i2,i3,2)
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt1=x(i1,i2,i3,1)*cy
           tt2=x(i1,i2,i3,2)*cy
           do l=lowfil,lupfil
              j=modulo(i2+l,n2+1)
              tt1=tt1+x(i1,j,i3,1)*fil(1,l,2)-x(i1,j,i3,2)*fil(2,l,2)
              tt2=tt2+x(i1,j,i3,2)*fil(1,l,2)+x(i1,j,i3,1)*fil(2,l,2)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2

           kstrt2=kstrt2+tt1*x(i1,i2,i3,1)+tt2*x(i1,i2,i3,2)
        enddo
     enddo
     
  enddo
!$omp enddo

!$omp do 
  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt1=x(i1,i2,i3,1)*cz
           tt2=x(i1,i2,i3,2)*cz
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt1=tt1+x(i1,i2,j,1)*fil(1,l,3)-x(i1,i2,j,2)*fil(2,l,3)
              tt2=tt2+x(i1,i2,j,2)*fil(1,l,3)+x(i1,i2,j,1)*fil(2,l,3)
           enddo
           y(i1,i2,i3,1)=y(i1,i2,i3,1)+tt1
           y(i1,i2,i3,2)=y(i1,i2,i3,2)+tt2

           kstrt3=kstrt3+tt1*x(i1,i2,i3,1)+tt2*x(i1,i2,i3,2)
        enddo
     enddo
  enddo
!$omp enddo

!$omp critical
  kstrten(1)=kstrten(1)+kstrt1
  kstrten(2)=kstrten(2)+kstrt2
  kstrten(3)=kstrten(3)+kstrt3
!$omp end critical

!$omp end parallel  
END SUBROUTINE convolut_kinetic_per_T_k_notranspose

