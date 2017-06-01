!> @file
!!  Non-optimized convolution routines for kinetic operator
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Applies the modified kinetic energy operator onto x to get y. 
!! Works for the slab BC.
!! Modified kinetic energy operator:
!! A=-1/2 exp(Ikr) Delta exp(-Ikr)+C
!! where k=(k1,k2,k3); r=(x,y,z)
subroutine convolut_kinetic_slab_c_k(n1,n2,n3,hgrid,x,y,c_in,k1,k2,k3)

  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c_in,k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt1,tt2,c
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,lowfil:lupfil,3) :: fil  

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

  do i=1,14
     fil(1,-i,:)=fil(1,i,:)
  enddo

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
!$omp parallel default (private) shared(n1,n2,n3,x,y,c,fil)
!$omp do

  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt1=x(1,i1,i2,i3)*c
           tt2=x(2,i1,i2,i3)*c
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt1=tt1+x(1,j,i2,i3)*fil(1,l,1)-x(2,j,i2,i3)*fil(2,l,1)
              tt2=tt2+x(2,j,i2,i3)*fil(1,l,1)+x(1,j,i2,i3)*fil(2,l,1)
           enddo
           y(1,i1,i2,i3)=tt1
           y(2,i1,i2,i3)=tt2
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt1=0._wp
           tt2=0._wp
!           do l=lowfil,lupfil
!              j=modulo(i2+l,n2+1)
           do l=max(lowfil,-i2),min(lupfil,n2-i2)
           j=i2+l 
              tt1=tt1+x(1,i1,j,i3)*fil(1,l,2)-x(2,i1,j,i3)*fil(2,l,2)
              tt2=tt2+x(2,i1,j,i3)*fil(1,l,2)+x(1,i1,j,i3)*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
        enddo
     enddo
     
  enddo
!$omp enddo
!$omp do

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt1=0._wp
           tt2=0._wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt1=tt1+x(1,i1,i2,j)*fil(1,l,3)-x(2,i1,i2,j)*fil(2,l,3)
              tt2=tt2+x(2,i1,i2,j)*fil(1,l,3)+x(1,i1,i2,j)*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
        enddo
     enddo
  enddo
!$omp enddo
!$omp end parallel  
END SUBROUTINE convolut_kinetic_slab_c_k


!> Applies the modified kinetic energy operator onto x to get y. 
!! Computes the kinetic energy too.
!! Works for the slab BC.
!! Modified kinetic energy operator:
!! A=-1/2 exp(Ikr) Delta exp(-Ikr)
!! where k=(k1,k2,k3); r=(x,y,z)
subroutine convolut_kinetic_slab_T_k(n1,n2,n3,hgrid,x,y,ener,k1,k2,k3)

  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::k1,k2,k3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(2,0:n1,0:n2,0:n3), intent(inout) :: y
  real(wp), intent(out) :: ener
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt1,tt2,c,enerp
  real(wp), dimension(3) :: scale,scale1
  real(wp), dimension(2,lowfil:lupfil,3) :: fil  

  scale (:)=real(-.5_gp/hgrid(:)**2,wp)
  scale1(1)=real(k1/hgrid(1),wp)
  scale1(2)=real(k2/hgrid(2),wp)
  scale1(3)=real(k3/hgrid(3),wp)
  c=.5_wp*(k1*k1+k2*k2+k3*k3)

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

  ener=0._wp
!$omp parallel default (private) shared(x,y,ener,fil,c,n1,n2,n3)
  enerp=0._wp

!$omp do 

  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt1=x(1,i1,i2,i3)*c
           tt2=x(2,i1,i2,i3)*c
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt1=tt1+x(1,j,i2,i3)*fil(1,l,1)-x(2,j,i2,i3)*fil(2,l,1)
              tt2=tt2+x(2,j,i2,i3)*fil(1,l,1)+x(1,j,i2,i3)*fil(2,l,1)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
         enerp=enerp+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt1=0._wp
           tt2=0._wp
           !do l=lowfil,lupfil
           !   j=modulo(i2+l,n2+1)
           do l=max(lowfil,-i2),min(lupfil,n2-i2)
           j=i2+l 
              tt1=tt1+x(1,i1,j,i3)*fil(1,l,2)-x(2,i1,j,i3)*fil(2,l,2)
              tt2=tt2+x(2,i1,j,i3)*fil(1,l,2)+x(1,i1,j,i3)*fil(2,l,2)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
         enerp=enerp+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
     
  enddo
!$omp enddo

!$omp do 
  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt1=0._wp
           tt2=0._wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt1=tt1+x(1,i1,i2,j)*fil(1,l,3)-x(2,i1,i2,j)*fil(2,l,3)
              tt2=tt2+x(2,i1,i2,j)*fil(1,l,3)+x(1,i1,i2,j)*fil(2,l,3)
           enddo
           y(1,i1,i2,i3)=y(1,i1,i2,i3)+tt1
           y(2,i1,i2,i3)=y(2,i1,i2,i3)+tt2
         enerp=enerp+tt1*x(1,i1,i2,i3)+tt2*x(2,i1,i2,i3)
        enddo
     enddo
  enddo
!$omp enddo
!$omp critical
ener=ener+enerp
!$omp end critical
!  ener=ener*.5_wp
!$omp end parallel  
END SUBROUTINE convolut_kinetic_slab_T_k
