!> @file
!! Simple routines of convolutions
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> A analysis wavelet transformation where the size of the data is forced to shrink
!! The input array y is not overwritten
subroutine analyse_slab(n1,n2,n3,ww,y,x)
  implicit real(kind=8) (a-h,o-z)
  real*8 x(0:n1,2,0:n2,2,0:n3,2)
  real*8 y (0:2*n1+1,-7:2*n2+8,0:2*n3+1)
  real*8 ww(0:2*n1+1,-7:2*n2+8,0:2*n3+1)

  ! i1,I2,i3 -> I2,i3,i1
  nt=(2*n2+16)*(2*n3+2)
  call  ana_rot_per(n1,nt,y,ww)
  ! I2,i3,i1 -> i3,i1,i2
  nt=(2*n3+2)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,y,x)

  return
END SUBROUTINE analyse_slab


!> A synthesis wavelet transformation where the size of the data is allowed to grow
!! The input array x is not overwritten
subroutine synthese_slab(n1,n2,n3,ww,x,y)
  implicit real(kind=8) (a-h,o-z)
  real*8 x(0:n1,2,0:n2,2,0:n3,2)
  real*8 y (0:2*n1+1,-7:2*n2+8,0:2*n3+1)
  real*8 ww(0:2*n1+1,-7:2*n2+8,0:2*n3+1)

  ! i1,i2,i3 -> i2,i3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_per(n1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,I2
  nt=(2*n3+2)*(2*n1+2)
  call  syn_rot_grow(n2,nt,y,ww)
  ! i3,i1,I2  -> i1,I2,i3
  nt=(2*n1+2)*(2*n2+16)
  call  syn_rot_per(n3,nt,ww,y)

END SUBROUTINE synthese_slab


!> A analysis wavelet transformation where the size of the data is forced to shrink
!! The input array y is overwritten
subroutine analyse_slab_self(n1,n2,n3,y,x)
  implicit none
  integer,intent(in)::n1,n2,n3
  real*8,dimension((2*n1+2)*(2*n2+16)*(2*n3+2))::x,y
  integer nt

  ! i1,I2,i3 -> I2,i3,i1
  nt=(2*n2+16)*(2*n3+2)
  call  ana_rot_per(n1,nt,y,x)
  ! I2,i3,i1 -> i3,i1,i2
  nt=(2*n3+2)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,x,y)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,y,x)

END SUBROUTINE analyse_slab_self

subroutine synthese_slab_self(n1,n2,n3,x,y)
  ! A synthesis wavelet transformation where the size of the data is allowed to grow
  ! The input array x is overwritten
  implicit none
  integer,intent(in)::n1,n2,n3
  real*8,dimension((2*n1+2)*(2*n2+16)*(2*n3+2))::x,y
  integer nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_per(n1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,I2
  nt=(2*n3+2)*(2*n1+2)
  call  syn_rot_grow(n2,nt,y,x)
  ! i3,i1,I2  -> i1,I2,i3
  nt=(2*n1+2)*(2*n2+16)
  call  syn_rot_per(n3,nt,x,y)

END SUBROUTINE synthese_slab_self


!> Applies the magic filter matrix in slabwise BC ( no transposition)
!! The input array x is overwritten
!! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_slab_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  !n(c) character(len=*), parameter :: subname='convolut_magic_n_per'
  !n(c) integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat,i_stat,i_all
  !n(c) real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
  !n(c) real(kind=4) filCUDA(lowfil:lupfil) !array of filters to be passed to CUDA interface
  !n(c) data filCUDA / &
  !n(c)     8.4334247333529341094733325815816e-7_4,&
  !n(c)     -0.1290557201342060969516786758559028e-4_4,&
  !n(c)     0.8762984476210559564689161894116397e-4_4,&
  !n(c)     -0.30158038132690463167163703826169879e-3_4,&
  !n(c)     0.174723713672993903449447812749852942e-2_4,&
  !n(c)     -0.942047030201080385922711540948195075e-2_4,&
  !n(c)     0.2373821463724942397566389712597274535e-1_4,&
  !n(c)     0.612625895831207982195380597e-1_4,&
  !n(c)     0.9940415697834003993178616713_4,&
  !n(c)     -0.604895289196983516002834636e-1_4, &
  !n(c)     -0.2103025160930381434955489412839065067e-1_4,&
  !n(c)     0.1337263414854794752733423467013220997e-1_4,&
  !n(c)     -0.344128144493493857280881509686821861e-2_4,&
  !n(c)     0.49443227688689919192282259476750972e-3_4,&
  !n(c)     -0.5185986881173432922848639136911487e-4_4,&
  !n(c)     2.72734492911979659657715313017228e-6_4 /


  if (.not. GPUconv) then !traditional CPU computation

     !  (i1,i2*i3) -> (i2*i3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_n_per(n1,ndat,x,y)
     !  (i2,i3*i1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+1)
     call convrot_grow(n2,ndat,y,x)
     !  (i3,i1*I2) -> (i1*I2,i3)
     ndat=(n1+1)*(n2+16)
     call convrot_n_per(n3,ndat,x,y)

  else
     stop 'the GPU part is not yet written'
  end if
END SUBROUTINE convolut_magic_n_slab_self


!> Applies the magic filter matrix in periodic BC ( no transposition)
!! The input array x is not overwritten
!! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_slab(n1,n2,n3,x,y,ww)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  character(len=*), parameter :: subname='convolut_magic_n_per'
  !n(c) integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat,i_stat,i_all
  real(wp), dimension(0:n1,-7:n2+8,0:n3):: ww ! work array
  !n(c) real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
  !n(c) real(kind=4) filCUDA(lowfil:lupfil) !array of filters to be passed to CUDA interface
  !n(c) data filCUDA / &
  !n(c)     8.4334247333529341094733325815816e-7_4,&
  !n(c)     -0.1290557201342060969516786758559028e-4_4,&
  !n(c)     0.8762984476210559564689161894116397e-4_4,&
  !n(c)     -0.30158038132690463167163703826169879e-3_4,&
  !n(c)     0.174723713672993903449447812749852942e-2_4,&
  !n(c)     -0.942047030201080385922711540948195075e-2_4,&
  !n(c)     0.2373821463724942397566389712597274535e-1_4,&
  !n(c)     0.612625895831207982195380597e-1_4,&
  !n(c)     0.9940415697834003993178616713_4,&
  !n(c)     -0.604895289196983516002834636e-1_4, &
  !n(c)     -0.2103025160930381434955489412839065067e-1_4,&
  !n(c)     0.1337263414854794752733423467013220997e-1_4,&
  !n(c)     -0.344128144493493857280881509686821861e-2_4,&
  !n(c)     0.49443227688689919192282259476750972e-3_4,&
  !n(c)     -0.5185986881173432922848639136911487e-4_4,&
  !n(c)     2.72734492911979659657715313017228e-6_4 /


  if (.not. GPUconv) then !traditional CPU computation

     !  (i1,i2*i3) -> (i2*i3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_n_per(n1,ndat,x,y)
     !  (i2,i3*i1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+1)
     call convrot_grow(n2,ndat,y,ww)
     !  (i3,i1*I2) -> (i1*I2,i3)
     ndat=(n1+1)*(n2+16)
     call convrot_n_per(n3,ndat,ww,y)

  else
      stop 'the GPU part is not yet written'
  end if
END SUBROUTINE convolut_magic_n_slab


!> Applies the magic filter matrix transposed in periodic BC 
!! The input array x is overwritten
!! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_t_slab_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  !n(c) character(len=*), parameter :: subname='convolut_magic_t_per'
  integer, parameter :: lowfil=-7,lupfil=8
  integer :: ndat,i_stat,i_all
  !n(c) real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
  !n(c) real(kind=4) filCUDA(lowfil:lupfil) !array of filters to be passed to CUDA interface
  !n(c) data filCUDA / &
  !n(c)     2.72734492911979659657715313017228e-6_4,&
  !n(c)     -0.5185986881173432922848639136911487e-4_4,&
  !n(c)     0.49443227688689919192282259476750972e-3_4,&
  !n(c)     -0.344128144493493857280881509686821861e-2_4,&
  !n(c)     0.1337263414854794752733423467013220997e-1_4,&
  !n(c)     -0.2103025160930381434955489412839065067e-1_4,&
  !n(c)     -0.604895289196983516002834636e-1_4,&
  !n(c)     0.9940415697834003993178616713_4,&
  !n(c)     0.612625895831207982195380597e-1_4,&
  !n(c)     0.2373821463724942397566389712597274535e-1_4,&
  !n(c)     -0.942047030201080385922711540948195075e-2_4,&
  !n(c)     0.174723713672993903449447812749852942e-2_4,&
  !n(c)     -0.30158038132690463167163703826169879e-3_4,&
  !n(c)     0.8762984476210559564689161894116397e-4_4,&
  !n(c)     -0.1290557201342060969516786758559028e-4_4,&
  !n(c)     8.4334247333529341094733325815816e-7_4 /

  
  if (.not. GPUconv) then

     !  (i1,I2*i3) -> (I2*i3,i1)
     ndat=(n2+16)*(n3+1)
     call convrot_t_per(n1,ndat,x,y)
     !  (I2,i3*i1) -> (i3*i1,i2)
     ndat=(n3+1)*(n1+1)
     call convrot_shrink(n2,ndat,y,x)
     !  (i3,i1*i2) -> (i1*i2,i3)
     ndat=(n1+1)*(n2+1)
     call convrot_t_per(n3,ndat,x,y)

  else
      stop 'the GPU part is not yet written'
  end if

END SUBROUTINE convolut_magic_t_slab_self


!>  Applies the kinetic energy operator onto x to get y. Works for periodic BC
subroutine convolut_kinetic_slab_c(n1,n2,n3,hgrid,x,y,c)
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
!          do l=lowfil,lupfil
!             j=modulo(i2+l,n2+1)
!             tt=tt+x(i1,j   ,i3)*fil(l,2)
           do l=max(lowfil,-i2),min(lupfil,n2-i2)
              tt=tt+x(i1,i2+l ,i3)*fil(l,2)
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
  
END SUBROUTINE convolut_kinetic_slab_c


!>  Applies the kinetic energy operator onto x to get y. Works for periodic BC
subroutine convolut_kinetic_slab_T(n1,n2,n3,hgrid,x,y,ekin)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
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
!          do l=lowfil,lupfil
!              j=modulo(i2+l,n2+1)
!              tt=tt+x(i1,j   ,i3)*fil(l,2)
           do l=max(lowfil,-i2),min(lupfil,n2-i2)
              tt=tt+x(i1,i2+l ,i3)*fil(l,2)
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
  
END SUBROUTINE convolut_kinetic_slab_T
