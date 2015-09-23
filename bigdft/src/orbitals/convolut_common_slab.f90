!> @file
!!  Common convolutions
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>   A analysis wavelet transformation where the size of the data is forced to shrink
!!   The input array y is overwritten
subroutine analyse_shrink(n1,n2,n3,ww,y,x)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8), intent(inout) :: ww
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(inout) :: y
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: x
  !local variables
  integer :: nt

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+16)*(2*n3+16)
  call  ana_rot_shrink(n1,nt,y,ww)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+16)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_shrink(n3,nt,y,x)

END SUBROUTINE analyse_shrink


!>   A synthesis wavelet transformation where the size of the data is allowed to grow
!!   The input array x is not overwritten
subroutine synthese_grow(n1,n2,n3,ww,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(in) :: x
  real(wp), dimension(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8), intent(inout) :: ww
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,I1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,ww)
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_grow(n3,nt,ww,y)

END SUBROUTINE synthese_grow


!>   An analysis wavelet transformation where the size of the data is forced to shrink
!!   The input array y is not overwritten
subroutine analyse_slab(n1,n2,n3,ww,y,x)
  implicit none
  !Arguments
  integer, intent(in) :: n1,n2,n3
  real(kind=8) :: x(0:n1,2,0:n2,2,0:n3,2)
  real(kind=8) :: y (0:2*n1+1,-7:2*n2+8,0:2*n3+1)
  real(kind=8) :: ww(0:2*n1+1,-7:2*n2+8,0:2*n3+1)
  !Local variables
  integer :: nt

  ! i1,I2,i3 -> I2,i3,i1
  nt=(2*n2+16)*(2*n3+2)
  call  ana_rot_per(n1,nt,y,ww)
  ! I2,i3,i1 -> i3,i1,i2
  nt=(2*n3+2)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,y,x)

END SUBROUTINE analyse_slab


!>   A synthesis wavelet transformation where the size of the data is allowed to grow
!!   The input array x is not overwritten
subroutine synthese_slab(n1,n2,n3,ww,x,y)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: x(0:n1,2,0:n2,2,0:n3,2)
  real(kind=8) :: y (0:2*n1+1,-7:2*n2+8,0:2*n3+1)
  real(kind=8) :: ww(0:2*n1+1,-7:2*n2+8,0:2*n3+1)

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


!>  A analysis wavelet transformation where the size of the data is forced to shrink
!!  The input array y is overwritten
subroutine analyse_slab_self(n1,n2,n3,y,x)
  implicit none
  integer,intent(in)::n1,n2,n3
  real(kind=8),dimension((2*n1+2)*(2*n2+16)*(2*n3+2))::x,y
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


!>   A synthesis wavelet transformation where the size of the data is allowed to grow
!!   The input array x is overwritten
subroutine synthese_slab_self(n1,n2,n3,x,y)
  implicit none
  integer,intent(in)::n1,n2,n3
  real(kind=8),dimension((2*n1+2)*(2*n2+16)*(2*n3+2))::x,y
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


!>   Applies the magic filter matrix in slabwise BC ( no transposition)
!!   The input array x is overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_slab_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  !n(c) character(len=*), parameter :: subname='convolut_magic_n_per'
  !n(c) integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat
  !real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
  !n(c) real(kind=4) filCUDA(lowfil:lupfil) !array of filters to be passed to CUDA interface
  !n(c) data filCUDA / &
  !n(c)      8.4334247333529341094733325815816e-7_4,&
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

     !  (i1,i2*i3) -> (i2*i3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_n_per(n1,ndat,x,y)
     !  (i2,i3*i1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+1)
     call convrot_grow(n2,ndat,y,x)
     !  (i3,i1*I2) -> (i1*I2,i3)
     ndat=(n1+1)*(n2+16)
     call convrot_n_per(n3,ndat,x,y)

END SUBROUTINE convolut_magic_n_slab_self


!>   A synthesis wavelet transformation where the size of the data is allowed to grow
!!   The input array x is overwritten
subroutine synthese_free_self(n1,n2,n3,x,y)
  implicit none
  integer,intent(in)::n1,n2,n3
  real(kind=8),dimension((2*n1+16)*(2*n2+16)*(2*n3+16))::x,y
  integer nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,x)
  ! i3,i1,I2  -> i1,I2,i3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_grow(n3,nt,x,y)

END SUBROUTINE synthese_free_self


!>   Applies the magic filter matrix in slabwise BC ( no transposition)
!!   The input array x is overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_free_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:n1+8,-7:n2+8,-7:n3+8), intent(inout) :: x
  real(wp), dimension(-7:n1+8,-7:n2+8,-7:n3+8), intent(inout) :: y
  !local variables
  !n(c) character(len=*), parameter :: subname='convolut_magic_n_per'
  !n(c) integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat
  !real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
  !n(c) real(kind=4) filCUDA(lowfil:lupfil) !array of filters to be passed to CUDA interface
  !n(c) data filCUDA / &
  !n(c)      8.4334247333529341094733325815816e-7_4,&
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

     !  (i1,i2*i3) -> (i2*i3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_grow(n1,ndat,x,y)
     !  (i2,i3*i1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+16)
     call convrot_grow(n2,ndat,y,x)
     !  (i3,i1*I2) -> (i1*I2,i3)
     ndat=(n1+16)*(n2+16)
     call convrot_grow(n3,ndat,x,y)
END SUBROUTINE convolut_magic_n_free_self


!>   Applies the magic filter matrix in periodic BC ( no transposition)
!!   The input array x is not overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_slab(n1,n2,n3,x,y,ww)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: y
  !local variables
  !n(c) character(len=*), parameter :: subname='convolut_magic_n_per'
  !n(c) integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat
  real(wp), dimension(0:n1,-7:n2+8,0:n3):: ww ! work array
  !real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
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


     !  (i1,i2*i3) -> (i2*i3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_n_per(n1,ndat,x,y)
     !  (i2,i3*i1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+1)
     call convrot_grow(n2,ndat,y,ww)
     !  (i3,i1*I2) -> (i1*I2,i3)
     ndat=(n1+1)*(n2+16)
     call convrot_n_per(n3,ndat,ww,y)
END SUBROUTINE convolut_magic_n_slab


!>   Applies the magic filter matrix transposed in periodic BC 
!!   The input array x is overwritten
!!   this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_t_slab_self(n1,n2,n3,x,y)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,-7:n2+8,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  !n(c) character(len=*), parameter :: subname='convolut_magic_t_per'
  !n(c) integer, parameter :: lowfil=-7,lupfil=8
  integer :: ndat
  !real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
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

     !  (i1,I2*i3) -> (I2*i3,i1)
     ndat=(n2+16)*(n3+1)
     call convrot_t_per(n1,ndat,x,y)
     !  (I2,i3*i1) -> (i3*i1,i2)
     ndat=(n3+1)*(n1+1)
     call convrot_shrink(n2,ndat,y,x)
     !  (i3,i1*i2) -> (i1*i2,i3)
     ndat=(n1+1)*(n2+1)
     call convrot_t_per(n3,ndat,x,y)


END SUBROUTINE convolut_magic_t_slab_self


!>   Applies the kinetic energy operator onto x to get y. Works for periodic BC
subroutine convolut_kinetic_slab_sdc(n1,n2,n3,x,y,cprecr,modul1,modul3,a,b,c,e)
  use module_defs, only: wp,gp
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in) :: cprecr
  real(wp), dimension(8,0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(8,0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,l,j
  real(wp) :: tt111,tt112,tt121,tt122,tt211,tt212,tt221,tt222
  integer,intent(in) :: modul1(lowfil:n1+lupfil)
  integer,intent(in) :: modul3(lowfil:n3+lupfil)
  real(gp),intent(in) :: a(lowfil:lupfil,3)
  real(gp),intent(in) :: b(lowfil:lupfil,3)
  real(gp),intent(in) :: c(lowfil:lupfil,3)
  real(gp),intent(in) :: e(lowfil:lupfil,3)
  
!  real(gp) :: tel
!  integer :: ncount0,ncount1,ncount2,ncount_rate,ncount_max
!  integer :: mflop1,mflop3
  
  ! filter length:29
  ! 8: wavelets+scfunction
  ! 4: flops for one iteration
!  mflop1=(n1+1)*(n2+1)*(n3+1)*29*8*4*2 ! convolution in the x and y direction
!  mflop3=(n1+1)*(n2+1)*(n3+1)*29*8*4   ! convolution in the z       direction
!dee
!  open(unit=97,file='time_check',status='unknown')
!  call system_clock(ncount0,ncount_rate,ncount_max)
  
!dee
!$omp parallel default(private) &
!$omp shared (n1,n2,n3,x,y,cprecr,a,b,c,e,modul1,modul3) 


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
           do l=max(lowfil,-i2),min(lupfil,n2-i2) ! so that 0 =< j =< n2
              j=i2+l
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
!$omp enddo
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
!$omp enddo 
!$omp end parallel
!  call system_clock(ncount2,ncount_rate,ncount_max)
!  tel=dble(ncount2-ncount0)/dble(ncount_rate)
!  write(97,'(a40,1x,e10.3,1x,f6.1)') 'z:',tel,1.d-6*mflop3/tel
!  close(97)

END SUBROUTINE convolut_kinetic_slab_sdc


!> BigDFT/prepare_sdc_slab
subroutine prepare_sdc_slab(n1,n3,modul1,modul3,a,b,c,e,hx,hy,hz)
  use module_defs, only: wp,gp
  implicit none
  integer,intent(in)::n1,n3
  real(gp),intent(in)::hx,hy,hz
  
  integer, parameter :: lowfil=-14,lupfil=14
  
  integer,intent(out)::modul1(lowfil:n1+lupfil)
  integer,intent(out)::modul3(lowfil:n3+lupfil)
  real(gp),intent(out)::a(lowfil:lupfil,3)
  real(gp),intent(out)::b(lowfil:lupfil,3)
  real(gp),intent(out)::c(lowfil:lupfil,3)
  real(gp),intent(out)::e(lowfil:lupfil,3)
  
  real(gp)::hgrid(3)
  integer::i
  real(gp)::scale(3)

  call fill_mod_arr(modul1,lowfil,n1+lupfil,n1+1)
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
END SUBROUTINE prepare_sdc_slab
