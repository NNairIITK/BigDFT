!> @file
!!  Routine usig kernels
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> hits the input array x with the kernel
!! @f$ ((-1/2\Delta+C)_{ij})^{-1} @f$
subroutine hit_with_kernel(x,z,kern_k1,kern_k2,kern_k3,n1,n2,n3,nd1,nd2,nd3,c)  
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,nd1,nd2,nd3
  real(gp),intent(in)::kern_k1(0:n1)
  real(gp),intent(in)::kern_k2(0:n2)
  real(gp),intent(in)::kern_k3(0:n3)
  real(gp),intent(in)::c

  real(wp),intent(inout)::x(0:n1,0:n2,0:n3)! input/output

  real(wp)::z(2,0:(n1+1)/2,0:n2,0:n3)! work array
  real(gp) tt,fac_n
  integer i1,i2,i3,isign,inzee
!***for fft:************************************************************************************
  include 'fftw3.f'
  integer(8)::plan_f,plan_b
!***********************************************************************************************
  
  fac_n=(n1+1)*(n2+1)*(n3+1)
  call dfftw_plan_dft_r2c_3d(plan_f,n1+1,n2+1,n3+1,x,z,fftw_estimate)
     call dfftw_execute(plan_f)
  call dfftw_destroy_plan(plan_f)
  ! hit the fourier transform of x with the kernel
  do i3=0,n3
    do i2=0,n2
      do i1=0,(n1+1)/2
        tt=(kern_k1(i1)+kern_k2(i2)+kern_k3(i3)+c)*fac_n
        z(:,i1,i2,i3)=z(:,i1,i2,i3)/tt
      enddo
    enddo
  enddo

  call dfftw_plan_dft_c2r_3d(plan_b,n1+1,n2+1,n3+1,z,x,fftw_estimate)
     call dfftw_execute(plan_b)
  call dfftw_destroy_plan(plan_b)

END SUBROUTINE hit_with_kernel


!> construct the kernel (-1/2 d^2/dx^2)_{ij}
!! at a real space grid with grid size hgrid
!! and then fourier transform it to momentum space
subroutine make_kernel(n1,hgrid,kern)
use module_base
implicit none
integer,intent(in)::n1
real(gp),intent(in)::hgrid
real(gp),intent(out)::kern(0:n1)

include 'fftw3.f'
integer*8::plan_f
real(gp),allocatable::z(:,:),zin(:,:)
real(gp) dmax,tt
integer i

integer,parameter::lowfil=-14,lupfil=14
real(gp) scale
real(gp)::fil(lowfil:lupfil)

scale=-.5_gp/hgrid**2

! second derivative filters for Daubechies 16
fil(0)=   -3.5536922899131901941296809374_gp*scale
fil(1)=    2.2191465938911163898794546405_gp*scale
fil(2)=   -0.6156141465570069496314853949_gp*scale
fil(3)=    0.2371780582153805636239247476_gp*scale
fil(4)=   -0.0822663999742123340987663521_gp*scale
fil(5)=    0.02207029188482255523789911295638968409_gp*scale
fil(6)=   -0.409765689342633823899327051188315485e-2_gp*scale
fil(7)=    0.45167920287502235349480037639758496e-3_gp*scale
fil(8)=   -0.2398228524507599670405555359023135e-4_gp*scale
fil(9)=    2.0904234952920365957922889447361e-6_gp*scale
fil(10)=  -3.7230763047369275848791496973044e-7_gp*scale
fil(11)=  -1.05857055496741470373494132287e-8_gp*scale
fil(12)=  -5.813879830282540547959250667e-11_gp*scale
fil(13)=   2.70800493626319438269856689037647576e-13_gp*scale
fil(14)=  -6.924474940639200152025730585882e-18_gp*scale

do i=1,14
  fil(-i)=fil(i)
enddo

allocate(z(2,0:n1))
allocate(zin(2,0:n1))

! construct the kernel in real space
zin=0._gp
zin(1,0)=fil(0)
do i=1,14
  zin(1,i)     =fil(i)
  zin(1,n1+1-i)=fil(i)
enddo

call dfftw_plan_dft_1d(plan_f,n1+1,zin,z,FFTW_FORWARD ,FFTW_ESTIMATE)
call dfftw_execute(plan_f)
call dfftw_destroy_plan(plan_f)

do i=0,n1
  kern(i)=z(1,i)
enddo

deallocate(z,zin)
END SUBROUTINE make_kernel
