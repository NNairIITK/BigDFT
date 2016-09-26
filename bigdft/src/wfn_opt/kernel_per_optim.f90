!> @file
!!  OPtimzed routines using kernels
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Hit the Fourier transform of x with the kernel. At the same time, transform the array
!!  from the form z3 (where only half of values of i3 are stored)
!!  to the form z1   (where only half of values of i1 are stored)
!!  The latter thing could be done separately by the subroutine z3_to_z1 that is contained
!!  in FFT_back, but then the code would be slower.
subroutine hit_with_kernel_fac(x,z1,z3,kern_k1,kern_k2,kern_k3,n1,n2,n3,nd1,nd2,nd3,&
     n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,c,fac)
  ! hits the input array x with the kernel
  ! @f$ ((-1/2\Delta+C)_{ij})^{-1} @f$
  use module_base
  implicit none
  integer,intent(in) :: n1,n2,n3,nd1,nd2,nd3
  integer,intent(in) :: n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b
  real(gp),intent(in) :: kern_k1(n1)
  real(gp),intent(in) :: kern_k2(n2)
  real(gp),intent(in) :: kern_k3(n3)
  real(gp),intent(in) :: c,fac

  real(wp),intent(inout)::x(n1,n2,n3)! input/output

  real(wp) :: z1(2,nd1b,nd2,nd3,2)! work array
  real(wp) :: z3(2,nd1,nd2,nd3f,2)! work array
  real(gp) :: tt
  integer :: i1,i2,i3,inzee

  ! fft the input array x:

  call FFT_for(n1,n2,n3,n1f,n3f,nd1,nd2,nd3,nd1f,nd3f,x,z1,z3,inzee)
  z1=0.d0 !to be put in the loop below
  !$omp parallel default (private) shared(z1,z3,kern_k1,kern_k2,kern_k3,c,fac)&
  !$omp shared(n1b,n3f,inzee,n1,n2,n3)

  ! i3=1: then z1 is contained in z3 
  !$omp do 
  do i2=1,n2
     do i1=1,n1b
        tt=fac/(kern_k1(i1)+kern_k2(i2)+kern_k3(1)+c)
        z1(1,i1,i2,1,inzee)=z3(1,i1,i2,1,inzee)*tt
        z1(2,i1,i2,1,inzee)=z3(2,i1,i2,1,inzee)*tt
     enddo
  enddo
  !$omp enddo

  !$omp do
  do i3=2,n3f
     ! i2=1
     ! i1=1
     tt=fac/(kern_k1(1)+kern_k2(1)+kern_k3(i3)+c)
     z1(1,1,1,i3,inzee)=z3(1,1,1,i3,inzee)*tt
     z1(2,1,1,i3,inzee)=z3(2,1,1,i3,inzee)*tt

     z1(1,1,1,n3+2-i3,inzee)=z3(1,1,1,i3,inzee)*tt
     z1(2,1,1,n3+2-i3,inzee)=-z3(2,1,1,i3,inzee)*tt

     ! i2=1
     do i1=2,n1b
        tt=fac/(kern_k1(i1)+kern_k2(1)+kern_k3(i3)+c)
        z1(1,i1,1,i3,inzee)=z3(1,i1,1,i3,inzee)*tt
        z1(2,i1,1,i3,inzee)=z3(2,i1,1,i3,inzee)*tt

        z1(1,i1,1,n3+2-i3,inzee)= z3(1,n1+2-i1,1,i3,inzee)*tt
        z1(2,i1,1,n3+2-i3,inzee)=-z3(2,n1+2-i1,1,i3,inzee)*tt
     enddo

     do i2=2,n2
        ! i1=1
        tt=fac/(kern_k1(1)+kern_k2(i2)+kern_k3(i3)+c)
        z1(1,1,i2,i3,inzee)=z3(1,1,i2,i3,inzee)*tt
        z1(2,1,i2,i3,inzee)=z3(2,1,i2,i3,inzee)*tt

        z1(1,1,i2,n3+2-i3,inzee)= z3(1,1,n2+2-i2,i3,inzee)*tt
        z1(2,1,i2,n3+2-i3,inzee)=-z3(2,1,n2+2-i2,i3,inzee)*tt

        do i1=2,n1b
           tt=fac/(kern_k1(i1)+kern_k2(i2)+kern_k3(i3)+c)
           z1(1,i1,i2,i3,inzee)=z3(1,i1,i2,i3,inzee)*tt
           z1(2,i1,i2,i3,inzee)=z3(2,i1,i2,i3,inzee)*tt

           z1(1,i1,i2,n3+2-i3,inzee)= z3(1,n1+2-i1,n2+2-i2,i3,inzee)*tt
           z1(2,i1,i2,n3+2-i3,inzee)=-z3(2,n1+2-i1,n2+2-i2,i3,inzee)*tt
        enddo
     enddo
  enddo
  !$omp enddo

  !$omp end parallel

  call FFT_back(n1,n2,n3,n1b,n3f,n3b,nd1,nd2,nd3,nd1b,nd3f,nd3b,x,z1,z3,inzee)

END SUBROUTINE hit_with_kernel_fac


!>  Hits the input array x with the kernel @f$((-1/2\Delta+C)_{ij})^{-1}@f$
subroutine hit_with_kernel(x,z1,z3,kern_k1,kern_k2,kern_k3,n1,n2,n3,nd1,nd2,nd3,&
  n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,c)
  use module_base
  implicit none
! Arguments
  integer,intent(in) :: n1,n2,n3,nd1,nd2,nd3
  integer,intent(in) :: n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b
  real(gp),intent(in) :: kern_k1(n1)
  real(gp),intent(in) :: kern_k2(n2)
  real(gp),intent(in) :: kern_k3(n3)
  real(gp),intent(in) :: c
  real(wp),intent(inout) :: x(n1,n2,n3)! input/output
!Local variables
  real(wp) :: z1(2,nd1b,nd2,nd3,2)! work array
  real(wp) :: z3(2,nd1,nd2,nd3f,2)! work array
  real(gp) :: tt
  integer :: i1,i2,i3,inzee

! fft the input array x:
  !call to_zero(2*nd1b*nd2*nd3*2,z1(1,1,1,1,1))
  !call to_zero(2*nd1*nd2*nd3f*2,z3(1,1,1,1,1))
  call FFT_for(n1,n2,n3,n1f,n3f,nd1,nd2,nd3,nd1f,nd3f,x,z1,z3,inzee)
  call f_zero(z1) !to be included in the loop below
! hit the Fourier transform of x with the kernel. At the same time, transform the array
! from the form z3 (where only half of values of i3 are stored)
! to the form z1   (where only half of values of i1 are stored)
! The latter thing could be done separately by the subroutine z3_to_z1 that is contained
! in FFT_back, but then the code would be slower.

  !$omp parallel default (private) shared(z1,z3,kern_k1,kern_k2,kern_k3,c)&
  !$omp shared(n1b,n3f,inzee,n1,n2,n3)

  ! i3=1: then z1 is contained in z3 


  !$omp do 
  do i2=1,n2
    do i1=1,n1b
      tt=1._gp/(kern_k1(i1)+kern_k2(i2)+kern_k3(1)+c)
      z1(1,i1,i2,1,inzee)=z3(1,i1,i2,1,inzee)*tt
      z1(2,i1,i2,1,inzee)=z3(2,i1,i2,1,inzee)*tt
    enddo
  enddo  
  !$omp enddo

  !$omp do
  do i3=2,n3f
    ! i2=1
    ! i1=1
    tt=1._gp/(kern_k1(1)+kern_k2(1)+kern_k3(i3)+c)
    z1(1,1,1,i3,inzee)=z3(1,1,1,i3,inzee)*tt
    z1(2,1,1,i3,inzee)=z3(2,1,1,i3,inzee)*tt

    z1(1,1,1,n3+2-i3,inzee)=z3(1,1,1,i3,inzee)*tt
    z1(2,1,1,n3+2-i3,inzee)=-z3(2,1,1,i3,inzee)*tt

    ! i2=1
    do i1=2,n1b  
      tt=1._gp/(kern_k1(i1)+kern_k2(1)+kern_k3(i3)+c)
      z1(1,i1,1,i3,inzee)=z3(1,i1,1,i3,inzee)*tt
      z1(2,i1,1,i3,inzee)=z3(2,i1,1,i3,inzee)*tt

      z1(1,i1,1,n3+2-i3,inzee)= z3(1,n1+2-i1,1,i3,inzee)*tt
      z1(2,i1,1,n3+2-i3,inzee)=-z3(2,n1+2-i1,1,i3,inzee)*tt
    enddo

    do i2=2,n2
      ! i1=1
      tt=1._gp/(kern_k1(1)+kern_k2(i2)+kern_k3(i3)+c)
      z1(1,1,i2,i3,inzee)=z3(1,1,i2,i3,inzee)*tt
      z1(2,1,i2,i3,inzee)=z3(2,1,i2,i3,inzee)*tt

      z1(1,1,i2,n3+2-i3,inzee)= z3(1,1,n2+2-i2,i3,inzee)*tt
      z1(2,1,i2,n3+2-i3,inzee)=-z3(2,1,n2+2-i2,i3,inzee)*tt

      do i1=2,n1b
        tt=1._gp/(kern_k1(i1)+kern_k2(i2)+kern_k3(i3)+c)
        z1(1,i1,i2,i3,inzee)=z3(1,i1,i2,i3,inzee)*tt
        z1(2,i1,i2,i3,inzee)=z3(2,i1,i2,i3,inzee)*tt

        z1(1,i1,i2,n3+2-i3,inzee)= z3(1,n1+2-i1,n2+2-i2,i3,inzee)*tt
        z1(2,i1,i2,n3+2-i3,inzee)=-z3(2,n1+2-i1,n2+2-i2,i3,inzee)*tt
      enddo
    enddo
  enddo
  !$omp enddo

  !$omp end parallel

  call FFT_back(n1,n2,n3,n1b,n3f,n3b,nd1,nd2,nd3,nd1b,nd3f,nd3b,x,z1,z3,inzee)

END SUBROUTINE hit_with_kernel


!> Construct the kernel (-1/2 d^2/dx^2)_{ij}
!! at a real space grid with grid size hgrid
!! and then fourier transform it to momentum space
subroutine make_kernel(n1,hgrid,kern)
  use module_fft_sg
  use module_base
  implicit none
  integer,intent(in)::n1
  real(gp),intent(in)::hgrid

  real(gp),intent(out)::kern(0:n1)

!***********************************************************************************************
  integer :: now(7),after(7),before(7)
  integer :: isign=1
  real(gp),allocatable :: trig(:,:),z(:,:,:)
  integer :: inzee,i,nd1,ic
!***********************************************************************************************


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

! construct the kernel in real space

!***********************************************************************************************
! fourier transform the kernel

  nd1=n1+2
  allocate(trig(2,nd1))
  allocate(z(2,nd1,2))

  inzee=1
  z=0.0_gp
  z(1,1,inzee)=fil(0)

  if (nd1 < 14) then
     call f_err_throw('ERROR: dimension too little dimension n=' // trim(yaml_toa(nd1)) // ' n1=' // trim(yaml_toa(n1)), &
          &         err_name='BIGDFT_INPUT_VARIABLES_ERROR')
     !write(*,*)'ERROR: dimension too little dimension n',nd1,n1
     !stop
  end if

  do i=1,14
    z(1,i+1,inzee)=fil(i)
  enddo

  do i=1,14
    z(1,n1+2-i,inzee)=z(1,n1+2-i,inzee)+fil(i)
  enddo

  isign=1
  call ctrig_sg(n1+1,nd1,trig,after,before,now,isign,ic)
  do i=1,ic
    call fftstp_sg(1,1,nd1,1,nd1,z(1,1,inzee),z(1,1,3-inzee),nd1,trig,after(i),now(i),before(i),isign)
    inzee=3-inzee
  enddo

  do i=0,n1
    kern(i)=z(1,i+1,inzee)
  enddo

  deallocate(trig,z)
END SUBROUTINE make_kernel
