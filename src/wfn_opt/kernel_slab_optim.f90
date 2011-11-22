!> @file
!!   Optimized routines using kernels
!! @author
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Hits the input array x with the kernel
!! @f$((-1/2\Delta+C)_{ij})^{-1}@f$
subroutine hit_with_kernel_slab(x,zx,kern_k1,kern_k3,n1,n2,n3,c,hgrid)
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3
  real(gp),intent(in)::kern_k1(0:n1)
  real(gp),intent(in)::kern_k3(0:n3)
  real(gp),intent(in)::c
  real(gp),intent(in)::hgrid

  real(wp),intent(inout)::x(0:n1,0:n2,0:n3)! input and output
  real(wp)::zx(2,0:(n1+1)/2,0:n2,0:n3)! work array
  integer ::  i1,i2,i3,isign,inzee,i,j,l
  integer :: nd1,nd2,nd3,ntrig

  !*******************************************************************************************
  ! for 2-dimensional FFT:
  real(wp),allocatable::trig(:,:),z(:,:,:,:),zw(:)
  integer,parameter::ncache=4*1024

  real(wp) t1,t2,tela
  integer count1,count2,count_rate,count_max

  nd1=n1+2
  nd2=n2+2
  nd3=n3+2

  ! fourier transform the x


  !allocations to be reformulated
  ntrig=max(n1,n2,n3)
  allocate(trig(2,ntrig))

  call cpu_time(t1)
  call system_clock(count1,count_rate,count_max)      

  call forward_fft(n1,n2,n3,nd1,nd3,x,zx,ntrig,trig)

  call segment_invert(n1,n2,n3,kern_k1,kern_k3,c,zx,hgrid)

  call backward_fft(n1,n2,n3,nd1,nd3,x,zx,ntrig,trig)

  call cpu_time(t2)
  call system_clock(count2,count_rate,count_max)      
  tela=(count2-count1)/real(count_rate)
  !write(*,*) 'Time (CPU,ELA)  (sec):' ,t2-t1,tela

  !deallocations to be reformulated
  deallocate(trig)
END SUBROUTINE hit_with_kernel_slab


subroutine forward_fft(n1,n2,n3,nd1,nd3,x,zx,ntrig,trig)
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,nd1,nd3,ntrig
  real(wp),intent(in)::x(0:n1,0:n2,0:n3)
  real (wp),intent(in):: trig(2,ntrig)
  real(wp),intent(inout)::zx(2,0:(n1+1)/2,0:n2,0:n3)! work array
  integer,parameter::ncache=4*1024
  real(wp),allocatable::zw(:)
  real(wp),allocatable::z(:,:,:,:)
  integer::i1,i2,i3,isign,inzee,i,j,l
!$omp parallel default(private) shared(n2,n3,n1,x,nd1,nd3,zx,zw)
  allocate(z(2,nd1,nd3,2))
  allocate(zw(ncache+1))
  isign=1
  inzee=1
  ! Fourier transform along x:
!$omp do
    do i2=0,n2
       do i3=0,n3
          do i1=0,n1
             z(1,i1+1,i3+1,inzee)=x(i1,i2,i3)
             z(2,i1+1,i3+1,inzee)=0.0_wp
          enddo
       enddo
       call fft2d(n1+1,n3+1,nd1,nd3,z,isign,inzee,zw,ncache)
       do i3=0,n3
          do i1=0,(n1+1)/2
             zx(1,i1,i2,i3)=z(1,i1+1,i3+1,inzee)
             zx(2,i1,i2,i3)=z(2,i1+1,i3+1,inzee)
          enddo
       enddo
    enddo
!$omp enddo
  deallocate(z,zw)
!$omp end parallel

END SUBROUTINE forward_fft


subroutine backward_fft(n1,n2,n3,nd1,nd3,x,zx,ntrig,trig)
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,nd1,nd3,ntrig
  real(wp),intent(inout)::x(0:n1,0:n2,0:n3)! input and output
  real(wp),intent(in)::zx(2,0:(n1+1)/2,0:n2,0:n3)! work array
  real (wp),intent(in):: trig(2,ntrig)
  integer,parameter::ncache=4*1024
  real(wp),allocatable::z(:,:,:,:)
  real(wp),allocatable::zw(:)
  integer ::i1,i2,i3,isign,inzee,i,j,l

!$omp parallel default(private) shared(n2,n3,n1,x,nd1,nd3,zx,zw)
  allocate(z(2,nd1,nd3,2))
  allocate(zw(ncache+1))

  isign=-1
  inzee=1
!$omp do
  do i2=0,n2
     ! i3=0
     ! zx(i1,i2,0)=zx*(n1+1-i1,i2,0) for i1 != 0
     z(1,1,1,inzee)=zx(1,0,i2,0)
     z(2,1,1,inzee)=zx(2,0,i2,0)
     do i1=1,(n1+1)/2
        z(1,i1+1,1,inzee)=zx(1,i1,i2,0)
        z(2,i1+1,1,inzee)=zx(2,i1,i2,0)

        z(1,n1+2-i1,1,inzee)=zx(1,i1,i2,0)
        z(2,n1+2-i1,1,inzee)=-zx(2,i1,i2,0)
     enddo

     ! nonzero i3
     ! zx(i1,i2,i3)=zx*(n1+1-i1,i2,n3+1-i3) for i1 != 0
     do i3=1,n3

        z(1,1,i3+1,inzee)=zx(1,0,i2,i3)
        z(2,1,i3+1,inzee)=zx(2,0,i2,i3)
        do i1=1,(n1+1)/2
           z(1,i1+1,i3+1,inzee)=zx(1,i1,i2,i3)
           z(2,i1+1,i3+1,inzee)=zx(2,i1,i2,i3)

           z(1,n1+2-i1,n3+2-i3,inzee)=zx(1,i1,i2,i3)
           z(2,n1+2-i1,n3+2-i3,inzee)=-zx(2,i1,i2,i3)
        enddo

     enddo

     call fft2d(n1+1,n3+1,nd1,nd3,z,isign,inzee,zw,ncache)
     do i3=0,n3
        do i1=0,n1
           x(i1,i2,i3)=z(1,i1+1,i3+1,inzee)
        enddo
     enddo
  enddo
!$omp enddo
  deallocate(z,zw)
!$omp end parallel
END SUBROUTINE backward_fft
