!> @file
!!   Routines using kernels with FFTW
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> hits the input array x with the kernel
!! @f$((-1/2\Delta+C)_{ij})^{-1}@f$
subroutine hit_with_kernel_slab(x,z  ,kern_k1,kern_k3,n1,n2,n3,c,hgrid)    
    implicit none
    integer,intent(in)::n1,n2,n3
    real*8,intent(in)::kern_k1(0:n1)
    real*8,intent(in)::kern_k3(0:n3)
    real*8,intent(in)::c,hgrid

    real*8,intent(inout)::x(0:n1,0:n2,0:n3)! input and output

    real*8::z(2,0:(n1+1)/2,0:n2,0:n3)! work array 

    integer i1,i2,i3
!***for fft:************************************************************************************
    include 'fftw3.f'
    integer(8)::plan_f,plan_b
!***********************************************************************************************
    real*8,allocatable:: x2(:,:),y2(:,:,:)
    real*8 tt,t0,t1,t2,t3


    allocate(x2(0:n1,0:n3),y2(2,0:(n1+1)/2,0:n3))

! fft the input array x:

    call cpu_time(t0)

    call dfftw_plan_dft_r2c_2d(plan_f,n1+1,n3+1,x2,y2,fftw_estimate)

    do i2=0,n2
        do i3=0,n3
            do i1=0,n1
                x2(i1,i3)=x(i1,i2,i3)
            enddo
        enddo

        call dfftw_execute(plan_f)

        do i3=0,n3
            do i1=0,(n1+1)/2
                z(1,i1,i2,i3)=y2(1,i1,i3)
                z(2,i1,i2,i3)=y2(2,i1,i3)
            enddo
        enddo
    enddo
    call dfftw_destroy_plan(plan_f)

    call cpu_time(t1)

    call segment_invert(n1,n2,n3,kern_k1,kern_k3,c,z,hgrid)

    call cpu_time(t2)
    !! fourier transform x back; the result is y

    call dfftw_plan_dft_c2r_2d(plan_b,n1+1,n3+1,y2,x2,fftw_estimate)
    do i2=0,n2
        do i3=0,n3
            do i1=0,(n1+1)/2
                y2(1,i1,i3)=z(1,i1,i2,i3)
                y2(2,i1,i3)=z(2,i1,i2,i3)
            enddo
        enddo

        call dfftw_execute(plan_b)

        do i3=0,n3
            do i1=0,n1
                x(i1,i2,i3)=x2(i1,i3)
            enddo
        enddo
    enddo

    call dfftw_destroy_plan(plan_b)
    call cpu_time(t3)
    write(*,*)'total cpu time:',t3-t0,'sec'
    write(*,*)'inversion time:',t2-t1,'sec'

    deallocate(x2,y2)

END SUBROUTINE hit_with_kernel_slab
