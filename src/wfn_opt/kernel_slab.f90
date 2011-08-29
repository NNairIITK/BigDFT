!> @file
!!  Routines using kernels
!! @author
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
 

!>   Hits the input array x with the kernel
!!   @f$((-1/2\Delta+C)_{ij})^{-1}@f$
!! See the optimized version (hit_kernel_slab_optim) and
!!     non-optimized version (hit_kernal_slab_simple)
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

    call slab_invert(n1,n2,n3,z,kern_k1,kern_k3,c,hgrid)

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


subroutine slab_invert(n1,n2,n3,z,kern_k1,kern_k3,c,hgrid)
    implicit none
    integer,intent(in)::n1,n2,n3
    real*8,intent(in)::kern_k1(0:n1)
    real*8,intent(in)::kern_k3(0:n3)
    real*8,intent(in)::c,hgrid

    real*8 z(2,0:(n1+1)/2,0:n2,0:n3)
    !*******************************************************************************************
    ! for the matrix inversion with lapack:
    real*8,allocatable,dimension(:,:)::ab,amat,amat1
    real*8,allocatable::b(:,:)
    !     .. Scalar Arguments ..
    INTEGER::INFO, Kd, LDAB, LDB, NRHS=2,n
    integer,allocatable::ipiv(:)
    integer i1,i2,i3
    integer j,l,i
    real*8 ct

    include 'd2s16.inc'

    fil=fil*(n1+1)*(n3+1) ! include the factor from the Fourier transform

    ! initialize the variables for matrix inversion
    n=n2+1
    kd=lupfil
    ldab=kd+1!
    ldb=n2+1
    allocate(ab(ldab,n))
    allocate(b(0:n2,2))
    allocate(ipiv(n2+1))

    !! hit the fourier transform of x with the kernel
    do i3=0,n3
        do i1=0,(n1+1)/2
            ! include the factor from the Fourier transform
            ct=(kern_k1(i1)+kern_k3(i3)+c)*(n1+1)*(n3+1)
    
            ! ab has to be reinitialized each time
            ! since it is overwritten in the course of  dgbsv
            do j=1,n
                do i=max(1,j-lupfil),j
                    ab(kd+1+i-j,j)=fil(i-j)
                enddo
                ab(kd+1,j)=fil(0)+ct
            enddo
    
            do i2=0,n2
                b(i2,1)=z(1,i1,i2,i3)
                b(i2,2)=z(2,i1,i2,i3)
            enddo
    
            call DPBSV( 'U', N, KD, NRHS, AB, LDAB, B, LDB, INFO )
    
            do i2=0,n2
                z(1,i1,i2,i3)=b(i2,1)
                z(2,i1,i2,i3)=b(i2,2)
            enddo
    
        enddo
    enddo

    deallocate(b,ipiv)

END SUBROUTINE slab_invert


!>  Construct the kernel (-1/2 d^2/dx^2)_{ij}
!!  at a real space grid with grid size hgrid
!!  and then fourier transform it to momentum space
subroutine make_kernel_slab(n1,hgrid,kern)
implicit none
integer,intent(in)::n1
real*8,intent(in)::hgrid
real*8,intent(out)::kern(0:n1)

include 'fftw3.f'
integer*8::plan_f
real*8,allocatable::z(:,:),zin(:,:)
real*8 dmax,tt
integer i
include 'd2s16.inc'

allocate(z(2,0:n1))
allocate(zin(2,0:n1))

! construct the kernel in real space
zin=0.d0
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
END SUBROUTINE make_kernel_slab
