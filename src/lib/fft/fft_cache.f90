!!****p* BigDFT/fft_check
!!
!! DESCRIPTION
!!   3-dimensional complex-complex FFT routine: 
!!   When compared to the best vendor implementations on RISC architectures 
!!   it gives close to optimal performance (perhaps loosing 20 percent in speed)
!!   and it is significanly faster than many not so good vendor implementations 
!!   as well as other portable FFT's. 
!!   On all vector machines tested so far (Cray, NEC, Fujitsu) is 
!!   was significantly faster than the vendor routines
!!
!! The theoretical background is described in :
!! 1) S. Goedecker: Rotating a three-dimensional array in optimal
!! positions for vector processing: Case study for a three-dimensional Fast
!! Fourier Transform, Comp. Phys. Commun. \underline{76}, 294 (1993)
!! Citing of this reference is greatly appreciated if the routines are used 
!! for scientific work.
!!
!! COPYRIGHT
!!   Copyright (C) Stefan Goedecker, CEA Grenoble, 2002, Basel University, 2009
!!   This file is distributed under the terms of the
!!   GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! SYNOPSIS
!!   Presumably good compiler flags:
!!   IBM, serial power 2: xlf -qarch=pwr2 -O2 -qmaxmem=-1
!!   with OpenMP: IBM: xlf_r -qfree -O4 -qarch=pwr3 -qtune=pwr3 -qsmp=omp -qmaxmem=-1 ; 
!!                     a.out
!!   DEC: f90 -O3 -arch ev67 -pipeline
!!   with OpenMP: DEC: f90 -O3 -arch ev67 -pipeline -omp -lelan ; 
!!                     prun -N1 -c4 a.out
!!
!! SOURCE
!!
program fft_cache

   use module_fft_sg
   implicit none
   integer, parameter :: iunit=21
   integer :: i,j,n1,ntime
   real(kind=8) :: t1,t2
   logical :: ok_test
   !write(*,'(a)') 'FFT test: (n1,n2,n3)'
   open(unit=iunit,file="fft_cache.dat",status="unknown",position="append")
   do i=1,ndata
      call cpu_time(t1)
      n1 = i_data(i)
      ntime=2*(nfft_max/n1+1)
      do j=1,ntime
         !call do_fft(n1, 3, 3,j)
         call do_fft(n1, 3, 3, 0, ok_test)
         if (.not.ok_test) STOP 'fft_cache: One test has failed!'
      end do
      call cpu_time(t2)
      !write(unit=6,fmt=*) j,t2-t1
      write(unit=iunit,fmt=*) ncache,n1,t2-t1
   end do
   write(unit=iunit,fmt=*) 
   close(unit=iunit)

contains

   subroutine do_fft(n1,n2,n3,itime,ok_test)

      implicit none

      ! dimension parameters
      integer, intent(in) :: n1, n2, n3, itime
      logical, intent(out) :: ok_test
      ! Local variables
      integer :: count1,count2,count_rate,count_max,i,inzee,i_sign
      real(kind=8) :: ttm,tta,t1,t2,tela,time,flops
      ! parameters for FFT
      integer :: nd1, nd2, nd3
      ! general array
      real(kind=8), allocatable :: zin(:,:)
      ! arrays for FFT 
      real(kind=8), allocatable :: z(:,:,:)
      character(len=10) :: message

      nd1=n1+1
      nd2=n2+1
      nd3=n3+1
      if (itime == 1) write(6,'(i7,3(i6))',advance='no') ncache,n1,n2,n3

      ! Allocations
      allocate(zin(2,n1*n2*n3))
      allocate(z(2,nd1*nd2*nd3,2))

      do i=1,nd1*nd2*nd3
         z(1,i,1)=0.d0
         z(2,i,1)=0.d0
         z(1,i,2)=0.d0
         z(2,i,2)=0.d0
      end do

        call init(n1,n2,n3,nd1,nd2,nd3,zin,z)

        i_sign=-1
        inzee=1
        call fft(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)

        call cpu_time(t1)
        call system_clock(count1,count_rate,count_max)      

        i_sign=1
        call fft(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)

        call cpu_time(t2)
        call system_clock(count2,count_rate,count_max)      
        time=(t2-t1)
        tela=(count2-count1)/real(count_rate,kind=8)

        call vgl(n1,n2,n3,nd1,nd2,nd3,z(1,1,inzee), &
                       n1,n2,n3,zin,1.d0/real(n1*n2*n3,kind=8),tta,ttm)
        if (ttm.gt.1.d-10) then
           message = 'Failed'
           ok_test = .False.
        else
           message = 'Succeeded'
           ok_test = .True.
        end if
        flops=5*n1*n2*n3*log(1.d0*n1*n2*n3)/log(2.d0)
        !write(6,'(a,2(x,e11.4),x,i4)')  'Time (CPU,ELA) per FFT call (sec):' ,time,tela
        !write(6,*) 'Estimated floating point operations per FFT call',flops
        !write(6,*) 'CPU MFlops',1.d-6*flops/time
        if (itime == 1) write(6,'(1x,a,2(1pg9.2),1x,a)',advance='no') 'Backw<>Forw:ttm=,tta=',ttm,tta,message

        ! De-allocations
        deallocate(z)
        deallocate(zin)

   end subroutine do_fft

   subroutine init(n1,n2,n3,nd1,nd2,nd3,zin,z)
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3
      real*8 :: zin(2,n1,n2,n3),z(2,nd1,nd2,nd3)
      !Local variables
      integer :: i1,i2,i3
      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               zin(1,i1,i2,i3) = cos(1.23d0*real(i1*111 + i2*11 + i3,kind=8))
               zin(2,i1,i2,i3) = sin(3.21d0*real(i3*111 + i2*11 + i1,kind=8))
               z(1,i1,i2,i3) = zin(1,i1,i2,i3) 
               z(2,i1,i2,i3) = zin(2,i1,i2,i3) 
            end do
         end do
      end do
   end subroutine init


   subroutine vgl(n1,n2,n3,nd1,nd2,nd3,x,md1,md2,md3,y,scale,tta,ttm)
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3
      real(kind=8), intent(in) :: x(2,nd1,nd2,nd3),y(2,md1,md2,md3)
      real(kind=8), intent(in) :: scale
      !Local variables
      real(kind=8) :: ttm,tta,ttr,tti
      integer :: i1,i2,i3
      ttm=0.d0
      tta=0.d0
      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               ttr=abs(x(1,i1,i2,i3)*scale-y(1,i1,i2,i3))/abs(y(1,i1,i2,i3))
               tti=abs(x(2,i1,i2,i3)*scale-y(2,i1,i2,i3))/abs(y(2,i1,i2,i3))
               ttm=max(ttr,tti,ttm)
               tta=tta+ttr+tti
            end do
         end do
      end do
      tta=tta/(n1*n2*n3)
   end subroutine vgl

end program fft_cache
!!***

