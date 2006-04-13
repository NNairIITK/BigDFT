

   subroutine test_kernel(n1,n2,n3,hgrid,pot_ion,kernel,rhopot)
   implicit none
   integer :: n1,n2,n3
   real*8 :: a_gauss,a2
   real*8 :: rhotot,hgrid,shft1,shft2,shft3,ehart,eexcu,vexcu
   real*8, dimension(n1,n2,n3) :: rhopot,kernel,pot_ion
   real*8 :: pi,x1,x2,x3,r,r2,factor,derf,max_diff,diff,tt
   integer :: i1,i2,i3,ii1,ii2,ii3

   a_gauss=4.d0*hgrid
   a2 = a_gauss**2
   
! Shift the center of the Gaussian away from central grid point to break symmetries
   shft1=1.3d0*hgrid
   shft2=0.5d0*hgrid
   shft3=0.1d0*hgrid
   
   !Initialisation
   pi = 4.d0*atan(1.d0)
   !Normalisation
   factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
   !Gaussian function
   rhotot=0.d0
   do i3=1,n3
      x3 = hgrid*dble(i3-n3/2)-shft3
      do i2=1,n2
         x2 = hgrid*dble(i2-n2/2)-shft2
         do i1=1,n1
            x1 = hgrid*dble(i1-n1/2)-shft1
            r2 = x1*x1+x2*x2+x3*x3
            rhopot(i1,i2,i3) = factor*exp(-r2/a2)
            rhotot=rhotot+rhopot(i1,i2,i3)
         end do
      end do
   end do
   rhotot=rhotot*hgrid**3

!! Plot values along x axis
!   open(unit=11,file='rho.dat')
!   do i3=1,n3
!      x1 =                            - shft1
!      x2 =                            - shft2
!      x3 = hgrid*dble(i3-n3/2) - shft3
!      r=sqrt(x1**2+x2**2+x3**2)
!      write(unit=11,fmt="(e10.3,e12.5,2(e21.14),e9.2,2(e12.5))") &
!           r, rhopot(n1/2,n2/2,i3)
!   end do
!   close(unit=11)
   
   !Calculate potential using Poisson Solver
   write(*,*) 'testing poisson solver'
   call PSolver_Kernel(n1,n2,n3,hgrid,kernel,.false.,pot_ion,rhopot,ehart,eexcu,vexcu)
   
!! Plot values along x axis
!   open(unit=11,file='pot.dat')
!   do i3=1,n3
!      x1 =                            - shft1
!      x2 =                            - shft2
!      x3 = hgrid*dble(i3-n3/2) - shft3
!      r=sqrt(x1**2+x2**2+x3**2)
!      if (r == 0.d0) then
!         !limit_{x -> 0} erf(x/x) = 2/sqrt(pi)
!         factor = 2.d0/(sqrt(pi)*a_gauss)
!         tt=0.d0
!      else
!         factor = derf(r/a_gauss)/r
!         tt=abs(1.d0/r)
!      end if
!      write(unit=11,fmt="(e10.3,3(e21.14),e9.2,2(e12.5))") &
!           r, rhopot(n1/2,n2/2,i3),factor,tt
!   end do
!   close(unit=11)



! Global error 
   max_diff = 0.d0
   do i1=1,n1
      x1 = hgrid*dble(i1-n1/2) - shft1
      do i2=1,n2
         x2 = hgrid*dble(i2-n2/2) - shft2
         do i3=1,n3
            x3 = hgrid*dble(i3-n3/2) - shft3
            r=sqrt(x1**2+x2**2+x3**2)
             if (r == 0.d0) then
                !limit_{x -> 0} erf(x/x) = 2/sqrt(pi)
                factor = 2.d0/(sqrt(pi)*a_gauss)
             else
                factor = derf(r/a_gauss)/r
             end if
            diff=abs(rhopot(i1,i2,i3)-factor)
            if (diff.gt.max_diff) then
               max_diff=diff
               ii1=i1
               ii2=i2
               ii3=i3
            endif
         end do
      end do
   end do
   
   write(*,*) 'Testing Poisson Solver for a_gauss=',a_gauss
   write(*,'(a,f7.2,1x,e10.3,1x,e10.3)') 'hgridh,Deltarho,max_diff',hgrid,rhotot-1.d0,max_diff
   write(*,*) 'Max diff at : ',ii1,ii2,ii3
   write(*,*) 'Poisson Solver test finished'

end subroutine test_kernel




!!****h* Poisson_Solver/PSolver_Kernel
!! NAME
!!   PSolver
!!
!! FUNCTION
!!    Solves Poisson's equation by applying a kernel
!! replaces the charge density contained in rhopot by the Hartree stored as well in rhopot.
!! If xc_on is true it also adds the XC potential and ionic potential pot_ion
!!
!! SYNOPSIS
!!    Poisson solver applying a kernel and using Fourier transform for the convolution.
!!    We double the size of the mesh.
!! WARNING
!!    For the use of FFT routine
!!           inzee=1: first part of Z is data (output) array, 
!!                    second part work array
!!           inzee=2: first part of Z is work array, second part data array
!!                real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
!!                imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
!!           inzee on output is in general different from inzee on input
!!
!! AUTHOR
!!    Thierry Deutsch
!! COPYRIGHT
!!    Copyright (C) 2005 CEA
!! CREATION DATE
!!    13/07/2005
!!
!! MODIFICATION HISTORY
!!    SG added XC part
!!
!! SOURCE
!!
subroutine PSolver_Kernel(n01,n02,n03,hgrid,kernel,xc_on,pot_ion,rhopot,ehart,eexcu,vexcu)
   implicit none
   !Arguments
   logical xc_on
   integer :: n01,n02,n03
   real*8 :: hgrid
   real*8, dimension(n01,n02,n03) :: kernel
   real*8, dimension(n01,n02,n03) :: rhopot,pot_ion
   
   !Local variables
   real*8, dimension(:,:,:), allocatable :: zarray,karray
   real*8 :: factor,a,b,c,d,ehart,eexcu,vexcu
   integer :: i,n1,n2,n3,nd1,nd2,nd3
   integer :: inzee,inkee,i_sign,modulo
   integer :: iprint = 1

   !Initialisation
   call fourier_dim(2*n01,n1)
   call fourier_dim(2*n02,n2)
   call fourier_dim(2*n03,n3)
   
   nd1 = n1 + modulo(n1+1,2)
   nd2 = n2 + modulo(n2+1,2)
   nd3 = n3 + modulo(n3+1,2)

   if (iprint.eq.1) then
   write(*,*) 'initial  FFT dimension',2*n01,2*n02,2*n03
   write(*,*) 'adjusted FFT dimension',n1,n2,n3
   iprint=0
   endif

   !Allocation
   allocate(zarray(2,nd1*nd2*nd3,2))
   allocate(karray(2,nd1*nd2*nd3,2))
   
   !Set zarray
   inzee=1
   call zarray_in(n01,n02,n03,nd1,nd2,nd3,rhopot,zarray(1,1,inzee))
   !Set karray
   inkee=1
   call karray_in(n01,n02,n03,n1,n2,n3,nd1,nd2,nd3,kernel,karray(1,1,inkee))
   
   !FFT
!   print *,"Do a 3D FFT for the density"
   i_sign=1
   call fft(n1,n2,n3,nd1,nd2,nd3,zarray,i_sign,inzee)
   
!   print *,"Do a 3D FFT for the kernel"
   call fft(n1,n2,n3,nd1,nd2,nd3,karray,i_sign,inkee)
   
!   print *, "Apply the kernel"
   do i=1,nd1*nd2*nd3
      a=zarray(1,i,inzee)
      b=zarray(2,i,inzee)
      c=karray(1,i,inkee)
      d=karray(2,i,inkee)
      zarray(1,i,inzee) = a*c - b*d
      zarray(2,i,inzee) = a*d + b*c
   end do
   
   !Inverse FFT
   i_sign=-1
!   print *,"Do a 3D inverse FFT"
   call fft(n1,n2,n3,nd1,nd2,nd3,zarray,i_sign,inzee)
   !Recolt the result
   !We have to multiply by a factor
   factor = hgrid**3/(n1*n2*n3)
  if (xc_on) then
!   print *,"Add XC and ionic potential"
    call excpotu(2,n01,n02,n03,nd1,nd2,nd3,hgrid,factor,zarray(1,1,inzee),pot_ion,rhopot,ehart,eexcu,vexcu)
  else
!   Calling this routine gives only the Hartree potential
    print *,"PSolver does not add XC and ionic potential"
    call zarray_out(n01,n02,n03,nd1,nd2,nd3,zarray(1,1,inzee),factor,rhopot)
    ehart=0.d0 ; eexcu=0.d0 ; vexcu=0.d0
  endif
   
   !De-allocation
   deallocate(zarray)
   deallocate(karray)
   
end subroutine PSolver_Kernel


!!****h* Poisson_Solver/SPSolver_Kernel
!! NAME
!!   SPSolver
!!
!! FUNCTION
!!    Solves Poisson's equation by applying a kernel
!! replaces the charge density contained in rhopot by the Hartree stored as well in rhopot.
!! If xc_on is true it also adds the XC potential and ionic potential pot_ion
!!
!! SYNOPSIS
!!    Poisson solver applying a kernel and using Fourier transform for the convolution.
!!    We double the size of the mesh.
!! WARNING
!!    For the use of FFT routine
!!           inzee=1: first part of Z is data (output) array, 
!!                    second part work array
!!           inzee=2: first part of Z is work array, second part data array
!!                real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
!!                imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
!!           inzee on output is in general different from inzee on input
!!
!! AUTHOR
!!    Thierry Deutsch
!! COPYRIGHT
!!    Copyright (C) 2005 CEA
!! CREATION DATE
!!    13/07/2005
!!
!! MODIFICATION HISTORY
!!    SG added XC part
!!
!! SOURCE
!!
subroutine SPSolver_Kernel(n01,n02,n03,hgrid,kernel,xc_on,pot_ion,rhopot,ehart,eexcu,vexcu)
   implicit none
   !Arguments
   logical xc_on
   integer :: n01,n02,n03
   real*8 :: hgrid
   real*8, dimension(n01/2+8,n02/2+8,n03/2+8) :: kernel
   real*8, dimension(n01,n02,n03) :: rhopot,pot_ion
   
   !Local variables
   real*8, dimension(:,:,:), allocatable :: zarray,karray,srhopot,brhopot
   real*8 :: factor,a,b,c,d,ehart,eexcu,vexcu
   integer :: i,n1,n2,n3,nd1,nd2,nd3,ns1,ns2,ns3
   integer :: inzee,inkee,i_sign,modulo
   integer :: iprint = 1

!  smooth rhopt to make it a smaller array
   ns1=n01/2
   ns2=n02/2
   ns3=n03/2
   allocate(srhopot(ns1+8,ns2+8,ns3+8))
   call restrict(ns1,ns2,ns3,rhopot,srhopot)

   !Initialisation
   call fourier_dim(2*(ns1+8),n1)
   call fourier_dim(2*(ns2+8),n2)
   call fourier_dim(2*(ns3+8),n3)
   
   nd1 = n1 + modulo(n1+1,2)
   nd2 = n2 + modulo(n2+1,2)
   nd3 = n3 + modulo(n3+1,2)

   if (iprint.eq.1) then
   write(*,*) 'ns1+8,ns2+8,ns3+8 ',ns1+8,ns2+8,ns3+8
   write(*,*) 'initial  FFT dimension',2*(ns1+8),2*(ns2+8),2*(ns3+8)
   write(*,*) 'adjusted FFT dimension',n1,n2,n3
   write(*,*) 'leading FFT dimension',nd1,nd2,nd3
   iprint=0
   endif

   !Allocation
   allocate(zarray(2,nd1*nd2*nd3,2))
   allocate(karray(2,nd1*nd2*nd3,2))
   
   !Set zarray
   inzee=1
   call zarray_in(ns1+8,ns2+8,ns3+8,nd1,nd2,nd3,srhopot,zarray(1,1,inzee))
   !Set karray
   inkee=1
   call karray_in(ns1+8,ns2+8,ns3+8,n1,n2,n3,nd1,nd2,nd3,kernel,karray(1,1,inkee))
   
   !FFT
!   print *,"Do a 3D FFT for the density"
   i_sign=1
   call fft(n1,n2,n3,nd1,nd2,nd3,zarray,i_sign,inzee)
   
!   print *,"Do a 3D FFT for the kernel"
   call fft(n1,n2,n3,nd1,nd2,nd3,karray,i_sign,inkee)
   
!   print *, "Apply the kernel"
   do i=1,nd1*nd2*nd3
      a=zarray(1,i,inzee)
      b=zarray(2,i,inzee)
      c=karray(1,i,inkee)
      d=karray(2,i,inkee)
      zarray(1,i,inzee) = a*c - b*d
      zarray(2,i,inzee) = a*d + b*c
   end do
   

   !Inverse FFT
   i_sign=-1
!   print *,"Do a 3D inverse FFT"
   call fft(n1,n2,n3,nd1,nd2,nd3,zarray,i_sign,inzee)
   !Recolt the result
   !We have to multiply by a factor
   factor = (2.d0*hgrid)**3/(n1*n2*n3)
    call zarray_out(ns1+8,ns2+8,ns3+8,nd1,nd2,nd3,zarray(1,1,inzee),factor,srhopot)
   deallocate(zarray)
   deallocate(karray)


  if (xc_on) then
!   print *,"Add XC and ionic potential"
    allocate(brhopot(n01,n02,n03))
    call extrapolate(ns1,ns2,ns3,srhopot,brhopot)
    call excpotu(1,n01,n02,n03,n01,n02,n03,hgrid,1.d0,brhopot,pot_ion,rhopot,ehart,eexcu,vexcu)
    deallocate(brhopot)
  else
!   Calling this routine gives only the Hartree potential
    print *,"PSolver does not add XC and ionic potential"
    call extrapolate(ns1,ns2,ns3,srhopot,rhopot)
    ehart=0.d0 ; eexcu=0.d0 ; vexcu=0.d0
  endif
   
   !De-allocation
   deallocate(srhopot)
   
end subroutine SPSolver_Kernel


subroutine zarray_in(n01,n02,n03,nd1,nd2,nd3,density,zarray)
   implicit none
   !Arguments
   integer :: n01,n02,n03,nd1,nd2,nd3,i1,i2,i3
   real*8, dimension(n01,n02,n03) :: density
   real*8, dimension(2,nd1,nd2,nd3) :: zarray
   !Body
     do 300, i3=1,n03

     do 200,i2=1,n02
     do 100, i1=1,n01
     zarray(1,i1,i2,i3) = density(i1,i2,i3)
     zarray(2,i1,i2,i3) = 0.d0
100  continue
     do 101, i1 = n01+1, nd1
     zarray(1,i1,i2,i3) = 0.d0
     zarray(2,i1,i2,i3) = 0.d0
101  continue
200  continue

     do 210,i2=n02+1,nd2
     do 110, i1=1,nd1
     zarray(1,i1,i2,i3) = 0.d0
     zarray(2,i1,i2,i3) = 0.d0
110  continue
210  continue

300  continue

     do 333,i3=n03+1,nd3
     do 333,i2=1,nd2
     do 333,i1=1,nd1
     zarray(1,i1,i2,i3) = 0.d0
     zarray(2,i1,i2,i3) = 0.d0
333  continue

end subroutine zarray_in


subroutine karray_in(n01,n02,n03,n1,n2,n3,nd1,nd2,nd3,kernel,karray)
! this routine can be elimiated once we have a real FFT
   implicit none
   !Arguments
   integer :: n01,n02,n03,n1,n2,n3,nd1,nd2,nd3
   real*8, dimension(n01,n02,n03) :: kernel
   real*8, dimension(2,nd1,nd2,nd3) :: karray
   !Local variables
   integer :: i1,i2,i3
   
   !Body
   do i3=1,nd3
   do i2=1,nd2
   do i1=1,nd1
   karray(1,i1,i2,i3) = 0.d0
   karray(2,i1,i2,i3) = 0.d0
   enddo ; enddo ; enddo
   do i3=1,n03
      do i2=1,n02
         do i1=1,n01
            karray(1,i1,i2,i3) = kernel(i1,i2,i3)
         end do
         do i1=2,n01-1
            karray(1,n1-i1+2,i2,i3) = karray(1,i1,i2,i3)
            karray(2,n1-i1+2,i2,i3) = karray(2,i1,i2,i3)
         end do
      end do
      do i2=2,n02-1
         do i1=1,nd1
            karray(1,i1,n2-i2+2,i3) = karray(1,i1,i2,i3)
            karray(2,i1,n2-i2+2,i3) = karray(2,i1,i2,i3)
         end do
      end do
   end do
   do i3=2,n03-1
      do i2=1,nd2
         do i1=1,nd1
            karray(1,i1,i2,n3-i3+2) = karray(1,i1,i2,i3)
            karray(2,i1,i2,n3-i3+2) = karray(2,i1,i2,i3)
         end do
      end do
   end do
end subroutine karray_in
   
subroutine zarray_out(n1,n2,n3,nd1,nd2,nd3,zarray,factor,potential)
   implicit none
   !Arguments
   integer :: n1,n2,n3,nd1,nd2,nd3,i1,i2,i3
   real*8, dimension(n1,n2,n3) :: potential
   real*8, dimension(2,nd1,nd2,nd3) :: zarray
   real*8 :: factor
   !Body 
   do i3=1,n3
   do i2=1,n2
   do i1=1,n1
   potential(i1,i2,i3) = factor*zarray(1,i1,i2,i3)
   enddo ; enddo ; enddo
end subroutine zarray_out


	subroutine excpotu(ncmplx,n1,n2,n3,nd1,nd2,nd3,hgrid,factor,zarray,pot_ion,rhopot,ehart,eexcu,vexcu)
! Takes the Hartree potential as contained in zarray from the Poisson solver 
! and sums the Exc and Hartree potential in the array rhopot
! Calculates also the Hartree energy ehart and XC energy eexcu
	implicit real*8 (a-h,o-z)
      parameter (a0u=.4581652932831429d0, a1u=2.217058676663745d0,  & 
                 a2u=0.7405551735357053d0,a3u=0.01968227878617998d0)
      parameter (b1u=1.0d0,               b2u=4.504130959426697d0,  & 
                 b3u=1.110667363742916d0, b4u=0.02359291751427506d0)
      parameter  & 
      (c1u=4.d0*a0u*b1u/3.0d0,  c2u=5.0d0*a0u*b2u/3.0d0+a1u*b1u,  & 
              c3u=2.0d0*a0u*b3u+4.0d0*a1u*b2u/3.0d0+2.0d0*a2u*b1u/3.0d0,  & 
      c4u=7.0d0*a0u*b4u/3.0d0+5.0d0*a1u*b3u/3.0d0+a2u*b2u+a3u*b1u/3.0d0,  & 
              c5u=2.0d0*a1u*b4u+4.0d0*a2u*b3u/3.0d0+2.0d0*a3u*b2u/3.0d0,  & 
              c6u=5.0d0*a2u*b4u/3.0d0+a3u*b3u,c7u=4.0d0*a3u*b4u/3.0d0)
      parameter (rsfac=.6203504908994000d0)
!      parameter(eps2=1.d-28)
      parameter(thirdm=-1.d0/3.d0)

      dimension rhopot(n1,n2,n3),zarray(ncmplx,nd1,nd2,nd3),pot_ion(n1,n2,n3)

	eexcu=0.d0
	vexcu=0.d0
	ehart=0.d0
!	x1=1.d0
!	ic=0
	do 205,i3=1,n3
	do 205,i2=1,n2
	do 205,i1=1,n1

	rhou1=rhopot(i1,i2,i3)

      if (rhou1.lt.1.d-20) then
        h1=factor*zarray(1,i1,i2,i3)
        ehart=ehart+rhou1*h1
        rhopot(i1,i2,i3)=h1+pot_ion(i1,i2,i3)
      else

!        rsu1=rsfac*rhou1**(thirdm)
        rsu1=rsfac*exp(thirdm*log(rhou1))
	
!10	continue
!	d1=x1-rhou1/x1**2
!	ic=ic+1
!	x1=x1-.333333333333333d0*d1
!	if (d1**2.gt.eps2) goto 10
!	rsu1=rsfac/x1

	topu1=a2u+rsu1*a3u
	topu1=a1u+rsu1*topu1
	topu1=a0u+rsu1*topu1
	dtopu1=c6u+rsu1*c7u
	dtopu1=c5u+rsu1*dtopu1
	dtopu1=c4u+rsu1*dtopu1
	dtopu1=c3u+rsu1*dtopu1
	dtopu1=c2u+rsu1*dtopu1
	dtopu1=c1u+rsu1*dtopu1
	dtopu1=-rsu1*dtopu1
	botu1=b3u+rsu1*b4u
	botu1=b2u+rsu1*botu1
	botu1=b1u+rsu1*botu1
	botu1=rsu1*botu1
	t1=1.d0/botu1
        epsxcu1=-topu1*t1

	eexcu=eexcu+epsxcu1*rhou1
	vexcu=vexcu+(dtopu1*t1*t1)*rhou1
        h1=factor*zarray(1,i1,i2,i3)
        p1=pot_ion(i1,i2,i3)
        ehart=ehart+rhou1*h1
        rhopot(i1,i2,i3)=h1+(dtopu1*t1*t1)+p1

      endif

205	continue
	eexcu=eexcu*hgrid**3
	vexcu=vexcu*hgrid**3
	ehart=.5d0*ehart*hgrid**3

!	write(6,*) 'ehart,eexcu,vexcu',ehart,eexcu,vexcu
!	write(6,*) 'average iterations for root in excpotu',2.d0*ic/(n1*n2*n3)

	return
	end subroutine excpotu


subroutine fourier_dim(n,n_next)
   implicit none
   !Arguments
   integer :: i,n,n_next
   
   integer, parameter :: ndata = 82
   integer, dimension(ndata) :: idata 
              data  ( idata(i),i=1,82 ) /   &
               3,    4,   5,     6,    8,    9,   12,   15,   16,   18, &
              20,   24,   25,   27,   30,   32,   36,   40,   45,   48, &
              54,   60,   64,   72,   75,   80,   81,   90,   96,  100, &
             108,  120,  125,  128,  135,  144,  150,  160,  162,  180, &
             192,  200,  216,  225,  240,  243,  256,  270,  288,  300, &
             320,  324,  360,  375,  384,  400,  405,  432,  450,  480, &
             486,  500,  512,  540,  576,  600,  625,  640,  648,  675, &
             720,  729,  750,  768,  800,  810,  864,  900,  960,  972, &
            1000,  1024 /
         
   loop_data: do i=1,ndata
      if (n <= idata(i)) then
         n_next = idata(i)
         return
      end if
   end do loop_data
   write(unit=*,fmt=*) "fourier_dim: ",n," is bigger than ",idata(ndata)
   stop
end subroutine fourier_dim

