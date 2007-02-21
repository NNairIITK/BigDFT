!!****f* BigDFT/PSolver_Kernel
!! NAME
!!   PSolver_Kernel
!!
!! FUNCTION
!!    Solver of Poisson equation applying a kernel
!!
!! SYNOPSIS
!!    Poisson solver applying a kernel and 
!!    using Fourier transform for the convolution.
!!    rhopot : input  -> the density
!!             output -> the Hartree potential + pot_ion
!!    The potential pot_ion is ADDED in the array rhopot.
!!    Calculate also the Hartree potential
!!
!!    Replaces the charge density contained in rhopot 
!!    by the Hartree stored as well in rhopot.
!!    The XC potential is chosen from the value of ixc, 
!!    by following the same rules as ABINIT.
!!    If ixc is 0, it also adds the XC potential and
!!    ionic potential pot_ion
!!
!!    We double the size of the mesh except in one dimension
!!    in order to use the property of the density to be real.
!! WARNING
!!    For the use of FFT routine
!!        inzee=1: first part of Z is data (output) array, 
!!                 second part work array
!!        inzee=2: first part of Z is work array, second part data array
!!                 real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
!!                 imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
!!        inzee on output is in general different from inzee on input
!!
!! AUTHOR
!!    Thierry Deutsch, Luigi Genovese
!! COPYRIGHT
!!    Copyright (C) 2005 CEA
!! CREATION DATE
!!    13/07/2005
!!
!! MODIFICATION HISTORY
!!    12/2005 Kernel stored into memory
!!    12/2005 Real Kernel FFT and use less memory
!!
!! SOURCE
!!
subroutine PSolver_Kernel(n01,n02,n03,nfft1,nfft2,nfft3, &
     hgrid,karray,ixc,pot_ion,rhopot,ehartree,eexcu,vexcu)
   implicit none
   !Arguments
   integer, intent(in)  :: n01,n02,n03,nfft1,nfft2,nfft3,ixc
   real*8, intent(in) :: hgrid
   !logical, intent(in) :: xc_on
   real*8, intent(in), dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1) :: karray
   real*8, intent(in), dimension(n01,n02,n03) :: pot_ion
   real*8, intent(inout), dimension(n01,n02,n03) :: rhopot
   real*8, intent(out) :: ehartree,eexcu,vexcu
   !Local variables
   real*8, dimension(:,:,:), allocatable :: zarray
   real*8 :: factor,exc,vxc,eht
   integer :: n1,n2,n3,nd1,nd2,nd3,n1h,nd1h
   integer :: inzee,i_sign,i_allocated

   call timing(0,'PSolv_comput  ','ON')

   !Dimension of the FFT
   call dimensions_FFT(n01,n02,n03,n1,n2,n3)
   !Half size of nd1
   n1h=n1/2
   nd1 = n1 + modulo(n1+1,2)
   nd2 = n2 + modulo(n2+1,2)
   nd3 = n3 + modulo(n3+1,2)
   nd1h=(nd1+1)/2
   !Allocations
   i_allocated=0
   allocate(zarray(2,nd1h*nd2*nd3,2),stat=i_allocated)
   if (i_allocated /= 0) then
      print *,"PSolver_Kernel:Problem of memory allocation"
      stop
   end if
   !Set zarray
   call zarray_in(n01,n02,n03,nd1h,nd2,nd3,rhopot,zarray)

   !FFT
   !print *,"Do a 3D HalFFT for the density"
   i_sign=1
   inzee=1
   call fft(n1h,n2,n3,nd1h,nd2,nd3,zarray,i_sign,inzee)
 
   !print *, "Apply the kernel"
   call kernel_application(n1,n2,n3,nd1h,nd2,nd3,nfft1,nfft2,nfft3,zarray,karray,inzee)

   !Inverse FFT
   i_sign=-1
   !print *,"Do a 3D inverse HalFFT"
   call fft(n1h,n2,n3,nd1h,nd2,nd3,zarray,i_sign,inzee)
 
   !Recollect the result
   !We have to multiply by a factor
   factor = hgrid**3/(n1*n2*n3)


   if (ixc /= 0) then
      call timing(0,'PSolv_comput  ','OF')

      call timing(0,'Exchangecorr  ','ON')
      call xc_energy(n01,n02,n03,2*nd1h,nd2,nd3,ixc,factor,hgrid,rhopot,&
           pot_ion,zarray(1,1,inzee),ehartree,eexcu,vexcu)
      call timing(0,'Exchangecorr  ','OF')
      write(*,'(1x,a,2(f20.12),3x,a,f20.12)') "The xc energies are",eexcu,vexcu,"Hartree",ehartree

   else
      ! Calling this routine gives only the Hartree potential
      call zarray_out(n01,n02,n03,nd1h,nd2,nd3,&
           rhopot,zarray(1,1,inzee),factor,hgrid,ehartree)
      write(*,'(1x,a,f20.12)') "The Hartree energy is",ehartree
      eexcu=0.d0
      vexcu=0.d0
         
      call timing(0,'PSolv_comput  ','OF')

   endif

   !De-allocations
   deallocate(zarray)
end subroutine PSolver_Kernel
!!***

!!****f* BigDFT/kernel_application
!! NAME
!!   kernel_application
!!
!! FUNCTION
!!    Multiply the FFT of the density by the FFT of the kernel
!!
!! SYNOPSIS
!!    zarray(:,:,:,:,inzee) : IN -> FFT of the density with the x dimension divided by two
!!                            (HalFFT), OUT -> FFT of the potential
!!    karray                : kernel FFT (real, 1/8 of the total grid)
!!    n1h,n2,n3             : dimension of the FFT grid for zarray
!!    nd1h,nd2,nd3          : dimensions of zarray
!!    nfft1,nfft2,nfft3     : original FFT grid dimensions, to be used for karray dimensions
!!
!! WARNING
!!    We use all the zarray vector, storing the auxiliary part using ouzee=3-inzee
!!    All the loop are unrolled such to avoid different conditions
!!    the "min" functions are substituted by kink computed with absolute values
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    March 2006
!!
!! SOURCE
!!
subroutine kernel_application(n1,n2,n3,nd1h,nd2,nd3,nfft1,nfft2,nfft3,zarray,karray,inzee)
   implicit none
   !Arguments
   integer, intent(in)  :: n1,n2,n3,nd1h,nd2,nd3,nfft1,nfft2,nfft3,inzee
   real*8, intent(in), dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1) :: karray
   real*8, intent(inout), dimension(2,nd1h,nd2,nd3,2) :: zarray
   !Local variables
   real*8, dimension(:), allocatable :: cos_array,sin_array
   real*8 :: a,b,c,d,pi2,g1,cp,sp
   real*8 :: rfe,ife,rfo,ifo,rk,ik,rk2,ik2,re,ro,ie,io,rhk,ihk
   integer :: i1,i2,i3,j1,j2,j3,i_allocated,i_stat,ouzee,n1h,n2h,n3h
   integer :: si1,si2,si3

   !Body
   n1h=n1/2
   n2h=n2/2
   n3h=n3/2
   !Allocations
   i_allocated=0
   allocate(cos_array(n1h+1),stat=i_stat)
   i_allocated=i_allocated+i_stat
   allocate(sin_array(n1h+1),stat=i_stat)
   i_allocated=i_allocated+i_stat
   if (i_allocated /= 0) then
      print *,"kernel_application:Problem of memory allocation"
      stop
   end if

   pi2=8.d0*datan(1.d0)
   pi2=pi2/real(n1,kind=8)
   do i1=1,n1h+1
      cos_array(i1)=dcos(pi2*(i1-1))
      sin_array(i1)=-dsin(pi2*(i1-1))
   end do
   
   ouzee=3-inzee


       
!--------------------------------------------!
!--- Starting reconstruction half -> full ---!
!--------------------------------------------!   

   !-------------Case i3 = 1
   i3=1
   j3=1
   si3=1
   
   !-------------Case i2 = 1, i3 = 1
   i2=1
   j2=1
   si2=1
   
   !Case i1 == 1
   i1=1
   si1=1
   a=zarray(1,i1,i2,i3,inzee)
   b=zarray(2,i1,i2,i3,inzee)
   c=zarray(1,si1,si2,si3,inzee)
   d=zarray(2,si1,si2,si3,inzee)
   rfe=.5d0*(a+c)
   ife=.5d0*(b-d)
   rfo=.5d0*(a-c)
   ifo=.5d0*(b+d) 
   cp=cos_array(i1)
   sp=sin_array(i1)
   rk=rfe+cp*ifo-sp*rfo
   ik=ife-cp*rfo-sp*ifo
   g1=karray(i1,j2,j3)
   rk2=rk*g1
   ik2=ik*g1
   
   zarray(1,1,i2,i3,ouzee) = rk2
   zarray(2,1,i2,i3,ouzee) = ik2
   
   !Case i1=2,n1h
   do i1=2,n1h
      si1=n1h+2-i1
      
      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1
      
      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
   end do
   
   !Case i1=n1h+1
   i1=n1h+1
   si1=n1h+2-i1
   
   a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
   b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
   c=zarray(1,si1,si2,si3,inzee)
   d=zarray(2,si1,si2,si3,inzee)
   rfe=.5d0*(a+c)
   ife=.5d0*(b-d)
   rfo=.5d0*(a-c)
   ifo=.5d0*(b+d) 
   cp=cos_array(i1)
   sp=sin_array(i1)
   rk=rfe+cp*ifo-sp*rfo
   ik=ife-cp*rfo-sp*ifo
   g1=karray(i1,j2,j3)
   rk2=rk*g1
   ik2=ik*g1
   
   zarray(1,i1,i2,i3,ouzee) = rk2
   zarray(2,i1,i2,i3,ouzee) = ik2
   !-------------END case i2 = 1 , i3=1
   
   !case i2 >=2
   do i2=2,n2
      j2=n2h+1-abs(n2h+1-i2)
      si2=n2+2-i2 !if i2 /=1, otherwise si2=1
      
      !Case i1 == 1
      i1=1
      si1=1
      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1
      
      zarray(1,1,i2,i3,ouzee) = rk2
      zarray(2,1,i2,i3,ouzee) = ik2
      
      !Case i1=2,n1h
      do i1=2,n1h
         si1=n1h+2-i1
         
         a=zarray(1,i1,i2,i3,inzee)
         b=zarray(2,i1,i2,i3,inzee)
         c=zarray(1,si1,si2,si3,inzee)
         d=zarray(2,si1,si2,si3,inzee)
         rfe=.5d0*(a+c)
         ife=.5d0*(b-d)
         rfo=.5d0*(a-c)
         ifo=.5d0*(b+d) 
         cp=cos_array(i1)
         sp=sin_array(i1)
         rk=rfe+cp*ifo-sp*rfo
         ik=ife-cp*rfo-sp*ifo
         g1=karray(i1,j2,j3)
         rk2=rk*g1
         ik2=ik*g1
         
         zarray(1,i1,i2,i3,ouzee) = rk2
         zarray(2,i1,i2,i3,ouzee) = ik2
      end do
      
      !Case i1=n1h+1
      i1=n1h+1
      si1=n1h+2-i1
      
      a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
      b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1
      
      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
   end do
   !-------------END Case i3 = 1

   !case i3 >=2
   do i3=2,n3
      j3=n3h+1-abs(n3h+1-i3)
      si3=n3+2-i3 !if i3 /=1, otherwise si3=1

      !-------------Case i2 = 1
      i2=1
      j2=1
      si2=1
      
      !Case i1 == 1
      i1=1
      si1=1
      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1
      
      zarray(1,1,i2,i3,ouzee) = rk2
      zarray(2,1,i2,i3,ouzee) = ik2
      
      !Case i1=2,n1h
      do i1=2,n1h
         si1=n1h+2-i1
         
         a=zarray(1,i1,i2,i3,inzee)
         b=zarray(2,i1,i2,i3,inzee)
         c=zarray(1,si1,si2,si3,inzee)
         d=zarray(2,si1,si2,si3,inzee)
         rfe=.5d0*(a+c)
         ife=.5d0*(b-d)
         rfo=.5d0*(a-c)
         ifo=.5d0*(b+d) 
         cp=cos_array(i1)
         sp=sin_array(i1)
         rk=rfe+cp*ifo-sp*rfo
         ik=ife-cp*rfo-sp*ifo
         g1=karray(i1,j2,j3)
         rk2=rk*g1
         ik2=ik*g1
         
         zarray(1,i1,i2,i3,ouzee) = rk2
         zarray(2,i1,i2,i3,ouzee) = ik2
      end do
      
      !Case i1=n1h+1
      i1=n1h+1
      si1=n1h+2-i1
      
      a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
      b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1
      
      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
      !-------------END case i2 = 1      
         
      !case i2 >=2
      do i2=2,n2
         j2=n2h+1-abs(n2h+1-i2)
         si2=n2+2-i2 !if i2 /=1, otherwise si2=1

         !Case i1 == 1
         i1=1
         si1=1
         a=zarray(1,i1,i2,i3,inzee)
         b=zarray(2,i1,i2,i3,inzee)
         c=zarray(1,si1,si2,si3,inzee)
         d=zarray(2,si1,si2,si3,inzee)
         rfe=.5d0*(a+c)
         ife=.5d0*(b-d)
         rfo=.5d0*(a-c)
         ifo=.5d0*(b+d) 
         cp=cos_array(i1)
         sp=sin_array(i1)
         rk=rfe+cp*ifo-sp*rfo
         ik=ife-cp*rfo-sp*ifo
         g1=karray(i1,j2,j3)
         rk2=rk*g1
         ik2=ik*g1
         
         zarray(1,1,i2,i3,ouzee) = rk2
         zarray(2,1,i2,i3,ouzee) = ik2
         
         !Case i1=2,n1h
         do i1=2,n1h
            si1=n1h+2-i1
            
            a=zarray(1,i1,i2,i3,inzee)
            b=zarray(2,i1,i2,i3,inzee)
            c=zarray(1,si1,si2,si3,inzee)
            d=zarray(2,si1,si2,si3,inzee)
            rfe=.5d0*(a+c)
            ife=.5d0*(b-d)
            rfo=.5d0*(a-c)
            ifo=.5d0*(b+d) 
            cp=cos_array(i1)
            sp=sin_array(i1)
            rk=rfe+cp*ifo-sp*rfo
            ik=ife-cp*rfo-sp*ifo
            g1=karray(i1,j2,j3)
            rk2=rk*g1
            ik2=ik*g1
            
            zarray(1,i1,i2,i3,ouzee) = rk2
            zarray(2,i1,i2,i3,ouzee) = ik2
         end do
         
         !Case i1=n1h+1
         i1=n1h+1
         si1=n1h+2-i1
         
         a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
         b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
         c=zarray(1,si1,si2,si3,inzee)
         d=zarray(2,si1,si2,si3,inzee)
         rfe=.5d0*(a+c)
         ife=.5d0*(b-d)
         rfo=.5d0*(a-c)
         ifo=.5d0*(b+d) 
         cp=cos_array(i1)
         sp=sin_array(i1)
         rk=rfe+cp*ifo-sp*rfo
         ik=ife-cp*rfo-sp*ifo
         g1=karray(i1,j2,j3)
         rk2=rk*g1
         ik2=ik*g1
         
         zarray(1,i1,i2,i3,ouzee) = rk2
         zarray(2,i1,i2,i3,ouzee) = ik2
      end do

   end do


!--------------------------------------------!
!--- Starting reconstruction full -> half ---!
!--------------------------------------------!   

   !case i3=1
   i3=1
   j3=1
   !case i2=1
   i2=1
   j2=1
   do i1=1,n1h
      j1=n1h+2-i1
      
      a=zarray(1,i1,i2,i3,ouzee)
      b=zarray(2,i1,i2,i3,ouzee)
      c=zarray(1,j1,j2,j3,ouzee)
      d=-zarray(2,j1,j2,j3,ouzee)
      cp=cos_array(i1)
      sp=sin_array(i1)
      re=(a+c)
      ie=(b+d)
      ro=(a-c)*cp-(b-d)*sp
      io=(a-c)*sp+(b-d)*cp
      rhk=re-io 
      ihk=ie+ro
      
      zarray(1,i1,i2,i3,inzee)=rhk
      zarray(2,i1,i2,i3,inzee)=ihk
   end do
   !case i2 >= 2
   do i2=2,n2
      j2=nd2+1-i2
      do i1=1,n1h
         j1=n1h+2-i1
         
         a=zarray(1,i1,i2,i3,ouzee)
         b=zarray(2,i1,i2,i3,ouzee)
         c=zarray(1,j1,j2,j3,ouzee)
         d=-zarray(2,j1,j2,j3,ouzee)
         cp=cos_array(i1)
         sp=sin_array(i1)
         re=(a+c)
         ie=(b+d)
         ro=(a-c)*cp-(b-d)*sp
         io=(a-c)*sp+(b-d)*cp
         rhk=re-io 
         ihk=ie+ro
         
         zarray(1,i1,i2,i3,inzee)=rhk
         zarray(2,i1,i2,i3,inzee)=ihk
      end do
   end do
   
   
   !case i3 >=2
   do i3=2,n3
      j3=nd3+1-i3
      !case i2=1
      i2=1
      j2=1
      do i1=1,n1h
         j1=n1h+2-i1
         
         a=zarray(1,i1,i2,i3,ouzee)
         b=zarray(2,i1,i2,i3,ouzee)
         c=zarray(1,j1,j2,j3,ouzee)
         d=-zarray(2,j1,j2,j3,ouzee)
         cp=cos_array(i1)
         sp=sin_array(i1)
         re=(a+c)
         ie=(b+d)
         ro=(a-c)*cp-(b-d)*sp
         io=(a-c)*sp+(b-d)*cp
         rhk=re-io 
         ihk=ie+ro
         
         zarray(1,i1,i2,i3,inzee)=rhk
         zarray(2,i1,i2,i3,inzee)=ihk
      end do
      !case i2 >= 2
      do i2=2,n2
         j2=nd2+1-i2
         do i1=1,n1h
            j1=n1h+2-i1

            a=zarray(1,i1,i2,i3,ouzee)
            b=zarray(2,i1,i2,i3,ouzee)
            c=zarray(1,j1,j2,j3,ouzee)
            d=-zarray(2,j1,j2,j3,ouzee)
            cp=cos_array(i1)
            sp=sin_array(i1)
            re=(a+c)
            ie=(b+d)
            ro=(a-c)*cp-(b-d)*sp
            io=(a-c)*sp+(b-d)*cp
            rhk=re-io 
            ihk=ie+ro

            zarray(1,i1,i2,i3,inzee)=rhk
            zarray(2,i1,i2,i3,inzee)=ihk
         end do
      end do

     end do

   !De-allocations
   deallocate(cos_array)
   deallocate(sin_array)

 end subroutine kernel_application



!!****f* BigDFT/norm_ind
!! NAME
!!   norm_ind
!!
!! FUNCTION
!!   Index in zarray
!!
!! SOURCE
!!
subroutine norm_ind(nd1,nd2,nd3,i1,i2,i3,ind)
  implicit none
  !Arguments
  integer :: nd1,nd2,nd3,i1,i2,i3
  integer :: ind
  !Local variables
  integer :: a1,a2,a3
  if ( i1 == nd1 ) then
     a1=1
  else
     a1=i1
  end if
  if ( i2 == nd2 ) then
     a2=1
  else
     a2=i2
  end if
  if ( i3 == nd3 ) then
     a3=1
  else
     a3=i3
  end if
  ind=a1+nd1*(a2-1)+nd1*nd2*(a3-1)
end subroutine norm_ind
!!***


!!****f* BigDFT/symm_ind
!! NAME
!!   symm_ind
!!
!! FUNCTION
!!   Index in zarray for -g vector
!!
!! SOURCE
!!
subroutine symm_ind(nd1,nd2,nd3,i1,i2,i3,ind)
  implicit none
  !Arguments
  integer :: nd1,nd2,nd3,i1,i2,i3
  integer :: ind
  !Local variables
  integer ::  a1,a2,a3
  if (i1 /= 1) then 
     a1=nd1+1-i1
  else
     a1=i1
  end if
  if (i2 /= 1) then 
     a2=nd2+1-i2
  else
     a2=i2
  end if
  if (i3 /= 1) then 
     a3=nd3+1-i3
  else
     a3=i3
  end if
  ind=a1+nd1*(a2-1)+nd1*nd2*(a3-1)
end subroutine symm_ind
!!***


!!****f* BigDFT/symm_ind3
!! NAME
!!   symm_ind3
!!
!! FUNCTION
!!   From index of g, index of -g
!!
!! SOURCE
!!
subroutine symm_ind3(nd1,nd2,nd3,i1,i2,i3,a1,a2,a3)
  implicit none
  !Arguments
  integer :: nd1,nd2,nd3,i1,i2,i3
  integer ::a1,a2,a3
  !Local variables
  if (i1 /= 1) then 
     a1=nd1+1-i1
  else
     a1=i1
  end if
  if (i2 /= 1) then 
     a2=nd2+1-i2
  else
     a2=i2
  end if
  if (i3 /= 1) then 
     a3=nd3+1-i3
  else
     a3=i3
  end if
end subroutine symm_ind3
!!***


!!****f* BigDFT/zarray_in
!! NAME
!!   zarray_in
!!
!! FUNCTION
!!   Put the density into zarray
!!
!! SOURCE
!!
subroutine zarray_in(n01,n02,n03,nd1,nd2,nd3,density,zarray)
   implicit none
   !Arguments
   integer :: n01,n02,n03,nd1,nd2,nd3
   real*8, dimension(n01,n02,n03) :: density
   real*8, dimension(2,nd1,nd2,nd3) :: zarray
   !Local variables
   integer :: i1,i2,i3,n01h,nd1hm,nd3hm,nd2hm
   !Half the size of n01
   n01h=n01/2
   nd1hm=(nd1-1)/2
   nd2hm=(nd2-1)/2
   nd3hm=(nd3-1)/2
   !Set to zero
   do i3=1,nd3
   do i2=1,nd2
   do i1=1,nd1
   zarray(1,i1,i2,i3) = 0.d0
   zarray(2,i1,i2,i3) = 0.d0
   enddo ; enddo ; enddo
   !Set zarray
   do i3=1,n03
     do i2=1,n02
       do i1=1,n01h
         zarray(1,i1+nd1hm,i2+nd2hm,i3+nd3hm) = density(2*i1-1,i2,i3)
         zarray(2,i1+nd1hm,i2+nd2hm,i3+nd3hm) = density(2*i1,i2,i3)
      end do
    end do
  end do
  if(modulo(n01,2) == 1) then
      do i3=1,n03
         do i2=1,n02
            zarray(1,n01h+1+nd1hm,i2+nd2hm,i3+nd3hm) = density(n01,i2,i3)
         end do
      end do
  end if
end subroutine zarray_in
!!***


!!****f* BigDFT/zarray_out
!! NAME
!!   zarray_out
!!
!! FUNCTION
!!   Set the potential (rhopot) from zarray
!!   Calculate the Hartree energy.
!!
!! SOURCE
!!
subroutine zarray_out(n01,n02,n03,nd1,nd2,nd3,&
     rhopot,zarray,factor,hgrid,ehartree)
  implicit none
  !Arguments
  integer :: n01,n02,n03,nd1,nd2,nd3
  real*8, dimension(n01,n02,n03) :: rhopot
  !Convert zarray(2,nd1,nd2,nd3) -> zarray(2*nd1,nd2,nd3)
  !to use i1=1,n01 instead of i1=1,n1h + special case for modulo(n01,2)
  real*8, dimension(2*nd1,nd2,nd3) :: zarray
  real*8 :: factor,hgrid
  real*8 :: ehartree
  !Local variables
  real*8 :: pot1
  integer :: i1,i2,i3
  !
  ehartree=0.d0
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           pot1 = factor*zarray(i1,i2,i3)
           ehartree = ehartree + pot1 * rhopot(i1,i2,i3)
           rhopot(i1,i2,i3) = pot1
        end do
     end do
  end do
  !Double counting and integration step
  ehartree=0.5d0*ehartree*hgrid**3
end subroutine zarray_out
!!***


!!****f* BigDFT/excpotu
!! NAME
!!   excpotu
!!
!! FUNCTION
!!   Takes the Hartree potential as contained in pot_ion from the
!!   Poisson solver and sums the Exc and Hartree potential in the array rhopot
!!   Calculates also the XC energy eexcu and XC potential energy vexcu
!!
!! SYNOPSIS
!!  rhopot(n01,n02,n03)   (inout) input density, output potential
!!  pot_ion(n01,n02,n03)  (in)    ionic potential
!!  zarray(nd1_2,nd2,nd3) (in)    nd1_2=2*nd1
!!        Use nd1_2 for the call of SPSolver_Kernel
!!
!! AUTHOR
!!    S. Goedecker
!!
!! SOURCE
!!
subroutine excpotu(n01,n02,n03,nd1_2,nd2,nd3,&
     rhopot,pot_ion,zarray,factor,hgrid,ehartree,eexcu,vexcu)
  implicit none
  !Arguments
  integer :: n01,n02,n03,nd1_2,nd2,nd3
  real*8, dimension(n01,n02,n03) :: rhopot
  real*8, dimension(n01,n02,n03) :: pot_ion
  !Convert zarray(2,nd1,nd2,nd3) -> zarray(2*nd1,nd2,nd3)
  !to use i1=1,n01 instead of i1=1,n1h + special case for modulo(n1,2)
  real*8, dimension(nd1_2,nd2,nd3) :: zarray
  real*8 :: factor,hgrid
  real*8 :: ehartree,eexcu,vexcu
  !Local variables
  real*8, parameter :: &
       a0u=.4581652932831429d0, &
       a1u=2.217058676663745d0, & 
       a2u=0.7405551735357053d0,&
       a3u=0.01968227878617998d0
  real*8, parameter :: &
       b1u=1.0d0, &
       b2u=4.504130959426697d0, & 
       b3u=1.110667363742916d0, &
       b4u=0.02359291751427506d0
  real*8, parameter :: & 
       c1u=4.d0*a0u*b1u/3.0d0, &
       c2u=5.0d0*a0u*b2u/3.0d0+a1u*b1u, & 
       c3u=2.0d0*a0u*b3u+4.0d0*a1u*b2u/3.0d0+2.0d0*a2u*b1u/3.0d0,  & 
       c4u=7.0d0*a0u*b4u/3.0d0+5.0d0*a1u*b3u/3.0d0+a2u*b2u+a3u*b1u/3.0d0,  & 
       c5u=2.0d0*a1u*b4u+4.0d0*a2u*b3u/3.0d0+2.0d0*a3u*b2u/3.0d0,  & 
       c6u=5.0d0*a2u*b4u/3.0d0+a3u*b3u, &
       c7u=4.0d0*a3u*b4u/3.0d0
  real*8, parameter :: rsfac=.6203504908994000d0
  ! real*8, parameter :: eps2=1.d-28
  real*8, parameter :: thirdm=-1.d0/3.d0
  integer :: i3,i2,i1
  real*8 :: rhou1,pot1,rsu1,topu1,dtopu1,botu1,t1
  real*8 :: epsxcu1,p1
  !Body
  eexcu=0.d0
  vexcu=0.d0
  ehartree=0.d0
  !	x1=1.d0
  !	ic=0
  do i3=1,n03
     do i2=1,n02
	do i1=1,n01
           rhou1=rhopot(i1,i2,i3)
           if (rhou1 < 1.d-20) then
              pot1=factor*zarray(i1,i2,i3)
              ehartree=ehartree+rhou1*pot1
              rhopot(i1,i2,i3)=pot1+pot_ion(i1,i2,i3)
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
              pot1=factor*zarray(i1,i2,i3)
              p1=pot_ion(i1,i2,i3)
              ehartree=ehartree+rhou1*pot1
              rhopot(i1,i2,i3)=pot1+(dtopu1*t1*t1)+p1
              
           endif
           
        end do
     end do
  end do
  eexcu=eexcu*hgrid**3
  vexcu=vexcu*hgrid**3
  ehartree=0.5d0*ehartree*hgrid**3
  write(6,*) 'ehartree,eexcu,vexcu',ehartree,eexcu,vexcu
  !write(6,*) 'average iterations for root in excpotu',2.d0*ic/(n01*n02*n03)
end subroutine excpotu
!!***


!!****f* BigDFT/check_symmetry
!! NAME
!!   check_symmetry
!!
!! FUNCTION
!!   Check the symmetry of zarray
!!
!! SOURCE
!!
subroutine check_symmetry(nd1,nd2,nd3,zarray,inzee)
  implicit none
  !Arguments
  real*8, dimension(2,nd1*nd2*nd3,2) :: zarray
  !Local variables 
  integer :: i1,i2,i3,nd1,nd2,nd3,ind1,ind2,inzee,f1,f2,f3
  f1=nd1
  f2=nd2
  f3=nd3
  print *,"Checking proper symmetry..."
  do i3=1,f3
     do i2=1,f2
        do i1=1,f1
           call norm_ind(nd1,nd2,nd3,i1,i2,i3,ind1)
           call symm_ind(nd1,nd2,nd3,i1,i2,i3,ind2)
           if(abs(zarray(1,ind1,inzee)-zarray(1,ind2,inzee)) <= 1d-10 .and. &
                abs(zarray(2,ind1,inzee)+zarray(2,ind2,inzee)) <= 1d-10 ) then
           else
              print *,"no symmetry -> reality",&
                   i1,i2,i3,nd1,nd2,nd3,&
                   zarray(1,ind1,inzee),zarray(1,ind2,inzee) 
              stop
           end if
        end do
     end do
  end do
  print *,"...ok."
end subroutine check_symmetry
!!***


!!****f* BigDFT/test_kernel
!! NAME
!!   test_kernel
!!
!! FUNCTION
!!   Test the kernel
!!
!! SOURCE
!!
subroutine test_kernel(n01,n02,n03,nfft1,nfft2,nfft3,&
     hgrid,karray,pot_ion,rhopot)
  implicit none
  !Arguments
  integer :: n01,n02,n03,nfft1,nfft2,nfft3
  real*8 :: hgrid
  real*8, dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1) :: karray
  real*8, dimension(n01,n02,n03) :: pot_ion
  real*8, dimension(n01,n02,n03) :: rhopot
  !Local variables
  real*8 :: a_gauss,a2
  real*8 :: rhotot,shft1,shft2,shft3,ehart,eexcu,vexcu
  real*8 :: pi,x1,x2,x3,r,r2,factor,derf,max_diff,diff,tt
  integer :: i1,i2,i3,ii1,ii2,ii3
  
  a_gauss=4.d0*hgrid
  a2 = a_gauss**2

  write(*,*) 'test_kernel, dim kernel',nfft1/2+1,nfft2/2+1,nfft3/2+1
  
  !Shift the center of the Gaussian 
  !away from central grid point to break symmetries
  !shft1=1.3d0*hgrid
  !shft2=0.5d0*hgrid
  !shft3=0.1d0*hgrid
  shft1=0.d0
  shft2=0.d0
  shft3=0.d0
  
  !Initialisation
  pi = 4.d0*atan(1.d0)
  !Normalisation
  factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
  !Gaussian function
  rhotot=0.d0
  do i3=1,n03
     x3 = hgrid*(i3-n03/2)-shft3
     do i2=1,n02
        x2 = hgrid*(i2-n02/2)-shft2
        do i1=1,n01
           x1 = hgrid*(i1-n01/2)-shft1
           r2 = x1*x1+x2*x2+x3*x3
           rhopot(i1,i2,i3) = factor*exp(-r2/a2)
           rhotot=rhotot+rhopot(i1,i2,i3)
        end do
     end do
  end do
  rhotot=rhotot*hgrid**3
  
  !! Plot values along x axis
  !   open(unit=11,file='rho.dat')
  !   do i3=1,n03
  !      x1 =                            - shft1
  !      x2 =                            - shft2
  !      x3 = hgrid*(i3-n03/2) - shft3
  !      r=sqrt(x1**2+x2**2+x3**2)
  !      write(unit=11,fmt="(e10.3,e12.5,2(e21.14),e9.2,2(e12.5))") &
  !           r, rhopot(n01/2,n02/2,i3)
  !   end do
  !   close(unit=11)
  
   !Calculate potential using Poisson Solver
      write(*,*) 'testing poisson solver'
      call PSolver_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,&
           hgrid,karray,0,pot_ion,rhopot,ehart,eexcu,vexcu)
   
   !! Plot values along x axis
   !   open(unit=11,file='pot.dat')
   !   do i3=1,n03
   !      x1 =                            - shft1
   !      x2 =                            - shft2
   !      x3 = hgrid*(i3-n03/2) - shft3
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
   !           r, rhopot(n01/2,n02/2,i3),factor,tt
   !   end do
   !   close(unit=11)
   
   
   
   ! Global error 
   max_diff = 0.d0
   do i3=1,n03
      x3 = hgrid*(i3-n03/2) - shft3
      do i2=1,n02
         x2 = hgrid*(i2-n02/2) - shft2
         do i1=1,n01
            x1 = hgrid*(i1-n01/2) - shft1
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
   write(*,'(1x,a,f7.2,1x,e10.3,1x,e10.3)') &
        'hgridh,Deltarho,max_diff',hgrid,rhotot-1.d0,max_diff
   write(*,*) 'Max diff at : ',ii1,ii2,ii3
   write(*,*) 'Poisson Solver test finished'

end subroutine test_kernel
!!***

