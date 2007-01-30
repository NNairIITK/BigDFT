!!****h* BigDFT/PSolver
!! NAME
!!   PSolver
!!
!! FUNCTION
!!    Program test for solver Poisson
!!    Laplacian V = 4pi rho
!!    May work either in parallel or in serial case
!!
!! AUTHOR
!!    Thierry Deutsch
!!
!! COPYRIGHT
!!    Copyright (C) 2005 CEA
!! CREATION DATE
!!    06/07/2005
!!
!! MODIFICATION HISTORY
!!    15/12/2005 Change calls of routine to be compatible with Stefan's code
!!
!! SOURCE
!!
program PSolver
  implicit none
  include 'mpif.h'
  integer, parameter :: i_test=1
  !Order of interpolating scaling function
  integer, parameter :: itype_scf=8
  real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
  !Error function
  real(kind=8) :: derf
  !Length of the box
  real(kind=8), parameter :: acell = 10.d0
  character(len=50) :: name,chain
  real(kind=8), dimension(:,:,:), allocatable :: density, rhopot, karray
  real(kind=8), dimension(:,:,:), allocatable :: pot_ion
  real(kind=8) :: pi,x1,x2,x3,r,r2,factor,hgrid,max_diff
  real(kind=8) :: ehartree,eexcu,vexcu
  integer :: n01,n02,n03,m1,m2,m3,md1,md2,md3,nd1,nd2,nd3,n1,n2,n3
  integer :: i1,i2,i3,i1_max,i2_max,i3_max,n1k,n2k,n3k,iproc,nproc,ierr
  integer :: n_cell,i_allocated,l1,nsp1,nsp2,nsp3,nfft1,nfft2,nfft3
  
  !Use arguments
  call getarg(1,chain)
  read(unit=chain,fmt=*) n01
  call getarg(2,chain)
  read(unit=chain,fmt=*) n02
  call getarg(3,chain)
  read(unit=chain,fmt=*) n03
  

  print *,"PSolver: ",n01,n02,n03


  !Step size
  n_cell = min(n01,n02,n03)
  hgrid = acell/real(n_cell,kind=8)
  pi = 4.d0*atan(1.d0)
  
  !Allocations
  !Density
  allocate(density(n01,n02,n03))
  !Density then potential
  allocate(rhopot(n01,n02,n03))
  !Ionic potential (not used here)
  allocate(pot_ion(n01,n02,n03))
  pot_ion(:,:,:) = 0.d0


  !Initialisation
  if (i_test == 1) then
     !Normalisation
     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
     !gaussian function
     do i3=1,n03
        x3 = hgrid*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hgrid*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hgrid*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3
              density(i1,i2,i3) = factor*exp(-r2/a2)
           end do
        end do
     end do
  else
     !second derivative of a gaussian: (-6/a2+4r^2/a^4)exp(-r2/a4)
     factor = 1.d0/(4.d0*pi*a2)
     do i3=1,n03
        x3 = hgrid*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hgrid*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hgrid*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3
              density(i1,i2,i3) = factor*(6.d0-4.d0*r2/a2)*exp(-r2/a2)
           end do
        end do
     end do
  end if
  
  !Test PSolver for this gaussian (copy into rhopot)
  rhopot(:,:,:) = density(:,:,:)

  !Build the Kernel for full size
  !Fullsize
  call Dimensions_FFT(n01,n02,n03,nfft1,nfft2,nfft3)
  n1k=nfft1/2+1
  n2k=nfft2/2+1
  n3k=nfft3/2+1
  allocate(karray(n1k,n2k,n3k),stat=i_allocated)
  if (i_allocated /= 0) then
     print *,"Problem of memory allocation"
     stop
  end if
  call Build_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,hgrid,itype_scf,karray)
!!$  write(unit=*,fmt="(1x,3('-'),a,30('-'))") "test_kernel (fullsize -- start)"
!!$  call test_kernel(n01,n02,n03,nfft1,nfft2,nfft3,hgrid,karray,pot_ion,rhopot)
!!$  write(unit=*,fmt="(1x,3('-'),a,30('-'))") "test_kernel (fullsize -- done)"
  
  call PSolver_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,&
       hgrid,karray,.false.,pot_ion,rhopot,ehartree,eexcu,vexcu)
  !end do avoid
   deallocate(karray)

  !Maximum difference
  max_diff = 0.d0
  do i3=1,n03
     x3 = hgrid*real(i3-n03/2,kind=8)
     do i2=1,n02 
        x2 = hgrid*real(i2-n02/2,kind=8)
        do i1=1,n01
           x1 = hgrid*real(i1-n01/2,kind=8)
           r2 = x1*x1+x2*x2+x3*x3
           r = sqrt(r2)
           if (i_test == 1) then
              !Potential from a gaussian
              if (r == 0.d0) then
                 factor = abs(rhopot(i1,i2,i3)-2.d0/(sqrt(pi)*a_gauss))
              else
                 factor = abs(rhopot(i1,i2,i3)-derf(r/a_gauss)/r)
              end if
           else
              !Potential from a second derivative of a gaussian
              factor = abs(rhopot(i1,i2,i3)-exp(-r2/a2))
           end if
           if (max_diff < factor) then
              max_diff = factor
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do

  if (i_test == 1) then
     write(unit=name,fmt="(a,i0,a,i0,a,i0,a)") &
          "gaussian-",n01,"-",n02,"-",n03,".dat"
  else
     write(unit=name,fmt="(a,i0,a,i0,a,i0,a)") &
          "gaussianD-",n01,"-",n02,"-",n03,".dat"
  end if
  !print *,"datafile: ",name
  open(unit=11,file=name)
  write(unit=11,fmt="(a,f10.3,3(a,i0))") &
       "#hgrid=",hgrid," n01=",n01," n02=",n02," n03=",n03
  write(unit=11,fmt="(a)") "# x density potential derf(x)/x diff kernel"
  do i1=1,n01
     x1 = hgrid*real(i1-n01/2,kind=8)
     if (i_test == 1) then
        if (i1 == n01/2) then
           !limit_{x -> 0} erf(x/x) = 2/sqrt(pi)
           factor = 2.d0/(sqrt(pi)*a_gauss)
        else
           factor = derf(x1/a_gauss)/x1
        end if
     else
        r2 = x1*x1
        factor = exp(-r2/a2)
     end if
     write(unit=11,fmt="(6(e24.15))") &
          x1, density(i1,n02/2,n03/2), rhopot(i1,n02/2,n03/2), &
          factor, abs(rhopot(i1,n02/2,n03/2)-factor)
  end do
  close(unit=11)
!!$
  write(*,*) 'Testing Poisson Solver for a_gauss=',a_gauss
  write(unit=*,fmt="(1x,a,1pe23.15)") "ehartree=",ehartree
  write(unit=*,fmt="(1x,2(a,1pe23.15))") "eexcu=",eexcu," vexcu=",vexcu
  write(unit=*,fmt="(1x,a,2(1pe23.15))") 'hgrid, Max diff: ',hgrid,max_diff
  write(*,*) 'Max diff at : ',i1_max,i2_max,i3_max
  write(99,'(4(i4),1(1x,e12.5))')nproc,n01,n02,n03,max_diff

  !De-allocations
  deallocate(density)
  deallocate(rhopot)
  deallocate(pot_ion)

end program PSolver
!!***
