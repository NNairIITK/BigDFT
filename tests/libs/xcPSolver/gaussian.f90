!> @file
!!  Use integral form for Poisson solver
!! @author
!!    Copyright (c) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!> Program testing new ideas like momentum-preserving gaussian integrals
program MP_gaussian
  use module_base
  use gaussians
  use yaml_output
  implicit none
  integer, parameter :: nmoms=16,nsteps=1
  integer :: npts,j,imoms,pow
  real(gp) :: hgrid,pgauss,x0,reference
  real(gp), dimension(0:nmoms,2) :: moments
  real(gp), dimension(3,2,0:nmoms) :: avgmaxmin
  real(gp), dimension(:), allocatable :: fj_phi,fj_coll,fj_phi2
  

  call f_lib_initialize()

  !number of points of the array
  npts=1000
  hgrid=1.0_gp
  pgauss=8.0e-3_dp*1.25_dp**(6*(6-1))
  x0=0.1_gp
  pow=1

  call yaml_map('sigma/h',sqrt(0.5_gp/pgauss)/hgrid)
  !plot function
  do j=-npts,npts
     if (pow /= 0) then
        write(16,*)j,exp(-pgauss*(j*hgrid-x0)**2)*((j*hgrid-x0)**pow)
     else
        write(16,*)j,exp(-pgauss*(j*hgrid-x0)**2)
     end if
  end do

  avgmaxmin=0.d0
  avgmaxmin(3,:,:)=1.d100

  !array where we have to write the value of the discretization
  fj_phi=f_malloc(-npts .to. npts,id='fj_phi')
  fj_coll=f_malloc(-npts .to. npts,id='fj_coll')

  call initialize_real_space_conversion() !initialize the work arrays needed to integrate with isf

  call evaluate_moments(nmoms,npts,hgrid,pgauss,pow,x0,fj_phi,fj_coll,moments)
!!$  !print moments value
!!$  do imoms=0,nmoms
!!$     reference=gauint0(pgauss,imoms+pow)
!!$     if (reference /=0.0_gp) then
!!$        call yaml_map('Mom No.'//trim(yaml_toa(imoms)),&
!!$             (moments(imoms,:)-reference)/reference,fmt='(1pe22.14)',advance='no')
!!$     else
!!$        call yaml_map('Mom No.'//trim(yaml_toa(imoms)),&
!!$             moments(imoms,:),fmt='(1pe22.14)',advance='no')
!!$     end if
!!$     call yaml_comment('Ref: '//trim(yaml_toa(reference,fmt='(1pe22.14)')))
!!$  end do

  !calculate average, maximum and minimum
  do j=1,2
     do imoms=0,nmoms
        avgmaxmin(1,j,imoms)=avgmaxmin(1,j,imoms)+moments(imoms,j)/real(nsteps,gp)
        avgmaxmin(2,j,imoms)=max(moments(imoms,j),avgmaxmin(2,j,imoms))
        avgmaxmin(3,j,imoms)=min(moments(imoms,j),avgmaxmin(3,j,imoms))
     end do
  end do
  

  call yaml_map('Results',reshape(avgmaxmin,(/6,nmoms+1/)))

  call finalize_real_space_conversion('Main program')
  
  call f_free(fj_phi,fj_coll)
  call f_lib_finalize()
end program MP_gaussian

!> classify the quality of a multipole extraction in both cases
subroutine evaluate_moments(nmoms,npts,hgrid,pgauss,pow,x0,fj_phi,fj_coll,moments)
  use module_base, only: gp
  use gaussians, only: scfdotf
  implicit none
  integer, intent(in) :: npts,pow,nmoms
  real(gp), intent(in) :: hgrid,pgauss,x0
  real(gp), dimension(0:nmoms,2), intent(out) :: moments
  real(gp), dimension(-npts:npts), intent(out) :: fj_phi,fj_coll
  !local variables
  integer :: j

  !use the elemental property of the scfdotf function
  fj_phi=scfdotf((/(j,j=-npts,npts)/),hgrid,pgauss,x0,pow)
  call moments_1d(2*npts+1,fj_phi,x0+hgrid*(npts+1),hgrid,nmoms,moments(0,1))

  !collocation array
  if (pow /=0) then
     fj_coll=(/(exp(-pgauss*(j*hgrid-x0)**2)*(j*hgrid-x0)**pow,j=-npts,npts)/)
  else
     fj_coll=(/(exp(-pgauss*(j*hgrid-x0)**2),j=-npts,npts)/)
  end if
  call moments_1d(2*npts+1,fj_coll,x0+hgrid*(npts+1),hgrid,nmoms,moments(0,2))

end subroutine evaluate_moments

!>calculate the moments of a array with respect to a reference point 
subroutine moments_1d(n,array,x0,h,nmoms,moments)
  use module_base, only:gp
  implicit none
  integer, intent(in) :: nmoms,n
  real(gp), intent(in) :: x0,h !<grid spacing
  real(gp), dimension(n), intent(in) :: array
  real(gp), dimension(0:nmoms), intent(out) :: moments
  !local variables
  integer :: j,k
  real(gp) :: x

  do j=0,nmoms
     moments(j)=0.0_gp
     do k=1,n
        x=real(k,gp)*h-x0
        moments(j)=moments(j)+x**j*array(k)
     end do
     moments(j)=moments(j)*h
  end do

end subroutine moments_1d

!!$
!!$!> takes a Gaussian of exponent pgauss an center x0 and discretize it on a grid of size 1 in units of sqrt(0.5*[pgauss])
!!$!! f(x)=fac*exp(-pgauss*(x-x0)**2)
!!$subroutine discretize_gaussian(nrange,fac,pgauss,x0,hgrid,filename)
!!$  use module_base
!!$  use gaussians
!!$  !not yet use dynamic_memory !let's start
!!$  implicit none
!!$  integer, intent(in) :: nrange
!!$  real(gp), intent(in) :: fac,pgauss,x0,hgrid
!!$  character(len=*), intent(in) :: filename
!!$  !local variables
!!$  integer, parameter :: itype_scf=16,npoints=2**6,nmoms=6
!!$  integer :: k,n_scf,n_range,j
!!$  real(gp) :: dx,x
!!$  real(gp), dimension(0:nmoms,3) :: moments !<to check all methods
!!$  real(dp), dimension(0:2048) :: fISF
!!$  real(gp), dimension(:), allocatable :: f_i,f_phi_i,f_phi_i2,x_scf,y_scf
!!$
!!$  !allocate the three gaussians arrays
!!$  allocate(f_i(-nrange:nrange))
!!$  allocate(f_phi_i(-nrange:nrange))
!!$  allocate(f_phi_i2(-nrange:nrange))
!!$
!!$  !Number of integration points = 2*itype_scf*n_points
!!$  n_scf=2*itype_scf*npoints
!!$
!!$  !Other allocations
!!$  allocate(x_scf(0:n_scf))
!!$  allocate(y_scf(0:n_scf))
!!$
!!$  !Build the scaling function
!!$  call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
!!$
!!$  !Step grid for the integration
!!$  dx = real(nrange,gp)/real(n_scf,gp)
!!$
!!$  !first, collocate the gaussian
!!$  do k=-nrange,nrange
!!$     f_i(k)=fac*exp(-pgauss*(hgrid*real(k,gp)-x0)**2)
!!$  end do
!!$
!!$  f_phi_i2=0.0_gp
!!$  f_phi_i=0.0_gp
!!$  !then, calculate its values integrated with SCF
!!$  !use the i2 array as a workspace
!!$  call my_gauss_conv_scf(itype_scf,pgauss,x0,hgrid,dx,nrange,n_scf,x_scf,y_scf,f_phi_i,f_phi_i2)
!!$  call dscal(2*nrange+1,fac,f_phi_i,1)
!!$  !Only itype=8,14,16,20,24,30,40,50,60,100
!!$  call four_isf_coeffs(itype_scf,fISF)
!!$  f_phi_i2=0.0_gp
!!$  !then calculate the values via the analytic method (for the moment only on the second half)
!!$  !this routine is not useful for performances of precision and speed
!!$  !call my_analytic_integral(hgrid*sqrt(pgauss),x0/hgrid,&
!!$  !     nrange,itype_scf,f_phi_i2,fISF,64)
!!$  !  do k=-nrange,-1
!!$  !     f_phi_i2(k)=f_phi_i2(-k)
!!$  !  end do
!!$  !call dscal(2*nrange+1,fac,f_phi_i2,1)
!!$
!!$  !use daubechies wavelets for expressing the function
!!$  !  call gauss_to_daub(hgrid,fac,x0,sqrt(0.5/pgauss),0,&!no err, errsuc
!!$  !       2*nrange+1,n_left,n_right,c,err_norm,&         !no err_wav. nmax instead of n_intvx
!!$  !       ww,nwork,.false.)                             !added work arrays ww with dimension nwork
!!$  !create the isf array in the gaussian module
!!$  call initialize_real_space_conversion()
!!$  do j=-nrange,nrange
!!$     f_phi_i2(j)=scfdotf(j,hgrid,pgauss,x0,0)
!!$     !print *,'j',j,f_phi_i2(j)
!!$  end do
!!$  !use the elemental property of the scfdotf function
!!$  f_phi_i2=scfdotf((/(j,j=-nrange,nrange)/),hgrid,pgauss,x0,0)
!!$  call finalize_real_space_conversion('discretize_gaussian')
!!$
!!$  !then print the results
!!$  moments=0.0_dp
!!$  !open(unit=200, file = trim(filename), status = 'unknown')
!!$  do j=0,nmoms
!!$     open(unit=200+j)
!!$     do k=-nrange,nrange
!!$        x=real(k,gp)*hgrid-x0
!!$        moments(j,1)=moments(j,1)+x**j*f_i(k)
!!$        moments(j,2)=moments(j,2)+x**j*f_phi_i(k)
!!$        moments(j,3)=moments(j,3)+x**j*f_phi_i2(k)
!!$        write(200+j,'(i4,3(1pe26.17e3))')k,x**j*f_i(k),x**j*f_phi_i(k),x**j*f_phi_i2(k)
!!$     end do
!!$     moments(j,:) = moments(j,:)*hgrid
!!$     print '(a,i3,8(1pe14.5))',trim(filename),j,&
!!$          pgauss,moments(j,:),gauint0(pgauss,j),moments(j,:)-gauint0(pgauss,j)
!!$     close(unit=200+j)
!!$  end do
!!$
!!$  deallocate(f_i,f_phi_i,f_phi_i2)
!!$
!!$end subroutine discretize_gaussian
!!$
!!$
!!$!> Write a routine which performs the integration of a function on a given grid
!!$subroutine my_gauss_conv_scf(itype_scf,pgauss,x0,hgrid,dx,n_range,n_scf,x_scf,y_scf,kernel_scf,work)
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: n_range,n_scf,itype_scf
!!$  real(dp), intent(in) :: pgauss,hgrid,dx,x0
!!$  real(dp), dimension(0:n_scf), intent(in) :: x_scf
!!$  real(dp), dimension(0:n_scf), intent(in) :: y_scf
!!$  real(dp), dimension(-n_range:n_range), intent(inout) :: work
!!$  real(dp), dimension(-n_range:n_range), intent(inout) :: kernel_scf
!!$  !local variables
!!$  real(dp), parameter :: p0_ref = 1.0_dp
!!$  integer :: n_iter,i_kern,i
!!$  real(dp) :: p0_cell,p0gauss,absci,kern
!!$
!!$  !Step grid for the integration
!!$  !dx = real(n_range,dp)/real(n_scf,dp)
!!$
!!$  !To have a correct integration
!!$  p0_cell = p0_ref/(hgrid*hgrid)
!!$
!!$  !write(*,*) 'p0_cell = ', p0_cell 
!!$
!!$  !We calculate the number of iterations to go from pgauss to p0_ref
!!$  n_iter = nint((log(pgauss) - log(p0_cell))/log(4.0_dp))
!!$
!!$  write(*,*) 'n_iter = ', n_iter
!!$
!!$  if (n_iter <= 0 .or. x0 /= 0.d0)then
!!$     n_iter = 0
!!$     p0gauss = pgauss
!!$  else
!!$     p0gauss = pgauss/4._dp**n_iter
!!$  end if
!!$
!!$  !  write(*,*) 'p0gauss =',  p0gauss
!!$
!!$  !Stupid integration
!!$  !Do the integration with the exponential centered in i_kern
!!$  kernel_scf(:) = 0.0_dp
!!$  do i_kern=-n_range,n_range
!!$     kern = 0.0_dp
!!$     do i=0,n_scf
!!$        absci = x_scf(i) - real(i_kern,dp)
!!$        absci = x0-absci*hgrid !sign has changed here
!!$        absci = absci*absci
!!$        kern = kern + y_scf(i)*dexp(-p0gauss*absci)
!!$     end do
!!$     kernel_scf(-i_kern) = kern*dx
!!$     !write(17,*)i_kern,kernel_scf(i_kern)
!!$     !kernel_scf(-i_kern) = kern*dx
!!$     !if (abs(kern) < 1.d-18) then
!!$     !   !Too small not useful to calculate
!!$     !   exit
!!$     !end if
!!$  end do
!!$
!!$  !  do i=0,n_scf
!!$  !     write(18,*)i,x_scf(i),y_scf(i)
!!$  !  end do
!!$
!!$  !Start the iteration to go from p0gauss to pgauss
!!$  call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,work)
!!$
!!$END SUBROUTINE my_gauss_conv_scf
!!$
!!$
!!$!> Here alpha corresponds to sqrt(alpha) in mathematica
!!$!! the final result is fwork(j+m)-fwork(j-m)
!!$subroutine my_analytic_integral(alpha,x0,ntot,m,fwork,fISF,argument_nf)
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: ntot,m
!!$  real(dp), intent(in) :: alpha,x0 !<x0 is the deviation wrt the grid spacing
!!$  real(dp), dimension(-ntot:ntot), intent(inout) :: fwork
!!$  real(dp), dimension(0:2048), intent(in) :: fISF
!!$  integer, intent(in) :: argument_nf 
!!$
!!$  !local variables
!!$  integer :: nf
!!$  real(dp), parameter :: pi=3.1415926535897932384_dp
!!$  logical :: flag,flag1,flag2
!!$  integer :: j,q,jz
!!$  real(dp) :: if,r1,r2,res,ypm,ymm,erfcpm,erfcmm,factor,re,ro,factorend
!!$
!!$  !write(*,*) fISF(1000), fISF16(1000)
!!$
!!$  !  if(present(argument_nf)) then 
!!$  nf=argument_nf
!!$  !  else
!!$  !     nf = 64 ! "default value"
!!$  !  endif
!!$
!!$  flag=.false.
!!$  factor=pi/real(2*m,dp)/alpha
!!$  factorend=sqrt(pi)/alpha/real(4*m,dp)*exp((alpha*x0)**2)
!!$  !fill work array
!!$  !the calculation for j=0 can be separated from the rest
!!$  !since it only contains error functions
!!$  !the values of the loop are not symmetric in principle
!!$  loop_nonzero: do j=-ntot,ntot
!!$     ypm=alpha*(real(j+m,dp)+x0)
!!$     ymm=alpha*(real(j-m,dp)+x0)
!!$     call derfcf(erfcpm,ypm) ! erfcpm = erfc(ypm)
!!$     call derfcf(erfcmm,ymm)
!!$
!!$     !assume nf even
!!$     res=0._dp
!!$     !reso=0._dp
!!$     do q=nf,2,-2
!!$        !the sign of q only influences the imaginary part
!!$        !so we can multiply by a factor of two
!!$
       ! call wofz_mod(alpha,m,q,-j-m,r1,if,flag1)
       ! call wofz_mod(alpha,m,q,-j+m,r2,if,flag2)
!!$        call GcplxInt(alpha,m,q,-j-m,-x0,r1,if,flag1)
       ! if (q==12 .and. j==15) then
       !    print *,'here',alpha,m,q,j,x0,r1,if,flag1
       !    stop
       ! end if
!!$        call GcplxInt(alpha,m,q,-j+m,-x0,r2,if,flag2)
!!$        re=r1-r2
!!$        flag=flag1 .or. flag2 .or. flag
!!$        !if (flag) then
!!$        !   print *,'here',r1,if,q,j,x0
!!$        !   stop 
!!$        !end if
       ! call wofz_mod(alpha,m,q-1,-j-m,r1,if,flag1)
       ! call wofz_mod(alpha,m,q-1,-j+m,r2,if,flag2)
!!$        call GcplxInt(alpha,m,q-1,-j-m,-x0,r1,if,flag1)
!!$        call GcplxInt(alpha,m,q-1,-j+m,-x0,r2,if,flag2)
!!$        ro=r1-r2
!!$        flag=flag1 .or. flag2 .or. flag
!!$        !if (flag) then
!!$        !   print *,'there',xo,y
!!$        !   stop 
!!$        !end if
!!$        !write(16,'(2(i4),6(1pe15.7))')j,q,re,ro,erfcmm-erfcpm
!!$        re=re*fISF(q)
!!$        ro=ro*fISF(q-1)
!!$        res=res+re-ro
!!$     end do
!!$     !q=0 
!!$     !fwork(j)=derf(y)+rese-reso
!!$     fwork(j)=erfcmm-erfcpm+2.0_dp*res!e-reso
!!$     fwork(j)=factorend*fwork(j)
!!$     !exit form the loop if it is close to zero
!!$     !     if (abs(fwork(j)) < 1.e-25_dp) exit loop_nonzero
!!$     !write(17,'(i4,8(1pe15.7))')j,derf(y),erfcsgn,rese,reso,derf(y)+rese-reso,&
!!$     !     -erfcsgn+rese-reso,erfcsgn+rese-reso
!!$  end do loop_nonzero
!!$
!!$  !check flag
!!$  if (flag) then
!!$     write (*,*)'value of alpha',alpha
!!$     stop 'problem occurred in wofz'
!!$  end if
!!$
!!$  !put to zero the rest
!!$  do jz=j+1,ntot
!!$     fwork(jz)=0.0_dp
!!$  end do
!!$
!!$END SUBROUTINE my_analytic_integral
!!$
