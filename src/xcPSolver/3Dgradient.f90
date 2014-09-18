!> @file
!! Routines to calculate gradients of the electronic densities for GGA
!! @author
!! Copyright (C) 2002-2011 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!> Calculates the White-Bird correction to the XC potential.
!! Since this correction strongly depends on the way of calculating the gradient
!! It is based on a finite difference calculation of the gradient, corrected at the border
!! (see the calc_gradient routine)
!! Works either in parallel or in serial, by proper adjustation of the arguments
!!
!! INPUTS and OUTPUT
!! @param f_i(n1,n2,n3,3,ispden) Three functions (depending on the fourth index) indicating the 
!!           derivative of rho*e_xc with respect to the three components of the gradient.
!!           Its dimensions are properly enlarged to allow another finite-difference calculation
!! @param wbl, wbr indicates the point starting from which we have to calculate the WB correction
!!           The correction is then added to the wb_vxc array in the proper points
!! @param n3grad is the effective dimension of interest
subroutine wb_correction(geocode,n1,n2,n3,n3grad,wbl,wbr,f_i,hx,hy,hz,nspden,&
     wb_vxc)
  use Poisson_Solver, only: dp
 implicit none
 !Arguments
 character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
 integer, intent(in) :: n1,n2,n3,n3grad,wbl,wbr,nspden
 real(dp), intent(in) :: hx,hy,hz
 real(dp), dimension(n1,n2,n3,3,nspden), intent(in) :: f_i
 real(dp), dimension(n1,n2,n3,nspden), intent(inout) :: wb_vxc
 !Local variables
 integer :: i1,i2,i3,ispden
 !Filters of finite difference derivative for order 4
 real(dp), parameter :: a1=0.8d0, a2=-0.2d0
 real(dp), parameter :: a3=0.038095238095238095238d0, a4=-0.0035714285714285714286d0
 real(dp) :: derx,dery,derz,c1,c2,c3,c4
 !Body

 !Some check
 if ((n3 /= n3grad + wbl + wbr -2)) then
    print *,'wb_correction:incompatibility of the dimensions, n3=',n3,'n3grad=',n3grad,&
         'wbl=',wbl,'wbr=',wbr
    stop
 end if
 if (geocode /='F' .and. wbl /= wbr ) then
    print *,'wb_correction:incompatibility of the dimensions, n3=',n3,'n3grad=',n3grad,&
         'wbl=',wbl,'wbr=',wbr,'geocode',geocode
    stop
 end if


 !coefficients to calculate the wb correction to the gradient
 c4=a4
 c3=c4+a3
 c2=c3+a2
 c1=c2+a1

 !loop over the spin components
 do ispden=1,nspden

    !loop over the different directions

    !x direction
    if (geocode == 'F') then
       do i3=wbl,n3grad+wbl-1
          do i2=1,n2
             derx=-c1*f_i(1,i2,i3,1,ispden)&
                  -c1*f_i(2,i2,i3,1,ispden)&
                  -c2*f_i(3,i2,i3,1,ispden)&
                  -c3*f_i(4,i2,i3,1,ispden)&
                  -c4*f_i(5,i2,i3,1,ispden)
             wb_vxc(1,i2,i3,ispden)=wb_vxc(1,i2,i3,ispden)+derx/hx

             derx=&
                  a1*f_i(1,i2,i3,1,ispden)&
                  -a1*f_i(3,i2,i3,1,ispden)-a2*f_i(4,i2,i3,1,ispden)&
                  -a3*f_i(5,i2,i3,1,ispden)-a4*f_i(6,i2,i3,1,ispden)
             wb_vxc(2,i2,i3,ispden)=wb_vxc(2,i2,i3,ispden)+derx/hx

             derx=&
                  a2*f_i(1,i2,i3,1,ispden)+a1*f_i(2,i2,i3,1,ispden)&
                  -a1*f_i(4,i2,i3,1,ispden)-a2*f_i(5,i2,i3,1,ispden)&
                  -a3*f_i(6,i2,i3,1,ispden)-a4*f_i(7,i2,i3,1,ispden)
             wb_vxc(3,i2,i3,ispden)=wb_vxc(3,i2,i3,ispden)+derx/hx

             derx=a3*f_i(1,i2,i3,1,ispden)&
                  +a2*f_i(2,i2,i3,1,ispden)+a1*f_i(3,i2,i3,1,ispden)&
                  -a1*f_i(5,i2,i3,1,ispden)-a2*f_i(6,i2,i3,1,ispden)&
                  -a3*f_i(7,i2,i3,1,ispden)-a4*f_i(8,i2,i3,1,ispden)
             wb_vxc(4,i2,i3,ispden)=wb_vxc(4,i2,i3,ispden)+derx/hx

             do i1=5,n1-4
                derx=-a1*(f_i(i1+1,i2,i3,1,ispden)-f_i(i1-1,i2,i3,1,ispden))&
                     -a2*(f_i(i1+2,i2,i3,1,ispden)-f_i(i1-2,i2,i3,1,ispden))&
                     -a3*(f_i(i1+3,i2,i3,1,ispden)-f_i(i1-3,i2,i3,1,ispden))&
                     -a4*(f_i(i1+4,i2,i3,1,ispden)-f_i(i1-4,i2,i3,1,ispden))
                wb_vxc(i1,i2,i3,ispden)=wb_vxc(i1,i2,i3,ispden)+derx/hx
             end do

             derx=-a1*(f_i(n1-2,i2,i3,1,ispden)-f_i(n1-4,i2,i3,1,ispden))&
                  -a2*(f_i(n1-1,i2,i3,1,ispden)-f_i(n1-5,i2,i3,1,ispden))&
                  -a3*(f_i(n1,i2,i3,1,ispden)-f_i(n1-6,i2,i3,1,ispden))&
                  -a4*(-f_i(n1-7,i2,i3,1,ispden))
             wb_vxc(n1-3,i2,i3,ispden)=wb_vxc(n1-3,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(n1-1,i2,i3,1,ispden)-f_i(n1-3,i2,i3,1,ispden))&
                  -a2*(f_i(n1,i2,i3,1,ispden)-f_i(n1-4,i2,i3,1,ispden))&
                  -a3*(-f_i(n1-5,i2,i3,1,ispden))&
                  -a4*(-f_i(n1-6,i2,i3,1,ispden))
             wb_vxc(n1-2,i2,i3,ispden)=wb_vxc(n1-2,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(n1,i2,i3,1,ispden)-f_i(n1-2,i2,i3,1,ispden))&
                  -a2*(-f_i(n1-3,i2,i3,1,ispden))&
                  -a3*(-f_i(n1-4,i2,i3,1,ispden))&
                  -a4*(-f_i(n1-5,i2,i3,1,ispden))
             wb_vxc(n1-1,i2,i3,ispden)=wb_vxc(n1-1,i2,i3,ispden)+derx/hx

             derx= c1*f_i(n1,i2,i3,1,ispden)&
                  +c1*f_i(n1-1,i2,i3,1,ispden)&
                  +c2*f_i(n1-2,i2,i3,1,ispden)&
                  +c3*f_i(n1-3,i2,i3,1,ispden)&
                  +c4*f_i(n1-4,i2,i3,1,ispden)
             wb_vxc(n1,i2,i3,ispden)=wb_vxc(n1,i2,i3,ispden)+derx/hx
          end do
       end do
    else !surface or periodic case is the same for x direction
       do i3=wbl,n3grad+wbl-1
          do i2=1,n2
             derx=-a1*(f_i(2,i2,i3,1,ispden)-f_i(n1,i2,i3,1,ispden))&
                  -a2*(f_i(3,i2,i3,1,ispden)-f_i(n1-1,i2,i3,1,ispden))&
                  -a3*(f_i(4,i2,i3,1,ispden)-f_i(n1-2,i2,i3,1,ispden))&
                  -a4*(f_i(5,i2,i3,1,ispden)-f_i(n1-3,i2,i3,1,ispden))
             wb_vxc(1,i2,i3,ispden)=wb_vxc(1,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(3,i2,i3,1,ispden)-f_i(1,i2,i3,1,ispden))&
                  -a2*(f_i(4,i2,i3,1,ispden)-f_i(n1,i2,i3,1,ispden))&
                  -a3*(f_i(5,i2,i3,1,ispden)-f_i(n1-1,i2,i3,1,ispden))&
                  -a4*(f_i(6,i2,i3,1,ispden)-f_i(n1-2,i2,i3,1,ispden))
             wb_vxc(2,i2,i3,ispden)=wb_vxc(2,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(4,i2,i3,1,ispden)-f_i(2,i2,i3,1,ispden))&
                  -a2*(f_i(5,i2,i3,1,ispden)-f_i(1,i2,i3,1,ispden))&
                  -a3*(f_i(6,i2,i3,1,ispden)-f_i(n1,i2,i3,1,ispden))&
                  -a4*(f_i(7,i2,i3,1,ispden)-f_i(n1-1,i2,i3,1,ispden))
             wb_vxc(3,i2,i3,ispden)=wb_vxc(3,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(5,i2,i3,1,ispden)-f_i(3,i2,i3,1,ispden))&
                  -a2*(f_i(6,i2,i3,1,ispden)-f_i(2,i2,i3,1,ispden))&
                  -a3*(f_i(7,i2,i3,1,ispden)-f_i(1,i2,i3,1,ispden))&
                  -a4*(f_i(8,i2,i3,1,ispden)-f_i(n1,i2,i3,1,ispden))
             wb_vxc(4,i2,i3,ispden)=wb_vxc(4,i2,i3,ispden)+derx/hx

             do i1=5,n1-4
                derx=-a1*(f_i(i1+1,i2,i3,1,ispden)-f_i(i1-1,i2,i3,1,ispden))&
                     -a2*(f_i(i1+2,i2,i3,1,ispden)-f_i(i1-2,i2,i3,1,ispden))&
                     -a3*(f_i(i1+3,i2,i3,1,ispden)-f_i(i1-3,i2,i3,1,ispden))&
                     -a4*(f_i(i1+4,i2,i3,1,ispden)-f_i(i1-4,i2,i3,1,ispden))
                wb_vxc(i1,i2,i3,ispden)=wb_vxc(i1,i2,i3,ispden)+derx/hx
             end do

             derx=-a1*(f_i(n1-2,i2,i3,1,ispden)-f_i(n1-4,i2,i3,1,ispden))&
                  -a2*(f_i(n1-1,i2,i3,1,ispden)-f_i(n1-5,i2,i3,1,ispden))&
                  -a3*(f_i(n1,i2,i3,1,ispden)-f_i(n1-6,i2,i3,1,ispden))&
                  -a4*(f_i(1,i2,i3,1,ispden)-f_i(n1-7,i2,i3,1,ispden))
             wb_vxc(n1-3,i2,i3,ispden)=wb_vxc(n1-3,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(n1-1,i2,i3,1,ispden)-f_i(n1-3,i2,i3,1,ispden))&
                  -a2*(f_i(n1,i2,i3,1,ispden)-f_i(n1-4,i2,i3,1,ispden))&
                  -a3*(f_i(1,i2,i3,1,ispden)-f_i(n1-5,i2,i3,1,ispden))&
                  -a4*(f_i(2,i2,i3,1,ispden)-f_i(n1-6,i2,i3,1,ispden))
             wb_vxc(n1-2,i2,i3,ispden)=wb_vxc(n1-2,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(n1,i2,i3,1,ispden)-f_i(n1-2,i2,i3,1,ispden))&
                  -a2*(f_i(1,i2,i3,1,ispden)-f_i(n1-3,i2,i3,1,ispden))&
                  -a3*(f_i(2,i2,i3,1,ispden)-f_i(n1-4,i2,i3,1,ispden))&
                  -a4*(f_i(3,i2,i3,1,ispden)-f_i(n1-5,i2,i3,1,ispden))
             wb_vxc(n1-1,i2,i3,ispden)=wb_vxc(n1-1,i2,i3,ispden)+derx/hx

             derx=-a1*(f_i(1,i2,i3,1,ispden)-f_i(n1-1,i2,i3,1,ispden))&
                  -a2*(f_i(2,i2,i3,1,ispden)-f_i(n1-2,i2,i3,1,ispden))&
                  -a3*(f_i(3,i2,i3,1,ispden)-f_i(n1-3,i2,i3,1,ispden))&
                  -a4*(f_i(4,i2,i3,1,ispden)-f_i(n1-4,i2,i3,1,ispden))
             wb_vxc(n1,i2,i3,ispden)=wb_vxc(n1,i2,i3,ispden)+derx/hx
          end do
       end do
    end if

    !y direction
    if (geocode /= 'P') then
       !only the periodic case is periodic in y
       do i3=wbl,n3grad+wbl-1
          do i1=1,n1
             dery=-c1*f_i(i1,1,i3,2,ispden)&
                  -c1*f_i(i1,2,i3,2,ispden)&
                  -c2*f_i(i1,3,i3,2,ispden)&
                  -c3*f_i(i1,4,i3,2,ispden)&
                  -c4*f_i(i1,5,i3,2,ispden)
             wb_vxc(i1,1,i3,ispden)=wb_vxc(i1,1,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=&
                  a1*f_i(i1,1,i3,2,ispden)&
                  -a1*f_i(i1,3,i3,2,ispden)-a2*f_i(i1,4,i3,2,ispden)&
                  -a3*f_i(i1,5,i3,2,ispden)-a4*f_i(i1,6,i3,2,ispden)
             wb_vxc(i1,2,i3,ispden)=wb_vxc(i1,2,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=&
                  a2*f_i(i1,1,i3,2,ispden)+a1*f_i(i1,2,i3,2,ispden)&
                  -a1*f_i(i1,4,i3,2,ispden)-a2*f_i(i1,5,i3,2,ispden)&
                  -a3*f_i(i1,6,i3,2,ispden)-a4*f_i(i1,7,i3,2,ispden)
             wb_vxc(i1,3,i3,ispden)=wb_vxc(i1,3,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=a3*f_i(i1,1,i3,2,ispden)&
                  +a2*f_i(i1,2,i3,2,ispden)+a1*f_i(i1,3,i3,2,ispden)&
                  -a1*f_i(i1,5,i3,2,ispden)-a2*f_i(i1,6,i3,2,ispden)&
                  -a3*f_i(i1,7,i3,2,ispden)-a4*f_i(i1,8,i3,2,ispden)
             wb_vxc(i1,4,i3,ispden)=wb_vxc(i1,4,i3,ispden)+dery/hy
          end do
          do i2=5,n2-4
             do i1=1,n1
                dery=-a1*(f_i(i1,i2+1,i3,2,ispden)-f_i(i1,i2-1,i3,2,ispden))&
                     -a2*(f_i(i1,i2+2,i3,2,ispden)-f_i(i1,i2-2,i3,2,ispden))&
                     -a3*(f_i(i1,i2+3,i3,2,ispden)-f_i(i1,i2-3,i3,2,ispden))&
                     -a4*(f_i(i1,i2+4,i3,2,ispden)-f_i(i1,i2-4,i3,2,ispden))
                wb_vxc(i1,i2,i3,ispden)=wb_vxc(i1,i2,i3,ispden)+dery/hy
             end do
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,n2-2,i3,2,ispden)-f_i(i1,n2-4,i3,2,ispden))&
                  -a2*(f_i(i1,n2-1,i3,2,ispden)-f_i(i1,n2-5,i3,2,ispden))&
                  -a3*(f_i(i1,n2,i3,2,ispden)-f_i(i1,n2-6,i3,2,ispden))&
                  -a4*(-f_i(i1,n2-7,i3,2,ispden))
             wb_vxc(i1,n2-3,i3,ispden)=wb_vxc(i1,n2-3,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,n2-1,i3,2,ispden)-f_i(i1,n2-3,i3,2,ispden))&
                  -a2*(f_i(i1,n2,i3,2,ispden)-f_i(i1,n2-4,i3,2,ispden))&
                  -a3*(-f_i(i1,n2-5,i3,2,ispden))&
                  -a4*(-f_i(i1,n2-6,i3,2,ispden))
             wb_vxc(i1,n2-2,i3,ispden)=wb_vxc(i1,n2-2,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,n2,i3,2,ispden)-f_i(i1,n2-2,i3,2,ispden))&
                  -a2*(-f_i(i1,n2-3,i3,2,ispden))&
                  -a3*(-f_i(i1,n2-4,i3,2,ispden))&
                  -a4*(-f_i(i1,n2-5,i3,2,ispden))
             wb_vxc(i1,n2-1,i3,ispden)=wb_vxc(i1,n2-1,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery= c1*f_i(i1,n2,i3,2,ispden)&
                  +c1*f_i(i1,n2-1,i3,2,ispden)&
                  +c2*f_i(i1,n2-2,i3,2,ispden)&
                  +c3*f_i(i1,n2-3,i3,2,ispden)&
                  +c4*f_i(i1,n2-4,i3,2,ispden)
             wb_vxc(i1,n2,i3,ispden)=wb_vxc(i1,n2,i3,ispden)+dery/hy
          end do
       end do
    else
       do i3=wbl,n3grad+wbl-1
          do i1=1,n1
             dery=-a1*(f_i(i1,2,i3,2,ispden)-f_i(i1,n2,i3,2,ispden))&
                  -a2*(f_i(i1,3,i3,2,ispden)-f_i(i1,n2-1,i3,2,ispden))&
                  -a3*(f_i(i1,4,i3,2,ispden)-f_i(i1,n2-2,i3,2,ispden))&
                  -a4*(f_i(i1,5,i3,2,ispden)-f_i(i1,n2-3,i3,2,ispden))
             wb_vxc(i1,1,i3,ispden)=wb_vxc(i1,1,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,3,i3,2,ispden)-f_i(i1,1,i3,2,ispden))&
                  -a2*(f_i(i1,4,i3,2,ispden)-f_i(i1,n2,i3,2,ispden))&
                  -a3*(f_i(i1,5,i3,2,ispden)-f_i(i1,n2-1,i3,2,ispden))&
                  -a4*(f_i(i1,6,i3,2,ispden)-f_i(i1,n2-2,i3,2,ispden))
             wb_vxc(i1,2,i3,ispden)=wb_vxc(i1,2,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,4,i3,2,ispden)-f_i(i1,2,i3,2,ispden))&
                  -a2*(f_i(i1,5,i3,2,ispden)-f_i(i1,1,i3,2,ispden))&
                  -a3*(f_i(i1,6,i3,2,ispden)-f_i(i1,n2,i3,2,ispden))&
                  -a4*(f_i(i1,7,i3,2,ispden)-f_i(i1,n2-1,i3,2,ispden))
             wb_vxc(i1,3,i3,ispden)=wb_vxc(i1,3,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,5,i3,2,ispden)-f_i(i1,3,i3,2,ispden))&
                  -a2*(f_i(i1,6,i3,2,ispden)-f_i(i1,2,i3,2,ispden))&
                  -a3*(f_i(i1,7,i3,2,ispden)-f_i(i1,1,i3,2,ispden))&
                  -a4*(f_i(i1,8,i3,2,ispden)-f_i(i1,n2,i3,2,ispden))
             wb_vxc(i1,4,i3,ispden)=wb_vxc(i1,4,i3,ispden)+dery/hy
          end do
          do i2=5,n2-4
             do i1=1,n1
                dery=-a1*(f_i(i1,i2+1,i3,2,ispden)-f_i(i1,i2-1,i3,2,ispden))&
                     -a2*(f_i(i1,i2+2,i3,2,ispden)-f_i(i1,i2-2,i3,2,ispden))&
                     -a3*(f_i(i1,i2+3,i3,2,ispden)-f_i(i1,i2-3,i3,2,ispden))&
                     -a4*(f_i(i1,i2+4,i3,2,ispden)-f_i(i1,i2-4,i3,2,ispden))
                wb_vxc(i1,i2,i3,ispden)=wb_vxc(i1,i2,i3,ispden)+dery/hy
             end do
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,n2-2,i3,2,ispden)-f_i(i1,n2-4,i3,2,ispden))&
                  -a2*(f_i(i1,n2-1,i3,2,ispden)-f_i(i1,n2-5,i3,2,ispden))&
                  -a3*(f_i(i1,n2,i3,2,ispden)-f_i(i1,n2-6,i3,2,ispden))&
                  -a4*(f_i(i1,1,i3,2,ispden)-f_i(i1,n2-7,i3,2,ispden))
             wb_vxc(i1,n2-3,i3,ispden)=wb_vxc(i1,n2-3,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,n2-1,i3,2,ispden)-f_i(i1,n2-3,i3,2,ispden))&
                  -a2*(f_i(i1,n2,i3,2,ispden)-f_i(i1,n2-4,i3,2,ispden))&
                  -a3*(f_i(i1,1,i3,2,ispden)-f_i(i1,n2-5,i3,2,ispden))&
                  -a4*(f_i(i1,2,i3,2,ispden)-f_i(i1,n2-6,i3,2,ispden))
             wb_vxc(i1,n2-2,i3,ispden)=wb_vxc(i1,n2-2,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,n2,i3,2,ispden)-f_i(i1,n2-2,i3,2,ispden))&
                  -a2*(f_i(i1,1,i3,2,ispden)-f_i(i1,n2-3,i3,2,ispden))&
                  -a3*(f_i(i1,2,i3,2,ispden)-f_i(i1,n2-4,i3,2,ispden))&
                  -a4*(f_i(i1,3,i3,2,ispden)-f_i(i1,n2-5,i3,2,ispden))
             wb_vxc(i1,n2-1,i3,ispden)=wb_vxc(i1,n2-1,i3,ispden)+dery/hy
          end do
          do i1=1,n1
             dery=-a1*(f_i(i1,1,i3,2,ispden)-f_i(i1,n2-1,i3,2,ispden))&
                  -a2*(f_i(i1,2,i3,2,ispden)-f_i(i1,n2-2,i3,2,ispden))&
                  -a3*(f_i(i1,3,i3,2,ispden)-f_i(i1,n2-3,i3,2,ispden))&
                  -a4*(f_i(i1,4,i3,2,ispden)-f_i(i1,n2-4,i3,2,ispden))
             wb_vxc(i1,n2,i3,ispden)=wb_vxc(i1,n2,i3,ispden)+dery/hy
          end do
       end do
    end if

    !z direction
    if (geocode == 'F') then
       if(wbl <= 1) then
          do i2=1,n2
             do i1=1,n1
                derz=-c1*f_i(i1,i2,1,3,ispden)&
                     -c1*f_i(i1,i2,2,3,ispden)&
                     -c2*f_i(i1,i2,3,3,ispden)&
                     -c3*f_i(i1,i2,4,3,ispden)&
                     -c4*f_i(i1,i2,5,3,ispden)
                wb_vxc(i1,i2,1,ispden)=wb_vxc(i1,i2,1,ispden)+derz/hz
             end do
          end do
       end if
       if(wbl <= 2 .and. n3 > 5 ) then
          do i2=1,n2
             do i1=1,n1
                derz=&
                     a1*f_i(i1,i2,1,3,ispden)&
                     -a1*f_i(i1,i2,3,3,ispden)-a2*f_i(i1,i2,4,3,ispden)&
                     -a3*f_i(i1,i2,5,3,ispden)-a4*f_i(i1,i2,6,3,ispden)
                wb_vxc(i1,i2,2,ispden)=wb_vxc(i1,i2,2,ispden)+derz/hz
             end do
          end do
       end if
       if(wbl <= 3 .and. n3 > 6 ) then
          do i2=1,n2
             do i1=1,n1
                derz=&
                     a2*f_i(i1,i2,1,3,ispden)+a1*f_i(i1,i2,2,3,ispden)&
                     -a1*f_i(i1,i2,4,3,ispden)-a2*f_i(i1,i2,5,3,ispden)&
                     -a3*f_i(i1,i2,6,3,ispden)-a4*f_i(i1,i2,7,3,ispden)
                wb_vxc(i1,i2,3,ispden)=wb_vxc(i1,i2,3,ispden)+derz/hz
             end do
          end do
       end if
       if(wbl <= 4 .and. n3 > 7 ) then
          do i2=1,n2
             do i1=1,n1
                derz=a3*f_i(i1,i2,1,3,ispden)&
                     +a2*f_i(i1,i2,2,3,ispden)+a1*f_i(i1,i2,3,3,ispden)&
                     -a1*f_i(i1,i2,5,3,ispden)-a2*f_i(i1,i2,6,3,ispden)&
                     -a3*f_i(i1,i2,7,3,ispden)-a4*f_i(i1,i2,8,3,ispden)
                wb_vxc(i1,i2,4,ispden)=wb_vxc(i1,i2,4,ispden)+derz/hz
             end do
          end do
       end if
       do i3=5,n3-4
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,i3+1,3,ispden)-f_i(i1,i2,i3-1,3,ispden))&
                     -a2*(f_i(i1,i2,i3+2,3,ispden)-f_i(i1,i2,i3-2,3,ispden))&
                     -a3*(f_i(i1,i2,i3+3,3,ispden)-f_i(i1,i2,i3-3,3,ispden))&
                     -a4*(f_i(i1,i2,i3+4,3,ispden)-f_i(i1,i2,i3-4,3,ispden))
                wb_vxc(i1,i2,i3,ispden)=wb_vxc(i1,i2,i3,ispden)+derz/hz
             end do
          end do
       end do
       if (wbr <=4 .and. n3 > 7 ) then
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,n3-2,3,ispden)-f_i(i1,i2,n3-4,3,ispden))&
                     -a2*(f_i(i1,i2,n3-1,3,ispden)-f_i(i1,i2,n3-5,3,ispden))&
                     -a3*(f_i(i1,i2,n3,3,ispden)-f_i(i1,i2,n3-6,3,ispden))&
                     -a4*(-f_i(i1,i2,n3-7,3,ispden))
                wb_vxc(i1,i2,n3-3,ispden)=wb_vxc(i1,i2,n3-3,ispden)+derz/hz
             end do
          end do
       end if
       if(wbr <= 3 .and. n3 > 6) then
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,n3-1,3,ispden)-f_i(i1,i2,n3-3,3,ispden))&
                     -a2*(f_i(i1,i2,n3,3,ispden)-f_i(i1,i2,n3-4,3,ispden))&
                     -a3*(-f_i(i1,i2,n3-5,3,ispden))&
                     -a4*(-f_i(i1,i2,n3-6,3,ispden))
                wb_vxc(i1,i2,n3-2,ispden)=wb_vxc(i1,i2,n3-2,ispden)+derz/hz
             end do
          end do
       end if
       if(wbr <= 2 .and. n3 > 5) then
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,n3,3,ispden)-f_i(i1,i2,n3-2,3,ispden))&
                     -a2*(-f_i(i1,i2,n3-3,3,ispden))&
                     -a3*(-f_i(i1,i2,n3-4,3,ispden))&
                     -a4*(-f_i(i1,i2,n3-5,3,ispden))
                wb_vxc(i1,i2,n3-1,ispden)=wb_vxc(i1,i2,n3-1,ispden)+derz/hz
             end do
          end do
       end if
       if(wbr <= 1) then
          do i2=1,n2
             do i1=1,n1
                derz= c1*f_i(i1,i2,n3,3,ispden)&
                     +c1*f_i(i1,i2,n3-1,3,ispden)&
                     +c2*f_i(i1,i2,n3-2,3,ispden)&
                     +c3*f_i(i1,i2,n3-3,3,ispden)&
                     +c4*f_i(i1,i2,n3-4,3,ispden)
                wb_vxc(i1,i2,n3,ispden)=wb_vxc(i1,i2,n3,ispden)+derz/hz
             end do
          end do
       end if
    else
       if(wbl == 1) then
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,2,3,ispden)-f_i(i1,i2,n3,3,ispden))&
                     -a2*(f_i(i1,i2,3,3,ispden)-f_i(i1,i2,n3-1,3,ispden))&
                     -a3*(f_i(i1,i2,4,3,ispden)-f_i(i1,i2,n3-2,3,ispden))&
                     -a4*(f_i(i1,i2,5,3,ispden)-f_i(i1,i2,n3-3,3,ispden))
                wb_vxc(i1,i2,1,ispden)=wb_vxc(i1,i2,1,ispden)+derz/hz
             end do
          end do
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,3,3,ispden)-f_i(i1,i2,1,3,ispden))&
                     -a2*(f_i(i1,i2,4,3,ispden)-f_i(i1,i2,n3,3,ispden))&
                     -a3*(f_i(i1,i2,5,3,ispden)-f_i(i1,i2,n3-1,3,ispden))&
                     -a4*(f_i(i1,i2,6,3,ispden)-f_i(i1,i2,n3-2,3,ispden))
                wb_vxc(i1,i2,2,ispden)=wb_vxc(i1,i2,2,ispden)+derz/hz
             end do
          end do
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,4,3,ispden)-f_i(i1,i2,2,3,ispden))&
                     -a2*(f_i(i1,i2,5,3,ispden)-f_i(i1,i2,1,3,ispden))&
                     -a3*(f_i(i1,i2,6,3,ispden)-f_i(i1,i2,n3,3,ispden))&
                     -a4*(f_i(i1,i2,7,3,ispden)-f_i(i1,i2,n3-1,3,ispden))
                wb_vxc(i1,i2,3,ispden)=wb_vxc(i1,i2,3,ispden)+derz/hz
             end do
          end do
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,5,3,ispden)-f_i(i1,i2,3,3,ispden))&
                     -a2*(f_i(i1,i2,6,3,ispden)-f_i(i1,i2,2,3,ispden))&
                     -a3*(f_i(i1,i2,7,3,ispden)-f_i(i1,i2,1,3,ispden))&
                     -a4*(f_i(i1,i2,8,3,ispden)-f_i(i1,i2,n3,3,ispden))
                wb_vxc(i1,i2,4,ispden)=wb_vxc(i1,i2,4,ispden)+derz/hz
             end do
          end do
!!!          i3s=5
!!!          i3e=n3-4
       else
!!!          i3s=1
!!!          i3e=n3
       end if
       do i3=5,n3-4 
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,i3+1,3,ispden)-f_i(i1,i2,i3-1,3,ispden))&
                     -a2*(f_i(i1,i2,i3+2,3,ispden)-f_i(i1,i2,i3-2,3,ispden))&
                     -a3*(f_i(i1,i2,i3+3,3,ispden)-f_i(i1,i2,i3-3,3,ispden))&
                     -a4*(f_i(i1,i2,i3+4,3,ispden)-f_i(i1,i2,i3-4,3,ispden))
                wb_vxc(i1,i2,i3,ispden)=wb_vxc(i1,i2,i3,ispden)+derz/hz
             end do
          end do
       end do
       !wbr or wbl is the same in periodic BC
       if (wbr == 1) then
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,n3-2,3,ispden)-f_i(i1,i2,n3-4,3,ispden))&
                     -a2*(f_i(i1,i2,n3-1,3,ispden)-f_i(i1,i2,n3-5,3,ispden))&
                     -a3*(f_i(i1,i2,n3,3,ispden)-f_i(i1,i2,n3-6,3,ispden))&
                     -a4*(f_i(i1,i2,1,3,ispden)-f_i(i1,i2,n3-7,3,ispden))
                wb_vxc(i1,i2,n3-3,ispden)=wb_vxc(i1,i2,n3-3,ispden)+derz/hz
             end do
          end do
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,n3-1,3,ispden)-f_i(i1,i2,n3-3,3,ispden))&
                     -a2*(f_i(i1,i2,n3,3,ispden)-f_i(i1,i2,n3-4,3,ispden))&
                     -a3*(f_i(i1,i2,1,3,ispden)-f_i(i1,i2,n3-5,3,ispden))&
                     -a4*(f_i(i1,i2,2,3,ispden)-f_i(i1,i2,n3-6,3,ispden))
                wb_vxc(i1,i2,n3-2,ispden)=wb_vxc(i1,i2,n3-2,ispden)+derz/hz
             end do
          end do
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,n3,3,ispden)-f_i(i1,i2,n3-2,3,ispden))&
                     -a2*(f_i(i1,i2,1,3,ispden)-f_i(i1,i2,n3-3,3,ispden))&
                     -a3*(f_i(i1,i2,2,3,ispden)-f_i(i1,i2,n3-4,3,ispden))&
                     -a4*(f_i(i1,i2,3,3,ispden)-f_i(i1,i2,n3-5,3,ispden))
                wb_vxc(i1,i2,n3-1,ispden)=wb_vxc(i1,i2,n3-1,ispden)+derz/hz
             end do
          end do
          do i2=1,n2
             do i1=1,n1
                derz=-a1*(f_i(i1,i2,1,3,ispden)-f_i(i1,i2,n3-1,3,ispden))&
                     -a2*(f_i(i1,i2,2,3,ispden)-f_i(i1,i2,n3-2,3,ispden))&
                     -a3*(f_i(i1,i2,3,3,ispden)-f_i(i1,i2,n3-3,3,ispden))&
                     -a4*(f_i(i1,i2,4,3,ispden)-f_i(i1,i2,n3-4,3,ispden))
                wb_vxc(i1,i2,n3,ispden)=wb_vxc(i1,i2,n3,ispden)+derz/hz
             end do
          end do
       end if
    end if
    !end of the loop over spin components
 end do

END SUBROUTINE wb_correction


!> Calculates the finite difference gradient White-Bird correction to the XC potential.
!! The gradient in point x is calculated by taking four point before and after x.
!! The lack of points near the border is solved by ideally prolungating the input
!! function outside the grid, by assigning the same value of the border point.
!! To implement this operation an auxiliary array, bigger than the original grid, is defined.
!! This routines works either in parallel or in serial by proper adjustation of the arguments.
!!
!! INPUTS 
!! @param rhoinp(n1,n2,n3,nspden) input function
!!           the dimension of rhoinp are more than needed because in the parallel case we must
!!           take into account the points before and after the selected interval, 
!!           of effective (third) dimension n3grad.
!! @param hx,hy,hz grid spacing in the three directions
!! 
!! OUTPUT
!! @param gradient(n1,n2,n3grad,2*nspden-1,0:3) the square modulus of the gradient (index 0) 
!!           and the three components.
!!           In the spin-polarized case,
!!           Following ABINIT conventions:\n
!!           - gradient(i1,i2,i3,1,0:3) is the gradient of the spin up density
!!           - gradient(i1,i2,i3,2,0:3) is the gradient of the spin down density
!!           - gradient(i1,i2,i3,3,0:3) is the gradient of the total density
!!           of course gradient(i1,i2,i3,3,0) /= gradient(i1,i2,i3,2,0) + gradient(i1,i2,i3,2,0)
!!           since the modulus is non-linear
!!           The drivexc routine uses only the first part of this array (component 0)
!!           The rest of the array is defined for later use in the wb postprocessing routine.
subroutine calc_gradient(geocode,n1,n2,n3,n3grad,deltaleft,deltaright,rhoinp,nspden,hx,hy,hz,&
     gradient,rhocore)
  use Poisson_Solver, only: dp, wp
  use memory_profiling
  use wrapper_linalg
  use dynamic_memory
 implicit none
 !Arguments
 character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
 integer, intent(in) :: n1,n2,n3,n3grad,deltaleft,deltaright,nspden
 real(dp), intent(in) :: hx,hy,hz
 real(dp), dimension(n1,n2,n3,nspden), intent(inout) :: rhoinp
 real(dp), dimension(n1,n2,n3grad,2*nspden-1,0:3), intent(out) :: gradient
 real(dp), dimension(:,:,:,:), pointer :: rhocore
 !Local variables
 character(len=*), parameter :: subname='calc_gradient'
 integer :: i1,i2,i3,j3,ispden
 !filters of finite difference derivative for order 4
 real(dp), parameter :: a1=0.8d0, a2=-0.2d0
 real(dp), parameter :: a3=0.038095238095238095238d0, a4=-0.0035714285714285714286d0
 real(dp) :: derx,dery,derz,drcx,drcy,drcz
 real(dp), dimension(:,:,:), allocatable :: density
 !Body

 !some check
 if (n3 /= n3grad + deltaleft + deltaright) then
    print *,'calc_gradient:incompatibility of the dimensions, n3=',n3,'n3grad=',n3grad,&
         'deltaleft=',deltaleft,'deltaright=',deltaright
    stop
 end if

 if (geocode /= 'F' .and. deltaleft/=deltaright) then
    print *,'calc_gradient:incompatibility of the dimensions, n3=',n3,'n3grad=',n3grad,&
         'deltaleft=',deltaleft,'deltaright=',deltaright,'geocode=',geocode
    stop
 end if

 !let us initialize the larger vector to calculate the gradient
 density = f_malloc((/ n1+8 , n2+8 , n3grad+8/),id='density')

 do ispden=1,nspden !loop over up/dw densities
    select case(geocode)
    case('F')
       do i3=1,4-deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3)=rhoinp(i1-4,i2-4,1,ispden)
             end do
          end do
       end do
       do i3=1,deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
       do i3=deltaleft+1,n3grad+deltaleft
          do i2=1,4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,1,i3,ispden)
             end do
          end do
          do i2=5,n2+4
             do i1=1,4
                density(i1,i2,i3+4-deltaleft)=rhoinp(1,i2-4,i3,ispden)
             end do
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
             do i1=n1+5,n1+8
                density(i1,i2,i3+4-deltaleft)=rhoinp(n1,i2-4,i3,ispden)
             end do
          end do
          do i2=n2+5,n2+8
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,n2,i3,ispden)
             end do
          end do
       end do
       do i3=n3grad+deltaleft+1,n3
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
       do i3=1,4-deltaright
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+n3+4-deltaleft)=rhoinp(i1-4,i2-4,n3,ispden)
             end do
          end do
       end do
    case('P')
       !periodic BC, deltaleft/right are zero only if the array is not cut
       do i3=1,4-deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3)=rhoinp(i1-4,i2-4,n3-4+i3,ispden)
             end do
          end do
       end do
       !normal case, the array rhoinp contains all the necessary
       do i3=1,deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
       do i3=deltaleft+1,n3grad+deltaleft
          do i2=1,4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,n2-4+i2,i3,ispden)
             end do
          end do
          do i2=5,n2+4
             do i1=1,4
                density(i1,i2,i3+4-deltaleft)=rhoinp(n1-4+i1,i2-4,i3,ispden)
             end do
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
             do i1=n1+5,n1+8
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-n1-4,i2-4,i3,ispden)
             end do
          end do
          do i2=n2+5,n2+8
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-n2-4,i3,ispden)
             end do
          end do
       end do
       do i3=n3grad+deltaleft+1,n3
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
       do i3=1,4-deltaright
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+n3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
    case('S')
       !surface case, periodic in x and z while isolated for y
       do i3=1,4-deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3)=rhoinp(i1-4,i2-4,n3-4+i3,ispden)
             end do
          end do
       end do
       do i3=1,deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
       do i3=deltaleft+1,n3grad+deltaleft
          do i2=1,4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,1,i3,ispden)
             end do
          end do
          do i2=5,n2+4
             do i1=1,4
                density(i1,i2,i3+4-deltaleft)=rhoinp(n1-4+i1,i2-4,i3,ispden)
             end do
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
             do i1=n1+5,n1+8
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-n1-4,i2-4,i3,ispden)
             end do
          end do
          do i2=n2+5,n2+8
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,n2,i3,ispden)
             end do
          end do
       end do
       do i3=n3grad+deltaleft+1,n3
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
       do i3=1,4-deltaright
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+n3+4-deltaleft)=rhoinp(i1-4,i2-4,i3,ispden)
             end do
          end do
       end do
    case default
       write(*,*)' Wrong geocode for gradient calculation:',geocode
       stop
    end select

    if (associated(rhocore)) then
       !calculating the gradient by using the auxiliary array
       !$omp parallel do default(shared) private(i3,i2,i1,derx,dery,derz,drcx,drcy,drcz)
       do i3=5,n3grad+4 
          do i2=5,n2+4
             do i1=5,n1+4
                !gradient in the x direction
                derx=a1*(density(i1+1,i2,i3)-density(i1-1,i2,i3))&
                     +a2*(density(i1+2,i2,i3)-density(i1-2,i2,i3))&
                     +a3*(density(i1+3,i2,i3)-density(i1-3,i2,i3))&
                     +a4*(density(i1+4,i2,i3)-density(i1-4,i2,i3))
                !gradient in the y direction
                dery=a1*(density(i1,i2+1,i3)-density(i1,i2-1,i3))&
                     +a2*(density(i1,i2+2,i3)-density(i1,i2-2,i3))&
                     +a3*(density(i1,i2+3,i3)-density(i1,i2-3,i3))&
                     +a4*(density(i1,i2+4,i3)-density(i1,i2-4,i3))
                !gradient in the z direction
                derz=a1*(density(i1,i2,i3+1)-density(i1,i2,i3-1))&
                     +a2*(density(i1,i2,i3+2)-density(i1,i2,i3-2))&
                     +a3*(density(i1,i2,i3+3)-density(i1,i2,i3-3))&
                     +a4*(density(i1,i2,i3+4)-density(i1,i2,i3-4))
                !rhocore gradients
                drcx=rhocore(i1-4,i2-4,i3-4+deltaleft,2)
                drcy=rhocore(i1-4,i2-4,i3-4+deltaleft,3)
                drcz=rhocore(i1-4,i2-4,i3-4+deltaleft,4)
                !square modulus
                gradient(i1-4,i2-4,i3-4,ispden,0)=(derx/hx)**2+(dery/hy)**2+(derz/hz)**2+&
                     2.0_dp*((derx/hx)*drcx+(dery/hy)*drcy+(derz/hz)*drcz)+&
                     drcx**2+drcy**2+drcz**2
                !different components
                gradient(i1-4,i2-4,i3-4,ispden,1)=derx/hx+drcx
                gradient(i1-4,i2-4,i3-4,ispden,2)=dery/hy+drcy
                gradient(i1-4,i2-4,i3-4,ispden,3)=derz/hz+drcz
             end do
          end do
       end do
       !$end parallel do
    else
       !$omp parallel do default(shared) private(i3,i2,i1,derx,dery,derz)
       do i3=5,n3grad+4 
          do i2=5,n2+4
             do i1=5,n1+4
                !gradient in the x direction
                derx=a1*(density(i1+1,i2,i3)-density(i1-1,i2,i3))&
                     +a2*(density(i1+2,i2,i3)-density(i1-2,i2,i3))&
                     +a3*(density(i1+3,i2,i3)-density(i1-3,i2,i3))&
                     +a4*(density(i1+4,i2,i3)-density(i1-4,i2,i3))
                !gradient in the y direction
                dery=a1*(density(i1,i2+1,i3)-density(i1,i2-1,i3))&
                     +a2*(density(i1,i2+2,i3)-density(i1,i2-2,i3))&
                     +a3*(density(i1,i2+3,i3)-density(i1,i2-3,i3))&
                     +a4*(density(i1,i2+4,i3)-density(i1,i2-4,i3))
                !gradient in the z direction
                derz=a1*(density(i1,i2,i3+1)-density(i1,i2,i3-1))&
                     +a2*(density(i1,i2,i3+2)-density(i1,i2,i3-2))&
                     +a3*(density(i1,i2,i3+3)-density(i1,i2,i3-3))&
                     +a4*(density(i1,i2,i3+4)-density(i1,i2,i3-4))
                !square modulus
                gradient(i1-4,i2-4,i3-4,ispden,0)=(derx/hx)**2+(dery/hy)**2+(derz/hz)**2
                !different components
                gradient(i1-4,i2-4,i3-4,ispden,1)=derx/hx
                gradient(i1-4,i2-4,i3-4,ispden,2)=dery/hy
                gradient(i1-4,i2-4,i3-4,ispden,3)=derz/hz
             end do
          end do
       end do
       !$end parallel do
    end if
 end do
 
 !Once again for total density 
 !(adding up component only since density already contains dw component from above)
 if(nspden==2) then
    select case(geocode)
    case('F')
       do i3=1,4-deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3)=density(i1,i2,i3)+rhoinp(i1-4,i2-4,1,1)
             end do
          end do
       end do
       do i3=1,deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
       do i3=deltaleft+1,n3grad+deltaleft
          do i2=1,4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,1,i3,1)
             end do
          end do
          do i2=5,n2+4
             do i1=1,4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(1,i2-4,i3,1)
             end do
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
             do i1=n1+5,n1+8
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(n1,i2-4,i3,1)
             end do
          end do
          do i2=n2+5,n2+8
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,n2,i3,1)
             end do
          end do
       end do
       do i3=n3grad+deltaleft+1,n3
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
       do i3=1,4-deltaright
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+n3+4-deltaleft)=density(i1,i2,i3+n3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,n3,1)
             end do
          end do
       end do
    case('P')
       do i3=1,4-deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3)=density(i1,i2,i3)+rhoinp(i1-4,i2-4,n3-4+i3,1)
             end do
          end do
       end do
       do i3=1,deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
       do i3=deltaleft+1,n3grad+deltaleft
          do i2=1,4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,n2-4+i2,i3,1)
             end do
          end do
          do i2=5,n2+4
             do i1=1,4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(n1-4+i1,i2-4,i3,1)
             end do
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
             do i1=n1+5,n1+8
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-n1-4,i2-4,i3,1)
             end do
          end do
          do i2=n2+5,n2+8
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-n2-4,i3,1)
             end do
          end do
       end do
       do i3=n3grad+deltaleft+1,n3
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
       do i3=1,4-deltaright
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+n3+4-deltaleft)=density(i1,i2,i3+n3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
    case('S')
       do i3=1,4-deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3)=density(i1,i2,i3)+rhoinp(i1-4,i2-4,n3-4+i3,1)
             end do
          end do
       end do
       do i3=1,deltaleft
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
       do i3=deltaleft+1,n3grad+deltaleft
          do i2=1,4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,1,i3,1)
             end do
          end do
          do i2=5,n2+4
             do i1=1,4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(n1-4+i1,i2-4,i3,1)
             end do
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
             do i1=n1+5,n1+8
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-n1-4,i2-4,i3,1)
             end do
          end do
          do i2=n2+5,n2+8
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,n2,i3,1)
             end do
          end do
       end do
       do i3=n3grad+deltaleft+1,n3
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+4-deltaleft)=density(i1,i2,i3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
       do i3=1,4-deltaright
          do i2=5,n2+4
             do i1=5,n1+4
                density(i1,i2,i3+n3+4-deltaleft)=density(i1,i2,i3+n3+4-deltaleft)+&
                     rhoinp(i1-4,i2-4,i3,1)
             end do
          end do
       end do
    end select
    
    if (associated(rhocore)) then
       !calculating the gradient by using the auxiliary array
       !$omp parallel do default(shared) private(i3,i2,i1,derx,dery,derz,drcx,drcy,drcz)
       do i3=5,n3grad+4 
          do i2=5,n2+4
             do i1=5,n1+4
                !gradient in the x direction
                derx=a1*(density(i1+1,i2,i3)-density(i1-1,i2,i3))&
                     +a2*(density(i1+2,i2,i3)-density(i1-2,i2,i3))&
                     +a3*(density(i1+3,i2,i3)-density(i1-3,i2,i3))&
                     +a4*(density(i1+4,i2,i3)-density(i1-4,i2,i3))
                !gradient in the y direction
                dery=a1*(density(i1,i2+1,i3)-density(i1,i2-1,i3))&
                     +a2*(density(i1,i2+2,i3)-density(i1,i2-2,i3))&
                     +a3*(density(i1,i2+3,i3)-density(i1,i2-3,i3))&
                     +a4*(density(i1,i2+4,i3)-density(i1,i2-4,i3))
                !gradient in the z direction
                derz=a1*(density(i1,i2,i3+1)-density(i1,i2,i3-1))&
                     +a2*(density(i1,i2,i3+2)-density(i1,i2,i3-2))&
                     +a3*(density(i1,i2,i3+3)-density(i1,i2,i3-3))&
                     +a4*(density(i1,i2,i3+4)-density(i1,i2,i3-4))
                !rhocore gradients
                drcx=rhocore(i1-4,i2-4,i3-4+deltaleft,2)
                drcy=rhocore(i1-4,i2-4,i3-4+deltaleft,3)
                drcz=rhocore(i1-4,i2-4,i3-4+deltaleft,4)

                !square modulus
                gradient(i1-4,i2-4,i3-4,3,0)=(derx/hx)**2+(dery/hy)**2+(derz/hz)**2+&
                     4.0_dp*((derx/hx)*drcx+(dery/hy)*drcy+(derz/hz)*drcz)+&
                     4.0_dp*(drcx**2+drcy**2+drcz**2)
                !different components
                gradient(i1-4,i2-4,i3-4,3,1)=derx/hx+2.0_dp*drcx
                gradient(i1-4,i2-4,i3-4,3,2)=dery/hy+2.0_dp*drcy
                gradient(i1-4,i2-4,i3-4,3,3)=derz/hz+2.0_dp*drcz
             end do
          end do
       end do
       !$end parallel do
    else
       !calculating the gradient by using the auxiliary array
       !$omp parallel do default(shared) private(i3,i2,i1,derx,dery,derz)
       do i3=5,n3grad+4 
          do i2=5,n2+4
             do i1=5,n1+4
                !gradient in the x direction
                derx=a1*(density(i1+1,i2,i3)-density(i1-1,i2,i3))&
                     +a2*(density(i1+2,i2,i3)-density(i1-2,i2,i3))&
                     +a3*(density(i1+3,i2,i3)-density(i1-3,i2,i3))&
                     +a4*(density(i1+4,i2,i3)-density(i1-4,i2,i3))
                !gradient in the y direction
                dery=a1*(density(i1,i2+1,i3)-density(i1,i2-1,i3))&
                     +a2*(density(i1,i2+2,i3)-density(i1,i2-2,i3))&
                     +a3*(density(i1,i2+3,i3)-density(i1,i2-3,i3))&
                     +a4*(density(i1,i2+4,i3)-density(i1,i2-4,i3))
                !gradient in the z direction
                derz=a1*(density(i1,i2,i3+1)-density(i1,i2,i3-1))&
                     +a2*(density(i1,i2,i3+2)-density(i1,i2,i3-2))&
                     +a3*(density(i1,i2,i3+3)-density(i1,i2,i3-3))&
                     +a4*(density(i1,i2,i3+4)-density(i1,i2,i3-4))
                !square modulus
                gradient(i1-4,i2-4,i3-4,3,0)=(derx/hx)**2+(dery/hy)**2+(derz/hz)**2
                !different components
                gradient(i1-4,i2-4,i3-4,3,1)=derx/hx
                gradient(i1-4,i2-4,i3-4,3,2)=dery/hy
                gradient(i1-4,i2-4,i3-4,3,3)=derz/hz
             end do
          end do
       end do
       !$end parallel do
    end if

    !if n3 /= n3grad (which can appear only in parallel) translate the density to be contiguous
    !such as it can be passed to the ABINIT XC routines (drivexc)
    if (n3/=n3grad) then
       j3=0
       do i3=n3grad+deltaleft+1,n3 !we have deltaright points
          j3=j3+1
          do i2=1,n2
             do i1=1,n1
                rhoinp(i1,i2,i3,1)=rhoinp(i1,i2,deltaleft+j3,2)!+0.5_dp*rhocore(1,1,deltaleft+j3,1)
             end do
          end do
       end do
       do i3=1,n3grad-deltaright
          j3=j3+1
          do i2=1,n2
             do i1=1,n1
                rhoinp(i1,i2,i3,2)=rhoinp(i1,i2,deltaleft+j3,2)!+0.5_dp*rhocore(1,1,deltaleft+j3,1)
             end do
          end do
       end do
    end if
 end if

 !add the core density 
 if (associated(rhocore)) then! too particular, should not enter here in this case .and. n1*n2*(n3grad+deltaright).gt.0) then
    call axpy(n1*n2*(n3grad+deltaleft),0.5_wp,rhocore(1,1,1,1),1,rhoinp(1,1,1,1),1)
    if (nspden==2) then
       if (n3grad /= n3) then
          if (deltaright > 0) then
             call axpy(n1*n2*(n3grad+deltaright),0.5_wp,rhocore(1,1,deltaleft+1,1),1,&
                  rhoinp(1,1,n3grad+deltaleft+1,1),1)
          else
             call axpy(n1*n2*(n3grad+deltaright),0.5_wp,rhocore(1,1,deltaleft+1,1),1,&
                  rhoinp(1,1,1,2),1)
          end if
       else
          call axpy(n1*n2*n3,0.5_wp,rhocore(1,1,1,1),1,&
               rhoinp(1,1,1,2),1)
       end if
    end if
 end if


  call f_free(density)

END SUBROUTINE calc_gradient
