!> @file
!! atomic program for generating and optimizing hgh pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by Matthias Krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Generates logarithmic grids
 subroutine radgrid(nrad,rr,rw,rd,a_grd,b_grd,rmax)
   implicit none
   !Arguments
   integer, intent(inout) :: nrad
   real(kind=8), intent(in) :: a_grd, b_grd, rmax
   real(kind=8), dimension(nrad), intent(out) :: rr !< radial grid: rr(i)=a_grd*b_grd**(i-1)-c_grd
   real(kind=8), dimension(nrad), intent(out) :: rw !< weights for radial integration (dr/di)
   real(kind=8), dimension(nrad), intent(out) :: rd !< rd() di/dr  (/4pir?)
   !Local variables
   real(kind=8), parameter :: fourpi=16.d0*atan(1.d0)
   integer :: i
   do i=1,nrad
      rr(i)=a_grd*exp(b_grd*(i-1))
      rw(i)=b_grd*rr(i)
      rd(i)=1.d0/rw(i)
      rw(i)=rw(i)*fourpi*rr(i)**2
      if (rr(i) > rmax) then
         nrad=i-1
         ! modify weights at en point for improved accuracy
         rw(1)=rw(1)*17.d0/48.d0
         rw(2)=rw(2)*59.d0/48.d0
         rw(3)=rw(3)*43.d0/48.d0
         rw(4)=rw(4)*49.d0/48.d0
         !finish
         return
      end if
   end do
   write(6,'(1x,a)') 'rmax too large, stopped in radgrid'
   stop
   
end subroutine radgrid
