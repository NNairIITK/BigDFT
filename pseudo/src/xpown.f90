!> @file
!! Part of the pseudo program (pseudopotential generation)
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


!> Calculates r^n, n=int
!! This routine avoids problem with CRAY
real(kind=8) function xpown(r,n)
   implicit none
   !Arguments
   real(kind=8), intent(in) :: r
   integer, intent(in) :: n
   !Local variables
   integer :: i
   if (n == 0) then
      xpown=1.d0
      return
   endif
   if (r == 0.d0) then
      xpown=0.d0
      return
   endif
   xpown=1.d0
   do i=1,abs(n)
      xpown=xpown*r
   enddo
   if (n < 0) xpown=1.d0/xpown

end function xpown
