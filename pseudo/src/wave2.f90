!> @file
!! atomic program for generating and optimizing hgh pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


real(kind=8) function wave2(ng,ll,psi,expxpr,r,k,nint)
   !Arguments
   integer, intent(in) :: ng,ll,k,nint
   real(kind=8), dimension(0:ng), intent(in) :: psi
   real(kind=8), dimension(0:ng,nint), intent(in) :: expxpr
   real(kind=8), intent(in) :: r
   !Local variables
   real(kind=8) :: DDOT

   wave2=DDOT(ng+1,psi,1,expxpr(0,k),1)
   ! cwh
   if (ll.gt.0) wave2=wave2*r**ll

end function wave2
