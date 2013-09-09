!> @file
!! atomic program for generating and optimizing hgh pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/


!> Determine NPSP
subroutine detnp(nn,r,rad0,npsp)
   implicit none
   !Arguments
   integer, intent(in) :: nn
   real(kind=8), dimension(nn), intent(in) :: r
   real(kind=8), intent(in) :: rad0
   integer, intent(out) :: npsp
   !Local variables
   real(kind=8) :: rmin
   integer :: i
   rmin=1.d10
   do i=1,nn
      if (abs(r(i)-rad0) < rmin) then
      rmin=abs(r(i)-rad0)
      npsp=i
      end if
   end do
end subroutine detnp
