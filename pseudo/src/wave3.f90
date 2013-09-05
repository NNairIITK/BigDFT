!> @file
!! atomic program for generating and optimizing hgh pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/


subroutine wave3(ng,ll,xp,psi,expxpr,r,k,nint,w,dw,ddw)
   implicit none
   !Arguments
   integer, intent(in) :: ng,ll,k,nint
   real(kind=8), dimension(0:ng), intent(in) :: xp,psi
   real(kind=8), dimension(0:ng,nint), intent(in) :: expxpr
   real(kind=8), intent(in) :: r
   real(kind=8), intent(out) :: w,dw,ddw
   !Local variables
   real(kind=8) :: byr,byrr,rpll,uu1,uu2,uu3,uu4,uu5,uu6,tt1,tt2,tt3,tt4
   real(kind=8) :: xpown
   integer :: i

   w=0.d0
   dw=0.d0
   ddw=0.d0
   byr  = 1.d0/r
   byrr = byr/r
   rpll = xpown(r,ll)
   uu1  = -2.0d0*r*rpll
   uu2  = rpll*ll*byr
   uu3  = -2.0d0*rpll
   uu4  = +4.d0*rpll*r*r
   uu5  = -4.d0*rpll*ll
   uu6  = rpll*ll*(ll-1)*byrr
   do i=0,ng
      tt1=expxpr(i,k)
      tt2= psi(i)*tt1
      w=w + tt2
      tt3=tt2*xp(i)
      dw = dw +uu1*tt3 +uu2*tt2 
      tt4=tt3*xp(i)
      ddw=ddw +uu3*tt3+uu4*tt4+uu5*tt3+uu6*tt2
   enddo
   ! cwh
   if (ll.gt.0) w=w*rpll
   
end subroutine wave3



