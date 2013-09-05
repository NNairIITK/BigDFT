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
