!> Calculates gamma(x), fuer halb/ganzzahlige positive x
real(kind=8) function gamma(x)
   implicit none
   !Argument
   real(kind=8), intent(in) :: x
   !Local variables
   real(kind=8), parameter :: sqrtpi=sqrt(4.d0*atan(1.0d0))
   integer :: n,i
   real(kind=8) :: xx,gam

   ! print*,'entered gamma with x=',x
   if (x.le.0) then 
      write(6,*) 'stopped in gamma x<=0'
      stop
   endif
   if ( (x-int(x)) == 0.d0) then
      n=int(x)
      if (n == 0) then 
         write(6,*) 'stopped in gamma x=0'
         stop
      else
         gam=1
         do i=1,n-1
            gam=gam*i
         enddo
      endif
   else
      xx=x-0.5d0
      if ( (xx-int(xx)) /= 0.d0) then
         write(6,*) 'stopped in gamma x<>n+1/2'
         write(6,*) 'x=',x
         stop
      endif
      n=int(xx)
      gam=sqrtpi
      do i=1,2*n-1,2
         gam=gam*i
      enddo
      if (n > 0) gam=gam/2**n
   endif
   gamma = gam
end function gamma
