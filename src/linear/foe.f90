! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft(A,B,N,cc,ef,fscale,tmprtr)
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: A, B, ef, fscale, tmprtr
  integer,intent(in) :: n
  real(8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  real(kind=8) :: bma, bpa, y, arg, fac, tt, erfcc
  real(kind=8),dimension(1000) :: cf
  real(kind=8),parameter :: pi=4.d0*atan(1.d0)

  if (n>1000) stop 'chebft'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  do k=1,n
      y=cos(pi*(k-0.5D0)*(1.d0/n))
      arg=y*bma+bpa
      if (tmprtr.eq.0.d0) then
          cf(k)=.5d0*erfcc((arg-ef)*(1.d0/fscale))
      else
          cf(k)=1.d0/(1.d0+exp( (arg-ef)*(1.d0/tmprtr) ) )
      end if
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do  k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.D0/n)))
      end do
      cc(j)=fac*tt
  end do

end subroutine chebft



! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft2(a,b,n,cc)
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(kind=8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  real(kind=8),parameter :: pi=4.d0*atan(1.d0)
  real(kind=8) :: tt, y, arg, fac, bma, bpa
  real(kind=8),dimension(1000) :: cf

  if (n>1000) stop 'chebft2'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  ! 3 gives broder safety zone than 4
  !tt=3.0d0*n/(B-A)
  tt=4.d0*n/(b-a)
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      cf(k)=exp((arg-b)*tt)
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do
end subroutine chebft2





! Calculates the error function complement with an error of less than 1.2E-7
  function erfcc(x)
  implicit none

  ! Calling arguments
  real(8),intent(in) :: x
  real(8) :: erfcc

  ! Local variables
  real(8) :: z, t

  z=abs(x)
  t=1.d0/(1.+0.5d0*z)
  erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
        & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
        & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
  if (x.lt.0.) erfcc=2.D0-erfcc

end function erfcc
