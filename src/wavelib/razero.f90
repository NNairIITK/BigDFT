!>   Set to zero an array x(n)
!!
subroutine razero(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n)
  !Local variables
  integer :: i
  do i=1,n
     x(i)=0.d0
  end do
END SUBROUTINE razero



!>   Set to zero an array x(n): omp version of razero
!!
subroutine omp_razero(n,x)
  use module_base
  implicit none
  !Arguments
  integer, intent(in) :: n
  real(wp), intent(out) :: x(n)
  !Local variables
  integer :: i

  !$omp do
  do i=1,n
     x(i)=0._wp
  end do
  !$omp enddo
END SUBROUTINE omp_razero



!>   Set to 10^-20 an array x(n)
!!
subroutine tenminustwenty(n,x,nproc)
  implicit none
! Arguments
  integer :: n,nproc
  real(kind=8) :: x(n)
! Local variables
  integer :: i
  do i=1,n
     x(i)=1.d-20/real(nproc,kind=8)
  end do
END SUBROUTINE tenminustwenty



!>   To be used in the following function.
!!
module randomData
  implicit none

  integer, parameter :: ntab=32

  logical :: start = .true.
  integer :: iy = 0
  integer, dimension(NTAB) :: iv
end module randomData



!>   Random Number generator from Numerical Recipes
!!   To be used for reproducibility of the results
!!
!!
function builtin_rand(idum)
  use randomData, only : ntab, iy, iv, start

  implicit none

  integer, intent(inout) :: idum
  real(kind=4) :: builtin_rand
  !local variables
  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ndiv=1+(im-1)/ntab
  real(kind=4), parameter :: am=1.e0/im,eps=1.2e-7,rnmx=1.-eps
  integer :: j,k

  if (start) then
     iv(:) = 0
     start = .false.
  end if
  if (idum <= 0.or. iy == 0) then
     idum=max(-idum,1)
     do j=ntab+8,1,-1
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum < 0) idum=idum+im
        if (j <= ntab) iv(j)=idum
     end do
     iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if (idum <= 0) idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  builtin_rand=min(am*iy,rnmx)
END FUNCTION builtin_rand

