!!****f* BigDFT/razero
!! FUNCTION
!!   Set to zero an array x(n)
!!
!! SOURCE
!!
subroutine razero(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n)
  !Local variables
  integer :: i
!$omp parallel default(private) shared(n,x)
!$omp do 
  do i=1,n
     x(i)=0.d0
  end do
!$omp enddo
!$omp end parallel
END SUBROUTINE razero
!!***


!!****f* BigDFT/omp_razero
!! FUNCTION
!!   Set to zero an array x(n): omp version of razero
!!
!! SOURCE
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
!!***


!!****f* BigDFT/tenminustwenty
!! FUNCTION
!!   Set to 10^-20 an array x(n)
!!
!! SOURCE
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
!!***



!!****f* BigDFT/builtin_rand
!! FUNCTION
!!   Random Number generator from Numerical Recipes
!!   To be used for reproductability of the results
!!
!! SOURCE
!!
function builtin_rand(idum)
  implicit none
  integer, intent(inout) :: idum
  real(kind=4) :: builtin_rand
  !local variables
  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab
  real(kind=4), parameter :: am=1.e0/real(im,kind=4),eps=1.2e-7,rnmx=1.-eps
  integer :: j,k
  integer , save :: iy
  integer , dimension(NTAB), save :: iv
  data iv /NTAB*0/, iy /0/

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
end function builtin_rand
      
