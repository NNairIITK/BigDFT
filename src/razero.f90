!!****f* BigDFT/razero
!! NAME
!!   razero
!!
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
  do i=1,n
     x(i)=0.d0
  end do
end subroutine razero
!!***

!!****f* BigDFT/tenminustwenty
!! NAME
!!   tenminustwenty
!!
!! FUNCTION
!!   Set to 10^-20 an array x(n)
!!
!! SOURCE
!!
subroutine tenminustwenty(n,x,nproc)
  implicit none
! Arguments
  integer :: n,nproc
  real*8 :: x(n)
! Local variables
  integer :: i
  do i=1,n
     x(i)=1.d-20/real(nproc,kind=8)
  end do
END SUBROUTINE tenminustwenty
!!***
