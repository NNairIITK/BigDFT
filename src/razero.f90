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
