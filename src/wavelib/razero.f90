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
