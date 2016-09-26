! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module sort

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A, id)
  real(kind=8),dimension(:),intent(inout) :: A
  integer,dimension(:),intent(inout) :: id
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, id, iq)
     call QsortC(A(:iq-1), id(:iq-1))
     call QsortC(A(iq:), id(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, id, marker)
  real(kind=8),dimension(:),intent(inout) :: A
  integer,dimension(:),intent(inout) :: id
  integer,intent(out) :: marker
  integer :: i, j, itemp, ii
  real(kind=8) :: temp
  real(kind=8) :: x      ! pivot point

  ii = floor(0.5d0*real(size(A),kind=8))
  x = A(ii)
  !x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        itemp = id(i)
        id(i) = id(j)
        id(j) = itemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module sort

!!program sortdriver
!!  use qsort_c_module
!!  implicit none
!!  integer,parameter :: r = 10
!!  real(kind=8),dimension(1:r) :: myarray = &        ! (1:r)
!!     (/0, 50, 20, 25, 90, 10, 5, 1, 99, 75/)
!!  integer,dimension(1:r) :: id = (/1,2,3,4,5,6,7,8,9,10/)
!!  print *, "myarray is ", myarray
!!  call QsortC(myarray, id)
!!  print *, "sorted array is ", myarray
!!  print *, "sorted array is ", id
!!end program sortdriver
