!> Module used in the builtin_rand function (see razero.f90)
module randomData
  implicit none

  integer, parameter :: ntab=32

  logical :: start = .true.
  integer :: iy = 0
  integer, dimension(NTAB) :: iv
end module randomData
