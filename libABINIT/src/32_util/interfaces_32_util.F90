!!****m* ABINIT/interfaces_32_util
!! NAME
!! interfaces_32_util
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/32_util
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!!
!! SOURCE

module interfaces_32_util

 implicit none

interface
 subroutine canon9(num,red,shift)
  use defs_basis
  implicit none
  real(dp),intent(in) :: num
  real(dp),intent(out) :: red
  real(dp),intent(out) :: shift
 end subroutine canon9
end interface

interface
 subroutine derfc(derfc_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derfc_yy
  real(dp),intent(in) :: yy
 end subroutine derfc
end interface

interface
 subroutine derf_ab(derf_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(in) :: yy
  real(dp),intent(out) :: derf_yy
 end subroutine derf_ab
end interface

interface
 subroutine mati3inv(mm,mit)
  implicit none
  integer,intent(out) :: mit(3,3)
  integer,intent(in) :: mm(3,3)
 end subroutine mati3inv
end interface

interface
 subroutine matr3inv(aa,ait)
  use defs_basis
  implicit none
  real(dp),intent(in) :: aa(3,3)
  real(dp),intent(out) :: ait(3,3)
 end subroutine matr3inv
end interface

interface
 subroutine matrginv(a,lda,n)
  use defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  real(dp),intent(inout) :: a(lda,n)
 end subroutine matrginv
end interface

interface
 subroutine wrap2_pmhalf(num,red,shift)
  use defs_basis
  implicit none
  real(dp),intent(in) :: num
  real(dp),intent(out) :: red
  real(dp),intent(out) :: shift
 end subroutine wrap2_pmhalf
end interface

end module interfaces_32_util
!!***
