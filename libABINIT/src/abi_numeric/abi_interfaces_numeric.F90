!!****m* ABINIT/abi_interfaces_numeric
!! NAME
!! abi_interfaces_numeric
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/abi_numeric
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

module abi_interfaces_numeric

 implicit none

interface
 subroutine abi_derfcf(derfc_yy,yy)
  use abi_defs_basis
  implicit none
  real(dp),intent(out) :: derfc_yy
  real(dp),intent(in) :: yy
 end subroutine abi_derfcf
end interface

interface
 subroutine abi_derf_ab(derf_yy,yy)
  use abi_defs_basis
  implicit none
  real(dp),intent(in) :: yy
  real(dp),intent(out) :: derf_yy
 end subroutine abi_derf_ab
end interface

interface
 subroutine abi_sort_dp(n,list,iperm,tol)
  implicit none
  integer :: n
  double precision :: tol
  integer :: iperm(n)
  double precision :: list(n)
 end subroutine abi_sort_dp
end interface

interface
 subroutine abi_sort_int(n,list,iperm)
  implicit none
  integer :: n
  integer :: iperm(n)
  integer :: list(n)
 end subroutine abi_sort_int
end interface

interface
 subroutine abi_mati3inv(mm,mit)
  implicit none
  integer,intent(out) :: mit(3,3)
  integer,intent(in) :: mm(3,3)
 end subroutine abi_mati3inv
end interface

interface
 subroutine abi_matr3inv(aa,ait)
  use abi_defs_basis
  implicit none
  real(dp),intent(in) :: aa(3,3)
  real(dp),intent(out) :: ait(3,3)
 end subroutine abi_matr3inv
end interface

interface
 subroutine abi_matrginv(a,lda,n)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  real(dp),intent(inout) :: a(lda,n)
 end subroutine abi_matrginv
end interface

interface
 function abi_uniformrandom(seed) 
  implicit none
  integer :: seed
  double precision :: abi_uniformrandom
 end function abi_uniformrandom
end interface

interface
 subroutine wrap2_pmhalf(num,red,shift)
  use abi_defs_basis
  implicit none
  real(dp),intent(in) :: num
  real(dp),intent(out) :: red
  real(dp),intent(out) :: shift
 end subroutine wrap2_pmhalf
end interface

end module abi_interfaces_numeric
!!***
