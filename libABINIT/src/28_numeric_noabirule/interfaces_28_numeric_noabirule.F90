!!****m* ABINIT/interfaces_28_numeric_noabirule
!! NAME
!! interfaces_28_numeric_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/28_numeric_noabirule
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

module interfaces_28_numeric_noabirule

 implicit none

interface
 subroutine sort_dp(n,list,iperm,tol)
  implicit none
  integer :: n
  double precision :: tol
  integer :: iperm(n)
  double precision :: list(n)
 end subroutine sort_dp
end interface

interface
 subroutine sort_int(n,list,iperm)
  implicit none
  integer :: n
  integer :: iperm(n)
  integer :: list(n)
 end subroutine sort_int
end interface

interface
 subroutine splfit(arg,derfun,fun,ider,newarg,newfun,numarg,numnew)
  implicit none
  integer, intent(in) :: ider
  integer, intent(in) :: numarg
  integer, intent(in) :: numnew
  double precision, intent(in) :: arg(numarg)
  double precision, intent(out) :: derfun(numnew)
  double precision, intent(in) :: fun(numarg,2)
  double precision, intent(in) :: newarg(numnew)
  double precision, intent(out) :: newfun(numnew)
 end subroutine splfit
end interface

interface
 subroutine spline( t, y, n, ybcbeg, ybcend, ypp )
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: ybcbeg
  double precision, intent(in) :: ybcend
  double precision, intent(in) :: t(n)
  double precision, intent(in) :: y(n)
  double precision, intent(out) :: ypp(n)
 end subroutine spline
end interface

interface
 subroutine spline_c( nomega_lo, nomega_li, omega_lo, omega_li, splined_li, tospline_lo)
  use abi_defs_basis
  implicit none
  integer, intent(in) :: nomega_li
  integer, intent(in) :: nomega_lo
  real(dp), intent(in) :: omega_li(nomega_li)
  real(dp), intent(in) :: omega_lo(nomega_lo)
  complex(dpc), intent(out) :: splined_li(nomega_li)
  complex(dpc), intent(in) :: tospline_lo(nomega_lo)
 end subroutine spline_c
end interface

interface
 subroutine spline_complex( t, y, n, ybcbeg, ybcend, ypp )
  use abi_defs_basis
  implicit none
  integer, intent(in) :: n
  complex(dpc), intent(in) :: ybcbeg
  complex(dpc), intent(in) :: ybcend
  real(dp), intent(in) :: t(n)
  complex(dpc), intent(in) :: y(n)
  complex(dpc), intent(out) :: ypp(n)
 end subroutine spline_complex
end interface

interface
 subroutine splint (nspline,xspline,yspline,ysplin2,&  
  &  nfit,xfit,yfit)
  implicit none
  integer, intent(in) :: nfit
  integer, intent(in) :: nspline
  double precision, intent(in) :: xfit(nfit)
  double precision, intent(in) :: xspline(nspline)
  double precision, intent(out) :: yfit(nfit)
  double precision, intent(in) :: ysplin2(nspline)
  double precision, intent(in) :: yspline(nspline)
 end subroutine splint
end interface

interface
 subroutine splint_complex (nspline,xspline,yspline,ysplin2,&  
  &  nfit,xfit,yfit)
  use abi_defs_basis
  implicit none
  integer, intent(in) :: nfit
  integer, intent(in) :: nspline
  real(dp), intent(in) :: xfit(nfit)
  real(dp), intent(in) :: xspline(nspline)
  complex(dpc), intent(out) :: yfit(nfit)
  complex(dpc), intent(in) :: ysplin2(nspline)
  complex(dpc), intent(in) :: yspline(nspline)
 end subroutine splint_complex
end interface

interface
 function uniformrandom(seed) 
  implicit none
  integer :: seed
  double precision :: uniformrandom
 end function uniformrandom
end interface

end module interfaces_28_numeric_noabirule
!!***
