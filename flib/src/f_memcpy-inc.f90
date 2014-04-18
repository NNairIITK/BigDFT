!> @file
!! Include fortran file for memcpy interfaces
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine f_memcpy_i0(dest,src,n)
  implicit none
  integer, intent(in) :: n !<nelems
  integer, intent(inout) :: dest !<destination buffer address
  integer, intent(inout) :: src !<source buffer address
  !local variables
  integer :: n1,n2
  n1=n
  n2=n
  include 'f_memcpy-base-inc.f90'
end subroutine f_memcpy_i0

subroutine f_memcpy_i1(dest,src)
  implicit none
  integer, dimension(:), intent(inout) :: dest !<destination buffer
  integer, dimension(:), intent(in) :: src !<source buffer 
  !local variables
  integer :: n1,n2
  n1=size(dest)
  n2=size(src)
  include 'f_memcpy-base-inc.f90'
end subroutine f_memcpy_i1

subroutine f_memcpy_d0(dest,src,n)
  implicit none
  integer, intent(in) :: n !<nelems
  double precision, intent(inout) :: dest !<destination buffer address
  double precision, intent(inout) :: src !<source buffer address
  !local variables
  integer :: n1,n2
  n1=n
  n2=n
  include 'f_memcpy-base-inc.f90'
end subroutine f_memcpy_d0


subroutine f_memcpy_r0(dest,src,n)
  implicit none
  integer, intent(in) :: n !<nelems
  real, intent(inout) :: dest !<destination buffer address
  real, intent(inout) :: src !<source buffer address
  !local variables
  integer :: n1,n2
  n1=n
  n2=n
  include 'f_memcpy-base-inc.f90'
end subroutine f_memcpy_r0

subroutine f_memcpy_l0(dest,src,n)
  implicit none
  integer, intent(in) :: n !<nelems
  logical, intent(inout) :: dest !<destination buffer address
  logical, intent(inout) :: src !<source buffer address
  !local variables
  integer :: n1,n2
  n1=n
  n2=n
  include 'f_memcpy-base-inc.f90'
end subroutine f_memcpy_l0
