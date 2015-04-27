!> @file
!!  Define a module to extend numerical functions
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module numerics
  implicit none
  private

  interface safe_exp
     module procedure safe_dexp
  end interface safe_exp

  public :: safe_exp

  contains

    !> fpe-free way of calling exp.
    !! Crop the results to zero in the case of underflow
    pure function safe_dexp(x,extra_crop_order,underflow) result(ex)
      implicit none
      !> argument of the exponential function
      double precision, intent(in) :: x
      !> determine the value under which the result is assumed to be zero
      double precision, intent(in), optional :: underflow
      !> further restrict the valid range of the function
      !! by the order of magnitude indicated.
      !! Useful when the function has to be multiplied by extra terms
      !! The default is log of epsilon**2
      integer, intent(in), optional :: extra_crop_order
      double precision :: ex
      !local variables
      !> if the exponent is bigger than this value, the result is tiny(1.0)
      double precision, parameter :: mn_expo=-708.396418532264d0 ! = log(tiny(1.d0))
      !> if the exponent is lower than this value, the result is huge(1.0)
      double precision, parameter :: mx_expo=709.78271289338397d0 ! = log(huge(1.d0))
      !> the value of the cropping
      double precision, parameter :: crop_expo=72.0873067782343d0 ! = -2*log(epsilon(1.d0))
      double precision :: crop,mn,mx

      if (x==0.d0) then
         ex=1.d0
         return
      end if
      crop=crop_expo
      if (present(extra_crop_order)) crop=real(extra_crop_order,kind=8)
      mn=mn_expo+crop
      mx=mx_expo-crop
      if (present(underflow)) mn=log(abs(underflow))
      if (x > mn .and. x< mx) then
         ex=exp(x)
      else if (x <= mn) then
         ex=0.d0
      else
         ex=exp(mx)
      end if

    end function safe_dexp

    !> give a function which takes into account overflows and underflows even in the gaussian arguments
    pure function safe_gaussian(x0,x,alpha) result(gau)
      implicit none
      double precision, intent(in) :: x0 !< gaussian center
      double precision, intent(in) :: x !< argument
      !double precision, intent(in), optional :: sigma !<standard deviation
      double precision, intent(in) :: alpha !< exponent
      double precision :: gau
      !local variables
      !> if the sqrt is bigger than this value, the result is tiny(1.0)
      double precision, parameter :: mn_sqrt= sqrt(tiny(1.d0))
      !> if the sqrt is lower than this value, the result is huge(1.0)
      double precision, parameter :: mx_sqrt= sqrt(huge(1.d0))

      double precision :: gau_arg,xd

      !evaluate in safe way gau_arg
      xd=abs(x-x0) !assume that this is legal
      if (xd > mn_sqrt .and. xd< mx_sqrt) then
         xd=xd*xd
         gau_arg=-alpha*xd
         !if everything goes fine
         gau=safe_exp(gau_arg)
      else if (x <= mn_sqrt) then
         gau=1.d0
      else
         gau=0.d0
      end if
    end function safe_gaussian

end module numerics
