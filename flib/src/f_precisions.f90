!> @file
!! Define precisions in portable format for future usage
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_precisions
  implicit none

  public

  !for reals and complex, to be verified if supported
  integer, parameter :: f_simple = selected_real_kind(6, 37)
  integer, parameter :: f_double = selected_real_kind(15, 307)
  integer, parameter :: f_quadruple = selected_real_kind(33, 4931)

  !for integers to be verified if they are supported, especially long
  integer, parameter :: f_short=selected_int_kind(4)
  integer, parameter :: f_integer=selected_int_kind(8)
  integer, parameter :: f_long=selected_int_kind(16)

  !logicals to be done also, and tested against bits and bytes with f_loc
  !integer, parameter :: bit=0 !not supported
  integer, parameter :: f_byte=1

  !> kind of a address, long integer on most machines. Might be detected at configure time
  integer, parameter :: f_address = f_long

  !> portable carriage return, contains both CR for unix and DOS
  character(len=*), parameter :: f_cr=char(13)//char(10)

  !>characters for the equality and the association of the parameters
  integer, parameter, private :: num_size=4
  character(len=num_size), parameter, private :: c_0='zero'
  character(len=num_size), parameter, private :: c_1='one '

  !>type for the definition of a parameter
  type, public :: f_parameter
     character(len=num_size) :: val
  end type f_parameter

  type(f_parameter), parameter :: f_0=f_parameter(c_0)
  type(f_parameter), parameter :: f_1=f_parameter(c_1)

  !bytes representation of true and false
  logical(f_byte), parameter :: f_T=.true._f_byte
  logical(f_byte), parameter :: f_F=.false._f_byte

  !> function to localize the address of anything
  integer(f_address), external :: f_loc

  interface assignment(=)
     module procedure set_d,set_i,set_li
  end interface assignment(=)

  private :: set_d,set_i,set_li
  
  contains

    elemental subroutine set_d(val,par)
      implicit none
      real(f_double), intent(out) :: val
      type(f_parameter), intent(in) :: par
      select case(par%val)
      case(c_0)
         val=0.0_f_double
      case(c_1)
         val=1.0_f_double
      end select
    end subroutine set_d

    elemental subroutine set_i(val,par)
      implicit none
      integer(f_integer), intent(out) :: val
      type(f_parameter), intent(in) :: par
      select case(par%val)
      case(c_0)
         val=int(0,f_integer)
      case(c_1)
         val=int(1,f_integer)
      end select
    end subroutine set_i

    elemental subroutine set_li(val,par)
      implicit none
      integer(f_long), intent(out) :: val
      type(f_parameter), intent(in) :: par
      select case(par%val)
      case(c_0)
         val=int(0,f_long)
      case(c_1)
         val=int(1,f_long)
      end select
    end subroutine set_li


end module f_precisions
