!> @file
!! Routines to deal with the address of objects or external functions
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Module used by the module to manage the memory allocations
!! needed to pass to the C routines the correct address
!! in order to take the address of the metadata
module metadata_interfaces

  implicit none

  private

  integer, parameter :: longsize=20              !<could be lower
  character(len=*), parameter :: fmtlong='(i20)' !< conversion of long integer

  interface
     subroutine geti1(array,iadd)
       implicit none
       integer, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti1

     subroutine geti2(array,iadd)
       implicit none
       integer, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti2

     subroutine geti3(array,iadd)
       implicit none
       integer, dimension(:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti3

     subroutine geti4(array,iadd)
       implicit none
       integer, dimension(:,:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti4

     !template for character arrays of variable length
     include 'getadd-c-inc.f90'

     subroutine getl1(array,iadd)
       implicit none
       logical, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getl1

     subroutine getl2(array,iadd)
       implicit none
       logical, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getl2

     subroutine getdp1(array,iadd)
       implicit none
       double precision, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp1

     subroutine getdp2(array,iadd)
       implicit none
       double precision, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp2

     subroutine getdp3(array,iadd)
       implicit none
       double precision, dimension(:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp3

     subroutine getdp4(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp4

     subroutine getdp1ptr(array,iadd)
       implicit none
       double precision, dimension(:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp1ptr

     subroutine getdp2ptr(array,iadd)
       implicit none
       double precision, dimension(:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp2ptr

     subroutine getdp3ptr(array,iadd)
       implicit none
       double precision, dimension(:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp3ptr

     subroutine getdp4ptr(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp4ptr

     subroutine getdp5ptr(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp5ptr

     subroutine geti1ptr(array,iadd)
       implicit none
       integer, dimension(:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti1ptr

     subroutine geti2ptr(array,iadd)
       implicit none
       integer, dimension(:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti2ptr

  end interface

interface pad_array
  module procedure pad_i1,pad_i2,pad_i3,pad_i4
  module procedure pad_c1
  module procedure pad_l1,pad_l2
  module procedure pad_dp1,pad_dp2,pad_dp3,pad_dp4,pad_dp5
end interface

public :: pad_array,geti1,geti2,geti3,geti4
public :: getc1
!, getc1_2, getc1_3, getc1_4, getc1_5, getc1_6, getc1_7, getc1_8, getc1_9, &
!!$     getc1_10, getc1_11, getc1_12, getc1_13, getc1_14, getc1_15, getc1_16, getc1_17, getc1_18, &
!!$     getc1_19, getc1_20, getc1_21, getc1_22, getc1_23, getc1_24, getc1_25, getc1_26, getc1_27, &
!!$     getc1_28, getc1_29, getc1_30, getc1_31, getc1_32, getc1_33, getc1_34, getc1_35, getc1_36, &
!!$     getc1_37, getc1_38, getc1_39, getc1_40
public :: getl1,getl2
public :: getdp1,getdp2,getdp3,getdp4!,getlongaddress
public :: getdp1ptr,getdp2ptr,getdp3ptr,getdp4ptr,getdp5ptr,geti1ptr,geti2ptr
public :: address_toi,long_toa

contains

  subroutine pad_i1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    integer, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_i1

  subroutine pad_i2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    integer, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_i2

  subroutine pad_i3(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(3), intent(in) :: shp
    integer, dimension(shp(1),shp(2),shp(3)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,product(shp),product(shp(1:2))*(shp(3)+ndebug))

  end subroutine pad_i3

  subroutine pad_i4(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(4), intent(in) :: shp
    integer, dimension(shp(1),shp(2),shp(3),shp(4)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,product(shp),product(shp(1:3))*(shp(4)+ndebug))

  end subroutine pad_i4

  subroutine pad_c1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    character(len=*), dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_character(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_c1

  subroutine pad_l1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    logical, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_logical(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_l1

  subroutine pad_l2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    logical, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_logical(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_l2

  subroutine pad_dp1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    double precision, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_dp1

  subroutine pad_dp2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    double precision, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_dp2

  subroutine pad_dp3(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(3), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:2))*(shp(3)+ndebug))

  end subroutine pad_dp3

  subroutine pad_dp4(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(4), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3),shp(4)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:3))*(shp(4)+ndebug))

  end subroutine pad_dp4

  subroutine pad_dp5(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(5), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3),shp(4),shp(5)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:4))*(shp(5)+ndebug))

  end subroutine pad_dp5

  subroutine pad_double(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    double precision, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) call razero(ndim_tot,array)
    do i=ndim_tot+1,ndim_extra
       array(i)=d_nan()
    end do
  end subroutine pad_double

  subroutine pad_simple(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    real, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i,i_nan
    real :: r_nan1
    equivalence (r_nan1,i_nan)

    if (init) call razero_simple(ndim_tot,array)
    do i=ndim_tot+1,ndim_extra
       array(i)=r_nan()
    end do
  end subroutine pad_simple

  subroutine pad_logical(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    logical, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) then
       do i=1,ndim_tot
          array(i)=.false.
       end do
    end if
    do i=ndim_tot+1,ndim_extra
       array(i)=.true.
    end do
  end subroutine pad_logical

  subroutine pad_integer(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    integer, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) call razero_integer(ndim_tot,array)
    do i=ndim_tot+1,ndim_extra
       array(i)= 2147483647 !i_nan
    end do
  end subroutine pad_integer

  subroutine pad_character(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    character(len=*), dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) then
       do i=1,ndim_tot
          array(i)=repeat(' ',len(array(1)))
       end do
    end if
    do i=ndim_tot+1,ndim_extra
       array(i)=repeat('X',len(array(1)))
    end do
  end subroutine pad_character

  !> Function which specify NaN according to IEEE specifications
  function d_nan()
   implicit none
   double precision :: d_nan
   !local variables
   double precision :: dnan
   integer, dimension(2) :: inan
   equivalence (dnan, inan)
   ! This first assignment is for big-endian machines
   inan(1) = 2147483647
   ! The second assignment is for little-endian machines
   inan(2) = 2147483647
   d_nan = dnan
  end function d_nan

  !> Function which specify NaN according to IEEE specifications
  function r_nan()
   implicit none
   real :: r_nan
   !local variables
   real :: rnan
   integer :: inan
   equivalence (rnan, inan)
   inan = 2147483647
   r_nan = rnan
  end function r_nan

  function exa_toi(a)
    character(len=1), intent(in) :: a
    integer :: exa_toi
    select case(a)
    case('0')
       exa_toi=0
    case('1')
       exa_toi=1
    case('2')
       exa_toi=2
    case('3')
       exa_toi=3
    case('4')
       exa_toi=4
    case('5')
       exa_toi=5
    case('6')
       exa_toi=6
    case('7')
       exa_toi=7
    case('8')
       exa_toi=8
    case('9')
       exa_toi=9
    case('a')
       exa_toi=10
    case('b')
       exa_toi=11
    case('c')
       exa_toi=12
    case('d')
       exa_toi=13
    case('e')
       exa_toi=14
    case('f')
       exa_toi=15
    case default
       !raise an error
       print *,'a: ',a
       stop 'undefined value'
    end select
    
  end function exa_toi

  function address_toi(address)
    character(len=*), intent(in) ::  address
    integer(kind=8) :: address_toi
    !local variables
    integer :: i,l
    integer(kind=8) :: j
    character(len=1) :: a

    l=len_trim(address)
    address_toi=0
    do i=l-2,1,-1
       a=address(i+2:i+2)
       j=int(16**(l-2-i),kind=8)*int(exa_toi(a),kind=8)
       !print *,'i,a',i,a,exa_toi(a)
       address_toi=address_toi+j
    end do
    
  end function address_toi

  pure function long_toa(iadd)
    use yaml_strings
    implicit none 
    integer(kind=8), intent(in) :: iadd
    character(len=longsize) :: long_toa
    
    long_toa=adjustl(yaml_toa(iadd,fmt=fmtlong))

  end function long_toa

end module metadata_interfaces


!> Routine to call an external routine with an integer argument
subroutine call_external(routine,args)
  implicit none
  external :: routine                  !< Routine to be called
  integer(kind=8), intent(in) :: args  !< Argument of the called routine


  print *,'calling external, args',args
  
  if (args==0) then
     call routine()
  else
     call routine(args)
  end if
end subroutine call_external


!> Call the external routine with no argument
!! to be generalized to the case where extra arguments are needed
subroutine call_external_f(routine)!,args)
  implicit none
  external :: routine                  !< Routine to be called
!  integer(kind=8), intent(in) :: args


!  print *,'calling external, args',args
  
!  if (args==0) then

     call routine()
!  else
!     call routine(args)
!  end if

end subroutine call_external_f


!> Function which identify the address of the scalar object
!! associated to a unknown quantity
function f_loc(routine)
  implicit none
  external :: routine       !< Object
  integer(kind=8) :: f_loc  !< Address of the object routine

  call getlongaddress(routine,f_loc)

end function f_loc

