!> @file
!! Manage dynamic memory allocation
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module used by the module to manage the memory allocations
module metadata_interfaces
  implicit none

  private

  integer, parameter :: longsize=20 !<could be lower
  character(len=*), parameter :: fmtlong='(i20)' !< conversion of long integer

  interface
     subroutine geti1(array,iadd)
       implicit none
       integer, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti1

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
  module procedure pad_i1,pad_i2,pad_dp1,pad_dp2,pad_dp3,pad_dp4,pad_dp5
end interface

public :: pad_array,geti1,getdp1,getdp2,getdp3,getdp4!,getlongaddress
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

  function long_toa(iadd)
    use yaml_strings
    implicit none 
    integer(kind=8), intent(in) :: iadd
    character(len=longsize) :: long_toa
    
    long_toa=adjustl(yaml_toa(iadd,fmt=fmtlong))

  end function long_toa

end module metadata_interfaces


!> Module used to manage memory allocations and de-allocations
module dynamic_memory
  !use m_profiling, except => ndebug, and=> d_nan, also=> r_nan
  use memory_profiling, except => ndebug
  use dictionaries, info_length => max_field_length
  implicit none

  private 

  integer, parameter :: namelen=32          !< length of the character variables
  integer, parameter :: error_string_len=80 !< length of error string
  integer, parameter :: ndebug=0            !< size of debug parameters
  !> errorcodes
  character(len=error_string_len) :: lasterror=repeat(' ',len(lasterror))
  integer, parameter :: SUCCESS                = 0
  integer, parameter :: INVALID_RANK           = -1979
  integer :: ierror=SUCCESS

  logical :: profile_initialized=.false.  !< global variables for initialization
  !>dictionaries needed for profiling storage
  type(dictionary), pointer :: dict_global
  type(dictionary), pointer :: dict_routine           
  type(dictionary), pointer :: dict_calling_sequence
  type(dictionary), pointer :: dict_codepoint=>null() !< save variable which says where we are in the code
  logical :: routine_opened=.false.                   !< global variable (can be stored in dictionaries)
  character(len=namelen) :: present_routine=repeat(' ',namelen)
!  character(len=namelen) :: last_opened_routine=repeat(' ',namelen)

  !> parameters for defitions of internal dictionary
  character(len=*), parameter :: arrayid='Array Id'
  character(len=*), parameter :: routineid='Allocating Routine Id'
  character(len=*), parameter :: sizeid='Size (Bytes)'
  character(len=*), parameter :: metadatadd='Address of metadata'
  character(len=*), parameter :: firstadd='Address of first element'
  character(len=*), parameter :: processid='Process Id'

  !error codes
  integer :: ERR_ALLOCATE
  integer :: ERR_DEALLOCATE
  integer :: ERR_MEMLIMIT

  !> Structure needed to allocate an allocatable array
  type, public :: malloc_information_all
     logical :: pin         !< flag to control the pinning of the address
     logical :: try         !< raise an exception
     logical :: put_to_zero !<initialize to zero after allocation
     integer :: rank        !< rank of the array
     integer, dimension(7) :: shape,lbounds,ubounds
     integer(kind=8) :: metadata_add      !<physical address of the fortran metadata
     character(len=namelen) :: array_id   !< label the array
     character(len=namelen) :: routine_id !<label the routine
  end type malloc_information_all

  !> Structure needed to allocate a pointer
  type, public :: malloc_information_ptr
     logical :: ptr         !< just to make the structures different, to see if needed
     logical :: pin         !< flag to control the pinning of the address
     logical :: try         !<raise an exception
     logical :: put_to_zero !<initialize to zero after allocation
     integer :: rank        !< rank of the pointer
     integer, dimension(7) :: shape,lbounds,ubounds
     integer(kind=8) :: metadata_add      !<physical address of the fortran metadata
     character(len=namelen) :: array_id   !< label the array
     character(len=namelen) :: routine_id !<label the routine
  end type malloc_information_ptr


  type, public :: array_bounds
     integer :: nlow  !<lower bounds
     integer :: nhigh !<higher bounds
  end type array_bounds

  interface assignment(=)
     module procedure i1_all,d1_all,d2_all,d3_all,d4_all
     module procedure d1_ptr,d2_ptr,d3_ptr,d4_ptr,d5_ptr
     module procedure i1_ptr,i2_ptr
  end interface

  interface operator(.to.)
     module procedure bounds
  end interface

  interface f_free
     module procedure i1_all_free,d1_all_free,d2_all_free,d1_all_free_multi,d3_all_free,d4_all_free
  end interface

  interface f_free_ptr
     module procedure d1_ptr_free,d2_ptr_free,d3_ptr_free,d4_ptr_free,d5_ptr_free
     module procedure i1_ptr_free,i2_ptr_free
  end interface


!!$  interface pad_with_nan
!!$     module procedure i_padding,dp_padding,c_padding,l_padding,sp_padding,dp_padding2,dp_padding3
!!$  end interface

  interface f_malloc
     module procedure f_malloc,f_malloc_simple,f_malloc_bounds,f_malloc_bound
  end interface

  interface f_malloc0
     module procedure f_malloc0,f_malloc0_simple,f_malloc0_bounds,f_malloc0_bound
  end interface

  interface f_malloc_ptr
     module procedure f_malloc_ptr,f_malloc_ptr_simple,f_malloc_ptr_bounds,f_malloc_ptr_bound
  end interface

  interface f_malloc0_ptr
     module procedure f_malloc0_ptr,f_malloc0_ptr_simple,f_malloc0_ptr_bounds,f_malloc0_ptr_bound
  end interface


  !> Public routines
  public :: f_malloc,f_malloc0,f_malloc_ptr,f_malloc0_ptr,f_malloc_dump_status
  public :: f_free,f_free_ptr
  public :: f_routine,f_release_routine,f_set_status,f_finalize
  public :: assignment(=),operator(.to.)

contains

  function bounds(nlow,nhigh)
    implicit none
    integer, intent(in) :: nlow,nhigh
    type(array_bounds) :: bounds

    bounds%nlow=nlow
    bounds%nhigh=nhigh
  end function bounds

  function malloc_information_ptr_null() result(m)
    implicit none
    type(malloc_information_ptr) :: m
    !local variables
    integer :: i

    m%ptr=.true.
    m%pin=.false.
    m%try=.false.
    m%put_to_zero=.false.
    m%metadata_add=0
    m%rank=1
    m%shape=0
    m%lbounds=1
    m%ubounds=0
    do i=1,namelen
       m%array_id(i:i)=' '
       m%routine_id(i:i)=' '
    end do
  end function malloc_information_ptr_null

  function malloc_information_all_null() result(m)
    implicit none
    type(malloc_information_all) :: m
    !local variables
    integer :: i

    m%pin=.false.
    m%try=.false.
    m%put_to_zero=.false.
    m%metadata_add=0
    m%rank=1
    m%shape=0
    m%lbounds=1
    m%ubounds=0
    do i=1,namelen
       m%array_id(i:i)=' '
       m%routine_id(i:i)=' '
    end do

  end function malloc_information_all_null

  !for rank-1 arrays
  function f_malloc_simple(size,id,routine_id,try) result(m)
    integer, intent(in) :: size
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt

    m=malloc_information_all_null()
    m%rank=1
    m%shape(1)=size
    m%ubounds(1)=m%shape(1)

    include 'f_malloc-inc.f90'

  end function f_malloc_simple

  !for rank-1 arrays
  function f_malloc_bound(bounds,id,routine_id,try) result(m)
    type(array_bounds), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt

    m=malloc_information_all_null()
    m%rank=1
    m%lbounds(1)=bounds%nlow
    m%ubounds(1)=bounds%nhigh
    m%shape(1)=m%ubounds(1)-m%lbounds(1)+1

    include 'f_malloc-inc.f90'

  end function f_malloc_bound

  !define the allocation information for  arrays of different rank
  function f_malloc_bounds(bounds,id,routine_id,try) result(m)
    implicit none
    type(array_bounds), dimension(:), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_all_null()
    
    m%rank=size(bounds)
    do i=1,m%rank
       m%lbounds(i)=bounds(i)%nlow
       m%ubounds(i)=bounds(i)%nhigh
       m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
       end do

    include 'f_malloc-inc.f90'

  end function f_malloc_bounds


  !define the allocation information for  arrays of different rank
  function f_malloc(shape,id,routine_id,lbounds,ubounds,try) result(m)
    implicit none
    integer, dimension(:), intent(in), optional :: shape,lbounds,ubounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_all_null()

    include 'f_malloc-extra-inc.f90'
    include 'f_malloc-inc.f90'

  end function f_malloc


  !>for rank-1 arrays
  function f_malloc0_simple(size,id,routine_id,try) result(m)
    implicit none
    integer, intent(in) :: size
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt

    m=malloc_information_all_null()
    m%put_to_zero=.true.
    m%rank=1
    m%shape(1)=size
    m%ubounds(1)=m%shape(1)

    include 'f_malloc-inc.f90'

  end function f_malloc0_simple

  !>for rank-1 arrays
  function f_malloc0_bound(bounds,id,routine_id,try) result(m)
    implicit none
    type(array_bounds), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt

    m=malloc_information_all_null()
    m%put_to_zero=.true.
    m%rank=1
    m%lbounds(1)=bounds%nlow
    m%ubounds(1)=bounds%nhigh
    m%shape(1)=m%ubounds(1)-m%lbounds(1)+1

    include 'f_malloc-inc.f90'
  end function f_malloc0_bound

  function f_malloc0_bounds(bounds,id,routine_id,try) result(m)
    implicit none
    type(array_bounds), dimension(:), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_all_null()

    m%put_to_zero=.true.

    m%rank=size(bounds)
    do i=1,m%rank
       m%lbounds(i)=bounds(i)%nlow
       m%ubounds(i)=bounds(i)%nhigh
       m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
    end do

    include 'f_malloc-inc.f90'
  end function f_malloc0_bounds

  !> define the allocation information for  arrays of different rank
  function f_malloc0(shape,id,routine_id,lbounds,ubounds,try) result(m)
    implicit none
    integer, dimension(:), intent(in), optional :: shape,lbounds,ubounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_all) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_all_null()

    m%put_to_zero=.true.

    include 'f_malloc-extra-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc0

  !> for rank-1 arrays
  function f_malloc_ptr_simple(size,id,routine_id,try) result(m)
    integer, intent(in) :: size
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt

    m=malloc_information_ptr_null()
    m%rank=1
    m%shape(1)=size
    m%ubounds(1)=m%shape(1)

    include 'f_malloc-inc.f90'
  end function f_malloc_ptr_simple

  !for rank-1 arrays
  function f_malloc_ptr_bound(bounds,id,routine_id,try) result(m)
    type(array_bounds), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt

    m=malloc_information_ptr_null()
    m%rank=1
    m%lbounds(1)=bounds%nlow
    m%ubounds(1)=bounds%nhigh
    m%shape(1)=m%ubounds(1)-m%lbounds(1)+1

    include 'f_malloc-inc.f90'
  end function f_malloc_ptr_bound

    !define the allocation information for  arrays of different rank
  function f_malloc_ptr_bounds(bounds,id,routine_id,try) result(m)
    implicit none
    type(array_bounds), dimension(:), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_ptr_null()
    
    m%rank=size(bounds)
    do i=1,m%rank
       m%lbounds(i)=bounds(i)%nlow
       m%ubounds(i)=bounds(i)%nhigh
       m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
    end do

    include 'f_malloc-inc.f90'
  end function f_malloc_ptr_bounds


  !define the allocation information for  arrays of different rank
  function f_malloc_ptr(shape,id,routine_id,lbounds,ubounds,try) result(m)
    implicit none
    integer, dimension(:), intent(in), optional :: shape,lbounds,ubounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_ptr_null()

    include 'f_malloc-extra-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_ptr

  !for rank-1 arrays
  function f_malloc0_ptr_simple(size,id,routine_id,try) result(m)
    integer, intent(in) :: size
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt

    m=malloc_information_ptr_null()
    m%put_to_zero=.true.
    m%rank=1
    m%shape(1)=size
    m%ubounds(1)=m%shape(1)

    include 'f_malloc-inc.f90'
  end function f_malloc0_ptr_simple

  !for rank-1 arrays
  function f_malloc0_ptr_bound(bounds,id,routine_id,try) result(m)
    type(array_bounds), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt

    m=malloc_information_ptr_null()
    m%put_to_zero=.true.
    m%rank=1
    m%lbounds(1)=bounds%nlow
    m%ubounds(1)=bounds%nhigh
    m%shape(1)=m%ubounds(1)-m%lbounds(1)+1

    include 'f_malloc-inc.f90'
  end function f_malloc0_ptr_bound

  !define the allocation information for  arrays of different rank
  function f_malloc0_ptr_bounds(bounds,id,routine_id,try) result(m)
    implicit none
    type(array_bounds), dimension(:), intent(in) :: bounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_ptr_null()

    m%put_to_zero=.true.
    m%rank=size(bounds)
    do i=1,m%rank
       m%lbounds(i)=bounds(i)%nlow
       m%ubounds(i)=bounds(i)%nhigh
       m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
    end do

    include 'f_malloc-inc.f90'
  end function f_malloc0_ptr_bounds

  !define the allocation information for  arrays of different rank
  function f_malloc0_ptr(shape,id,routine_id,lbounds,ubounds,try) result(m)
    implicit none
    integer, dimension(:), intent(in), optional :: shape,lbounds,ubounds
    logical, intent(in), optional :: try
    character(len=*), intent(in), optional :: id,routine_id
    type(malloc_information_ptr) :: m
    !local variables
    integer :: lgt,i
    m=malloc_information_ptr_null()
    m%put_to_zero=.true.

    include 'f_malloc-extra-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc0_ptr

  !> This routine adds the corresponding subprogram name to the dictionary
  !! and prepend the dictionary to the global info dictionary
  !! if it is called more than once for the same name it has no effect
  subroutine f_routine(id)
    implicit none
    character(len=*), intent(in), optional :: id
    !local variables
    integer :: lgt

    if (.not. present(id)) return !no effect

    if (trim(present_routine) /= trim(id)) then
       if(associated(dict_routine)) then
          call prepend(dict_global,dict_routine)
          nullify(dict_routine)
       end if
       !this means that the previous routine has not been closed
       if (routine_opened) then
          call open_routine(dict_codepoint)
!          last_opened_routine=present_routine
       end if
       routine_opened=.true.
       call add(dict_codepoint,trim(id))

       present_routine=repeat(' ',namelen)
       lgt=min(len(id),namelen)
       present_routine(1:lgt)=id(1:lgt)

    end if
  end subroutine f_routine

  !> Close a previously opened routine
  subroutine f_release_routine()
!    use yaml_output
    implicit none
    if (associated(dict_routine)) then
       call prepend(dict_global,dict_routine)
       nullify(dict_routine)
    end if
    call close_routine(dict_codepoint,.not. routine_opened)!trim(dict_key(dict_codepoint)))
!!$    call yaml_open_map('Codepoint after closing')
!!$    call yaml_map('Potential Reference Routine',trim(dict_key(dict_codepoint)))
!!$      call yaml_dict_dump(dict_codepoint)
!!$    call yaml_close_map()
    present_routine=trim(dict_key(dict_codepoint))
    !last_opened_routine=trim(dict_key(dict_codepoint))!repeat(' ',namelen)
    routine_opened=.false.
  end subroutine f_release_routine

  subroutine open_routine(dict)
    implicit none
    type(dictionary), pointer :: dict
    !local variables
    integer :: ival
    character(len=info_length) :: routinename
    type(dictionary), pointer :: dict_tmp

    !now imagine that a new routine is created
    ival=dict_len(dict)-1
    routinename=dict//ival

    !call yaml_map('The routine which has to be converted is',trim(routinename))

    call pop(dict,ival)

    dict_tmp=>dict//ival//trim(routinename)
    dict => dict_tmp
    nullify(dict_tmp)

  end subroutine open_routine

  subroutine close_routine(dict,jump_up)
!    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(in) :: jump_up
    !character(len=*), intent(in) :: name
    !local variables

    !integer :: ival
    type(dictionary), pointer :: dict_tmp

    if (.not. associated(dict)) stop 'ERROR, routine not associated' 

!!$    !jump_up=(trim(present_routine) /= trim(name))
!!$    jump_up=(trim(last_opened_routine) /= trim(name))
!!$
!!$    call yaml_open_map('Test of the example')
!!$    call yaml_map('Name',trim(name))
!!$    call yaml_map('Last opened routine',trim(last_opened_routine))
!!$    call yaml_map('Present Routine',present_routine)
!!$    call yaml_map('Willing to jump up',jump_up)
!!$    call yaml_close_map()

    if (jump_up) then
       !now the routine has to be closed
       !we should jump at the upper level
       dict_tmp=>dict%parent 
       if (associated(dict_tmp%parent)) then
          nullify(dict)
          !this might be null if we are at the topmost level
          dict=>dict_tmp%parent
       end if
       nullify(dict_tmp)
    end if

  end subroutine close_routine

  !> Decide the error messages associated to the dynamic memory
  subroutine malloc_errors()
    use error_handling
    implicit none
    
    call f_err_define(err_name='ERR_ALLOCATE',err_msg='Allocation error',err_id=ERR_ALLOCATE,&
         err_action='Control the order of the allocation of if the memory limit has been reached')
    call f_err_define(err_name='ERR_DEALLOCATE',err_msg='Dellocation error',err_id=ERR_DEALLOCATE,&
         err_action='Control the order of the allocation of if the memory limit has been reached')
    call f_err_define(err_name='ERR_MEMLIMIT',err_msg='Memory limit reached',err_id=ERR_MEMLIMIT,&
         err_action='Control the size of the arrays needed for this run with bigdft-tool program')

    
  end subroutine malloc_errors


  !> Initialize the library
  subroutine f_set_status(memory_limit,output_level,logfile_name,unit,iproc)
    use yaml_output, only: yaml_date_and_time_toa
    implicit none
    !Arguments
    character(len=*), intent(in), optional :: logfile_name   !< Name of the logfile
    real(kind=4), intent(in), optional :: memory_limit       !< Memory limit
    integer, intent(in), optional :: output_level            !< Level of output for memocc
                                                             !! 0 no file, 1 light, 2 full
    integer, intent(in), optional :: unit                    !< Indicate file unit for the output
    integer, intent(in), optional :: iproc                   !< Process Id (used to dump, by default one 0)

    if (.not. profile_initialized) then
       profile_initialized=.true.
       call malloc_errors()
       !initalize the dictionary with the allocation information
       nullify(dict_routine)
       call dict_init(dict_global)
       call set(dict_global//'Timestamp of Profile initialization',trim(yaml_date_and_time_toa()))
       !Process Id (used to dump)
       call set(dict_global//processid,0)
       call dict_init(dict_calling_sequence)
       !in principle the calling sequence starts from the main
       dict_codepoint => dict_calling_sequence//'Calling sequence of Main program'
    end if

    if (present(memory_limit)) call memocc_set_memory_limit(memory_limit)

    if (present(output_level)) call memocc_set_state(output_level)

    if (present(unit)) call memocc_set_stdout(unit)

    if (present(logfile_name)) call memocc_set_filename(logfile_name)
       
    if (present(iproc)) call set(dict_global//processid,iproc)
  end subroutine f_set_status

  !> Finalize f_malloc (Display status)
  subroutine f_finalize(dump)
    use yaml_output, only: yaml_warning,yaml_open_map,yaml_close_map,yaml_dict_dump,yaml_get_default_stream
    implicit none
    !Arguments
    logical, intent(in), optional :: dump !< Dump always information, 
                                          !! otherwise only for Process Id == 0 and errors
    !local variables
    integer :: pid
    logical :: dump_status
    !integer :: unt
    !put the last values in the dictionary if not freed
    if (associated(dict_routine)) then
       !call yaml_get_default_stream(unt)
       !call yaml_stream_attributes(unit=unt)
       !call yaml_warning('Not all the arrays have been freed: memory leaks are possible')
       call prepend(dict_global,dict_routine)
       nullify(dict_routine)
       !      end if
       !      if (.false.) then !residual memory to be defined
    end if

    if (present(dump)) then
       dump_status=.true.
    else 
       pid = dict_global//processid
       if (pid == 0) then
          dump_status=.true.
       else
          dump_status=.false.
       end if
       !Print if error
       if (dict_size(dict_global) == 2) dump_status=.false.
       !print *,'dict_size',dict_size(dict_global)
       !print *,'dict_len',dict_len(dict_global)
    end if
    if (dump_status) then
       call yaml_open_map('Status of the memory at finalization')
       !call yaml_dict_dump(dict_global)
       call dump_leaked_memory(dict_global)
       call yaml_close_map()
    end if
       call dict_free(dict_global)
    !    call yaml_open_map('Calling sequence')
    !    call yaml_dict_dump(dict_calling_sequence)
    !    call yaml_close_map()
       call dict_free(dict_calling_sequence)

       call memocc_report()
!    end if
    profile_initialized=.false.
    present_routine=repeat(' ',namelen)
    routine_opened=.false.
  end subroutine f_finalize

  !> Check error of allocations or deallocations
  subroutine check_for_errors(ierror,try)
    use yaml_output!, only: yaml_warning,yaml_open_map,yaml_map,yaml_close_map,yaml_dict_dump,yaml_get_default_stream
    implicit none
    logical, intent(in) :: try
    integer, intent(in) :: ierror
    !local variables
    integer :: unt

    !recuperate possible error
    if (ierror /= INVALID_RANK) lasterror='Fortran (de)allocation problem'

    !raise exception if exit
    if (ierror /= SUCCESS) then
       if (try) then
          return
       else
          call yaml_open_map('(de)allocation error exiting')
             call yaml_map('Error code',ierror)
             call yaml_map('Last error',lasterror)
          call yaml_close_map()
          !write(*,*)'(de)allocation error, exiting. Error code:',ierror
          !write(*,*)'last error:',lasterror
          call yaml_get_default_stream(unt)
          if (associated(dict_routine)) then
             call yaml_open_map('Status of the routine before exiting')
             !call yaml_dict_dump(dict_routine)
             call dump_leaked_memory(dict_routine)
             call yaml_close_map()
          end if
          call yaml_dict_dump(dict_calling_sequence)
          stop
       end if
    end if

  end subroutine check_for_errors


  !> Dump all allocations
  subroutine dump_leaked_memory(dict)
    use metadata_interfaces, only: address_toi
     use yaml_output
     implicit none
     type(dictionary), pointer, intent(in) :: dict
     !Local variables
     type(dictionary), pointer :: dict_ptr!, dict_tmp
     character(len=256) :: array_id
     dict_ptr => dict_next(dict)
     do while(associated(dict_ptr))
!!$        call yaml_map('Key has been found',has_key(dict_ptr,trim(arrayid)))
!!$        call try()
!!$        array_id = dict_ptr//arrayid
!!$        call close_try()
        if (has_key(dict_ptr,trim(arrayid))) then
           array_id = dict_ptr//arrayid
           call yaml_open_map(trim(array_id))
           call yaml_dict_dump(dict_ptr)
           call yaml_map(metadatadd,trim(dict_key(dict_ptr)))
!!$           call yaml_map('Long integer conversion',&
!!$                address_toi(trim(adjustl(dict_key(dict_ptr)))))
           call yaml_close_map()
        else
           call yaml_open_map(trim(dict_key(dict_ptr)))
           call yaml_dict_dump(dict_ptr)
           call yaml_close_map()
        end if
!!$           PRINT *,'ERROR',TRY_ERROR()
!!$           call yaml_map('key',dict_key(dict_ptr))
!!$           print *,'nelems: ',dict_ptr%data%nelems
!!$           print *,'key: "',trim(dict_ptr%data%key),'"'
!!$           print *,'value: "',trim(dict_ptr%data%value),'"'
!!$           call yaml_dict_dump(dict_ptr%parent)
!!$           print *,'aa'
        dict_ptr=>dict_next(dict_ptr)
     end do
  end subroutine dump_leaked_memory


  !> Use the address of the allocated pointer to profile the deallocation
  subroutine profile_deallocation(ierr,ilsize,address)
    use yaml_output!, only: yaml_warning,yaml_open_map,yaml_close_map,yaml_dict_dump,yaml_map
    implicit none
    integer, intent(in) :: ierr
    integer(kind=8), intent(in) :: ilsize
    character(len=*), intent(in) :: address
    !local variables
    logical :: use_global
    character(len=namelen) :: array_id,routine_id
    integer(kind=8) :: jlsize
    type(dictionary), pointer :: dict_add!,dict_tmp

    call check_for_errors(ierr,.false.)
    !search in the dictionaries the address
    dict_add=>find_key(dict_routine,trim(address))
    if (.not. associated(dict_add)) then
       dict_add=>find_key(dict_global,trim(address))
       if (.not. associated(dict_add)) then
          print *,'address:',trim(address)
          call f_malloc_dump_status()
          call yaml_dict_dump(dict_calling_sequence)
          stop 'profile deallocations: address not present'
       else
          use_global=.true.
       end if
    else
       use_global=.false.
    end if
    !the global dictionary should be used instead
    array_id=dict_add//arrayid
    routine_id=dict_add//routineid
    jlsize=dict_add//sizeid

    call memocc(ierr,-int(ilsize),trim(array_id),trim(routine_id))
!!$    call yaml_map('Deallocating',trim(array_id))
!!$    call yaml_map('Use global',use_global)
!!$    call yaml_map('Associated global',associated(dict_global))
    if (use_global) then
       !call yaml_dict_dump(dict_global)
       call pop(dict_global,trim(address))
    else
       call pop(dict_routine,trim(address))
    end if

!!$    call yaml_comment('Test',hfill='-')
!!$!    call yaml_map('Associated',associated(dict_tmp))
!!$    call yaml_map('Associated routine',associated(dict_routine))
!!$    !call yaml_dict_dump(dict_tmp)

  end subroutine profile_deallocation


  subroutine profile_allocation(ierr,iadd,address,sizeof,m)
    implicit none
    integer, intent(in) :: ierr,sizeof
    integer(kind=8), intent(in) :: iadd
    character(len=*), intent(in) :: address
    type(malloc_information_all), intent(in) :: m
    !local variables
    integer :: i
    integer(kind=8) :: ilsize
    type(dictionary), pointer :: dict_tmp

    if (.not. associated(dict_routine)) then
       call dict_init(dict_routine)
    end if
    !size
    ilsize=int(sizeof,kind=8)
    do i=1,m%rank
       ilsize=ilsize*int(m%shape(i),kind=8)
    end do
    !create the dictionary array
    !add the array to the routine
    call dict_array(m%routine_id,m%array_id,ilsize,dict_tmp)
    call set(dict_routine//trim(address),dict_tmp)
    call check_for_errors(ierr,m%try)
    call memocc(ierr,product(m%shape(1:m%rank))*sizeof,m%array_id,m%routine_id)
  contains

    subroutine dict_array(routine_id,array_id,size,dict_tmp)
      implicit none
      character(len=*), intent(in) :: array_id,routine_id
      integer(kind=8), intent(in) :: size !< in bytes
      type(dictionary), pointer :: dict_tmp
      nullify(dict_tmp)
      call dict_init(dict_tmp)
      call set(dict_tmp//arrayid,trim(array_id))
      call set(dict_tmp//sizeid,size)
      call set(dict_tmp//routineid,trim(routine_id))
      call set(dict_tmp//metadatadd,iadd)

    end subroutine dict_array

  end subroutine profile_allocation

  subroutine f_malloc_dump_status()
    use yaml_output
    implicit none
!!$    if (associated(dict_routine)) call yaml_dict_dump(dict_routine)
!!$    call yaml_dict_dump(dict_global)
    call yaml_newline()
!    call yaml_map('Present routine',trim(present_routine))
!    call yaml_open_map(Routine dictionary')
    call yaml_dict_dump(dict_calling_sequence)
    if (associated(dict_routine)) then
       call yaml_open_map('Routine dictionary')
       call dump_leaked_memory(dict_routine)
       call yaml_close_map()
    end if
    call yaml_open_map('Global dictionary')
    !call yaml_dict_dump(dict_global)
    call dump_leaked_memory(dict_global)
    call yaml_close_map()

  end subroutine f_malloc_dump_status

  subroutine i1_all(array,m)
    use metadata_interfaces, metadata_address => geti1
    implicit none
    type(malloc_information_all), intent(in) :: m
    integer, dimension(:), allocatable, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd
    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine i1_all

  subroutine i1_all_free(array)
    use metadata_interfaces, metadata_address => geti1
    implicit none
    integer, dimension(:), allocatable, intent(inout) :: array
    include 'deallocate-inc.f90' 
  end subroutine i1_all_free


  subroutine d1_all(array,m)
    use metadata_interfaces, metadata_address => getdp1
    implicit none
    type(malloc_information_all), intent(in) :: m
    double precision, dimension(:), allocatable, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd
    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine d1_all

  subroutine d1_all_free(array)
    use metadata_interfaces, metadata_address => getdp1
    implicit none
    double precision, dimension(:), allocatable, intent(inout) :: array
    include 'deallocate-inc.f90' 
  end subroutine d1_all_free


  subroutine d2_all(array,m)
    use metadata_interfaces, metadata_address => getdp2
    implicit none
    type(malloc_information_all), intent(in) :: m
    double precision, dimension(:,:), allocatable, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd

    allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine d2_all

  subroutine d2_all_free(array)
    use metadata_interfaces, metadata_address => getdp2
    implicit none
    double precision, dimension(:,:), allocatable, intent(inout) :: array
    include 'deallocate-inc.f90' 
  end subroutine d2_all_free


  subroutine d3_all(array,m)
    use metadata_interfaces, metadata_address => getdp3
    implicit none
    type(malloc_information_all), intent(in) :: m
    double precision, dimension(:,:,:), allocatable, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd

    allocate(array(m%lbounds(1):m%ubounds(1),&
         m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine d3_all

  subroutine d3_all_free(array)
    use metadata_interfaces, metadata_address => getdp3
    implicit none
    double precision, dimension(:,:,:), allocatable, intent(inout) :: array
    include 'deallocate-inc.f90' 
  end subroutine d3_all_free


  subroutine d4_all(array,m)
    use metadata_interfaces, metadata_address => getdp4
    implicit none
    type(malloc_information_all), intent(in) :: m
    double precision, dimension(:,:,:,:), allocatable, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd

    allocate(array(m%lbounds(1):m%ubounds(1),&
         m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
         m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine d4_all

  subroutine d4_all_free(array)
    use metadata_interfaces, metadata_address => getdp4
    implicit none
    double precision, dimension(:,:,:,:), allocatable, intent(inout) :: array
    include 'deallocate-inc.f90' 
  end subroutine d4_all_free

  subroutine d1_all_free_multi(arrayA,arrayB,arrayC,arrayD,arrayE,arrayF,arrayG,arrayH)
    implicit none
    double precision, dimension(:), allocatable, intent(inout) :: arrayA
    double precision, dimension(:), allocatable, intent(inout) :: arrayB
    double precision, dimension(:), allocatable, optional, intent(inout) :: arrayC
    double precision, dimension(:), allocatable, optional, intent(inout) :: arrayD
    double precision, dimension(:), allocatable, optional, intent(inout) :: arrayE
    double precision, dimension(:), allocatable, optional, intent(inout) :: arrayF
    double precision, dimension(:), allocatable, optional, intent(inout) :: arrayG
    double precision, dimension(:), allocatable, optional, intent(inout) :: arrayH

    call d1_all_free(arrayA)
    call d1_all_free(arrayB)
    if (present(arrayC)) then
       call d1_all_free(arrayC)
    end if
    if (present(arrayD)) then
       call d1_all_free(arrayD)
    end if
    if (present(arrayE)) then
       call d1_all_free(arrayE)
    end if
    if (present(arrayF)) then
       call d1_all_free(arrayF)
    end if
    if (present(arrayG)) then
       call d1_all_free(arrayG)
    end if
    if (present(arrayH)) then
       call d1_all_free(arrayH)
    end if
  end subroutine d1_all_free_multi

  !pointers
  subroutine d1_ptr(array,m)
    use metadata_interfaces, metadata_address => getdp1ptr
    implicit none
    type(malloc_information_ptr), intent(in) :: m
    double precision, dimension(:), pointer, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd

    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)

    include 'allocate-inc.f90'
  end subroutine d1_ptr

  subroutine d1_ptr_free(array)
    use metadata_interfaces, metadata_address => getdp1ptr
    implicit none
    double precision, dimension(:), pointer, intent(inout) :: array
    include 'deallocate-inc.f90'
    nullify(array)
  end subroutine d1_ptr_free

  subroutine i1_ptr(array,m)
    use metadata_interfaces, metadata_address => geti1ptr
    implicit none
    type(malloc_information_ptr), intent(in) :: m
    integer, dimension(:), pointer, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd

    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)

    include 'allocate-inc.f90'
  end subroutine i1_ptr

  subroutine i1_ptr_free(array)
    use metadata_interfaces, metadata_address => geti1ptr
    implicit none
    integer, dimension(:), pointer, intent(inout) :: array
    include 'deallocate-inc.f90'
    nullify(array)
  end subroutine i1_ptr_free


  subroutine d2_ptr(array,m)
    use metadata_interfaces, metadata_address => getdp2ptr
    implicit none
    type(malloc_information_ptr), intent(in) :: m
    double precision, dimension(:,:), pointer, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd
    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)

    include 'allocate-inc.f90'
  end subroutine d2_ptr

  subroutine d2_ptr_free(array)
    use metadata_interfaces, metadata_address => getdp2ptr
    implicit none
    double precision, dimension(:,:), pointer, intent(inout) :: array
    include 'deallocate-inc.f90'
    nullify(array)
  end subroutine d2_ptr_free

  subroutine i2_ptr(array,m)
    use metadata_interfaces, metadata_address => geti2ptr
    implicit none
    type(malloc_information_ptr), intent(in) :: m
    integer, dimension(:,:), pointer, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd

    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)

    include 'allocate-inc.f90'
  end subroutine i2_ptr

  subroutine i2_ptr_free(array)
    use metadata_interfaces, metadata_address => geti2ptr
    implicit none
    integer, dimension(:,:), pointer, intent(inout) :: array
    include 'deallocate-inc.f90'
    nullify(array)
  end subroutine i2_ptr_free


  subroutine d3_ptr(array,m)
    use metadata_interfaces, metadata_address => getdp3ptr
    implicit none
    type(malloc_information_ptr), intent(in) :: m
    double precision, dimension(:,:,:), pointer, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd
    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
         m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine d3_ptr

  subroutine d3_ptr_free(array)
    use metadata_interfaces, metadata_address => getdp3ptr
    implicit none
    double precision, dimension(:,:,:), pointer, intent(inout) :: array
    include 'deallocate-inc.f90'
    nullify(array)
  end subroutine d3_ptr_free

  subroutine d4_ptr(array,m)
    use metadata_interfaces, metadata_address => getdp4ptr
    implicit none
    type(malloc_information_ptr), intent(in) :: m
    double precision, dimension(:,:,:,:), pointer, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd
    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
         m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine d4_ptr

  subroutine d4_ptr_free(array)
    use metadata_interfaces, metadata_address => getdp4ptr
    implicit none
    double precision, dimension(:,:,:,:), pointer, intent(inout) :: array
    include 'deallocate-inc.f90'
    nullify(array)
  end subroutine d4_ptr_free

  subroutine d5_ptr(array,m)
    use metadata_interfaces, metadata_address => getdp5ptr
    implicit none
    type(malloc_information_ptr), intent(in) :: m
    double precision, dimension(:,:,:,:,:), pointer, intent(inout) :: array
    !local variables
    integer(kind=8) :: iadd
    !allocate the array
    allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
         m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4),&
         m%lbounds(5):m%ubounds(5)+ndebug),stat=ierror)
    include 'allocate-inc.f90'
  end subroutine d5_ptr

  subroutine d5_ptr_free(array)
    use metadata_interfaces, metadata_address => getdp5ptr
    implicit none
    double precision, dimension(:,:,:,:,:), pointer, intent(inout) :: array
    include 'deallocate-inc.f90'
    nullify(array)
  end subroutine d5_ptr_free


end module dynamic_memory
