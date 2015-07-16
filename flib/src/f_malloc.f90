!> @file
!! Manage dynamic memory allocation control structures
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!> Module used to manage memory allocations control structures.
!! This module has to be intended as a submodule of dynamic_memory module
module module_f_malloc

  use dictionaries, only: f_loc,f_err_throw,f_err_raise

  !>global parameter of the module telling if the profile has to be activated
  !! this parameter can be modified only by dynamic memory module
  integer(kind=4), parameter :: f_malloc_namelen=32          !< length of the character variables
  integer, parameter :: max_rank=7          !< maximum rank in fortran

  !to be initialized in the dynamic_memory module
  integer, save :: ERR_INVALID_MALLOC

  logical, save, public :: f_malloc_default_profiling=.true.
  character(len=f_malloc_namelen), save, public :: f_malloc_routine_name=repeat(' ',f_malloc_namelen)

  !> Structure needed to allocate an allocatable array
  type, public :: malloc_information_all
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the array
     integer(kind=4), dimension(max_rank) :: shape   !< shape of the structure 
     integer(kind=4), dimension(max_rank) :: lbounds !< lower bounds
     integer(kind=4), dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
     
  end type malloc_information_all

  !> Structure needed to allocate an allocatable array of string of implicit length (for non-2003 compilers)
  type, public :: malloc_information_str_all
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the array
     integer(kind=4) :: len                  !< length of the character
     integer(kind=4), dimension(max_rank) :: shape   !< shape of the structure 
     integer(kind=4), dimension(max_rank) :: lbounds !< lower bounds
     integer(kind=4), dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
  end type malloc_information_str_all

  !> Structure needed to allocate a pointer
  type, public :: malloc_information_ptr
     logical :: ptr                          !< just to make the structures different, to see if needed
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the pointer
     integer(kind=4), dimension(max_rank) :: shape   !< shape of the structure 
     integer(kind=4), dimension(max_rank) :: lbounds !< lower bounds
     integer(kind=4), dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
  end type malloc_information_ptr

  !> Structure needed to allocate a pointer of string of implicit length (for non-2003 complilers)
  type, public :: malloc_information_str_ptr
     logical :: ptr                          !< just to make the structures different, to see if needed
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the pointer
     integer(kind=4) :: len                  !< length of the character
     integer(kind=4), dimension(max_rank) :: shape   !< shape of the structure 
     integer(kind=4), dimension(max_rank) :: lbounds !< lower bounds
     integer(kind=4), dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
  end type malloc_information_str_ptr

  type, public :: array_bounds
     integer(kind=4) :: nlow  !<lower bounds
     integer(kind=4) :: nhigh !<higher bounds
  end type array_bounds

  interface operator(.to.)
     module procedure f_array_bounds_intint, f_array_bounds_longlong
     module procedure f_array_bounds_intlong, f_array_bounds_longint
  end interface

  interface nullify_malloc_information
     module procedure nullify_malloc_information_all
     module procedure nullify_malloc_information_ptr
     module procedure nullify_malloc_information_str_all
     module procedure nullify_malloc_information_str_ptr
  end interface 

  !> Fake structure needed to document common arguments of the module
  type, private :: doc
     !> integer indicating the size of the array. Can be specified only if the target 
     !! array has to be of rank one.
     integer :: size
     !> identification of the allocation. Usually called as the name of the allocatable variable
     character(len=1) :: id
  end type doc

  !> Structure to perform heap allocation. Can be used to allocate allocatable arrays of any intrinsic
  !! type, kind and rank.
  !! The assignment of an allocatable array to this structure will allocate the array according to the
  !! specifications. It will also store the information about the allocation in an iternal database
  !! so that memory leaks and allocation problems can be more easily detected.
  !! @param mapname  @copydoc doc::size
  !! @param id       @copydoc doc::id
  interface f_malloc
     module procedure f_malloc,f_malloci_simple
     module procedure f_mallocli_simple
     module procedure f_malloc_bounds,f_malloc_bound
     !here also the procedures for the copying of arrays have to be defined
     module procedure f_malloc_i2,f_malloc_d2
     module procedure f_malloc_d1,f_malloc_i3
     module procedure f_malloc_d3,f_malloc_d4
     module procedure f_malloc_r1,f_malloc_r2,f_malloc_r3,f_malloc_r4
  end interface

  interface f_malloc0
     module procedure f_malloc0
     !!module procedure f_mallocli0
     module procedure f_malloci0_simple
     module procedure f_mallocli0_simple
     module procedure f_malloc0_bounds,f_malloc0_bound
  end interface

  interface f_malloc_ptr
     module procedure f_malloc_ptr,f_malloci_ptr_simple
     module procedure f_mallocli_ptr_simple
     module procedure f_malloc_ptr_bounds,f_malloc_ptr_bound
     !module procedure f_malloc_ptr_i2,f_malloc_ptr_i3
     module procedure f_malloc_ptr_i2_sp,f_malloc_ptr_i3_sp
     module procedure f_malloc_ptr_i4_sp
     !module procedure f_malloc_ptr_d1,f_malloc_ptr_d2
     !module procedure f_malloc_ptr_d3,f_malloc_ptr_d4
     module procedure f_malloc_ptr_d1_sp,f_malloc_ptr_d2_sp
     module procedure f_malloc_ptr_d3_sp,f_malloc_ptr_d4_sp
  end interface

  interface f_malloc0_ptr
     module procedure f_malloc0_ptr
     !!module procedure f_mallocli0_ptr
     module procedure f_malloci0_ptr_simple
     module procedure f_mallocli0_ptr_simple
     module procedure f_malloc0_ptr_bounds,f_malloc0_ptr_bound
  end interface

  interface f_malloc_str
     module procedure f_malloci_str,f_malloci_str_simple,f_mallocli_str_simple
     module procedure f_malloc_istr_bounds,f_malloc_istr_bound
     !!module procedure f_mallocli_str
     !!module procedure f_malloc_istr_simple
     !!module procedure f_malloc_listr_bounds
     !!module procedure f_malloc_listr_bound
  end interface

  interface f_malloc0_str
     module procedure f_malloci0_str,f_malloci0_str_simple,f_mallocli0_str_simple
     module procedure f_malloc0_istr_bounds,f_malloc0_istr_bound
     !!module procedure f_mallocli0_str
     !!module procedure f_malloc0_listr_bounds,f_malloc0_listr_bound
  end interface

  interface f_malloc_str_ptr
     module procedure f_malloci_str_ptr
     !!module procedure f_mallocli_str_ptr
     module procedure f_malloci_str_ptr_simple
     module procedure f_mallocli_str_ptr_simple
     module procedure f_malloc_str_ptr_bound
     module procedure f_malloc_str_ptr_bounds
  end interface

  interface f_malloc0_str_ptr
     module procedure f_malloci0_str_ptr
     !!module procedure f_mallocli0_str_ptr
     module procedure f_malloci0_str_ptr_simple
     module procedure f_malloc0_str_ptr_bounds,f_malloc0_str_ptr_bound
     !!module procedure f_malloc0_listr_ptr
     module procedure f_mallocli0_str_ptr_simple
  end interface

contains
  
  elemental pure function f_array_bounds_intint(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(kind=4), intent(in) :: nlow
    integer(kind=4), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,kind=4)
    f_array_bounds%nhigh=int(nhigh,kind=4)
  end function f_array_bounds_intint

  elemental pure function f_array_bounds_longlong(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(kind=8), intent(in) :: nlow
    integer(kind=8), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,kind=4)
    f_array_bounds%nhigh=int(nhigh,kind=4)
  end function f_array_bounds_longlong

  elemental pure function f_array_bounds_intlong(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(kind=4), intent(in) :: nlow
    integer(kind=8), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,kind=4)
    f_array_bounds%nhigh=int(nhigh,kind=4)
  end function f_array_bounds_intlong

  elemental pure function f_array_bounds_longint(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(kind=8), intent(in) :: nlow
    integer(kind=4), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,kind=4)
    f_array_bounds%nhigh=int(nhigh,kind=4)
  end function f_array_bounds_longint

  pure subroutine nullify_malloc_information_all(m)
    implicit none
    type(malloc_information_all), intent(out) :: m
    include 'f_malloc-null-inc.f90'
  end subroutine nullify_malloc_information_all

  pure subroutine nullify_malloc_information_ptr(m)
    implicit none
    type(malloc_information_ptr), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%ptr=.true.
  end subroutine nullify_malloc_information_ptr

  pure subroutine nullify_malloc_information_str_all(m)
    implicit none
    type(malloc_information_str_all), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%len=0
  end subroutine nullify_malloc_information_str_all
  
  pure subroutine nullify_malloc_information_str_ptr(m)
    implicit none
    type(malloc_information_str_ptr), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%len=0
    m%ptr=.true.
  end subroutine nullify_malloc_information_str_ptr
  

  !---routines for low-level dynamic memory handling
  
  !> For rank-1 arrays
  pure function f_malloci_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloci-simple-inc.f90'
  end function f_malloci_simple
  pure function f_mallocli_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_mallocli-simple-inc.f90'
  end function f_mallocli_simple
  !> For rank-1 arrays
  pure function f_malloci0_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloci-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_malloci0_simple
  pure function f_mallocli0_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_mallocli-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_mallocli0_simple
  !> For rank-1 arrays
  pure function f_malloci_ptr_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloci-simple-inc.f90'
  end function f_malloci_ptr_simple
  pure function f_mallocli_ptr_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_mallocli-simple-inc.f90'
  end function f_mallocli_ptr_simple
  !> For rank-1 arrays
  pure function f_malloci0_ptr_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloci-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_malloci0_ptr_simple
  pure function f_mallocli0_ptr_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_mallocli-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_mallocli0_ptr_simple
  !> For rank-1 arrays
  pure function f_malloci_str_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloci_str_simple
  pure function f_mallocli_str_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=int(length,kind=4)
  end function f_mallocli_str_simple
  !> For rank-1 arrays
  pure function f_malloci0_str_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloci0_str_simple
  pure function f_mallocli0_str_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_mallocli0_str_simple
  !> For rank-1 arrays
  pure function f_malloci_str_ptr_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloci_str_ptr_simple
  pure function f_mallocli_str_ptr_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=int(length,kind=4)
  end function f_mallocli_str_ptr_simple
  !> For rank-1 arrays
  pure function f_malloci0_str_ptr_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloci0_str_ptr_simple
  pure function f_mallocli0_str_ptr_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_mallocli0_str_ptr_simple


  !> For rank-1 arrays, with bounds
  pure function f_malloc_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bound-inc.f90'
  end function f_malloc_bound
  !>for rank-1 arrays, with boundaries
  pure function f_malloc0_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bound-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_bound
  !> For rank-1 arrays
  pure function f_malloc_ptr_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bound-inc.f90'
  end function f_malloc_ptr_bound
  !> For rank-1 arrays
  pure function f_malloc0_ptr_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bound-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr_bound
  !> For rank-1 arrays, with bounds
  pure function f_malloc_istr_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloc_istr_bound
  !!pure function f_malloc_listr_bound(length,bounds,id,routine_id,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bound-inc.f90'
  !!  m%len=int(length,kind=4)
  !!end function f_malloc_listr_bound
  !>for rank-1 arrays, with boundaries
  pure function f_malloc0_istr_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloc0_istr_bound
  !!pure function f_malloc0_listr_bound(length,bounds,id,routine_id,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bound-inc.f90'
  !!  m%len=int(length,kind=4)
  !!  m%put_to_zero=.true.
  !!end function f_malloc0_listr_bound
  !> For rank-1 arrays
  pure function f_malloc_str_ptr_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloc_str_ptr_bound
  !> For rank-1 arrays
  pure function f_malloc0_str_ptr_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr_bound
  pure function f_malloc0_listr_ptr_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloc0_listr_ptr_bound


  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bounds-inc.f90'
  end function f_malloc_bounds
  pure function f_malloc0_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bounds-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_ptr_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bounds-inc.f90'
  end function f_malloc_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc0_ptr_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bounds-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_istr_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloc_istr_bounds
  !!pure function f_malloc_listr_bounds(length,bounds,id,routine_id,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bounds-inc.f90'
  !!  m%len=int(length,kind=4)
  !!end function f_malloc_listr_bounds
  pure function f_malloc0_istr_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloc0_istr_bounds
  !!pure function f_malloc0_listr_bounds(length,bounds,id,routine_id,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bounds-inc.f90'
  !!  m%len=int(length,kind=4)
  !!  m%put_to_zero=.true.
  !!end function f_malloc0_listr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_str_ptr_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloc_str_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc0_str_ptr_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr_bounds

  !> Define the allocation information for arrays of different rank
  function f_malloc(sizes,id,routine_id,lbounds,ubounds,profile,src,src_ptr) result(m)
    implicit none
    !the integer array src is here added to avoid problems in resolving the ambiguity with f_malloc_src
    integer, dimension(:), intent(in), optional :: src
    integer, dimension(:), pointer, intent(in), optional :: src_ptr
    type(malloc_information_all) :: m
    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
    !local variables
    integer :: i
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-check-inc.f90'
    include 'f_malloc-extra-inc.f90'
  end function f_malloc
  !> define the allocation information for  arrays of different rank
  function f_malloc0(sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-total-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0
  !!function f_mallocli0(sizes,id,routine_id,lbounds,ubounds,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_all) :: m
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%put_to_zero=.true.
  !!end function f_mallocli0
  !> Define the allocation information for  arrays of different rank
  function f_malloc_ptr(sizes,id,routine_id,lbounds,ubounds,profile,src,src_ptr) result(m)
    implicit none
    !the integer array src is here added to avoid problems in resolving the ambiguity
    integer, dimension(:), intent(in), optional :: src
    integer, dimension(:), pointer, intent(in), optional :: src_ptr
    type(malloc_information_ptr) :: m
    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
    !local variables
    integer :: i
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-check-inc.f90'
    include 'f_malloc-extra-inc.f90'
  end function f_malloc_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloc0_ptr(sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-total-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr
  !!function f_mallocli0_ptr(sizes,id,routine_id,lbounds,ubounds,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_ptr) :: m
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%put_to_zero=.true.
  !!end function f_mallocli0_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloci_str(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-total-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloci_str
  !!function f_mallocli_str(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%len=int(length,kind=4)
  !!end function f_mallocli_str
  !> define the allocation information for  arrays of different rank
  function f_malloci0_str(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-total-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloci0_str
  !!function f_mallocli0_str(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%len=int(length,kind=4)
  !!  m%put_to_zero=.true.
  !!end function f_mallocli0_str
  !> Define the allocation information for  arrays of different rank
  function f_malloci_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-total-inc.f90'
    m%len=int(length,kind=4)
  end function f_malloci_str_ptr
  !!function f_mallocli_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_ptr) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%len=int(length,kind=4)
  !!end function f_mallocli_str_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloci0_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer(kind=4), intent(in) :: length
    include 'f_malloci-total-inc.f90'
    m%len=int(length,kind=4)
    m%put_to_zero=.true.
  end function f_malloci0_str_ptr
  !!function f_mallocli0_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
  !!  implicit none
  !!  type(malloc_information_str_ptr) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%len=int(length,kind=4)
  !!  m%put_to_zero=.true.
  !!end function f_mallocli0_str_ptr

  function f_malloc_d1(src,lbounds,ubounds,id,routine_id,profile) result(m)
    implicit none
    double precision, dimension(:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d1

  function f_malloc_d2(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    double precision, dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d2

  function f_malloc_d3(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    double precision, dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d3

  function f_malloc_d4(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    double precision, dimension(:,:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d4

  function f_malloc_r1(src,lbounds,ubounds,id,routine_id,profile) result(m)
    implicit none
    real, dimension(:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r1

  function f_malloc_r2(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    real, dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r2

  function f_malloc_r3(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    real, dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r3

  function f_malloc_r4(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    real, dimension(:,:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r4

  function f_malloc_i2(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    integer, dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_i2

  function f_malloc_i3(src,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    integer, dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_i3

!!$  function f_malloc_ptr_i2(src,id,routine_id,sizes,lbounds,ubounds,profile) result(m)
!!$    implicit none
!!$    integer, dimension(:,:), intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_i2

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_i2_sp(src_ptr,id,routine_id,profile) result(m)
    implicit none
    integer, dimension(:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_i2_sp

!!$  function f_malloc_ptr_i3(src,id,routine_id,sizes,lbounds,ubounds,profile) result(m)
!!$    implicit none
!!$    integer, dimension(:,:,:), intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_i3

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_i3_sp(src_ptr,id,routine_id,profile) result(m)
    implicit none
    integer, dimension(:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_i3_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_i4_sp(src_ptr,id,routine_id,profile) result(m)
    implicit none
    integer, dimension(:,:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_i4_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d1_sp(src_ptr,id,routine_id,profile) result(m)
    implicit none
    double precision, dimension(:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d1_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d2_sp(src_ptr,id,routine_id,profile) result(m)
    implicit none
    double precision, dimension(:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d2_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d3_sp(src_ptr,id,routine_id,profile) result(m)
    implicit none
    double precision, dimension(:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d3_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d4_sp(src_ptr,id,routine_id,profile) result(m)
    implicit none
    double precision, dimension(:,:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d4_sp
  
!!$  function f_malloc_ptr_d1(id,src,routine_id,sizes,lbounds,ubounds,profile) result(m)
!!$    implicit none
!!$    double precision, dimension(:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d1

!!$  function f_malloc_ptr_d2(id,src,routine_id,sizes,lbounds,ubounds,profile) result(m)
!!$    implicit none
!!$    double precision, dimension(:,:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d2
!!$
!!$  function f_malloc_ptr_d3(id,src,routine_id,sizes,lbounds,ubounds,profile) result(m)
!!$    implicit none
!!$    double precision, dimension(:,:,:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d3
!!$
!!$  function f_malloc_ptr_d4(id,src,routine_id,sizes,lbounds,ubounds,profile) result(m)
!!$    implicit none
!!$    double precision, dimension(:,:,:,:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d4
  
end module module_f_malloc
