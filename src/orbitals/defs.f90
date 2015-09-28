!> @file
!!  File defining parameters for BigDFT package (processed by the build system)
!! @author
!!    Copyright (C) 2008-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!! @warning
!!   THIS FILE IS PROCESSED BY THE BUILD SYSTEM.

!> Modules which contains the low level definitions, as well as some profiling procedures
module module_defs
  use wrapper_MPI
  use f_precisions

  implicit none  

  private

  include 'configure.inc' !< Include variables set from configure.

  integer, public :: verbose=2    !< Verbosity of the output, control the level of writing (normal by default)

  ! General precision, density and the wavefunctions types
  integer, parameter, public :: gp=f_double!kind(1.0d0)  !< general-type precision
  integer, parameter, public :: dp=f_double!kind(1.0d0)  !< density-type precision
  integer, parameter, public :: wp=f_double!kind(1.0d0)  !< wavefunction-type precision
  integer, parameter, public :: tp=f_double!kind(1.0d0)  !< DIIS precision (single in this context, if double is only for non-regression)

  !> Define type of data for MPI
  integer, parameter, public :: mpidtypw=MPI_DOUBLE_PRECISION
  integer, parameter, public :: mpidtypd=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
  !integer, parameter :: mpidtypw=MPI_REAL,mpidtypd=MPI_REAL !in case of single precision

!!$  !> Flag for GPU computing, if CUDA libraries are present
!!$  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
!!$  !! this value can be changed in the read_input_variables routine
!!$  logical :: GPUconv=.false.,GPUshare=.true.

  !> Flag for GPU computing, if OpenCL libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical, parameter, public :: ASYNCconv=.true.

  !> Logical parameter for the projectors application strategy (true for distributed way)
  !! if the projector allocation passes the memorylimit this is switched to true
  !! inside localize_projectors routines
  logical, public :: DistProjApply=.true. !<then copied as a element of the nlpsp structure

  !> experimental variables to test the add of new functionalities
  logical :: experimental_modulebase_var_onlyfion=.false.

  !> Evergreen
  real(dp), parameter, public :: pi_param=3.141592653589793238462643383279502884197_dp

  !> Error codes, to be documented little by little
  integer, save, public :: BIGDFT_RUNTIME_ERROR                   !< Error during runtime
  integer, save, public :: BIGDFT_MPI_ERROR                       !< See error definitions below
  integer, save, public :: BIGDFT_LINALG_ERROR                    !< To be moved to linalg wrappers
  integer, save, public :: BIGDFT_INPUT_VARIABLES_ERROR           !< Problems in parsing or in consistency of input variables
  integer, save, public :: BIGDFT_INPUT_FILE_ERROR                !< The file does not exist!

  !> Code constants.
  !real(gp), parameter :: UNINITIALISED = -123456789._gp

  private :: f_double,f_simple,f_long,f_short,f_integer

  !interface for uninitialized variable
  interface UNINITIALIZED
     module procedure uninitialized_dbl,uninitialized_int,uninitialized_real,uninitialized_long, uninitialized_logical
  end interface

  public :: UNINITIALIZED

  contains

    pure function uninitialized_int(one) 
      implicit none
      integer(kind = 4), intent(in) :: one
      integer(kind = 4) :: uninitialized_int
      integer :: foo
      foo = kind(one)
      uninitialized_int=-123456789
    end function uninitialized_int

    pure function uninitialized_long(one) 
      implicit none
      integer(kind = 8), intent(in) :: one
      integer(kind = 8) :: uninitialized_long
      integer :: foo
      foo = kind(one)
      uninitialized_long=-123456789
    end function uninitialized_long

    pure function uninitialized_real(one) 
      implicit none
      real(kind=4), intent(in) :: one
      real(kind=4) :: uninitialized_real
      integer :: foo
      foo = kind(one)
      uninitialized_real=-123456789.e0
    end function uninitialized_real

    pure function uninitialized_dbl(one) 
      implicit none
      real(kind=8), intent(in) :: one
      real(kind=8) :: uninitialized_dbl
      integer :: foo
      foo = kind(one)
      uninitialized_dbl=-123456789.d0
    end function uninitialized_dbl

    pure function uninitialized_logical(one) 
      implicit none
      logical, intent(in) :: one
      logical :: uninitialized_logical
      integer :: foo
      foo = kind(one)
      uninitialized_logical=.false.
    end function uninitialized_logical

end module module_defs
