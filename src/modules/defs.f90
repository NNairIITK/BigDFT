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
  use wrapper_linalg

  implicit none  

  include 'configure.inc' !< Include variables set from configure.

  integer :: verbose=2    !< Verbosity of the output, control the level of writing (minimal by default)

  ! General precision, density and the wavefunctions types
  integer, parameter :: gp=kind(1.0d0)  !< general-type precision
  integer, parameter :: dp=kind(1.0d0)  !< density-type precision
  integer, parameter :: wp=kind(1.0d0)  !< wavefunction-type precision
  integer, parameter :: tp=kind(1.0d0)  !< DIIS precision (single in this context, if double is only for non-regression)

  !> Define type of data for MPI
  integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypd=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
  !integer, parameter :: mpidtypw=MPI_REAL,mpidtypd=MPI_REAL !in case of single precision

  !> Flag for GPU computing, if CUDA libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical :: GPUconv=.false.,GPUshare=.true.

  !> Flag for GPU computing, if OpenCL libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical, parameter :: ASYNCconv=.true.

  !> Logical parameter for the projectors application strategy (true for distributed way)
  !! if the projector allocation passes the memorylimit this is switched to true
  !! inside localize_projectors routines
  logical :: DistProjApply=.true. !<then copied as a element of the nlpsp structure

  !> experimental variables to test the add of new functionalities
  logical :: experimental_modulebase_var_onlyfion=.false.

  type(mpi_environment) :: bigdft_mpi !< Contains all data needed for MPI processes

  !> Physical constants.
  real(gp), parameter :: Bohr_Ang = 0.52917721092_gp                    !< 1 AU in angstroem
  real(gp), parameter :: Ha_cmm1=219474.6313705_gp                      !< 1 Hartree, in cm^-1 (from abinit 5.7.x)
  real(gp), parameter :: Ha_eV=27.21138505_gp                           !< 1 Hartree, in eV
  real(gp), parameter :: eV_Ha=3.674932379e-2_gp                        !< 1 ev, in Hartree
  real(gp), parameter :: Ha_K=315774.65_gp                              !< 1 Hartree, in Kelvin
  real(gp), parameter :: Ha_THz=6579.683920722_gp                       !< 1 Hartree, in THz
  real(gp), parameter :: Ha_J=4.35974394d-18                            !< 1 Hartree, in J
  real(gp), parameter :: e_Cb=1.602176487d-19                           !< minus the electron charge, in Coulomb
  real(gp), parameter :: kb_HaK=8.617343d-5/Ha_eV                       !< Boltzmann constant in Ha/K
  real(gp), parameter :: amu_emass=1.660538782e-27_gp/9.10938215e-31_gp !< 1 atomic mass unit, in electronic mass
  real(gp), parameter :: AU_GPa=29421.010901602753_gp                   !< 1 Ha/Bohr^3 in GPa
  real(gp), parameter :: Radian_Degree = 57.29577951308232087679_gp     !< 1 radian in degrees
  real(gp), parameter :: eVAng_HaBohr = Bohr_Ang*eV_Ha                  !< convert forces from eV/Angstroem to hartree/bohr
  real(gp), parameter :: kcalMol_Ha = 0.001593601437458137_gp        !< from kcal_th/mol to hartree
                                                                     !!(thermochemical calorie used in amber: 1cal_th=4.184J)
                                                                     !!also see: http://archive.ambermd.org/201009/0039.html
  real(gp), parameter :: kcalMolAng_HaBohr =0.0008432975639921999_gp !<convert forces from kcal_th/mol/angstrom to hartree/bohr


  !> Evergreen
  real(dp), parameter :: pi_param=3.141592653589793238462643383279502884197_dp

  !> Error codes, to be documented little by little
  integer, save :: BIGDFT_RUNTIME_ERROR                   !< Error during runtime
  integer, save :: BIGDFT_MPI_ERROR                       !< See error definitions below
  integer, save :: BIGDFT_LINALG_ERROR                    !< To be moved to linalg wrappers
  integer, save :: BIGDFT_INPUT_VARIABLES_ERROR           !< Problems in parsing or in consistency of input variables
  integer, save :: BIGDFT_INPUT_FILE_ERROR                !< The file does not exist!

  !> Code constants.
  !real(gp), parameter :: UNINITIALISED = -123456789._gp

  !interface for uninitialized variable
  interface UNINITIALIZED
     module procedure uninitialized_dbl,uninitialized_int,uninitialized_real,uninitialized_long, uninitialized_logical
  end interface

  interface safe_exp
     module procedure safe_dexp
  end interface safe_exp

!!$  interface
!!$     subroutine bigdft_utils_flush(unit)
!!$       integer, intent(in) :: unit
!!$     end subroutine bigdft_utils_flush
!!$  end interface

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

end module module_defs
