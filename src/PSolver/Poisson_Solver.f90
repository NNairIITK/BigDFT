!!****h* BigDFT/Poisson_Solver
!! NAME
!!    Poisson_Solver
!!
!! FUNCTION
!!    The module of the Poisson Solver.
!!    It must be used in the parent routine. 
!!
!! USAGE
!!    In the main routine in which the Poisson Solver is called
!!    1) The Poisson kernel must be declared as a pointer, then the 
!!       routine createKernel can be called. On exit, the kernel will be allocated and
!!       ready to use. See the documentation of the createKernel routine for more details
!!    2) The correct sizes for allocating the density/potential and the pot_ion arrays
!!       are given from the routine PS_dim4allocation (see routine documentation for details).
!!       Its argument MUST be in agreement with the arguments of the PSolver routine. 
!!       WARNING: No cross-check of the arguments is performed!
!!    3) The PSolver routine can then be called. On exit, the Hartree potential is computed
!!       and summed (following ixc value) to XC and external potential. 
!!       The input array is overwritten. Again, see routine documentation for details.
!!
!! WARNING
!!    This module REQUIRE the module of XC functional from ABINIT, defs_xc, which
!!    require defs_basis and defs_datatypes. 
!!    Such routines are provided inside the abinit directory of this bundle.
!!    They are based on XC functionals present in ABINIT 5.x
!!    If you want to use this Poisson Solver without the XC functionals, you can comment out
!!    the XC part in the PSolver routine
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
module Poisson_Solver

  private

  !calculate the allocation dimensions
  public :: PS_dim4allocation
  !routine that creates the kernel
  public :: createKernel
  !calculate the poisson solver
  public :: PSolver
  !calculate the allocation dimensions
  public :: P_FFT_dimensions, S_FFT_dimensions, F_FFT_dimensions

contains

  include 'PSolver_launch.f90'
  include 'Build_Kernel.f90'
  include 'PSolver_Base.f90'
  include 'xcenergy.f90'
  include '3Dgradient.f90'
  include 'fft3d.f90'
  include 'scaling_function.f90'

end module Poisson_Solver
