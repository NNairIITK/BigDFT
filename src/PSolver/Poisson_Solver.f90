module Poisson_Solver

  private

  !calculate the allocation dimensions
  public :: PS_dim4allocation
  !routine that creates the kernel
  public :: createKernel
  !calculate the poisson solver
  public :: PSolver

contains

  include 'PSolver_launch.f90'
  include 'Build_Kernel.f90'
  include 'PSolver_Base.f90'
  include 'xcenergy.f90'
  include '3Dgradient.f90'
  include 'fft3d.f90'
  include 'scaling_function.f90'

end module Poisson_Solver
