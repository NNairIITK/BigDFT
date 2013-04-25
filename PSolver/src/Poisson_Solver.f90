!> @file
!!    Define the module for Poisson Solver
!!
!! @author
!!    Luigi Genovese (February 2007)
!!    PSolverNC added by Anders Bergman, March 2008
!!    Copyright (C) 2002-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> The module of the Poisson Solver.
!!    It must be used in the parent routine. 
!! USAGE
!!    In the main routine in which the Poisson Solver is called
!!    -# The Poisson kernel must be declared as a pointer, then the 
!!       routine createKernel can be called. On exit, the kernel will be allocated and
!!       ready to use. See the documentation of the createKernel routine for more details
!!    -# The correct sizes for allocating the density/potential and the pot_ion arrays
!!       are given from the routine PS_dim4allocation (see routine documentation for details).
!!       Its argument MUST be in agreement with the arguments of the PSolver routine. 
!!       WARNING: No cross-check of the arguments is performed!
!!    -# The PSolver routine can then be called. On exit, the Hartree potential is computed
!!       and summed (following ixc value) to XC and external potential. 
!!       The input array is overwritten. Again, see routine documentation for details.
!!    -# QUICK INSTRUCTION FOR THE IMPATIENT:If you want to use the Poisson Solver in the 
!!       "ordinary" way, for a grid of dimensions nx,ny,nz and grid spacings hx,hy,hz, 
!!       just create the Kernel with
!!           call createKernel(geocode,nx,ny,nz,hx,hy,hz,14,0,1,kernel)
!!       where kernel is a pointer as described above; 
!!       geocode is 'F','S' or 'P' for Free, Surfaces of Periodic BC respectively.
!!       (Beware that for Surfaces BC the isolated direction is y!)
!!       After that you can calculate the potential with
!!           call PSolver(geocode,'G',0,1,nx,ny,nz,0,hx,hy,hz,&
!!                rhopot,kernel,fake_arr,energy,fake_exc,fake_vxc,0.d0,.false.,1)
!!       where:
!!         @param rhopot    is the density on input, and the electrostatic potential on output
!!                   (an array of dimension(nx,ny,nz))
!!         @param energy    is the result of @f$ 1/2 \int dx rho(x) potA @f$(x)
!!         @param fake_arr  is an array of dimension(1), untouched
!!         @param fake_*xc  values of the XC energies, automatically zero in that case
!!
!!       Any other changment of the arguments require reading of the documentation.
!!       See documentations of the Public routines
!! @warning
!!    This module REQUIRES the module of XC functional from ABINIT, defs_xc, which
!!    require defs_basis and defs_datatypes. 
!!    Such routines are provided inside the abinit directory of this bundle.
!!    They are based on XC functionals present in ABINIT 5.x
!!    If you want to use this Poisson Solver without the XC functionals, you can comment out
!!    the XC part in the PSolver routine
!!    Search for
!!

module Poisson_Solver
  use wrapper_linalg
  use wrapper_MPI
  ! TO BE REMOVED with f_malloc
  use m_profiling
  ! TO BE REMOVED with f_malloc

  implicit none

  private

  ! General precision, density and the potential types
  integer, parameter, public :: gp=kind(1.0d0)  !< general-type precision
  integer, parameter, public :: dp=kind(1.0d0)  !< density-type precision
  integer, parameter, public :: wp=kind(1.0d0)  !< potential-type precision
  ! Associated MPI precisions.
  integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypd=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION

  include 'configure.inc'

  !> Defines the fundamental structure for the kernel
  type, public :: coulomb_operator
     !variables with physical meaning
     integer :: itype_scf !< order of the ISF family to be used
     real(gp) :: mu !< inverse screening length for the Helmholtz Eq. (Poisson Eq. -> mu=0)
     character(len=1) :: geocode !< Code for the boundary conditions
     integer, dimension(3) :: ndims !< dimension of the box of the density
     real(gp), dimension(3) :: hgrids !<grid spacings in each direction
     real(gp), dimension(3) :: angrad !< angles in radiants between each of the axis
     real(dp), dimension(:), pointer :: kernel !< kernel of the Poisson Solver
     real(dp) :: work1_GPU,work2_GPU,k_GPU
     integer, dimension(5) :: plan
     integer, dimension(3) :: geo
     !variables with computational meaning
!     integer :: iproc_world !iproc in the general communicator
!     integer :: iproc,nproc
!     integer :: mpi_comm
     type(mpi_environment) :: mpi_env
     integer :: igpu !< control the usage of the GPU
     integer :: initCufftPlan
     integer :: keepGPUmemory
  end type coulomb_operator

  !calculate the allocation dimensions
  public :: PS_dim4allocation, PS_getVersion
  !routine that creates the kernel
  public :: pkernel_init, pkernel_set, pkernel_free
  !calculate the poisson solver
  public :: H_potential 
  !calculate the allocation dimensions
  public :: P_FFT_dimensions, S_FFT_dimensions, F_FFT_dimensions, W_FFT_dimensions, xc_dimensions

contains

  function pkernel_null() result(k)
    type(coulomb_operator) :: k
    k%itype_scf=0
    k%geocode='F'
    k%mu=0.0_gp
    k%ndims=(/0,0,0/)
    k%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    k%angrad=(/0.0_gp,0.0_gp,0.0_gp/)
    nullify(k%kernel)
    k%work1_GPU=0.d0
    k%work2_GPU=0.d0
    k%k_GPU=0.d0
    k%plan=(/0,0,0,0,0/)
    k%geo=(/0,0,0/)
    k%mpi_env=mpi_environment_null()
    k%igpu=0
  end function pkernel_null

  function PS_getVersion() result(str)
    character(len = 128) :: str

    write(str, "(A)") package_version
  end function PS_getVersion

  include 'PSolver_Main.f90'
  include 'createKernel.f90'
end module Poisson_Solver
