!> @file
!!    Modulefile for the definition of the basic structures
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module PStypes
  use f_enums
  use wrapper_MPI
  use PSbase
  use environment, only: cavity_data,cavity_default
  implicit none

  private
  
  !>Defines the internal information for application of the FFT between the kernel and the 
  !!density
  type, public :: FFT_metadata
     integer :: m1,m2,m3 !<original real dimension, with m2 in z direction and and m3 in y direction
     integer :: n1,n2,n3 !<dimension of the FFT operation, taking into account the zero padding if needed
     integer :: md1,md2,md3 !< Dimension of the real unpadded space, 
     !!md2 is further enlarged to be a multiple of number of processes
     integer :: nd1,nd2,nd3 !<fourier dimensions for which the kernel is injective,
     !!                formally 1/8 of the fourier grid. Here the dimension nd3 is
     !!                enlarged to be a multiple of nproc
     integer :: istart,iend,n3p !<start, endpoints and number of planes of the given processor
     real(dp) :: scal !<factor to rescale the solver such that the FFT is unitary, divided by 4pi
  end type FFT_metadata

  !> define the work arrays needed for the treatment of the 
  !!Poisson Equation. Not all of them are allocated, the actual memory usage 
  !!depends on the treatment
  type, public :: PS_workarrays
     !> logaritmic derivative of the dielectric function,
     !! to be used in the case of Polarization Iteration method
     real(dp), dimension(:,:,:,:), pointer :: dlogeps
     !> inverse of the dielectric function
     !! in the case of Polarization Iteration method
     !! inverse of the square root of epsilon
     !! in the case of the Preconditioned Conjugate Gradient
     real(dp), dimension(:,:), pointer :: oneoeps
     !> correction term, given in terms of the multiplicative factor of nabla*eps*nabla
     !! to be used for Preconditioned Conjugate Gradient 
     real(dp), dimension(:,:), pointer :: corr
     !> inner rigid cavity to be integrated in the sccs method to avoit inner
     !! cavity discontinuity due to near-zero edens near atoms
     real(dp), dimension(:,:), pointer :: epsinnersccs
     real(dp), dimension(:,:,:), pointer :: zf
     !> input guess vectors to be preserved for future use
     real(dp), dimension(:,:), pointer :: rho_old,res_old,pot_old
     !> Polarization charge vector for print purpose only.
     real(dp), dimension(:,:), pointer :: pol_charge
     integer(f_address) :: work1_GPU,work2_GPU,k_GPU !<addresses for the GPU memory 
     integer(f_address) :: p_GPU,q_GPU,r_GPU,x_GPU,z_GPU,oneoeps_GPU,corr_GPU!<addresses for the GPU memory 
     !> GPU scalars. Event if they are scalars of course their address is needed
     integer(f_address) :: alpha_GPU, beta_GPU, kappa_GPU, beta0_GPU
  end type PS_workarrays

  !> Defines the fundamental structure for the kernel
  type, public :: coulomb_operator
     !variables with physical meaning
     integer :: itype_scf             !< Order of the ISF family to be used
     real(gp) :: mu                   !< Inverse screening length for the Helmholtz Eq. (Poisson Eq. -> mu=0)
     !> geocode is used in all the code to specify the boundary conditions (BC) the problem:
     !!          - 'F' free BC, isolated systems.
     !!                The program calculates the solution as if the given density is
     !!                "alone" in R^3 space.
     !!          - 'S' surface BC, isolated in y direction, periodic in xz plane                
     !!                The given density is supposed to be periodic in the xz plane,
     !!                so the dimensions in these direction mus be compatible with the FFT
     !!                Beware of the fact that the isolated direction is y!
     !!          - 'P' periodic BC.
     !!                The density is supposed to be periodic in all the three directions,
     !!                then all the dimensions must be compatible with the FFT.
     !!                No need for setting up the kernel (in principle for Plane Waves)
     !!          - 'W' Wires BC.
     !!                The density is supposed to be periodic in z direction, 
     !!                which has to be compatible with the FFT.
     !!          - 'H' Helmholtz Equation Solver
     character(len=1) :: geocode
     !> method of embedding in the environment
     !!          - 'VAC' Poisson Equation in vacuum. Default case.
     !!          - 'PCG' Generalized Poisson Equation, Preconditioned Conjugate Gradient
     !!          - 'PI'  Generalized Poisson Equation, Polarization Iteration method
     !character(len=3) :: method 
     !! this represents the information for the equation and the algorithm to be solved
     !! this enumerator contains the algorithm and has the attribute associated to the 
     !! type of cavity to be used
     type(f_enumerator) :: method
     integer, dimension(3) :: ndims   !< dimension of the box of the density
     real(gp), dimension(3) :: hgrids !<grid spacings in each direction
     real(gp), dimension(3) :: angrad !< angles in radiants between each of the axis
     type(cavity_data) :: cavity !< description of the cavity for the dielectric medium
     real(dp), dimension(:), pointer :: kernel !< kernel of the Poisson Solver
     integer, dimension(5) :: plan
     integer, dimension(3) :: geo
     !>workarrays for the application of the Solver. Might have different
     !!memory footprints dependently of the treatment.
     type(PS_workarrays) :: w
     !variables with computational meaning
     type(mpi_environment) :: mpi_env !< complete environment for the POisson Solver
     type(mpi_environment) :: inplane_mpi,part_mpi !<mpi_environment for internal ini-plane parallelization
     type(FFT_metadata) :: grid !<dimensions of the FFT grid associated to this kernel
     integer :: igpu !< control the usage of the GPU
     integer :: gpuPCGRed !< control if GPU can be used for PCG reductions
     integer :: initCufftPlan
     integer :: keepGPUmemory
     integer :: keepzf
     !parameters for the iterative methods
     !> Order of accuracy for derivatives into ApplyLaplace subroutine = Total number of points at left and right of the x0 where we want to calculate the derivative.
     integer :: nord
     integer :: max_iter !< maximum number of convergence iterations
     real(dp) :: minres !< convergence criterion for the iteration
     real(dp) :: PI_eta !<parameter for the update of PI iteration

     integer, dimension(:), pointer :: counts !<array needed to gather the information of the poisson solver
     integer, dimension(:), pointer :: displs !<array needed to gather the information of the poisson solver
     integer, dimension(:), pointer :: rhocounts !<array needed to gather the information of the poisson solver on multiple gpus
     integer, dimension(:), pointer :: rhodispls !<array needed to gather the information of the poisson solver on multiple gpus
  end type coulomb_operator

  !> define the energy terms for the Poisson Operator and Generalized applications
  type, public :: PSolver_energies
     !> hartree energy, defined as the @f$\int \rho(\mathbf{r}) V(\mathbf{r}) \mathrm d r $ @f, with @f$\rho$@f being the 
     !! input density and @f$V@f the potential defined by this density
     !! in the case when rho_ion is passed, the electrostatic contribution is only filled
     real(gp) :: hartree 
     !> electrostatic energy, defined as the hartree energy but with @f$\rho$@f and @f$V$@f coming from @f$\rho + \rho_\text{ion}$$f
     !! the hartree energy can be obtained by subtraction with the potential energy terms
     real(gp) :: elec 
     !> Energetic term coming from the @f$\int \rho V_extra$@f, in the case of a @f$\rho$@f-dependent cavity.
     !! clearly this term is calculated only if the potential is corrected. When the cavity is fixed, the eVextra is zero.
     real(gp) :: eVextra
     !> surface and volume term
     real(gp) :: cavitation
     !> stress tensor, to be calculated when calculate_strten is .true.
     real(gp), dimension(6) :: strten
  end type PSolver_energies

  !>Datatype defining the mode for the running of the Poisson solver
  type, public :: PSolver_options
     !> @copydoc poisson_solver::doc::datacode
     character(len=1) :: datacode
     !> integer variable setting the verbosity, from silent (0) to high
     integer :: verbosity_level
     !> add pot_ion to the final potential
     logical :: add_pot_ion
     !> add rho_ion to the initial density potential
     logical :: add_rho_ion
     !> if .true., and the cavity is set to 'sccs' attribute, then the epsilon is updated according to rhov
     logical :: update_cavity
     !> if .true. calculate the stress tensor components.
     !! The stress tensor calculation are so far validated only for orthorhombic cells and vacuum treatments
     logical :: calculate_strten
     !> Use the input guess procedure for the solution in the case of a electrostatic medium
     !! this value is ignored in the case of a vacuum solver
     logical :: use_input_guess
     !> Total integral on the supercell of the final potential on output
     !! clearly meaningful only for Fully periodic BC, ignored in the other cases
     real(gp) :: potential_integral
  end type PSolver_options

  public :: pkernel_null,PSolver_energies_null,pkernel_free

contains

  pure function PSolver_energies_null() result(e)
    implicit none
    type(PSolver_energies) :: e
    e%hartree    =0.0_gp
    e%elec       =0.0_gp
    e%eVextra    =0.0_gp
    e%cavitation =0.0_gp
    e%strten     =0.0_gp
  end function PSolver_energies_null

  pure function FFT_metadata_null() result(d)
    implicit none
    type(FFT_metadata) :: d
    d%m1=0
    d%m2=0
    d%m3=0
    d%n1=0
    d%n2=0
    d%n3=0
    d%md1=0
    d%md2=0
    d%md3=0
    d%nd1=0
    d%nd2=0
    d%nd3=0
    d%istart=0
    d%iend=0
    d%n3p=0
    d%scal=0.0_dp
  end function FFT_metadata_null

  pure subroutine nullify_work_arrays(w)
    use f_utils, only: f_zero
    implicit none
    type(PS_workarrays), intent(out) :: w
    nullify(w%dlogeps)
    nullify(w%oneoeps)
    nullify(w%corr)
    nullify(w%epsinnersccs)
    nullify(w%pol_charge)
    nullify(w%rho_old)
    nullify(w%pot_old)
    nullify(w%res_old)
    nullify(w%zf)
    call f_zero(w%work1_GPU)
    call f_zero(w%work2_GPU)
    call f_zero(w%k_GPU)
    call f_zero(w%p_GPU)
    call f_zero(w%q_GPU)
    call f_zero(w%r_GPU)
    call f_zero(w%x_GPU)
    call f_zero(w%z_GPU)
    call f_zero(w%oneoeps_GPU)
    call f_zero(w%corr_GPU)
    call f_zero(w%alpha_GPU)
    call f_zero(w%beta_GPU)
    call f_zero(w%kappa_GPU)
    call f_zero(w%beta0_GPU)
  end subroutine nullify_work_arrays

  pure function pkernel_null() result(k)
    implicit none
    type(coulomb_operator) :: k
    k%itype_scf=0
    k%geocode='F'
    call nullify_f_enum(k%method)
    k%cavity=cavity_default()
    k%mu=0.0_gp
    k%ndims=(/0,0,0/)
    k%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    k%angrad=(/0.0_gp,0.0_gp,0.0_gp/)
    nullify(k%kernel)
    k%plan=(/0,0,0,0,0/)
    k%geo=(/0,0,0/)
    call nullify_work_arrays(k%w)
    k%mpi_env=mpi_environment_null()
    k%inplane_mpi=mpi_environment_null()
    k%part_mpi=mpi_environment_null()
    k%grid=FFT_metadata_null()
    k%igpu=0
    k%initCufftPlan=0
    k%keepGPUmemory=1
    k%keepzf=1
    k%nord=0
    k%max_iter=0
    k%PI_eta=0.0_dp
    k%minres=0.0_dp
    nullify(k%counts)
    nullify(k%displs)
  end function pkernel_null

  subroutine free_PS_workarrays(iproc,igpu,keepzf,gpuPCGred,keepGPUmemory,w)
    use dynamic_memory
    integer, intent(in) :: keepzf,gpuPCGred,keepGPUmemory,igpu,iproc
    type(PS_workarrays), intent(inout) :: w
    call f_free_ptr(w%dlogeps)
    call f_free_ptr(w%oneoeps)
    call f_free_ptr(w%corr)
    call f_free_ptr(w%epsinnersccs)
    call f_free_ptr(w%pol_charge)
    call f_free_ptr(w%pot_old)
    call f_free_ptr(w%rho_old)
    call f_free_ptr(w%res_old)
    if(keepzf == 1) call f_free_ptr(w%zf)
    if (gpuPCGRed == 1) then
       if (keepGPUmemory == 1) then
          call cudafree(w%z_GPU)
          call cudafree(w%r_GPU)
          call cudafree(w%oneoeps_GPU)
          call cudafree(w%p_GPU)
          call cudafree(w%q_GPU)
          call cudafree(w%x_GPU)
          call cudafree(w%corr_GPU)
          call cudafree(w%alpha_GPU)
          call cudafree(w%beta_GPU)
          call cudafree(w%beta0_GPU)
          call cudafree(w%kappa_GPU)
       end if
    end if
    if (igpu == 1) then
       if (iproc == 0) then
          if (keepGPUmemory == 1) then
             call cudafree(w%work1_GPU)
             call cudafree(w%work2_GPU)
          endif
          call cudafree(w%k_GPU)
       endif
    end if

  end subroutine free_PS_workarrays

  !> Free memory used by the kernel operation
  !! @ingroup PSOLVER
  subroutine pkernel_free(kernel)
    use dynamic_memory
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    integer :: i_stat
    call f_free_ptr(kernel%kernel)
    call f_free_ptr(kernel%counts)
    call f_free_ptr(kernel%displs)
    call free_PS_workarrays(kernel%mpi_env%iproc,kernel%igpu,&
         kernel%keepzf,kernel%gpuPCGred,kernel%keepGPUmemory,kernel%w)
    !free GPU data
    if (kernel%igpu == 1) then
       if (kernel%mpi_env%iproc == 0) then
          call f_free_ptr(kernel%rhocounts)
          call f_free_ptr(kernel%rhodispls)
          if (kernel%initCufftPlan == 1) then
             call cufftDestroy(kernel%plan(1))
             call cufftDestroy(kernel%plan(2))
             call cufftDestroy(kernel%plan(3))
             call cufftDestroy(kernel%plan(4))
          endif
       endif
    end if

    !cannot yet free the communicators of the poisson kernel
    !lack of handling of the reference counter

  end subroutine pkernel_free


end module PStypes
