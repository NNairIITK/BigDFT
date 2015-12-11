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
  use dynamic_memory
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
     !> dielectric function epsilon, continuous and differentiable in the whole
     !domain, to be used in the case of Preconditioned Conjugate Gradient method
     real(dp), dimension(:,:), pointer :: eps
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
     !>work array needed to store the zero-padded part of the 
     !!density and the potential. Cannot be used as-is, it must be
     !!copied back into a distributed array with the 
     !!finalize_hartree_results routine.
     real(dp), dimension(:,:,:), pointer :: zf
     !> input guess vectors to be preserved for future use
     !!or work arrays, might be of variable dimension
     !!(either full of distributed)
     real(dp), dimension(:,:), pointer :: pot
     real(dp), dimension(:,:), pointer :: rho
     !> Polarization charge vector for print purpose only.
     real(dp), dimension(:,:), pointer :: rho_pol
     !> arrays for the execution of the PCG algorithm
     real(dp), dimension(:,:), pointer :: res,z,p,q

     integer(f_address) :: work1_GPU,work2_GPU,rho_GPU,pot_ion_GPU,k_GPU !<addresses for the GPU memory 
     integer(f_address) :: p_GPU,q_GPU,r_GPU,x_GPU,z_GPU,oneoeps_GPU,corr_GPU!<addresses for the GPU memory 
     !> GPU scalars. Event if they are scalars of course their address is needed
     integer(f_address) :: alpha_GPU, beta_GPU, kappa_GPU, beta0_GPU, eexctX_GPU, reduc_GPU, ehart_GPU
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
     integer :: stay_on_gpu
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
     !> if .true., and the cavity is set to 'sccs' attribute, then the epsilon is updated according to rhov
     logical :: update_cavity
     !> if .true. calculate the stress tensor components.
     !! The stress tensor calculation are so far validated only for orthorhombic cells and vacuum treatments
     logical :: calculate_strten
     !> Use the input guess procedure for the solution in the case of a electrostatic medium
     !! this value is ignored in the case of a vacuum solver
     logical :: use_input_guess
     !> Trigger the calculation of the electrostatic contribution only
     !! if .true., the code only calculates the electrostatic contribution
     !! and the cavitation terms are neglected
     !> extract the polarization charge and the dielectric function, to be used for plotting purposes
     logical :: cavity_info
     logical :: only_electrostatic
     !> Total integral on the supercell of the final potential on output
     !! clearly meaningful only for Fully periodic BC, ignored in the other cases
     real(gp) :: potential_integral
  end type PSolver_options

  public :: pkernel_null,PSolver_energies_null,pkernel_free,pkernel_allocate_cavity
  public :: pkernel_set_epsilon,PS_allocate_cavity_workarrays,build_cavity_from_rho
  public :: ps_allocate_lowlevel_workarrays,PSolver_options_null
  public :: release_PS_workarrays,PS_release_lowlevel_workarrays

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

  pure function PSolver_options_null() result(o)
    implicit none
    type(PSolver_options) :: o

    o%datacode           ='X'
    o%verbosity_level    =0
    o%update_cavity      =.false.
    o%calculate_strten   =.false.
    o%use_input_guess    =.false.
    o%cavity_info        =.false.
    o%only_electrostatic =.true.
    o%potential_integral =0.0_gp
   end function PSolver_options_null


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
    nullify(w%rho_pol)
    nullify(w%rho)
    nullify(w%pot)
    nullify(w%res)
    nullify(w%zf)
    nullify(w%z)
    nullify(w%p)
    nullify(w%q)
    nullify(w%eps)
    call f_zero(w%work1_GPU)
    call f_zero(w%work2_GPU)
    call f_zero(w%rho_GPU)
    call f_zero(w%pot_ion_GPU)
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
    call f_zero(w%eexctX_GPU)
    call f_zero(w%ehart_GPU)
    call f_zero(w%reduc_GPU)
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
    call nullify_mpi_environment(k%mpi_env)
    call nullify_mpi_environment(k%inplane_mpi)
    call nullify_mpi_environment(k%part_mpi)
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
    integer, intent(in) :: keepzf,gpuPCGred,keepGPUmemory,igpu,iproc
    type(PS_workarrays), intent(inout) :: w
    call f_free_ptr(w%eps)
    call f_free_ptr(w%dlogeps)
    call f_free_ptr(w%oneoeps)
    call f_free_ptr(w%corr)
    call f_free_ptr(w%epsinnersccs)
    call f_free_ptr(w%rho_pol)
    call f_free_ptr(w%pot)
    call f_free_ptr(w%rho)
    call f_free_ptr(w%res)
    call f_free_ptr(w%z)
    call f_free_ptr(w%p)
    call f_free_ptr(w%q)
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
          call cudafree(w%ehart_GPU)
          call cudafree(w%eexctX_GPU)
          call cudafree(w%reduc_GPU)
       end if
    end if
    if (igpu == 1) then
       if (iproc == 0) then
	  call cudadestroystream()
          if (keepGPUmemory == 1) then
             call cudafree(w%work1_GPU)
             call cudafree(w%work2_GPU)
             call cudafree(w%rho_GPU)
             call cudafree(w%pot_ion_GPU)
          endif
          call cudafree(w%k_GPU)
       endif
    end if

  end subroutine free_PS_workarrays


  subroutine release_PS_workarrays(keepzf,w,use_input_guess)
    implicit none
    integer, intent(in) :: keepzf
    type(PS_workarrays), intent(inout) :: w
    logical, intent(in) :: use_input_guess
    if(keepzf /= 1) call f_free_ptr(w%zf)
!    call f_free_ptr(w%pot)
    if (.not. use_input_guess) call f_free_ptr(w%pot)
  end subroutine release_PS_workarrays

  !> Free memory used by the kernel operation
  !! @ingroup PSOLVER
  subroutine pkernel_free(kernel)
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
    call release_mpi_environment(kernel%inplane_mpi)
    call release_mpi_environment(kernel%part_mpi)
    call release_mpi_environment(kernel%mpi_env)
  end subroutine pkernel_free

  !> allocate the workarrays needed to perform the 
  !! GPS operation
  subroutine PS_allocate_lowlevel_workarrays(cudasolver,use_input_guess,rho,kernel)
    use f_utils, only: f_zero
    use wrapper_linalg, only: axpy
    implicit none
    logical, intent(in) :: cudasolver,use_input_guess
    type(coulomb_operator), intent(inout) :: kernel
    real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(in) :: rho !< initial rho, needed for PCG
    !local variables
    integer :: n1,n23

    !we need to reallocate the zf array with the right size when called with stress_tensor and gpu
    if(kernel%keepzf == 1) then
       if(kernel%igpu==1 .and. .not. cudasolver) then !LG: what means that?
          call f_free_ptr(kernel%w%zf)
          kernel%w%zf = f_malloc_ptr([kernel%grid%md1, kernel%grid%md3, &
               2*kernel%grid%md2/kernel%mpi_env%nproc],id='zf')
       end if
    else
       kernel%w%zf = f_malloc_ptr([kernel%grid%md1, kernel%grid%md3, &
            2*kernel%grid%md2/kernel%mpi_env%nproc],id='zf')
    end if
    n23=kernel%grid%m3*kernel%grid%n3p
    n1=kernel%grid%m1

    select case(trim(str(kernel%method)))
    case('PCG')
       if (use_input_guess .and. &
            all([associated(kernel%w%res),associated(kernel%w%pot)])) then
          !call axpy(n1*n23,1.0_gp,rho(1,1),1,kernel%w%res(1,1),1)
          call f_memcpy(src=rho,dest=kernel%w%res)
       else
          !allocate if it is the first time
          if (associated(kernel%w%pot)) then
             call f_zero(kernel%w%pot)
          else
             kernel%w%pot=f_malloc0_ptr([n1,n23],id='pot')
          end if
          if (.not. associated(kernel%w%res)) &
               kernel%w%res=f_malloc_ptr([n1,n23],id='res')
          call f_memcpy(src=rho,dest=kernel%w%res)
       end if
       kernel%w%q=f_malloc0_ptr([n1,n23],id='q')
       kernel%w%p=f_malloc0_ptr([n1,n23],id='p')
       kernel%w%z=f_malloc_ptr([n1,n23],id='z')
       kernel%w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol') !>>>>>>>>>>here the switch
    case('PI')
       if (use_input_guess .and. &
            associated(kernel%w%pot)) then
       else
          !allocate if it is the first time
          if (associated(kernel%w%pot)) then
             call f_zero(kernel%w%pot)
          else
       kernel%w%pot=f_malloc_ptr([kernel%ndims(1),kernel%ndims(2)*kernel%ndims(3)],id='pot')
          end if
       end if

       !kernel%w%pot=f_malloc_ptr([kernel%ndims(1),kernel%ndims(2)*kernel%ndims(3)],id='pot')
       kernel%w%rho=f_malloc0_ptr([kernel%ndims(1),kernel%ndims(2)*kernel%ndims(3)],id='rho')
       kernel%w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol') !>>>>>>>>>>here the switch
    end select

  end subroutine PS_allocate_lowlevel_workarrays

  !> this is useful to deallocate useless space and to 
  !! also perform extra treatment for the inputguess
  subroutine PS_release_lowlevel_workarrays(cavity_info,use_input_guess,rho,kernel)
    use wrapper_linalg, only: axpy
    implicit none
    logical, intent(in) :: cavity_info,use_input_guess
    type(coulomb_operator), intent(inout) :: kernel
    real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(in) :: rho !< initial rho, needed for PCG

    select case(trim(str(kernel%method)))
    case('PCG')
!       if (use_input_guess) then
!          !preserve the previous values for the input guess
!          call axpy(size(kernel%w%res),-1.0_gp,rho(1,1),1,kernel%w%res(1,1),1)
!       else
!          call f_free_ptr(kernel%w%res)
!       end if
       if (.not. use_input_guess) call f_free_ptr(kernel%w%res)
       call f_free_ptr(kernel%w%q)
       call f_free_ptr(kernel%w%p)
       call f_free_ptr(kernel%w%z)
       call f_free_ptr(kernel%w%rho_pol) !>>>>>>>>>>here the switch
    case('PI')
       call f_free_ptr(kernel%w%rho)
       call f_free_ptr(kernel%w%rho_pol) !>>>>>>>>>>here the switch
    end select

  end subroutine PS_release_lowlevel_workarrays

  !> allocate the workarrays for the first initialization
  !! their allocation depends on the treatment which we are going to
  !! apply
  subroutine PS_allocate_cavity_workarrays(n1,n23,ndims,method,w)
    use environment, only: PS_SCCS_ENUM
    use dynamic_memory
    implicit none
    integer, intent(in) :: n1,n23
    integer, dimension(3), intent(in) :: ndims
    type(f_enumerator), intent(in) :: method
    type(PS_workarrays), intent(inout) :: w

    select case(trim(str(method)))
    case('PCG')
       !w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol') !>>>>>>>>>>here the switch
       w%eps=f_malloc_ptr([n1,n23],id='eps')
       w%corr=f_malloc_ptr([n1,n23],id='corr')
       w%oneoeps=f_malloc_ptr([n1,n23],id='oneosqrteps')
       w%epsinnersccs=f_malloc_ptr([n1,n23],id='epsinnersccs')
    case('PI')
       !w%rho_pol=f_malloc_ptr([n1,n23],id='rho_pol') !>>>>>>>>>>here the switch
       w%dlogeps=f_malloc_ptr([3,ndims(1),ndims(2),ndims(3)],id='dlogeps')
       w%oneoeps=f_malloc_ptr([n1,n23],id='oneoeps')
       w%epsinnersccs=f_malloc_ptr([n1,n23],id='epsinnersccs')
       !w%epsinnersccs=f_malloc_ptr([n1,ndims(2)*ndims(3)],id='epsinnersccs')
    end select
   end subroutine PS_allocate_cavity_workarrays

  !> create the memory space needed to store the arrays for the 
  !! description of the cavity
  subroutine pkernel_allocate_cavity(kernel,vacuum)
    use environment, only: PS_SCCS_ENUM,vacuum_eps
    use f_utils, only: f_zero
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    logical, intent(in), optional :: vacuum !<if .true. the cavity is allocated as no cavity exists, i.e. only vacuum
    !local variables
    integer :: n1,n23,i1,i23

    n1=kernel%ndims(1)
    n23=kernel%ndims(2)*kernel%grid%n3p
!!$    call PS_allocate_cavity_workarrays(n1,n23,kernel%ndims,&
!!$         kernel%method,kernel%w)
    if (present(vacuum)) then
       if (vacuum) then
          select case(trim(str(kernel%method)))
          case('PCG')
             call f_zero(kernel%w%corr)
             do i23=1,n23
                do i1=1,n1
                   kernel%w%eps(i1,i23)=vacuum_eps
                   kernel%w%oneoeps(i1,i23)=1.0_dp/sqrt(vacuum_eps)
                end do
             end do
          case('PI')
             call f_zero(kernel%w%dlogeps)
             do i23=1,n23
                do i1=1,n1
                   kernel%w%oneoeps(i1,i23)=1.0_dp/vacuum_eps
                end do
             end do
          end select
          if (kernel%method .hasattr. PS_SCCS_ENUM) &
               call f_zero(kernel%w%epsinnersccs)
       end if
    end if

  end subroutine pkernel_allocate_cavity

  !> set the epsilon in the pkernel structure as a function of the seteps variable.
  !! This routine has to be called each time the dielectric function has to be set
  subroutine pkernel_set_epsilon(kernel,eps,dlogeps,oneoeps,oneosqrteps,corr)
    use yaml_strings
    use dynamic_memory
    use FDder
    use dictionaries, only: f_err_throw
    implicit none
    !> Poisson Solver kernel
    type(coulomb_operator), intent(inout) :: kernel
    !> dielectric function. Needed for non VAC methods, given in full dimensions
    real(dp), dimension(:,:,:), intent(in), optional :: eps
    !> logarithmic derivative of epsilon. Needed for PCG method.
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:,:), intent(in), optional :: dlogeps
    !> inverse of epsilon. Needed for PI method.
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:), intent(in), optional :: oneoeps
    !> inverse square root of epsilon. Needed for PCG method.
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:), intent(in), optional :: oneosqrteps
    !> correction term of the Generalized Laplacian
    !! if absent, it will be calculated from the array of epsilon
    real(dp), dimension(:,:,:), intent(in), optional :: corr
    !local variables
    integer :: n1,n23,i3s,i23,i3,i2,i1
    real(kind=8) :: pi
    real(dp), dimension(:,:,:), allocatable :: de2,ddeps
    real(dp), dimension(:,:,:,:), allocatable :: deps

    pi=4.0_dp*atan(1.0_dp)
    if (present(corr)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(corr) /= kernel%ndims)) &
            call f_err_throw('Error in the dimensions of the array corr,'//&
            trim(yaml_toa(shape(corr))))
    end if
    if (present(eps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(eps) /= kernel%ndims)) &
            call f_err_throw('Error in the dimensions of the array epsilon,'//&
            trim(yaml_toa(shape(eps))))
    end if
    if (present(oneoeps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(oneoeps) /= kernel%ndims)) &
            call f_err_throw('Error in the dimensions of the array oneoeps,'//&
            trim(yaml_toa(shape(oneoeps))))
    end if
    if (present(oneosqrteps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(oneosqrteps) /= kernel%ndims)) &
            call f_err_throw('Error in the dimensions of the array oneosqrteps,'//&
            trim(yaml_toa(shape(oneosqrteps))))
    end if
    if (present(dlogeps)) then
       !check the dimensions (for the moment no parallelism)
       if (any(shape(dlogeps) /= &
            [3,kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)])) &
            call f_err_throw('Error in the dimensions of the array dlogeps,'//&
            trim(yaml_toa(shape(dlogeps))))
    end if

    !store the arrays needed for the method
    !the stored arrays are of rank two to collapse indices for
    !omp parallelism
    n1=kernel%ndims(1)
    n23=kernel%ndims(2)*kernel%grid%n3p
    !starting point in third direction
    i3s=kernel%grid%istart+1
    if (kernel%grid%n3p==0) i3s=1
    select case(trim(str(kernel%method)))
    case('PCG')
       !check the dimensions of the associated arrays
       if (all([associated(kernel%w%corr),associated(kernel%w%oneoeps)])) then
          !then check the shapes
          if (any(shape(kernel%w%oneoeps) /= [n1,n23])) &
               call f_err_throw('Incorrect shape of oneoeps')
          if (any(shape(kernel%w%corr) /= [n1,n23])) &
               call f_err_throw('Incorrect shape of corr')
          if (present(corr)) then
             call f_memcpy(n=n1*n23,src=corr(1,1,i3s),dest=kernel%w%corr)
          else if (present(eps)) then
!!$        !allocate work arrays
             deps=f_malloc([kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),3],id='deps')
             de2 =f_malloc(kernel%ndims,id='de2')
             ddeps=f_malloc(kernel%ndims,id='ddeps')

             call nabla_u_and_square(kernel%geocode,kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
                  eps,deps,de2,kernel%nord,kernel%hgrids)

             call div_u_i(kernel%geocode,kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
                  deps,ddeps,kernel%nord,kernel%hgrids)
             i23=1
             do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
                do i2=1,kernel%ndims(2)
                   do i1=1,kernel%ndims(1)
                      kernel%w%corr(i1,i23)=(-0.125d0/pi)*&
                           (0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
                   end do
                   i23=i23+1
                end do
             end do
             call f_free(deps)
             call f_free(ddeps)
             call f_free(de2)
          else
             call f_err_throw('For method "PCG" the arrays corr or epsilon should be present')   
          end if
          if (present(oneosqrteps)) then
             call f_memcpy(n=n1*n23,src=oneosqrteps(1,1,i3s),&
                  dest=kernel%w%oneoeps)
          else if (present(eps)) then
             i23=1
             do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
                do i2=1,kernel%ndims(2)
                   do i1=1,kernel%ndims(1)
                      kernel%w%oneoeps(i1,i23)=1.0_dp/sqrt(eps(i1,i2,i3))
                   end do
                   i23=i23+1
                end do
             end do
          else
             call f_err_throw('For method "PCG" the arrays oneosqrteps or epsilon should be present')
          end if
          if (present(eps)) then
             call f_memcpy(n=n1*n23,src=eps(1,1,i3s),&
                  dest=kernel%w%eps)
          else
             call f_err_throw('For method "PCG" the arrays eps should be present')
          end if
       else
          call f_err_throw('For method "PCG" the arrays oneosqrteps'//&
               ' and corr have to be associated, call PS_allocate_cavity_workarrays')
       end if
    case('PI')
       if (all([associated(kernel%w%dlogeps),associated(kernel%w%oneoeps)])) then
          !then check the shapes
          if (any(shape(kernel%w%oneoeps) /= [n1,n23])) &
               call f_err_throw('Incorrect shape of oneoeps')
          if (any(shape(kernel%w%dlogeps) /= [3,kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)])) &
               call f_err_throw('Incorrect shape of dlogeps')
          if (present(dlogeps)) then
             call f_memcpy(src=dlogeps,dest=kernel%w%dlogeps)
          else if (present(eps)) then
             !allocate arrays
             deps=f_malloc([kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),3],id='deps')
             call nabla_u(kernel%geocode,kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
                  eps,deps,kernel%nord,kernel%hgrids)
             do i3=1,kernel%ndims(3)
                do i2=1,kernel%ndims(2)
                   do i1=1,kernel%ndims(1)
                      !switch and create the logarithmic derivative of epsilon
                      kernel%w%dlogeps(1,i1,i2,i3)=deps(i1,i2,i3,1)/eps(i1,i2,i3)
                      kernel%w%dlogeps(2,i1,i2,i3)=deps(i1,i2,i3,2)/eps(i1,i2,i3)
                      kernel%w%dlogeps(3,i1,i2,i3)=deps(i1,i2,i3,3)/eps(i1,i2,i3)
                   end do
                end do
             end do
             call f_free(deps)
          else
             call f_err_throw('For method "PI" the arrays dlogeps or epsilon should be present')
          end if
          if (present(oneoeps)) then
             call f_memcpy(n=n1*n23,src=oneoeps(1,1,i3s),&
                  dest=kernel%w%oneoeps)
          else if (present(eps)) then
             i23=1
             do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
                do i2=1,kernel%ndims(2)
                   do i1=1,kernel%ndims(1)
                      kernel%w%oneoeps(i1,i23)=1.0_dp/eps(i1,i2,i3)
                   end do
                   i23=i23+1
                end do
             end do
          else
             call f_err_throw('For method "PI" the arrays oneoeps or epsilon should be present')
          end if
       else
          call f_err_throw('For method "PI" the arrays oneoeps '//&
               'and dlogeps have to be associated, call PS_allocate_cavity_workarrays')
       end if
    end select

  end subroutine pkernel_set_epsilon

  subroutine build_cavity_from_rho(rho,nabla2_rho,delta_rho,cc_rho,kernel,&
       depsdrho,dsurfdrho,IntSur,IntVol)
    use environment
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(in) :: rho,nabla2_rho,delta_rho,cc_rho
    !> functional derivative of the sc epsilon with respect to 
    !! the electronic density, in distributed memory
    real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(out) :: depsdrho
    !> functional derivative of the surface integral with respect to 
    !! the electronic density, in distributed memory
    real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(out) :: dsurfdrho
    real(dp), intent(out) :: IntSur,IntVol
    !local variables
    real(dp), parameter :: innervalue = 0.9d0 !to be defined differently
    integer :: n01,n02,n03,i3s,i1,i2,i3,i23
    real(dp) :: rh,d2,d,dd,de,epsm1

    IntSur=0.d0
    IntVol=0.d0

    n01=kernel%ndims(1)
    n02=kernel%ndims(2)
    n03=kernel%ndims(3)
    !starting point in third direction
    i3s=kernel%grid%istart+1
    epsm1=(kernel%cavity%epsilon0-1.0_gp)
    !now fill the pkernel arrays according the the chosen method
    select case(trim(str(kernel%method)))
    case('PCG')
       !in PCG we only need corr, oneosqrtepsilon
       i23=1
       do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
          !do i3=1,n03
          do i2=1,n02
             do i1=1,n01
                if (kernel%w%epsinnersccs(i1,i23).gt.innervalue) then 
                   kernel%w%eps(i1,i23)=1.d0 !eps(i1,i2,i3)
                   kernel%w%oneoeps(i1,i23)=1.d0 !oneosqrteps(i1,i2,i3)
                   kernel%w%corr(i1,i23)=0.d0 !corr(i1,i2,i3)
                   depsdrho(i1,i23)=0.d0
                   dsurfdrho(i1,i23)=0.d0
                else
                   rh=rho(i1,i23)
                   d2=nabla2_rho(i1,i23)
                   d=sqrt(d2)
                   dd = delta_rho(i1,i23)
                   de=epsprime(rh,kernel%cavity)
                   depsdrho(i1,i23)=de
                   kernel%w%eps(i1,i23)=eps(rh,kernel%cavity)
                   kernel%w%oneoeps(i1,i23)=oneosqrteps(rh,kernel%cavity)
                   kernel%w%corr(i1,i23)=corr_term(rh,d2,dd,kernel%cavity)
                   dsurfdrho(i1,i23)=surf_term(rh,d2,dd,cc_rho(i1,i23),kernel%cavity)/epsm1
                   !evaluate surfaces and volume integrals
                   IntSur=IntSur + de*d/epsm1
                   IntVol=IntVol + (kernel%cavity%epsilon0-eps(rh,kernel%cavity))/epsm1
                end if
             end do
             i23=i23+1
          end do
       end do
    case('PI')
       !for PI we need  dlogeps,oneoeps
       !first oneovereps
       i23=1
       do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
          do i2=1,n02
             do i1=1,n01
                if (kernel%w%epsinnersccs(i1,i23).gt.innervalue) then ! Check for inner sccs cavity value to fix as vacuum
                   kernel%w%oneoeps(i1,i23)=1.d0 !oneoeps(i1,i2,i3)
                   depsdrho(i1,i23)=0.d0
                   dsurfdrho(i1,i23)=0.d0
                else
                   rh=rho(i1,i23)
                   d2=nabla2_rho(i1,i23)
                   d=sqrt(d2)
                   dd = delta_rho(i1,i23)
                   de=epsprime(rh,kernel%cavity)
                   depsdrho(i1,i23)=de
                   kernel%w%oneoeps(i1,i23)=oneoeps(rh,kernel%cavity) 
                   dsurfdrho(i1,i23)=surf_term(rh,d2,dd,cc_rho(i1,i23),kernel%cavity)/epsm1

                   !evaluate surfaces and volume integrals
                   IntSur=IntSur + de*d/epsm1
                   IntVol=IntVol + (kernel%cavity%epsilon0-eps(rh,kernel%cavity))/epsm1
                end if
             end do
             i23=i23+1
          end do
       end do
    end select

    IntSur=IntSur*product(kernel%hgrids)
    IntVol=IntVol*product(kernel%hgrids)


  end subroutine build_cavity_from_rho



end module PStypes
