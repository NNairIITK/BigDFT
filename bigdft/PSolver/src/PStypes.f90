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
  use f_input_file, only: ATTRS
  implicit none

  private

  character(len=*), parameter :: KERNEL_VARIABLES        = 'kernel'
  character(len=*), parameter :: SCREENING               = 'screening' 
  character(len=*), parameter :: ISF_ORDER               = 'isf_order' 
  character(len=*), parameter :: STRESS_TENSOR           = 'stress_tensor' 
  character(len=*), parameter :: ENVIRONMENT_VARIABLES   = 'environment'
  character(len=*), parameter :: CAVITY_KEY              = 'cavity' 
  character(len=*), parameter :: EPSILON_KEY             = 'epsilon' 
  character(len=*), parameter :: EDENSMAXMIN             = 'edensmaxmin' 
  character(len=*), parameter :: DELTA_KEY               = 'delta' 
  character(len=*), parameter :: CAVITATION              = 'cavitation' 
  character(len=*), parameter :: GAMMAS_KEY              = 'gammaS' 
  character(len=*), parameter :: ALPHAS_KEY              = 'alphaS' 
  character(len=*), parameter :: BETAV_KEY               = 'betaV' 
  character(len=*), parameter :: GPS_ALGORITHM           = 'gps_algorithm'
  character(len=*), parameter :: PI_ETA                  = 'pi_eta' 
  character(len=*), parameter :: INPUT_GUESS             = 'input_guess' 
  character(len=*), parameter :: FD_ORDER                = 'fd_order' 
  character(len=*), parameter :: ITERMAX                 = 'itermax' 
  character(len=*), parameter :: MINRES                  = 'minres' 
  character(len=*), parameter, public :: SETUP_VARIABLES         = 'setup'
  character(len=*), parameter :: ACCEL                   = 'accel' 
  character(len=*), parameter :: KEEP_GPU_MEMORY         = 'keep_gpu_memory' 
  character(len=*), parameter :: TASKGROUP_SIZE_KEY      = 'taskgroup_size' 
  character(len=*), parameter :: GLOBAL_DATA             = 'global_data' 
  character(len=*), parameter, public :: VERBOSITY               = 'verbose' 
  character(len=*), parameter :: OUTPUT                  = 'output' 
  character(len=*), parameter :: DICT_COMPLETED          = '__dict_has_been_checked__'//ATTRS

  
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
     integer(f_address) :: alpha_GPU, beta_GPU, kappa_GPU, beta0_GPU
  end type PS_workarrays

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
     type(PSolver_options) :: opt !<Datatype controlling the operations of the solver
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

  public :: pkernel_null,PSolver_energies_null,pkernel_free,pkernel_allocate_cavity
  public :: pkernel_set_epsilon,PS_allocate_cavity_workarrays,build_cavity_from_rho
  public :: ps_allocate_lowlevel_workarrays,PSolver_options_null,PS_input_dict
  public :: release_PS_workarrays,PS_release_lowlevel_workarrays,PS_set_options,pkernel_init

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
  end subroutine nullify_work_arrays

  pure function pkernel_null() result(k)
    use environment, only : PS_VAC_ENUM
    implicit none
    type(coulomb_operator) :: k
    k%itype_scf=0
    k%geocode='F'
    call nullify_f_enum(k%method)
    k%cavity=cavity_default()
    k%opt=PSolver_options_null()
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
       end if
    end if
    if (igpu == 1) then
       if (iproc == 0) then
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

  !> Initialization of the Poisson kernel
  !! @ingroup PSOLVER
  function pkernel_init(iproc,nproc,dict,geocode,ndims,hgrids,angrad,mpi_env) result(kernel)
    use yaml_output
    use dictionaries
    use numerics
    use wrapper_MPI
    use f_enums
    use environment
    use f_input_file, only: input_file_dump
    implicit none
    integer, intent(in) :: iproc      !< Proc Id
    integer, intent(in) :: nproc      !< Number of processes
    type(dictionary), pointer :: dict !< dictionary of the input variables
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids
    real(gp), dimension(3), intent(in), optional :: angrad
    type(mpi_environment), intent(in), optional :: mpi_env
    type(coulomb_operator) :: kernel
    !local variables
    integer :: nthreads,group_size,taskgroup_size
    !$ integer :: omp_get_max_threads

    !nullification
    kernel=pkernel_null()

    !geocode and ISF family
    kernel%geocode=geocode
    !dimensions and grid spacings
    kernel%ndims=ndims
    kernel%hgrids=hgrids

    if (present(angrad)) then
       kernel%angrad=angrad
    else
       kernel%angrad=onehalf* [ pi, pi, pi ]
    end if

    !new treatment for the kernel input variables
    kernel%method=PS_VAC_ENUM
    if (DICT_COMPLETED .notin. dict) call PS_input_dict(dict) !complete the dictionary

    call PS_fill_variables(kernel,kernel%opt,dict) !fill the structure with basic results

    kernel%keepzf=1

    !import the mpi_environment if present
    if (present(mpi_env)) then
       call copy_mpi_environment(src=mpi_env,dest=kernel%mpi_env)
    else

       !specialized treatment
       taskgroup_size=dict//SETUP_VARIABLES//TASKGROUP_SIZE_KEY

       group_size=nproc
       !if the taskgroup size is not a divisor of nproc do not create taskgroups
       if (nproc >1 .and. taskgroup_size > 0 .and. taskgroup_size < nproc .and.&
            mod(nproc,taskgroup_size)==0) then
          group_size=taskgroup_size
       end if
       call mpi_environment_set(kernel%mpi_env,iproc,nproc,mpiworld(),group_size)
    end if

    !gpu can be used only for one nproc
    if (nproc > 1) kernel%igpu=0

    !-------------------
    nthreads=0
    if (kernel%mpi_env%iproc == 0 .and. kernel%mpi_env%igroup == 0 .and. kernel%opt%verbosity_level==1) then
       if (kernel%mu==0.0_gp) then 
          call yaml_comment('Kernel Initialization',hfill='-')
          call yaml_mapping_open('Poisson Kernel Initialization')
       else
          call yaml_mapping_open('Helmholtz Kernel Initialization')
          call yaml_map('Screening Length (AU)',1.0_gp/kernel%mu,fmt='(g25.17)')
       end if
       !we might also perform an input_file dump if needed
       call input_file_dump(dict)
       !$ nthreads = omp_get_max_threads()
       call yaml_map('MPI tasks',kernel%mpi_env%nproc)
       if (nthreads /=0) call yaml_map('OpenMP threads per MPI task',nthreads)
       if (kernel%igpu==1) call yaml_map('Kernel copied on GPU',.true.)
       if (kernel%method /= 'VAC') call yaml_map('Iterative method for Generalised Equation',str(kernel%method))
       if (kernel%method .hasattr. PS_RIGID_ENUM) call yaml_map('Cavity determination','rigid')
       if (kernel%method .hasattr. PS_SCCS_ENUM) call yaml_map('Cavity determination','sccs')
       call yaml_mapping_close() !kernel
    end if

  end function pkernel_init

  !>modifies the options of the poisson solver to switch certain options
  subroutine PS_set_options(kernel,global_data,calculate_strten,verbose,&
       update_cavity,use_input_guess,cavity_info,cavitation_terms,&
       potential_integral)
    implicit none
    type(coulomb_operator), intent(inout) :: kernel
    logical, intent(in), optional :: global_data,calculate_strten,verbose
    logical, intent(in), optional :: update_cavity,use_input_guess
    logical, intent(in), optional :: cavity_info,cavitation_terms
    real(gp), intent(in), optional :: potential_integral

    if (present(global_data     )) then
       if (global_data) then
          kernel%opt%datacode='G'
       else
          kernel%opt%datacode='D'
       end if
    end if
    if (present(calculate_strten)) kernel%opt%calculate_strten=calculate_strten
    if (present(verbose         )) then
       if (verbose) then
          kernel%opt%verbosity_level=1
       else
          kernel%opt%verbosity_level=0
       end if
    end if
    if (present(update_cavity   )) kernel%opt%update_cavity=update_cavity
    if (present(use_input_guess )) kernel%opt%use_input_guess=use_input_guess
    if (present(cavity_info     )) kernel%opt%cavity_info=cavity_info
    if (present(cavitation_terms)) kernel%opt%only_electrostatic=.not. cavitation_terms
    if (present(potential_integral)) kernel%opt%potential_integral =potential_integral

  end subroutine PS_set_options


  !>routine to fill the input variables of the kernel
  subroutine PS_input_dict(dict,dict_minimal)
    use dictionaries
    use f_input_file
    use yaml_parse
    implicit none
    !>input dictionary, a copy of the user input, to be filled
    !!with all the variables on exit
    type(dictionary), pointer :: dict
    type(dictionary), pointer, optional :: dict_minimal
    !local variables
    integer(f_integer) :: params_size
    !integer(kind = 8) :: cbuf_add !< address of c buffer
    character, dimension(:), allocatable :: params
    type(dictionary), pointer :: parameters
    type(dictionary), pointer :: parsed_parameters
    type(dictionary), pointer :: profiles
    type(dictionary), pointer :: nested,asis

    call f_routine(id='PS_input_dict')

    nullify(parameters,parsed_parameters,profiles)

    !alternative filling of parameters from hard-coded source file
    !call getstaticinputdef(cbuf_add,params_size)
    call getpsinputdefsize(params_size)
    !allocate array
    params=f_malloc_str(1,params_size,id='params')
    !fill it and parse dictionary
    call getpsinputdef(params)

    call yaml_parse_from_char_array(parsed_parameters,params)
    !there is only one document in the input variables specifications
    parameters=>parsed_parameters//0
    profiles => parsed_parameters//1
    call f_free_str(1,params)

    call input_file_complete(parameters,dict,imports=profiles)

    if (present(dict_minimal)) then
       nullify(nested,asis)
       call input_file_minimal(parameters,dict,dict_minimal,nested,asis)
    end if

    if (associated(parsed_parameters)) then
       call dict_free(parsed_parameters)
       nullify(parameters)
       nullify(profiles)
    else
       call dict_free(parameters)
    end if

    !write in the dictionary that it has been completed
    call set(dict//DICT_COMPLETED,.true.)

    call f_release_routine()

  end subroutine PS_input_dict

  subroutine PS_fill_variables(k,opt,dict)
    use dictionaries
    implicit none
    type(coulomb_operator), intent(inout) :: k
    type(PSolver_options), intent(inout) :: opt
    type(dictionary), pointer :: dict
    !local variables
    type(dictionary), pointer :: lvl,var

    ! Transfer dict values into input_variables structure.
    lvl => dict_iter(dict)
    do while(associated(lvl))
       var => dict_iter(lvl)
       do while(associated(var))
          call PS_input_fill(k,opt,dict_key(lvl),var)
          var => dict_next(var)
       end do
       lvl => dict_next(lvl)
    end do

  end subroutine PS_fill_variables

  !> Set the dictionary from the input variables
  subroutine PS_input_fill(k,opt, level, val)
    use PSbase
    use environment
    use yaml_output, only: yaml_warning
    use dictionaries
    use numerics
    implicit none
    type(coulomb_operator), intent(inout) :: k
    type(PSolver_options), intent(inout) :: opt
    type(dictionary), pointer :: val
    character(len = *), intent(in) :: level
    !local variables
    logical :: dummy_l
    real(gp) :: dummy_d
    integer, dimension(2) :: dummy_int !<to use as filling for input variables
    real(gp), dimension(2) :: dummy_gp !< to fill the input variables
    logical, dimension(2) :: dummy_log !< to fill the input variables
    character(len=256) :: dummy_char
    character(len = max_field_length) :: strn
    integer :: i, ipos

    if (index(dict_key(val), "_attributes") > 0) return

    select case(trim(level))
    case(KERNEL_VARIABLES)
       select case (trim(dict_key(val)))
       case(SCREENING)
          k%mu=val
       case(ISF_ORDER)
          k%itype_scf=val
       case(STRESS_TENSOR)
          opt%calculate_strten=val
       case DEFAULT
          if (k%mpi_env%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (ENVIRONMENT_VARIABLES)
       select case (trim(dict_key(val)))
       case (CAVITY_KEY)
          strn=val

          select case(trim(strn))
          case('vacuum')
             call f_enum_attr(k%method,PS_NONE_ENUM)
          case('rigid')
             call f_enum_attr(k%method,PS_RIGID_ENUM)
          case('sccs')   
             call f_enum_attr(k%method,PS_SCCS_ENUM)
          end select
       case (EPSILON_KEY)
          k%cavity%epsilon0=val
       case (EDENSMAXMIN)
          dummy_gp=val
          k%cavity%edensmin=dummy_gp(1)
          k%cavity%edensmax=dummy_gp(2)
       case (DELTA_KEY)
          dummy_d=val
          ! Divided by 4 because both rigid cavities are 4*delta spread 
          k%cavity%delta=0.25_gp*dummy_d
       case (CAVITATION)
          dummy_l=val
          opt%only_electrostatic=.not. dummy_l
       case (GAMMAS_KEY)
          dummy_d=val
          k%cavity%gammaS=dummy_d*SurfAU
       case (ALPHAS_KEY)
          dummy_d=val
          k%cavity%alphaS=dummy_d*SurfAU
       case (BETAV_KEY)
          dummy_d=val
          k%cavity%betaV=dummy_d/AU_GPa
       case (GPS_ALGORITHM)
          strn=val
          select case(trim(strn))
          case('PI')
             call f_enum_update(dest=k%method,src=PS_PI_ENUM)
          case('PCG')
             call f_enum_update(dest=k%method,src=PS_PCG_ENUM)
          end select
       case (PI_ETA)
          k%PI_eta=val
       case (INPUT_GUESS)
          opt%use_input_guess=val
       case (FD_ORDER)
          k%nord=val
       case (ITERMAX)
          k%max_iter=val
       case (MINRES)
          k%minres=val
       case DEFAULT
          if (k%mpi_env%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (SETUP_VARIABLES)
       select case (trim(dict_key(val)))       
       case (ACCEL)
          dummy_l=val
          if (dummy_l) then
             k%igpu=1
          else
             k%igpu=0
          end if
       case (KEEP_GPU_MEMORY)
          dummy_l=val
          if (dummy_l) then
             k%keepGPUmemory=1
          else
             k%keepGPUmemory=0
          end if
       case (TASKGROUP_SIZE_KEY)

       case (GLOBAL_DATA)
          dummy_l=val
          if (dummy_l) then
             opt%datacode='G'
          else
             opt%datacode='D'
          end if
       case (VERBOSITY)
          dummy_l=val
          if (dummy_l) then
             opt%verbosity_level=1
          else
             opt%verbosity_level=0
          end if
       case (OUTPUT)
          !for the moment no treatment, to be added
       case DEFAULT
          if (k%mpi_env%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case DEFAULT
    end select
  END SUBROUTINE PS_input_fill


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