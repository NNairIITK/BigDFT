!> @file
!!    Modulefile for the definition of the basic structures
!!
!! @author
!!    B. Schaefer, L. Genovese
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module SPREDtypes
  use f_enums
  use wrapper_MPI
  use SPREDbase
  use dynamic_memory
  use f_input_file, only: ATTRS
  implicit none

  private

  character(len=*), parameter :: DICT_COMPLETED          = '__dict_has_been_checked__'//ATTRS

  
  !>Datatype defining the inputs variables for SPRED
  type, public :: SPRED_inputs
    !> spred input variables
    type(f_enumerator) :: fp_method
    integer :: natx_sphere       !< number of atoms in each sphere (for periodic fingerprint)
    integer :: angmom       !< angular momentum of gaussian orbitals for overlap matrix fingerprints (both periodic and free BC)
  end type SPRED_inputs


  public :: SPRED_input_dict

contains
  subroutine SPRED_init(inputs)
    implicit none
  end subroutine SPRED_init 

  !>routine to fill the input variables of the kernel
  subroutine SPRED_input_dict(dict,dict_minimal)
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

  end subroutine SPRED_input_dict

  subroutine PS_fill_variables(k,inputs,dict)
    use dictionaries
    implicit none
    type(coulomb_operator), intent(inout) :: k
    type(SPRED_inputs), intent(inout) :: inputs
    type(dictionary), pointer :: dict
    !local variables
    type(dictionary), pointer :: lvl,var

    ! Transfer dict values into input_variables structure.
    lvl => dict_iter(dict)
    do while(associated(lvl))
       var => dict_iter(lvl)
       do while(associated(var))
          call PS_input_fill(k,inputs,dict_key(lvl),var)
          var => dict_next(var)
       end do
       lvl => dict_next(lvl)
    end do

  end subroutine PS_fill_variables

  !> Set the dictionary from the input variables
  subroutine PS_input_fill(k,inputs, level, val)
    use PSbase
    use environment
    use yaml_output, only: yaml_warning
    use dictionaries
    use numerics
    implicit none
    type(coulomb_operator), intent(inout) :: k
    type(SPRED_inputs), intent(inout) :: inputs
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
          inputs%calculate_strten=val
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
          inputs%only_electrostatic=.not. dummy_l
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
          inputs%use_input_guess=val
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
          strn=val
          select case(trim(strn))
          case('CUDA')
             k%igpu=1
          case('none')
             k%igpu=0
          end select
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
             inputs%datacode='G'
          else
             inputs%datacode='D'
          end if
       case (VERBOSITY)
          dummy_l=val
          if (dummy_l) then
             inputs%verbosity_level=1
          else
             inputs%verbosity_level=0
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


end module SPREDtypes
