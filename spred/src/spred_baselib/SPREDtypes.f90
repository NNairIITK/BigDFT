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

  !> Datatype defining the inputs variables for SPRED
  type, public :: SPRED_inputs
    !> fingerprint input variables
    type(f_enumerator) :: fp_method
    integer :: fp_natx_sphere       !< number of atoms in each sphere (for periodic fingerprint)
    integer :: fp_angmom       !< angular momentum of gaussian orbitals for overlap matrix fingerprints (both periodic and free BC)


    !>global input variables
    integer :: glbl_nwrite
  end type SPRED_inputs

  public :: SPRED_read_uinp


contains

  !> read user input
  subroutine SPRED_read_uinp(radical,inputs,mpi_env_)
    use SPREDbase
    use dictionaries
    use yaml_output
    use f_input_file, only: input_file_dump
    implicit none
    !parameter
    character(len = *), intent(in) :: radical    !< The name of the run. use "input" if empty
    type(SPRED_inputs), intent(out) :: inputs
    type(mpi_environment), optional, intent(in) :: mpi_env_ !< Environment of the reading. Used for broadcasting the result
    !internal
    logical :: exists_user
    type(mpi_environment) :: mpi_env !< Environment of the reading. Used for broadcasting the result
    character(len = max_field_length) :: fname
    type(dictionary), pointer :: dict

    call f_routine(id='SPRED_read_uinp')

    if(present(mpi_env_))then
        mpi_env=mpi_env_
    else
        mpi_env%iproc=0
        mpi_env%nproc=1
    endif

    call dict_init(dict)

    ! fill dictionary with hard coded default values
    call SPRED_init_input_dict(dict)

    !Now read a user input yaml file
    ! try radical.yaml. If radical is empty, use input.yaml
    if (len_trim(radical) == 0) then
       fname(1:len(fname)) = "input.yaml"
    else
       fname(1:len(fname)) = trim(radical) // ".yaml"
    end if
    inquire(file = trim(fname), exist = exists_user)
    if (exists_user)  call merge_input_file_to_dict(dict, trim(fname),mpi_env)

    if (mpi_env%iproc == 0) then
    call input_file_dump(dict)
    endif

    !now transfer the information from the dictionary to the SPRED input data
    !structure
    call SPRED_fill_variables(inputs,dict)

    ! We put a barrier here to be sure that non master proc will be stopped
    ! by any issue on the master proc.
    if(present(mpi_env_))  call mpibarrier(comm=mpi_env%mpi_comm)

    call dict_free(dict)
    call f_release_routine()
  end subroutine SPRED_read_uinp

  subroutine merge_input_file_to_dict(dict, fname, mpi_env)
    use SPREDbase
    use dictionaries
    !use yaml_output, only :yaml_map
    use yaml_parse, only: yaml_parse_from_char_array
    use yaml_output
    implicit none
    !Arguments
    type(dictionary), pointer :: dict            !< Dictionary of the input files. Should be initialized on entry
    character(len = *), intent(in) :: fname      !< Name of the file where the dictionary has to be read from
    type(mpi_environment), intent(in) :: mpi_env !< Environment of the reading. Used for broadcasting the result
    !local variables
    integer(kind = 8) :: cbuf, cbuf_len
    integer :: ierr
    character(len = max_field_length) :: val
    character, dimension(:), allocatable :: fbuf
    type(dictionary), pointer :: udict
    external :: getFileContent,copyCBuffer,freeCBuffer

    call f_routine(id='merge_input_file_to_dict')
    if (mpi_env%iproc == 0) then
       call getFileContent(cbuf, cbuf_len, fname, len_trim(fname))
    end if

    if (mpi_env%nproc > 1) call mpibcast(cbuf_len,comm=mpi_env%mpi_comm)
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')

    if (mpi_env%iproc == 0) then
       call copyCBuffer(fbuf, cbuf, cbuf_len)
       call freeCBuffer(cbuf)
!!       if (mpi_env%nproc > 1 .and. cbuf_len > 0) &
!!            & call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER, 0,
!!            mpi_env%mpi_comm, ierr)
!!    else
!!       if (cbuf_len > 0) call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER,
!!       0, mpi_env%mpi_comm, ierr)
    end if

    !this call can be replaced with the size of the character array
    if (mpi_env%nproc > 1) call mpibcast(fbuf,comm=mpi_env%mpi_comm)

    call f_err_open_try()
    call yaml_parse_from_char_array(udict, fbuf)
    call f_free_str(1,fbuf)
    ! Handle with possible partial dictionary.
    if (dict_len(udict) > 0) then
       call dict_update(dict, udict // 0)
    end if
    call dict_free(udict)
    ierr = 0
    if (f_err_check()) ierr = f_get_last_error(val)
    !call f_dump_all_errors()
    call f_err_close_try()

    if (ierr /= 0) call f_err_throw(err_id = ierr, err_msg = val)
    call f_release_routine()


  END SUBROUTINE merge_input_file_to_dict


  !> initializes the input dictionary
  !! by populating it from hard-coded source file
  subroutine SPRED_init_input_dict(dict)
    use dictionaries
    use f_input_file
    use yaml_parse
    implicit none
    !>input dictionary, a copy of the user input, to be filled
    !!with all the variables on exit
    type(dictionary), pointer :: dict
    !local variables
    integer(f_integer) :: params_size
    !integer(kind = 8) :: cbuf_add !< address of c buffer
    character, dimension(:), allocatable :: params
    type(dictionary), pointer :: parameters
    type(dictionary), pointer :: parsed_parameters
    type(dictionary), pointer :: profiles
    type(dictionary), pointer :: nested,asis

    call f_routine(id='SPRED_init_input_dict')

    nullify(parameters,parsed_parameters,profiles)

    !alternative filling of parameters from hard-coded source file
    !call getstaticinputdef(cbuf_add,params_size)
    call getspredinputdefsize(params_size)
    !allocate array
    params=f_malloc_str(1,params_size,id='params')
    !fill it and parse dictionary
    call getspredinputdef(params)

    call yaml_parse_from_char_array(parsed_parameters,params)
    !there is only one document in the input variables specifications
    parameters=>parsed_parameters//0
    profiles => parsed_parameters//1
    call f_free_str(1,params)

    call input_file_complete(parameters,dict,imports=profiles)

!    if (present(dict_minimal)) then
!       nullify(nested,asis)
!       call input_file_minimal(parameters,dict,dict_minimal,nested,asis)
!    end if

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

  end subroutine SPRED_init_input_dict

  subroutine SPRED_fill_variables(inputs,dict)
    use dictionaries
    implicit none
    type(SPRED_inputs), intent(out) :: inputs
    type(dictionary), pointer :: dict
    !local variables
    type(dictionary), pointer :: lvl,var

    ! Transfer dict values into input_variables structure.
    lvl => dict_iter(dict)
    do while(associated(lvl))
       var => dict_iter(lvl)
       do while(associated(var))
          call SPRED_input_fill(inputs,dict_key(lvl),var)
          var => dict_next(var)
       end do
       lvl => dict_next(lvl)
    end do

  end subroutine SPRED_fill_variables

  !> Set the dictionary from the input variables
  subroutine SPRED_input_fill(inputs, level, val)
    use SPREDbase
    use SPRED_public_keys
    use SPRED_public_enums
    use yaml_output, only: yaml_warning
    use dictionaries
    implicit none
    type(SPRED_inputs), intent(inout) :: inputs
    type(dictionary), pointer :: val
    character(len = *), intent(in) :: level
    !local variables
    real(gp) :: dummy_d
    character(len = max_field_length) :: str
    integer :: ipos

    if (index(dict_key(val), "_attributes") > 0) return

    select case(trim(level))
    case(FP_VARIABLES)
       select case (trim(dict_key(val)))
       case(FP_METHOD)
          str=val
          select case(trim(str))
          case('OMF')
             inputs%fp_method=OMF_FP_METHOD
          case('OMP')
             inputs%fp_method=OMP_FP_METHOD
          case('OMPOLD')
             inputs%fp_method=OMPOLD_FP_METHOD
          case('OMSOLD')
             inputs%fp_method=OMSOLD_FP_METHOD
          end select
       case(FP_NATX_SPHERE)
          inputs%fp_natx_sphere=val
       case(FP_ANGMOM)
          inputs%fp_angmom=val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case(GLBL_VARIABLES)
       select case (trim(dict_key(val)))
       case(GLBL_NWRITE)
          inputs%glbl_nwrite=val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case DEFAULT
          call yaml_warning("unknown input section '" // trim(level) // "'")
    end select
  END SUBROUTINE SPRED_input_fill


end module SPREDtypes
