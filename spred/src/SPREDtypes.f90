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
    !> fingerprint input variables
    type(f_enumerator) :: fp_method
    integer :: fp_natx_sphere       !< number of atoms in each sphere (for periodic fingerprint)
    integer :: fp_angmom       !< angular momentum of gaussian orbitals for overlap matrix fingerprints (both periodic and free BC)
  end type SPRED_inputs


  public :: SPRED_input_dict

contains
  subroutine SPRED_init(inputs)
    implicit none
  end subroutine SPRED_init 

  subroutine merge_input_file_to_dict(dict, fname)
    use module_base
    !use yaml_output, only :yaml_map
    use yaml_parse, only: yaml_parse_from_char_array
    use yaml_output
    implicit none
    !Arguments
    type(dictionary), pointer :: dict            !< Dictionary of the input files. Should be initialized on entry
    character(len = *), intent(in) :: fname      !< Name of the file where the dictionary has to be read from 
    !local variables
    integer(kind = 8) :: cbuf, cbuf_len
    integer :: ierr
    character(len = max_field_length) :: val
    character, dimension(:), allocatable :: fbuf
    type(dictionary), pointer :: udict
    external :: getFileContent,copyCBuffer,freeCBuffer


    call f_routine(id='merge_input_file_to_dict')
       call getFileContent(cbuf, cbuf_len, fname, len_trim(fname))

    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')

       call copyCBuffer(fbuf, cbuf, cbuf_len)
       call freeCBuffer(cbuf)
    end if


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


  !> This function returns a dictionary with all the input variables of a SPRED run filled.
  !! This dictionary is constructed from an updated version of the input variables dictionary
  !! following the input files as defined by the user
  subroutine read_input_dict_from_files(radical,dict)
    use SPRED_base
    !use module_input_dicts, only: merge_input_file_to_dict
    use input_old_text_format
    !use yaml_output
    implicit none
    character(len = *), intent(in) :: radical    !< The name of the run. use "input" if empty
    type(dictionary), pointer :: dict            !< Input dictionary, has to be nullified at input
    !local variables
    integer :: ierr
    logical :: exists_default, exists_user
    character(len = max_field_length) :: fname

    call f_routine(id='read_input_dict_from_files')


    ! We try first default.yaml
    inquire(file = "default.yaml", exist = exists_default)
    if (exists_default) call merge_input_file_to_dict(dict, "default.yaml")
    ! We try then radical.yaml
    if (len_trim(radical) == 0 .or. trim(radical) == LOGFILE) then
       fname(1:len(fname)) = "input.yaml"
    else
       fname(1:len(fname)) = trim(radical) // ".yaml"
    end if
    inquire(file = trim(fname), exist = exists_user)
    if (exists_user) call merge_input_file_to_dict(dict, trim(fname))

    call f_release_routine()
  end subroutine read_input_dict_from_files

  !>routine to fill the input variables of spred
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

    call f_routine(id='SPRED_input_dict')

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

  subroutine SPRED_fill_variables(k,inputs,dict)
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
          call SPRED_input_fill(k,inputs,dict_key(lvl),var)
          var => dict_next(var)
       end do
       lvl => dict_next(lvl)
    end do

  end subroutine SPRED_fill_variables

  !> Set the dictionary from the input variables
  subroutine SPRED_input_fill(inputs, level, val)
    use SPREDbase
    use yaml_output, only: yaml_warning
    use dictionaries
    implicit none
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
    case(FP_VARIABLES)
       select case (trim(dict_key(val)))
       case(FP_METHOD)
          inputs%fp_method=val
       case(FP_NATX_SPHERE)
          inputs%fp_natx_sphere=val
       case(FP_ANGMOM)
          inputs%fp_angmom=val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case DEFAULT
    end select
  END SUBROUTINE SPRED_input_fill


end module SPREDtypes
