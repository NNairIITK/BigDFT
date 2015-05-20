!> @file
!!  Module to store all dictionary keys of the input files.
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Define all static strings to store input variables
module module_input_keys
  use dictionaries
  use public_keys
  implicit none
  
  private

  !public :: input_keys_init, input_keys_finalize
  public :: input_keys_fill_all, input_keys_dump

  type(dictionary), pointer :: parameters=>null()
  type(dictionary), pointer :: parsed_parameters=>null()


contains

  subroutine input_keys_init()
    use yaml_output
    use dynamic_memory
    use yaml_parse
    implicit none
    !local variables
    integer :: params_size
    !integer(kind = 8) :: cbuf_add !< address of c buffer
    character, dimension(:), allocatable :: params

    !alternative filling of parameters from hard-coded source file
    !call getstaticinputdef(cbuf_add,params_size)
    call getinputdefsize(params_size)
    !allocate array
    params=f_malloc_str(1,params_size,id='params')
    !fill it and parse dictionary
    !print *,'after', f_loc(params),f_loc(params(1)),'shape',shape(params),params_size
    !print *,'cbuf_add',cbuf_add
    call getinputdef(params)
    !write(*,*)'here definition'
    !write(*,'('//trim(yaml_toa(params_size))//'a)')params

    !call copycbuffer(params,cbuf_add,params_size)
    !print *,'there',params_size
    call yaml_parse_from_char_array(parsed_parameters,params)
    !there is only one document in the input variables specifications
    parameters=>parsed_parameters//0
    call f_free_str(1,params)

    !call yaml_dict_dump(parameters, comment_key = COMMENT)
    
!!$    !in the case the errors have not been initialized before
!!$    call input_keys_errors()

  END SUBROUTINE input_keys_init


  subroutine input_keys_finalize()
    use dictionaries
    implicit none

    if (associated(parsed_parameters)) then
       call dict_free(parsed_parameters)
       nullify(parameters)
    else
       call dict_free(parameters)
    end if
  END SUBROUTINE input_keys_finalize

  !> Fill all the input keys into dict
  subroutine input_keys_fill_all(dict,dict_minimal)
    use dynamic_memory
    use module_defs, only: gp, pi_param
    use f_input_file
    use public_keys
    use yaml_strings, only: operator(.eqv.)
    !use yaml_output
    implicit none
    type(dictionary), pointer :: dict,dict_minimal
    !local variables
    type(dictionary), pointer :: as_is,nested,no_check
    character(max_field_length) :: meth, prof
    real(gp) :: dtmax_, betax_
    logical :: user_defined

    if (f_err_raise(.not. associated(dict),'The input dictionary has to be associated',&
         err_name='BIGDFT_RUNTIME_ERROR')) return

    call f_routine(id='input_keys_fill_all')

    ! Overriding the default for isolated system
    if (POSINP .in. dict) then
       if (.not.has_key(dict//POSINP,ASTRUCT_CELL) .and. .not. has_key(dict//DFT_VARIABLES,DISABLE_SYM)) then
          call set(dict // DFT_VARIABLES // DISABLE_SYM,.true.)
       end if
    end if
    nested=>list_new(.item. LIN_BASIS_PARAMS)

    ! Check and complete dictionary.
    call input_keys_init()


    call input_file_complete(parameters,dict,nocheck=nested)

!!$    call input_keys_fill(parameters,dict, PERF_VARIABLES)
!!$    call input_keys_fill(parameters,dict, MODE_VARIABLES)
!!$    call input_keys_fill(parameters,dict, DFT_VARIABLES)
!!$    call input_keys_fill(parameters,dict, KPT_VARIABLES)
!!$    call input_keys_fill(parameters,dict, GEOPT_VARIABLES)
!!$    call input_keys_fill(parameters,dict, MIX_VARIABLES)
!!$    call input_keys_fill(parameters,dict, SIC_VARIABLES)
!!$    call input_keys_fill(parameters,dict, TDDFT_VARIABLES)
!!$    call input_keys_fill(parameters,dict, LIN_GENERAL)
!!$    call input_keys_fill(parameters,dict, LIN_BASIS)
!!$    call input_keys_fill(parameters,dict, LIN_KERNEL)
!!$    call input_keys_fill(parameters,dict, LIN_BASIS_PARAMS, check=.false.)

    !create a shortened dictionary which will be associated to the given run
    !call input_minimal(dict,dict_minimal)
    as_is =>list_new(.item. FRAG_VARIABLES,.item. IG_OCCUPATION, .item. POSINP, .item. OCCUPATION) 
    call input_file_minimal(parameters,dict,dict_minimal,nested,as_is)
    call dict_free(nested,as_is)

    ! Additional treatments.
    meth = dict // GEOPT_VARIABLES // GEOPT_METHOD
    if (trim(meth) .eqv. "FIRE") then
       !prof = input_keys_get_source(dict // GEOPT_VARIABLES, DTMAX, user_defined)
       !if (trim(prof) == DEFAULT .and. .not. user_defined) then
       if (input_value_is_default(dict // GEOPT_VARIABLES, DTMAX)) then
          betax_ = dict // GEOPT_VARIABLES // BETAX
          call set(dict // GEOPT_VARIABLES // DTMAX, 0.25 * pi_param * sqrt(betax_), fmt = "(F7.4)")
       end if
       !prof = input_keys_get_source(dict // GEOPT_VARIABLES, DTINIT, user_defined)
       !if (trim(prof) == DEFAULT .and. .not. user_defined) then
       if (input_value_is_default(dict // GEOPT_VARIABLES, DTINIT)) then
          dtmax_ = dict // GEOPT_VARIABLES // DTMAX
          call set(dict // GEOPT_VARIABLES // DTINIT, 0.5 * dtmax_, fmt = "(F7.4)")
       end if
    end if

    call input_keys_finalize()

    call f_release_routine()
  end subroutine input_keys_fill_all

  !> takes the posinp filename from the dictionary. Starting point is dict//POSINP
  subroutine astruct_dict_get_source(dict, source)
    use public_keys, only: POSINP_SOURCE
    use f_utils, only: f_zero
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(inout) :: source !<preserve previous value if present
    !local variables
    type(dictionary), pointer :: dict_tmp

    call f_zero(source)
    dict_tmp=dict .get. ASTRUCT_PROPERTIES
    source=dict_tmp .get. POSINP_SOURCE

!!$    write(source, "(A)") ""
!!$    if (has_key(dict, ASTRUCT_PROPERTIES)) then
!!$       if (has_key(dict // ASTRUCT_PROPERTIES, POSINP_SOURCE)) &
!!$            & source = dict_value(dict // ASTRUCT_PROPERTIES // POSINP_SOURCE)
!!$    end if
  end subroutine astruct_dict_get_source

  !> Dump the dictionary of the input variables.
  !! Should dump only the keys relative to the input variables and
  !! print out warnings for the ignored keys
  subroutine input_keys_dump(dict, userOnly)
    use yaml_strings, only: f_strcpy
    use yaml_output
    use f_input_file, only: input_file_dump
    use dynamic_memory
    use public_keys, only: POSINP,ASTRUCT_PROPERTIES,ASTRUCT_POSITIONS
    
    implicit none
    type(dictionary), pointer :: dict   !< Dictionary to dump
    logical, intent(in), optional :: userOnly

    !local variables
    integer, parameter :: natoms_dump=500
    integer :: i, dlen, skeys,natoms
    character(max_field_length), dimension(:), allocatable :: keys
    character(max_field_length) ::  sourcefile
    logical :: userOnly_
    type(dictionary), pointer :: tmp

    call f_routine(id='input_keys_dump')

    !new mechanism, to see if it works

    userOnly_ = .false.
    if (present(userOnly)) userOnly_ = userOnly

    !create nodump_list
    tmp => list_new(.item. POSINP)
    call input_file_dump(dict, useronly=userOnly_,nodump_list=tmp)

    call dict_free(tmp)
    !then treat separately the posinp list
    call f_strcpy(src='not provided',dest=sourcefile)
    tmp=dict .get. POSINP
    if ( .not. associated(tmp)) return
    !check the number of atoms
    natoms=-1
    if (ASTRUCT_POSITIONS .in. tmp) natoms=dict_len(tmp // ASTRUCT_POSITIONS)
    if (natoms > natoms_dump) then
       call astruct_dict_get_source(tmp, sourcefile)
       call yaml_map(POSINP,sourcefile)
    else
       call yaml_mapping_open(POSINP)
       call input_file_dump(TMP, useronly=userOnly_,msg='Atomic positions')
       call yaml_mapping_close()
    end if
    nullify(tmp)


!!$    call yaml_comment("Input parameters", hfill = "-")
!!$    !TEST (the first dictionary has no key)
!!$    !if (.not. associated(dict%parent)) then
!!$    if (associated(dict%child)) then
!!$       if (dict_len(dict) >= 1) then
!!$          ! List case.
!!$          dlen = dict_len(dict)
!!$          do i = 0,  dlen- 1, 1
!!$             call input_variable_dump(dict // i,userOnly_)
!!$          end do
!!$       else
!!$          ! Dictionary case
!!$          allocate(keys(dict_size(dict)))
!!$          keys = dict_keys(dict)
!!$          skeys = size(keys)
!!$          do i = 1, skeys
!!$             if (POSINP == trim(keys(i))) then
!!$                call f_strcpy(src='not provided',dest=sourcefile)
!!$                tmp=>dict//POSINP
!!$                !check the number of atoms
!!$                natoms=-1
!!$                if (ASTRUCT_POSITIONS .in. tmp) natoms=dict_len(tmp // ASTRUCT_POSITIONS)
!!$                if (natoms > natoms_dump) then
!!$                   call astruct_dict_get_source(tmp, sourcefile)
!!$                   call yaml_map(POSINP,sourcefile)
!!$                else
!!$                   call input_variable_dump(dict // keys(i),userOnly_)
!!$                end if
!!$             else
!!$                call input_variable_dump(dict // keys(i),userOnly_)
!!$             end if
!!$          end do
!!$          deallocate(keys)
!!$       end if
!!$    else
!!$       call yaml_scalar(dict%data%value)
!!$    end if

    call f_release_routine()

  end subroutine input_keys_dump

end module module_input_keys
