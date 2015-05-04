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

  public :: input_keys_init, input_keys_finalize
  public :: input_keys_fill_all, input_keys_dump
!  public :: input_keys_equal,input_keys_get_profiles

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


!!!!  subroutine input_keys_dump_def(fname, file)
!!!!    use dictionaries
!!!!    use yaml_output
!!!!    use f_utils, only: f_get_free_unit
!!!!    implicit none
!!!!    character(len = *), intent(in) :: fname
!!!!    character(len = *), intent(in), optional :: file !< Subsection of the input to be printed (old input.file)
!!!!    !local variables
!!!!    integer, parameter :: unt
!!!!    integer :: ierr !, iunit_def
!!!!    !integer :: iunit_def
!!!!
!!!!    unt=f_get_free_unit(789159) !< To be sure is not opened
!!!!    
!!!!    ! Switch YAML output stream (not needed anymore)
!!!!    !call yaml_get_default_stream(iunit_def)
!!!!    call yaml_set_stream(unit=unt,filename=trim(fname),tabbing=0,record_length=100,istat=ierr,setdefault=.false.)
!!!!
!!!!    if (ierr == 0) then
!!!!       call input_keys_init()
!!!!       if (present(file)) then
!!!!          if (file .in. parameters) then
!!!!             call yaml_dict_dump(parameters // file,unit=unt)
!!!!          else
!!!!             call yaml_dict_dump(parameters,unit=unt)
!!!!          end if
!!!!       else
!!!!          call yaml_dict_dump(parameters,unit=unt)
!!!!       end if
!!!!       call input_keys_finalize()
!!!!    end if
!!!!    call yaml_close_stream(unit=unt)
!!!!    ! Set back normal YAML output (not needed anymore)
!!!!    !call yaml_set_default_stream(iunit_def,ierr)
!!!!  end subroutine input_keys_dump_def
!!!!
!!!!
!!!!  !> Get for each keys available profiles.
!!!!  function input_keys_get_profiles(file)
!!!!    use dictionaries
!!!!    implicit none
!!!!    character(len = *), intent(in) :: file
!!!!    type(dictionary), pointer :: input_keys_get_profiles
!!!!
!!!!    type(dictionary), pointer :: p
!!!!    integer :: i, skeys
!!!!    character(max_field_length), dimension(:), allocatable :: keys
!!!!
!!!!    call input_keys_init()
!!!!    
!!!!    call dict_init(p)
!!!!    
!!!!    if (has_key(parameters, file)) then
!!!!       call vars(p, parameters // file)
!!!!    else
!!!!       allocate(keys(dict_size(parameters)))
!!!!       keys = dict_keys(parameters)
!!!!       skeys = size(keys)
!!!!       do i = 1, skeys, 1
!!!!          call vars(p // keys(i), parameters // keys(i))
!!!!       end do
!!!!       deallocate(keys)
!!!!    end if
!!!!
!!!!    call input_keys_finalize()
!!!!
!!!!    input_keys_get_profiles => p
!!!!
!!!!  contains
!!!!
!!!!    subroutine vars(dict, ref)
!!!!      use dictionaries
!!!!      implicit none
!!!!      type(dictionary), pointer :: dict, ref
!!!!
!!!!      integer :: i, svar
!!!!      character(max_field_length), dimension(:), allocatable :: var
!!!!
!!!!      allocate(var(dict_size(ref)))
!!!!      var = dict_keys(ref)
!!!!      !to avoid problems on BG/Q
!!!!      svar = size(var)
!!!!      do i = 1, svar, 1
!!!!         call generate(dict, var(i), ref // var(i))
!!!!      end do
!!!!      deallocate(var)
!!!!      if (dict_size(dict) == 0) then
!!!!         call set(dict, "no profile")
!!!!      end if
!!!!    end subroutine vars
!!!!
!!!!    subroutine generate(dict, key, ref)
!!!!      use dictionaries
!!!!      implicit none
!!!!      type(dictionary), pointer :: dict, ref
!!!!      character(max_field_length), intent(in) :: key
!!!!
!!!!      integer :: i, skeys
!!!!      character(max_field_length), dimension(:), allocatable :: keys
!!!!
!!!!      allocate(keys(dict_size(ref)))
!!!!      keys = dict_keys(ref)
!!!!      !to avoid problems on BG/Q
!!!!      skeys = size(keys)
!!!!      do i = 1, skeys, 1
!!!!         if (trim(keys(i)) /= COMMENT .and. &
!!!!              & trim(keys(i)) /= COND .and. &
!!!!              & trim(keys(i)) /= RANGE .and. &
!!!!              & trim(keys(i)) /= PROF_KEY .and. &
!!!!              & trim(keys(i)) /= EXCLUSIVE .and. &
!!!!              & trim(keys(i)) /= DEFAULT) then
!!!!            call add(dict // key // "profiles", "'" // trim(keys(i)) // "'")
!!!!            !call dict_copy(dict // key // "profiles" // keys(i), ref // keys(i))
!!!!!!$         else if (trim(keys(i)) == EXCLUSIVE) then
!!!!!!$            ! Add exclusive values.
!!!!!!$            call dict_copy(dict // key // "allowed values", ref // EXCLUSIVE)
!!!!!!$         else if (trim(keys(i)) == RANGE) then
!!!!!!$            ! Add range definition
!!!!!!$            call dict_copy(dict // key // "within range", ref // RANGE)
!!!!!!$         else if (trim(keys(i)) == DEFAULT) then
!!!!!!$            ! Add range definition
!!!!!!$            call dict_copy(dict // key // "default value", ref // DEFAULT)
!!!!         end if
!!!!      end do
!!!!      deallocate(keys)
!!!!    end subroutine generate
!!!!  END FUNCTION input_keys_get_profiles
!!!!
!!!!
!!!!  !> Compare two strings (case-insensitive). Blanks are relevant!
!!!!  function input_keys_equal(stra,strb)
!!!!    implicit none
!!!!    character(len=*), intent(in) :: stra,strb
!!!!    logical :: input_keys_equal
!!!!    !Local variables
!!!!    integer :: i,ica,icb,ila,ilb,ilength
!!!!    ila=len(stra)
!!!!    ilb=len(strb)
!!!!    ilength=min(ila,ilb)
!!!!    ica=ichar(stra(1:1))
!!!!    icb=ichar(strb(1:1))
!!!!    input_keys_equal=(modulo(ica-icb,32) == 0) .and. (ila==ilb)
!!!!    do i=2,ilength
!!!!       ica=ichar(stra(i:i))
!!!!       icb=ichar(strb(i:i))
!!!!       input_keys_equal=input_keys_equal .and. &
!!!!            &   (modulo(ica-icb,32) == 0)
!!!!       if (.not. input_keys_equal) exit
!!!!    end do
!!!!  END FUNCTION input_keys_equal
!!!!


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
    type(dictionary), pointer :: as_is,nested
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

    ! Check and complete dictionary.
    call input_keys_init()

    !check for some fields that the user did not specify any unsupported key
    !now this mechanism is embedded in each of the checks
    !call input_keys_control(dict,DFT_VARIABLES)

    call input_keys_fill(parameters,dict, PERF_VARIABLES)
    call input_keys_fill(parameters,dict, MODE_VARIABLES)
    call input_keys_fill(parameters,dict, DFT_VARIABLES)
    call input_keys_fill(parameters,dict, KPT_VARIABLES)
    call input_keys_fill(parameters,dict, GEOPT_VARIABLES)
    call input_keys_fill(parameters,dict, MIX_VARIABLES)
    call input_keys_fill(parameters,dict, SIC_VARIABLES)
    call input_keys_fill(parameters,dict, TDDFT_VARIABLES)
    call input_keys_fill(parameters,dict, LIN_GENERAL)
    call input_keys_fill(parameters,dict, LIN_BASIS)
    call input_keys_fill(parameters,dict, LIN_KERNEL)
    call input_keys_fill(parameters,dict, LIN_BASIS_PARAMS, check=.false.)

    !create a shortened dictionary which will be associated to the given run
    !call input_minimal(dict,dict_minimal)
    nested=>list_new(.item. LIN_BASIS_PARAMS)
    as_is =>list_new(.item. FRAG_VARIABLES,.item. IG_OCCUPATION, .item. POSINP, .item. OCCUPATION) 
    call input_file_minimal(parameters,dict,dict_minimal,nested,as_is)
    call dict_free(nested)
    call dict_free(as_is)

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

!!$  subroutine input_control_callback()
!!$    use yaml_output
!!$    use dictionaries
!!$    implicit none
!!$    call f_err_severe()
!!$  end subroutine input_control_callback

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
    use f_input_file, only: input_variable_dump
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

    userOnly_ = .false.
    if (present(userOnly)) userOnly_ = userOnly

    call yaml_comment("Input parameters", hfill = "-")
    !TEST (the first dictionary has no key)
    !if (.not. associated(dict%parent)) then
    if (associated(dict%child)) then
       if (dict_len(dict) >= 1) then
          ! List case.
          dlen = dict_len(dict)
          do i = 0,  dlen- 1, 1
             call input_variable_dump(dict // i,userOnly_)
          end do
       else
          ! Dictionary case
          allocate(keys(dict_size(dict)))
          keys = dict_keys(dict)
          skeys = size(keys)
          do i = 1, skeys
             if (POSINP == trim(keys(i))) then
                call f_strcpy(src='not provided',dest=sourcefile)
                tmp=>dict//POSINP
                !check the number of atoms
                natoms=-1
                if (ASTRUCT_POSITIONS .in. tmp) natoms=dict_len(tmp // ASTRUCT_POSITIONS)
                if (natoms > natoms_dump) then
                   call astruct_dict_get_source(tmp, sourcefile)
                   call yaml_map(POSINP,sourcefile)
                else
                   call input_variable_dump(dict // keys(i),userOnly_)
                end if
             else
                call input_variable_dump(dict // keys(i),userOnly_)
             end if
          end do
          deallocate(keys)
       end if
    else
       call yaml_scalar(dict%data%value)
    end if

    call f_release_routine()

!!$  contains
!!$
!!$    recursive subroutine dict_dump_(dict)
!!$      use yaml_output
!!$      use dictionaries
!!$      implicit none
!!$      type(dictionary), pointer :: dict
!!$
!!$      logical :: flow, userDef
!!$      integer :: i, dlen, skeys
!!$      type(dictionary), pointer :: parent, attr, iter
!!$      character(max_field_length) :: descr, tag, prof, output
!!$      character(max_field_length), dimension(:), allocatable :: keys
!!$
!!$      if (index(dict%data%key, ATTRS) > 0) return
!!$
!!$      descr = " "
!!$      tag = " "
!!$      prof = " "
!!$      userDef = .false.
!!$      parent => dict%parent
!!$      if (associated(parent)) then
!!$         if (has_key(parent, trim(dict%data%key) // ATTRS)) then
!!$            attr => parent // (trim(dict%data%key) // ATTRS)
!!$            if (has_key(attr, COMMENT)) descr = attr // COMMENT
!!$            !if (has_key(attr, PROF_KEY)) tag = attr // PROF_KEY
!!$            if (has_key(attr, PROF_KEY)) prof = attr // PROF_KEY
!!$            if (has_key(attr, USER_KEY)) userDef = attr // USER_KEY
!!$         end if
!!$      end if
!!$      
!!$      if (dict_len(dict) > 0) then
!!$         ! List case.
!!$         if (userOnly_ .and. .not.userDef .and. trim(dict%data%key) /= "") return
!!$
!!$         flow = (.not.associated(dict%child%child))
!!$         if (.not.flow .and. trim(descr) /= "") then
!!$            call yaml_sequence_open(trim(dict%data%key), tag = tag, advance = "no")
!!$            call yaml_comment(trim(descr), tabbing = 50)
!!$         else
!!$            call yaml_sequence_open(trim(dict%data%key), tag = tag, flow=flow)
!!$         end if
!!$         dlen = dict_len(dict)
!!$         do i = 0, dlen - 1, 1
!!$            call yaml_sequence("", advance = "no")
!!$            call dict_dump_(dict // i)
!!$         end do
!!$         if (flow .and. trim(descr) /= "") then
!!$            call yaml_sequence_close(advance = "no")
!!$            call yaml_comment(trim(descr), tabbing = 50)
!!$         else
!!$            call yaml_sequence_close()
!!$         end if
!!$      else if (dict_size(dict) > 0) then
!!$         ! Dictionary case
!!$         if (userOnly_ .and. .not.userDef) return
!!$
!!$         if (len_trim(dict%data%key) > 0) &
!!$              & call yaml_mapping_open(trim(dict%data%key),flow=.false.)
!!$         iter => dict_next(dict)
!!$         allocate(keys(dict_size(dict)))
!!$         keys = dict_keys(dict)
!!$         skeys = size(keys)
!!$         do i = 1, skeys, 1
!!$            call dict_dump_(dict // keys(i))
!!$         end do
!!$         deallocate(keys)
!!$         if (len_trim(dict%data%key) > 0) call yaml_mapping_close()
!!$      else if (associated(dict)) then
!!$         ! Leaf case.
!!$         if (dict%data%item >= 0) then
!!$            ! List entry
!!$            call yaml_sequence(trim(dict%data%value))
!!$         else
!!$            ! Dictionary entry
!!$            if (userOnly_ .and. .not.userDef) return
!!$
!!$            if (userOnly_ .and. trim(prof) /= "") then
!!$               output = prof
!!$            else
!!$               output = dict%data%value
!!$            end if
!!$            if (trim(descr) /= "") then
!!$               call yaml_map(trim(dict%data%key), trim(output), tag = tag, advance = "no")
!!$               call yaml_comment(trim(descr), tabbing = 50)
!!$            else
!!$               call yaml_map(trim(dict%data%key), trim(output), tag = tag)
!!$            end if
!!$         end if
!!$      end if
!!$
!!$    end subroutine dict_dump_
  end subroutine input_keys_dump

end module module_input_keys
