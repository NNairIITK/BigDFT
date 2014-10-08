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
  use dictionaries, only: dictionary
  use public_keys
  implicit none
  
  private

  !> Error ids for this module.
  integer, parameter :: ERR_UNDEF=1
  integer, public :: INPUT_VAR_NOT_IN_LIST = ERR_UNDEF
  integer, public :: INPUT_VAR_NOT_IN_RANGE = ERR_UNDEF
  integer, public :: INPUT_VAR_ILLEGAL = ERR_UNDEF
  type(dictionary), pointer :: failed_exclusive

!!$  character(len = *), parameter :: RANGE = "__range__", EXCLUSIVE = "__exclusive__"
!!$  character(len = *), parameter :: DEFAULT = "default", COMMENT = "__comment__"
!!$  character(len = *), parameter :: COND = "__condition__", WHEN = "__when__"
!!$  character(len = *), parameter :: MASTER_KEY = "__master_key__"
  character(len = *), parameter :: ATTRS = "_attributes"
  character(len = *), parameter :: PROF_KEY = "PROFILE_FROM"
  character(len = *), parameter :: USER_KEY = "USER_DEFINED"

  character(len = *), parameter :: COMMENT = "COMMENT"
  character(len = *), parameter :: DESCRIPTION = "DESCRIPTION"
  character(len = *), parameter :: RANGE = "RANGE"
  character(len = *), parameter :: EXCLUSIVE = "EXCLUSIVE"
  character(len = *), parameter :: DEFAULT = "default"
  character(len = *), parameter :: COND = "CONDITION"
  character(len = *), parameter :: WHEN = "WHEN"
  character(len = *), parameter :: MASTER_KEY = "MASTER_KEY"

  public :: input_keys_init, input_keys_finalize
  public :: input_keys_fill_all, input_keys_dump
  public :: input_keys_equal, input_keys_get_source, input_keys_dump_def
  public :: input_keys_get_profiles,input_keys_validate
  public :: input_keys_errors

  type(dictionary), pointer :: parameters=>null()
  type(dictionary), pointer :: parsed_parameters=>null()


contains


  !> Callback routine when an error occurs
  subroutine abort_excl()
    use yaml_output
    use dictionaries
    implicit none

    call f_dump_last_error()
    call yaml_mapping_open("allowed values")
    call yaml_dict_dump(failed_exclusive)
    call yaml_mapping_close()
    call f_err_severe()
  end subroutine abort_excl


  !> Callback routine for illegal input variables
  subroutine warn_illegal()
    implicit none
    
  end subroutine warn_illegal


  subroutine input_keys_init()
    use yaml_output
    use dictionaries
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
    
    !in the case the errors have not been initialized before
    call input_keys_errors()

  END SUBROUTINE input_keys_init

  !> define the errors of the module
  subroutine input_keys_errors()
    use dictionaries, only: f_err_define
    implicit none

    if (INPUT_VAR_NOT_IN_LIST == ERR_UNDEF) then
       call f_err_define(err_name='INPUT_VAR_NOT_IN_LIST',&
            err_msg='given value not in allowed list.',&
            err_action='choose a value from the list below.',&
            err_id=INPUT_VAR_NOT_IN_LIST,callback=abort_excl)
    end if
    if (INPUT_VAR_NOT_IN_RANGE == ERR_UNDEF) then
       call f_err_define(err_name='INPUT_VAR_NOT_IN_RANGE',&
            err_msg='given value not in allowed range.',&
            err_action='adjust the given value.',&
            err_id=INPUT_VAR_NOT_IN_RANGE)
    end if
    if (INPUT_VAR_ILLEGAL == ERR_UNDEF) then
       call f_err_define(err_name='INPUT_VAR_ILLEGAL',&
            err_msg='provided variable is not allowed in this context.',&
            err_action='correct or remove the input variable.',&
            err_id=INPUT_VAR_ILLEGAL,callback=warn_illegal)
    end if

    
  end subroutine input_keys_errors

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


  subroutine input_keys_dump_def(fname, file)
    use dictionaries
    use yaml_output
    implicit none
    character(len = *), intent(in) :: fname
    character(len = *), intent(in), optional :: file !< Subsection of the input to be printed (old input.file)
    !local variables
    integer, parameter :: unt=789159 !< To be sure is not opened
    integer :: ierr !, iunit_def
    !integer :: iunit_def

    ! Switch YAML output stream (not needed anymore)
    !call yaml_get_default_stream(iunit_def)
    call yaml_set_stream(unit=unt,filename=trim(fname),tabbing=0,record_length=100,istat=ierr,setdefault=.false.)

    if (ierr == 0) then
       call input_keys_init()
       if (present(file)) then
          if (has_key(parameters, file)) then
             call yaml_dict_dump(parameters // file,unit=unt)
          else
             call yaml_dict_dump(parameters,unit=unt)
          end if
       else
          call yaml_dict_dump(parameters,unit=unt)
       end if
       call input_keys_finalize()
    end if
    call yaml_close_stream(unit=unt)
    ! Set back normal YAML output (not needed anymore)
    !call yaml_set_default_stream(iunit_def,ierr)
  end subroutine input_keys_dump_def


  !> Get for each keys available profiles.
  function input_keys_get_profiles(file)
    use dictionaries
    implicit none
    character(len = *), intent(in) :: file
    type(dictionary), pointer :: input_keys_get_profiles

    type(dictionary), pointer :: p
    integer :: i, skeys
    character(max_field_length), dimension(:), allocatable :: keys

    call input_keys_init()
    
    call dict_init(p)
    
    if (has_key(parameters, file)) then
       call vars(p, parameters // file)
    else
       allocate(keys(dict_size(parameters)))
       keys = dict_keys(parameters)
       skeys = size(keys)
       do i = 1, skeys, 1
          call vars(p // keys(i), parameters // keys(i))
       end do
       deallocate(keys)
    end if

    call input_keys_finalize()

    input_keys_get_profiles => p

  contains

    subroutine vars(dict, ref)
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict, ref

      integer :: i, svar
      character(max_field_length), dimension(:), allocatable :: var

      allocate(var(dict_size(ref)))
      var = dict_keys(ref)
      !to avoid problems on BG/Q
      svar = size(var)
      do i = 1, svar, 1
         call generate(dict, var(i), ref // var(i))
      end do
      deallocate(var)
      if (dict_size(dict) == 0) then
         call set(dict, "no profile")
      end if
    end subroutine vars

    subroutine generate(dict, key, ref)
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict, ref
      character(max_field_length), intent(in) :: key

      integer :: i, skeys
      character(max_field_length), dimension(:), allocatable :: keys

      allocate(keys(dict_size(ref)))
      keys = dict_keys(ref)
      !to avoid problems on BG/Q
      skeys = size(keys)
      do i = 1, skeys, 1
         if (trim(keys(i)) /= COMMENT .and. &
              & trim(keys(i)) /= COND .and. &
              & trim(keys(i)) /= RANGE .and. &
              & trim(keys(i)) /= PROF_KEY .and. &
              & trim(keys(i)) /= EXCLUSIVE .and. &
              & trim(keys(i)) /= DEFAULT) then
            call add(dict // key // "profiles", "'" // trim(keys(i)) // "'")
            !call dict_copy(dict // key // "profiles" // keys(i), ref // keys(i))
!!$         else if (trim(keys(i)) == EXCLUSIVE) then
!!$            ! Add exclusive values.
!!$            call dict_copy(dict // key // "allowed values", ref // EXCLUSIVE)
!!$         else if (trim(keys(i)) == RANGE) then
!!$            ! Add range definition
!!$            call dict_copy(dict // key // "within range", ref // RANGE)
!!$         else if (trim(keys(i)) == DEFAULT) then
!!$            ! Add range definition
!!$            call dict_copy(dict // key // "default value", ref // DEFAULT)
         end if
      end do
      deallocate(keys)
    end subroutine generate
  END FUNCTION input_keys_get_profiles


  !> Compare two strings (case-insensitive). Blanks are relevant!
  function input_keys_equal(stra,strb)
    implicit none
    character(len=*), intent(in) :: stra,strb
    logical :: input_keys_equal
    !Local variables
    integer :: i,ica,icb,ila,ilb,ilength
    ila=len(stra)
    ilb=len(strb)
    ilength=min(ila,ilb)
    ica=ichar(stra(1:1))
    icb=ichar(strb(1:1))
    input_keys_equal=(modulo(ica-icb,32) == 0) .and. (ila==ilb)
    do i=2,ilength
       ica=ichar(stra(i:i))
       icb=ichar(strb(i:i))
       input_keys_equal=input_keys_equal .and. &
            &   (modulo(ica-icb,32) == 0)
       if (.not. input_keys_equal) exit
    end do
  END FUNCTION input_keys_equal


  function input_keys_get_source(dict, key, user_defined)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: key
    logical, intent(out) :: user_defined

    character(len = max_field_length) :: input_keys_get_source

    user_defined = .false.
    input_keys_get_source(1:max_field_length) = DEFAULT
    if (has_key(dict, trim(key) // ATTRS)) then
       if (has_key(dict // (trim(key) // ATTRS), USER_KEY)) &
            & user_defined = dict // (trim(key) // ATTRS) // USER_KEY
       if (has_key(dict // (trim(key) // ATTRS), PROF_KEY)) &
            & input_keys_get_source = dict // (trim(key) // ATTRS) // PROF_KEY
    end if
  end function input_keys_get_source

  !> this routine controls that the keys which are defined in the 
  !! input dictionary are all valid.
  !! in case there are some keys which are different, raise an error
  subroutine input_keys_validate(dict)
    use dictionaries
    use yaml_output
    use module_base, only: bigdft_mpi
    implicit none
    type(dictionary), pointer :: dict
    !local variables
    logical :: found
    type(dictionary), pointer :: valid_entries,valid_patterns
    type(dictionary), pointer :: iter,invalid_entries,iter2
    

    !> fill the list of valid entries
    valid_entries=>list_new([&
         .item. RADICAL_NAME, &
         .item. POSINP,&  
         .item. PERF_VARIABLES,&  
         .item. DFT_VARIABLES,&   
         .item. KPT_VARIABLES,&   
         .item. GEOPT_VARIABLES,& 
         .item. MIX_VARIABLES,&   
         .item. SIC_VARIABLES,&   
         .item. TDDFT_VARIABLES,& 
         .item. LIN_GENERAL,&     
         .item. LIN_BASIS,&       
         .item. LIN_KERNEL,&      
         .item. LIN_BASIS_PARAMS,&
         .item. OCCUPATION,&
         .item. IG_OCCUPATION,&
         .item. FRAG_VARIABLES])
    
    !then the list of vaid patterns
    valid_patterns=>list_new(&
         .item. 'psppar' &
         )

    call dict_init(invalid_entries)
    !for any of the keys of the dictionary iterate to find if it is allowed
    iter=>dict_iter(dict)
    do while(associated(iter))
       if ((valid_entries .index. dict_key(iter)) < 0) then
          found=.false.
          iter2=>dict_iter(valid_patterns)
          !check also if the key contains the allowed patterns
          find_patterns: do while(associated(iter2))
             if (index(dict_key(iter),trim(dict_value(iter2))) > 0) then
                found=.true.
                exit find_patterns
             end if
             iter2=>dict_next(iter2)
          end do find_patterns
          if (.not. found) call add(invalid_entries,dict_key(iter))
       end if
       iter=>dict_next(iter)
    end do

    if (dict_len(invalid_entries) > 0) then
       if (bigdft_mpi%iproc==0) then
          call yaml_map('Allowed keys',valid_entries)
          call yaml_map('Allowed key patterns',valid_patterns)
          call yaml_map('Invalid entries of the input dictionary',invalid_entries)
       end if
       call f_err_throw('The input dictionary contains invalid entries,'//&
            ' check above the valid entries',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    end if

    call dict_free(invalid_entries)
    call dict_free(valid_entries)
    call dict_free(valid_patterns)

  end subroutine input_keys_validate

  !> Fill all the input keys into dict
  subroutine input_keys_fill_all(dict,dict_minimal)
    use dictionaries
    use dynamic_memory
    use module_defs, only: gp, pi_param
    use yaml_output
    use public_keys
    !use yaml_output
    implicit none
    type(dictionary), pointer :: dict,dict_minimal

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

    call input_keys_fill(dict, PERF_VARIABLES)
    call input_keys_fill(dict, DFT_VARIABLES)
    call input_keys_fill(dict, KPT_VARIABLES)
    call input_keys_fill(dict, GEOPT_VARIABLES)
    call input_keys_fill(dict, MIX_VARIABLES)
    call input_keys_fill(dict, SIC_VARIABLES)
    call input_keys_fill(dict, TDDFT_VARIABLES)
    call input_keys_fill(dict, LIN_GENERAL)
    call input_keys_fill(dict, LIN_BASIS)
    call input_keys_fill(dict, LIN_KERNEL)
    call input_keys_fill(dict, LIN_BASIS_PARAMS, check=.false.)

    !create a shortened dictionary which will be associated to the given run
    call input_minimal(dict,dict_minimal)

    ! Additional treatments.
    meth = dict // GEOPT_VARIABLES // GEOPT_METHOD
    if (input_keys_equal(trim(meth), "FIRE")) then
       prof = input_keys_get_source(dict // GEOPT_VARIABLES, DTMAX, user_defined)
       if (trim(prof) == DEFAULT .and. .not. user_defined) then
          betax_ = dict // GEOPT_VARIABLES // BETAX
          call set(dict // GEOPT_VARIABLES // DTMAX, 0.25 * pi_param * sqrt(betax_), fmt = "(F7.4)")
       end if
       prof = input_keys_get_source(dict // GEOPT_VARIABLES, DTINIT, user_defined)
       if (trim(prof) == DEFAULT .and. .not. user_defined) then
          dtmax_ = dict // GEOPT_VARIABLES // DTMAX
          call set(dict // GEOPT_VARIABLES // DTINIT, 0.5 * dtmax_, fmt = "(F7.4)")
       end if
    end if

    call input_keys_finalize()

    call f_release_routine()
  end subroutine input_keys_fill_all


  !> This routine is used to create a minimal dictionary which can be used at the place 
  !! of the one provided as an indication on the understood variables
  subroutine input_minimal(dict,minimal)
    use dictionaries
    use yaml_output
    implicit none
    type(dictionary), pointer, intent(in) :: dict
    type(dictionary), pointer, intent(out) :: minimal
    !local variables
    type(dictionary), pointer :: dict_tmp,min_cat,as_is
    character(len=max_field_length) :: category
    logical :: cat_found

    nullify(minimal)

    !recursively search into the reference input variables

    dict_tmp => dict_iter(parameters)

    do while(associated(dict_tmp))
       !for any of the keys of parameters look at the corresponding value of the dictionary
       category=dict_key(dict_tmp)
       !call yaml_map('dict category',parameters//category)
       !print *,'category',trim(category),has_key(dict,category)
       !call yaml_map('dict category',dict_tmp)
       if (has_key(dict,category)) then
          call minimal_category(dict_tmp,dict//category,min_cat)
          if (associated(min_cat)) then
             if (.not. associated(minimal)) call dict_init(minimal)
             call set(minimal//category,min_cat)
          end if
       end if
!stop
       dict_tmp => dict_next(dict_tmp)
    end do

    !then add other information to the minimal dictionary which is associated
    !to specific system parameters
    !! basis set
    if (LIN_BASIS_PARAMS .in. dict) then
      dict_tmp => dict_iter(dict//LIN_BASIS_PARAMS)
      do while(associated(dict_tmp))
       category=dict_key(dict_tmp)
       !Pb with stack (Cray - ftn 05/2015), solved with the cat_found temporary variable
       cat_found = category .in. parameters//LIN_BASIS_PARAMS
       if (.not. cat_found .and. index(category,ATTRS) == 0 ) then
          !verify that no parameters correspond to default values
          call minimal_category(parameters//LIN_BASIS_PARAMS,dict_tmp,min_cat)
          if (associated(min_cat)) then
             !call dict_copy(minimal//LIN_BASIS_PARAMS//category,dict_tmp)
             call set(minimal//LIN_BASIS_PARAMS//category,min_cat)
          end if
       end if
       dict_tmp => dict_next(dict_tmp)
      end do
   end if

   as_is=>list_new(.item. FRAG_VARIABLES,.item. IG_OCCUPATION, .item. POSINP, .item. OCCUPATION) 

   !fragment dictionary has to be copied as-is
   !other variables have to follow the same treatment
   dict_tmp => dict_iter(as_is)
   do while(associated(dict_tmp))
      category=dict_value(dict_tmp)
      if (category .in. dict) &
           call dict_copy(minimal//category,dict//category)
      dict_tmp => dict_next(dict_tmp)
   end do
   call dict_free(as_is)

    contains
      
      subroutine minimal_category(vars,input,minim)
        implicit none
        type(dictionary), pointer :: vars,input,minim
        !local variables
        logical :: profile_found
        character(len=max_field_length) :: def_var,var_prof,prof_var
        type(dictionary), pointer :: defvar,var
        nullify(minim)

        var=>dict_iter(vars)
!        call yaml_map('var dict',var)

        do while(associated(var))
           def_var=dict_key(var)
!           print *,'here2 ',trim(def_var),has_key(input,def_var)
           
           !search if the input data have values among the profiles
           if (has_key(input,def_var)) then
              profile_found=.false.
              !see if the dictionary has the PROF_KEY in thier possibilities
              prof_var(1:len(prof_var))=' '
              if (has_key(var,PROF_KEY)) then
                 if (has_key(input,dict_value(var//PROF_KEY))) &
                    prof_var=dict_value(input//dict_value(var//PROF_KEY))
              end if

              defvar => dict_iter(var)
              !              call yaml_map('var dict inside',defvar)
              check_profile: do while(associated(defvar))
                 !exclude keys for definition of the variable
                 var_prof=dict_key(defvar)
!              call yaml_map('key',var_prof)
                 if (trim(var_prof) /= COMMENT .and. trim(var_prof) /= COND .and.&
                      trim(var_prof) /= RANGE .and. trim(var_prof) /= PROF_KEY .and. &
                      trim(var_prof) /= EXCLUSIVE) then
                    !check if some profile meets desired values
!call yaml_map('defvar',defvar)
!call yaml_map('input',input//def_var)
!call yaml_map('result',defvar == input//def_var)
!call yaml_map('var_prof',var_prof)
!call yaml_map('test',trim(var_prof) /= DEFAULT .and. var_prof /= prof_var)
!print *,'key',def_var
                   profile_found= (defvar == input//def_var) 
                   if (profile_found) then
                       if (trim(var_prof) /= DEFAULT .and. var_prof /= prof_var) then
                          if (.not. associated(minim)) call dict_init(minim)
                          call set(minim//def_var,var_prof)
                       end if
                       exit check_profile
                    end if
                 end if
                 defvar => dict_next(defvar)
              end do check_profile
              !the key has not been found among the profiles, therefore it should be entered as is
              if (.not. profile_found .and. len_trim(dict_value(input//def_var))/=0) then
                 if (.not. associated(minim)) call dict_init(minim)
                 !clean the list items if the dictionary is a list with all the items identical
                 defvar => dict_iter(input//def_var)
                 var_prof=repeat(' ',len(var_prof))
                 if (dict_len(input // def_var)==0) nullify(defvar)
                 compact_list: do while(associated(defvar))
                    !if scalar, retrieve the value, otherwise exit
                    if (dict_size(defvar) == 0 .and. dict_len(defvar)==0) then
                       prof_var=defvar
                    else
                       var_prof=repeat(' ',len(var_prof))
                       exit compact_list
                    end if
                    !check if all the values of the list are equal to the first one
                    if (len_trim(var_prof) == 0) then
                       var_prof=prof_var
                    else
                       !check if all the values are OK, otherwise exit at first failure
                       if (var_prof /= prof_var) then
                         var_prof=repeat(' ',len(var_prof))
                         exit compact_list
                      end if
                    end if
                    defvar => dict_next(defvar)
                 end do compact_list
                 !if the dictionary is not a one-level list or if it is a list with different values
                 ! copy it as-is
                 if (len_trim(var_prof) == 0) then
                    call dict_copy(minim//def_var,input//def_var)
                 else !otherwise put the scalar value associated
                    call set(minim//def_var,var_prof)
                 end if
              end if
           end if
           var => dict_next(var)
        end do
      end subroutine minimal_category
      
    end subroutine input_minimal


  subroutine input_keys_fill(dict, file,check)
    use dictionaries
    use dynamic_memory
    use yaml_output
    implicit none
    !Arguments
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file
    logical, intent(in), optional :: check
    !Local variables
    !integer :: i
    logical :: user, hasUserDef,docheck
    type(dictionary), pointer :: ref,ref_iter
    !character(len=max_field_length), dimension(:), allocatable :: keys
    
    docheck=.true.
    if (present(check)) docheck=check

!    call f_routine(id='input_keys_fill')
    if (docheck) call input_keys_control(dict,file)

    ref => parameters // file

    ref_iter => dict_iter(parameters // file)
    hasUserDef = .false.
    do while(associated(ref_iter))
       if (trim(dict_key(ref_iter)) /= DESCRIPTION) then
          call input_keys_set(user, dict // file, file, dict_key(ref_iter))
          hasUserDef = (hasUserDef .or. user)
       end if
       ref_iter=> dict_next(ref_iter)
    end do

!    keys = f_malloc_str(max_field_length,dict_size(ref),id='keys')
!    keys = dict_keys(ref)
!    hasUserDef = .false.
!    do i=1,size(keys)
!       call input_keys_set(user, dict // file, file, keys(i))
!        hasUserDef = (hasUserDef .or. user)
!    end do
!    call set(dict // (trim(file) // ATTRS) // USER_KEY, hasUserDef)
!
!    call f_free_str(max_field_length, keys)
!    call f_release_routine()
  END SUBROUTINE input_keys_fill


  !> control if all the keys which are defined in a given field are associated with a true input variable
  subroutine input_keys_control(dict,file)
    use dictionaries
    use yaml_output, only: yaml_map,yaml_toa,yaml_warning
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file
    !local variables
    type(dictionary), pointer :: dict_tmp,ref,dict_err,dict_it
    
    ref=> parameters // file
    !parse all the keys of the dictionary
    dict_tmp=>dict_iter(dict//file)
    do while(associated(dict_tmp))
       if (.not. (dict_key(dict_tmp) .in. ref) .and. &
            & index(dict_key(dict_tmp), ATTRS) == 0) then
    !      call yaml_map('Allowed keys',dict_keys(ref))
          !even in a f_err_open_try section this error is assumed to be fatal
          !for the moment. A mechanism to downgrade its gravity should be
          !provided; -> LG 30/06/14: should work in try environment at present
          dict_it=>dict_iter(ref)
          call dict_init(dict_err)
          do while(associated(dict_it))
             if (trim(dict_key(dict_it)) /= DESCRIPTION) &
                  call add(dict_err,dict_key(dict_it))
             dict_it=>dict_next(dict_it)
          end do
          !dict_err=>list_new(.item. dict_keys(ref))
          call yaml_warning('Input file, section "'//file//&
            '"; invalid key "'//trim(dict_key(dict_tmp))//'".')
          call yaml_map('Allowed keys',dict_err)
          call dict_free(dict_err)
          call f_err_throw('An invalid key ('//trim(dict_key(dict_tmp))&
              //') has been found in section "'&
                //file//'". Check above the allowed keys.' ,&
            err_id=INPUT_VAR_ILLEGAL,callback=f_err_severe)
       end if
       dict_tmp=> dict_next(dict_tmp)
    end do
  end subroutine input_keys_control


  subroutine input_control_callback()
    use yaml_output
    use dictionaries
    implicit none
    call f_err_severe()
  end subroutine input_control_callback


  subroutine input_keys_set(userDef, dict, file, key)
    use dictionaries
    use yaml_output
    use dynamic_memory
    implicit none
    logical, intent(out) :: userDef
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file, key

    integer :: i, skeys
    type(dictionary), pointer :: ref
    character(len = max_field_length) :: val, profile_
    character(len = max_field_length), dimension(:), allocatable :: keys
    double precision, dimension(2) :: rg
    logical :: found

!    call f_routine(id='input_keys_set')

    ref => parameters // file // key

    profile_(1:max_field_length) = " "
    if (trim(profile_) == "") profile_(1:max_field_length) = DEFAULT
    
    userDef = (has_key(dict, key))
    if (userDef) then
       ! Key should be present only for some unmet conditions.
       if (.not.set_(dict, ref)) then
          !to see if the f_release_routine has to be controlled automatically
         ! call f_release_routine() !to be called before raising the error
!!$       if (f_err_raise(.not.set_(dict, ref), err_id = INPUT_VAR_ILLEGAL, &
!!$            & err_msg = trim(file) // "/" // trim(key) // " is not allowed in this context.")) then
          call f_err_throw(err_id = INPUT_VAR_ILLEGAL, &
               & err_msg = trim(file) // "/" // trim(key) // " is not allowed in this context.")
          return
       end if
       val = dict // key

       ! There is already a value in dict.
       if (has_key(ref, val)) then
          ! The key was asking for a profile, we copy it.
          call dict_copy(dict // key, ref // val)
          profile_ = val
       else
          ! The value is user-defined, we validate it with info from dict from.
          if (has_key(ref, RANGE)) then
             rg(1) = ref // RANGE // 0
             rg(2) = ref // RANGE // 1
             call validate(dict // key, key, rg)
          else if (has_key(ref, EXCLUSIVE)) then
             failed_exclusive => ref // EXCLUSIVE
             allocate(keys(dict_size(failed_exclusive)))
             keys = dict_keys(failed_exclusive)
             found = .false.
             skeys = size(keys)
             do i = 1, skeys, 1
                found = input_keys_equal(trim(val), trim(keys(i)))
                if (found) exit
             end do
             deallocate(keys)
             if (f_err_raise(.not. found, err_id = INPUT_VAR_NOT_IN_LIST, &
                  & err_msg = trim(key) // " = '" // trim(val) // "' is not allowed.")) return
             nullify(failed_exclusive)
          end if
       end if
    else
       ! Key should be present only for some unmet conditions.
       if (.not.set_(dict, ref)) then
!          call f_release_routine()
          return
       end if
       ! Hard-coded profile from key.
       if (has_key(ref, PROF_KEY)) then
          val = ref // PROF_KEY
          if (has_key(dict, val)) then
             profile_ = dict // val
          end if
       end if

       ! There is no value in dict, we take it from ref.
       if (.not. has_key(ref, profile_)) profile_ = DEFAULT
       call dict_copy(dict // key, ref // profile_)
    end if

    ! Copy the comment.
    if (has_key(ref, COMMENT)) &
         & call dict_copy(dict // (trim(key) // ATTRS) // COMMENT, ref // COMMENT)
    ! Save the source.
    if (userDef) &
         & call set(dict // (trim(key) // ATTRS) // USER_KEY, .true.)
    if (profile_ /= DEFAULT) &
         & call set(dict // (trim(key) // ATTRS) // PROF_KEY, profile_)

!    call f_release_routine()
  contains

    function set_(dict, ref)
      implicit none
      type(dictionary), pointer :: dict, ref
      logical :: set_

      integer :: j, dlen
      type(dictionary), pointer :: tmp
      character(max_field_length) :: mkey, val_master, val_when

      set_ = .true.
      if (has_key(ref, COND)) then
         mkey = ref // COND // MASTER_KEY
         if (.not. has_key(dict, mkey)) then
            set_ = .false.
            return
         end if
         val_master = dict // mkey
         set_ = .false.
         tmp => ref // COND // WHEN
         do j = 0, dict_len(tmp) - 1, 1
            val_when = tmp // j
            set_ = set_ .or. &
                 & (input_keys_equal(trim(val_master), trim(val_when)))
         end do
      end if
    end function set_

    recursive subroutine validate(dict, key, rg)
      implicit none
      type(dictionary), pointer :: dict
      character(len = *), intent(in) :: key
      double precision, dimension(2), intent(in) :: rg
      
      character(len = max_field_length) :: val
      character(max_field_length), dimension(:), allocatable :: keys
      double precision :: var
      integer :: dlen, skeys

      if (associated(dict%child)) then
         if (dict_len(dict) >= 1) then
            ! List case.
            dlen = dict_len(dict)
            do i = 0, dlen - 1, 1
               call validate(dict // i, key, rg)
            end do            
         else
            ! Dictionary case
            allocate(keys(dict_size(dict)))
            keys = dict_keys(dict)
            skeys = size(keys)
            do i = 1, skeys, 1
               call validate(dict // keys(i), key, rg)
            end do
            deallocate(keys)
         end if
      else
         var = dict
         if (var < rg(1) .or. var > rg(2)) then
            val = dict
            call f_err_throw(err_id = INPUT_VAR_NOT_IN_RANGE, &
              & err_msg = trim(key) // " = '" // trim(val) // &
              & "' not in range.")
            return
         end if
      end if
    end subroutine validate
  END SUBROUTINE input_keys_set


  !> Dump the dictionary of the input variables.
  !! Should dump only the keys relative to the iunput variables and
  !! print out warnings for the ignored keys
  subroutine input_keys_dump(dict, userOnly)
    use yaml_output
    use dictionaries
    use dynamic_memory
    implicit none
    type(dictionary), pointer :: dict   !< Dictionary to dump
    logical, intent(in), optional :: userOnly

    !local variables
    integer :: i, dlen, skeys
    character(max_field_length), dimension(:), allocatable :: keys
    logical :: userOnly_

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
             call dict_dump_(dict // i)
          end do
       else
          ! Dictionary case
          allocate(keys(dict_size(dict)))
          keys = dict_keys(dict)
          skeys = size(keys)
          do i = 1, skeys, 1
             call dict_dump_(dict // keys(i))
          end do
          deallocate(keys)
       end if
    else
       call yaml_scalar(dict%data%value)
    end if

    call f_release_routine()

  contains

    recursive subroutine dict_dump_(dict)
      use yaml_output
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict

      logical :: flow, userDef
      integer :: i, dlen, skeys
      type(dictionary), pointer :: parent, attr, iter
      character(max_field_length) :: descr, tag, prof, output
      character(max_field_length), dimension(:), allocatable :: keys

      if (index(dict%data%key, ATTRS) > 0) return

      descr = " "
      tag = " "
      prof = " "
      userDef = .false.
      parent => dict%parent
      if (associated(parent)) then
         if (has_key(parent, trim(dict%data%key) // ATTRS)) then
            attr => parent // (trim(dict%data%key) // ATTRS)
            if (has_key(attr, COMMENT)) descr = attr // COMMENT
            !if (has_key(attr, PROF_KEY)) tag = attr // PROF_KEY
            if (has_key(attr, PROF_KEY)) prof = attr // PROF_KEY
            if (has_key(attr, USER_KEY)) userDef = attr // USER_KEY
         end if
      end if
      
      if (dict_len(dict) > 0) then
         ! List case.
         if (userOnly_ .and. .not.userDef .and. trim(dict%data%key) /= "") return

         flow = (.not.associated(dict%child%child))
         if (.not.flow .and. trim(descr) /= "") then
            call yaml_sequence_open(trim(dict%data%key), tag = tag, advance = "no")
            call yaml_comment(trim(descr), tabbing = 50)
         else
            call yaml_sequence_open(trim(dict%data%key), tag = tag, flow=flow)
         end if
         dlen = dict_len(dict)
         do i = 0, dlen - 1, 1
            call yaml_sequence("", advance = "no")
            call dict_dump_(dict // i)
         end do
         if (flow .and. trim(descr) /= "") then
            call yaml_sequence_close(advance = "no")
            call yaml_comment(trim(descr), tabbing = 50)
         else
            call yaml_sequence_close()
         end if
      else if (dict_size(dict) > 0) then
         ! Dictionary case
         if (userOnly_ .and. .not.userDef) return

         if (len_trim(dict%data%key) > 0) &
              & call yaml_mapping_open(trim(dict%data%key),flow=.false.)
         iter => dict_next(dict)
         allocate(keys(dict_size(dict)))
         keys = dict_keys(dict)
         skeys = size(keys)
         do i = 1, skeys, 1
            call dict_dump_(dict // keys(i))
         end do
         deallocate(keys)
         if (len_trim(dict%data%key) > 0) call yaml_mapping_close()
      else if (associated(dict)) then
         ! Leaf case.
         if (dict%data%item >= 0) then
            ! List entry
            call yaml_sequence(trim(dict%data%value))
         else
            ! Dictionary entry
            if (userOnly_ .and. .not.userDef) return

            if (userOnly_ .and. trim(prof) /= "") then
               output = prof
            else
               output = dict%data%value
            end if
            if (trim(descr) /= "") then
               call yaml_map(trim(dict%data%key), trim(output), tag = tag, advance = "no")
               call yaml_comment(trim(descr), tabbing = 50)
            else
               call yaml_map(trim(dict%data%key), trim(output), tag = tag)
            end if
         end if
      end if

    end subroutine dict_dump_
  end subroutine input_keys_dump

end module module_input_keys
