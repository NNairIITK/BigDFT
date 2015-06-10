!> @file
!>  Modules which contains all the interfaces to parse input dictionary.
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Modules which contains all interfaces to parse input dictionary.
module module_input_dicts
  use public_keys
  use dictionaries
  implicit none

  private

  ! Update a dictionary from a input file
  public :: merge_input_file_to_dict

  public :: occupation_set_from_dict
  public :: neb_set_from_dict

  ! Dictionaries from files (old formats).
  public :: occupation_data_file_merge_to_dict
  public :: dict_set_run_properties,dict_get_run_properties,dict_run_new,bigdft_options
  public :: set_dict_run_file,create_log_file,dict_run_validate

  !> Keys of a run dict. All private, use get_run_prop() and set_run_prop() to change them.
  character(len = *), parameter :: RADICAL_NAME = "radical"
  character(len = *), parameter :: INPUT_NAME   = "input_file"
  character(len = *), parameter :: OUTDIR       = "outdir"
  character(len = *), parameter :: LOGFILE      = "logfile"
  character(len = *), parameter :: USE_FILES    = "run_from_files"


contains

  subroutine bigdft_options(parser)
    use yaml_parse
    use dictionaries, only: dict_new,operator(.is.)
    implicit none
    type(yaml_cl_parse), intent(inout) :: parser

    call yaml_cl_parse_option(parser,'name','None',&
         'name of the run','n',&
         dict_new('Usage' .is. &
         'Name of the run. When <name> is given, input files like <name>.* are used. '//&
         'The file "default.yaml" set the default values. If name is given as a list in yaml format, '//&
         'this is interpreted as a list of runs as if a runs-file has been given.',&
         'Allowed values' .is. &
         'String value.'),first_option=.true.)

    call yaml_cl_parse_option(parser,OUTDIR,'.',&
         'output directory','d',&
         dict_new('Usage' .is. &
         'Set the directory where all the output files have to be written.',&
         'Allowed values' .is. &
         'String value, indicating the path of the directory. If the last subdirectory is not existing, it will be created'))

    call yaml_cl_parse_option(parser,LOGFILE,'No',&
         'create logfile','l',&
         dict_new('Usage' .is. &
         'When "Yes", write the result of the run in file "log.yaml" or "log-<name>.yaml" if the run has a specified name.',&
         'Allowed values' .is. &
         'Boolean (yaml syntax). Automatically set to true when using runs-file or output directory different from "."'))

    call yaml_cl_parse_option(parser,'runs-file','None',&
         'list_posinp filename','r',&
         dict_new('Usage' .is. &
         'File containing the list of the run ids which have to be launched independently (list in yaml format). '//&
         'The option runs-file is not compatible with the --name option.',&
         'Allowed values' .is. &
         'String value. Should be associated to a existing filename'),&
         conflicts='[name]')

    call yaml_cl_parse_option(parser,'taskgroup-size','None',&
         'mpi_groupsize (number of MPI runs for a single instance of BigDFT)','t',&
         dict_new('Usage' .is. &
         'Indicates the number of mpi tasks associated to a single instance of BigDFT run',&
         'Allowed values' .is. &
         'Integer value. Disabled if not a divisor of the total No. of MPI tasks.'))

  end subroutine bigdft_options


  !> default run properties
  subroutine dict_run_new(run)
    implicit none
    type(dictionary), pointer :: run

    run => dict_new(RADICAL_NAME .is. ' ')
  end subroutine dict_run_new

  subroutine set_dict_run_file(run_id,options)
    implicit none
    character(len=*), intent(in) :: run_id
    type(dictionary), pointer, intent(inout) :: options
    !local variables
    logical :: lval
    character(len=max_field_length) :: val
    type(dictionary), pointer :: drun

    call dict_run_new(drun)
    if (trim(run_id) /= "input" ) then
       call dict_set_run_properties(drun, run_id = trim(run_id))
    end if
    if (OUTDIR .in. options) then
       val = options // OUTDIR
       call dict_set_run_properties(drun, outdir_id = trim(val))
    end if
    if (LOGFILE .in. options) then
       lval = options // LOGFILE
       call dict_set_run_properties(drun, log_to_disk = lval)
    end if
    call set(drun // USE_FILES, .true.)
    call add(options//'BigDFT', drun)
  end subroutine set_dict_run_file


  !> set the parameters of the run 
  subroutine dict_set_run_properties(run,run_id,input_id,posinp_id, &
       & outdir_id,log_to_disk,run_from_files)
    use public_keys, only: POSINP
    implicit none
    type(dictionary), pointer :: run !< nullified if not initialized
    character(len=*), intent(in), optional :: run_id !< Radical of the run
    character(len=*), intent(in), optional :: input_id, posinp_id !< Input file name and posinp file name.
    character(len=*), intent(in), optional :: outdir_id !< Output directory (automatically add a trailing "/" if not any
    logical, intent(in), optional :: log_to_disk !< Write logfile to disk instead of screen.
    logical, intent(in), optional :: run_from_files !< Run_objects should be initialised from files.

    integer :: lgt

    if (.not. associated(run)) call dict_run_new(run)

    if (present(run_id)) then
       if (trim(run_id) /= "input") then
          call set(run // RADICAL_NAME, trim(run_id))
       else
          call set(run // RADICAL_NAME, " ")
       end if
    end if
    if (present(input_id)) call set(run // INPUT_NAME, trim(input_id))
    if (present(posinp_id)) call set(run // POSINP, trim(posinp_id))
    if (present(outdir_id)) then
       lgt = len_trim(outdir_id)
       if (outdir_id(lgt:lgt) == "/") then
          call set(run // OUTDIR, trim(outdir_id))
       else
          call set(run // OUTDIR, trim(outdir_id) // "/")
       end if
    end if
    if (present(log_to_disk)) call set(run // LOGFILE, log_to_disk)
    if (present(run_from_files)) call set(run // USE_FILES, run_from_files)

  end subroutine dict_set_run_properties

  !> get the parameters of the run 
  subroutine dict_get_run_properties(run,run_id,input_id,posinp_id,naming_id, &
       & outdir_id,log_to_disk,run_from_files)
    use public_keys, only: POSINP
    use f_utils, only: f_zero
    use yaml_strings, only: f_strcpy
    implicit none
    type(dictionary), pointer :: run
    character(len=*), intent(out), optional :: run_id, naming_id
    character(len=*), intent(out), optional :: input_id, posinp_id
    character(len=*), intent(inout), optional :: outdir_id
    logical, intent(inout), optional :: log_to_disk
    logical, intent(inout), optional :: run_from_files

    if (present(input_id)) then
       if (INPUT_NAME .in. run) then
          input_id = run // INPUT_NAME
       else if (RADICAL_NAME .in. run) then
          input_id = run // RADICAL_NAME
       else
          call f_zero(input_id)
          !input_id = " "
       end if
       if (len_trim(input_id) == 0) call f_strcpy(src=&
            "input" // trim(run_id_toa()),dest=input_id)
    end if
    if (present(posinp_id)) then 
       if (POSINP .in. run) then
          posinp_id = run // POSINP
       else if (RADICAL_NAME .in. run) then
          posinp_id = run // RADICAL_NAME
       else
          call f_zero(posinp_id)
          !posinp_id = " "
       end if
       if (len_trim(posinp_id) == 0) call f_strcpy(src=&
            "posinp" // trim(run_id_toa()),dest=posinp_id)
    end if
    if (present(naming_id)) then
       if (RADICAL_NAME .in. run) then
          naming_id = run // RADICAL_NAME
       else
          call f_zero(naming_id)
          !naming_id = " "
       end if
       if (len_trim(naming_id) == 0) then
          !naming_id = trim(run_id_toa())
          call f_strcpy(src=trim(run_id_toa()),dest=naming_id)
       else
          call f_strcpy(src="-" // trim(naming_id),dest=naming_id)
          !naming_id = "-" // trim(naming_id)
       end if
    end if
    if (present(run_id)) then
       if (RADICAL_NAME .in. run) then
          run_id = run // RADICAL_NAME
       else
          call f_zero(run_id)
          !run_id = " "
       end if
    end if
    if (present(outdir_id) .and. has_key(run, OUTDIR)) outdir_id = run // OUTDIR
    if (present(log_to_disk) .and. has_key(run, LOGFILE)) log_to_disk = run // LOGFILE
    if (present(run_from_files) .and. has_key(run, USE_FILES)) run_from_files = run // USE_FILES

  end subroutine dict_get_run_properties

  function run_id_toa()
    use yaml_output, only: yaml_toa
    use module_base, only: bigdft_mpi
    implicit none
    character(len=20) :: run_id_toa

    run_id_toa=repeat(' ',len(run_id_toa))

    if (bigdft_mpi%ngroup>1) then
       run_id_toa=adjustl(trim(yaml_toa(bigdft_mpi%igroup,fmt='(i15)')))
    end if

  end function run_id_toa

  !> this routine controls that the keys which are defined in the 
  !! input dictionary are all valid.
  !! in case there are some keys which are different, raise an error
  subroutine dict_run_validate(dict)
    use dictionaries
    use yaml_output
    use yaml_strings, only: operator(.eqv.)
    use module_base, only: bigdft_mpi
    use public_keys, only: POSINP, PERF_VARIABLES, DFT_VARIABLES, KPT_VARIABLES, &
         & GEOPT_VARIABLES, MIX_VARIABLES, SIC_VARIABLES, TDDFT_VARIABLES, LIN_GENERAL, &
         & LIN_BASIS, LIN_KERNEL, LIN_BASIS_PARAMS, OCCUPATION, IG_OCCUPATION, FRAG_VARIABLES, &
         & MODE_VARIABLES
    implicit none
    type(dictionary), pointer :: dict
    !local variables
    logical :: found,loginput
    type(dictionary), pointer :: valid_entries,valid_patterns
    type(dictionary), pointer :: iter,invalid_entries,iter2


    !> fill the list of valid entries
    valid_entries=>list_new([&
         .item. OUTDIR,&
         .item. RADICAL_NAME,&
         .item. USE_FILES,&
         .item. INPUT_NAME,&
         .item. LOGFILE,&
         .item. POSINP,&
         .item. MODE_VARIABLES,&
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

    !if the input file is a logfile, this will be only a warning
    loginput=.false.

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
          !if the key "Version number"  is present this means that we are parsing a logfile
          !which is tolerated
          loginput = (trim(dict_key(iter)) .eqv. 'Version Number' ) .or. loginput
       end if
       iter=>dict_next(iter)
    end do

    if (dict_len(invalid_entries) > 0 .and. .not. loginput) then
       if (bigdft_mpi%iproc==0) then
          call yaml_map('Allowed keys',valid_entries)
          call yaml_map('Allowed key patterns',valid_patterns)
          call yaml_map('Invalid entries of the input dictionary',invalid_entries)
       end if
       call f_err_throw('The input dictionary contains invalid entries,'//&
            ' check above the valid entries',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    !else if (loginput) then
    !      call yaml_warning('This input file has been created from a logfile')
    end if

    call dict_free(invalid_entries)
    call dict_free(valid_entries)
    call dict_free(valid_patterns)

  end subroutine dict_run_validate


  subroutine create_log_file(dict,dict_from_files)
    use module_base, enum_int => f_int
    use yaml_strings
    use yaml_output
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(out) :: dict_from_files !<identifies if the dictionary comes from files
    !local variables
    integer :: ierror,lgt,unit_log,ierrr
    integer(kind=4) :: ierr
    character(len = max_field_length) :: writing_directory, run_name
    character(len=500) :: logfilename,path
    integer :: iproc_node, nproc_node
    logical :: log_to_disk

    ! Get user input writing_directory.
    writing_directory = "."
    call dict_get_run_properties(dict, outdir_id = writing_directory)
    ! Create writing_directory and parents if needed and broadcast everything.
    if (trim(writing_directory) /= '.') then
       !path=repeat(' ',len(path))
       call f_zero(path)
       !add the output directory in the directory name
       if (bigdft_mpi%iproc == 0 .and. trim(writing_directory) /= '.') then
          call f_mkdir(writing_directory,path)
!!$          call getdir(writing_directory,&
!!$               int(len_trim(writing_directory),kind=4),path,int(len(path),kind=4),ierr)
!!$          if (ierr /= 0) then
!!$             write(*,*) "ERROR: cannot create writing directory '"&
!!$                  //trim(writing_directory) // "'."
!!$             call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
!!$          end if
       end if
!!$       call MPI_BCAST(path,len(path),MPI_CHARACTER,0,bigdft_mpi%mpi_comm,ierr)
!!$       lgt=min(len(writing_directory),len(path))
!!$       writing_directory(1:lgt)=path(1:lgt)
       call mpibcast(path,comm=bigdft_mpi%mpi_comm)
       call f_strcpy(src=path,dest=writing_directory)
    end if
    ! Add trailing slash if missing.
    lgt = len_trim(writing_directory)
    if (writing_directory(lgt:lgt) /= "/") &
         & writing_directory(min(lgt+1, len(writing_directory)):min(lgt+1, len(writing_directory))) = "/"

    ! Test if logging on disk is required.
    log_to_disk = (bigdft_mpi%ngroup > 1)
    call dict_get_run_properties(dict, log_to_disk = log_to_disk) !< May overwrite with user choice

    ! Save modified infos in dict.
    call dict_set_run_properties(dict, outdir_id = writing_directory, log_to_disk = log_to_disk)

    ! Now, create the logfile if needed.
    if (bigdft_mpi%iproc == 0) then
       if (log_to_disk) then
          ! Get Create log file name.
          call dict_get_run_properties(dict, naming_id = run_name)
          logfilename = "log" // trim(run_name) // ".yaml"
          path = trim(writing_directory)//trim(logfilename)
          call yaml_map('<BigDFT> log of the run will be written in logfile',path,unit=6)
          ! Check if logfile is already connected.
          call yaml_stream_connected(trim(path), unit_log, ierrr)
          if (ierrr /= 0) then
             ! Move possible existing log file.
             call ensure_log_file(trim(writing_directory), trim(logfilename), ierr)
             if (ierr /= 0) call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
             ! Close active stream and logfile if any. (TO BE MOVED IN RUN_UPDATE TO AVOID CLOSURE OF UPLEVEL INSTANCE)
             call yaml_get_default_stream(unit_log)
             if (unit_log /= 6) call yaml_close_stream(unit_log, ierrr)
             !Create stream and logfile
             call yaml_set_stream(filename=trim(path),record_length=92,istat=ierrr)
             !create that only if the stream is not already present, otherwise print a warning
             if (ierrr == 0) then
                call yaml_get_default_stream(unit_log)
                call old_input_stdout(unit_log)
                !call input_set_stdout(unit=unit_log)
             else
                call yaml_warning('Logfile '//trim(path)//' cannot be created, stream already present. Ignoring...')
             end if
          else
             call yaml_release_document(unit_log)
             call yaml_set_default_stream(unit_log, ierrr)
          end if ! Logfile already connected
       else
          !use stdout, do not crash if unit is present
          call yaml_set_stream(record_length=92,istat=ierrr)
       end if ! Need to create a named logfile.

       !start writing on logfile
       call yaml_new_document()
       !welcome screen
       call print_logo()
    end if ! Logfile is created by master proc only

    if (bigdft_mpi%nproc >1) call processor_id_per_node(bigdft_mpi%iproc,bigdft_mpi%nproc,iproc_node,nproc_node)

    if (bigdft_mpi%iproc==0) then
       if (bigdft_mpi%nproc >1) call yaml_map('MPI tasks of root process node',nproc_node)
       call print_configure_options()
    end if

    dict_from_files = .false.
    if (USE_FILES .in. dict) dict_from_files = dict // USE_FILES

  END SUBROUTINE create_log_file

  !> Routine to read YAML input files and create input dictionary.
  !! Update the input dictionary with the result of yaml_parse
  subroutine merge_input_file_to_dict(dict, fname, mpi_env)
    use module_base
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
       !if (mpi_env%nproc > 1) &
       !     & call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
    !else
       !call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
    end if

    if (mpi_env%nproc > 1) call mpibcast(cbuf_len,comm=mpi_env%mpi_comm)
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')

    if (mpi_env%iproc == 0) then
       call copyCBuffer(fbuf, cbuf, cbuf_len)
       call freeCBuffer(cbuf)
!       if (mpi_env%nproc > 1 .and. cbuf_len > 0) &
!            & call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
!    else
!       if (cbuf_len > 0) call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
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

  !subroutine aocc_to_dict(dict, nspin, noncoll, nstart, aocc, nelecmax, lmax, nsccode)
  !  use module_defs, only: gp
  !  use dictionaries
  !  implicit none
  !  integer, intent(in) :: nelecmax, lmax, nsccode, nspin, noncoll, nstart
  !  type(dictionary), pointer :: dict
  !  real(gp), dimension(nelecmax), intent(in) :: aocc

  !  type(dictionary), pointer :: val
  !  character(len = 4) :: key
  !  integer :: l, inl, nl, iocc, sccode, nsc, i
  !  character(len = 1), dimension(4), parameter :: lname = (/ "s", "p", "d", "f" /)

  !  call dict_init(dict)

  !  sccode = nsccode
  !  iocc=0
  !  do l = 1, lmax
  !     iocc=iocc+1
  !     ! Get number of shells for this channel 
  !     !(to be corrected, the rule is not the same)
  !     nl = int(aocc(iocc))
  !     ! Get number of semi cores for this channel
  !     nsc = modulo(sccode, 4)
  !     sccode = sccode / 4
  !     if (nl == 0) cycle
  !     do inl = 1, nl, 1
  !        if (inl <= nsc) then
  !           write(key, "(A1,I1,A1,A1)") "(", nstart + inl, lname(l), ")"
  !        else
  !           write(key, "(I1, A1)") nstart + inl, lname(l)
  !        end if
  !        call dict_init(val)
  !        do i = 1, nspin * noncoll * (2 * l - 1), 1
  !           iocc=iocc+1
  !           call add(val, aocc(iocc))
  !        end do
  !        call set(dict // key, val)
  !     end do
  !  end do
  !end subroutine aocc_to_dict


  subroutine occupation_set_from_dict(dict, key, norbu, norbd, occup, &
       & nkpts, nspin, norbsempty, qelec_up, qelec_down, norb_max)
    use module_defs, only: gp
    use dynamic_memory
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: key
    real(gp), dimension(:), pointer :: occup
    integer, intent(in) :: nkpts, nspin, norbsempty, norb_max
    real(gp), intent(in) :: qelec_up, qelec_down
    integer, intent(out) :: norbu, norbd

    integer :: norb
    integer :: ikpt,ne_up,ne_dwn
    type(dictionary), pointer :: occup_src
    character(len = 12) :: kpt_key
    call f_routine(id='occupation_set_from_dict')
    ! Default case.
    !integer approximation of the number of electrons
    ne_up=int_elec(qelec_up)
    !the same for down case
    ne_dwn=int_elec(qelec_down)

    if (nspin == 1) then
       !norb  = min((nelec_up + 1) / 2, norb_max)
       norb  = min((ne_up + 1) / 2, norb_max)
       norbu = norb
    else
       !norb = min(nelec_up + nelec_down, 2 * norb_max)
       norb = min(ne_up + ne_dwn, 2 * norb_max)
       if (nspin == 2) then
          !norbu = min(nelec_up, norb_max)
          norbu = min(ne_up, norb_max)
       else
          !norbu = min(nelec_up, 2 * norb_max)
          norbu = min(ne_up, 2 * norb_max)
       end if
    end if
    norbd = norb - norbu
!!$    write(*,*) nelec_up, nelec_down, norbsempty, norb_max
!!$    write(*,*) norbu, norbd, norb
!!$    stop
    ! Modify the default with occupation
    occup_src = dict .get. key
    !nullify(occup_src)
    !if (has_key(dict, key)) then
    if (associated(occup_src)) then
       !occup_src => dict //key
       ! Occupation is provided.
       if (has_key(occup_src, "K point 1")) then
          call count_for_kpt(occup_src // "K point 1")
       else if (nkpts == 1) then
          call count_for_kpt(occup_src)
       end if
       do ikpt = 2, nkpts, 1
          write(kpt_key, "(A)") "K point" // trim(yaml_toa(ikpt, fmt = "(I0)"))
          if (has_key(occup_src, kpt_key)) call count_for_kpt(occup_src // kpt_key)
       end do
    else if (norbsempty > 0) then
       !value of empty orbitals up and down, needed to fill occupation numbers
       if (nspin == 4 .or. nspin == 1) then
          norbu = norbu + min(norbsempty, norb_max - norbu)
       else if (nspin == 2) then
          norbu = norbu + min(norbsempty, norb_max - norbu)
          norbd = norbd + min(norbsempty, norb_max - norbd)
       end if
    end if

    ! Summarize and check.
    norb = norbu + norbd
    if (((nspin == 1 .or. nspin == 2) .and. (norbu > norb_max .or. norbd > norb_max)) &
         & .or. (nspin == 4 .and. (norbu > 2 * norb_max .or. norbd > 0))) then
       call yaml_warning('Total number of orbitals (found ' // trim(yaml_toa(norb)) &
            & // ') exceeds the available input guess orbitals (being ' &
            & // trim(yaml_toa(norb_max)) // ').')
       stop
    end if

    ! Allocate occupation accordingly.
    occup = f_malloc_ptr(norb * nkpts, id = "occup", routine_id = "occupation_set_from_dict")
    ! Setup occupation
    if (nspin==1) then
       do ikpt = 1, nkpts, 1
          call fill_default((ikpt - 1) * norb, 2, qelec_up, norb)
          if (associated(occup_src)) then
             write(kpt_key, "(A)") "K point" // trim(yaml_toa(ikpt, fmt = "(I0)"))
             if (ikpt == 0 .and. .not. has_key(occup_src, kpt_key)) then
                call fill_for_kpt((ikpt - 1) * norb, occup_src)
             else
                call fill_for_kpt((ikpt - 1) * norb, occup_src // kpt_key)
             end if
          end if
       end do
    else
       do ikpt = 0, nkpts - 1, 1
          call fill_default(ikpt * norb, 1, qelec_up, norbu)
          call fill_default(ikpt * norb + norbu, 1, qelec_down, norbd)
          if (associated(occup_src)) then
             write(kpt_key, "(A)") "K point" // trim(yaml_toa(ikpt, fmt = "(I0)"))
             if (ikpt == 0 .and. .not. has_key(occup_src, kpt_key)) then
                call fill_for_kpt((ikpt - 1) * norb, occup_src // "up")
                call fill_for_kpt((ikpt - 1) * norb + norbu, occup_src // "down")
             else
                call fill_for_kpt((ikpt - 1) * norb, occup_src // kpt_key // "up")
                call fill_for_kpt((ikpt - 1) * norb + norbu, occup_src // kpt_key // "down")
             end if
          end if
       end do
    end if

    !Check if sum(occup)=nelec
    if (abs(sum(occup) / nkpts - (qelec_up + qelec_down))>1.e-6_gp) then
       call f_err_throw('The total number of electrons ' &
            & // trim(yaml_toa(sum(occup) / nkpts,fmt='(f13.6)')) &
            & // ' is not equal to' // trim(yaml_toa(qelec_up + qelec_down)),&
            err_name='BIGDFT_INPUT_FILE_ERROR')
    end if

    call f_release_routine()

  contains

    pure function int_elec(qelec) result(ne)
      implicit none
      real(gp), intent(in) :: qelec
      integer :: ne
      
      ne=nint(qelec)
      !if we have an excess of electrons, add one orbital
      if (qelec - real(ne,gp) > 1.e-12_gp) ne=ne+1
    end function int_elec

    subroutine count_for_kpt(occ)
      implicit none
      type(dictionary), pointer :: occ
      
      if (nspin == 2) then
         if (.not. has_key(occ, "up") .or. &
              & .not. has_key(occ, "down")) stop "missing up or down"
         call count_orbs(norbu, occ // "up")
         call count_orbs(norbd, occ // "down")
      else
         call count_orbs(norbu, occ)
      end if
    end subroutine count_for_kpt

    subroutine count_orbs(n, occ)
      implicit none
      type(dictionary), pointer :: occ
      integer, intent(inout) :: n
      
      type(dictionary), pointer :: it
      character(len = max_field_length) :: key
      integer :: iorb

      it => dict_iter(occ)
      do while(associated(it))
         key = dict_key(it)
         read(key(index(key, " ") + 1:), *) iorb
         n = max(n, iorb)
         it => dict_next(it)
      end do
    end subroutine count_orbs

    subroutine fill_default(isorb, nfill, qelec, norb)
      implicit none
      integer, intent(in) :: isorb, nfill, norb
      real(gp), intent(in) :: qelec
      !local variables
      integer :: iorb, ne
      real(gp) :: rit,rnt

      rnt=0.0_gp
      !ne = (nelec + 1) / nfill
      ne = (int_elec(qelec) + 1) / nfill
      do iorb=isorb + 1, isorb + min(ne, norb)
         rit=min(real(nfill,gp),qelec-rnt)
         occup(iorb)=rit
         rnt=rnt+rit
      enddo
      do iorb=isorb+min(ne, norb)+1,isorb+norb
         occup(iorb)=0._gp
      end do
    end subroutine fill_default

    subroutine fill_for_kpt(isorb, occ)
      implicit none
      integer, intent(in) :: isorb
      type(dictionary), pointer :: occ

      type(dictionary), pointer :: it
      character(len = max_field_length) :: key
      integer :: iorb

      it => dict_iter(occ)
      do while(associated(it))
         key = dict_key(it)
         read(key(index(key, " ") + 1:), *) iorb
         occup(isorb + iorb) = it
         it => dict_next(it)
      end do
    end subroutine fill_for_kpt
  end subroutine occupation_set_from_dict


  subroutine occupation_data_file_merge_to_dict(dict, key, filename)
    use module_defs, only: gp, UNINITIALIZED
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key

    logical :: exists
    integer :: ierror, ntu, ntd, nt, i, iorb, lline, lstring
    character(len = 100) :: line, string
    type(dictionary), pointer :: valu, vald
    
    inquire(file = filename, exist = exists)
    if (.not. exists) return

    open(unit=91,file=filename,status='old',iostat=ierror)
    !Check the open statement
    if (ierror /= 0) then
       call yaml_warning('Failed to open the existing file '// trim(filename))
       stop
    end if

    !The first line gives the number of orbitals
    read(unit=91,fmt='(a100)') line

    read(line,fmt=*,iostat=ierror) ntu, ntd
    if (ierror /= 0) then
       !The first line gives the number of orbitals
       ntd = 0
       read(line,fmt=*,iostat=ierror) ntu
       if (ierror /=0) stop 'ERROR: reading the number of orbitals.'
    end if

    call dict_init(valu)
    if (ntd > 0) call dict_init(vald)

    nt = 1
    do
       read(unit=91,fmt='(a100)',iostat=ierror) line
       if (ierror /= 0) then
          exit
       end if
       !Transform the line in case there are slashes (to ease the parsing)
       lline = len(line)
       do i=1,lline
          if (line(i:i) == '/') then
             line(i:i) = ':'
          end if
       end do
       read(line,*,iostat=ierror) iorb,string
       if (ierror /= 0) then
          exit
       end if
       !Transform back the ':' into '/'
       lstring = len(string)
       do i=1,lstring
          if (string(i:i) == ':') then
             string(i:i) = '/'
          end if
       end do
       nt=nt+1

       if (iorb<0 .or. iorb>ntu + ntd) then
          !if (iproc==0) then
          write(*,'(1x,a,i0,a)') 'ERROR in line ',nt,' of the file "[name].occ"'
          write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
          !end if
          stop
       else
          if (iorb <= ntu) then
             call set(valu // ("Orbital" // trim(yaml_toa(iorb, fmt = "(I0)"))), string)
          else
             call set(vald // ("Orbital" // trim(yaml_toa(iorb, fmt = "(I0)"))), string)
          end if
       end if
    end do

    close(unit = 91)

    if (ntd > 0) then
       call set(dict // key // "K point 1" // "up", valu)
       call set(dict // key // "K point 1" // "down", vald)
    else
       call set(dict // key // "K point 1", valu)
    end if

    call set(dict // key // SOURCE_KEY, filename)

  end subroutine occupation_data_file_merge_to_dict

  subroutine neb_set_from_dict(dict, opt, climbing_, imax, nimg_, &
       & cv, tol, ds_, kmin, kmax, temp_, damp_, meth)
    use module_defs, only: gp
    !use module_input_keys
    use public_keys
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(out) :: opt, climbing_
    integer, intent(out) :: imax, nimg_
    real(kind = gp), intent(out) :: cv, tol, ds_, damp_, kmin, kmax, temp_
    character(len = max_field_length), intent(out) :: meth

    if (.not. has_key(dict, GEOPT_VARIABLES)) return
    if (trim(dict_value(dict // GEOPT_VARIABLES // GEOPT_METHOD)) /= "NEB") return

    opt       = dict // GEOPT_VARIABLES // EXTREMA_OPT
    climbing_ = dict // GEOPT_VARIABLES // NEB_CLIMBING
    imax      = dict // GEOPT_VARIABLES // NCOUNT_CLUSTER_X
    nimg_     = dict // GEOPT_VARIABLES // NIMG
    cv        = dict // GEOPT_VARIABLES // FORCEMAX
    tol       = dict // GEOPT_VARIABLES // FIX_TOL
    ds_       = dict // GEOPT_VARIABLES // BETAX
    kmin      = dict // GEOPT_VARIABLES // SPRINGS_K // 0
    kmax      = dict // GEOPT_VARIABLES // SPRINGS_K // 1
    meth      = dict // GEOPT_VARIABLES // NEB_METHOD
    if (has_key(dict // GEOPT_VARIABLES, TEMP)) temp_ = dict // GEOPT_VARIABLES // TEMP
    if (has_key(dict // GEOPT_VARIABLES, NEB_DAMP)) damp_ = dict // GEOPT_VARIABLES // NEB_DAMP
  end subroutine neb_set_from_dict
end module module_input_dicts

!> Module reading the old format (before 1.7) for the input
module input_old_text_format
  use yaml_strings, only: operator(.eqv.)
  use public_enums
  implicit none
  private
  integer, parameter :: WF_N_FORMAT      = 4
  character(len = 12), dimension(0:WF_N_FORMAT-1), parameter :: wf_format_names = &
       (/ "none        ", &
       "plain text  ", &
       "Fortran bin.", &
       "ETSF        " /)


  !> All possible values of input psi (determination of the input guess)
  integer, dimension(12), parameter :: input_psi_values = &
       (/ INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_CP2K, &
       INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, INPUT_PSI_DISK_WVL, &
       INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS, &
       INPUT_PSI_LINEAR_AO, INPUT_PSI_DISK_LINEAR, INPUT_PSI_MEMORY_LINEAR /)


  public :: read_input_dict_from_files!input_from_old_text_format
contains

  function input_psi_names(id)
    use public_enums
    integer, intent(in) :: id
    character(len = 14) :: input_psi_names

    select case(id)
    case(INPUT_PSI_EMPTY)
       write(input_psi_names, "(A)") "empty"
    case(INPUT_PSI_RANDOM)
       write(input_psi_names, "(A)") "random"
    case(INPUT_PSI_CP2K)
       write(input_psi_names, "(A)") "CP2K"
    case(INPUT_PSI_LCAO)
       write(input_psi_names, "(A)") "LCAO"
    case(INPUT_PSI_MEMORY_WVL)
       write(input_psi_names, "(A)") "wvl. in mem."
    case(INPUT_PSI_DISK_WVL)
       write(input_psi_names, "(A)") "wvl. on disk"
    case(INPUT_PSI_LCAO_GAUSS)
       write(input_psi_names, "(A)") "LCAO + gauss."
    case(INPUT_PSI_MEMORY_GAUSS)
       write(input_psi_names, "(A)") "gauss. in mem."
    case(INPUT_PSI_DISK_GAUSS)
       write(input_psi_names, "(A)") "gauss. on disk"
    case(INPUT_PSI_LINEAR_AO)
       write(input_psi_names, "(A)") "Linear AO"
    case(INPUT_PSI_MEMORY_LINEAR)
       write(input_psi_names, "(A)") "Linear restart"
    case(INPUT_PSI_DISK_LINEAR)
       write(input_psi_names, "(A)") "Linear on disk"
    case default
       write(input_psi_names, "(A)") "Error"
    end select
  end function input_psi_names


  subroutine input_psi_help()
    integer :: i

    write(*, "(1x,A)") "Available values of inputPsiId are:"
    do i = 1, size(input_psi_values)
       write(*, "(1x,A,I5,A,A)") " | ", input_psi_values(i), &
            & " - ", input_psi_names(input_psi_values(i))
    end do
  end subroutine input_psi_help

  subroutine output_wf_format_help()
    integer :: i

    write(*, "(1x,A)") "Available values of output_wf are:"
    do i = 0, size(wf_format_names) - 1
       write(*, "(1x,A,I5,A,A)") " | ", i, &
            & " - ", wf_format_names(i)
    end do
  end subroutine output_wf_format_help

  !> This function returns a dictionary with all the input variables of a BigDFT run filled.
  !! This dictionary is constructed from a updated version of the input variables dictionary
  !! following the input files as defined  by the user
  subroutine read_input_dict_from_files(radical,mpi_env,dict)
    use module_base
    use wrapper_MPI
    use module_input_dicts, only: merge_input_file_to_dict
    !use yaml_output
    use dynamic_memory
    implicit none
    character(len = *), intent(in) :: radical    !< The name of the run. use "input" if empty
    type(mpi_environment), intent(in) :: mpi_env !< The environment where the variables have to be updated
    type(dictionary), pointer :: dict            !< Input dictionary, has to be nullified at input
    !local variables
    integer :: ierr
    logical :: exists_default, exists_user
    character(len = max_field_length) :: fname

    call f_routine(id='read_input_dict_from_files')


    ! We try first default.yaml
    inquire(file = "default.yaml", exist = exists_default)
    if (exists_default) call merge_input_file_to_dict(dict, "default.yaml", mpi_env)

    ! We try then radical.yaml
    if (len_trim(radical) == 0) then
       fname(1:len(fname)) = "input.yaml"
    else
       fname(1:len(fname)) = trim(radical) // ".yaml"
    end if
    inquire(file = trim(fname), exist = exists_user)
    if (exists_user) call merge_input_file_to_dict(dict, trim(fname), mpi_env)

    ! We fallback on the old text format (to be eliminated in the future)
    if (.not.exists_default .and. .not. exists_user) then
       call input_from_old_text_format(radical,mpi_env,dict)
    end if

    ! We put a barrier here to be sure that non master proc will be stopped
    ! by any issue on the master proc.
    call mpibarrier(comm=mpi_env%mpi_comm)

    call f_release_routine()
  end subroutine read_input_dict_from_files

  subroutine input_from_old_text_format(radical,mpi_env,dict)
    use public_keys
    use wrapper_MPI
    use dictionaries
    use yaml_output
    implicit none
    character(len = *), intent(in) :: radical    !< The name of the run. use "input" if empty
    type(mpi_environment), intent(in) :: mpi_env !< The environment where the variables have to be updated
    type(dictionary), pointer :: dict            !< Input dictionary
    !local variables
    character(len = 100) :: f0
    type(dictionary), pointer :: vals
    
    ! Parse all files.
    call set_inputfile(f0, radical, PERF_VARIABLES)
    nullify(vals)
    call read_perf_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//PERF_VARIABLES, vals)

    call set_inputfile(f0, radical, DFT_VARIABLES)
    nullify(vals)
    call read_dft_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//DFT_VARIABLES, vals)

    call set_inputfile(f0, radical, KPT_VARIABLES)
    nullify(vals)
    call read_kpt_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//KPT_VARIABLES, vals)

    call set_inputfile(f0, radical, GEOPT_VARIABLES)
    nullify(vals)
    call read_geopt_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//GEOPT_VARIABLES, vals)

    call set_inputfile(f0, radical, MIX_VARIABLES)
    nullify(vals)
    call read_mix_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//MIX_VARIABLES, vals)

    call set_inputfile(f0, radical, SIC_VARIABLES)
    nullify(vals)
    call read_sic_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//SIC_VARIABLES, vals)

    call set_inputfile(f0, radical, TDDFT_VARIABLES)
    nullify(vals)
    call read_tddft_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//TDDFT_VARIABLES, vals)

    !call set_inputfile(f0, radical, 'lin')
    !call read_lin_and_frag_from_text_format(mpi_env%iproc,dict,trim(radical)) !as it also reads fragment

    call set_inputfile(f0, radical, 'neb')
    nullify(vals)
    call read_neb_from_text_format(mpi_env%iproc,vals, trim(f0))
    if (associated(vals)) call set(dict//GEOPT_VARIABLES, vals)

    if (mpi_env%iproc==0) then
       call yaml_warning('Input files read in the old format.'//&
            'Use the input_minimal.yaml file to switch to new format. '//&
            'In future versions this will be deprecated')
    end if

  end subroutine input_from_old_text_format

  !> Set and check the input file
  !! if radical is empty verify if the file input.ext exists. 
  !! otherwise search for radical.ext
  !! if the so defined file is not existing, then filename becomes default.ext
  subroutine set_inputfile(filename, radical, ext)
    implicit none
    character(len = *), intent(in) :: radical, ext
    character(len = 100), intent(out) :: filename

    logical :: exists

    write(filename, "(A)") ""
    if (trim(radical) == "") then
       write(filename, "(A,A,A)") "input", ".", trim(ext)
    else
       write(filename, "(A,A,A)") trim(radical), ".", trim(ext)
    end if

    inquire(file=trim(filename),exist=exists)
    if (.not. exists .and. (trim(radical) /= "input" .and. trim(radical) /= "")) &
         & write(filename, "(A,A,A)") "default", ".", trim(ext)
  end subroutine set_inputfile

  subroutine read_dft_from_text_format(iproc,dict,filename)
    use module_base
    use module_input
    use public_keys
    use dictionaries
    !  use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iproc
    !local variables
    logical :: exists
    integer :: ierror
    real(gp), dimension(2), parameter :: hgrid_rng=(/0.0_gp,2.0_gp/)
    real(gp), dimension(2), parameter :: xrmult_rng=(/0.0_gp,100.0_gp/)

    logical :: dummy_bool
    integer :: dummy_int
    real(gp) :: dummy_real
    real(gp), dimension(3) :: dummy_real3

    !dft parameters, needed for the SCF part
    call input_set_file(iproc,(iproc == 0),trim(filename),exists, DFT_VARIABLES)
    !if (exists) in%files = in%files + INPUTS_DFT
    !call the variable, its default value, the line ends if there is a comment
    if (.not. exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    !grid spacings
    call input_var(dummy_real3(1),'0.45',dict//HGRIDS//0,ranges=hgrid_rng)
    call input_var(dummy_real3(2),'0.45',dict//HGRIDS//1,ranges=hgrid_rng)
    call input_var(dummy_real3(3),'0.45',dict//HGRIDS//2,ranges=hgrid_rng,&
         comment='hx,hy,hz: grid spacing in the three directions')

    !coarse and fine radii around atoms
    call input_var(dummy_real,'5.0',dict//RMULT//0,ranges=xrmult_rng)
    call input_var(dummy_real,'8.0',dict//RMULT//1,ranges=xrmult_rng,&
         comment='c(f)rmult: c(f)rmult*radii_cf(:,1(2))=coarse(fine) atom-based radius')

    !XC functional (ABINIT XC codes)
    call input_var(dummy_int,'1',dict//IXC,comment='ixc: exchange-correlation parameter (LDA=1,PBE=11)')

    !charge and electric field
    call input_var(dummy_int,'0',dict//NCHARGE,ranges=(/-500,500/))

    call input_var(dummy_real3(1),'0.',dict//ELECFIELD//0)
    call input_var(dummy_real3(2),'0.',dict//ELECFIELD//1)
    call input_var(dummy_real3(3),'0.',dict//ELECFIELD//2,&
         comment='charge of the system, Electric field (Ex,Ey,Ez)')

    !spin and polarization
    call input_var(dummy_int,'1',dict//NSPIN,exclusive=(/1,2,4/))
    call input_var(dummy_int,'0',dict//MPOL,comment='nspin=1 non-spin polarization, mpol=total magnetic moment')

    !convergence parameters
    call input_var(dummy_real,'1.e-4',dict//GNRM_CV,ranges=(/1.e-20_gp,1.0_gp/),&
         comment='gnrm_cv: convergence criterion gradient')
    call input_var(dummy_int,'50',dict//ITERMAX,ranges=(/0,10000/))
    call input_var(dummy_int,'1',dict//NREPMAX,ranges=(/0,1000/),&
         comment='itermax,nrepmax: max. # of wfn. opt. steps and of re-diag. runs')

    !convergence parameters
    call input_var(dummy_int,'6',dict//NCONG,ranges=(/0,20/))
    call input_var(dummy_int,'6',dict//IDSX,ranges=(/0,15/),&
         comment='ncong, idsx: # of CG it. for preconditioning eq., wfn. diis history')

    !dispersion parameter
    call input_var(dummy_int,'0',dict//DISPERSION,ranges=(/0,5/),&
         comment='dispersion correction potential (values 1,2,3,4,5), 0=none')

    ! Now the variables which are to be used only for the last run
    call input_var(dummy_int,'0',dict//INPUTPSIID,&
         exclusive=(/-2,-1,0,2,10,12,13,100,101,102/),input_iostat=ierror)
    ! Validate inputPsiId value (Can be added via error handling exception)
    if (ierror /=0 .and. iproc == 0) then
       write( *,'(1x,a,I0,a)')'ERROR: illegal value of inputPsiId (', dummy_int, ').'
       call input_psi_help()
       call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierror)
    end if

    call input_var(dummy_int,'0',dict//OUTPUT_WF,exclusive=(/0,1,2,3/),input_iostat=ierror)
    ! Validate output_wf value.
    if (ierror /=0 .and. iproc == 0) then
       write( *,'(1x,a,I0,a)')'ERROR: illegal value of output_wf (', dummy_int, ').'
       call output_wf_format_help()
       call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierror)
    end if

    call input_var(dummy_int,'0',dict//OUTPUT_DENSPOT,exclusive=(/0,1,2,10,11,12,20,21,22/),&
         comment='InputPsiId, output_wf, output_denspot')

    ! Tail treatment.
    call input_var(dummy_real,'0.0',dict//RBUF,ranges=(/0.0_gp,10.0_gp/))
    call input_var(dummy_int,'30',dict//NCONGT,ranges=(/1,50/),&
         comment='rbuf, ncongt: length of the tail (AU),# tail CG iterations')

    !davidson treatment
    call input_var(dummy_int,'0',dict//NORBV,ranges=(/-9999,9999/))
    call input_var(dummy_int,'0',dict//NVIRT)
    call input_var(dummy_int,'0',dict//NPLOT,&
         comment='Davidson subspace dim., # of opt. orbs, # of plotted orbs')

    ! Line to disable automatic behaviours (currently only symmetries).
    call input_var(dummy_bool,'F',dict//DISABLE_SYM,comment='disable the symmetry detection')

    call input_free(.false.)

  end subroutine read_dft_from_text_format


  !> Read the input variables needed for the geometry optmization
  !! Every argument should be considered as mandatory
  subroutine read_geopt_from_text_format(iproc,dict,filename)
    use module_base
    use module_input
    use public_keys
    use dictionaries
    implicit none
    integer, intent(in) :: iproc
    character(len=*), intent(in) :: filename
    type(dictionary), pointer :: dict
    !local variables
    character(len=*), parameter :: subname='read_geopt_from_text_format'
    integer :: i
    logical :: exists

    character(len = 5) :: dummy_str
    integer :: dummy_int, ionmov_
    real(gp) :: dummy_real

    !geometry input parameters
    call input_set_file(iproc,(iproc == 0),trim(filename),exists,GEOPT_VARIABLES)  
    !if (exists) in%files = in%files + INPUTS_GEOPT
    !call the variable, its default value, the line ends if there is a comment
    if (.not. exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    call input_var(dummy_str,"BFGS",dict // GEOPT_METHOD, comment = "")
    !call set(dict // GEOPT_METHOD, dummy_str)
    call input_var(dummy_int,'1',dict // NCOUNT_CLUSTER_X,comment="")
    !call set(dict // NCOUNT_CLUSTER_X, dummy_int)

    call input_var(dummy_real,'1.0',dict // FRAC_FLUCT)
    !call set(dict // FRAC_FLUCT, dummy_real, fmt = "(E8.2)")
    call input_var(dummy_real,'0.0',dict // FORCEMAX,comment="")
    !call set(dict // FORCEMAX, dummy_real, fmt = "(E8.2)")
    call input_var(dummy_real,'0.0',dict // RANDDIS,comment="")
    !call set(dict // RANDDIS, dummy_real, fmt = "(E8.2)")

    if (trim(dummy_str) .eqv. "AB6MD") then
       call input_var(ionmov_,'6',dict // IONMOV,comment="")
       !call set(dict // IONMOV, ionmov_)
       call input_var(dummy_real,'20.670689',dict // DTION,comment="")
       !call set(dict // DTION, dummy_real)
       if (ionmov_ == 6) then
          call input_var(dummy_real,'300',dict // MDITEMP,comment="")
          !call set(dict // MDITEMP, dummy_real)
       elseif (ionmov_ > 7) then
          call input_var(dummy_real,'300',dict // MDITEMP)
          !call set(dict // MDITEMP, dummy_real)
          call input_var(dummy_real,'300',dict // MDFTEMP,comment="")
          !call set(dict // MDFTEMP, dummy_real)
       end if

       if (ionmov_ == 8) then
          call input_var(dummy_real,'1.e5',dict // NOSEINERT,comment="")
          !call set(dict // NOSEINERT, dummy_real)
       else if (ionmov_ == 9) then
          call input_var(dummy_real,'1.e-3',dict // FRICTION,comment="")
          !call set(dict // FRICTION, dummy_real)
          call input_var(dummy_real,'1.e4',dict // MDWALL,comment="")
          !call set(dict // MDWALL, dummy_real)
       else if (ionmov_ == 13) then
          !here no dictionary
          call input_var(dummy_int,'0',ranges=(/0,100/),comment="")
          do i=1,dummy_int-1
             call input_var(dummy_real,'0.0',dict // QMASS // (i-1))
             !call set(dict // QMASS // (i-1), dummy_real)
          end do
          if (dummy_int > 0) then
             call input_var(dummy_real,'0.0',dict // QMASS // (dummy_int-1),comment="")
             !call set(dict // QMASS // (dummy_int-1), dummy_real)
          end if
          call input_var(dummy_real,'10',dict // BMASS)
          !call set(dict // BMASS, dummy_real)
          call input_var(dummy_real,'1.0',dict // VMASS,comment="")
          !call set(dict // VMASS, dummy_real)
       end if
    else if (trim(dummy_str) .eqv. "DIIS") then
       call input_var(dummy_real,'2.0',dict // BETAX)
       !call set(dict // BETAX, dummy_real, fmt = "(F6.3)")
       call input_var(dummy_int,'4',dict // HISTORY,comment="")
       !call set(dict // HISTORY, dummy_int)
    else
       call input_var(dummy_real,'4.0',dict // BETAX,comment="")
       !call set(dict // BETAX, dummy_real, fmt = "(F6.3)")
    end if
    if (trim(dummy_str) .eqv. "FIRE") then
       call input_var(dummy_real,'0.75',dict // DTINIT)
       !call set(dict // DTINIT, dummy_real, fmt = "(F6.3)")
       call input_var(dummy_real, '1.5',dict // DTMAX,comment="")
       !call set(dict // DTMAX, dummy_real, fmt = "(F6.3)")
    endif

    call input_free(.false.)

  END SUBROUTINE read_geopt_from_text_format

  !> Read the input variables needed for the geometry optmization
  !!    Every argument should be considered as mandatory
  subroutine read_mix_from_text_format(iproc,dict,filename)
    use module_base
    use module_input
    use public_keys
    use dictionaries
    implicit none
    !Arguments
    integer, intent(in) :: iproc
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: filename
    !local variables
    !n(c) character(len=*), parameter :: subname='mix_input_variables'
    logical :: exists
    integer :: dummy_int
    real(gp) :: dummy_real

    !Mix parameters, needed for the SCF poart with Davidson
    call input_set_file(iproc,(iproc == 0),trim(filename),exists,MIX_VARIABLES)
    !if (exists) in%files = in%files + INPUTS_MIX
    !call the variable, its default value, the line ends if there is a comment
    if (.not.exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    !Controls the self-consistency: 0 direct minimisation otherwise ABINIT convention
    call input_var(dummy_int,'0',dict // ISCF,comment="")
    !call set(dict // ISCF, dummy_int)
    call input_var(dummy_int,'1',dict // ITRPMAX,comment="")
    !call set(dict // ITRPMAX, dummy_int)
    call input_var(dummy_real,'1.e-4',dict // RPNRM_CV,comment="")
    !call set(dict // RPNRM_CV, dummy_real, fmt = "(E8.1)")
    call input_var(dummy_int,'0',dict // NORBSEMPTY)
    !call set(dict // NORBSEMPTY, dummy_int)
    call input_var(dummy_real,'0.0',dict // TEL) 
    !call set(dict // TEL, dummy_real, fmt = "(E9.2)")
    call input_var(dummy_int,'1',dict // OCCOPT,comment="")
    !call set(dict // OCCOPT, dummy_int)
    call input_var(dummy_real,'0.0',dict // ALPHAMIX)
    !call set(dict // ALPHAMIX, dummy_real, fmt = "(F6.3)")
    call input_var(dummy_real,'2.0',dict // ALPHADIIS,comment="")
    !call set(dict // ALPHADIIS, dummy_real, fmt = "(F6.3)")

    call input_free(.false.)
  END SUBROUTINE read_mix_from_text_format

  !> Read Self-Interaction Correction (SIC) input parameters
  subroutine read_sic_from_text_format(iproc,dict,filename)
    use module_input
    use public_keys
    use dictionaries
    implicit none
    integer, intent(in) :: iproc
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: filename
    !local variables
    logical :: exists
    !n(c) character(len=*), parameter :: subname='sic_input_variables'
    double precision :: dummy_real
    character(len = 4) :: dummy_str

    !Self-Interaction Correction input parameters
    call input_set_file(iproc,(iproc == 0),trim(filename),exists,'SIC Parameters')  
    !if (exists) in%files = in%files + INPUTS_SIC
    if (.not.exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    call input_var(dummy_str,'NONE',dict // SIC_APPROACH,comment='')
    !call set(dict // SIC_APPROACH, dummy_str)
    call input_var(dummy_real,'0.0',dict // SIC_ALPHA,comment='')
    !call set(dict // SIC_ALPHA, dummy_real, fmt = "(E8.2)")
    
    if (trim(dummy_str) .eqv. 'NK') then
       call input_var(dummy_real,'0.0',dict // SIC_FREF,comment='')
       !call set(dict // SIC_FREF, dummy_real, fmt = "(E8.2)")
    end if

    call input_free(.false.)
  END SUBROUTINE read_sic_from_text_format

  subroutine read_tddft_from_text_format(iproc,dict,filename)
    use module_input
    use public_keys
    use dictionaries
    implicit none
    integer, intent(in) :: iproc
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: filename
    !local variables
    logical :: exists
    !n(c) character(len=*), parameter :: subname='tddft_input_variables'
    character(len = 4) :: dummy_str

    !TD-DFT parameters
    call input_set_file(iproc,(iproc == 0),trim(filename),exists,'TD-DFT Parameters')  
    !if (exists) in%files = in%files + INPUTS_TDDFT
    !call the variable, its default value, the line ends if there is a comment
    if (.not. exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    call input_var(dummy_str,"NONE",dict // TDDFT_APPROACH,comment="")
    !call set(dict // TDDFT_APPROACH, dummy_str)

    call input_free(.false.)

  END SUBROUTINE read_tddft_from_text_format

  subroutine read_kpt_from_text_format(iproc,dict,filename)
    use module_base
    use dictionaries
    use module_input
    use public_keys
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iproc
    type(dictionary), pointer :: dict
    !local variables
    logical :: exists
    character(len=*), parameter :: subname='read_kpt_from_text_format'

    integer :: dummy_int, nseg, i, ierror
    integer, dimension(3) :: dummy_int3
    real(gp) :: dummy_real
    real(gp), dimension(3) :: dummy_real3
    character(len = max_field_length) :: dummy_str

    !kpt parameters, needed for the SCF part
    call input_set_file(iproc,(iproc == 0),trim(filename),exists, KPT_VARIABLES)
    !if (exists) in%files = in%files + INPUTS_KPT
    !call the variable, its default value, the line ends if there is a comment
    if (.not. exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    !if the file does exist, we fill up the dictionary.
    call input_var(dummy_str, 'manual',dict//KPT_METHOD, comment='K-point sampling method')
    !call set(dict//KPT_METHOD, trim(dummy_str))

    if (trim(dummy_str) .eqv.  'auto') then
       call input_var(dummy_real,'0.0',dict//KPTRLEN, comment='Equivalent length of K-space resolution (Bohr)')
       !call set(dict//KPTRLEN, dummy_real)
    else if (trim(dummy_str) .eqv. 'mpgrid') then
       !take the points of Monckorst-pack grid
       call input_var(dummy_int3(1),'1',dict//NGKPT//0)
       call input_var(dummy_int3(2),'1',dict//NGKPT//1)
       call input_var(dummy_int3(3),'1',dict//NGKPT//2, comment='No. of Monkhorst-Pack grid points')
       !call set(dict//NGKPT//0, dummy_int3(1))
       !call set(dict//NGKPT//1, dummy_int3(2))
       !call set(dict//NGKPT//2, dummy_int3(3))
       !shift
       !no dict here
       call input_var(dummy_int,'1',ranges=(/1,8/),comment='No. of different shifts')
       !read the shifts
       do i=1,dummy_int
          call input_var(dummy_real3(1),'0.',dict//SHIFTK//(i-1)//0)
          call input_var(dummy_real3(2),'0.',dict//SHIFTK//(i-1)//1)
          call input_var(dummy_real3(3),'0.',dict//SHIFTK//(i-1)//2,comment=' ')
          !call set(dict//SHIFTK//(i-1)//0, dummy_real3(1), fmt = "(F6.4)")
          !call set(dict//SHIFTK//(i-1)//1, dummy_real3(2), fmt = "(F6.4)")
          !call set(dict//SHIFTK//(i-1)//2, dummy_real3(3), fmt = "(F6.4)")
       end do
    else if (trim(dummy_str) .eqv. 'manual') then
       call input_var(dummy_int,'1',ranges=(/1,10000/),&
            comment='Number of K-points')
       do i=1,dummy_int
          call input_var(dummy_real3(1),'0.',dict//KPT//(i-1)//0)
          call input_var(dummy_real3(2),'0.',dict//KPT//(i-1)//1)
          call input_var(dummy_real3(3),'0.',dict//KPT//(i-1)//2)
          !call set(dict//KPT//(i-1)//0, dummy_real3(1), fmt = "(F6.4)")
          !call set(dict//KPT//(i-1)//1, dummy_real3(2), fmt = "(F6.4)")
          !call set(dict//KPT//(i-1)//2, dummy_real3(3), fmt = "(F6.4)")
          call input_var(dummy_real,'1.',dict//WKPT//(i-1),comment='K-pt coords, K-pt weigth')
          !call set(dict//WKPT//(i-1), dummy_real, fmt = "(F6.4)")
       end do
    end if

    ! Now read the band structure definition. do it only if the file exists
    !no dictionary here
    call input_var(dummy_str,'bands',comment='For doing band structure calculation',&
         input_iostat=ierror)
    call set(dict//BANDS, (ierror==0))
    if (ierror==0) then
       call input_var(nseg,'1',ranges=(/1,1000/),&
            comment='# of segments of the BZ path')
       !number of points for each segment, parallel granularity
       do i=1,nseg
          call input_var(dummy_int,'1',dict//ISEG)
          !call set(dict//ISEG, dummy_int)
       end do
       call input_var(dummy_int,'1',dict//NGRANULARITY,&
            comment='points for each segment, # of points done for each group')
       !call set(dict//NGRANULARITY, dummy_int)

       call input_var(dummy_real3(1),'0.',dict//KPTV//0//0)
       call input_var(dummy_real3(2),'0.',dict//KPTV//0//1)
       call input_var(dummy_real3(3),'0.',dict//KPTV//0//2,comment=' ')
!       call set(dict//KPTV//0//0, dummy_real3(1))
!       call set(dict//KPTV//0//1, dummy_real3(2))
!       call set(dict//KPTV//0//2, dummy_real3(3))
       do i=1,nseg
          call input_var(dummy_real3(1),'0.5',dict//KPTV//(i-1)//0)
          call input_var(dummy_real3(2),'0.5',dict//KPTV//(i-1)//1)
          call input_var(dummy_real3(3),'0.5',dict//KPTV//(i-1)//2,comment=' ')
          !call set(dict//KPTV//(i-1)//0, dummy_real3(1))
          !call set(dict//KPTV//(i-1)//1, dummy_real3(2))
          !call set(dict//KPTV//(i-1)//2, dummy_real3(3))
       end do

       !read an optional line to see if there is a file associated
       !no dict for the moment
       call input_var(dummy_str,' ',&
            comment=' ',input_iostat=ierror)
       if (ierror == 0) then
          !since a file for the local potential is already given, do not perform ground state calculation
          call set(dict//BAND_STRUCTURE_FILENAME, dummy_str)
       end if
    end if

    call input_free(.false.)

  end subroutine read_kpt_from_text_format


  !> Read the input variables which can be used for performances
  subroutine read_perf_from_text_format(iproc,dict,filename)
    use module_input
    use dictionaries
    use public_keys
    implicit none
    character(len=*), intent(in) :: filename
    type(dictionary), pointer :: dict
    integer, intent(in) :: iproc
    !local variables
    !n(c) character(len=*), parameter :: subname='perf_input_variables'
    logical :: exists, dummy_bool
    integer :: dummy_int, blocks(2)
    double precision :: dummy_real
    character(len = 7) :: dummy_str

    call input_set_file(iproc, (iproc == 0), filename, exists, PERF_VARIABLES)
    !if (exists) in%files = in%files + INPUTS_PERF
    if (.not. exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    call input_var("debug", .false., "Debug option", dummy_bool)
    call set(dict // DEBUG, dummy_bool)
    call input_var("fftcache", 8*1024, "Cache size for the FFT", dummy_int)
    call set(dict // FFTCACHE, dummy_int)
    call input_var("accel", "NO", "Acceleration", dummy_str)
    call set(dict // ACCEL, dummy_str)

    !determine desired OCL platform which is used for acceleration
    call input_var("OCL_platform"," ", "Chosen OCL platform", dummy_str)
    call set(dict // OCL_PLATFORM, dummy_str)
    call input_var("OCL_devices"," ", "Chosen OCL devices", dummy_str)
    call set(dict // OCL_DEVICES, dummy_str)

    !!@TODO to relocate
    call input_var("blas", .false., "CUBLAS acceleration", dummy_bool)
    call set(dict // BLAS, dummy_bool)
    call input_var("projrad", 15.0d0, "Radius ", dummy_real)
    call set(dict // PROJRAD, dummy_real, fmt = "(F6.3)")
    call input_var("exctxpar", "OP2P", "Exact exchange parallelisation scheme", dummy_str)
    call set(dict // EXCTXPAR, dummy_str)
    call input_var("ig_diag", .true.,"Input guess", dummy_bool)
    call set(dict // IG_DIAG, dummy_bool)
    call input_var("ig_norbp", 5, "Input guess: ", dummy_int)
    call set(dict // IG_NORBP, dummy_int)
    call input_var("ig_blocks", (/ 300, 800 /), "Input guess: ", blocks)
    call set(dict // IG_BLOCKS // 0, blocks(1))
    call set(dict // IG_BLOCKS // 1, blocks(2))
    call input_var("ig_tol", 1d-4, "Input guess: Tolerance criterion", dummy_real)
    call set(dict // IG_TOL, dummy_real, fmt = "(E8.1)")
    call input_var("methortho", 0, "Orthogonalisation ", dummy_int)
    call set(dict // METHORTHO, dummy_int)
    call input_var("rho_commun", "DEF","Density communication scheme (DBL, RSC, MIX)",dummy_str)
    call set(dict // RHO_COMMUN, dummy_str)
    call input_var("psolver_groupsize",0, "Size of ", dummy_int)
    call set(dict // PSOLVER_GROUPSIZE, dummy_int)
    call input_var("psolver_accel",0, "Acceleration ", dummy_int)
    call set(dict // PSOLVER_ACCEL, dummy_int)
    call input_var("unblock_comms", "OFF", "Overlap Com)",dummy_str)
    call set(dict // UNBLOCK_COMMS, dummy_str)
    call input_var("linear", 'OFF', "Linear Input Guess approach",dummy_str)
    call set(dict // LINEAR, dummy_str)
    call input_var("tolsym", 1d-8, "Tolerance for symmetry detection",dummy_real)
    call set(dict // TOLSYM, dummy_real, fmt = "(E8.1)")
    call input_var("signaling", .false., "Expose calculation results on Network",dummy_bool)
    call set(dict // SIGNALING, dummy_bool)
    call input_var("signalTimeout", 0, "Time out on startup for signal connection",dummy_int)  
    call set(dict // SIGNALTIMEOUT, dummy_int)
    call input_var("domain", "", "Domain to add to the hostname to find the IP", dummy_str)
    call set(dict // DOMAIN, dummy_str)
    call input_var("inguess_geopt", 0,"0= wavelet input ",dummy_int)
    call set(dict // INGUESS_GEOPT, dummy_int)
    call input_var("store_index", .true., "Linear scaling: store ", dummy_bool)
    call set(dict // STORE_INDEX, dummy_bool)
    !verbosity of the output
    call input_var("verbosity", 2, "Verbosity of the output 0=low, 2=high",dummy_int)
    call set(dict // VERBOSITY, dummy_int)

    !If false, apply the projectors in the once-and-for-all scheme, otherwise on-the-fly
    call input_var("psp_onfly", .true., "Calculate the PSP projectors on the fly (less memory)",dummy_bool)
    call set(dict // PSP_ONFLY, dummy_bool)

    !If true, preserve the multipole of the ionic part (local potential) projecting on delta instead of ISF
    call input_var("multipole_preserving", .false., "Preserve multipole moment of the ionic charge",dummy_bool)
    call set(dict // MULTIPOLE_PRESERVING, dummy_bool)
    call input_var("mp_isf", 16, "Interpolating scaling function for the multipole preserving option",dummy_int)
    call set(dict // MP_ISF, dummy_int)

    !block size for pdsyev/pdsygv, pdgemm (negative -> sequential)
    call input_var("pdsyev_blocksize",-8,"SCALAPACK linear scaling blocksize",dummy_int) !ranges=(/-100,1000/)
    call set(dict // PDSYEV_BLOCKSIZE, dummy_int)
    call input_var("pdgemm_blocksize",-8,"SCALAPACK linear scaling blocksize",dummy_int) !ranges=(/-100,1000/)
    call set(dict // PDGEMM_BLOCKSIZE, dummy_int)

    !max number of process uses for pdsyev/pdsygv, pdgemm
    call input_var("maxproc_pdsyev",4,"SCALAPACK linear scaling max num procs",dummy_int) !ranges=(/1,100000/)
    call set(dict // MAXPROC_PDSYEV, dummy_int)
    call input_var("maxproc_pdgemm",4,"SCALAPACK linear scaling max num procs",dummy_int) !ranges=(/1,100000/)
    call set(dict // MAXPROC_PDGEMM, dummy_int)

    !FOE: if the determinant of the interpolation matrix to find the Fermi energy
    !is smaller than this value, switch from cubic to linear interpolation.
    call input_var("ef_interpol_det",1.d-20,"FOE: max ",dummy_real)
    call set(dict // EF_INTERPOL_DET, dummy_real, fmt = "(E9.2)")
    call input_var("ef_interpol_chargediff",10.d0,"FOE: max ",dummy_real)
    call set(dict // EF_INTERPOL_CHARGEDIFF, dummy_real, fmt = "(E9.2)")

    !determines whether a mixing step shall be preformed after the input guess !(linear version)
    call input_var("mixing_after_inputguess",1,"mixing after inguess (0/1/2)",dummy_int)
    call set(dict // MIXING_AFTER_INPUTGUESS, dummy_int)

    !determines whether the input guess support functions are orthogonalized iteratively (T) or in the standard way (F)
    call input_var("iterative_orthogonalization",.false.," orbitals",dummy_bool)
    call set(dict // ITERATIVE_ORTHOGONALIZATION, dummy_bool)

    call input_var("check_sumrho", 2, (/0,1,2/), "linear sumrho: 0=no check, 1=light check, 2=full check", dummy_int)
    call set(dict // CHECK_SUMRHO, dummy_int)

    call input_var("check_overlap", 2, (/0,1,2/), "linear overlap: 0=no check, 1=light check, 2=full check", dummy_int)
    call set(dict // CHECK_OVERLAP, dummy_int)

    call input_var("experimental_mode", .false., "linear scaling: activate the experimental mode", dummy_bool)
    call set(dict // EXPERIMENTAL_MODE, dummy_bool)

    call input_var("write_orbitals", 0, "(LS): write KS orbitals for cubic restart (0: no, 1: wvl, 2: wvl+isf)", dummy_int)
    call set(dict // WRITE_ORBITALS, dummy_int)

    call input_var("explicit_locregcenters", .false., "linear scaling: explicitely specify localization centers", dummy_bool)
    call set(dict // EXPLICIT_LOCREGCENTERS, dummy_bool)

    call input_var("calculate_KS_residue", .true., "linear scaling: calculate Kohn-Sham residue", dummy_bool)
    call set(dict // CALCULATE_KS_RESIDUE, dummy_bool)

    call input_var("intermediate_forces", .false., "linear scaling: calculate intermediate forces", dummy_bool)
    call set(dict // INTERMEDIATE_FORCES, dummy_bool)

    call input_var("kappa_conv", 0.1d0, "exit kappa for extended input guess (experimental mode)", dummy_real)
    call set(dict // KAPPA_CONV, dummy_real)

    call input_var("evbounds_nsatur", 3, "number of FOE cycles before the eigenvalue bounds are shrinked", dummy_int)
    call set(dict // EVBOUNDS_NSATUR, dummy_int)

    call input_var("evboundsshrink_nsatur", 4, "maximal number of unsuccessful eigenvalue bounds shrinkings", dummy_int)
    call set(dict // EVBOUNDSSHRINK_NSATUR, dummy_int)

    call input_var("method_updatekernel", 0, (/0,1,2/), "K update (sup fun opt) (0: purific., 1: FOE, 2: renorm.)", dummy_int)
    call set(dict // METHOD_UPDATEKERNEL, dummy_int)

    call input_var("purification_quickreturn", .true., "linear scaling: quick return in purification", dummy_bool)
    call set(dict // PURIFICATION_QUICKRETURN, dummy_bool)

    call input_var("adjust_FOE_temperature", .true., "dynamic adjustment of FOE error function decay length", dummy_bool)
    call set(dict // ADJUST_FOE_TEMPERATURE, dummy_bool)

    call input_var("calculate_gap", .false., "calculate the HOMO LUMO gap", dummy_bool)
    call set(dict // CALCULATE_GAP, dummy_bool)

    call input_var("loewdin_charge_analysis", .false., "perform a Loewdin charge analysis at the end", dummy_bool)
    call set(dict // LOEWDIN_CHARGE_ANALYSIS, dummy_bool)

    call input_var("check_matrix_compression", .true., "perform a check of the matrix compression routines", dummy_bool)
    call set(dict // CHECK_MATRIX_COMPRESSION, dummy_bool)

    call input_var("correction_co_contra", .true., "correction covariant / contravariant gradient", dummy_bool)
    call set(dict // CORRECTION_CO_CONTRA, dummy_bool)

    call input_var("fscale_lowerbound", 5.d-3, "lower bound for the error function decay length", dummy_real)
    call set(dict // FSCALE_LOWERBOUND, dummy_real)

    call input_var("fscale_upperbound", 5.d-2, "upper bound for the error function decay length", dummy_real)
    call set(dict // FSCALE_UPPERBOUND, dummy_real)

    call input_var("imethod_overlap", 1, (/1,2/), "lin scaling method to calculate overlap matrix (1:old, 2:new)", dummy_int)
    call set(dict // IMETHOD_OVERLAP, dummy_int)

    call input_var("enable_matrix_taskgroups", .true., "enable matrix taskgroups", dummy_bool)
    call set(dict // ENABLE_MATRIX_TASKGROUPS, dummy_bool)

    call input_var("hamapp_radius_incr", 8, "radius enlargement for Ham application", dummy_int)
    call set(dict // HAMAPP_RADIUS_INCR, dummy_int)

    call input_var("adjust_kernel_iterations", .true., "addaptive ajustment of the number of kernel iterations", dummy_bool)
    call set(dict // ADJUST_KERNEL_ITERATIONS, dummy_bool)

    call input_var("wf_extent_analysis", .false., "extent analysis of the support functions / KS orbitals", dummy_bool)
    call set(dict // WF_EXTENT_ANALYSIS, dummy_bool)

    call input_free(.false.)

  END SUBROUTINE read_perf_from_text_format


  !> Read the linear input variables
  subroutine read_lin_and_frag_from_text_format(iproc,dict,run_name)
    use module_base
    use module_input
    use public_keys
    implicit none
    character(len=*), intent(in) :: run_name
    type(dictionary), pointer :: dict
    integer, intent(in) :: iproc
    !local variables
    !n(c) character(len=*), parameter :: subname='perf_input_variables'
    logical :: exists, dummy_bool,frag_bool
    integer :: dummy_int,ios
    double precision :: dummy_real
    character(len=256) :: comments,dummy_char,filename
    type(dictionary), pointer :: dict_basis

    !call f_err_throw('For the linear version the input parameters must be read in the .yaml format, &
    !    &the old version is deprecated', err_name='BIGDFT_INPUT_VARIABLES_ERROR')

    filename=repeat(' ',len(filename))
    call set_inputfile(filename, trim(run_name),    "lin")
    ! This name seems to be too long..
    !call input_set_file(iproc,(iproc == 0),trim(filename),exists,'Parameters for Localized basis generation (O(N) approach)')
    call input_set_file(iproc,(iproc == 0),trim(filename),exists,'Parameters for O(N) approach')
!    call input_set_file(iproc, (iproc == 0), filename, exists, LIN_GENERAL)
!    call input_set_file(iproc, (iproc == 0), filename, exists, LIN_BASIS)
!    call input_set_file(iproc, (iproc == 0), filename, exists, LIN_KERNEL)
    !if (exists) in%files = in%files + INPUTS_PERF
    if (.not. exists) then
       call input_free(.false.)
       return
    end if

    if (.not. associated(dict)) call dict_init(dict)

    ! General variables #######################################################

    comments='number of accuracy levels: either 2 (for low/high accuracy) or 1 (for hybrid mode)'
    call input_var(dummy_int,'2',ranges=(/1,2/),comment=comments)
    call dict_set(dict//LIN_GENERAL//HYBRID,dummy_int==1)

    ! number of iterations
    comments = 'outer loop iterations (low, high)'
    call input_var(dummy_int,'15',dict//LIN_GENERAL//NIT//0,ranges=(/0,100000/))
    call input_var(dummy_int,'1',dict//LIN_GENERAL//NIT//1,ranges=(/0,100000/),comment=comments)

    comments = 'basis iterations (low, high)'
    call input_var(dummy_int,'12',dict//LIN_BASIS//NIT//0,ranges=(/0,100000/))
    call input_var(dummy_int,'50',dict//LIN_BASIS//NIT//1,ranges=(/0,100000/),comment=comments)

    comments = 'kernel iterations (low, high) - directmin only'
    call input_var(dummy_int,'1',dict//LIN_KERNEL//NSTEP//0,ranges=(/0,1000/))
    call input_var(dummy_int,'1',dict//LIN_KERNEL//NSTEP//1,ranges=(/0,1000/),comment=comments)

    comments = 'density iterations (low, high)'
    call input_var(dummy_int,'15',dict//LIN_KERNEL//NIT//0,ranges=(/0,1000/))
    call input_var(dummy_int,'15',dict//LIN_KERNEL//NIT//1,ranges=(/0,1000/),comment=comments)

    ! DIIS history lengths
    comments = 'DIIS history for basis (low, high)'
    call input_var(dummy_int,'5',dict//LIN_BASIS//IDSX//0,ranges=(/0,100/))
    call input_var(dummy_int,'0',dict//LIN_BASIS//IDSX//1,ranges=(/0,100/),comment=comments)

    comments = 'DIIS history for kernel (low, high) - directmin only'
    call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX_COEFF//0,ranges=(/0,100/))
    call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX_COEFF//1,ranges=(/0,100/),comment=comments)

    comments = 'DIIS history for density mixing (low, high)'
    call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX//0,ranges=(/0,100/))
    call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX//1,ranges=(/0,100/),comment=comments)

    ! mixing parameters
    comments = 'density mixing parameter (low, high)'
    call input_var(dummy_real,'.5d0',dict//LIN_KERNEL//ALPHAMIX//0,ranges=(/0.d0,1.d0/))
    call input_var(dummy_real,'.5d0',dict//LIN_KERNEL//ALPHAMIX//1,ranges=(/0.d0,1.d0/),comment=comments)

    ! Convergence criteria
    comments = 'outer loop convergence (low, high)'
    call input_var(dummy_real,'1.d-8' ,dict//LIN_GENERAL//RPNRM_CV//0,ranges=(/0.d0,1.d0/))
    call input_var(dummy_real,'1.d-12',dict//LIN_GENERAL//RPNRM_CV//1,ranges=(/0.d0,1.d0/),comment=comments)

    comments = 'basis convergence (low, high) ; early stop TMB optimization, dynamic gnrm, activate dyn (exp. mode only)'
    call input_var(dummy_real,'1.d-3',dict//LIN_BASIS//GNRM_CV//0,ranges=(/0.0_gp,1.0_gp/))
    call input_var(dummy_real,'1.d-5',dict//LIN_BASIS//GNRM_CV//1,ranges=(/0.0_gp,1.0_gp/))
    call input_var(dummy_real,'1.d-4',dict//LIN_BASIS//DELTAE_CV,ranges=(/0.0_gp,1.0_gp/))
    call input_var(dummy_real,'1.d-4',dict//LIN_BASIS//GNRM_DYN,ranges=(/0.0_gp,1.0_gp/))
    call input_var(dummy_real,'1.d-3',dict//LIN_BASIS//MIN_GNRM_FOR_DYNAMIC,ranges=(/1.d-7,1.0_gp/),comment=comments)

    comments = 'factor to reduce the confinement. Only used for hybrid mode.'
    call input_var(dummy_real,'0.5d0',dict//LIN_GENERAL//CONF_DAMPING,ranges=(/-1.d100,1.d0/),comment=comments)

    comments = 'kernel convergence (low, high) - directmin only'
    call input_var(dummy_real,'75.d-5',dict//LIN_KERNEL//GNRM_CV_COEFF//0,ranges=(/0.d0,1.d0/))
    call input_var(dummy_real,'1.d-5',dict//LIN_KERNEL//GNRM_CV_COEFF//1,ranges=(/0.d0,1.d0/),comment=comments)

    comments = 'density convergence (low, high)'
    call input_var(dummy_real,'1.d-13',dict//LIN_KERNEL//RPNRM_CV//0,ranges=(/0.d0,1.d0/))
    call input_var(dummy_real,'1.d-13',dict//LIN_KERNEL//RPNRM_CV//1,ranges=(/0.d0,1.d0/),comment=comments)

    comments = 'convergence criterion on density to fix TMBS'
    call input_var(dummy_real,'1.d-10',dict//LIN_BASIS//FIX_BASIS,ranges=(/1.d-14,1.d-6/),comment=comments)
    !call input_var(in%lin%support_functions_converged,'1.d-10',ranges=(/0.d0,1.d0/),comment=comments)

    comments='mixing method: 100 (direct minimization), 101 (simple dens mixing), 102 (simple pot mixing), 103 (FOE)'
    call input_var(dummy_int,'100',ranges=(/100,103/),comment=comments)
    select case(dummy_int)
    case(100)
       call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'DIRMIN')
    case(101) 
       call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'DIAG')
       call dict_set(dict//LIN_KERNEL//MIXING_METHOD,'DEN')
    case(102)      
       call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'DIAG')
       call dict_set(dict//LIN_KERNEL//MIXING_METHOD,'POT')
    case(103)
       call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'FOE')
    end select

    comments = 'initial step size for basis optimization (DIIS, SD)' ! DELETE ONE
    call input_var(dummy_real,'1.d0',dict//LIN_BASIS//ALPHA_DIIS,ranges=(/0.0_gp,10.0_gp/))
    call input_var(dummy_real,'1.d0',dict//LIN_BASIS//ALPHA_SD,ranges=(/0.0_gp,10.0_gp/),comment=comments)

    comments = 'initial step size for kernel update (SD), curve fitting for alpha update - directmin only'
    call input_var(dummy_real,'1.d0',dict//LIN_KERNEL//ALPHA_SD_COEFF,ranges=(/0.0_gp,10.0_gp/))
    call input_var(dummy_bool,'F',dict//LIN_KERNEL//ALPHA_FIT_COEFF,comment=comments)

    comments = 'lower and upper bound for the eigenvalue spectrum (FOE). Will be adjusted automatically if chosen too small'
    call input_var(dummy_real,'-.5d0',dict//LIN_KERNEL//EVAL_RANGE_FOE//0,ranges=(/-10.d0,-1.d-10/))
    call input_var(dummy_real,'-.5d0',dict//LIN_KERNEL//EVAL_RANGE_FOE//1,ranges=(/1.d-10,10.d0/),comment=comments)

    !comments='number of iterations in the preconditioner, order of Taylor approximations'
    comments='number of iterations in the preconditioner'
    call input_var(dummy_int,'5',dict//LIN_BASIS//NSTEP_PREC,ranges=(/1,100/),comment=comments)
    !!call input_var(dummy_int,'1',dict//LIN_GENERAL//TAYLOR_ORDER,ranges=(/-100,100/),comment=comments)
    !call input_var(in%lin%order_taylor,'1',ranges=(/1,100/),comment=comments)

    comments = '0-> exact Loewdin, 1-> taylor expansion; &
               &in orthoconstraint: correction for non-orthogonality (0) or no correction (1)'
    call input_var(dummy_int,'1',dict//LIN_GENERAL//TAYLOR_ORDER,ranges=(/-100,10000/))
    call input_var(dummy_int,'1',dict//LIN_BASIS//CORRECTION_ORTHOCONSTRAINT,comment=comments)
    !call input_var(in%lin%correctionOrthoconstraint,'1',ranges=(/0,1/),comment=comments)

    comments='fscale: length scale over which complementary error function decays from 1 to 0'
    call input_var(dummy_real,'1.d-2',dict//LIN_KERNEL//FSCALE_FOE,ranges=(/0.d0,1.d0/),comment=comments)

    !plot basis functions: true or false
    comments='Output support functions (i: wvl, i+10: wvl+isf): i=0 No, i=1 formatted, i=2 Fortran bin, i=3 ETSF ;'//&
             'calculate dipole ; pulay correction (old and new); diagonalization at the end (dmin, FOE)'
    call input_var(dummy_int,'0',dict//LIN_GENERAL//OUTPUT_WF,exclusive=(/0,1,2,3,10,11,12,13/))
    call input_var(dummy_bool,'F',dict//LIN_GENERAL//CALC_DIPOLE)
    call input_var(dummy_bool,'T',dict//LIN_GENERAL//CALC_PULAY//0)
    call input_var(dummy_bool,'F',dict//LIN_GENERAL//CALC_PULAY//1)

!    in%lin%pulay_correction=dummy_bool
!    call input_var(in%lin%new_pulay_correction,'F')
    call input_var(dummy_bool,'F',dict//LIN_GENERAL//SUBSPACE_DIAG,comment=comments)

  !fragment calculation and transfer integrals: true or false
  comments='fragment calculation; calculate transfer_integrals; constrained DFT calculation; extra states to optimize (dmin only)'
  !these should becode dummy variables to build dictionary
  !!call input_var(in%lin%fragment_calculation,'F')
  !!call input_var(in%lin%calc_transfer_integrals,'F')
  !!call input_var(in%lin%constrained_dft,'F')
  !!call input_var(in%lin%extra_states,'0',ranges=(/0,10000/),comment=comments)
  call input_var(frag_bool,'F')
  !this variable makes sense only if fragments are specified
  call input_var(dummy_bool,'F')
  if (frag_bool) call dict_set(dict//FRAG_VARIABLES//TRANSFER_INTEGRALS,dummy_bool)
  
  call input_var(dummy_bool,'F') !constrained DFT, obtained via the charges
  call input_var(dummy_int,'0',dict//LIN_GENERAL//EXTRA_STATES,&
       ranges=(/0,10000/),comment=comments)

  ! Now read in the parameters specific for each atom type.
  comments = 'Atom name, number of basis functions per atom, prefactor for confinement potential,'//&
       'localization radius, kernel cutoff, kernel cutoff FOE'
  read_basis: do !while(itype <= atoms%astruct%ntypes) 
  !!   if (exists) then
        call input_var(dummy_char,'C',input_iostat=ios)
        if (ios /= 0) exit read_basis
     dict_basis=>dict//LIN_BASIS_PARAMS//trim(dummy_char)
     call input_var(dummy_int,'1',dict_basis//NBASIS,ranges=(/1,100/))
     call input_var(dummy_real,'1.2d-2',dict_basis//AO_CONFINEMENT,&
          ranges=(/-1.0_gp,1.0_gp/))
     call input_var(dummy_real,'1.2d-2',dict_basis//CONFINEMENT//0,&
          ranges=(/-1.0_gp,1.0_gp/))
     call input_var(dummy_real,'5.d-5',dict_basis//CONFINEMENT//1,&
          ranges=(/-1.0_gp,1.0_gp/))
     call input_var(dummy_real,'10.d0',dict_basis//RLOC//0,&
          ranges=(/1.0_gp,10000.0_gp/))
     call input_var(dummy_real,'10.d0',dict_basis//RLOC//1,&
          ranges=(/1.0_gp,10000.0_gp/))
     call input_var(dummy_real,'12.d0',dict_basis//RLOC_KERNEL,&
          ranges=(/1.0_gp,10000.0_gp/))
     call input_var(dummy_real,'20.d0',dict_basis//RLOC_KERNEL_FOE,&
          ranges=(/1.0_gp,10000.0_gp/),comment=comments)

  !!   if (.not. exists) exit read_basis !default has been filled
  end do read_basis

    call input_free(.false.)

    !read extensively the file and build the temporary variable
    !from which the input dictionary is updated
    if (frag_bool) then
       filename=repeat(' ',len(filename))
       call set_inputfile(filename, run_name,   "frag")
       call fragment_input_variables_from_text_format(iproc,.false.,&
            trim(filename),frag_bool,dict//FRAG_VARIABLES)
    end if

  END SUBROUTINE read_lin_and_frag_from_text_format


  subroutine read_neb_from_text_format(iproc,dict,filename)
    use module_base
    use module_input
    use public_keys
    use dictionaries
    implicit none
    character(len=*), intent(in) :: filename
    type(dictionary), pointer :: dict
    integer, intent(in) :: iproc

    INTEGER :: num_of_images
    CHARACTER (LEN=20) :: minimization_scheme
    logical :: climbing, optimization, restart, exists
    integer :: max_iterations
    real(gp) :: convergence, damp, k_min, k_max, ds, temp_req, tolerance
    CHARACTER (LEN=80) :: first_config, last_config, job_name, scratch_dir

    NAMELIST /NEB/ scratch_dir,         &
         climbing,            &
         optimization,        &
         minimization_scheme, &
         damp,                &
         temp_req,            &
         k_max, k_min,        &
         ds,                  &
         max_iterations,      &
         tolerance,           &
         convergence,         &
         num_of_images,       &
         restart,             & ! not used
         job_name,            & ! not used
         first_config,        & ! not used
         last_config            ! not used

    inquire(file=trim(filename),exist=exists)
    if (.not. exists) return

    open(unit = 123, file = trim(filename), action = "read")
    READ(123 , NML=NEB )
    close(123)

    if (.not. associated(dict)) call dict_init(dict)

    call set(dict // GEOPT_METHOD, "NEB")
    call set(dict // NEB_CLIMBING, climbing)
    call set(dict // EXTREMA_OPT, optimization)
    call set(dict // NEB_METHOD, minimization_scheme)
    if (trim(minimization_scheme) == 'damped-verlet') call set(dict // NEB_DAMP, damp)
    call set(dict // SPRINGS_K // 0, k_min)
    call set(dict // SPRINGS_K // 1, k_max)
    if (trim(minimization_scheme) == 'sim-annealing') call set(dict // TEMP, temp_req)
    call set(dict // BETAX, ds)
    call set(dict // NCOUNT_CLUSTER_X, max_iterations)
    call set(dict // FIX_TOL, tolerance)
    call set(dict // FORCEMAX, convergence)
    call set(dict // NIMG, num_of_images)

  end subroutine read_neb_from_text_format

  !> Read fragment input parameters
  subroutine fragment_input_variables_from_text_format(iproc,dump,filename,shouldexist,dict)
    use module_base
    use fragment_base
    use module_input
    use yaml_output, only: yaml_toa,yaml_map
    implicit none
    logical, intent(in) :: shouldexist
    integer, intent(in) :: iproc
    character(len=*), intent(in) :: filename
    type(dictionary), pointer :: dict
    logical, intent(in) :: dump
    !local variables
    !character(len=*), parameter :: subname='fragment_input_variables'
    logical :: exists
    character(len=256) :: comments
    integer :: ifrag, frag_num
    real(gp) :: charge
    type(fragmentInputParameters) :: frag
    type(dictionary), pointer :: dict_frag

    !Linear input parameters
    call input_set_file(iproc,dump,trim(filename),exists,'Fragment Parameters') 

    if (.not. exists .and. shouldexist) then ! we should be doing a fragment calculation, so this is a problem
       call f_err_throw("The file 'input.frag' is missing and fragment calculation was specified",&
            err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       call input_free(.false.)
       return
    end if

    call nullifyInputFragParameters(frag)

    !example of interpreted fragment yaml file
    !!frag:
    !!  transfer_integrals: Yes
    !!  frag_name1: [1, ... , 3, 5, ... , 8]
    !!  frag_name2: [9, 10, 13, 16, ... ,18]
    !!  frag_name3: [11, 12, 15]
    !!  constrained_dft:
    !!    Fragment No. 9: +1
    !!    Fragment No. 15: -1
    !!

    ! number of reference fragments
    comments='# number of fragments in reference system, number of fragments in current system'
    call input_var(frag%nfrag_ref,'1',ranges=(/1,100000/))
    call input_var(frag%nfrag,'1',ranges=(/1,100000/),comment=comments)

    ! Allocate fragment pointers
    call allocateInputFragArrays(frag)

    ! ADD A SENSIBLE DEFAULT AND ALLOW FOR USER NOT TO SPECIFY FRAGMENT NAMES
    comments = '#  reference fragment number i, fragment label'
    do ifrag=1,frag%nfrag_ref
       call input_var(frag_num,'1',ranges=(/1,frag%nfrag_ref/))
       if (frag_num/=ifrag) then
          call f_err_throw("The file 'input.frag'  has an error when specifying"//&
               " the reference fragments",err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       end if
       call input_var(frag%label(frag_num),' ',comment=comments)
       frag%label(frag_num)=trim(frag%label(frag_num))
       ! keep dirname blank if this isn't a fragment calculation
       if (len(trim(frag%label(frag_num)))>=1) then
          frag%dirname(frag_num)='data-'//trim(frag%label(frag_num))//'/'
       else
          frag%dirname(frag_num)=''
       end if
    end do

    comments = '# fragment number j, reference fragment i this corresponds to, charge on this fragment'
    do ifrag=1,frag%nfrag
       call input_var(frag_num,'1',ranges=(/1,frag%nfrag/))
       if (frag_num/=ifrag) then
          call f_err_throw("The file 'input.frag'  has an error when specifying"//&
               " the system fragments",err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       end if
       call input_var(frag%frag_index(frag_num),'1',ranges=(/0,100000/))
       call input_var(charge,'0.d0',ranges=(/-500.d0,500.d0/),comment=comments)
       frag%charge(frag_num)=charge
       !call input_var(frag%charge(frag_num),'1',ranges=(/-500,500/),comment=comments)
    end do

    call input_free(.false.)

    call dict_from_frag(frag,dict_frag)

    !call yaml_map('Fragment dictionary',dict_frag)

    call dict_update(dict,dict_frag)
    call dict_free(dict_frag)
    call deallocateInputFragArrays(frag)

  END SUBROUTINE fragment_input_variables_from_text_format

end module input_old_text_format
