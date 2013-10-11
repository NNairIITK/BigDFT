!> @file
!!  Routines to read and print input variables
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> this function returns a dictionary with all the input variables of a BigDFT run filled
!! this dictionary is constructed from a updated version of the input variables dictionary
!! following the input files as defined  by the user
function read_input_dict_from_files(radical, mpi_env) result(dict)
  use dictionaries
  use wrapper_MPI
  use module_input_keys
  use module_interfaces, only: merge_input_file_to_dict
  use input_old_text_format
  implicit none
  character(len = *), intent(in) :: radical !< the name of the run. use "input" if empty
  type(mpi_environment), intent(in) :: mpi_env !< the environment where the variables have to be updated

  integer :: ierr
  type(dictionary), pointer :: dict
  logical :: exists_default, exists_user
  character(len = max_field_length) :: fname
  character(len = 100) :: f0

  ! Handle error with master proc only.
  if (mpi_env%iproc > 0) call f_err_set_callback(f_err_ignore)
  
  nullify(dict)
  ! We try first default.yaml
  inquire(file = "default.yaml", exist = exists_default)
  if (exists_default) call merge_input_file_to_dict(dict, "default.yaml", mpi_env)

  ! We try then radical.yaml
  if (len_trim(radical) == 0) then
     fname = "input.yaml"
  else
     fname(1:max_field_length) = trim(radical) // ".yaml"
  end if
  inquire(file = trim(fname), exist = exists_user)
  if (exists_user) call merge_input_file_to_dict(dict, trim(fname), mpi_env)

  ! We fallback on the old text format (to be eliminated in the future)
  if (.not.exists_default .and. .not. exists_user) then
     ! Parse all files.
     call dict_init(dict)
     call set_inputfile(f0, radical, PERF_VARIABLES)
     call read_perf_from_text_format(mpi_env%iproc,dict//PERF_VARIABLES, trim(f0))
     call set_inputfile(f0, radical, DFT_VARIABLES)
     call read_dft_from_text_format(mpi_env%iproc,dict//DFT_VARIABLES, trim(f0))
     call set_inputfile(f0, radical, KPT_VARIABLES)
     call read_kpt_from_text_format(mpi_env%iproc,dict//KPT_VARIABLES, trim(f0))
     call set_inputfile(f0, radical, GEOPT_VARIABLES)
     call read_geopt_from_text_format(mpi_env%iproc,dict//GEOPT_VARIABLES, trim(f0))
     call set_inputfile(f0, radical, MIX_VARIABLES)
     call read_mix_from_text_format(mpi_env%iproc,dict//MIX_VARIABLES, trim(f0))
     call set_inputfile(f0, radical, SIC_VARIABLES)
     call read_sic_from_text_format(mpi_env%iproc,dict//SIC_VARIABLES, trim(f0))
     call set_inputfile(f0, radical, TDDFT_VARIABLES)
     call read_tddft_from_text_format(mpi_env%iproc,dict//TDDFT_VARIABLES, trim(f0))
  else
     ! We add an overloading input.perf (for automatic test purposes).
     ! This will be changed in far future when only YAML input will be allowed.
     call set_inputfile(f0, radical, PERF_VARIABLES)
     call read_perf_from_text_format(mpi_env%iproc,dict//PERF_VARIABLES, trim(f0))
  end if

  if (mpi_env%iproc > 0) call f_err_severe_restore()

  ! We put a barrier here to be sure that non master proc will be stop
  ! by any issue on the master proc.
  call mpi_barrier(mpi_env%mpi_comm, ierr)
end function read_input_dict_from_files

!> Routine to read YAML input files and create input dictionary.
subroutine merge_input_file_to_dict(dict, fname, mpi_env)
  use module_input_keys
  use dictionaries
  use yaml_parse
  use wrapper_MPI
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: fname
  type(mpi_environment), intent(in) :: mpi_env

  integer(kind = 8) :: cbuf, cbuf_len
  integer :: ierr
  character(len = max_field_length) :: val
  character, dimension(:), allocatable :: fbuf
  type(dictionary), pointer :: udict
  
  if (mpi_env%iproc == 0) then
     call getFileContent(cbuf, cbuf_len, fname, len_trim(fname))
     if (mpi_env%nproc > 1) &
          & call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
  else
     call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
  end if
  allocate(fbuf(cbuf_len))
  if (mpi_env%iproc == 0) then
     call copyCBuffer(fbuf(1), cbuf, cbuf_len)
     call freeCBuffer(cbuf)
     if (mpi_env%nproc > 1) &
          & call mpi_bcast(fbuf(1), cbuf_len, MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
  else
     call mpi_bcast(fbuf(1), cbuf_len, MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
  end if

  call f_err_open_try()
  call yaml_parse_from_char_array(udict, fbuf)

  ! Handle with possible partial dictionary.
  deallocate(fbuf)
  call dict_update(dict, udict // 0)
  call dict_free(udict)
  
  ierr = 0
  if (f_err_check()) ierr = f_get_last_error(val)
  call f_err_close_try()
  !in the present implementation f_err_check is not cleaned after the close of the try
  if (ierr /= 0) call f_err_throw(err_id = ierr, err_msg = val)
end subroutine merge_input_file_to_dict

!> Fill the input_variables structure with the information
!! contained in the dictionary dict
!! the dictionary should be completes to fill all the information
subroutine inputs_from_dict(in, atoms, dict, dump)
  use module_types
  use module_defs
  use yaml_output
  use module_interfaces, except => inputs_from_dict
  use dictionaries
  use module_input_keys
  implicit none
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  type(dictionary), pointer :: dict
  logical, intent(in) :: dump

  !type(dictionary), pointer :: profs
  integer :: ierr

  call default_input_variables(in)

  ! To avoid race conditions where procs create the default file and other test its
  ! presence, we put a barrier here.
  if (bigdft_mpi%nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)

  ! Analyse the input dictionary and transfer it to in.
  call input_keys_fill_all(dict)
  call perf_input_analyse(bigdft_mpi%iproc, in, dict//PERF_VARIABLES)
  call dft_input_analyse(bigdft_mpi%iproc, in, dict//DFT_VARIABLES)
  call geopt_input_analyse(bigdft_mpi%iproc, in, dict//GEOPT_VARIABLES)
  call mix_input_analyse(bigdft_mpi%iproc, in, dict//MIX_VARIABLES)
  call sic_input_analyse(bigdft_mpi%iproc, in, dict//SIC_VARIABLES, in%ixc)
  call tddft_input_analyse(bigdft_mpi%iproc, in, dict//TDDFT_VARIABLES)

  ! Shake atoms, if required.
  call astruct_set_displacement(atoms%astruct, in%randdis)
  if (bigdft_mpi%nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
  ! Update atoms with symmetry information
  call astruct_set_symmetries(atoms%astruct, in%disableSym, in%symTol, in%elecfield)

  call kpt_input_analyse(bigdft_mpi%iproc, in, dict//KPT_VARIABLES, &
       & atoms%astruct%sym, atoms%astruct%geocode, atoms%astruct%cell_dim)
  
  if (bigdft_mpi%iproc == 0 .and. dump) call input_keys_dump(dict)

  if (in%gen_nkpt > 1 .and. in%gaussian_help) then
     if (bigdft_mpi%iproc==0) call yaml_warning('Gaussian projection is not implemented with k-point support')
     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
  end if
  
  if(in%inputpsiid==100 .or. in%inputpsiid==101 .or. in%inputpsiid==102) &
      DistProjApply=.true.
  if(in%linear /= INPUT_IG_OFF .and. in%linear /= INPUT_IG_LIG) then
     !only on the fly calculation
     DistProjApply=.true.
  end if

  ! Stop the code if it is trying to run GPU with spin=4
  if (in%nspin == 4 .and. (GPUconv .or. OCLconv)) then
     if (bigdft_mpi%iproc==0) call yaml_warning('GPU calculation not implemented with non-collinear spin')
     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
  end if

  ! This should be moved to init_atomic_values, but we can't since
  ! input.lin needs some arrays allocated here.
  if (.not. associated(atoms%nzatom)) then
     call allocate_atoms_nat(atoms, "inputs_from_dict")
     call allocate_atoms_ntypes(atoms, "inputs_from_dict")
  end if

  ! Linear scaling (if given)
  !in%lin%fragment_calculation=.false. ! to make sure that if we're not doing a linear calculation we don't read fragment information
  call lin_input_variables_new(bigdft_mpi%iproc,dump .and. (in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR), trim(in%file_lin),in,atoms)

  ! Fragment information (if given)
  call fragment_input_variables(bigdft_mpi%iproc,dump .and. (in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR).and.in%lin%fragment_calculation,trim(in%file_frag),in,atoms)

!!$  ! Stop code for unproper input variables combination.
!!$  if (in%ncount_cluster_x > 0 .and. .not. in%disableSym .and. atoms%geocode == 'S') then
!!$     if (bigdft_mpi%iproc==0) then
!!$        write(*,'(1x,a)') 'Change "F" into "T" in the last line of "input.dft"'   
!!$        write(*,'(1x,a)') 'Forces are not implemented with symmetry support, disable symmetry please (T)'
!!$     end if
!!$     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
!!$  end if

!!$  if (bigdft_mpi%iproc == 0) then
!!$     profs => input_keys_get_profiles("")
!!$     call yaml_dict_dump(profs)
!!$     call dict_free(profs)
!!$  end if

  !check whether a directory name should be associated for the data storage
  call check_for_data_writing_directory(bigdft_mpi%iproc,in)
  
  !check if an error has been found and raise an exception to be handled
  if (f_err_check()) then
     call f_err_throw('Error in reading input variables from dictionary',&
          err_name='BIGDFT_INPUT_VARIABLES_ERROR')
  end if

end subroutine inputs_from_dict

!> Check the directory of data (create if not present)
subroutine check_for_data_writing_directory(iproc,in)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(inout) :: in
  !local variables
  logical :: shouldwrite

  if (iproc==0) call yaml_comment('|',hfill='-')

  !initialize directory name
  shouldwrite=.false.

  shouldwrite=shouldwrite .or. &
       in%output_wf_format /= WF_FORMAT_NONE .or. &    !write wavefunctions
       in%output_denspot /= output_denspot_NONE .or. & !write output density
       in%ncount_cluster_x > 1 .or. &                  !write posouts or posmds
       in%inputPsiId == 2 .or. &                       !have wavefunctions to read
       in%inputPsiId == 12 .or.  &                     !read in gaussian basis
       in%gaussian_help .or. &                         !Mulliken and local density of states
       in%writing_directory /= '.' .or. &              !have an explicit local output directory
       bigdft_mpi%ngroup > 1   .or. &                  !taskgroups have been inserted
       in%lin%plotBasisFunctions > 0 .or. &            !dumping of basis functions for locreg runs
       in%inputPsiId == 102                            !reading of basis functions

  !here you can check whether the etsf format is compiled

  if (shouldwrite) then
     call create_dir_output(iproc, in)
     if (iproc==0) call yaml_map('Data Writing directory',trim(in%dir_output))
  else
     if (iproc==0) call yaml_map('Data Writing directory','./')
     in%dir_output=repeat(' ',len(in%dir_output))
  end if
END SUBROUTINE check_for_data_writing_directory

subroutine create_dir_output(iproc, in)
  use yaml_output
  use module_types
  use module_base
  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(inout) :: in

  character(len=100) :: dirname
  integer :: i_stat,ierror,ierr

  ! Create a directory to put the files in.
  dirname=repeat(' ',len(dirname))
  if (iproc == 0) then
     call getdir(in%dir_output, len_trim(in%dir_output), dirname, 100, i_stat)
     if (i_stat /= 0) then
        call yaml_warning("Cannot create output directory '" // trim(in%dir_output) // "'.")
        call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
     end if
  end if
  call MPI_BCAST(dirname,len(dirname),MPI_CHARACTER,0,bigdft_mpi%mpi_comm,ierr)
  in%dir_output=dirname
END SUBROUTINE create_dir_output

!> Set default values for input variables
subroutine default_input_variables(in)
  use module_base
  use module_types
  use dictionaries
  implicit none

  type(input_variables), intent(inout) :: in

  ! Default values.
  in%output_wf_format = WF_FORMAT_NONE
  in%output_denspot_format = output_denspot_FORMAT_CUBE
  nullify(in%gen_kpt)
  nullify(in%gen_wkpt)
  nullify(in%kptv)
  nullify(in%nkptsv_group)
  ! Default abscalc variables
  call abscalc_input_variables_default(in)
  ! Default frequencies variables
  call frequencies_input_variables_default(in)
  ! Default values for geopt.
  call geopt_input_variables_default(in) 
  ! Default values for mixing procedure
  call mix_input_variables_default(in) 
  ! Default values for tddft
  call tddft_input_variables_default(in)
  !Default for Self-Interaction Correction variables
  call sic_input_variables_default(in)
  ! Default for signaling
  in%gmainloop = 0.d0
  ! Default for lin.
  nullify(in%lin%potentialPrefac_lowaccuracy)
  nullify(in%lin%potentialPrefac_highaccuracy)
  nullify(in%lin%potentialPrefac_ao)
  nullify(in%lin%norbsPerType)
  nullify(in%lin%locrad)
  nullify(in%lin%locrad_lowaccuracy)
  nullify(in%lin%locrad_highaccuracy)
  nullify(in%lin%locrad_type)
  nullify(in%lin%kernel_cutoff)
  !nullify(in%frag%frag_info)
  nullify(in%frag%label)
  nullify(in%frag%dirname)
  nullify(in%frag%frag_index)
  nullify(in%frag%charge)
END SUBROUTINE default_input_variables

!> Assign default values for mixing variables
subroutine mix_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in

  !mixing treatement (hard-coded values)
  in%iscf=0
  in%itrpmax=1
  in%alphamix=0.0_gp
  in%rpnrm_cv=1.e-4_gp
  in%gnrm_startmix=0.0_gp
  in%norbsempty=0
  in%Tel=0.0_gp
  in%occopt=SMEARING_DIST_ERF
  in%alphadiis=2.d0

END SUBROUTINE mix_input_variables_default

!> Assign default values for GEOPT variables
subroutine geopt_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in

  !put some fake values for the geometry optimsation case
  in%geopt_approach='SDCG'
  in%ncount_cluster_x=0
  in%frac_fluct=1.0_gp
  in%forcemax=0.0_gp
  in%randdis=0.0_gp
  in%betax=2.0_gp
  in%history = 1
  in%wfn_history = 1
  in%ionmov = -1
  in%dtion = 0.0_gp
  in%strtarget(:)=0.0_gp
  nullify(in%qmass)

END SUBROUTINE geopt_input_variables_default

!> Assign default values for self-interaction correction variables
subroutine sic_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in

  in%SIC%approach='NONE'
  in%SIC%alpha=0.0_gp
  in%SIC%fref=0.0_gp

END SUBROUTINE sic_input_variables_default




!> Assign default values for TDDFT variables
subroutine tddft_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in

  in%tddft_approach='NONE'

END SUBROUTINE tddft_input_variables_default

  subroutine allocateInputFragArrays(input_frag)
    use module_types
    implicit none
  
    ! Calling arguments
    type(fragmentInputParameters),intent(inout) :: input_frag
  
    ! Local variables
    integer :: i_stat
    character(len=*),parameter :: subname='allocateInputFragArrays'
 
    allocate(input_frag%frag_index(input_frag%nfrag), stat=i_stat)
    call memocc(i_stat, input_frag%frag_index, 'input_frag%frag_index', subname)

    allocate(input_frag%charge(input_frag%nfrag), stat=i_stat)
    call memocc(i_stat, input_frag%charge, 'input_frag%charge', subname)

    !allocate(input_frag%frag_info(input_frag%nfrag_ref,2), stat=i_stat)
    !call memocc(i_stat, input_frag%frag_info, 'input_frag%frag_info', subname)

    allocate(input_frag%label(input_frag%nfrag_ref), stat=i_stat)
    call memocc(i_stat, input_frag%label, 'input_frag%label', subname)

    allocate(input_frag%dirname(input_frag%nfrag_ref), stat=i_stat)
    call memocc(i_stat, input_frag%dirname, 'input_frag%dirname', subname)

  end subroutine allocateInputFragArrays

  subroutine deallocateInputFragArrays(input_frag)
    use module_types
    implicit none
  
    ! Calling arguments
    type(fragmentInputParameters),intent(inout) :: input_frag
  
    ! Local variables
    integer :: i_stat,i_all
    character(len=*),parameter :: subname='deallocateInputFragArrays'
 
    !if(associated(input_frag%frag_info)) then
    !  i_all = -product(shape(input_frag%frag_info))*kind(input_frag%frag_info)
    !  deallocate(input_frag%frag_info,stat=i_stat)
    !  call memocc(i_stat,i_all,'input_frag%frag_info',subname)
    !  nullify(input_frag%frag_info)
    !end if 
 
    if(associated(input_frag%frag_index)) then
      i_all = -product(shape(input_frag%frag_index))*kind(input_frag%frag_index)
      deallocate(input_frag%frag_index,stat=i_stat)
      call memocc(i_stat,i_all,'input_frag%frag_index',subname)
      nullify(input_frag%frag_index)
    end if 
 
    if(associated(input_frag%charge)) then
      i_all = -product(shape(input_frag%charge))*kind(input_frag%charge)
      deallocate(input_frag%charge,stat=i_stat)
      call memocc(i_stat,i_all,'input_frag%charge',subname)
      nullify(input_frag%charge)
    end if 

    if(associated(input_frag%label)) then
      i_all = -product(shape(input_frag%label))*kind(input_frag%label)
      deallocate(input_frag%label,stat=i_stat)
      call memocc(i_stat,i_all,'input_frag%label',subname)
      nullify(input_frag%label)
    end if 

    if(associated(input_frag%dirname)) then
      i_all = -product(shape(input_frag%dirname))*kind(input_frag%dirname)
      deallocate(input_frag%dirname,stat=i_stat)
      call memocc(i_stat,i_all,'input_frag%dirname',subname)
      nullify(input_frag%dirname)
    end if 

  end subroutine deallocateInputFragArrays


  subroutine nullifyInputFragParameters(input_frag)
    use module_types
    implicit none

    ! Calling arguments
    type(fragmentInputParameters),intent(inout) :: input_frag

    nullify(input_frag%frag_index)
    nullify(input_frag%charge)
    !nullify(input_frag%frag_info)
    nullify(input_frag%label)
    nullify(input_frag%dirname)

  end subroutine nullifyInputFragParameters


subroutine create_log_file(iproc,inputs)

  use module_base
  use module_types
  use module_input
  use yaml_strings
  use yaml_output

  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(inout) :: inputs
  !local variables
  integer :: ierr,ierror,lgt
  logical :: exists
  character(len=500) :: logfile,logfile_old,logfile_dir

  logfile=repeat(' ',len(logfile))
  logfile_old=repeat(' ',len(logfile_old))
  logfile_dir=repeat(' ',len(logfile_dir))
  !open the logfile if needed, and set stdout
  !if (trim(in%writing_directory) /= '.') then
  if (.true.) then
     !add the output directory in the directory name
     if (iproc == 0 .and. trim(inputs%writing_directory) /= '.') then
        call getdir(inputs%writing_directory,&
             len_trim(inputs%writing_directory),logfile,len(logfile),ierr)
        if (ierr /= 0) then
           write(*,*) "ERROR: cannot create writing directory '"&
                //trim(inputs%writing_directory) // "'."
           call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
        end if
     end if
     call MPI_BCAST(logfile,len(logfile),MPI_CHARACTER,0,bigdft_mpi%mpi_comm,ierr)
     lgt=min(len(inputs%writing_directory),len(logfile))
     inputs%writing_directory(1:lgt)=logfile(1:lgt)
     lgt=0
     call buffer_string(inputs%dir_output,len(inputs%dir_output),&
          trim(logfile),lgt,back=.true.)
     if (iproc ==0) then
        logfile=repeat(' ',len(logfile))
        if (len_trim(inputs%run_name) >0) then
!           logfile='log-'//trim(inputs%run_name)//trim(bigdft_run_id_toa())//'.yaml'
           logfile='log-'//trim(inputs%run_name)//'.yaml'
        else
           logfile='log'//trim(bigdft_run_id_toa())//'.yaml'
        end if
        !inquire for the existence of a logfile
        call yaml_map('<BigDFT> log of the run will be written in logfile',&
             trim(inputs%writing_directory)//trim(logfile),unit=6)
        inquire(file=trim(inputs%writing_directory)//trim(logfile),exist=exists)
        if (exists) then
           logfile_old=trim(inputs%writing_directory)//'logfiles'
           call getdir(logfile_old,&
                len_trim(logfile_old),logfile_dir,len(logfile_dir),ierr)
           if (ierr /= 0) then
              write(*,*) "ERROR: cannot create writing directory '" //trim(logfile_dir) // "'."
              call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
           end if
           logfile_old=trim(logfile_dir)//trim(logfile)
           logfile=trim(inputs%writing_directory)//trim(logfile)
           !change the name of the existing logfile
           lgt=index(logfile_old,'.yaml')
           call buffer_string(logfile_old,len(logfile_old),&
                trim(adjustl(yaml_time_toa()))//'.yaml',lgt)
           call movefile(trim(logfile),len_trim(logfile),trim(logfile_old),len_trim(logfile_old),ierr)
           if (ierr /= 0) then
              write(*,*) "ERROR: cannot move logfile '"//trim(logfile)
              write(*,*) '                      into '//trim(logfile_old)// "'."
              call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
           end if
           call yaml_map('<BigDFT> Logfile existing, renamed into',&
                trim(logfile_old),unit=6)

        else
           logfile=trim(inputs%writing_directory)//trim(logfile)
        end if
        !Create stream and logfile
        call yaml_set_stream(unit=70,filename=trim(logfile),record_length=92,istat=ierr)
        !create that only if the stream is not already present, otherwise print a warning
        if (ierr == 0) then
           call input_set_stdout(unit=70)
           if (len_trim(inputs%run_name) == 0) then
              call f_malloc_set_status(unit=70, &
                   & logfile_name='malloc' // trim(bigdft_run_id_toa()) // '.prc')
           else
              call f_malloc_set_status(unit=70, &
                   & logfile_name='malloc-' // trim(inputs%run_name) // '.prc')
           end if
           !call memocc_set_stdout(unit=70)
        else
           call yaml_warning('Logfile '//trim(logfile)//' cannot be created, stream already present. Ignoring...')
        end if
     end if
  else
     !use stdout, do not crash if unit is present
     if (iproc==0) call yaml_set_stream(record_length=92,istat=ierr)
  end if
    
END SUBROUTINE create_log_file


!>  Free all dynamically allocated memory from the kpt input file.
subroutine free_kpt_variables(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  character(len=*), parameter :: subname='free_kpt_variables'
  integer :: i_stat, i_all

  if (associated(in%gen_kpt)) then
     i_all=-product(shape(in%gen_kpt))*kind(in%gen_kpt)
     deallocate(in%gen_kpt,stat=i_stat)
     call memocc(i_stat,i_all,'in%gen_kpt',subname)
  end if
  if (associated(in%gen_wkpt)) then
     i_all=-product(shape(in%gen_wkpt))*kind(in%gen_wkpt)
     deallocate(in%gen_wkpt,stat=i_stat)
     call memocc(i_stat,i_all,'in%gen_wkpt',subname)
  end if
  if (associated(in%kptv)) then
     i_all=-product(shape(in%kptv))*kind(in%kptv)
     deallocate(in%kptv,stat=i_stat)
     call memocc(i_stat,i_all,'in%kptv',subname)
  end if
  if (associated(in%nkptsv_group)) then
     i_all=-product(shape(in%nkptsv_group))*kind(in%nkptsv_group)
     deallocate(in%nkptsv_group,stat=i_stat)
     call memocc(i_stat,i_all,'in%nkptsv_group',subname)
  end if
  nullify(in%gen_kpt)
  nullify(in%gen_wkpt)
  nullify(in%kptv)
  nullify(in%nkptsv_group)
end subroutine free_kpt_variables

!>  Free all dynamically allocated memory from the geopt input file.
subroutine free_geopt_variables(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  character(len=*), parameter :: subname='free_geopt_variables'
  integer :: i_stat, i_all

  if (associated(in%qmass)) then
     i_all=-product(shape(in%qmass))*kind(in%qmass)
     deallocate(in%qmass,stat=i_stat)
     call memocc(i_stat,i_all,'in%qmass',subname)
  end if
  nullify(in%qmass)
end subroutine free_geopt_variables

!>  Free all dynamically allocated memory from the input variable structure.
subroutine free_input_variables(in)
  use module_base
  use module_types
  use module_xc
  implicit none
  type(input_variables), intent(inout) :: in
  character(len=*), parameter :: subname='free_input_variables'

!!$  if(in%linear /= INPUT_IG_OFF .and. in%linear /= INPUT_IG_LIG) &
!!$       & call deallocateBasicArraysInput(in%lin)

  call free_geopt_variables(in)
  call free_kpt_variables(in)
  call deallocateBasicArraysInput(in%lin)
  call deallocateInputFragArrays(in%frag)

  ! Free the libXC stuff if necessary, related to the choice of in%ixc.
  call xc_end()

!!$  if (associated(in%Gabs_coeffs) ) then
!!$     i_all=-product(shape(in%Gabs_coeffs))*kind(in%Gabs_coeffs)
!!$     deallocate(in%Gabs_coeffs,stat=i_stat)
!!$     call memocc(i_stat,i_all,'in%Gabs_coeffs',subname)
!!$  end if

  ! Stop the signaling stuff.
  !Destroy C wrappers on Fortran objects,
  ! and stop the GMainLoop.
  if (in%gmainloop /= 0.d0) then
     call bigdft_signals_free(in%gmainloop)
  end if
END SUBROUTINE free_input_variables


!> Assign default values for ABSCALC variables
subroutine abscalc_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(out) :: in

  in%c_absorbtion=.false.
  in%potshortcut=0
  in%iat_absorber=0
  in%abscalc_bottomshift=0
  in%abscalc_S_do_cg=.false.
  in%abscalc_Sinv_do_cg=.false.
END SUBROUTINE abscalc_input_variables_default

!> Assign default values for frequencies variables
!!    freq_alpha: frequencies step for finite difference = alpha*hx, alpha*hy, alpha*hz
!!    freq_order; order of the finite difference (2 or 3 i.e. 2 or 4 points)
!!    freq_method: 1 - systematic moves of atoms over each direction
subroutine frequencies_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(out) :: in

  in%freq_alpha=1.d0/real(64,kind(1.d0))
  in%freq_order=2
  in%freq_method=1

END SUBROUTINE frequencies_input_variables_default

subroutine dft_input_analyse(iproc, in, dict_dft)
  use module_base
  use module_types
  use module_input
  use dictionaries
  use module_xc
  use yaml_output
  use module_input_keys
  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(inout) :: in
  type(dictionary), pointer :: dict_dft

  !grid spacings (here the profiles can be used is we already read PSPs)
  in%hx = dict_dft//HGRIDS//0
  in%hy = dict_dft//HGRIDS//1
  in%hz = dict_dft//HGRIDS//2

  !coarse and fine radii around atoms
  in%crmult = dict_dft//RMULT//0
  in%frmult = dict_dft//RMULT//1

  !XC functional (ABINIT XC codes)
  in%ixc = dict_dft//IXC

  !charge and electric field
  in%ncharge = dict_dft//NCHARGE
  in%elecfield(1) = dict_dft//ELECFIELD//0
  in%elecfield(2) = dict_dft//ELECFIELD//1
  in%elecfield(3) = dict_dft//ELECFIELD//2

  !spin and polarization
  in%nspin = dict_dft//NSPIN
  in%mpol = dict_dft//MPOL

  ! Initialise XC calculation
  if (in%ixc < 0) then
     call xc_init(in%ixc, XC_MIXED, in%nspin)
  else
     call xc_init(in%ixc, XC_ABINIT, in%nspin)
  end if

  !convergence parameters
  in%gnrm_cv = dict_dft//GNRM_CV
  in%itermax = dict_dft//ITERMAX
  in%nrepmax = dict_dft//NREPMAX

  !convergence parameters
  in%ncong = dict_dft//NCONG
  in%idsx = dict_dft//IDSX
  !does not make sense a DIIS history longer than the number of iterations
  !only if the iscf is not particular
  in%idsx = min(in%idsx, in%itermax)

  !dispersion parameter
  in%dispersion = dict_dft//DISPERSION
    
  ! Now the variables which are to be used only for the last run
  in%inputPsiId = dict_dft//INPUTPSIID
  in%output_wf_format = dict_dft//OUTPUT_WF
  in%output_denspot = dict_dft//OUTPUT_DENSPOT

  !project however the wavefunction on gaussians if asking to write them on disk
  ! But not if we use linear scaling version (in%inputPsiId >= 100)
  in%gaussian_help=(in%inputPsiId >= 10 .and. in%inputPsiId < 100)

  !switch on the gaussian auxiliary treatment 
  !and the zero of the forces
  if (in%inputPsiId == 10) then
     in%inputPsiId = 0
  else if (in%inputPsiId == 13) then
     in%inputPsiId = 2
  end if
  ! Setup out grid parameters.
  if (in%output_denspot >= 0) then
     in%output_denspot_format = in%output_denspot / 10
  else
     in%output_denspot_format = output_denspot_FORMAT_CUBE
     in%output_denspot = abs(in%output_denspot)
  end if
  in%output_denspot = modulo(in%output_denspot, 10)

  ! Tail treatment.
  in%rbuf = dict_dft//RBUF
  in%ncongt = dict_dft//NCONGT

  !davidson treatment
  in%norbv = dict_dft//NORBV
  in%nvirt = dict_dft//NVIRT
  in%nplot = dict_dft//NPLOT

  ! Line to disable automatic behaviours (currently only symmetries).
  in%disableSym = dict_dft//DISABLE_SYM

  !define whether there should be a last_run after geometry optimization
  !also the mulliken charge population should be inserted
  if ((in%rbuf > 0.0_gp) .or. in%output_wf_format /= WF_FORMAT_NONE .or. &
       in%output_denspot /= output_denspot_NONE .or. in%norbv /= 0) then
     in%last_run=-1 !last run to be done depending of the external conditions
  else
     in%last_run=0
  end if
end subroutine dft_input_analyse

subroutine kpt_input_analyse(iproc, in, dict, sym, geocode, alat)
  use module_base
  use module_types
  use defs_basis
  use m_ab6_kpoints
  use yaml_output
  use module_input_keys
  use dictionaries
  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(inout) :: in
  type(dictionary), pointer :: dict
  type(symmetry_data), intent(in) :: sym
  character(len = 1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  real(gp), intent(in) :: alat(3)
  !local variables
  logical :: lstat
  character(len=*), parameter :: subname='kpt_input_analyse'
  integer :: i_stat,ierror,i,nshiftk, ngkpt_(3), ikpt, j, ncount, nseg, iseg_, ngranularity_
  real(gp) :: kptrlen_, shiftk_(3,8), norm, alat_(3)
  character(len = 6) :: method
  
  ! Set default values.
  in%gen_nkpt=1
  in%nkptv=0
  in%ngroups_kptv=1

  call free_kpt_variables(in)
  nullify(in%kptv, in%nkptsv_group)
  nullify(in%gen_kpt, in%gen_wkpt)

  method = dict // KPT_METHOD
  if (input_keys_equal(trim(method), 'auto')) then
     kptrlen_ = dict // KPTRLEN
     if (geocode == 'F') then
        in%gen_nkpt = 1
        allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
        call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
        in%gen_kpt = 0.
        allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
        call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
        in%gen_wkpt = 1.
     else
        call kpoints_get_auto_k_grid(sym%symObj, in%gen_nkpt, in%gen_kpt, in%gen_wkpt, &
             & kptrlen_, ierror)
        if (ierror /= AB6_NO_ERROR) then
           if (iproc==0) &
                & call yaml_warning("ERROR: cannot generate automatic k-point grid." // &
                & " Error code is " // trim(yaml_toa(ierror,fmt='(i0)')))
           stop
        end if
        !assumes that the allocation went through
        call memocc(0,in%gen_kpt,'in%gen_kpt',subname)
        call memocc(0,in%gen_wkpt,'in%gen_wkpt',subname)
     end if
  else if (input_keys_equal(trim(method), 'mpgrid')) then
     !take the points of Monkhorst-pack grid
     ngkpt_(1) = dict // NGKPT // 0
     ngkpt_(2) = dict // NGKPT // 1
     ngkpt_(3) = dict // NGKPT // 2
     if (geocode == 'S') ngkpt_(2) = 1
     !shift
     nshiftk = dict_len(dict//SHIFTK)
     !read the shifts
     shiftk_=0.0_gp
     do i=1,nshiftk
        shiftk_(1,i) = dict // SHIFTK // (i-1) // 0
        shiftk_(2,i) = dict // SHIFTK // (i-1) // 1
        shiftk_(3,i) = dict // SHIFTK // (i-1) // 2
     end do

     !control whether we are giving k-points to Free BC
     if (geocode == 'F') then
        if (iproc==0 .and. (maxval(ngkpt_) > 1 .or. maxval(abs(shiftk_)) > 0.)) &
             & call yaml_warning('Found input k-points with Free Boundary Conditions, reduce run to Gamma point')
        in%gen_nkpt = 1
        allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
        call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
        in%gen_kpt = 0.
        allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
        call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
        in%gen_wkpt = 1.
     else
        call kpoints_get_mp_k_grid(sym%symObj, in%gen_nkpt, in%gen_kpt, in%gen_wkpt, &
             & ngkpt_, nshiftk, shiftk_, ierror)
        if (ierror /= AB6_NO_ERROR) then
           if (iproc==0) &
                & call yaml_warning("ERROR: cannot generate MP k-point grid." // &
                & " Error code is " // trim(yaml_toa(ierror,fmt='(i0)')))
           stop
        end if
        !assumes that the allocation went through
        call memocc(0,in%gen_kpt,'in%gen_kpt',subname)
        call memocc(0,in%gen_wkpt,'in%gen_wkpt',subname)
     end if
  else if (input_keys_equal(trim(method), 'manual')) then
     in%gen_nkpt = max(1, dict_len(dict//KPT))
     if (geocode == 'F' .and. in%gen_nkpt > 1) then
        if (iproc==0) call yaml_warning('Found input k-points with Free Boundary Conditions, reduce run to Gamma point')
        in%gen_nkpt = 1
     end if
     allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
     call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
     allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
     call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
     norm=0.0_gp
     do i=1,in%gen_nkpt
        in%gen_kpt(1, i) = dict // KPT // (i-1) // 0
        in%gen_kpt(2, i) = dict // KPT // (i-1) // 1
        in%gen_kpt(3, i) = dict // KPT // (i-1) // 2
        if (geocode == 'S' .and. in%gen_kpt(2,i) /= 0.) then
           in%gen_kpt(2,i) = 0.
           if (iproc==0) call yaml_warning('Surface conditions, supressing k-points along y.')
        end if
        in%gen_wkpt(i) = dict // WKPT // (i-1)
        if (geocode == 'F') then
           in%gen_kpt = 0.
           in%gen_wkpt = 1.
        end if
        norm=norm+in%gen_wkpt(i)
     end do
     ! We normalise the weights.
     in%gen_wkpt(:)=in%gen_wkpt/norm
  else
     if (iproc==0) &
          & call yaml_warning("ERROR: wrong k-point sampling method (" // &
          & trim(method) // ").")
     stop
  end if

  ! Convert reduced coordinates into BZ coordinates.
  alat_ = alat
  if (geocode /= 'P') alat_(2) = 1.0_gp
  if (geocode == 'F') then
     alat_(1)=1.0_gp
     alat_(3)=1.0_gp
  end if
  do i = 1, in%gen_nkpt, 1
     in%gen_kpt(:, i) = in%gen_kpt(:, i) / alat_(:) * two_pi
  end do
 
  in%band_structure_filename=''
  lstat = dict // BANDS
  if (lstat) then
     !calculate the number of groups of for the band structure
     in%nkptv=1
     nseg = dict_len(dict // ISEG)
     do i=1,nseg
        iseg_ = dict // ISEG // (i-1)
        in%nkptv=in%nkptv+iseg_
     end do
     ngranularity_ = dict // NGRANULARITY

     in%ngroups_kptv=&
          ceiling(real(in%nkptv,gp)/real(ngranularity_,gp))

     allocate(in%nkptsv_group(in%ngroups_kptv+ndebug),stat=i_stat)
     call memocc(i_stat,in%nkptsv_group,'in%nkptsv_group',subname)

     ncount=0
     do i=1,in%ngroups_kptv-1
        !if ngranularity is bigger than nkptv  then ngroups is one
        in%nkptsv_group(i)=ngranularity_
        ncount=ncount+ngranularity_
     end do
     !put the rest in the last group
     in%nkptsv_group(in%ngroups_kptv)=in%nkptv-ncount

     allocate(in%kptv(3,in%nkptv+ndebug),stat=i_stat)
     call memocc(i_stat,in%kptv,'in%kptv',subname)

     ikpt = 0
     do i=1,nseg
        iseg_ = dict // ISEG // (i-1)
        ikpt=ikpt+iseg_
        in%kptv(1,ikpt) = dict // KPTV // (ikpt - 1) // 0
        in%kptv(2,ikpt) = dict // KPTV // (ikpt - 1) // 1
        in%kptv(3,ikpt) = dict // KPTV // (ikpt - 1) // 2
        !interpolate the values
        do j=ikpt-iseg_+1,ikpt-1
           in%kptv(:,j)=in%kptv(:,ikpt-iseg_) + &
                (in%kptv(:,ikpt)-in%kptv(:,ikpt-iseg_)) * &
                real(j-ikpt+iseg_,gp)/real(iseg_, gp)
        end do
     end do

     ! Convert reduced coordinates into BZ coordinates.
     do i = 1, in%nkptv, 1
        in%kptv(:, i) = in%kptv(:, i) / alat_(:) * two_pi
     end do

     if (has_key(dict, BAND_STRUCTURE_FILENAME)) then
        in%band_structure_filename = dict // BAND_STRUCTURE_FILENAME
        !since a file for the local potential is already given, do not perform ground state calculation
        if (iproc==0) then
           write(*,'(1x,a)')'Local Potential read from file, '//trim(in%band_structure_filename)//&
                ', do not optimise GS wavefunctions'
        end if
        in%nrepmax=0
        in%itermax=0
        in%itrpmax=0
        in%inputPsiId=-1000 !allocate empty wavefunctions
        in%output_denspot=0
     end if
  else
     in%nkptv = 0
     allocate(in%kptv(3,in%nkptv+ndebug),stat=i_stat)
     call memocc(i_stat,in%kptv,'in%kptv',subname)
  end if

  if (in%nkptv > 0 .and. geocode == 'F' .and. iproc == 0) &
       & call yaml_warning('Defining a k-point path in free boundary conditions.') 
END SUBROUTINE kpt_input_analyse

!> Read the input variables needed for the geometry optimisation
!! Every argument should be considered as mandatory
subroutine geopt_input_analyse(iproc,in,dict)
  use module_base
  use module_types
  use module_input_keys
  use dictionaries
  use yaml_output
  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(inout) :: in
  type(dictionary), pointer :: dict
  !local variables
  character(len=*), parameter :: subname='geopt_input_analyse'
  integer :: i_stat,i
  character(len = max_field_length) :: prof, meth
  real(gp) :: betax_, dtmax_

  ! Additional treatments.
  meth = dict // GEOPT_METHOD
  if (input_keys_equal(trim(meth), "FIRE")) then
     prof = input_keys_get_source(dict, DTMAX)
     if (trim(prof) == "default") then
        betax_ = dict // BETAX
        call set(dict // DTMAX, 0.25 * pi_param * sqrt(betax_), fmt = "(F7.4)")
     end if
     prof = input_keys_get_source(dict, DTINIT)
     if (trim(prof) == "default") then
        dtmax_ = dict // DTMAX
        call set(dict // DTINIT, 0.5 * dtmax_, fmt = "(F7.4)")
     end if
  end if
!!$  call yaml_dict_dump(dict)

  call free_geopt_variables(in)

  !target stress tensor
  in%strtarget(:)=0.0_gp

  !geometry input parameters
  in%geopt_approach = dict // GEOPT_METHOD
  in%ncount_cluster_x = dict // NCOUNT_CLUSTER_X
  !here the parsing of the wavefunction history should be added
  in%wfn_history=1

  in%frac_fluct = dict // FRAC_FLUCT
  in%forcemax = dict // FORCEMAX
  in%randdis = dict // RANDDIS

  if (input_keys_equal(trim(in%geopt_approach), "AB6MD")) then
     in%nnos=0
     in%ionmov = dict // IONMOV
     in%dtion = dict // DTION
     if (in%ionmov == 6) then
        in%mditemp = dict // MDITEMP
     elseif (in%ionmov > 7) then
        in%mditemp = dict // MDITEMP
        in%mdftemp = dict // MDFTEMP
     end if

     if (in%ionmov == 8) then
        in%noseinert = dict // NOSEINERT
     else if (in%ionmov == 9) then
        in%friction = dict // FRICTION
        in%mdwall = dict // MDWALL
     else if (in%ionmov == 13) then
        in%nnos = dict_len(dict // QMASS)
        allocate(in%qmass(in%nnos+ndebug),stat=i_stat)
        call memocc(i_stat,in%qmass,'in%qmass',subname)
        do i=1,in%nnos-1
           in%qmass(i) = dict_len(dict // QMASS // (i-1))
        end do
        in%bmass = dict // BMASS
        in%vmass = dict // VMASS
     end if

     if (in%ionmov /= 13) then
        !the allocation of this pointer should be done in any case
        allocate(in%qmass(in%nnos+ndebug),stat=i_stat)
        call memocc(i_stat,in%qmass,'in%qmass',subname)
     end if
  else if (input_keys_equal(trim(in%geopt_approach),"DIIS")) then
     in%betax = dict // BETAX
     in%history = dict // HISTORY
  else
     in%betax = dict // BETAX
  end if

  if (input_keys_equal(trim(in%geopt_approach),"FIRE")) then
     in%dtinit = dict // DTINIT
     in%dtmax = dict // DTMAX
  endif
END SUBROUTINE geopt_input_analyse

subroutine mix_input_analyse(iproc,in,dict)
  use module_base
  use module_types
  use module_input_keys
  use dictionaries
  implicit none
  !Arguments
  integer, intent(in) :: iproc
  type(dictionary), pointer :: dict
  type(input_variables), intent(inout) :: in
  !local variables
  !n(c) character(len=*), parameter :: subname='mix_input_variables'

  in%iscf = dict // ISCF
  
  in%itrpmax = dict // ITRPMAX
  in%rpnrm_cv = dict // RPNRM_CV

  in%norbsempty = dict // NORBSEMPTY
  in%Tel = dict // TEL
  in%occopt = dict // OCCOPT

  in%alphamix = dict // ALPHAMIX
  in%alphadiis = dict //ALPHADIIS

  !put the startmix if the mixing has to be done
  if (in%iscf >  SCF_KIND_DIRECT_MINIMIZATION) in%gnrm_startmix=1.e300_gp

END SUBROUTINE mix_input_analyse

subroutine sic_input_analyse(iproc,in,dict,ixc_)
  use module_base
  use module_types
  use module_input_keys
  use dictionaries
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: iproc
  type(dictionary), pointer :: dict
  type(input_variables), intent(inout) :: in
  integer, intent(in) :: ixc_
  !local variables
  !n(c) character(len=*), parameter :: subname='sic_input_variables'

  in%SIC%approach(1:len(in%SIC%approach)) = dict // SIC_APPROACH
  in%SIC%alpha = dict // SIC_ALPHA
  if (input_keys_equal(trim(in%SIC%approach), "NK")) in%SIC%fref = dict // SIC_FREF
  in%SIC%ixc = ixc_

END SUBROUTINE sic_input_analyse

subroutine tddft_input_analyse(iproc,in,dict)
  use module_base
  use module_types
  use module_input_keys
  use dictionaries
  implicit none
  integer, intent(in) :: iproc
  type(dictionary), pointer :: dict
  type(input_variables), intent(inout) :: in
  !local variables
  !n(c) character(len=*), parameter :: subname='tddft_input_variables'

  !TD-DFT parameters
  in%tddft_approach(1:len(in%tddft_approach)) = dict // TDDFT_APPROACH

END SUBROUTINE tddft_input_analyse

!> Read the input variables which can be used for performances
subroutine perf_input_analyse(iproc,in,dict)
  use module_base
  use module_types
  use module_input_keys
  use yaml_strings
  use yaml_output
  use dictionaries
  implicit none
  integer, intent(in) :: iproc
  type(dictionary), pointer :: dict
  type(input_variables), intent(inout) :: in
  !local variables
  !n(c) character(len=*), parameter :: subname='perf_input_variables'
  integer :: ierr,ipos,i,iproc_node,nproc_node
  character(len = 7) :: val

  ! Performance and setting-up options
  ! ----------------------------------
  in%debug = dict // DEBUG
  if (.not. in%debug) then
     call f_malloc_set_status(output_level=1)
     !call memocc_set_state(1)
  end if
  in%ncache_fft = dict // FFTCACHE
  call set_cache_size(in%ncache_fft)
  in%verbosity = dict // VERBOSITY
  if (in%verbosity == 0 ) then
     call f_malloc_set_status(output_level=0)
     !call memocc_set_state(0)
  end if
  in%writing_directory = dict // OUTDIR
  !here the logfile should be opened in the usual way, differentiating between 
  ! logfiles in case of multiple taskgroups
  if (trim(in%writing_directory) /= '.' .or. bigdft_mpi%ngroup > 1) then
     call create_log_file(iproc,in)
  else
     !use stdout, do not crash if unit is present
     if (iproc==0) call yaml_set_stream(record_length=92,istat=ierr)
  end if

  !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
  if (iproc==0) then
     !start writing on logfile
     call yaml_new_document()
     !welcome screen
     call print_logo()
  end if
  if (bigdft_mpi%nproc >1) call processor_id_per_node(bigdft_mpi%iproc,bigdft_mpi%nproc,iproc_node,nproc_node)
  if (iproc ==0) then
     if (bigdft_mpi%nproc >1) call yaml_map('MPI tasks of root process node',nproc_node)
     call print_configure_options()
  end if

  ! Miscellaneous options
  ! ---------------------
  in%symTol = dict // TOLSYM
  in%projrad = dict // PROJRAD
  in%exctxpar = dict // EXCTXPAR
  in%inguess_geopt = dict // INGUESS_GEOPT

  ! Material acceleration options
  ! -----------------------------
  in%matacc=material_acceleration_null()
  val = dict // ACCEL
  if (input_keys_equal(trim(val), "CUDAGPU")) then
     in%matacc%iacceleration = 1
  else if (input_keys_equal(trim(val), "OCLGPU")) then
     in%matacc%iacceleration = 2
  else if (input_keys_equal(trim(val), "OCLCPU")) then
     in%matacc%iacceleration = 3
  else if (input_keys_equal(trim(val), "OCLACC")) then
     in%matacc%iacceleration = 4
  else 
     in%matacc%iacceleration = 0
  end if
  !determine desired OCL platform which is used for acceleration
  in%matacc%OCL_platform = dict // OCL_PLATFORM
  ipos=min(len(in%matacc%OCL_platform),len(trim(in%matacc%OCL_platform))+1)
  do i=ipos,len(in%matacc%OCL_platform)
     in%matacc%OCL_platform(i:i)=achar(0)
  end do
  in%matacc%OCL_devices = dict // OCL_DEVICES
  ipos=min(len(in%matacc%OCL_devices),len(trim(in%matacc%OCL_devices))+1)
  do i=ipos,len(in%matacc%OCL_devices)
     in%matacc%OCL_devices(i:i)=achar(0)
  end do
  in%matacc%PSolver_igpu = dict // PSOLVER_ACCEL

  ! Signaling parameters
  in%signaling = dict // SIGNALING
  in%signalTimeout = dict // SIGNALTIMEOUT
  in%domain = dict // DOMAIN

  !!@TODO to relocate
  GPUblas = dict // BLAS
  DistProjApply = dict // PSP_ONFLY 

  in%orthpar%directDiag = dict // IG_DIAG
  in%orthpar%norbpInguess = dict // IG_NORBP
  in%orthpar%iguessTol = dict // IG_TOL
  in%orthpar%methOrtho = dict // METHORTHO
  !Block size used for the orthonormalization
  in%orthpar%bsLow = dict // IG_BLOCKS // 0
  in%orthpar%bsUp  = dict // IG_BLOCKS // 1
  !Check after collecting all values
  if(.not.in%orthpar%directDiag .or. in%orthpar%methOrtho==1) then 
     write(*,'(1x,a)') 'Input Guess: Block size used for the orthonormalization (ig_blocks)'
     if(in%orthpar%bsLow==in%orthpar%bsUp) then
        write(*,'(5x,a,i0)') 'Take block size specified by user: ',in%orthpar%bsLow
     else if(in%orthpar%bsLow<in%orthpar%bsUp) then
        write(*,'(5x,2(a,i0))') 'Choose block size automatically between ',in%orthpar%bsLow,' and ',in%orthpar%bsUp
     else
        write(*,'(1x,a)') "ERROR: invalid values of inputs%bsLow and inputs%bsUp. Change them in 'inputs.perf'!"
        call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
     end if
     write(*,'(5x,a)') 'This values will be adjusted if it is larger than the number of orbitals.'
  end if

  in%rho_commun = dict // RHO_COMMUN
  in%PSolver_groupsize = dict // PSOLVER_GROUPSIZE
  in%unblock_comms = dict // UNBLOCK_COMMS

  !Use Linear scaling methods
  val = dict // LINEAR
  if (input_keys_equal(trim(val), "LIG")) then
     in%linear = INPUT_IG_LIG
  else if (input_keys_equal(trim(val), "FUL")) then
     in%linear = INPUT_IG_FULL
  else if (input_keys_equal(trim(val), "TMO")) then
     in%linear = INPUT_IG_TMO
  else
     in%linear = INPUT_IG_OFF
  end if

  in%store_index = dict // STORE_INDEX

  !block size for pdsyev/pdsygv, pdgemm (negative -> sequential)
  in%lin%blocksize_pdsyev = dict // PDSYEV_BLOCKSIZE
  in%lin%blocksize_pdgemm = dict // PDGEMM_BLOCKSIZE
  !max number of process uses for pdsyev/pdsygv, pdgemm
  in%lin%nproc_pdsyev = dict // MAXPROC_PDSYEV
  in%lin%nproc_pdgemm = dict // MAXPROC_PDGEMM
  !FOE: if the determinant of the interpolation matrix to find the Fermi energy
  !is smaller than this value, switch from cubic to linear interpolation.
  in%lin%ef_interpol_det = dict // EF_INTERPOL_DET
  in%lin%ef_interpol_chargediff = dict // EF_INTERPOL_CHARGEDIFF
  !determines whether a mixing step shall be preformed after the input guess !(linear version)
  in%lin%mixing_after_inputguess = dict // MIXING_AFTER_INPUTGUESS
  !determines whether the input guess support functions are orthogonalized iteratively (T) or in the standard way (F)
  in%lin%iterative_orthogonalization = dict // ITERATIVE_ORTHOGONALIZATION

  in%check_sumrho = dict // CHECK_SUMRHO
!  call input_var("mpi_groupsize",0, "number of MPI processes for BigDFT run (0=nproc)", in%mpi_groupsize)
END SUBROUTINE perf_input_analyse

