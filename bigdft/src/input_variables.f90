!> @file
!!  Routines to read and print input variables
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> This function returns a dictionary with all the input variables of a BigDFT run filled.
!! This dictionary is constructed from a updated version of the input variables dictionary
!! following the input files as defined  by the user
subroutine read_input_dict_from_files(radical,mpi_env,dict)
  use dictionaries
  use wrapper_MPI
  !use module_input_keys
  use public_keys
  use module_input_dicts, only: merge_input_file_to_dict
  use input_old_text_format
  use yaml_output
  use dynamic_memory
  implicit none
  character(len = *), intent(in) :: radical    !< The name of the run. use "input" if empty
  type(mpi_environment), intent(in) :: mpi_env !< The environment where the variables have to be updated
  type(dictionary), pointer :: dict            !< Input dictionary, has to be nullified at input
  !local variables
  integer :: ierr
  logical :: exists_default, exists_user
  character(len = max_field_length) :: fname
  character(len = 100) :: f0

  call f_routine(id='read_input_dict_from_files')

!!$  if (f_err_raise(associated(dict),'The output dictionary should be nullified at input',&
!!$       err_name='BIGDFT_RUNTIME_ERROR')) return
!!$
!!$  nullify(dict) !this is however put in the case the dictionary comes undefined

!!$  call dict_init(dict)
  if (trim(radical) /= "" .and. trim(radical) /= "input") &
       & call set(dict // RADICAL_NAME, radical)

  ! Handle error with master proc only.
  !LG: modified, better to handle errors with all the 
  !! processors now that each of the cores has its own way of dumping 
  !! error codes
  !if (mpi_env%iproc > 0) call f_err_set_callback(f_err_ignore)

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
     ! Parse all files.
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
     !call set_inputfile(f0, radical, 'lin')
     call read_lin_and_frag_from_text_format(mpi_env%iproc,dict,trim(radical)) !as it also reads fragment

     call set_inputfile(f0, radical, 'neb')
     call read_neb_from_text_format(mpi_env%iproc,dict//GEOPT_VARIABLES, trim(f0))
  end if

  !LG modfication of errors (see above)
  !in case it should be restored the bigdft_severe shoudl be called instead
  !if (mpi_env%iproc > 0) call f_err_severe_restore()

  ! We put a barrier here to be sure that non master proc will be stopped
  ! by any issue on the master proc.
  call mpi_barrier(mpi_env%mpi_comm, ierr)

  call f_release_routine()
end subroutine read_input_dict_from_files


!> Fill the input_variables and atoms_data structures from the information
!! contained in the dictionary dict
!! the dictionary should be completes to fill all the information
subroutine inputs_from_dict(in, atoms, dict)
  use module_types
  use module_defs
  use yaml_output
  use module_interfaces, except => inputs_from_dict
  use dictionaries
  use module_input_keys
  use public_keys
  use module_input_dicts
  use dynamic_memory
  use module_xc
  use input_old_text_format, only: dict_from_frag
  use module_atoms, only: atoms_data,atoms_data_null,atomic_data_set_from_dict,check_atoms_positions
  use yaml_strings, only: f_strcpy
  use psp_projectors, only: PSPCODE_PAW
  use m_ab6_symmetry, only: symmetry_get_n_sym
  use interfaces_42_libpaw
  use bigdft_run, only: bigdft_run_id_toa
  implicit none
  !Arguments
  type(input_variables), intent(out) :: in
  type(atoms_data), intent(out) :: atoms
  type(dictionary), pointer :: dict
  !Local variables
  !type(dictionary), pointer :: profs, dict_frag
  logical :: found
  integer :: ierr, ityp, nelec_up, nelec_down, norb_max, jtype
  character(len = max_field_length) :: msg,filename
!  type(f_dict) :: dict
  type(dictionary), pointer :: dict_minimal, var, lvl

  integer, parameter :: pawlcutd = 10, pawlmix = 10, pawnphi = 13, pawntheta = 12, pawxcdev = 1
  integer, parameter :: xclevel = 1, usepotzero = 0
  integer :: nsym
  real(gp) :: gsqcut_shp

!  dict => dict//key

!  dict = dict//key

  call f_routine(id='inputs_from_dict')

  ! Atoms case.
  atoms = atoms_data_null()
  if (.not. has_key(dict, "posinp")) stop "missing posinp"
  call astruct_set_from_dict(dict // "posinp", atoms%astruct)

  ! Input variables case.
  call default_input_variables(in)

  !call yaml_map('Dictionary parsed',dict)

  ! Analyse the input dictionary and transfer it to in.
  call input_keys_validate(dict)

  ! extract also the minimal dictionary which is necessary to do this run
  call input_keys_fill_all(dict,dict_minimal)

  ! Transfer dict values into input_variables structure.
  lvl => dict_iter(dict)
  do while(associated(lvl))
     var => dict_iter(lvl)
     if (.not. associated(var)) then
        ! Scalar case.
        call input_set(in, trim(dict_key(lvl)), lvl)
     else
        do while(associated(var))
           call input_set(in, trim(dict_key(lvl)), var)
           var => dict_next(var)
        end do
     end if
     lvl => dict_next(lvl)
  end do

  call set_cache_size(in%ncache_fft)

  !status of the allocation verbosity and profiling
  !filename of the memory allocation status, if needed
  if (len_trim(in%run_name) == 0) then
     call f_strcpy(src=trim(in%writing_directory)//'/memstatus' // trim(bigdft_run_id_toa()) // '.yaml',&
          dest=filename)
  else
     call f_strcpy(src=trim(in%writing_directory)//'/memstatus-' // trim(in%run_name) // '.yaml',&
          dest=filename)
  end if

  if (.not. in%debug) then
     if (in%verbosity==3) then
        call f_malloc_set_status(output_level=1,&
             iproc=bigdft_mpi%iproc,logfile_name=filename)
     else
        call f_malloc_set_status(output_level=0,&
             iproc=bigdft_mpi%iproc)
     end if
  else
     call f_malloc_set_status(output_level=2,&
          iproc=bigdft_mpi%iproc,logfile_name=filename)
  end if

  call nullifyInputLinparameters(in%lin)
  call allocateBasicArraysInputLin(in%lin, atoms%astruct%ntypes)

  !First fill all the types by the default, then override by per-type values
  do jtype=1,atoms%astruct%ntypes
     var => dict_iter(dict//LIN_BASIS_PARAMS)
     do while(associated(var))
        call basis_params_set_dict(var,in%lin,jtype)
        var => dict_next(var)
     end do
     !then check if the objects exists in separate specifications
     if (has_key(dict//LIN_BASIS_PARAMS,trim(atoms%astruct%atomnames(jtype)))) then
        var => &
             dict_iter(dict//LIN_BASIS_PARAMS//trim(atoms%astruct%atomnames(jtype)))
     end if
     do while(associated(var))
        call basis_params_set_dict(var,in%lin,jtype)
        var => dict_next(var)
     end do
  end do

  call allocate_extra_lin_arrays(in%lin,atoms%astruct)

  ! Cross check values of input_variables.
  call input_analyze(in,atoms%astruct)

  ! Shake atoms, if required.
  call astruct_set_displacement(atoms%astruct, in%randdis)
  if (bigdft_mpi%nproc > 1) call mpibarrier(bigdft_mpi%mpi_comm)
  ! Update atoms with symmetry information
  call astruct_set_symmetries(atoms%astruct, in%disableSym, in%symTol, in%elecfield, in%nspin)

  call kpt_input_analyse(bigdft_mpi%iproc, in, dict//KPT_VARIABLES, &
       & atoms%astruct%sym, atoms%astruct%geocode, atoms%astruct%cell_dim)

  ! Add missing pseudo information.
  do ityp = 1, atoms%astruct%ntypes, 1
     call psp_dict_fill_all(dict, atoms%astruct%atomnames(ityp), in%ixc, &
          & in%projrad, in%crmult, in%frmult)
  end do

  ! Update atoms with pseudo information.
  call psp_dict_analyse(dict, atoms)
  call atomic_data_set_from_dict(dict,IG_OCCUPATION, atoms, in%nspin)

  ! Add multipole preserving information
  atoms%multipole_preserving = in%multipole_preserving

  ! Generate orbital occupation
  call read_n_orbitals(bigdft_mpi%iproc, nelec_up, nelec_down, norb_max, atoms, &
       & in%ncharge, in%nspin, in%mpol, in%norbsempty)
  if (norb_max == 0) norb_max = nelec_up + nelec_down ! electron gas case
  call occupation_set_from_dict(dict, OCCUPATION, &
       & in%gen_norbu, in%gen_norbd, in%gen_occup, &
       & in%gen_nkpt, in%nspin, in%norbsempty, nelec_up, nelec_down, norb_max)
  in%gen_norb = in%gen_norbu + in%gen_norbd

  ! Complement PAW initialisation.
  if (any(atoms%npspcode == PSPCODE_PAW)) then
     !gsqcut_shp = two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     gsqcut_shp = 2._gp * 2.2_gp / pi_param ** 2
     call symmetry_get_n_sym(atoms%astruct%sym%symObj, nsym, ierr)
     call pawinit(1, gsqcut_shp, pawlcutd, pawlmix, maxval(atoms%pawtab(:)%lmn_size) + 1, &
          & pawnphi, nsym, pawntheta, atoms%pawang, atoms%pawrad, 0, &
          & atoms%pawtab, pawxcdev, xclevel, usepotzero)
  end if
  
  if (in%gen_nkpt > 1 .and. in%gaussian_help) then
     if (bigdft_mpi%iproc==0) call yaml_warning('Gaussian projection is not implemented with k-point support')
     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
  end if
  
  if(in%inputpsiid == INPUT_PSI_LINEAR_AO .or. &
     in%inputpsiid == INPUT_PSI_MEMORY_LINEAR .or. &
     in%inputpsiid == INPUT_PSI_DISK_LINEAR) &
      DistProjApply=.true.
  if(in%linear /= INPUT_IG_OFF .and. in%linear /= INPUT_IG_LIG) then
     !only on the fly calculation
     DistProjApply=.true.
  end if

  !if other steps are supposed to be done leave the last_run to minus one
  !otherwise put it to one
  if (in%last_run == -1 .and. in%ncount_cluster_x <=1 .or. in%ncount_cluster_x <= 1) then
     in%last_run = 1
  end if

  ! Stop the code if it is trying to run GPU with spin=4
  if (in%nspin == 4 .and. in%matacc%iacceleration /= 0) then
     if (bigdft_mpi%iproc==0) call yaml_warning('GPU calculation not implemented with non-collinear spin')
     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
  end if

  !control atom positions
  call check_atoms_positions(atoms%astruct, (bigdft_mpi%iproc == 0))

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

  !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

  ! Warn for all INPUT_VAR_ILLEGAL errors.
  do while (f_err_pop(INPUT_VAR_ILLEGAL, add_msg = msg) /= 0)
     call yaml_warning(trim(msg))
  end do
  !check if an error has been found and raise an exception to be handled
  if (f_err_check()) then
     call f_err_throw('Error in reading input variables from dictionary',&
          err_name='BIGDFT_INPUT_VARIABLES_ERROR')
  end if

  ! not sure whether to actually make this an input variable or not so just set to false for now
  in%lin%diag_start=.false.

  !then fill also fragment variables
  in%lin%fragment_calculation=FRAG_VARIABLES .in. dict
  in%lin%calc_transfer_integrals=.false.
  in%lin%constrained_dft=.false.
  if (in%lin%fragment_calculation) then
     in%lin%constrained_dft=CONSTRAINED_DFT .in. dict // FRAG_VARIABLES
     found = TRANSFER_INTEGRALS .in. dict // FRAG_VARIABLES
     if (found) in%lin%calc_transfer_integrals=dict//FRAG_VARIABLES//TRANSFER_INTEGRALS
     call frag_from_dict(dict//FRAG_VARIABLES,in%frag)

!!$     ! again recheck
!!$     call dict_from_frag(in%frag,dict_frag)
!!$     call yaml_map('again',dict_frag)
!!$     call dict_free(dict_frag)
!!$     stop
  else
     call default_fragmentInputParameters(in%frag)
  end if
  if (bigdft_mpi%iproc==0) then
     call input_keys_dump(dict)
  end if

  !check whether a directory name should be associated for the data storage
  call check_for_data_writing_directory(bigdft_mpi%iproc,in)

  if (bigdft_mpi%iproc == 0)  call print_general_parameters(in,atoms)

  if (associated(dict_minimal) .and. bigdft_mpi%iproc == 0) then
     if (len_trim(in%run_name) == 0) then
        call f_strcpy(src='input_minimal.yaml',dest=filename)
     else
        call f_strcpy(src=trim(in%run_name)//'_minimal.yaml',dest=filename)
     end if

     call yaml_set_stream(unit=99971,filename=trim(in%writing_directory)//'/'//trim(filename),&
          record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
     if (ierr==0) then
        call yaml_comment('Minimal input file',hfill='-',unit=99971)
        call yaml_comment('This file indicates the minimal set of input variables which has to be given '//&
             'to perform the run. The code would produce the same output if this file is used as input.',unit=99971)
        call yaml_dict_dump(dict_minimal,unit=99971)
        call yaml_close_stream(unit=99971)
     else
        call yaml_warning('Failed to create input_minimal.yaml, error code='//trim(yaml_toa(ierr)))
     end if
  end if
  if (associated(dict_minimal)) call dict_free(dict_minimal)

  call f_release_routine()

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


!> Create the directory output
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
  
  in%refcnt=f_ref_new('inputs')

  in%matacc=material_acceleration_null()

  ! Default values.
  in%run_name = " "
  in%writing_directory = "."
  in%dir_output = "data"
  in%output_wf_format = WF_FORMAT_NONE
  in%output_denspot_format = output_denspot_FORMAT_CUBE
  nullify(in%gen_kpt)
  nullify(in%gen_wkpt)
  nullify(in%kptv)
  nullify(in%nkptsv_group)
  in%gen_norb = UNINITIALIZED(0)
  in%gen_norbu = UNINITIALIZED(0)
  in%gen_norbd = UNINITIALIZED(0)
  nullify(in%gen_occup)
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
  nullify(in%lin%kernel_cutoff_FOE)
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
  in%mditemp = UNINITIALIZED(in%mditemp)
  in%mdftemp = UNINITIALIZED(in%mdftemp)
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


!> Allocate the arrays for the input related to the fragment
subroutine allocateInputFragArrays(input_frag)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(fragmentInputParameters),intent(inout) :: input_frag

  ! Local variables
  character(len=*),parameter :: subname='allocateInputFragArrays'

  input_frag%frag_index = f_malloc_ptr(input_frag%nfrag,id='input_frag%frag_index')
  input_frag%charge = f_malloc_ptr(input_frag%nfrag,id='input_frag%charge')

  !allocate(input_frag%frag_info(input_frag%nfrag_ref,2), stat=i_stat)
  !call memocc(i_stat, input_frag%frag_info, 'input_frag%frag_info', subname)

  input_frag%label=f_malloc_str_ptr(len(input_frag%label),&
       input_frag%nfrag_ref,id='input_frag%label')
!!$  allocate(input_frag%label(input_frag%nfrag_ref), stat=i_stat)
!!$  call memocc(i_stat, input_frag%label, 'input_frag%label', subname)


  !f_malloc0_str_ptr should be used here
  input_frag%dirname=f_malloc_str_ptr(len(input_frag%dirname),&
       input_frag%nfrag_ref,id='input_frag%label')

!!$  allocate(input_frag%dirname(input_frag%nfrag_ref), stat=i_stat)
!!$  call memocc(i_stat, input_frag%dirname, 'input_frag%dirname', subname)

  !set the variables to their default value

end subroutine allocateInputFragArrays


!> Deallocate the arrays related to the input for the fragments
subroutine deallocateInputFragArrays(input_frag)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(fragmentInputParameters),intent(inout) :: input_frag

  ! Local variables
  character(len=*),parameter :: subname='deallocateInputFragArrays'

  !if(associated(input_frag%frag_info)) then
  !  i_all = -product(shape(input_frag%frag_info))*kind(input_frag%frag_info)
  !  deallocate(input_frag%frag_info,stat=i_stat)
  !  call memocc(i_stat,i_all,'input_frag%frag_info',subname)
  !  nullify(input_frag%frag_info)
  !end if 

  call f_free_ptr(input_frag%frag_index)
  
  
  call f_free_ptr(input_frag%charge)
  
  call f_free_str_ptr(len(input_frag%label),input_frag%label)
  call f_free_str_ptr(len(input_frag%dirname),input_frag%dirname)

!!$  if(associated(input_frag%label)) then
!!$     i_all = -product(shape(input_frag%label))*kind(input_frag%label)
!!$     deallocate(input_frag%label,stat=i_stat)
!!$     call memocc(i_stat,i_all,'input_frag%label',subname)
!!$     nullify(input_frag%label)
!!$  end if
!!$
!!$  if(associated(input_frag%dirname)) then
!!$     i_all = -product(shape(input_frag%dirname))*kind(input_frag%dirname)
!!$     deallocate(input_frag%dirname,stat=i_stat)
!!$     call memocc(i_stat,i_all,'input_frag%dirname',subname)
!!$     nullify(input_frag%dirname)
!!$  end if

end subroutine deallocateInputFragArrays

!> initialize fragment input parameters to their default value
subroutine default_fragmentInputParameters(frag)
  use module_types, only: fragmentInputParameters
  implicit none
  type(fragmentInputParameters),intent(out) :: frag

  !first nullify
  call nullifyInputFragParameters(frag)
  !then set defaults
  frag%nfrag_ref=1
  frag%nfrag=1
  !then allocate
  call allocateInputFragArrays(frag)
  !and fill to neutral values
  frag%label(1)=repeat(' ',len(frag%label))
  frag%dirname(1)=repeat(' ',len(frag%dirname))
  frag%frag_index(1)=1
  frag%charge(1)=0.0d0

end subroutine default_fragmentInputParameters

!> Nullify the parameters related to the fragments
subroutine nullifyInputFragParameters(input_frag)
  use module_types
  implicit none

  ! Calling arguments
  type(fragmentInputParameters),intent(inout) :: input_frag

  !default scalar variables
  nullify(input_frag%frag_index)
  nullify(input_frag%charge)
  !nullify(input_frag%frag_info)
  nullify(input_frag%label)
  nullify(input_frag%dirname)

end subroutine nullifyInputFragParameters


!> Creation of the log file (by default log.yaml)
!>  Free all dynamically allocated memory from the kpt input file.
subroutine free_kpt_variables(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  character(len=*), parameter :: subname='free_kpt_variables'

!!$  if (associated(in%gen_kpt)) then
!!$     i_all=-product(shape(in%gen_kpt))*kind(in%gen_kpt)
!!$     deallocate(in%gen_kpt,stat=i_stat)
!!$     call memocc(i_stat,i_all,'in%gen_kpt',subname)
!!$  end if
!!$  if (associated(in%gen_wkpt)) then
!!$     i_all=-product(shape(in%gen_wkpt))*kind(in%gen_wkpt)
!!$     deallocate(in%gen_wkpt,stat=i_stat)
!!$     call memocc(i_stat,i_all,'in%gen_wkpt',subname)
!!$  end if
  call f_free_ptr(in%gen_kpt)
  call f_free_ptr(in%gen_wkpt)
  call f_free_ptr(in%kptv)
  call f_free_ptr(in%nkptsv_group)
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
  ! integer :: i_stat, i_all

  if (associated(in%qmass)) then
     call f_free_ptr(in%qmass)
  end if
  nullify(in%qmass)
end subroutine free_geopt_variables


!>  Free all dynamically allocated memory from the input variable structure.
subroutine free_input_variables(in)
  use module_base
  use module_types
  use module_xc
  use dynamic_memory, only: f_free_ptr
  implicit none
  type(input_variables), intent(inout) :: in
  character(len=*), parameter :: subname='free_input_variables'

  !check if freeing is possible
  call f_ref_free(in%refcnt)

  call free_geopt_variables(in)
  call free_kpt_variables(in)
  call f_free_ptr(in%gen_occup)
  call deallocateBasicArraysInput(in%lin)
  call deallocateInputFragArrays(in%frag)

  ! Free the libXC stuff if necessary, related to the choice of in%ixc.
!!$  call xc_end(in%xcObj)

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


!> Cross check values of input_variables.
!! and change if necessary
subroutine input_analyze(in,astruct)
  use module_types, only: input_variables,output_denspot_FORMAT_CUBE, &
       output_denspot_NONE, WF_FORMAT_NONE, KERNELMODE_DIRMIN,&
       KERNELMODE_DIAG, KERNELMODE_FOE, MIXINGMODE_DENS, MIXINGMODE_POT, &
       LINEAR_DIRECT_MINIMIZATION, LINEAR_MIXDENS_SIMPLE, &
       LINEAR_MIXPOT_SIMPLE, LINEAR_FOE
  use module_atoms, only: atomic_structure
  use module_base
  use module_input_keys, only: input_keys_equal
  implicit none
  type(input_variables), intent(inout) :: in
  type(atomic_structure), intent(in) :: astruct

  integer :: ierr

  call f_routine(id='input_analyze')

  ! the PERF variables -----------------------------------------------------
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
  !@todo also the inputguess variable should be checked if BC are nonFree

  ! the DFT variables ------------------------------------------------------
  in%SIC%ixc = in%ixc

  in%idsx = min(in%idsx, in%itermax)

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

  !define whether there should be a last_run after geometry optimization
  !also the mulliken charge population should be inserted
  if ((in%rbuf > 0.0_gp) .or. in%output_wf_format /= WF_FORMAT_NONE .or. &
       in%output_denspot /= output_denspot_NONE .or. in%norbv /= 0) then
     in%last_run=-1 !last run to be done depending of the external conditions
  else
     in%last_run=0
  end if

  if (astruct%geocode == 'F' .or. astruct%nat == 0) then
     !Disable the symmetry
     in%disableSym = .true.
  end if
  
  ! the GEOPT variables ----------------------------------------------------
  !target stress tensor
  in%strtarget(:)=0.0_gp

  if (input_keys_equal(trim(in%geopt_approach), "AB6MD")) then
     if (in%ionmov /= 13) then
        in%nnos=0
        in%qmass = f_malloc_ptr(in%nnos, id = "in%qmass")
     end if
  end if

  ! Determine the SCF mode
  select case (in%lin%kernel_mode)
  case (KERNELMODE_DIRMIN)
      in%lin%scf_mode = LINEAR_DIRECT_MINIMIZATION
  case (KERNELMODE_DIAG)
      select case (in%lin%mixing_mode)
      case (MIXINGMODE_DENS)
          in%lin%scf_mode = LINEAR_MIXDENS_SIMPLE
      case (MIXINGMODE_POT)
          in%lin%scf_mode = LINEAR_MIXPOT_SIMPLE
      case default
          stop 'wrong value of in%lin%mixing_mode'
      end select
  case (KERNELMODE_FOE)
      in%lin%scf_mode = LINEAR_FOE
  case default
      stop 'wrong value of in%lin%kernel_mode'
  end select

  ! It is not possible to use both the old and the new Pulay correction at the same time
  if (in%lin%pulay_correction .and. in%lin%new_pulay_correction) then
     call f_err_throw('It is not possible to use both the old and the new Pulay correction at the same time!',&
          err_name='BIGDFT_INPUT_VARIABLES_ERROR')
  end if


  call f_release_routine()
END SUBROUTINE input_analyze


!> Analyse the kpt input and calculates k points if needed
subroutine kpt_input_analyse(iproc, in, dict, sym, geocode, alat)
  use module_base
  use module_types
  use module_atoms, only: symmetry_data
  use defs_basis
  use m_ab6_kpoints
  use yaml_output
  use module_input_keys, only: input_keys_equal
  use public_keys
  use dictionaries
  implicit none
  !Arguments
  integer, intent(in) :: iproc
  type(input_variables), intent(inout) :: in
  type(dictionary), pointer, intent(in) :: dict
  type(symmetry_data), intent(in) :: sym
  character(len = 1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  real(gp), dimension(3), intent(in) :: alat
  !local variables
  logical :: lstat
  character(len=*), parameter :: subname='kpt_input_analyse'
  integer :: ierror,i, nshiftk, ikpt, j, ncount, nseg, iseg_, ngranularity_
  integer, dimension(3) :: ngkpt_
  real(gp), dimension(3) :: alat_
  real(gp), dimension(3,8) :: shiftk_
  real(gp) :: kptrlen_, norm
  character(len = 6) :: method
  real(gp), dimension(:,:), pointer :: gen_kpt   !< K points coordinates
  real(gp), dimension(:), pointer :: gen_wkpt    !< Weights of k points
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
!!$        allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
!!$        in%gen_kpt = 0.
        in%gen_kpt=f_malloc0_ptr([3, in%gen_nkpt],id='gen_kpt')

!!$        allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
        in%gen_kpt=f_malloc_ptr(in%gen_nkpt,id='gen_wkpt')
                
        in%gen_wkpt = 1.
     else
        call kpoints_get_auto_k_grid(sym%symObj, in%gen_nkpt, gen_kpt, gen_wkpt, &
             & kptrlen_, ierror)
        if (ierror /= AB7_NO_ERROR) then
           if (iproc==0) &
                & call yaml_warning("ERROR: cannot generate automatic k-point grid." // &
                & " Error code is " // trim(yaml_toa(ierror,fmt='(i0)')))
           stop
        end if
        !assumes that the allocation went through (arrays allocated by abinit routines)
        in%gen_kpt=f_malloc_ptr(src=gen_kpt,id='gen_kpt')
        in%gen_wkpt=f_malloc_ptr(src=gen_wkpt,id='gen_wkpt')
        deallocate(gen_kpt,gen_wkpt)
!!$        call memocc(0,in%gen_kpt,'in%gen_kpt',subname)
!!$        call memocc(0,in%gen_wkpt,'in%gen_wkpt',subname)
     end if
  else if (input_keys_equal(trim(method), 'mpgrid')) then
     !take the points of Monkhorst-pack grid
     ngkpt_(1:3) = dict // NGKPT
     if (geocode == 'S') ngkpt_(2) = 1
     !shift
     nshiftk = dict_len(dict//SHIFTK)
     !read the shifts
     shiftk_=0.0_gp
     do i=1,nshiftk
        shiftk_(1:3,i) = dict // SHIFTK // (i-1)
     end do

     !control whether we are giving k-points to Free BC
     if (geocode == 'F') then
        if (iproc==0 .and. (maxval(ngkpt_) > 1 .or. maxval(abs(shiftk_)) > 0.)) &
             & call yaml_warning('Found input k-points with Free Boundary Conditions, reduce run to Gamma point')
        in%gen_nkpt = 1
        in%gen_kpt=f_malloc0_ptr([3, in%gen_nkpt],id='gen_kpt')
!!$        allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
!!$        in%gen_kpt = 0.
!!$        allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
        in%gen_kpt=f_malloc_ptr(in%gen_nkpt,id='gen_wkpt')
        in%gen_wkpt = 1.
     else
        call kpoints_get_mp_k_grid(sym%symObj, in%gen_nkpt, gen_kpt, gen_wkpt, &
             & ngkpt_, nshiftk, shiftk_, ierror)
        if (ierror /= AB7_NO_ERROR) then
           if (iproc==0) &
                & call yaml_warning("ERROR: cannot generate MP k-point grid." // &
                & " Error code is " // trim(yaml_toa(ierror,fmt='(i0)')))
           stop
        end if
        !assumes that the allocation went through 
        !(arrays allocated by abinit routines)
        in%gen_kpt=f_malloc_ptr(src=gen_kpt,id='gen_kpt')
        in%gen_wkpt=f_malloc_ptr(src=gen_wkpt,id='gen_wkpt')
        deallocate(gen_kpt,gen_wkpt)
!!$        call memocc(0,in%gen_kpt,'in%gen_kpt',subname)
!!$        call memocc(0,in%gen_wkpt,'in%gen_wkpt',subname)
     end if
  else if (input_keys_equal(trim(method), 'manual')) then
     in%gen_nkpt = max(1, dict_len(dict//KPT))
     if (geocode == 'F' .and. in%gen_nkpt > 1) then
        if (iproc==0) call yaml_warning('Found input k-points with Free Boundary Conditions, reduce run to Gamma point')
        in%gen_nkpt = 1
     end if
     in%gen_kpt=f_malloc_ptr([3, in%gen_nkpt],id='gen_kpt')
     in%gen_wkpt=f_malloc_ptr(in%gen_nkpt,id='gen_wkpt')

!!$     allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
!!$     call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
!!$     allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
!!$     call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
     norm=0.0_gp
     do i=1,in%gen_nkpt
        in%gen_kpt(1:3, i) = dict // KPT // (i-1)
        if (geocode == 'S' .and. in%gen_kpt(2,i) /= 0.) then
           in%gen_kpt(2,i) = 0.
           if (iproc==0) call yaml_warning('Surface conditions, suppressing k-points along y.')
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

     in%nkptsv_group = f_malloc_ptr(in%ngroups_kptv,id='in%nkptsv_group')

     ncount=0
     do i=1,in%ngroups_kptv-1
        !if ngranularity is bigger than nkptv  then ngroups is one
        in%nkptsv_group(i)=ngranularity_
        ncount=ncount+ngranularity_
     end do
     !put the rest in the last group
     in%nkptsv_group(in%ngroups_kptv)=in%nkptv-ncount

     in%kptv = f_malloc_ptr((/ 3, in%nkptv /),id='in%kptv')

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
     in%kptv = f_malloc_ptr((/ 3, in%nkptv /),id='in%kptv')
  end if

  if (in%nkptv > 0 .and. geocode == 'F' .and. iproc == 0) &
       & call yaml_warning('Defining a k-point path in free boundary conditions.') 

END SUBROUTINE kpt_input_analyse
