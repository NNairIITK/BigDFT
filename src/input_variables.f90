!!> @file
!!  Routines to read and print input variables
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Do all initialisation for all different files of BigDFT. 
!! Set default values if not any.
!! Initialize memocc
!! @todo
!!   Should be better for debug purpose to read input.perf before
subroutine bigdft_set_input(radical,posinp,in,atoms)
  use module_base
  use module_types
  use module_interfaces, except_this_one => bigdft_set_input
  use module_input_keys
  use yaml_output
  use dictionaries, only: dictionary
  implicit none

  !Arguments
  character(len=*),intent(in) :: posinp
  character(len=*),intent(in) :: radical
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(out) :: atoms

  character(len=*), parameter :: subname='bigdft_set_input'
  type(dictionary), pointer :: dict
!!$  logical :: exist_list

  atoms=atoms_null()
  ! Read atomic file
  call read_atomic_file(trim(posinp),bigdft_mpi%iproc,atoms%astruct)

  ! initialize mpi environment (this shouldn't be done twice)
!  call mpi_environment_set(bigdft_mpi,iproc,nproc,MPI_COMM_WORLD,nproc)
  dict => read_input_dict_from_files(trim(radical), bigdft_mpi)

  call standard_inputfile_names(in,trim(radical))
  call inputs_from_dict(in, atoms, dict, .true.)
  call dict_free(dict)

  ! Generate the description of input variables.
  if (bigdft_mpi%iproc == 0) then
     call input_keys_dump_def(trim(in%writing_directory) // "/input_help.yaml")
  end if

  ! Read associated pseudo files.
  call init_atomic_values((bigdft_mpi%iproc == 0), atoms, in%ixc)
  call read_atomic_variables(atoms, trim(in%file_igpop),in%nspin)

  ! Start the signaling loop in a thread if necessary.
  if (in%signaling .and. bigdft_mpi%iproc == 0) then
     call bigdft_signals_init(in%gmainloop, 2, in%domain, len(trim(in%domain)))
     call bigdft_signals_start(in%gmainloop, in%signalTimeout)
  end if

  if (bigdft_mpi%iproc == 0) then
     call print_general_parameters(in,atoms)
     !call write_input_parameters(inputs,atoms)
  end if

  !if other steps are supposed to be done leave the last_run to minus one
  !otherwise put it to one
  if (in%last_run == -1 .and. in%ncount_cluster_x <=1 .or. in%ncount_cluster_x <= 1) then
     in%last_run = 1
  end if

END SUBROUTINE bigdft_set_input

!> De-allocate the variable of type input_variables
subroutine bigdft_free_input(in)
  use module_base
  use module_types
  use yaml_output
  type(input_variables), intent(inout) :: in
  
  call free_input_variables(in)
  call f_lib_finalize()
  !free all yaml_streams active
  call yaml_close_all_streams()

end subroutine bigdft_free_input


!> Read the options in the command line using get_command statement
subroutine command_line_information(mpi_groupsize,posinp_file,run_id,ierr)
  use module_types
  implicit none
  integer, intent(out) :: mpi_groupsize
  character(len=*), intent(out) :: posinp_file !< file for list of radicals
  character(len=*), intent(out) :: run_id !< file for radical name
  integer, intent(out) :: ierr !< error code
  !local variables
  integer :: ncommands,icommands
  character(len=256) :: command

  ierr=BIGDFT_SUCCESS
  posinp_file=repeat(' ',len(posinp_file))
  run_id=repeat(' ',len(run_id))
  !traditional scheme
  !if (ncommands == 0) then
     run_id='input'
  !end if

  mpi_groupsize=0
  
  !first see how many arguments are present
  ncommands=COMMAND_ARGUMENT_COUNT()

  do icommands=1,ncommands
     command=repeat(' ',len(command))
     call get_command_argument(icommands,value=command,status=ierr)
     if (ierr /= 0) return
     !print *,'test',ncommands,icommands,command
     call find_command()
     if (ierr /= 0) return
  end do



contains

  subroutine find_command()
    implicit none
    integer :: ipos
    integer, external :: bigdft_error_ret

    if (index(command,'--taskgroup-size=') > 0) then
       if (mpi_groupsize /= 0) then
          ierr=bigdft_error_ret(BIGDFT_INVALID,'taskgroup size specified twice')
       end if
       ipos=index(command,'=')
       read(command(ipos+1:len(command)),*)mpi_groupsize
    else if (index(command,'--run-id=') > 0) then
       if (len_trim(run_id) > 0) then
          ierr=bigdft_error_ret(BIGDFT_INVALID,'run_id specified twice')
       end if
       ipos=index(command,'=')
       read(command(ipos+1:len(command)),*)run_id
    else if (index(command,'--runs-file=') > 0) then
       if (len_trim(posinp_file) > 0 .or. len_trim(run_id) >0) then
          ierr=bigdft_error_ret(BIGDFT_INVALID,'posinp_file specified twice or run_id already known')
       end if
       ipos=index(command,'=')
       read(command(ipos+1:len(command)),*)posinp_file
    else if (index(command,'--') > 0 .and. icommands==1) then
       !help screen
       call help_screen()
       stop
    else if (icommands==1) then
       read(command,*,iostat=ierr)run_id
    else
       call help_screen()
       stop
    end if
  end subroutine find_command

  subroutine help_screen()
    write(*,*)' Usage of the command line instruction'
    write(*,*)' --taskgroup-size=<mpi_groupsize>'
    write(*,*)' --runs-file=<list_posinp filename>'
    write(*,*)' --run-id=<name of the run>: it can be also specified as unique argument'
    write(*,*)' --help : prints this help screen'
  end subroutine help_screen


end subroutine command_line_information


!> Set and check the input file
subroutine set_inputfile(filename, radical, ext)
  implicit none
  character(len = 100), intent(out) :: filename
  character(len = *), intent(in) :: radical, ext
  
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


!> Define the name of the input files
subroutine standard_inputfile_names(in, radical)
  use module_types
  use module_base
  use yaml_output
  implicit none
  type(input_variables), intent(inout) :: in
  character(len = *), intent(in) :: radical

  !set prefix name of the run (as input by defaut for input.dft)
  in%run_name=repeat(' ',len(in%run_name))
  if (trim(radical) /= 'input') in%run_name=trim(radical)

  call set_inputfile(in%file_occnum, radical, "occ")
  call set_inputfile(in%file_igpop, radical,  "occup")
  call set_inputfile(in%file_lin, radical,    "lin")
  call set_inputfile(in%file_frag, radical,   "frag")

  if (trim(radical) == "input") then
     in%dir_output="data" // trim(bigdft_run_id_toa())
  else
     in%dir_output="data-"//trim(radical)!//trim(bigdft_run_id_toa())
  end if

  in%files = INPUTS_NONE
END SUBROUTINE standard_inputfile_names

function read_input_dict_from_files(radical, mpi_env) result(dict)
  use dictionaries
  use wrapper_MPI
  use module_input_keys
  use module_interfaces, only: merge_input_file_to_dict
  use input_old_text_format
  implicit none
  character(len = *), intent(in) :: radical
  type(mpi_environment), intent(in) :: mpi_env

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

  ! We fallback on the old text format
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

!> Read linear input parameters
subroutine lin_input_variables_new(iproc,dump,filename,in,atoms)
  use module_base
  use module_types
  use module_input
  implicit none
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: filename
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  logical, intent(in) :: dump
  !local variables
  logical :: exists
  character(len=*), parameter :: subname='lin_input_variables'
  character(len=256) :: comments
  logical,dimension(atoms%astruct%ntypes) :: parametersSpecified
  logical :: found
  character(len=20):: atomname
  integer :: itype, jtype, ios, ierr, iat, npt, iiorb, iorb, nlr, istat
  real(gp):: ppao, ppl, pph, lrl, lrh, kco
  real(gp),dimension(atoms%astruct%ntypes) :: locradType, locradType_lowaccur, locradType_highaccur

  !Linear input parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Linear Parameters')  
  if (exists) in%files = in%files + INPUTS_LIN

  ! number of accuracy levels: either 2 (for low/high accuracy) or 1 (for hybrid mode)
  comments='number of accuracy levels: either 2 (for low/high accuracy) or 1 (for hybrid mode)'
  call input_var(in%lin%nlevel_accuracy,'2',ranges=(/1,2/),comment=comments)

  ! number of iterations
  comments = 'outer loop iterations (low, high)'
  call input_var(in%lin%nit_lowaccuracy,'15',ranges=(/0,100000/))
  call input_var(in%lin%nit_highaccuracy,'1',ranges=(/0,100000/),comment=comments)

  comments = 'basis iterations (low, high)'
  call input_var(in%lin%nItBasis_lowaccuracy,'12',ranges=(/0,100000/))
  call input_var(in%lin%nItBasis_highaccuracy,'50',ranges=(/0,100000/),comment=comments)

  comments = 'kernel iterations (low, high) - directmin only'
  call input_var(in%lin%nItdmin_lowaccuracy,'1',ranges=(/0,1000/))
  call input_var(in%lin%nItdmin_highaccuracy,'1',ranges=(/0,1000/),comment=comments)

  comments = 'density iterations (low, high)'
  call input_var(in%lin%nItSCCWhenFixed_lowaccuracy,'15',ranges=(/0,1000/))
  call input_var(in%lin%nItSCCWhenFixed_highaccuracy,'15',ranges=(/0,1000/),comment=comments)

  ! DIIS history lengths
  comments = 'DIIS history for basis (low, high)'
  call input_var(in%lin%DIIS_hist_lowaccur,'5',ranges=(/0,100/))
  call input_var(in%lin%DIIS_hist_highaccur,'0',ranges=(/0,100/),comment=comments)

  comments = 'DIIS history for kernel (low, high) - directmin only'
  call input_var(in%lin%dmin_hist_lowaccuracy,'0',ranges=(/0,100/))
  call input_var(in%lin%dmin_hist_highaccuracy,'0',ranges=(/0,100/),comment=comments)

  comments = 'DIIS history for density mixing (low, high)'
  call input_var(in%lin%mixHist_lowaccuracy,'0',ranges=(/0,100/))
  call input_var(in%lin%mixHist_highaccuracy,'0',ranges=(/0,100/),comment=comments)

  ! mixing parameters
  comments = 'density mixing parameter (low, high)'
  call input_var(in%lin%alpha_mix_lowaccuracy,'.5d0',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%alpha_mix_highaccuracy,'.5d0',ranges=(/0.d0,1.d0/),comment=comments)

  ! Convergence criteria
  comments = 'outer loop convergence (low, high)'
  call input_var(in%lin%lowaccuracy_conv_crit,'1.d-8',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%highaccuracy_conv_crit,'1.d-12',ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'basis convergence (low, high)'
  call input_var(in%lin%convCrit_lowaccuracy,'1.d-3',ranges=(/0.0_gp,1.0_gp/))
  call input_var(in%lin%convCrit_highaccuracy,'1.d-5',ranges=(/0.0_gp,1.0_gp/),comment=comments)
  
  comments = 'multiplier to (exit one TMB optimization, fix TMB completely). Only used for hybrid mode'
  call input_var(in%lin%deltaenergy_multiplier_TMBexit,'1.d0',ranges=(/1.d-5,1.d1/))
  call input_var(in%lin%deltaenergy_multiplier_TMBfix,'1.d0',ranges=(/1.d-5,1.d1/),comment=comments)

  comments = 'factor to reduce the confinement. Only used for hybrid mode.'
  call input_var(in%lin%reduce_confinement_factor,'0.5d0',ranges=(/-1.d100,1.d0/),comment=comments)

  comments = 'kernel convergence (low, high) - directmin only'
  call input_var(in%lin%convCritdmin_lowaccuracy,'1.d-5',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%convCritdmin_highaccuracy,'1.d-5',ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'density convergence (low, high)'
  call input_var(in%lin%convCritMix_lowaccuracy,'1.d-13',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%convCritMix_highaccuracy,'1.d-13',ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'convergence criterion on density to fix TMBS'
  call input_var(in%lin%support_functions_converged,'1.d-10',ranges=(/0.d0,1.d0/),comment=comments)

  ! Miscellaneous
  comments='mixing method: 100 (direct minimization), 101 (simple dens mixing), 102 (simple pot mixing), 103 (FOE)'
  call input_var(in%lin%scf_mode,'100',ranges=(/100,103/),comment=comments)

  comments = 'initial step size for basis optimization (DIIS, SD)' ! DELETE ONE
  call input_var(in%lin%alphaDIIS,'1.d0',ranges=(/0.0_gp,10.0_gp/))
  call input_var(in%lin%alphaSD,'1.d0',ranges=(/0.0_gp,10.0_gp/),comment=comments)

  comments = 'initial step size for kernel update (SD), curve fitting for alpha update - directmin only'
  call input_var(in%lin%alphaSD_coeff,'1.d0',ranges=(/0.0_gp,10.0_gp/))
  call input_var(in%lin%curvefit_dmin,'F',comment=comments)

  comments = 'lower and upper bound for the eigenvalue spectrum (FOE). Will be adjusted automatically if chosen too small'
  call input_var(in%lin%evlow,'-.5d0',ranges=(/-10.d0,-1.d-10/))
  call input_var(in%lin%evhigh,'-.5d0',ranges=(/1.d-10,10.d0/),comment=comments)

  comments='number of iterations in the preconditioner'
  call input_var(in%lin%nItPrecond,'5',ranges=(/1,100/),comment=comments)
  
  comments = '0-> exact Loewdin, 1-> taylor expansion; &
             &in orthoconstraint: correction for non-orthogonality (0) or no correction (1)'
  call input_var(in%lin%methTransformOverlap,'1',ranges=(/-1,1/))
  call input_var(in%lin%correctionOrthoconstraint,'1',ranges=(/0,1/),comment=comments)

  comments='fscale: length scale over which complementary error function decays from 1 to 0'
  call input_var(in%lin%fscale,'1.d-2',ranges=(/0.d0,1.d0/),comment=comments)

  !plot basis functions: true or false
  comments='Output basis functions: 0 no output, 1 formatted output, 2 Fortran bin, 3 ETSF ;'//&
           'calculate dipole ; pulay correction; diagonalization at the end (dmin, FOE)'
  call input_var(in%lin%plotBasisFunctions,'0',ranges=(/0,3/))
  call input_var(in%lin%calc_dipole,'F')
  call input_var(in%lin%pulay_correction,'T')
  call input_var(in%lin%diag_end,'F',comment=comments)

  !fragment calculation and transfer integrals: true or false
  comments='fragment calculation; calculate transfer_integrals; constrained DFT calculation; extra states to optimize (dmin only)'
  call input_var(in%lin%fragment_calculation,'F')
  call input_var(in%lin%calc_transfer_integrals,'F')
  call input_var(in%lin%constrained_dft,'F')
  call input_var(in%lin%extra_states,'0',ranges=(/0,10000/),comment=comments)

  ! Allocate lin pointers and atoms%rloc
  call nullifyInputLinparameters(in%lin)
  call allocateBasicArraysInputLin(in%lin, atoms%astruct%ntypes)
  
  ! Now read in the parameters specific for each atom type.
  comments = 'Atom name, number of basis functions per atom, prefactor for confinement potential,'//&
             'localization radius, kernel cutoff'
  parametersSpecified=.false.
  itype = 1
  do
     !Check at the beginning to permit natom=0
     if (itype > atoms%astruct%ntypes) exit
     if (exists) then
        call input_var(atomname,'C',input_iostat=ios)
        if (ios /= 0) exit
     else
        call input_var(atomname,trim(atoms%astruct%atomnames(itype)))
        itype = itype + 1
     end if
     call input_var(npt,'1',ranges=(/1,100/),input_iostat=ios)
     call input_var(ppao,'1.2d-2',ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(ppl,'1.2d-2',ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(pph,'5.d-5',ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(lrl,'10.d0',ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios)
     call input_var(lrh,'10.d0',ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios)
     call input_var(kco,'20.d0',ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios,comment=comments)
     ! The reading was successful. Check whether this atom type is actually present.
     found=.false.
     do jtype=1,atoms%astruct%ntypes
        if(trim(atomname)==trim(atoms%astruct%atomnames(jtype))) then
           found=.true.
           parametersSpecified(jtype)=.true.
           in%lin%norbsPerType(jtype)=npt
           in%lin%potentialPrefac_ao(jtype)=ppao
           in%lin%potentialPrefac_lowaccuracy(jtype)=ppl
           in%lin%potentialPrefac_highaccuracy(jtype)=pph
           locradType(jtype)=lrl
           in%lin%locrad_type(jtype)=lrl
           locradType_lowaccur(jtype)=lrl
           locradType_highaccur(jtype)=lrh
           atoms%rloc(jtype,:)=locradType(jtype)
           in%lin%kernel_cutoff(jtype)=kco
        end if
     end do
     if(.not.found) then
        if(iproc==0 .and. dump) write(*,'(1x,3a)') "WARNING: you specified informations about the atomtype '",trim(atomname), &
             "', which is not present in the file containing the atomic coordinates."
     end if
  end do
  found  = .true.
  do jtype=1,atoms%astruct%ntypes
     found = found .and. parametersSpecified(jtype)
  end do
  if (.not. found) then
     ! The parameters were not specified for all atom types.
     if(iproc==0) then
        write(*,'(1x,a)',advance='no') "ERROR: the file 'input.lin' does not contain the parameters&
             & for the following atom types:"
        do jtype=1,atoms%astruct%ntypes
           if(.not.parametersSpecified(jtype)) write(*,'(1x,a)',advance='no') trim(atoms%astruct%atomnames(jtype))
        end do
     end if
     call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
     stop
  end if

  nlr=0
  do iat=1,atoms%astruct%nat
      itype=atoms%astruct%iatype(iat)
      nlr=nlr+in%lin%norbsPerType(itype)
  end do
  allocate(in%lin%locrad(nlr),stat=istat)
  call memocc(istat,in%lin%locrad,'in%lin%locrad',subname)
  allocate(in%lin%locrad_lowaccuracy(nlr),stat=istat)
  call memocc(istat,in%lin%locrad_lowaccuracy,'in%lin%locrad_lowaccuracy',subname)
  allocate(in%lin%locrad_highaccuracy(nlr),stat=istat)
  call memocc(istat,in%lin%locrad_highaccuracy,'in%lin%locrad_highaccuracy',subname)

  
  ! Assign the localization radius to each atom.
  iiorb=0
  do iat=1,atoms%astruct%nat
      itype=atoms%astruct%iatype(iat)
      do iorb=1,in%lin%norbsPerType(itype)
          iiorb=iiorb+1
          in%lin%locrad(iiorb)=locradType(itype)
          in%lin%locrad_lowaccuracy(iiorb)=locradType_lowaccur(itype)
          in%lin%locrad_highaccuracy(iiorb)=locradType_highaccur(itype)
      end do
  end do
  

  call input_free((iproc == 0) .and. dump)
END SUBROUTINE lin_input_variables_new


!> Assign default values for TDDFT variables
subroutine tddft_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in

  in%tddft_approach='NONE'

END SUBROUTINE tddft_input_variables_default

!> Read fragment input parameters
subroutine fragment_input_variables(iproc,dump,filename,in,atoms)
  use module_base
  use module_types
  use module_input
  implicit none
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: filename
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  logical, intent(in) :: dump
  !local variables
  logical :: exists
  character(len=*), parameter :: subname='fragment_input_variables'
  character(len=256) :: comments
  integer :: ifrag, frag_num, ierr

  !Linear input parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Fragment Parameters') 
  if (exists .and. dump) in%files = in%files + INPUTS_FRAG

  if (.not. exists .and. in%lin%fragment_calculation) then ! we should be doing a fragment calculation, so this is a problem
     write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' is missing and fragment calculation was specified"
     call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
     stop
  end if

  ! number of reference fragments
  comments='# number of fragments in reference system, number of fragments in current system'
  call input_var(in%frag%nfrag_ref,'1',ranges=(/1,100000/))
  call input_var(in%frag%nfrag,'1',ranges=(/1,100000/),comment=comments)
  
  ! Allocate fragment pointers
  call nullifyInputFragParameters(in%frag)
  call allocateInputFragArrays(in%frag)

  !comments = '# reference fragment number i, number of atoms in reference fragment i, '//&
  !           'number of atoms in corresponding environment'
  !do ifrag=1,in%frag%nfrag_ref
  !  call input_var(frag_num,'1',ranges=(/1,in%frag%nfrag_ref/))
  !  if (frag_num/=ifrag) then
  !      write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' has an error when specifying&
  !           & the reference fragments"
  !     call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  !     stop
  !  end if
  !  call input_var(in%frag%frag_info(frag_num,1),'1',ranges=(/1,100000/))
  !  call input_var(in%frag%frag_info(frag_num,2),'0',ranges=(/0,100000/),comment=comments)
  !end do

  ! ADD A SENSIBLE DEFAULT AND ALLOW FOR USER NOT TO SPECIFY FRAGMENT NAMES
  comments = '#  reference fragment number i, fragment label'
  do ifrag=1,in%frag%nfrag_ref
    call input_var(frag_num,'1',ranges=(/1,in%frag%nfrag_ref/))
    if (frag_num/=ifrag) then
        write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' has an error when specifying&
             & the reference fragments"
       call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
       stop
    end if
    call input_var(in%frag%label(frag_num),' ',comment=comments)
    in%frag%label(frag_num)=trim(in%frag%label(frag_num))
    ! keep dirname blank if this isn't a fragment calculation
    if (len(trim(in%frag%label(frag_num)))>=1) then
       in%frag%dirname(frag_num)='data-'//trim(in%frag%label(frag_num))//'/'
    else
       in%frag%dirname(frag_num)=''
    end if
  end do

  comments = '# fragment number j, reference fragment i this corresponds to, charge on this fragment'
  do ifrag=1,in%frag%nfrag
    call input_var(frag_num,'1',ranges=(/1,in%frag%nfrag/))
    if (frag_num/=ifrag) then
        write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' has an error when specifying&
             & the system fragments"
       call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
       stop
    end if
    call input_var(in%frag%frag_index(frag_num),'1',ranges=(/0,100000/))
    call input_var(in%frag%charge(frag_num),'1',ranges=(/-500,500/),comment=comments)
  end do

  call input_free((iproc == 0) .and. dump)

END SUBROUTINE fragment_input_variables

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


!> Read the input variables needed for the ABSCALC
!! Every argument should be considered as mandatory
subroutine abscalc_input_variables(iproc,filename,in)
  use module_base
  use module_types
  implicit none
  !Arguments
  type(input_variables), intent(inout) :: in
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  !Local variables
  integer, parameter :: iunit = 112
  integer :: ierror,iline, i

  character(len=*), parameter :: subname='abscalc_input_variables'
  integer :: i_stat

  ! Read the input variables.
  open(unit=iunit,file=filename,status='old')

  !line number, to control the input values
  iline=0

  !x-absorber treatment (in progress)

  read(iunit,*,iostat=ierror) in%iabscalc_type
  call check()


  read(iunit,*,iostat=ierror)  in%iat_absorber
  call check()
  read(iunit,*,iostat=ierror)  in%L_absorber
  call check()

  allocate(in%Gabs_coeffs(2*in%L_absorber +1+ndebug),stat=i_stat)
  call memocc(i_stat,in%Gabs_coeffs,'Gabs_coeffs',subname)

  read(iunit,*,iostat=ierror)  (in%Gabs_coeffs(i), i=1,2*in%L_absorber +1 )
  call check()

  read(iunit,*,iostat=ierror)  in%potshortcut
  call check()
  
  read(iunit,*,iostat=ierror)  in%nsteps
  call check()

  if( iand( in%potshortcut,4)>0) then
     read(iunit,'(a100)',iostat=ierror) in%extraOrbital
  end if
  
  read(iunit,*,iostat=ierror) in%abscalc_bottomshift
  if(ierror==0) then
  else
     in%abscalc_bottomshift=0
  endif

 

  read(iunit, '(a100)' ,iostat=ierror) in%xabs_res_prefix
  if(ierror==0) then
  else
     in%xabs_res_prefix=""
  endif


  read(iunit,*,iostat=ierror) in%abscalc_alterpot, in%abscalc_eqdiff 
  !!, &
  !!     in%abscalc_S_do_cg ,in%abscalc_Sinv_do_cg
  if(ierror==0) then
  else
     in%abscalc_alterpot=.false.
     in%abscalc_eqdiff =.false.
  endif



  in%c_absorbtion=.true.

  close(unit=iunit)

contains

  subroutine check()
    iline=iline+1
    if (ierror/=0) then
       if (iproc == 0) write(*,'(1x,a,a,a,i3)') &
            'Error while reading the file "',trim(filename),'", line=',iline
       stop
    end if
  END SUBROUTINE check

END SUBROUTINE abscalc_input_variables


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


!> Read the input variables needed for the frequencies calculation.
!! Every argument should be considered as mandatory.
subroutine frequencies_input_variables_new(iproc,dump,filename,in)
  use module_base
  use module_types
  use module_input
  implicit none
  !Arguments
  type(input_variables), intent(inout) :: in
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  logical, intent(in) :: dump
  !Local variables
  logical :: exists
  !n(c) integer, parameter :: iunit=111

  !Frequencies parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Frequencies Parameters')  
  if (exists) in%files = in%files + INPUTS_FREQ
  !call the variable, its default value, the line ends if there is a comment

  !Read in%freq_alpha (possible 1/64)
  call input_var(in%freq_alpha,'1/64',ranges=(/0.0_gp,1.0_gp/),&
       comment="Step size factor (alpha*hgrid)")
  !Read the order of finite difference scheme

  call input_var(in%freq_order,'2',exclusive=(/-1,1,2,3/),&
       comment="Order of the difference scheme")
  !Read the index of the method

  call input_var(in%freq_method,'1',exclusive=(/1/),&
       comment="Method used (only possible value=1)")
  call input_free((iproc == 0) .and. dump)

END SUBROUTINE frequencies_input_variables_new


!> Fill the arrays occup and spinsgn
!! if iunit /=0 this means that the file 'input.occ' does exist and it opens
subroutine occupation_input_variables(verb,iunit,nelec,norb,norbu,norbuempty,norbdempty,nspin,occup,spinsgn)
  use module_base
  use module_input
  use yaml_output
  implicit none
  ! Arguments
  logical, intent(in) :: verb
  integer, intent(in) :: nelec,nspin,norb,norbu,iunit,norbuempty,norbdempty
  real(gp), dimension(norb), intent(out) :: occup,spinsgn
  ! Local variables
  integer :: iorb,nt,ne,it,ierror,iorb1,i
  real(gp) :: rocc
  character(len=20) :: string
  character(len=100) :: line

  do iorb=1,norb
     spinsgn(iorb)=1.0_gp
  end do
  if (nspin/=1) then
     do iorb=1,norbu
        spinsgn(iorb)=1.0_gp
     end do
     do iorb=norbu+1,norb
        spinsgn(iorb)=-1.0_gp
     end do
  end if
  ! write(*,'(1x,a,5i4,30f6.2)')'Spins: ',norb,norbu,norbd,norbup,norbdp,(spinsgn(iorb),iorb=1,norb)

  ! First fill the occupation numbers by default
  nt=0
  if (nspin==1) then
     ne=(nelec+1)/2
     do iorb=1,ne
        it=min(2,nelec-nt)
        occup(iorb)=real(it,gp)
        nt=nt+it
     enddo
     do iorb=ne+1,norb
        occup(iorb)=0._gp
     end do
  else
     if (norbuempty+norbdempty == 0) then
        if (norb > nelec) then
           do iorb=1,min(norbu,norb/2+1)
              it=min(1,nelec-nt)
              occup(iorb)=real(it,gp)
              nt=nt+it
           enddo
           do iorb=min(norbu,norb/2+1)+1,norbu
              occup(iorb)=0.0_gp
           end do
           do iorb=norbu+1,norbu+min(norb-norbu,norb/2+1)
              it=min(1,nelec-nt)
              occup(iorb)=real(it,gp)
              nt=nt+it
           enddo
           do iorb=norbu+min(norb-norbu,norb/2+1)+1,norb
              occup(iorb)=0.0_gp
           end do
        else
           do iorb=1,norb
              occup(iorb)=1.0_gp
           end do
        end if
     else
        do iorb=1,norbu-norbuempty
           occup(iorb)=1.0_gp
        end do
        do iorb=norbu-norbuempty+1,norbu
           occup(iorb)=0.0_gp
        end do
        do iorb=1,norb-norbu-norbdempty
           occup(norbu+iorb)=1.0_gp
        end do
        do iorb=norb-norbu-norbdempty+1,norb-norbu
           occup(norbu+iorb)=0.0_gp
        end do
     end if
  end if
  ! Then read the file "input.occ" if does exist
  if (iunit /= 0) then
     nt=0
     do
        read(unit=iunit,fmt='(a100)',iostat=ierror) line
        if (ierror /= 0) then
           exit
        end if
        !Transform the line in case there are slashes (to ease the parsing)
        do i=1,len(line)
           if (line(i:i) == '/') then
              line(i:i) = ':'
           end if
        end do
        read(line,*,iostat=ierror) iorb,string
        call read_fraction_string(string,rocc,ierror) 
        if (ierror /= 0) then
           exit
        end if

        if (ierror/=0) then
           exit
        else
           nt=nt+1
           if (iorb<0 .or. iorb>norb) then
              !if (iproc==0) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "[name].occ"'
              write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
              !end if
              stop
           elseif (rocc<0._gp .or. rocc>2._gp) then
              !if (iproc==0) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "[name].occ"'
              write(*,'(10x,a,f5.2,a)') 'The occupation number ',rocc,' is not between 0. and 2.'
              !end if
              stop
           else
              occup(iorb)=rocc
           end if
        end if
     end do
     if (verb) then
        call yaml_comment('('//adjustl(trim(yaml_toa(nt)))//'lines read)')
        !write(*,'(1x,a,i0,a)') &
        !     'The occupation numbers are read from the file "[name].occ" (',nt,' lines read)'
     end if
     close(unit=iunit)

     if (nspin/=1) then
!!!        !Check if the polarisation is respected (mpol)
!!!        rup=sum(occup(1:norbu))
!!!        rdown=sum(occup(norbu+1:norb))
!!!        if (abs(rup-rdown-real(norbu-norbd,gp))>1.e-6_gp) then
!!!           if (iproc==0) then
!!!              write(*,'(1x,a,f13.6,a,i0)') 'From the file "input.occ", the polarization ',rup-rdown,&
!!!                             ' is not equal to ',norbu-norbd
!!!           end if
!!!           stop
!!!        end if
        !Fill spinsgn
        do iorb=1,norbu
           spinsgn(iorb)=1.0_gp
        end do
        do iorb=norbu+1,norb
           spinsgn(iorb)=-1.0_gp
        end do
     end if
  end if
  if (verb) then 
     call yaml_sequence(advance='no')
     call yaml_open_map('Occupation Numbers',flow=.true.)
     !write(*,'(1x,a,t28,i8)') 'Total Number of Orbitals',norb
     iorb1=1
     rocc=occup(1)
     do iorb=1,norb
        if (occup(iorb) /= rocc) then
           if (iorb1 == iorb-1) then
              call yaml_map('Orbital No.'//trim(yaml_toa(iorb1)),rocc,fmt='(f6.4)')
              !write(*,'(1x,a,i0,a,f6.4)') 'occup(',iorb1,')= ',rocc
           else
           call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
                adjustl(trim(yaml_toa(iorb-1))),rocc,fmt='(f6.4)')
           !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',iorb-1,')= ',rocc
           end if
           rocc=occup(iorb)
           iorb1=iorb
        end if
     enddo
     if (iorb1 == norb) then
        call yaml_map('Orbital No.'//trim(yaml_toa(norb)),occup(norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,f6.4)') 'occup(',norb,')= ',occup(norb)
     else
        call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
             adjustl(trim(yaml_toa(norb))),occup(norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',norb,')= ',occup(norb)
     end if
     call yaml_close_map()
  endif

  !Check if sum(occup)=nelec
  rocc=sum(occup)
  if (abs(rocc-real(nelec,gp))>1.e-6_gp) then
     call yaml_warning('ERROR in determining the occupation numbers: the total number of electrons ' &
        & // trim(yaml_toa(rocc,fmt='(f13.6)')) // ' is not equal to' // trim(yaml_toa(nelec)))
     !if (iproc==0) then
     !write(*,'(1x,a,f13.6,a,i0)') 'ERROR in determining the occupation numbers: the total number of electrons ',rocc,&
     !     ' is not equal to ',nelec
     !end if
     stop
  end if

END SUBROUTINE occupation_input_variables


module position_files
   implicit none
   contains
   subroutine directGetLine(line, ifile, eof)
      !Arguments
      integer, intent(in) :: ifile
      character(len=150), intent(out) :: line
      logical, intent(out) :: eof
      !Local variables
      integer :: i_stat

      eof = .false.
      read(ifile,'(a150)', iostat = i_stat) line
      if (i_stat /= 0) eof = .true.
   END SUBROUTINE directGetLine

   subroutine archiveGetLine(line, ifile, eof)
      !Arguments
      integer, intent(in) :: ifile
      character(len=150), intent(out) :: line
      logical, intent(out) :: eof
      !Local variables
      integer :: i_stat
      !The argument ifile is not used but it is used as argument routine
      !eof = .false.
      eof = (ifile /= ifile)
      call extractNextLine(line, i_stat)
      if (i_stat /= 0) eof = .true.
   END SUBROUTINE archiveGetLine
end module position_files

!> Read atomic file
subroutine read_atomic_file(file,iproc,astruct,status,comment,energy,fxyz)
   use module_base
   use module_types
   use module_interfaces, except_this_one => read_atomic_file
   use m_ab6_symmetry
   use position_files
   implicit none
   character(len=*), intent(in) :: file
   integer, intent(in) :: iproc
   type(atomic_structure), intent(inout) :: astruct
   integer, intent(out), optional :: status
   real(gp), intent(out), optional :: energy
   real(gp), dimension(:,:), pointer, optional :: fxyz
   character(len = *), intent(out), optional :: comment
   !Local variables
   character(len=*), parameter :: subname='read_atomic_file'
   integer :: l, extract, i_all, i_stat
   logical :: file_exists, archive
   character(len = 128) :: filename
   character(len = 15) :: arFile
   character(len = 6) :: ext
   real(gp) :: energy_
   real(gp), dimension(:,:), pointer :: fxyz_
   character(len = 1024) :: comment_

   file_exists = .false.
   archive = .false.
   if (present(status)) status = 0
   nullify(fxyz_)

   ! Extract from archive
   if (index(file, "posout_") == 1 .or. index(file, "posmd_") == 1) then
      write(arFile, "(A)") "posout.tar.bz2"
      if (index(file, "posmd_") == 1) write(arFile, "(A)") "posmd.tar.bz2"
      inquire(FILE = trim(arFile), EXIST = file_exists)
      if (file_exists) then
         !!$     call extractNextCompress(trim(arFile), len(trim(arFile)), &
         !!$          & trim(file), len(trim(file)), extract, ext)
         call openNextCompress(trim(arFile), len(trim(arFile)), &
         & trim(file), len(trim(file)), extract, ext)
         if (extract == 0) then
            write(*,*) "Can't find '", file, "' in archive."
            if (present(status)) then
               status = 1
               return
            else
               stop
            end if
         end if
         archive = .true.
         write(filename, "(A)") file//'.'//trim(ext)
         write(astruct%inputfile_format, "(A)") trim(ext)
      end if
   end if

   ! Test posinp.xyz
   if (.not. file_exists) then
      inquire(FILE = file//'.xyz', EXIST = file_exists)
      if (file_exists) then
         write(filename, "(A)") file//'.xyz'!"posinp.xyz"
         write(astruct%inputfile_format, "(A)") "xyz"
         open(unit=99,file=trim(filename),status='old')
      end if
   end if
   ! Test posinp.ascii
   if (.not. file_exists) then
      inquire(FILE = file//'.ascii', EXIST = file_exists)
      if (file_exists) then
         write(filename, "(A)") file//'.ascii'!"posinp.ascii"
         write(astruct%inputfile_format, "(A)") "ascii"
         open(unit=99,file=trim(filename),status='old')
      end if
   end if
   ! Test posinp.yaml
   if (.not. file_exists) then
      inquire(FILE = file//'.yaml', EXIST = file_exists)
      if (file_exists) then
         write(filename, "(A)") file//'.yaml'!"posinp.ascii"
         write(astruct%inputfile_format, "(A)") "yaml"
      end if
   end if
   ! Test the name directly
   if (.not. file_exists) then
      inquire(FILE = file, EXIST = file_exists)
      if (file_exists) then
         write(filename, "(A)") file
         l = len(file)
         if (file(l-3:l) == ".xyz") then
            write(astruct%inputfile_format, "(A)") "xyz"
         else if (file(l-5:l) == ".ascii") then
            write(astruct%inputfile_format, "(A)") "ascii"
         else if (file(l-4:l) == ".yaml") then
            write(astruct%inputfile_format, "(A)") "yaml"
         else
            write(*,*) "Atomic input file '" // trim(file) // "', format not recognised."
            write(*,*) " File should be *.yaml, *.ascii or *.xyz."
            if (present(status)) then
               status = 1
               return
            else
               stop
            end if
         end if
         if (trim(astruct%inputfile_format) /= "yaml") then
            open(unit=99,file=trim(filename),status='old')
         end if
      end if
   end if

   if (.not. file_exists) then
      if (present(status)) then
         status = 1
         return
      else
         write(*,*) "Atomic input file not found."
         write(*,*) " Files looked for were '"//file//".yaml', '"//file//".ascii', '"//file//".xyz' and '"//file//"'."
         stop 
      end if
   end if

   if (astruct%inputfile_format == "xyz") then
      !read atomic positions
      if (.not.archive) then
         call read_xyz_positions(iproc,99,astruct,comment_,energy_,fxyz_,directGetLine)
      else
         call read_xyz_positions(iproc,99,astruct,comment_,energy_,fxyz_,archiveGetLine)
      end if
   else if (astruct%inputfile_format == "ascii") then
      i_stat = iproc
      if (present(status)) i_stat = 1
      !read atomic positions
      if (.not.archive) then
         call read_ascii_positions(i_stat,99,astruct,comment_,energy_,fxyz_,directGetLine)
      else
         call read_ascii_positions(i_stat,99,astruct,comment_,energy_,fxyz_,archiveGetLine)
      end if
   else if (astruct%inputfile_format == "yaml") then
      !read atomic positions
      if (.not.archive) then
         call read_yaml_positions(trim(filename),astruct,comment_,energy_,fxyz_)
      else
         write(*,*) "Atomic input file in YAML not yet supported in archive file."
         stop
      end if
   end if

   !Check the number of atoms
   if (astruct%nat < 0) then
      if (present(status)) then
         status = 1
         return
      else
         write(*,'(1x,3a,i0,a)') "In the file '",trim(filename),&
              &  "', the number of atoms (",astruct%nat,") < 0 (should be >= 0)."
         stop 
      end if
   end if

   !control atom positions
   call check_atoms_positions(astruct,(iproc == 0))

   ! We delay the calculation of the symmetries.
!this should be already in the atoms_null routine
   astruct%sym=symm_null()
!   astruct%sym%symObj = -1
!   nullify(astruct%sym%irrzon)
!   nullify(astruct%sym%phnons)

   ! close open file.
   if (.not.archive .and. trim(astruct%inputfile_format) /= "yaml") then
      close(99)
      !!$  else
      !!$     call unlinkExtract(trim(filename), len(trim(filename)))
   end if
   
   ! We transfer optionals.
   if (present(energy)) then
      energy = energy_
   end if
   if (present(comment)) then
      write(comment, "(A)") comment_
   end if
   if (present(fxyz)) then
      fxyz => fxyz_
   else if (associated(fxyz_)) then
      i_all=-product(shape(fxyz_))*kind(fxyz_)
      deallocate(fxyz_,stat=i_stat)
      call memocc(i_stat,i_all,'fxyz_',subname)
   end if
END SUBROUTINE read_atomic_file

!> Write an atomic file
!Yaml output included
subroutine write_atomic_file(filename,energy,rxyz,atoms,comment,forces)
  use module_base
  use module_types
  use yaml_output
  implicit none
  character(len=*), intent(in) :: filename,comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3,atoms%astruct%nat), intent(in), optional :: forces
  !local variables
  character(len = 15) :: arFile
  integer :: iunit

  if (trim(filename) == "stdout") then
     iunit = 6
  else
     open(unit=9,file=trim(filename)//'.'//trim(atoms%astruct%inputfile_format))
     iunit = 9
  end if
  if (atoms%astruct%inputfile_format == "xyz") then
     call wtxyz(iunit,energy,rxyz,atoms,comment)
     if (present(forces)) call wtxyz_forces(9,forces,atoms)
  else if (atoms%astruct%inputfile_format == "ascii") then
     call wtascii(iunit,energy,rxyz,atoms,comment)
     if (present(forces)) call wtascii_forces(9,forces,atoms)
  else if (atoms%astruct%inputfile_format == 'yaml') then
     if (present(forces)) then
        call wtyaml(iunit,energy,rxyz,atoms,comment,.true.,forces)
     else
        call wtyaml(iunit,energy,rxyz,atoms,comment,.false.,rxyz)
     end if
  else
     write(*,*) "Error, unknown file format."
     stop
  end if
  if (trim(filename) /= "stdout") then
     close(unit=9)
  end if

  ! Add to archive
  if (index(filename, "posout_") == 1 .or. index(filename, "posmd_") == 1) then
     write(arFile, "(A)") "posout.tar.bz2"
     if (index(filename, "posmd_") == 1) write(arFile, "(A)") "posmd.tar.bz2"
     call addToCompress(trim(arFile), len(trim(arFile)), &
          & trim(filename)//'.'//trim(atoms%astruct%inputfile_format), &
          & len(trim(filename)//'.'//trim(atoms%astruct%inputfile_format)))
  end if
END SUBROUTINE write_atomic_file

!>Calculate the coefficient for moving atoms following the ifrztyp
subroutine frozen_alpha(ifrztyp,ixyz,alpha,alphai)
  use module_base
  implicit none
  integer, intent(in) :: ifrztyp,ixyz
  real(gp), intent(in) :: alpha
  real(gp), intent(out) :: alphai
  !local variables
  logical :: move_this_coordinate

  if (move_this_coordinate(ifrztyp,ixyz)) then
     alphai=alpha
  else
     alphai=0.0_gp
  end if
 
END SUBROUTINE frozen_alpha

!>Routine for moving atomic positions, takes into account the 
!!   frozen atoms and the size of the cell
!!   synopsis: rxyz=txyz+alpha*sxyz
!!   all the shift are inserted into the box if there are periodic directions
!!   if the atom are frozen they are not moved
subroutine atomic_axpy(atoms,txyz,alpha,sxyz,rxyz)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: alpha
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz

  do iat=1,atoms%astruct%nat
     !adjust the moving of the atoms following the frozen direction
     call frozen_alpha(atoms%astruct%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),3,alpha,alphaz)

     if (atoms%astruct%geocode == 'P') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(2,iat)=modulo(txyz(2,iat)+alphay*sxyz(2,iat),atoms%astruct%cell_dim(2))
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),atoms%astruct%cell_dim(3))
     else if (atoms%astruct%geocode == 'S') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),atoms%astruct%cell_dim(3))
     else
        rxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
     end if
  end do

END SUBROUTINE atomic_axpy


!>Routine for moving atomic positions, takes into account the 
!!   frozen atoms and the size of the cell
!!   synopsis: fxyz=txyz+alpha*sxyz
!!   update the forces taking into account the frozen atoms
!!   do not apply the modulo operation on forces 
subroutine atomic_axpy_forces(atoms,txyz,alpha,sxyz,fxyz)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: alpha
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: fxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz
  
  do iat=1,atoms%astruct%nat
     !adjust the moving of the forces following the frozen direction
     call frozen_alpha(atoms%astruct%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),3,alpha,alphaz)

     fxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
     fxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
     fxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
  end do
  
END SUBROUTINE atomic_axpy_forces


!>Calculate the scalar product between atomic positions by considering
!!   only non-blocked atoms
subroutine atomic_dot(atoms,x,y,scpr)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: x,y
  real(gp), intent(out) :: scpr
  !local variables
  integer :: iat
  real(gp) :: scpr1,scpr2,scpr3
  real(gp) :: alphax,alphay,alphaz

  scpr=0.0_gp

  do iat=1,atoms%astruct%nat
     call frozen_alpha(atoms%astruct%ifrztyp(iat),1,1.0_gp,alphax)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),2,1.0_gp,alphay)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),3,1.0_gp,alphaz)
     scpr1=alphax*x(1,iat)*y(1,iat)
     scpr2=alphay*x(2,iat)*y(2,iat)
     scpr3=alphaz*x(3,iat)*y(3,iat)
     scpr=scpr+scpr1+scpr2+scpr3
  end do
  
END SUBROUTINE atomic_dot


!>z=alpha*A*x + beta* y
subroutine atomic_gemv(atoms,m,alpha,A,x,beta,y,z)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: m
  real(gp), intent(in) :: alpha,beta
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: x
  real(gp), dimension(m), intent(in) :: y
  real(gp), dimension(m,3,atoms%astruct%nat), intent(in) :: A
  real(gp), dimension(m), intent(out) :: z
  !local variables
  integer :: iat,i,j
  real(gp) :: mv,alphai
  
  do i=1,m
     mv=0.0_gp
     do iat=1,atoms%astruct%nat
        do j=1,3
           call frozen_alpha(atoms%astruct%ifrztyp(iat),j,A(i,j,iat),alphai)
           mv=mv+alphai*x(j,iat)
        end do
     end do
     z(i)=alpha*mv+beta*y(i)
  end do

END SUBROUTINE atomic_gemv


!>  The function which controls all the moving positions
function move_this_coordinate(ifrztyp,ixyz)
  use module_base
  implicit none
  integer, intent(in) :: ixyz,ifrztyp
  logical :: move_this_coordinate
  
  move_this_coordinate= &
       ifrztyp == 0 .or. &
       (ifrztyp == 2 .and. ixyz /=2) .or. &
       (ifrztyp == 3 .and. ixyz ==2)
       
END FUNCTION move_this_coordinate


!> rxyz=txyz+alpha*sxyz
subroutine atomic_coordinate_axpy(atoms,ixyz,iat,t,alphas,r)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ixyz,iat
  real(gp), intent(in) :: t,alphas
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(out) :: r
  !local variables
  logical :: periodize
  real(gp) :: alat,alphai

  if (ixyz == 1) then
     alat=atoms%astruct%cell_dim(1)
  else if (ixyz == 2) then
     alat=atoms%astruct%cell_dim(2)
  else if (ixyz == 3) then
     alat=atoms%astruct%cell_dim(3)
  else
     alat = -1
     write(0,*) "Internal error"
     stop
  end if
  
  periodize= atoms%astruct%geocode == 'P' .or. &
       (atoms%astruct%geocode == 'S' .and. ixyz /= 2)

  call frozen_alpha(atoms%astruct%ifrztyp(iat),ixyz,alphas,alphai)

  if (periodize) then
     r=modulo(t+alphai,alat)
  else
     r=t+alphai
  end if

END SUBROUTINE atomic_coordinate_axpy


!> Initialization of acceleration (OpenCL)
subroutine init_material_acceleration(iproc,matacc,GPU)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in):: iproc
  type(material_acceleration), intent(in) :: matacc
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  integer :: iconv,iblas,initerror,ierror,useGPU,mproc,ierr,nproc_node

  if (matacc%iacceleration == 1) then
     call MPI_COMM_SIZE(bigdft_mpi%mpi_comm,mproc,ierr)
     !initialize the id_proc per node
     call processor_id_per_node(iproc,mproc,GPU%id_proc,nproc_node)
     call sg_init(GPUshare,useGPU,iproc,nproc_node,initerror)
     if (useGPU == 1) then
        iconv = 1
        iblas = 1
     else
        iconv = 0
        iblas = 0
     end if
     if (initerror == 1) then
        call yaml_warning('(iproc=' // trim(yaml_toa(iproc,fmt='(i0)')) // &
        &    ') S_GPU library init failed, aborting...')
        !write(*,'(1x,a)')'**** ERROR: S_GPU library init failed, aborting...'
        call MPI_ABORT(bigdft_mpi%mpi_comm,initerror,ierror)
     end if

     if (iconv == 1) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUconv=.true.
     end if
     if (iblas == 1) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUblas=.true.
     end if

     if (iproc == 0) then
        call yaml_map('Material acceleration','CUDA',advance='no')
        call yaml_comment('iproc=0')
       ! write(*,'(1x,a)') 'CUDA support activated (iproc=0)'
    end if

  else if (matacc%iacceleration >= 2) then
     ! OpenCL convolutions are activated
     ! use CUBLAS for the linear algebra for the moment
     if (.not. OCLconv) then
        call MPI_COMM_SIZE(bigdft_mpi%mpi_comm,mproc,ierr)
        !initialize the id_proc per node
        call processor_id_per_node(iproc,mproc,GPU%id_proc,nproc_node)
        !initialize the opencl context for any process in the node
        !call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
        !do jproc=0,mproc-1
        !   call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
        !   if (iproc == jproc) then
        !      print '(a,a,i4,i4)','Initializing for node: ',trim(nodename_local),iproc,GPU%id_proc
        call init_acceleration_OCL(matacc,GPU)
        !   end if
        !end do
        GPU%ndevices=min(GPU%ndevices,nproc_node)
        if (iproc == 0) then
           call yaml_map('Material acceleration','OpenCL',advance='no')
           call yaml_comment('iproc=0')
           call yaml_open_map('Number of OpenCL devices per node',flow=.true.)
           call yaml_map('used',trim(yaml_toa(min(GPU%ndevices,nproc_node),fmt='(i0)')))
           call yaml_map('available',trim(yaml_toa(GPU%ndevices,fmt='(i0)')))
           !write(*,'(1x,a,i5,i5)') 'OpenCL support activated, No. devices per node (used, available):',&
           !     min(GPU%ndevices,nproc_node),GPU%ndevices
           call yaml_close_map()
        end if
        !the number of devices is the min between the number of processes per node
        GPU%ndevices=min(GPU%ndevices,nproc_node)
        OCLconv=.true.
     end if

  else
     if (iproc == 0) then
        call yaml_map('Material acceleration',.false.,advance='no')
        call yaml_comment('iproc=0')
        ! write(*,'(1x,a)') 'No material acceleration (iproc=0)'
     end if
  end if

END SUBROUTINE init_material_acceleration


subroutine release_material_acceleration(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  
  if (GPUconv) then
     call sg_end()
  end if

  if (OCLconv) then
     call release_acceleration_OCL(GPU)
     OCLconv=.false.
  end if

END SUBROUTINE release_material_acceleration


!> Give the number of MPI processes per node (nproc_node) and before iproc (iproc_node)
subroutine processor_id_per_node(iproc,nproc,iproc_node,nproc_node)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(out) :: iproc_node,nproc_node
  !local variables
  character(len=*), parameter :: subname='processor_id_per_node'
  integer :: ierr,namelen,i_stat,i_all,jproc
  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename

  if (nproc == 1) then
     iproc_node=0
     nproc_node=1
  else
     allocate(nodename(0:nproc-1+ndebug),stat=i_stat)
     call memocc(i_stat,nodename,'nodename',subname)
     
     !initalise nodenames
     do jproc=0,nproc-1
        nodename(jproc)=repeat(' ',MPI_MAX_PROCESSOR_NAME)
     end do

     call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)

     !gather the result between all the process
     call MPI_ALLGATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
          nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
          bigdft_mpi%mpi_comm,ierr)

     !found the processors which belong to the same node
     !before the processor iproc
     iproc_node=0
     do jproc=0,iproc-1
        if (trim(nodename(jproc)) == trim(nodename(iproc))) then
           iproc_node=iproc_node+1
        end if
     end do
     nproc_node=iproc_node
     do jproc=iproc,nproc-1
        if (trim(nodename(jproc)) == trim(nodename(iproc))) then
           nproc_node=nproc_node+1
        end if
     end do
     
     i_all=-product(shape(nodename))*kind(nodename)
     deallocate(nodename,stat=i_stat)
     call memocc(i_stat,i_all,'nodename',subname)
  end if
END SUBROUTINE processor_id_per_node


!> this routine does the same operations as
!! read_atomic_file but uses inputs from memory
!! as input positions instead of inputs from file
!! Useful for QM/MM implementation of BigDFT-ART
!! @author Written by Laurent K Beland 2011 UdeM
subroutine initialize_atomic_file(iproc,atoms,rxyz)
  use module_base
  use module_types
  use module_interfaces, except_this_one => initialize_atomic_file
  use m_ab6_symmetry
  use yaml_output
  implicit none
  integer, intent(in) :: iproc
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(:,:), pointer :: rxyz
  !local variables
  character(len=*), parameter :: subname='initialize_atomic_file'
  integer :: i_stat
  integer :: iat,i,ierr

  allocate(atoms%amu(atoms%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%amu,'atoms%amu',subname)

  if (atoms%astruct%geocode=='S') then 
        atoms%astruct%cell_dim(2)=0.0_gp
  else if (atoms%astruct%geocode=='F') then !otherwise free bc    
        atoms%astruct%cell_dim(1)=0.0_gp
        atoms%astruct%cell_dim(2)=0.0_gp
        atoms%astruct%cell_dim(3)=0.0_gp
  else
        atoms%astruct%cell_dim(1)=0.0_gp
        atoms%astruct%cell_dim(2)=0.0_gp
        atoms%astruct%cell_dim(3)=0.0_gp
  end if

  !reduced coordinates are possible only with periodic units
  if (atoms%astruct%units == 'reduced' .and. atoms%astruct%geocode == 'F') then
     if (iproc==0) write(*,'(1x,a)')&
          'ERROR: Reduced coordinates are not allowed with isolated BC'
  end if

   !convert the values of the cell sizes in bohr
  if (atoms%astruct%units=='angstroem' .or. atoms%astruct%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     atoms%astruct%cell_dim(1)=atoms%astruct%cell_dim(1)/Bohr_Ang
     atoms%astruct%cell_dim(2)=atoms%astruct%cell_dim(2)/Bohr_Ang
     atoms%astruct%cell_dim(3)=atoms%astruct%cell_dim(3)/Bohr_Ang
  else if (atoms%astruct%units == 'reduced') then
     !assume that for reduced coordinates cell size is in bohr
     atoms%astruct%cell_dim(1)=real(atoms%astruct%cell_dim(1),gp)
     atoms%astruct%cell_dim(2)=real(atoms%astruct%cell_dim(2),gp)
     atoms%astruct%cell_dim(3)=real(atoms%astruct%cell_dim(3),gp)
  else
     call yaml_warning('Length units in input file unrecognized')
     call yaml_warning('recognized units are angstroem or atomic = bohr')
     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
  endif
  
  do iat=1,atoms%astruct%nat
     !xyz input file, allow extra information
     
     if (atoms%astruct%units == 'reduced') then !add treatment for reduced coordinates
        rxyz(1,iat)=modulo(rxyz(1,iat),1.0_gp)
        if (atoms%astruct%geocode == 'P') rxyz(2,iat)=modulo(rxyz(2,iat),1.0_gp)
        rxyz(3,iat)=modulo(rxyz(3,iat),1.0_gp)
     else if (atoms%astruct%geocode == 'P') then
        rxyz(1,iat)=modulo(rxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(2,iat)=modulo(rxyz(2,iat),atoms%astruct%cell_dim(2))
        rxyz(3,iat)=modulo(rxyz(3,iat),atoms%astruct%cell_dim(3))
     else if (atoms%astruct%geocode == 'S') then
        rxyz(1,iat)=modulo(rxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(3,iat)=modulo(rxyz(3,iat),atoms%astruct%cell_dim(3))
     end if
 
     if (atoms%astruct%units=='angstroem' .or. atoms%astruct%units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/Bohr_Ang
        enddo
     else if (atoms%astruct%units == 'reduced') then 
        rxyz(1,iat)=rxyz(1,iat)*atoms%astruct%cell_dim(1)
        if (atoms%astruct%geocode == 'P') rxyz(2,iat)=rxyz(2,iat)*atoms%astruct%cell_dim(2)
        rxyz(3,iat)=rxyz(3,iat)*atoms%astruct%cell_dim(3)
     endif
  enddo

  !control atom positions
  call check_atoms_positions(atoms,rxyz,(iproc == 0))

  ! We delay the calculation of the symmetries.
  atoms%astruct%sym%symObj = -1
  nullify(atoms%astruct%sym%irrzon)
  nullify(atoms%astruct%sym%phnons)

END SUBROUTINE initialize_atomic_file

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

  !grid spacings
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
!!$  call input_dict_var(in%nvirt, dict_dft//NVIRT, 0, ranges=(/0,abs(in%norbv)/))
!!$  call input_dict_var(in%nplot, dict_dft//NPLOT, 0, ranges=(/0,abs(in%norbv)/))

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
  character(len = 1), intent(in) :: geocode
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

!!$  call yaml_map('Error found',f_err_check())
!!$  if (f_err_check()) then
!!$     call f_dump_all_errors()
!!$  end if
!!$  stop

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

  ! Performence and setting-up options
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

