!> @file
!!  Define main module for using BigDFT as a blackbox
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG, DC)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module bigdft_run
  use module_defs, only: gp
  use dictionaries
  use module_types, only: input_variables,DFT_wavefunction,GPU_pointers,energy_terms
  use module_atoms, only: atoms_data

  private

  !>  Used to restart a new DFT calculation or to save information 
  !!  for post-treatment
  type, public :: restart_objects
     integer :: version !< 0=cubic, 100=linear
     integer :: n1,n2,n3,nat
     real(gp) :: hx_old,hy_old,hz_old
     real(gp), dimension(:,:), pointer :: rxyz_old,rxyz_new
     type(DFT_wavefunction) :: KSwfn !< Kohn-Sham wavefunctions
     type(DFT_wavefunction) :: tmb !<support functions for linear scaling
     type(GPU_pointers) :: GPU 
  end type restart_objects


  !> Public container to be used with call_bigdft().
  type, public :: run_objects
     type(dictionary), pointer :: user_inputs

     type(input_variables), pointer    :: inputs
     type(atoms_data), pointer         :: atoms
     type(restart_objects), pointer    :: rst
     real(gp), dimension(:,:), pointer :: radii_cf
  end type run_objects


  !> Used to store results of a DFT calculation.
  type, public :: DFT_global_output
     real(gp) :: energy, fnoise, pressure      !< Total energy, noise over forces and pressure
     type(energy_terms) :: energs              !< All energy terms
     integer :: fdim                           !< Dimension of allocated forces (second dimension)
     real(gp), dimension(:,:), pointer :: fxyz !< Atomic forces
     real(gp), dimension(6) :: strten          !< Stress Tensor
  end type DFT_global_output

  public :: init_global_output,deallocate_global_output,restart_objects_set_mat_acc
  public :: run_objects_free,copy_global_output,restart_objects_set_mode
  public :: run_objects_nullify,restart_objects_set_nat,restart_objects_new
  public :: run_objects_associate,run_objects_free_container,init_restart_objects
  public :: global_output_set_from_dict,free_restart_objects
  public :: run_objects_init,bigdft_init,bigdft_command_line_options,bigdft_nruns

!!$  ! interfaces of external routines 
!!$  interface
!!$     subroutine geopt(runObj,outs,nproc,iproc,ncount_bigdft)
!!$       use module_base
!!$       use module_types
!!$       implicit none
!!$       type(run_objects), intent(inout) :: runObj
!!$       type(DFT_global_output), intent(inout) :: outs
!!$       integer, intent(in) :: nproc,iproc
!!$       integer, intent(inout) :: ncount_bigdft
!!$     END SUBROUTINE geopt
!!$
!!$  end interface
  

  contains
    
    !> All in one routine to initialise and set-up restart objects.
    subroutine init_restart_objects(iproc,inputs,atoms,rst)
      implicit none
      !Arguments
      integer, intent(in) :: iproc
      type(input_variables), intent(in) :: inputs
      type(atoms_data), intent(in) :: atoms
      type(restart_objects), intent(out) :: rst

      call restart_objects_new(rst)
      call restart_objects_set_mode(rst, inputs%inputpsiid)
      call restart_objects_set_nat(rst, atoms%astruct%nat)
      call restart_objects_set_mat_acc(rst, iproc, inputs%matacc)
    END SUBROUTINE init_restart_objects

    !> Allocate and nullify restart objects
    pure subroutine restart_objects_new(rst)
      use module_defs, only: UNINITIALIZED
      use module_types, only: nullify_local_zone_descriptors,CUBIC_VERSION
      use locregs, only: nullify_locreg_descriptors
      implicit none
      !Arguments
      type(restart_objects), intent(out) :: rst

      ! Decide whether we use the cubic or the linear version
      rst%version = UNINITIALIZED(CUBIC_VERSION)

      !allocate pointers
      rst%nat = 0
      nullify(rst%rxyz_new)
      nullify(rst%rxyz_old)

      !nullify unallocated pointers
      rst%KSwfn%c_obj = 0
      nullify(rst%KSwfn%psi)
      nullify(rst%KSwfn%orbs%eval)

      nullify(rst%KSwfn%gaucoeffs)
      nullify(rst%KSwfn%oldpsis)

      call nullify_locreg_descriptors(rst%KSwfn%Lzd%Glr)
      nullify(rst%KSwfn%Lzd%Glr%wfd%keyglob)
      nullify(rst%KSwfn%Lzd%Glr%wfd%keygloc)
      nullify(rst%KSwfn%Lzd%Glr%wfd%keyvloc)
      nullify(rst%KSwfn%Lzd%Glr%wfd%keyvglob)

      nullify(rst%KSwfn%gbd%nshell)
      nullify(rst%KSwfn%gbd%ndoc)
      nullify(rst%KSwfn%gbd%nam)
      nullify(rst%KSwfn%gbd%xp)
      nullify(rst%KSwfn%gbd%psiat)
      nullify(rst%KSwfn%gbd%rxyz)

      !Nullify LZD for cubic version (new input guess)
      call nullify_local_zone_descriptors(rst%tmb%lzd)

      !Nullify GPU data
      rst%GPU%OCLconv=.false.
    END SUBROUTINE restart_objects_new


    pure subroutine restart_objects_set_mode(rst, inputpsiid)
      use module_types
      implicit none
      type(restart_objects), intent(inout) :: rst
      integer, intent(in) :: inputpsiid

      select case (inputpsiid)
      case (INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_CP2K, INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, &
           INPUT_PSI_DISK_WVL, INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS)
         rst%version = CUBIC_VERSION
      case (INPUT_PSI_LINEAR_AO, INPUT_PSI_MEMORY_LINEAR, INPUT_PSI_DISK_LINEAR)
         rst%version = LINEAR_VERSION
      end select
    END SUBROUTINE restart_objects_set_mode


    subroutine restart_objects_set_nat(rst, nat)
      use module_base
      implicit none
      !Arguments
      integer, intent(in) :: nat
      type(restart_objects), intent(inout) :: rst

      call f_free_ptr(rst%rxyz_old)
      call f_free_ptr(rst%rxyz_new)

      rst%nat = nat
      rst%rxyz_new = f_malloc_ptr((/ 3, nat /),id='rst%rxyz_new')
      rst%rxyz_old = f_malloc_ptr((/ 3, nat /),id='rst%rxyz_old')
    END SUBROUTINE restart_objects_set_nat

    subroutine restart_objects_set_mat_acc(rst, iproc, matacc)
      use module_types, only: material_acceleration
      implicit none
      !Arguments
      type(restart_objects), intent(inout) :: rst
      integer, intent(in) :: iproc
      type(material_acceleration), intent(in) :: matacc
      !initialise the acceleration strategy if required
      call init_material_acceleration(iproc,matacc,rst%GPU)
    END SUBROUTINE restart_objects_set_mat_acc


    !> De-Allocate restart_objects
    subroutine free_restart_objects(rst)
      use module_base
      use locregs
      use gaussians, only: deallocate_gwf
      use module_types
      implicit none
      type(restart_objects) :: rst
      !local variables
      integer :: istep

      if (rst%version == LINEAR_VERSION) then
         call destroy_DFT_wavefunction(rst%tmb)
      end if
      !always deallocate lzd for new input guess
      !call deallocate_lzd(rst%tmb%lzd)
      ! Modified by SM
      call deallocate_local_zone_descriptors(rst%tmb%lzd)

      call deallocate_locreg_descriptors(rst%KSwfn%Lzd%Glr)

      call f_free_ptr(rst%KSwfn%psi)
      call f_free_ptr(rst%KSwfn%orbs%eval)
      if (associated(rst%KSwfn%oldpsis)) then
         do istep=0,product(shape(rst%KSwfn%oldpsis))-1
            call old_wavefunction_free(rst%KSwfn%oldpsis(istep))
         end do
         deallocate(rst%KSwfn%oldpsis)
      end if
      call f_free_ptr(rst%rxyz_old)
      call f_free_ptr(rst%rxyz_new)

      !The gaussian basis descriptors are always allocated together
      !with the gaussian coefficients
      if (associated(rst%KSwfn%gbd%rxyz)) then
         nullify(rst%KSwfn%gbd%rxyz)
         call deallocate_gwf(rst%KSwfn%gbd)
      end if
      call f_free_ptr(rst%KSwfn%gaucoeffs)

      !finalise the material accelearion usage
      call release_material_acceleration(rst%GPU)

    END SUBROUTINE free_restart_objects


    !> Initialize the structure DFT_global_output
    subroutine nullify_global_output(outs)
      use module_defs, only: UNINITIALIZED
      use module_types, only: energy_terms_null
      implicit none
      type(DFT_global_output), intent(out) :: outs

      outs%energs=energy_terms_null()
      outs%fdim      = 0
      nullify(outs%fxyz)
      outs%energy    = UNINITIALIZED(1.0_gp)
      outs%fnoise    = UNINITIALIZED(1.0_gp)
      outs%pressure  = UNINITIALIZED(1.0_gp)
      outs%strten(:) = UNINITIALIZED(1.0_gp)
    END SUBROUTINE nullify_global_output


    subroutine init_global_output(outs, nat)
      use module_base
      use dynamic_memory
      implicit none
      type(DFT_global_output), intent(out) :: outs
      integer, intent(in) :: nat

      call nullify_global_output(outs)
      outs%fdim = nat
      outs%fxyz = f_malloc_ptr((/ 3, outs%fdim /),id='outs%fxyz')
      outs%fxyz(:,:) = UNINITIALIZED(1.0_gp)
    END SUBROUTINE init_global_output

    subroutine deallocate_global_output(outs, fxyz)
      use module_base
      implicit none
      type(DFT_global_output), intent(inout) :: outs
      real(gp), intent(out), optional :: fxyz

      if (associated(outs%fxyz)) then
         if (present(fxyz)) &
              call f_memcpy(src=outs%fxyz(1,1),dest=fxyz,n=3*outs%fdim)
         !call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz, 1)
         !end if
         call f_free_ptr(outs%fxyz)
      end if
    END SUBROUTINE deallocate_global_output

    !> Copies outsA to outsB
    subroutine copy_global_output(outsA,outsB)
      use module_base, only: f_err_throw,f_memcpy
      use yaml_strings, only: yaml_toa
      implicit none
      type(DFT_global_output), intent(in) :: outsA
      type(DFT_global_output), intent(inout) :: outsB
      integer :: i

      if(outsA%fdim /= outsB%fdim)then
         call f_err_throw("Error in copy_global_output: outsA and outsB have different sizes"//&
              trim(yaml_toa(outsA%fdim))//trim(yaml_toa(outsB%fdim)),&
              err_name='BIGDFT_RUNTIME_ERROR')
      endif
      !outsA%fdim == outsB%fdim so it does not have to be copied

      outsB%energy = outsA%energy
      outsB%fnoise = outsA%fnoise
      outsB%pressure = outsA%pressure
      !8.5.2014: outs%energs does not contain any pointers,
      !so we use intrinisc copy:
      outsB%energs = outsA%energs
      call f_memcpy(src=outsA%fxyz,dest=outsB%fxyz)
      !call vcopy(3 * outsB%fdim, outsA%fxyz(1,1), 1, outsB%fxyz(1,1), 1)
      call f_memcpy(src=outsA%strten,dest=outsB%strten)
!!$      do i=1,6
!!$         outsB%strten(i) = outsA%strten(i)
!!$      enddo
    end subroutine copy_global_output

    !> Associate to the structure run_objects, the input_variable structure and the atomic positions (atoms_data)
    subroutine run_objects_associate(runObj, inputs, atoms, rst, rxyz0)
      use module_base
      use module_types
      implicit none
      type(run_objects), intent(out) :: runObj
      type(input_variables), intent(in), target :: inputs
      type(atoms_data), intent(in), target :: atoms
      type(restart_objects), intent(in), target :: rst
      real(gp), intent(in), optional :: rxyz0

      call run_objects_free_container(runObj)
      runObj%atoms  => atoms
      runObj%inputs => inputs
      runObj%rst    => rst
      if (present(rxyz0)) then
         call vcopy(3 * atoms%astruct%nat, rxyz0, 1, runObj%atoms%astruct%rxyz(1,1), 1)
      end if

      runObj%radii_cf = f_malloc_ptr((/ runObj%atoms%astruct%ntypes, 3 /), id="runObj%radii_cf")
      call read_radii_variables(runObj%atoms, runObj%radii_cf, &
           & runObj%inputs%crmult, runObj%inputs%frmult, runObj%inputs%projrad)
    END SUBROUTINE run_objects_associate

!!    subroutine global_output_merge_to_dict(dict, outs, astruct)
!!      use module_defs, only: gp, UNINITIALIZED
!!      use module_atoms, only: atomic_structure
!!      use dictionaries
!!      implicit none
!!      type(dictionary), pointer :: dict
!!      type(DFT_global_output), intent(in) :: outs
!!      type(atomic_structure), intent(in) :: astruct
!!
!!      integer :: iat
!!      type(dictionary), pointer :: pos, fxyz
!!
!!      if (has_key(dict, GOUT_FORCES)) call dict_remove(dict, GOUT_FORCES)
!!      if (associated(outs%fxyz)) then
!!         pos => dict // GOUT_FORCES
!!         do iat=1,astruct%nat
!!!!$            call dict_init(fxyz)
!!!!$            call set(fxyz // astruct%atomnames(astruct%iatype(iat)) // 0, outs%fxyz(1, iat))
!!!!$            call set(fxyz // astruct%atomnames(astruct%iatype(iat)) // 1, outs%fxyz(2, iat))
!!!!$            call set(fxyz // astruct%atomnames(astruct%iatype(iat)) // 2, outs%fxyz(3, iat))
!!!!$            call add(pos, fxyz)
!!            call add(pos, dict_new(astruct%atomnames(astruct%iatype(iat)) .is. outs%fxyz(:,iat)))
!!         end do
!!      end if
!!
!!      if (has_key(dict, GOUT_ENERGY)) call dict_remove(dict, GOUT_ENERGY)
!!      if (outs%energy /= UNINITIALIZED(outs%energy)) &
!!           & call set(dict // GOUT_ENERGY, outs%energy)
!!
!!    end subroutine global_output_merge_to_dict

    subroutine global_output_set_from_dict(outs, dict)
      use dictionaries
      use module_input_dicts, only: GOUT_ENERGY,GOUT_FORCES
      implicit none
      type(dictionary), pointer :: dict
      type(DFT_global_output), intent(inout) :: outs

      integer :: i
      type(dictionary), pointer :: it,it0

      i=0
      it0 => dict_iter(dict .get. GOUT_FORCES)
      do while(associated(it0))
         !this will be done only once if the key exists
         if (.not. associated(outs%fxyz)) &
              call init_global_output(outs, dict_len(dict // GOUT_FORCES))

         i=i+1
         it => dict_iter(it0)
         find_forces: do while(associated(it))
            if (dict_len(it) == 3) then
               outs%fxyz(1:3, i) = it 
               exit find_forces
            end if
            it => dict_next(it)
         end do find_forces

         it0 => dict_next(it0)
      end do
      outs%energy = dict .get. GOUT_ENERGY

!!$      if (has_key(dict, GOUT_FORCES)) then
!!$         if (.not. associated(outs%fxyz)) &
!!$              & call init_global_output(outs, dict_len(dict // GOUT_FORCES))
!!$         do i = 1, outs%fdim, 1
!!$            it => dict_iter(dict // GOUT_FORCES // (i - 1))
!!$            find_forces: do while (associated(it))
!!$               if (dict_len(it) == 3) then
!!$                  outs%fxyz(1:3, i) = it 
!!$                  exit find_forces
!!$               end if
!!$               it => dict_next(it)
!!$            end do find_forces
!!$         end do
!!$      end if
!!$
!!$      if (has_key(dict, GOUT_ENERGY)) outs%energy = dict // GOUT_ENERGY
    end subroutine global_output_set_from_dict

    !> Routines to handle the argument objects of call_bigdft().
    subroutine run_objects_nullify(runObj)
      use module_types
      implicit none
      type(run_objects), intent(out) :: runObj

      nullify(runObj%user_inputs)
      nullify(runObj%inputs)
      nullify(runObj%atoms)
      nullify(runObj%rst)
      nullify(runObj%radii_cf)
    END SUBROUTINE run_objects_nullify


    !> Freed the run_objects structure
    subroutine run_objects_free(runObj)
      use module_types
      use module_base
      use dynamic_memory
      use yaml_output
      use dictionaries
      use  module_atoms, only: deallocate_atoms_data
      implicit none
      type(run_objects), intent(inout) :: runObj

      call dict_free(runObj%user_inputs)
      if (associated(runObj%rst)) then
         call free_restart_objects(runObj%rst)
         deallocate(runObj%rst)
      end if
      if (associated(runObj%atoms)) then
         call deallocate_atoms_data(runObj%atoms) 
         deallocate(runObj%atoms)
      end if
      if (associated(runObj%inputs)) then
         call free_input_variables(runObj%inputs)
         deallocate(runObj%inputs)
      end if
      call f_free_ptr(runObj%radii_cf)

    END SUBROUTINE run_objects_free

    !> Deallocate run_objects
    subroutine run_objects_free_container(runObj)
      use module_types
      use module_base
      use dynamic_memory
      use yaml_output
      implicit none
      type(run_objects), intent(inout) :: runObj

      ! User inputs are always owned by run objects.
      call dict_free(runObj%user_inputs)
      ! Radii_cf are always owned by run objects.
      call f_free_ptr(runObj%radii_cf)
      ! Currently do nothing except nullifying everything.
      call run_objects_nullify(runObj)
    END SUBROUTINE run_objects_free_container

    !> Read all input files and create the objects to run BigDFT
    subroutine run_objects_init(runObj, run_dict)
      use module_base, only: bigdft_mpi,dict_init
      use module_types
      use module_input_dicts, only: user_dict_from_files
      implicit none
      type(run_objects), intent(out) :: runObj
      type(dictionary), pointer :: run_dict
      !local variables
      character(len=max_field_length) :: radical, posinp

      radical = run_dict // 'name'
      posinp = run_dict // 'posinp' 

      call run_objects_nullify(runObj)

      ! Allocate persistent structures.
      allocate(runObj%rst)
      call restart_objects_new(runObj%rst)

      ! Generate input dictionary and parse it.
      call dict_init(runObj%user_inputs)
      call user_dict_from_files(runObj%user_inputs, radical, posinp, bigdft_mpi)
      call run_objects_parse(runObj)

      ! Start the signaling loop in a thread if necessary.
      if (runObj%inputs%signaling .and. bigdft_mpi%iproc == 0) then
         call bigdft_signals_init(runObj%inputs%gmainloop, 2, &
              & runObj%inputs%domain, len_trim(runObj%inputs%domain))
         call bigdft_signals_start(runObj%inputs%gmainloop, runObj%inputs%signalTimeout)
      end if
    END SUBROUTINE run_objects_init

    subroutine bigdft_init(options)
      use yaml_parse
      use dictionaries
      !use yaml_output, only: yaml_map
      use yaml_strings, only: f_strcpy,yaml_toa
      use module_defs, only: bigdft_mpi
      use module_input_dicts, only: merge_input_file_to_dict
      implicit none
      !> dictionary of the options of the run
      !! on entry, it contains the options for initializing
      !! on exit, it contains in the key "BigDFT", a list of the 
      !! dictionaries of each of the run that the local instance of BigDFT
      !! code has to execute.
      !! if this argument is not present, the code is only initialized
      !! in its normal mode: no taskgroups and default values of radical and posinp
      type(dictionary), pointer, optional :: options
      !local variables
      logical :: exist_list,posinp_name
      integer :: ierr,mpi_groupsize,iconfig
      character(len=max_field_length) :: posinp_id,run_id,err_msg
      integer, dimension(4) :: mpi_info
      type(dictionary), pointer :: dict_run,opts

      !coherence checks among the options (no disk access)

      !Initalize the global mpi environment
      call bigdft_mpi_init(ierr)
      !if (ierr /= MPI_SUCCESS) then
      !this part have to be included in mpi_init wrappers
      !   return
      !end if
      nullify(opts)
      if (present(options)) opts => options

      !taskgroup size
      mpi_groupsize=0
      mpi_groupsize=opts .get. 'taskgroup-size'

      !initialize the bigdft_mpi environment
      call bigdft_init_mpi_env(mpi_info, mpi_groupsize, ierr)
      !the error check has to be ierr

      !identify the list of the runs which are associated to the 
      !present processor
      nullify(dict_run)

      !logical flag telling that the input position name is still "posinp.*"
      !regardless of the fact that the input file is "input.*"
      posinp_name=.false.
      if ('name' .in. opts) then
         !check if the names are given as a list or as a scalar
         if (dict_len(opts) > 0) then
            call dict_copy(dict_run,opts//'name')
         else
            run_id = opts//'name'
         end if
      else
         call f_strcpy(src='input',dest=run_id)
         posinp_name=.true.
      end if
      !this is not possible if name option exists
      if ('runs-file' .in. opts) then
         posinp_name=.false.
         posinp_id = opts//'runs-file'
         !verify if file exists and if it is yaml compliant
         call f_file_exists(posinp_id,exist_list)
         if (exist_list) then
            call dict_init(dict_run)
            call f_err_open_try()
            call merge_input_file_to_dict(dict_run,posinp_id,bigdft_mpi)
            !verify if yaml_parsing found errors
            if (f_err_check()) ierr = f_get_last_error(err_msg)
            call f_err_close_try()
            !sanity check of the parsed file
            if (ierr /= 0) then
               call f_err_throw('Parsing error for runs-file "'//&
                    trim(posinp_id)//&
                    '", with message '//trim(err_msg),&
                    err_name='BIGDFT_INPUT_VARIABLES_ERROR')
            end if
            !check if the file is given by a list
            if (dict_len(dict_run) <= 0) then
               call f_err_throw('The runs-file "'//&
                    trim(posinp_id)//'" is not a list in yaml syntax',& 
                    err_name='BIGDFT_INPUT_VARIABLES_ERROR')
            end if
         else
            call f_err_throw('The runs-file specified ('//trim(posinp_id)//&
                 ') does not exists',err_name='BIGDFT_INPUT_FILE_ERROR')
         end if
         !here the run of the dicts has to be evaluated according to the taskgroups    
      else if (.not. associated(dict_run)) then
         call dict_init(dict_run)
         if (bigdft_mpi%ngroup == 1) then
            call add(dict_run,trim(run_id))
         else
            do iconfig=1,bigdft_mpi%ngroup
               call add(dict_run,trim(run_id)//&
                    trim(adjustl(yaml_toa(iconfig,fmt='(i3)'))))
            end do
         end if
      end if

      !call yaml_map('Dict of runs',dict_run)

      if (present(options)) then
         if (.not. associated(options)) call dict_init(options)
         !here the dict_run is given, and in each of the taskgroups a list of 
         !runs for BigDFT instances has to be given
         do iconfig=0,dict_len(dict_run)-1
            if (modulo(iconfig,bigdft_mpi%ngroup)==bigdft_mpi%igroup) then
               run_id=dict_run//iconfig
               if (posinp_name) then
                  if (dict_len(dict_run) == 1) then
                     call f_strcpy(src='posinp',dest=posinp_id)
                  else
                     call f_strcpy(src='posinp'//&
                          trim(adjustl(yaml_toa(iconfig,fmt='(i3)'))),&
                          dest=posinp_id)
                  end if
               else
                  posinp_id=run_id
               end if
               call add(options//'BigDFT',&
                    dict_new('name' .is. run_id,'posinp' .is. posinp_id))
            end if
         end do
      end if
      call dict_free(dict_run)

    end subroutine bigdft_init

    !>identify the options from command line
    !! and write the result in options dict
    subroutine bigdft_command_line_options(options)
      use yaml_parse
      use dictionaries
      implicit none
      !> dictionary of the options of the run
      !! on entry, it contains the options for initializing
      !! on exit, it contains in the key "BigDFT", a list of the 
      !! dictionaries of each of the run that the local instance of BigDFT
      !! code has to execute
      type(dictionary), pointer :: options
      !local variables
      type(yaml_cl_parse) :: parser !< command line parser

      !define command-line options
      parser=yaml_cl_parse_null()
      !between these lines, for another executable using BigDFT as a blackbox,
      !other command line options can be specified
      !then the bigdft options can be specified
      call bigdft_options(parser)
      !parse command line, and retrieve arguments
      call yaml_cl_parse_cmd_line(parser,args=options)
      !free command line parser information
      call yaml_cl_parse_free(parser)

    end subroutine bigdft_command_line_options

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

      call yaml_cl_parse_option(parser,'outdir','.',&
           'output directory','d',&
           dict_new('Usage' .is. &
           'Set the directory where all the output files have to be written.',&
           'Allowed values' .is. &
           'String value, indicating the path of the directory. If the last subdirectory is not existing, it will be created'))

      call yaml_cl_parse_option(parser,'logfile','No',&
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

    !> retrieve the number of runs for a given set of options
    !! gives 0 if the option dictionary is invalid 
    function bigdft_nruns(options)
      implicit none
      type(dictionary), pointer :: options !< filled options dictionary. bigdft_init has to be called
      integer :: bigdft_nruns
      type(dictionary), pointer :: runs
      runs = options .get. 'BigDFT'
      bigdft_nruns=dict_len(runs)
      if (bigdft_nruns < 0) bigdft_nruns=0
    end function bigdft_nruns

  
end module bigdft_run

!external wrapper temporary to make the code compiling with wrappers
subroutine run_objects_init_from_run_name(runObj, radical, posinp)
  use module_base
  use bigdft_run
  implicit none
  type(run_objects), intent(out) :: runObj
  character(len = *), intent(in) :: radical, posinp
  !local variables
  type(dictionary), pointer :: run_dict

  !create the ad-hoc dictionary run to wrap the module routine
  run_dict => dict_new('name' .is. radical, 'posinp' .is. posinp)
  
  call run_objects_init(runObj,run_dict)

  call dict_free(run_dict)
END SUBROUTINE run_objects_init_from_run_name


!> Routine to use BigDFT as a blackbox
subroutine call_bigdft(runObj,outs,nproc,iproc,infocode)
  use module_base
  use module_interfaces, only: write_atomic_file
  use yaml_output
  use module_types, only: INPUT_PSI_MEMORY_LINEAR,LINEAR_VERSION,INPUT_PSI_MEMORY_GAUSS, &
       INPUT_PSI_LCAO,INPUT_PSI_MEMORY_WVL,old_wavefunction_null,INPUT_PSI_LINEAR_AO,deallocate_wfd
  use bigdft_run
  !use communications_base
  implicit none
  integer, intent(in) :: iproc,nproc
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  integer, intent(inout) :: infocode
  !local variables
  character(len=*), parameter :: subname='call_bigdft'
  character(len=40) :: comment
  logical :: exists
  integer :: ierr,inputPsiId_orig,iat,iorb,istep
  real(gp) :: maxdiff
  external :: cluster,forces_via_finite_differences

  !put a barrier for all the processes
  call mpibarrier(bigdft_mpi%mpi_comm)
  call f_routine(id=subname)
  !Check the consistency between MPI processes of the atomic coordinates
  maxdiff=mpimaxdiff(runObj%atoms%astruct%rxyz,comm=bigdft_mpi%mpi_comm,bcast=.true.)
  !even when broadcasting results, there can be truncation errors
  !therefore increase the threshold for checking
  if (maxdiff > epsilon(1.0_gp)) then
     if (iproc==0) then
        call yaml_warning('Input positions not identical! '//&
             '(difference:'//trim(yaml_toa(maxdiff))//' ), broadcasting from master node.')
        call yaml_comment('If the code hangs here, this means that not all the tasks met the threshold')
        call yaml_flush_document()
     end if
        !the check=.true. is important here: it controls that each process
     !will participate in the broadcasting
     call mpibcast(runObj%atoms%astruct%rxyz,comm=bigdft_mpi%mpi_comm,&
          check=.true.)
  end if

  !fill the rxyz array with the positions
  !wrap the atoms in the periodic directions when needed
  select case(runObj%atoms%astruct%geocode)
  case('P')
     do iat=1,runObj%atoms%astruct%nat
        runObj%rst%rxyz_new(1,iat)=modulo(runObj%atoms%astruct%rxyz(1,iat),runObj%atoms%astruct%cell_dim(1))
        runObj%rst%rxyz_new(2,iat)=modulo(runObj%atoms%astruct%rxyz(2,iat),runObj%atoms%astruct%cell_dim(2))
        runObj%rst%rxyz_new(3,iat)=modulo(runObj%atoms%astruct%rxyz(3,iat),runObj%atoms%astruct%cell_dim(3))
     end do
  case('S')
     do iat=1,runObj%atoms%astruct%nat
        runObj%rst%rxyz_new(1,iat)=modulo(runObj%atoms%astruct%rxyz(1,iat),runObj%atoms%astruct%cell_dim(1))
        runObj%rst%rxyz_new(2,iat)=runObj%atoms%astruct%rxyz(2,iat)
        runObj%rst%rxyz_new(3,iat)=modulo(runObj%atoms%astruct%rxyz(3,iat),runObj%atoms%astruct%cell_dim(3))
     end do
  case('W')
     do iat=1,runObj%atoms%astruct%nat
        runObj%rst%rxyz_new(1,iat)=runObj%atoms%astruct%rxyz(1,iat)
        runObj%rst%rxyz_new(2,iat)=runObj%atoms%astruct%rxyz(2,iat)
        runObj%rst%rxyz_new(3,iat)=modulo(runObj%atoms%astruct%rxyz(3,iat),runObj%atoms%astruct%cell_dim(3))
     end do
  case('F')
     do iat=1,runObj%atoms%astruct%nat
        runObj%rst%rxyz_new(1,iat)=runObj%atoms%astruct%rxyz(1,iat)
        runObj%rst%rxyz_new(2,iat)=runObj%atoms%astruct%rxyz(2,iat)
        runObj%rst%rxyz_new(3,iat)=runObj%atoms%astruct%rxyz(3,iat)
     end do
  end select

  !assign the verbosity of the output
  !the verbose variables is defined in module_base
  verbose=runObj%inputs%verbosity

  ! Use the restart for the linear scaling version... probably to be modified.
  if(runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_WVL) then
     if (runObj%rst%version == LINEAR_VERSION) then
        runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_LINEAR
        if (any(runObj%inputs%lin%locrad_lowaccuracy /= &
             runObj%inputs%lin%locrad_highaccuracy)) then
           call f_err_throw('The radii for low and high accuracy must be the same '//&
                'when using the linear restart!',&
                err_name='BIGDFT_INPUT_VARIABLES_ERROR')
        end if
     end if
  end if
  inputPsiId_orig=runObj%inputs%inputPsiId

  loop_cluster: do
     !allocate history container if it has not been done
     if (runObj%inputs%wfn_history > 1  .and. .not. associated(runObj%rst%KSwfn%oldpsis)) then
        allocate(runObj%rst%KSwfn%oldpsis(0:runObj%inputs%wfn_history+1))
        runObj%rst%KSwfn%istep_history=0
        do istep=0,runObj%inputs%wfn_history+1
           runObj%rst%KSwfn%oldpsis(istep)=old_wavefunction_null() 
        end do
     end if

     if (runObj%inputs%inputPsiId == INPUT_PSI_LCAO .and. associated(runObj%rst%KSwfn%psi)) then
        call f_free_ptr(runObj%rst%KSwfn%psi)
        call f_free_ptr(runObj%rst%KSwfn%orbs%eval)
        call deallocate_wfd(runObj%rst%KSwfn%Lzd%Glr%wfd)
     end if

     !experimental, finite difference method for calculating forces on particular quantities
     call f_file_exists('input.finite_difference_forces',exists)
     if (exists) then
        runObj%inputs%last_run=1 !do the last_run things nonetheless
        runObj%inputs%inputPsiId = INPUT_PSI_LCAO !the first run always restart from IG
        !experimental_modulebase_var_onlyfion=.true. !put only ionic forces in the forces
     end if

     !Main routine calculating the KS orbitals
     call cluster(nproc,iproc,runObj%atoms,runObj%rst%rxyz_new, runObj%radii_cf, &
          outs%energy, outs%energs, outs%fxyz, outs%strten, outs%fnoise, outs%pressure,&
          runObj%rst%KSwfn,runObj%rst%tmb,&
          runObj%rst%rxyz_old,runObj%rst%hx_old,runObj%rst%hy_old,runObj%rst%hz_old,&
          runObj%inputs,runObj%rst%GPU,infocode)

     !save the new atomic positions in the rxyz_old array
     call f_memcpy(src=runObj%rst%rxyz_new,dest=runObj%rst%rxyz_old)
     !do iat=1,runObj%atoms%astruct%nat
     !   runObj%rst%rxyz_old(1,iat)=runObj%rst%rxyz_new(1,iat)
     !   runObj%rst%rxyz_old(2,iat)=runObj%rst%rxyz_new(2,iat)
     !   runObj%rst%rxyz_old(3,iat)=runObj%rst%rxyz_new(3,iat)
     !enddo

     if (exists) then
        call forces_via_finite_differences(iproc,nproc,runObj%atoms,runObj%inputs, &
             & outs%energy,outs%fxyz,outs%fnoise,runObj%rst,infocode)
     end if

     !Check infocode in function of the inputPsiId parameters
     !and change the strategy of input guess psi
     if (runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_WVL .and. infocode==2) then
        if (runObj%inputs%gaussian_help) then
           runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_GAUSS
        else
           runObj%inputs%inputPsiId = INPUT_PSI_LCAO
        end if
     else if (runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_LINEAR .and. infocode==2) then
        ! problems after restart for linear version
        runObj%inputs%inputPsiId = INPUT_PSI_LINEAR_AO
     else if ((runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_WVL .or. &
          runObj%inputs%inputPsiId == INPUT_PSI_LCAO) .and. infocode==1) then
        !runObj%inputs%inputPsiId=INPUT_PSI_LCAO !better to diagonalise than to restart an input guess
        runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_WVL
        !if (iproc==0) then
        !   call yaml_warning('Self-consistent cycle did not meet convergence criteria')
        !   write(*,*)&
        !        &   ' WARNING: Self-consistent cycle did not meet convergence criteria'
        !end if
        exit loop_cluster
     else if (runObj%inputs%inputPsiId == INPUT_PSI_LCAO .and. infocode==3) then
        if (iproc == 0) then
!!$               write( *,'(1x,a)')'Convergence error, cannot proceed.'
!!$               write( *,'(1x,a)')' writing positions in file posfail.xyz then exiting'
!!$               write(comment,'(a)')'UNCONVERGED WF '
           !call wtxyz('posfail',energy,rxyz,atoms,trim(comment))
           call write_atomic_file("posfail",outs%energy,runObj%rst%rxyz_new,&
                runObj%atoms%astruct%ixyz_int,runObj%atoms,'UNCONVERGED WF ')
        end if

        call f_free_ptr(runObj%rst%KSwfn%psi)
        call f_free_ptr(runObj%rst%KSwfn%orbs%eval)

        call deallocate_wfd(runObj%rst%KSwfn%Lzd%Glr%wfd)

!!$            !finalize memory counting (there are still at least positions and the forces allocated)
!!$            call memocc(0,0,'count','stop')
!!$
!!$            if (nproc > 1) call MPI_FINALIZE(ierr)

        !test if stderr works
        write(0,*)'unnormal end'
        call mpibarrier(bigdft_mpi%mpi_comm)
        call f_err_throw('Convergence error, cannot proceed. '//&
             'Writing positions in file posfail.xyz',err_name='BIGDFT_RUNTIME_ERROR')

     else
        exit loop_cluster
     end if

  end do loop_cluster

  !preserve the previous value
  runObj%inputs%inputPsiId=inputPsiId_orig

  !put a barrier for all the processes
  call f_release_routine()
  call mpibarrier(bigdft_mpi%mpi_comm)

END SUBROUTINE call_bigdft

subroutine run_objects_update(runObj, dict)
  use bigdft_run, only: run_objects
  use dictionaries, only: dictionary, dict_update
  implicit none
  type(run_objects), intent(inout) :: runObj
  type(dictionary), pointer :: dict

  ! We merge the previous dictionnary with new entries.
  call dict_update(runObj%user_inputs, dict)

  ! Parse new dictionnary.
  call run_objects_parse(runObj)
END SUBROUTINE run_objects_update

!> Parse the input dictionary and create all run_objects
subroutine run_objects_parse(runObj)
  use module_base, only: bigdft_mpi,f_err_throw
  use module_interfaces, only: atoms_new, inputs_new, inputs_from_dict, create_log_file
  use dynamic_memory
  use dictionaries
  use module_atoms, only: deallocate_atoms_data
  use bigdft_run
  implicit none
  type(run_objects), intent(inout) :: runObj
  character(len=*), parameter :: subname = "run_objects_parse"

  ! Free potential previous inputs and atoms.
  if (associated(runObj%atoms)) then
     call deallocate_atoms_data(runObj%atoms) 
     deallocate(runObj%atoms)
  end if
  ! Allocate atoms_data structure
  call atoms_new(runObj%atoms)
  if (associated(runObj%inputs)) then
     call free_input_variables(runObj%inputs)
     deallocate(runObj%inputs)
  end if
  !Allocation input_variables structure and initialize it with default values
  call inputs_new(runObj%inputs)

  ! Regenerate inputs and atoms.
  call inputs_from_dict(runObj%inputs, runObj%atoms, runObj%user_inputs)

  ! Number of atoms should not change.
  if (runObj%rst%nat > 0 .and. runObj%rst%nat /= runObj%atoms%astruct%nat) then
     call f_err_throw("The number of atoms changed!",err_name='BIGDFT_RUNTIME_ERROR')
  else if (runObj%rst%nat == 0) then
     call restart_objects_set_nat(runObj%rst, runObj%atoms%astruct%nat)
  end if
  call restart_objects_set_mode(runObj%rst, runObj%inputs%inputpsiid)
  if (associated(runObj%rst)) then
     call release_material_acceleration(runObj%rst%GPU)
  end if
  call restart_objects_set_mat_acc(runObj%rst, bigdft_mpi%iproc, runObj%inputs%matacc)

  ! Generate radii
  call f_free_ptr(runObj%radii_cf)

  runObj%radii_cf = f_malloc_ptr((/ runObj%atoms%astruct%ntypes, 3 /), id="runObj%radii_cf")
  call read_radii_variables(runObj%atoms, runObj%radii_cf, &
       & runObj%inputs%crmult, runObj%inputs%frmult, runObj%inputs%projrad)

END SUBROUTINE run_objects_parse


subroutine run_objects_system_setup(runObj, iproc, nproc, rxyz, shift, mem)
  use module_base, only: gp,f_memcpy
  use bigdft_run
  use module_types
  use module_fragments
  use module_interfaces, only: system_initialization
  use psp_projectors
  use communications_base, only: deallocate_comms
  implicit none
  type(run_objects), intent(inout) :: runObj
  integer, intent(in) :: iproc, nproc
  real(gp), dimension(3,runObj%atoms%astruct%nat), intent(out) :: rxyz
  real(gp), dimension(3), intent(out) :: shift
  type(memory_estimation), intent(out) :: mem

  integer :: inputpsi, input_wf_format
  type(DFT_PSP_projectors) :: nlpsp
  type(system_fragment), dimension(:), pointer :: ref_frags
  character(len = *), parameter :: subname = "run_objects_estimate_memory"

  ! Copy rxyz since system_size() will shift them.
!!$  allocate(rxyz(3,runObj%atoms%astruct%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,rxyz,'rxyz',subname)
  call f_memcpy(src=runObj%atoms%astruct%rxyz,dest=rxyz)
  !call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, rxyz(1,1), 1)

  call system_initialization(iproc, nproc, .true., inputpsi, input_wf_format, .true., &
       & runObj%inputs, runObj%atoms, rxyz, runObj%rst%GPU%OCLconv, runObj%rst%KSwfn%orbs, &
       & runObj%rst%tmb%npsidim_orbs, runObj%rst%tmb%npsidim_comp, &
       & runObj%rst%tmb%orbs, runObj%rst%KSwfn%Lzd, runObj%rst%tmb%Lzd, &
       & nlpsp, runObj%rst%KSwfn%comms, shift, runObj%radii_cf, &
       & ref_frags)
  call MemoryEstimator(nproc,runObj%inputs%idsx,runObj%rst%KSwfn%Lzd%Glr,&
       & runObj%rst%KSwfn%orbs%norb,runObj%rst%KSwfn%orbs%nspinor,&
       & runObj%rst%KSwfn%orbs%nkpts,nlpsp%nprojel,&
       & runObj%inputs%nspin,runObj%inputs%itrpmax,runObj%inputs%iscf,mem)

  ! De-allocations
  call deallocate_Lzd_except_Glr(runObj%rst%KSwfn%Lzd)
  call deallocate_comms(runObj%rst%KSwfn%comms)
  call deallocate_orbs(runObj%rst%KSwfn%orbs)
  call free_DFT_PSP_projectors(nlpsp)
  call deallocate_locreg_descriptors(runObj%rst%KSwfn%Lzd%Glr)
  call nullify_locreg_descriptors(runObj%rst%KSwfn%Lzd%Glr)
END SUBROUTINE run_objects_system_setup
