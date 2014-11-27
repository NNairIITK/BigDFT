!> @file
!!  Define main module for using BigDFT as a blackbox
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG, DC)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module handling the object for the runs of bigDFT (restart, output, ...)
module bigdft_run
  use module_defs, only: gp
  use dictionaries
  use module_types, only: input_variables,DFT_wavefunction,GPU_pointers,energy_terms
  use module_atoms, only: atoms_data
  use dynamic_memory, only: f_reference_counter,f_ref_new,f_ref,f_unref,&
       nullify_f_ref
  use f_utils
  private

  !>  Used to restart a new DFT calculation or to save information 
  !!  for post-treatment
  type, public :: restart_objects
     type(f_reference_counter) :: refcnt
     integer :: version !< 0=cubic, 100=linear
     integer :: n1,n2,n3,nat
     real(gp), dimension(:,:), pointer :: rxyz_old,rxyz_new
     type(DFT_wavefunction) :: KSwfn !< Kohn-Sham wavefunctions
     type(DFT_wavefunction) :: tmb !<support functions for linear scaling
     type(GPU_pointers) :: GPU 
  end type restart_objects

  !> Public container to be used with bigdft_state().
  type, public :: run_objects
     type(f_enumerator) :: run_mode
     !> user input specifications
     type(dictionary), pointer :: user_inputs
     !> structure of BigDFT input variables
     type(input_variables), pointer    :: inputs
     !> datatype describing the atomic system.
     type(atoms_data), pointer         :: atoms
     !> datatype describing the wavefunctions objects
     type(restart_objects), pointer    :: rst
  end type run_objects


  !> Used to store results of a DFT calculation.
  type, public :: state_properties
     real(gp) :: energy, fnoise, pressure      !< Total energy, noise over forces and pressure
     type(energy_terms) :: energs              !< All energy terms
     integer :: fdim                           !< Dimension of allocated forces (second dimension)
     real(gp), dimension(:,:), pointer :: fxyz !< Atomic forces
     real(gp), dimension(6) :: strten          !< Stress Tensor
  end type state_properties

  public :: init_state_properties,deallocate_state_properties,restart_objects_set_mat_acc
  public :: run_objects_free,copy_state_properties,restart_objects_set_mode
  public :: nullify_run_objects,restart_objects_set_nat,nullify_restart_objects
  public :: run_objects_associate,init_restart_objects,bigdft_set_rxyz
  public :: state_properties_set_from_dict,free_restart_objects,bigdft_get_rxyz_ptr
  public :: run_objects_init,bigdft_init,bigdft_command_line_options,bigdft_nruns
  public :: bigdft_nat,bigdft_state,free_run_objects,set_run_objects,bigdft_run_new
  public :: release_run_objects,bigdft_get_cell,bigdft_get_geocode,bigdft_get_run_properties
  public :: bigdft_get_astruct_ptr,bigdft_write_atomic_file,bigdft_set_run_properties
  public :: bigdft_norb,bigdft_get_eval,bigdft_run_id_toa,bigdft_get_rxyz
  public :: bigdft_dot,bigdft_nrm2

!!$  ! interfaces of external routines 
!!$  interface
!!$     subroutine geopt(runObj,outs,nproc,iproc,ncount_bigdft)
!!$       use module_base
!!$       use module_types
!!$       implicit none
!!$       type(run_objects), intent(inout) :: runObj
!!$       type(state_properties), intent(inout) :: outs
!!$       integer, intent(in) :: nproc,iproc
!!$       integer, intent(inout) :: ncount_bigdft
!!$     END SUBROUTINE geopt
!!$
!!$  end interface
  

  contains
    
    !> All in one routine to initialise and set-up restart objects.
    !! in case of previously initialized structure
    subroutine init_restart_objects(iproc,inputs,atoms,rst)
      implicit none
      !Arguments
      integer, intent(in) :: iproc
      type(input_variables), intent(in) :: inputs
      type(atoms_data), intent(in) :: atoms
      type(restart_objects), intent(inout) :: rst

      !call restart_objects_new(rst)
      ! Number of atoms should not change durung the calculation 
      if (rst%nat > 0 .and. rst%nat /= atoms%astruct%nat) then
         call f_err_throw("The number of atoms changed!",&
              err_name='BIGDFT_RUNTIME_ERROR')
      else if (rst%nat == 0) then
         !create reference counter
         rst%refcnt=f_ref_new('rst')
         call restart_objects_set_nat(rst, atoms%astruct%nat)
      end if
      call restart_objects_set_mode(rst, inputs%inputpsiid)
      call release_material_acceleration(rst%GPU)
      call restart_objects_set_mat_acc(rst,iproc, inputs%matacc)
    END SUBROUTINE init_restart_objects

    !> Allocate and nullify restart objects
    pure subroutine nullify_restart_objects(rst)
      use module_defs, only: UNINITIALIZED
      use module_types, only: nullify_local_zone_descriptors,CUBIC_VERSION,nullify_paw_objects
      use locregs, only: nullify_locreg_descriptors
      use gaussians, only: nullify_gaussian_basis
      implicit none
      !Arguments
      type(restart_objects), intent(out) :: rst

      call nullify_f_ref(rst%refcnt)
      ! Decide whether we use the cubic or the linear version
      rst%version = UNINITIALIZED(CUBIC_VERSION)

      !allocate pointers
      rst%nat = 0
      nullify(rst%rxyz_new)
      nullify(rst%rxyz_old)

      !nullify unallocated pointers
      !here whe should define the nullification of the whole 
      !DFT wavefunction structure
      rst%KSwfn%c_obj = 0
      nullify(rst%KSwfn%psi)
      nullify(rst%KSwfn%orbs%eval)
      call nullify_paw_objects(rst%KSwfn%paw)
      call nullify_paw_objects(rst%tmb%paw)

      nullify(rst%KSwfn%gaucoeffs)
      nullify(rst%KSwfn%oldpsis)

      call nullify_locreg_descriptors(rst%KSwfn%Lzd%Glr)
      call nullify_gaussian_basis(rst%KSwfn%gbd)

      !Nullify LZD for cubic version (new input guess)
      call nullify_local_zone_descriptors(rst%tmb%lzd)

      !Nullify GPU data
      rst%GPU%OCLconv=.false.
    END SUBROUTINE nullify_restart_objects


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

      !check if the object can be freed
      call f_ref_free(rst%refcnt)

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


    !> Initialize the structure state_properties
    subroutine nullify_state_properties(outs)
      use module_defs, only: UNINITIALIZED
      use module_types, only: energy_terms_null
      implicit none
      type(state_properties), intent(out) :: outs

      outs%energs=energy_terms_null()
      outs%fdim      = 0
      nullify(outs%fxyz)
      outs%energy    = UNINITIALIZED(1.0_gp)
      outs%fnoise    = UNINITIALIZED(1.0_gp)
      outs%pressure  = UNINITIALIZED(1.0_gp)
      outs%strten(:) = UNINITIALIZED(1.0_gp)
    END SUBROUTINE nullify_state_properties


    subroutine init_state_properties(outs, nat)
      use module_base
      use dynamic_memory
      implicit none
      type(state_properties), intent(out) :: outs
      integer, intent(in) :: nat

      call nullify_state_properties(outs)
      outs%fdim = nat
      outs%fxyz = f_malloc_ptr((/ 3, outs%fdim /),id='outs%fxyz')
      outs%fxyz(:,:) = UNINITIALIZED(1.0_gp)
    END SUBROUTINE init_state_properties

    subroutine deallocate_state_properties(outs, fxyz)
      use module_base
      implicit none
      type(state_properties), intent(inout) :: outs
      real(gp), intent(out), optional :: fxyz

      if (associated(outs%fxyz)) then
         if (present(fxyz)) &
              call f_memcpy(src=outs%fxyz(1,1),dest=fxyz,n=3*outs%fdim)
         !call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz, 1)
         !end if
         call f_free_ptr(outs%fxyz)
      end if
    END SUBROUTINE deallocate_state_properties

    !> Copies outsA to outsB
    !! outsB has to be allocated before
    subroutine copy_state_properties(outsA,outsB)
      use module_base, only: f_err_throw,f_memcpy
      use yaml_strings, only: yaml_toa
      implicit none
      type(state_properties), intent(in) :: outsA
      type(state_properties), intent(inout) :: outsB

      if(outsA%fdim /= outsB%fdim)then
         call f_err_throw("Error in copy_state_properties: outsA and outsB have different sizes"//&
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
      call f_memcpy(src=outsA%strten,dest=outsB%strten)
    end subroutine copy_state_properties

    !> Associate to the structure run_objects, the input_variable structure and the atomic positions (atoms_data)
    subroutine run_objects_associate(runObj, inputs, atoms, rst, rxyz0)
      use module_base
      use module_types
      implicit none
      type(run_objects), intent(inout) :: runObj
      type(input_variables), intent(in), target :: inputs
      type(atoms_data), intent(in), target :: atoms
      type(restart_objects), intent(in), target :: rst
      real(gp), intent(inout), optional :: rxyz0 !<fake intent(in)

      !call run_objects_free_container(runObj)
      !associate only meaningful objects
      if (f_ref_count(rst%refcnt) >= 0) then
         if (associated(runObj%rst)) call f_unref(runObj%rst%refcnt)
         runObj%rst    => rst
         call f_ref(runObj%rst%refcnt)
      end if

      if (f_ref_count(atoms%refcnt) >= 0) then
         if (associated(runObj%atoms)) call f_unref(runObj%atoms%refcnt)
         runObj%atoms  => atoms
         call f_ref(runObj%atoms%refcnt)
      end if

      if (f_ref_count(inputs%refcnt) >= 0) then
         if (associated(runObj%inputs)) call f_unref(runObj%inputs%refcnt)
         runObj%inputs => inputs
         call f_ref(runObj%inputs%refcnt)
      end if

      if (present(rxyz0)) then
         call bigdft_set_rxyz(runObj,rxyz_add=rxyz0)
         !call vcopy(3 * atoms%astruct%nat, rxyz0, 1, runObj%atoms%astruct%rxyz(1,1), 1)
      end if

    END SUBROUTINE run_objects_associate

    !> default run properties
    subroutine bigdft_run_new(run)
      implicit none
      type(dictionary), pointer :: run
      
      run => dict_new('name' .is. 'input','posinp' .is. 'posinp')
    end subroutine bigdft_run_new

    !> set the parameters of the run 
    subroutine bigdft_set_run_properties(run,run_id,posinp)
      implicit none
      type(dictionary), pointer :: run
      character(len=*), intent(in), optional :: run_id
      character(len=*), intent(in), optional :: posinp

      if (present(run_id)) call set(run // 'name',trim(run_id))
      if (present(posinp)) call set(run // 'posinp',trim(posinp))

    end subroutine bigdft_set_run_properties

    !> get the parameters of the run 
    subroutine bigdft_get_run_properties(run,run_id,posinp)
      implicit none
      type(dictionary), pointer :: run
      character(len=*), intent(out), optional :: run_id
      character(len=*), intent(out), optional :: posinp

      if (present(run_id)) run_id = run // 'name'
      if (present(posinp)) posinp = run // 'posinp'
      
    end subroutine bigdft_get_run_properties

    function bigdft_run_id_toa()
      use yaml_output, only: yaml_toa
      use module_base, only: bigdft_mpi
      implicit none
      character(len=20) :: bigdft_run_id_toa

      bigdft_run_id_toa=repeat(' ',len(bigdft_run_id_toa))

      if (bigdft_mpi%ngroup>1) then
         bigdft_run_id_toa=adjustl(trim(yaml_toa(bigdft_mpi%igroup,fmt='(i15)')))
      end if

    end function bigdft_run_id_toa


    !> copy the atom position in runObject into a workspace
    !! or retrieve the positions from a file
    subroutine bigdft_get_rxyz(runObj,filename,rxyz_add,rxyz)
      use dynamic_memory, only: f_memcpy
      use yaml_strings, only: yaml_toa
      use module_atoms, only: atomic_structure,nullify_atomic_structure,&
           set_astruct_from_file,deallocate_atomic_structure
      use module_base, only: bigdft_mpi
      implicit none
      !> run object from which the atomic positions have to be 
      !! retrieved
      type(run_objects), intent(inout), optional :: runObj
      !> filename from which the atomic positions have to be read
      !! alternative to runObj
      character(len=*), intent(in), optional :: filename
      !>starting position of the atomic position.
      !! the user is responsible to guarantee that the 
      !! correct memory space is allocated thereafter
      real(gp), intent(inout), optional :: rxyz_add
      !> array of correct size (3,bigdft_nat(runObj)) with the atomic positions
      real(gp), dimension(:,:), intent(out), optional :: rxyz
      !local variables
      integer :: n
      type(atomic_structure) :: astruct

      if (present(runObj) .eqv. present(filename)) then
         call f_err_throw('Error in bigdft_get_rxyz: runObj *xor* filename'//&
              'should be present',err_name='BIGDFT_RUNTIME_ERROR')
      end if

      if (present(runObj)) then
         n=3*bigdft_nat(runObj)

         if (present(rxyz_add) .eqv. present(rxyz)) then
            call f_err_throw('Error in bigdft_get_rxyz: rxyz_add *xor* rxyz'//&
                 'should be present',err_name='BIGDFT_RUNTIME_ERROR')
         end if

         if (present(rxyz_add)) then
            call f_memcpy(n=n,dest=rxyz_add,src=runObj%atoms%astruct%rxyz(1,1))
         else if (present(rxyz)) then
            if (n /= size(rxyz)) then
               call f_err_throw('Error in bigdft_get_rxyz: wrong size ('//&
                    trim(yaml_toa(n))//' /= '//trim(yaml_toa(size(rxyz)))//&
                    ')',err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call f_memcpy(dest=rxyz,src=runObj%atoms%astruct%rxyz)
         end if
      else if (present(filename)) then
         call nullify_atomic_structure(astruct)
         call set_astruct_from_file(trim(filename), bigdft_mpi%iproc, astruct)
         n=3*astruct%nat
         if (present(rxyz_add)) then
            call f_memcpy(n=n,dest=rxyz_add,src=astruct%rxyz(1,1))
         else if (present(rxyz)) then
            if (n /= size(rxyz)) then
               call f_err_throw('Error in bigdft_get_rxyz: wrong size ('//&
                    trim(yaml_toa(n))//' /= '//trim(yaml_toa(size(rxyz)))//&
                    ')',err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call f_memcpy(dest=rxyz,src=astruct%rxyz)
         end if

       call deallocate_atomic_structure(astruct)
       
      end if
    end subroutine bigdft_get_rxyz

    !> returns the pointer to the atomic positions of the run.
    !! it performas a shallwo copy therefore this routine is intended to 
    !! provide acces to position for reading.
    !! Use at own risk to modify the value of the atomic positions.
    !! it returns nullified pointer in the case runObj is not properly
    !! initialized
    function bigdft_get_rxyz_ptr(runObj) result(rxyz)
      implicit none
      type(run_objects), intent(in) :: runObj
      real(gp), dimension(:,:), pointer :: rxyz

      nullify(rxyz)
      if (associated(runObj%atoms)) rxyz => runObj%atoms%astruct%rxyz

    end function bigdft_get_rxyz_ptr

    function bigdft_get_astruct_ptr(runObj) result(astruct)
      use module_atoms, only: atomic_structure
      implicit none
      type(run_objects), intent(in), target :: runObj
      type(atomic_structure), pointer :: astruct

      nullify(astruct)
      if (associated(runObj%atoms)) astruct => runObj%atoms%astruct
    end function bigdft_get_astruct_ptr

    !> Routine to dump the atomic file in normal BigDFT mode
    subroutine bigdft_write_atomic_file(runObj,outs,filename,comment,&
         cwd_path)
      use module_base, only: bigdft_mpi
      use module_atoms, only: astruct_dump_to_file
      implicit none
      type(run_objects), intent(in) :: runObj
      type(state_properties), intent(in) :: outs
      character(len=*), intent(in) :: filename,comment
      !> when present and true, output the file in the main 
      !! working directory (where input files are present)
      !! and not in dir_output
      logical, intent(in), optional :: cwd_path
      !local variables
      logical :: indir

      !only master node dumps on disk
      if (bigdft_mpi%iproc /= 0) return

      indir=.false.
      if(present(cwd_path)) indir=cwd_path

      if (indir) then
         call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
              trim(filename),comment,&
              energy=outs%energy,forces=outs%fxyz)
      else
         call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
              trim(runObj%inputs%dir_output)//trim(filename),comment,&
              energy=outs%energy,forces=outs%fxyz)
      end if
      
    end subroutine bigdft_write_atomic_file

    !> Import positions for the run object from a given array
    subroutine bigdft_set_rxyz(runObj,rxyz_add,rxyz)
      use dynamic_memory, only: f_memcpy
      use yaml_strings, only: yaml_toa
      use yaml_output
      implicit none
      type(run_objects), intent(inout) :: runObj
      !>starting position of the atomic position.
      !! the user is responsible to guarantee that the 
      !! correct memory space is allocated thereafter
      real(gp), intent(inout), optional :: rxyz_add
      !> array of correct size (3,bigdft_nat(runObj)) with the atomic positions
      real(gp), dimension(:,:), intent(in), optional :: rxyz
      !local variables
      integer :: n

      n=3*bigdft_nat(runObj)

      if (present(rxyz_add) .eqv. present(rxyz)) then
         call f_err_throw('Error in bigdft_set_rxyz: rxyz_add *xor* rxyz'//&
              'should be present',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (present(rxyz_add)) then
         call f_memcpy(n=n,src=rxyz_add,dest=runObj%atoms%astruct%rxyz(1,1))
      else if (present(rxyz)) then
         if (n /= size(rxyz)) then
            call f_err_throw('Error in bigdft_set_rxyz: wrong size ('//&
                 trim(yaml_toa(n))//' /= '//trim(yaml_toa(size(rxyz)))//&
                 ')',err_name='BIGDFT_RUNTIME_ERROR')
         end if
         call f_memcpy(src=rxyz,dest=runObj%atoms%astruct%rxyz)
      end if

    end subroutine bigdft_set_rxyz

    subroutine state_properties_set_from_dict(outs, dict)
      use dictionaries
      use public_keys, only: GOUT_ENERGY,GOUT_FORCES
      implicit none
      type(dictionary), pointer :: dict
      type(state_properties), intent(inout) :: outs

      integer :: i
      type(dictionary), pointer :: it,it0

      i=0
      it0 => dict_iter(dict .get. GOUT_FORCES)
      do while(associated(it0))
         !this will be done only once if the key exists
         if (.not. associated(outs%fxyz)) &
              call init_state_properties(outs, dict_len(dict // GOUT_FORCES))

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
!!$              & call init_state_properties(outs, dict_len(dict // GOUT_FORCES))
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
    end subroutine state_properties_set_from_dict

    !> Routines to handle the argument objects of bigdft_state().
    pure subroutine nullify_run_objects(runObj)
      use module_types
      use f_utils, only: f_enumerator_null
      implicit none
      type(run_objects), intent(out) :: runObj
      runObj%run_mode=f_enumerator_null()
      nullify(runObj%user_inputs)
      nullify(runObj%inputs)
      nullify(runObj%atoms)
      nullify(runObj%rst)
    END SUBROUTINE nullify_run_objects
    
    !>release run_objects structure as a whole
    !! if the reference counter goes to zero, do not free
    !! the structure as other atoms and inputs may live somewhere
    !! freeing command has to be given explicitly until we are 
    !! sure that it would be illegal to have inputs and atoms living
    !! separately from run_objects
    !! should not give exception as release might be called indefinitely
    subroutine release_run_objects(runObj)
      use module_types
      use module_base
      use module_atoms, only: deallocate_atoms_data
      implicit none
      type(run_objects), intent(inout) :: runObj
      !local variables
      integer :: count

      if (associated(runObj%rst)) then
         call f_unref(runObj%rst%refcnt,count=count)
         if (count==0) then
            call free_restart_objects(runObj%rst)
         else
            nullify(runObj%rst)
         end if
      end if
      if (associated(runObj%atoms)) then
         call f_unref(runObj%atoms%refcnt,count=count)
         if (count==0) then
            call deallocate_atoms_data(runObj%atoms) 
         else
            nullify(runObj%atoms)
         end if
      end if
      if (associated(runObj%inputs)) then
         call f_unref(runObj%inputs%refcnt,count=count)
         if (count==0) then
            call free_input_variables(runObj%inputs)
         else
            nullify(runObj%inputs)
         end if
            
      end if
      call nullify_run_objects(runObj)
    end subroutine release_run_objects

    !> Free the run_objects structure, &
    !! if all the objects have been dereferenced
    subroutine free_run_objects(runObj)
      use module_types
      use module_base
      use module_atoms, only: deallocate_atoms_data
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
      call nullify_run_objects(runObj)
    END SUBROUTINE free_run_objects

    !> Read all input files and create the objects to run BigDFT
    subroutine run_objects_init(runObj,run_dict,source)
      use module_base, only: bigdft_mpi,dict_init
      use module_types
      use module_input_dicts, only: user_dict_from_files
      implicit none
      !> Object for BigDFT run. Has to be initialized by this routine in order to
      !! call bigdft main routine.
      type(run_objects), intent(out) :: runObj
      !> dictionary containing instructions for the run.
      !! identifies the instructions needed to initialize the run.
      !! essentially, run id, input positions id and so on.
      !! this dictionary is an element of the list that can be found 
      !! in options // 'BigDFT' after the call to bigdft_init routine
      !! if absent, the run objects need src argument to be initialized
      !! otherwise the runObj is a empty structure
      type(dictionary), pointer, optional :: run_dict
      !> template to perform a shallow copy of the arguments already inside.
      !! if run_dict is absent
      !! if present together with run_dict, it performs a shallow copy of the restart variables 
      !! given by it *except* the structures inputs and atoms, 
      !! which are provided by the informations given by run_dict
      type(run_objects), intent(in), optional :: source
      !local variables
      character(len=max_field_length) :: radical, posinp

      call nullify_run_objects(runObj)

      if (present(run_dict)) then
         radical = run_dict // 'name'
         posinp = run_dict // 'posinp'
         
         !here the control of the logfile can be inserted, driven by run_dict and not anymore by 
         ! user_inputs
 
         ! Generate input dictionary and parse it.
         call dict_init(runObj%user_inputs)
         call user_dict_from_files(runObj%user_inputs, radical, posinp, bigdft_mpi)

         !this will fill atoms and inputs
         call set_run_objects(runObj)

         !the user input is not needed anymore
         call dict_free(runObj%user_inputs)

         !decide what to do with restart
         if (present(source)) then
            if (associated(runObj%rst)) call f_unref(runObj%rst%refcnt)
            runObj%rst => source%rst
            call f_ref(runObj%rst%refcnt)
            !check the restart coherence
            ! Number of atoms should not change during the calculation 
            if (runObj%rst%nat > 0 .and. &
                 runObj%rst%nat /= runObj%atoms%astruct%nat) then
               call f_err_throw("The number of atoms changed!",&
                    err_name='BIGDFT_RUNTIME_ERROR')
            end if
         else
            ! Allocate persistent structures.
            allocate(runObj%rst)
            call nullify_restart_objects(runObj%rst)
            !init and update the restart objects
            call init_restart_objects(bigdft_mpi%iproc,&
                 runObj%inputs,runObj%atoms,runObj%rst)
         end if
         ! Start the signaling loop in a thread if necessary.
         if (runObj%inputs%signaling .and. bigdft_mpi%iproc == 0) then
            call bigdft_signals_init(runObj%inputs%gmainloop, 2, &
                 & runObj%inputs%domain, len_trim(runObj%inputs%domain))
            call bigdft_signals_start(runObj%inputs%gmainloop, runObj%inputs%signalTimeout)
         end if
         
      else if (present(source)) then
         call run_objects_associate(runObj,&
              source%inputs,source%atoms,source%rst)
      end if

    END SUBROUTINE run_objects_init

    subroutine bigdft_init(options)
      use yaml_parse
      use dictionaries
      !use yaml_output, only: yaml_map
      use yaml_strings, only: f_strcpy,yaml_toa
      use module_defs, only: bigdft_mpi
      use module_input_dicts, only: merge_input_file_to_dict
      use f_utils, only: f_file_exists
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

    !accessors for external programs
    !> Get the number of orbitals of the run in rst
    function bigdft_nat(runObj) result(nat)
      implicit none
      type(run_objects), intent(in) :: runObj !> BigDFT run structure
      integer :: nat !> Number of atoms

      if (associated(runObj%atoms)) then
         nat=runObj%atoms%astruct%nat
      else
         nat=-1
      end if
      if (f_err_raise(nat < 0 ,'Number of atoms uninitialized',&
           err_name='BIGDFT_RUNTIME_ERROR')) return

    end function bigdft_nat

    function bigdft_norb(runObj) result(norb)
      implicit none
      type(run_objects), intent(in) :: runObj
      integer :: norb

      norb=0
      if (associated(runObj%rst)) norb=runObj%rst%KSwfn%orbs%norb
      if (norb <= 0) call f_err_throw('Number of orbitals unitialized',&
           err_name='BIGDFT_RUNTIME_ERROR')
    end function bigdft_norb

    !> Fill the array eval with the number of orbitals of the last run
    !! the array eval should have been allocated with the correct size
    subroutine bigdft_get_eval(runObj,eval)
      use module_base, only: gp,f_memcpy
      implicit none
      type(run_objects), intent(in) :: runObj !< BigDFT run structure
      !> Buffer for eigenvectors. Should have at least dimension equal to 
      !! bigdft_norb(runObj)
      real(gp), dimension(*), intent(out) :: eval 
      !local variables
      integer :: norb

      norb=bigdft_norb(runObj)

      call f_memcpy(n=norb,src=runObj%rst%KSwfn%orbs%eval(1),dest=eval(1))
    end subroutine bigdft_get_eval


    function bigdft_get_geocode(runObj) result(geocode)
      implicit none
      type(run_objects), intent(in) :: runObj
      character :: geocode

      geocode=' '
      if (associated(runObj%atoms)) geocode=runObj%atoms%astruct%geocode
      if (all(geocode /= ['F','S','W','P'])) &
           call f_err_throw('Geometry code uninitialized',&
           err_name='BIGDFT_RUNTIME_ERROR')
            
    end function bigdft_get_geocode

    function bigdft_get_cell(runObj) result(cell)
      implicit none
      type(run_objects), intent(in) :: runObj
      real(gp), dimension(3) :: cell
      
      cell=0.0_gp
      if (associated(runObj%atoms)) then
         cell=runObj%atoms%astruct%cell_dim
      else
         call f_err_throw('Cell uninitialized',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end if
    end function bigdft_get_cell

!TODO: 1) define f_enumerator types
!!     2) regroup this routine in a QM_something routine
!!     3) Insert the if statements of mhgps in the bigdft_state objects
!!     4) Make everything compilable.

    !> Routine to use BigDFT as a blackbox
    subroutine bigdft_state(runObj,outs,infocode)
      use module_base
      use yaml_output
      use module_atoms, only: astruct_dump_to_file,rxyz_inside_box
      use module_types, only: INPUT_PSI_MEMORY_LINEAR,LINEAR_VERSION,INPUT_PSI_MEMORY_GAUSS, &
           INPUT_PSI_LCAO,INPUT_PSI_MEMORY_WVL,old_wavefunction_null,INPUT_PSI_LINEAR_AO,deallocate_wfd
      !use communications_base
      implicit none
      type(run_objects), intent(inout) :: runObj
      type(state_properties), intent(inout) :: outs
      integer, intent(inout) :: infocode
      !local variables
      character(len=*), parameter :: subname='bigdft_state'
      logical :: exists
      integer :: inputPsiId_orig,istep
      !integer :: iat
      real(gp) :: maxdiff
      external :: cluster,forces_via_finite_differences
      !put a barrier for all the processes
      call mpibarrier(bigdft_mpi%mpi_comm)

      call f_routine(id=subname)
      !Check the consistency between MPI processes of the atomic coordinates and broadcast them
      if (bigdft_mpi%nproc >1) then
         call mpibcast(runObj%atoms%astruct%rxyz,comm=bigdft_mpi%mpi_comm,&
              maxdiff=maxdiff)
         if (maxdiff > epsilon(1.0_gp)) then
            if (bigdft_mpi%iproc==0) then
               call yaml_warning('Input positions not identical! '//&
                    '(difference:'//trim(yaml_toa(maxdiff))//&
                    ' ), however broadcasting from master node.')
               call yaml_flush_document()
            end if
         end if
      end if
!!$      maxdiff=mpimaxdiff(runObj%atoms%astruct%rxyz,comm=bigdft_mpi%mpi_comm,bcast=.true.)
!!$      if (maxdiff > epsilon(1.0_gp)) then
!!$         if (bigdft_mpi%iproc==0) then
!!$            call yaml_warning('Input positions not identical! '//&
!!$                 '(difference:'//trim(yaml_toa(maxdiff))//' ), broadcasting from master node.')
!!$            call yaml_comment('If the code hangs here, this means that not all the tasks met the threshold')
!!$            call yaml_comment('This might be related to arithmetics in performing the comparison')
!!$            call yaml_flush_document()
!!$         end if
!!$         !the check=.true. is important here: it controls that each process
!!$         !will participate in the broadcasting
!!$         call mpibcast(runObj%atoms%astruct%rxyz,comm=bigdft_mpi%mpi_comm,&
!!$              check=.true.)
!!$      end if

      !fill the rxyz array with the positions
      !wrap the atoms in the periodic directions when needed
      call rxyz_inside_box(runObj%atoms%astruct,rxyz=runObj%rst%rxyz_new)

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
!write(*,*)'(BIGDFTbastian) debug after free memory,runObj%inputs%inputPsiId',runObj%inputs%inputPsiId,bigdft_mpi%iproc

         !backdoor for hacking, finite difference method for calculating forces on particular quantities
         call f_file_exists('input.finite_difference_forces',exists)
         if (exists) then
            runObj%inputs%last_run=1 !do the last_run things nonetheless
            runObj%inputs%inputPsiId = INPUT_PSI_LCAO !the first run always restart from IG
            !experimental_modulebase_var_onlyfion=.true. !put only ionic forces in the forces
         end if

         !Main routine calculating the KS orbitals
         call cluster(bigdft_mpi%nproc,bigdft_mpi%iproc,runObj%atoms,runObj%rst%rxyz_new, &
              outs%energy, outs%energs, outs%fxyz, outs%strten, outs%fnoise, outs%pressure,&
              runObj%rst%KSwfn,runObj%rst%tmb,&
              runObj%rst%rxyz_old,runObj%inputs,runObj%rst%GPU,infocode)

         !save the new atomic positions in the rxyz_old array
         call f_memcpy(src=runObj%rst%rxyz_new,dest=runObj%rst%rxyz_old)

         if (exists) then
            call forces_via_finite_differences(bigdft_mpi%iproc,bigdft_mpi%nproc,runObj%atoms,runObj%inputs, &
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
            runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_WVL
            exit loop_cluster
         else if (runObj%inputs%inputPsiId == INPUT_PSI_LCAO .and. infocode==3) then
            if (bigdft_mpi%iproc == 0) then
               call astruct_dump_to_file(runObj%atoms%astruct,"posfail",&
                    'UNCONVERGED WF ',outs%energy,runObj%rst%rxyz_new)
            end if

            call f_free_ptr(runObj%rst%KSwfn%psi)
            call f_free_ptr(runObj%rst%KSwfn%orbs%eval)

            call deallocate_wfd(runObj%rst%KSwfn%Lzd%Glr%wfd)

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

    END SUBROUTINE bigdft_state


    !> Parse the input dictionary and create all run_objects
    !! in particular this routine identifies the input and the atoms structure
    subroutine set_run_objects(runObj)
      use module_base, only: f_err_throw
      use module_interfaces, only: atoms_new, inputs_new, inputs_from_dict, create_log_file
      use module_atoms, only: deallocate_atoms_data
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

!!$  ! Number of atoms should not change durung the calculation 
!!$  if (runObj%rst%nat > 0 .and. runObj%rst%nat /= runObj%atoms%astruct%nat) then
!!$     call f_err_throw("The number of atoms changed!",err_name='BIGDFT_RUNTIME_ERROR')
!!$  else if (runObj%rst%nat == 0) then
!!$     call restart_objects_set_nat(runObj%rst, runObj%atoms%astruct%nat)
!!$  end if
!!$  call restart_objects_set_mode(runObj%rst, runObj%inputs%inputpsiid)
!!$  call release_material_acceleration(runObj%rst%GPU)
!!$  call restart_objects_set_mat_acc(runObj%rst, bigdft_mpi%iproc, runObj%inputs%matacc)

    END SUBROUTINE set_run_objects

    !> performs the scalar product between two
    !! set of atomic positions 
    !! The results is considered only in non-frozen directions
    function bigdft_dot(runObj,dx,dy_add,dy,dx_add) result(scpr)
      implicit none
      !> run_object bigdft structure
      type(run_objects), intent(in) :: runObj
      !> x vector of atomic positions. Has to be given in alternative to 
      !! dx_add
      real(gp), dimension(3,runObj%atoms%astruct%nat), intent(in), optional :: dx
      !> position of the first element of x vector. 
      !! Has to be given in alternative to dx_add
      real(gp), optional :: dx_add
      !> y vector of atomic positions. Has to be given in alternative to 
      !! dy_add
      real(gp), dimension(3,runObj%atoms%astruct%nat), intent(in), optional :: dy
      !> position of the first element of y vector. 
      !! Has to be given in alternative to dy_add
      real(gp), optional :: dy_add
      !> result of the ddot procedure
      real(gp) :: scpr
      external :: atomic_dot

      if (present(dx) .eqv. present(dx_add)) then
         call f_err_throw('Error in bigdft_dot: dx *xor* dx_add have to be present',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (present(dy) .eqv. present(dy_add)) then
         call f_err_throw('Error in bigdft_dot: dy *xor* dy_add have to be present',err_name='BIGDFT_RUNTIME_ERROR')
      end if

      if (present(dx)) then
         if (present(dy)) then
            call atomic_dot(bigdft_get_astruct_ptr(runObj),dx,dy,scpr)
         else
            call atomic_dot(bigdft_get_astruct_ptr(runObj),dx,dy_add,scpr)
         end if
      else
         if (present(dy)) then
            call atomic_dot(bigdft_get_astruct_ptr(runObj),dx_add,dy,scpr)
         else
            call atomic_dot(bigdft_get_astruct_ptr(runObj),dx_add,dy_add,scpr)
         end if
      end if
    end function bigdft_dot
  
    !> square root of bigdft_dot(runObj,dx,dx)
    function bigdft_nrm2(runObj,dx,dx_add) result(scpr)
      use yaml_output, only: yaml_toa
      implicit none
      !> run_object bigdft structure
      type(run_objects), intent(in) :: runObj
      !> x vector of atomic positions. Has to be given in alternative to 
      !! dx_add
      real(gp), dimension(3,runObj%atoms%astruct%nat), intent(in), optional :: dx
      !> position of the first element of x vector. 
      !! Has to be given in alternative to dx_add
      real(gp), optional :: dx_add
      !> result of the dnrm2 procedure
      real(gp) :: scpr

      if (present(dx) .eqv. present(dx_add)) then
         call f_err_throw('Error in bigdft_nrm2: dx *xor* dx_add have to be present',err_name='BIGDFT_RUNTIME_ERROR')
      end if

      if (present(dx)) then
         scpr=bigdft_dot(runObj,dx=dx,dy=dx)
      else
         scpr=bigdft_dot(runObj,dx_add=dx_add,dy_add=dx_add)
      end if
      if (scpr <=0.0_gp) call f_err_throw(&
           'Meaningless result in bigdft_nrm2! (<x|x>='//&
           trim(yaml_toa(scpr))//')',err_name='BIGDFT_RUNTIME_ERROR')
      
      scpr=sqrt(scpr)
    end function bigdft_nrm2

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

subroutine run_objects_update(runObj, dict)
  use module_base, only: bigdft_mpi
  use bigdft_run, only: run_objects,init_restart_objects,set_run_objects
  use dictionaries!, only: dictionary, dict_update,dict_copy,dict_free,dict_iter,dict_next
  use yaml_output
  implicit none
  type(run_objects), intent(inout) :: runObj
  type(dictionary), pointer :: dict
  !local variables
  type(dictionary), pointer :: item

  if (associated(runObj%user_inputs)) then
     item => dict_iter(dict)
     do while (associated(item))
        if (index(dict_key(item),'psppar') == 1) then
           call dict_copy(runObj%user_inputs//trim(dict_key(item)),item)
        end if
        item => dict_next(item)
     end do
  end if

  ! We merge the previous dictionary with new entries.
  call dict_update(runObj%user_inputs, dict)

  ! Parse new dictionary.
  call set_run_objects(runObj)

  !init and update the restart objects
  call init_restart_objects(bigdft_mpi%iproc,runObj%inputs,runObj%atoms,&
       runObj%rst)
END SUBROUTINE run_objects_update


!> this routine should be used in memguess executable also
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
       & nlpsp, runObj%rst%KSwfn%comms, shift, &
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
