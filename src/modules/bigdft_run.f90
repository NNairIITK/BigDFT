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
  use f_refcnts, only: f_reference_counter,f_ref_new,f_ref,f_unref,&
       nullify_f_ref,f_ref_free
  use f_utils
  use f_enums, f_str => str
  use module_input_dicts, only: bigdft_set_run_properties => dict_set_run_properties,&
       bigdft_get_run_properties => dict_get_run_properties
  use public_enums
  private

  !>  Used to restart a new DFT calculation or to save information
  !!  for post-treatment
  type, public :: QM_restart_objects
     type(f_reference_counter) :: refcnt
     integer :: version !< 0=cubic, 100=linear
     integer :: n1,n2,n3,nat
     real(gp), dimension(:,:), pointer :: rxyz_old,rxyz_new
     type(DFT_wavefunction) :: KSwfn !< Kohn-Sham wavefunctions
     type(DFT_wavefunction) :: tmb !<support functions for linear scaling
     type(GPU_pointers) :: GPU
  end type QM_restart_objects

  !>supplementary type in run_objects
  type, public :: MM_restart_objects
     type(f_reference_counter) :: refcnt
     !> array for temporary copy of atomic positions and forces
     real(gp), dimension(:,:), pointer :: rf_extra
     type(f_enumerator) :: run_mode !< run_mode for freeing the extra treatments
  end type MM_restart_objects

  !> Public container to be used with bigdft_state().
  type, public :: run_objects
     type(f_enumerator), pointer :: run_mode
     !> user input specifications
     type(dictionary), pointer :: user_inputs
     !> structure of BigDFT input variables
     type(input_variables), pointer    :: inputs
     !> datatype describing the atomic system.
     type(atoms_data), pointer         :: atoms
     !> datatype describing the wavefunctions objects
     type(QM_restart_objects), pointer    :: rst
     !> datatype describing extra MM information
     type(MM_restart_objects), pointer :: mm_rst
  end type run_objects

  !> Used to store results of a DFT calculation.
  type, public :: state_properties
     real(gp) :: energy, fnoise, pressure      !< Total energy, noise over forces and pressure
     type(energy_terms) :: energs              !< All energy terms
     integer :: fdim                           !< Dimension of allocated forces (second dimension)
     real(gp), dimension(:,:), pointer :: fxyz !< Atomic forces
     real(gp), dimension(6) :: strten          !< Stress Tensor
  end type state_properties

  public :: init_state_properties,deallocate_state_properties
  public :: run_objects_free,copy_state_properties
  public :: nullify_run_objects
  public :: run_objects_associate,bigdft_set_rxyz
  public :: state_properties_set_from_dict,bigdft_get_rxyz_ptr
  public :: run_objects_init,bigdft_init,bigdft_command_line_options,bigdft_nruns
  public :: init_QM_restart_objects,init_MM_restart_objects,set_run_objects,nullify_QM_restart_objects
  public :: nullify_MM_restart_objects
  public :: bigdft_nat,bigdft_state,free_run_objects
  public :: release_run_objects,bigdft_get_cell,bigdft_get_cell_ptr,bigdft_get_geocode,bigdft_get_run_properties
  public :: bigdft_get_units, bigdft_set_units
  public :: bigdft_get_astruct_ptr,bigdft_write_atomic_file,bigdft_set_run_properties
  public :: bigdft_norb,bigdft_get_eval,bigdft_run_id_toa,bigdft_get_rxyz
  public :: bigdft_dot,bigdft_nrm2
  public :: bigdft_get_input_policy
  public :: bigdft_set_input_policy

  !> Input policies
  integer,parameter,public :: INPUT_POLICY_SCRATCH = 10000 !< Start the calculation from scratch
  integer,parameter,public :: INPUT_POLICY_MEMORY  = 10001 !< Start the calculation from data in memory
  integer,parameter,public :: INPUT_POLICY_DISK    = 10002 !< Start the calculation from data on disk

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
  subroutine init_QM_restart_objects(iproc,inputs,atoms,rst)
    implicit none
    !Arguments
    integer, intent(in) :: iproc
    type(input_variables), intent(in) :: inputs
    type(atoms_data), intent(in) :: atoms
    type(QM_restart_objects), intent(inout) :: rst

    !call QM_restart_objects_new(rst)
    ! Number of atoms should not change during the calculation
    if (rst%nat > 0 .and. rst%nat /= atoms%astruct%nat) then
       call f_err_throw("The number of atoms changed!",&
            err_name='BIGDFT_RUNTIME_ERROR')
    else if (rst%nat == 0) then
       !create reference counter
       rst%refcnt=f_ref_new('rst')
       call QM_restart_objects_set_nat(rst, atoms%astruct%nat)
    end if
    call QM_restart_objects_set_mode(rst, inputs%inputpsiid)
    call release_material_acceleration(rst%GPU)
    call QM_restart_objects_set_mat_acc(rst,iproc, inputs%matacc)
  END SUBROUTINE init_QM_restart_objects

  pure subroutine nullify_MM_restart_objects(mm_rst)
    implicit none
    type(MM_restart_objects), intent(out) :: mm_rst
    call nullify_f_ref(mm_rst%refcnt)
    nullify(mm_rst%rf_extra)
    call nullify_f_enum(mm_rst%run_mode)
  end subroutine nullify_MM_restart_objects

  !> fill the run_mode with the input enumerator
  !! if the treatment requires some extra allocations, control
  !! that the arrays are present and with the good shape
  !! otherwise allocate them
  subroutine init_MM_restart_objects(runObj,mm_rst,nat,run_mode)
    use f_utils
    use dynamic_memory
    use public_enums
    use module_morse_bulk
    use module_tersoff
    use module_BornMayerHugginsTosiFumi
    use module_lj
    use module_lenosky_si
    use module_cp2k
    use yaml_output
    implicit none
    type(run_objects), intent(inout) :: runObj
    type(f_enumerator), intent(in) :: run_mode
    integer, intent(in) :: nat
    type(MM_restart_objects), intent(inout) :: mm_rst

    !then check if extra workspaces have to be allocated
    select case(trim(f_str(run_mode)))
    case('LENNARD_JONES_RUN_MODE')
       call nullify_MM_restart_objects(mm_rst)
       !create reference counter
       mm_rst%refcnt=f_ref_new('mm_rst')
       call init_lj(runObj%inputs%mm_paramset,&
            runObj%inputs%mm_paramfile,runObj%atoms%astruct%units)
    case('LENOSKY_SI_CLUSTERS_RUN_MODE')
       if (associated(mm_rst%rf_extra)) then
          if (size(mm_rst%rf_extra) == nat) then
             call f_zero(mm_rst%rf_extra)
          else
             call f_free_ptr(mm_rst%rf_extra)
             mm_rst%rf_extra=f_malloc0_ptr([3,nat],id='rf_extra')
          end if
       else
          call nullify_MM_restart_objects(mm_rst)
          !create reference counter
          mm_rst%refcnt=f_ref_new('mm_rst')
          mm_rst%rf_extra=f_malloc0_ptr([3,nat],id='rf_extra')
       end if
       call init_lensic(runObj%inputs%mm_paramset,&
            runObj%inputs%mm_paramfile,runObj%atoms%astruct%geocode,&
            runObj%atoms%astruct%units)
    case('LENOSKY_SI_BULK_RUN_MODE')
       if (associated(mm_rst%rf_extra)) then
          if (size(mm_rst%rf_extra) == nat) then
             call f_zero(mm_rst%rf_extra)
          else
             call f_free_ptr(mm_rst%rf_extra)
             mm_rst%rf_extra=f_malloc0_ptr([3,nat],id='rf_extra')
          end if
       else
          call nullify_MM_restart_objects(mm_rst)
          !create reference counter
          mm_rst%refcnt=f_ref_new('mm_rst')
          mm_rst%rf_extra=f_malloc0_ptr([3,nat],id='rf_extra')
       end if
    case('MORSE_SLAB_RUN_MODE')
       call nullify_MM_restart_objects(mm_rst)
       !create reference counter
       mm_rst%refcnt=f_ref_new('mm_rst')
       call init_morse_slab(runObj%inputs%mm_paramset,&
            runObj%inputs%mm_paramfile,runObj%atoms%astruct%geocode,&
            runObj%atoms%astruct%units)
    case('MORSE_BULK_RUN_MODE')
       call nullify_MM_restart_objects(mm_rst)
       !create reference counter
       mm_rst%refcnt=f_ref_new('mm_rst')
        call init_morse_bulk(runObj%inputs%mm_paramset,&
             runObj%inputs%mm_paramfile,runObj%atoms%astruct%geocode)
    case('TERSOFF_RUN_MODE')
       call nullify_MM_restart_objects(mm_rst)
       !create reference counter
       mm_rst%refcnt=f_ref_new('mm_rst')
       call init_tersoff(nat,runObj%atoms%astruct,runObj%inputs%mm_paramset,&
            runObj%inputs%mm_paramfile,runObj%atoms%astruct%geocode) 
    case('BMHTF_RUN_MODE')
       call nullify_MM_restart_objects(mm_rst)
       !create reference counter
       mm_rst%refcnt=f_ref_new('mm_rst')
       call init_bmhtf(nat,runObj%atoms%astruct,runObj%inputs%mm_paramset,&
            runObj%inputs%mm_paramfile,runObj%atoms%astruct%geocode) 
    case('AMBER_RUN_MODE')
       if (associated(mm_rst%rf_extra)) then
          if (size(mm_rst%rf_extra) == nat) then
             call f_zero(mm_rst%rf_extra)
          else
             call f_free_ptr(mm_rst%rf_extra)
             mm_rst%rf_extra=f_malloc0_ptr([3,nat],id='rf_extra')
          end if
       else
          call nullify_MM_restart_objects(mm_rst)
          !create reference counter
          mm_rst%refcnt=f_ref_new('mm_rst')
          mm_rst%rf_extra=f_malloc0_ptr([3,nat],id='rf_extra')
       endif
       call nab_init()
    case('CP2K_RUN_MODE')
       call nullify_MM_restart_objects(mm_rst)
       !create reference counter
       mm_rst%refcnt=f_ref_new('mm_rst')
       call init_cp2k(runObj%inputs%mm_paramfile,runObj%atoms%astruct%geocode)
    case default
       call nullify_MM_restart_objects(mm_rst)
       !create reference counter
       mm_rst%refcnt=f_ref_new('mm_rst')
    end select

    mm_rst%run_mode=run_mode
  end subroutine init_MM_restart_objects

  subroutine free_MM_restart_objects(mm_rst)
    use dynamic_memory
    use yaml_output
    use module_cp2k
    use module_BornMayerHugginsTosiFumi
    use f_enums, enum_int => int
    use yaml_strings
    implicit none
    type(MM_restart_objects), intent(inout) :: mm_rst
    !check if the object can be freed
    call f_ref_free(mm_rst%refcnt)
    call f_free_ptr(mm_rst%rf_extra)
    !free the extra variables
    select case(trim(f_str(mm_rst%run_mode)))
    case('BMHTF_RUN_MODE')
       call finalize_bmhtf()
    case('CP2K_RUN_MODE') ! CP2K run mode
       call finalize_cp2k()
    case('DFTBP_RUN_MODE') ! DFTB+ run mode
    case('LENNARD_JONES_RUN_MODE')
    case('MORSE_SLAB_RUN_MODE')
    case('MORSE_BULK_RUN_MODE')
    case('TERSOFF_RUN_MODE')
    case('LENOSKY_SI_CLUSTERS_RUN_MODE')
    case('LENOSKY_SI_BULK_RUN_MODE')
    case('AMBER_RUN_MODE')
    case('QM_RUN_MODE')
    case default
       call f_err_throw('Following method for evaluation of '//&
            'energies and forces is unknown: '+enum_int(mm_rst%run_mode),&
            err_name='BIGDFT_RUNTIME_ERROR')
    end select

    
  end subroutine free_MM_restart_objects

  !> Allocate and nullify restart objects
  pure subroutine nullify_QM_restart_objects(rst)
    use module_defs, only: UNINITIALIZED
    use module_types, only: nullify_local_zone_descriptors,nullify_paw_objects
    use locregs, only: nullify_locreg_descriptors
    use gaussians, only: nullify_gaussian_basis
    implicit none
    !Arguments
    type(QM_restart_objects), intent(out) :: rst

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
  END SUBROUTINE nullify_QM_restart_objects

  !pure 
  subroutine QM_restart_objects_set_mode(rst, inputpsiid)
    implicit none
    type(QM_restart_objects), intent(inout) :: rst
    type(f_enumerator), intent(in) :: inputpsiid !,integer

    if (inputpsiid .hasattr. 'CUBIC') then
       rst%version = CUBIC_VERSION
    else if (inputpsiid .hasattr. 'LINEAR') then
       rst%version = LINEAR_VERSION
    end if
  END SUBROUTINE QM_restart_objects_set_mode

  subroutine QM_restart_objects_set_nat(rst, nat)
    use module_base
    implicit none
    !Arguments
    integer, intent(in) :: nat
    type(QM_restart_objects), intent(inout) :: rst

    call f_free_ptr(rst%rxyz_old)
    call f_free_ptr(rst%rxyz_new)

    rst%nat = nat
    rst%rxyz_new = f_malloc_ptr((/ 3, nat /),id='rst%rxyz_new')
    rst%rxyz_old = f_malloc_ptr((/ 3, nat /),id='rst%rxyz_old')
  END SUBROUTINE QM_restart_objects_set_nat

  subroutine QM_restart_objects_set_mat_acc(rst, iproc, matacc)
    use module_input_keys, only: material_acceleration
    implicit none
    !Arguments
    type(QM_restart_objects), intent(inout) :: rst
    integer, intent(in) :: iproc
    type(material_acceleration), intent(in) :: matacc
    !initialise the acceleration strategy if required
    call init_material_acceleration(iproc,matacc,rst%GPU)
  END SUBROUTINE QM_restart_objects_set_mat_acc

  !> De-Allocate QM_restart_objects
  subroutine free_QM_restart_objects(rst)
    use module_base
    use locregs
    use gaussians, only: deallocate_gwf
    use module_types
    implicit none
    type(QM_restart_objects) :: rst
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

  END SUBROUTINE free_QM_restart_objects


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

  !>clean the outs object with empty (but meaningful)
  !! values so that the structure can be used for optimization
  subroutine clean_state_properties(outs)
    use module_types, only: energy_terms_null
    implicit none
    type(state_properties), intent(inout) :: outs

    outs%energs=energy_terms_null()
    call f_zero(outs%fxyz)
    outs%energy    = 0.0_gp
    outs%fnoise    = 0.0_gp
    outs%pressure  = 0.0_gp
    call f_zero(outs%strten)

  end subroutine clean_state_properties

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

  !> broadcast the state properties definitions among all the processors
  !! important to preserve coherency with respect the minimizers
  subroutine broadcast_state_properties(outs)
    use module_base
    use yaml_output
    implicit none
    type(state_properties), intent(inout) :: outs
    !local variables
    real(gp), dimension(:), allocatable :: data
    real(gp) :: maxdiff

    if (bigdft_mpi%nproc == 1) return

    !allocate datatransfer buffer: energy, pressure, fnoise, fxyz, strten
    data = f_malloc(3+size(outs%fxyz)+size(outs%strten),id='data')
    !fill data
    data(1)=outs%energy
    data(2)=outs%pressure
    data(3)=outs%fnoise
    if (size(outs%fxyz)>0) call f_memcpy(n=size(outs%fxyz),src=outs%fxyz(1,1),dest=data(4))
    call f_memcpy(n=size(outs%strten),src=outs%strten(1),dest=data(4+size(outs%fxyz)))

    call mpibcast(data,comm=bigdft_mpi%mpi_comm,&
         maxdiff=maxdiff)
    if (maxdiff > epsilon(1.0_gp)) then
       if (bigdft_mpi%iproc==0) then
          call yaml_warning('State properties not identical! '//&
               '(difference:'//trim(yaml_toa(maxdiff,fmt='(1pe15.5)'))//&
               ' ), broadcasting from master node.')
          call yaml_flush_document()
       end if
    end if

    !copy back in the structure
    outs%energy=data(1)
    outs%pressure=data(2)
    outs%fnoise=data(3)
    if (size(outs%fxyz)>0) call f_memcpy(n=size(outs%fxyz),dest=outs%fxyz(1,1),src=data(4))
    call f_memcpy(n=size(outs%strten),dest=outs%strten(1),src=data(4+size(outs%fxyz)))

    call f_free(data)
  end subroutine broadcast_state_properties


  !> Associate to the structure run_objects, the input_variable structure and the atomic positions (atoms_data)
  subroutine run_objects_associate(runObj, inputs, atoms, rst, mm_rst, rxyz0)
    use module_base
    use module_types
    implicit none
    type(run_objects), intent(inout) :: runObj
    type(input_variables), intent(in), target :: inputs
    type(atoms_data), intent(in), target :: atoms
    type(QM_restart_objects), intent(in), target :: rst
    type(MM_restart_objects), intent(in), target :: mm_rst
    real(gp), intent(inout), optional :: rxyz0 !<fake intent(in)

    !associate only meaningful objects
    if (f_ref_count(rst%refcnt) >= 0) then
       if (associated(runObj%rst)) call f_unref(runObj%rst%refcnt)
       runObj%rst    => rst
       call f_ref(runObj%rst%refcnt)
    end if

    if (f_ref_count(mm_rst%refcnt) >= 0) then
       if (associated(runObj%mm_rst)) call f_unref(runObj%mm_rst%refcnt)
       runObj%mm_rst    => mm_rst
       call f_ref(runObj%mm_rst%refcnt)
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
       runObj%run_mode => runObj%inputs%run_mode
    end if

    if (present(rxyz0)) then
       call bigdft_set_rxyz(runObj,rxyz_add=rxyz0)
    end if

  END SUBROUTINE run_objects_associate

  !> copy the atom position in runObject into a workspace
  !! or retrieve the positions from a file
  subroutine bigdft_get_rxyz(runObj,filename,rxyz_add,rxyz,energy,disableTrans)
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
    real(gp), intent(out), optional :: energy
    logical, intent(in), optional :: disableTrans
    !local variables
    integer :: n
    type(atomic_structure) :: astruct

    if (present(runObj) .eqv. present(filename)) then
       call f_err_throw('Error in bigdft_get_rxyz: runObj *xor* filename '//&
            'should be present',err_name='BIGDFT_RUNTIME_ERROR')
    elseif (present(runObj) .and. present(energy)) then
       call f_err_throw('Error in bigdft_get_rxyz: energy must not be present '//&
            'if runObj is present',err_name='BIGDFT_RUNTIME_ERROR')
    end if

    if (present(runObj)) then
       n=3*bigdft_nat(runObj)

       if (present(rxyz_add) .eqv. present(rxyz)) then
          call f_err_throw('Error in bigdft_get_rxyz: rxyz_add *xor* rxyz '//&
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
       call set_astruct_from_file(trim(filename), bigdft_mpi%iproc, &
            astruct,energy=energy,disableTrans=disableTrans)
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
  !! it performs a shallow copy therefore this routine is intended to
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

  end subroutine state_properties_set_from_dict

  !> Routines to handle the argument objects of bigdft_state().
  pure subroutine nullify_run_objects(runObj)
    use module_types
    implicit none
    type(run_objects), intent(out) :: runObj
    nullify(runObj%run_mode)
    nullify(runObj%user_inputs)
    nullify(runObj%inputs)
    nullify(runObj%atoms)
    nullify(runObj%rst)
    nullify(runObj%mm_rst)
  END SUBROUTINE nullify_run_objects

  !> Release run_objects structure as a whole
  !! if the reference counter goes to zero, do not free
  !! the structure as other atoms and inputs may live somewhere
  !! freeing command has to be given explicitly until we are
  !! sure that it would be illegal to have inputs and atoms living
  !! separately from run_objects
  !! should not give exception as release might be called indefinitely
  subroutine release_run_objects(runObj)
    use module_base
    use module_atoms, only: deallocate_atoms_data
    use yaml_output, only: yaml_sequence_close
    use module_input_keys, only: free_input_variables
    implicit none
    type(run_objects), intent(inout) :: runObj
    !local variables
    integer :: count
    if (associated(runObj%run_mode)) then
      if (bigdft_mpi%iproc==0 .and. runObj%run_mode /= 'QM_RUN_MODE')&
           call yaml_sequence_close()
    end if

    if (associated(runObj%rst)) then
       call f_unref(runObj%rst%refcnt,count=count)
       if (count==0) then
          call free_QM_restart_objects(runObj%rst)
       else
          nullify(runObj%rst)
       end if
    end if
    if (associated(runObj%mm_rst)) then
       call f_unref(runObj%mm_rst%refcnt,count=count)
       if (count==0) then
          call free_MM_restart_objects(runObj%mm_rst)
       else
          nullify(runObj%mm_rst)
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


  !> Free the run_objects structure,
  !! if all the objects have been dereferenced
  subroutine free_run_objects(runObj)
    use module_base
    use module_atoms, only: deallocate_atoms_data
    use yaml_output, only: yaml_sequence_close
    use module_input_keys, only: free_input_variables
    implicit none
    type(run_objects), intent(inout) :: runObj
    if (associated(runObj%run_mode)) then
      if (bigdft_mpi%iproc==0 .and. runObj%run_mode /= 'QM_RUN_MODE')&
           call yaml_sequence_close()
    end if

    call dict_free(runObj%user_inputs)
    if (associated(runObj%rst)) then
       call free_QM_restart_objects(runObj%rst)
       deallocate(runObj%rst)
    end if
    if (associated(runObj%mm_rst)) then
       call free_MM_restart_objects(runObj%mm_rst)
       deallocate(runObj%mm_rst)
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


  !> Parse the input dictionary and create all run_objects
  !! in particular this routine identifies the input and the atoms structure
  subroutine set_run_objects(runObj)
    use module_base, only: f_err_throw
    use module_interfaces, only: inputs_new, atoms_new
    use module_atoms, only: deallocate_atoms_data
    use module_input_dicts, only: dict_run_validate
    use module_input_keys, only: inputs_from_dict,free_input_variables
    use dynamic_memory
    use yaml_output
    implicit none
    type(run_objects), intent(inout) :: runObj
    character(len=*), parameter :: subname = "run_objects_parse"

    call f_routine(id='set_run_objects')

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
    call dict_run_validate(runObj%user_inputs)

    call inputs_from_dict(runObj%inputs, runObj%atoms, runObj%user_inputs)

    !associate the run_mode
    runObj%run_mode => runObj%inputs%run_mode

    call f_release_routine()

  END SUBROUTINE set_run_objects

  !> Read all input files and create the objects to run BigDFT
  subroutine run_objects_init(runObj,run_dict,source)
    use module_base, only: bigdft_mpi,dict_init
    use module_types
    use module_input_dicts, only: create_log_file
    use module_input_keys, only: user_dict_from_files
    use yaml_output
    use dynamic_memory
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
    logical :: dict_from_files
    character(len=max_field_length) :: radical, posinp_id

    call f_routine(id='run_objects_init')

    call nullify_run_objects(runObj)

    if (present(run_dict)) then
       !here the control of the logfile can be inserted, driven by run_dict and
       ! not anymore by user_inputs
       call create_log_file(run_dict,dict_from_files)
       if (dict_from_files) then
          ! Generate input dictionary.
          call dict_copy(runObj%user_inputs, run_dict)
          call bigdft_get_run_properties(run_dict, input_id = radical, posinp_id = posinp_id)
          call user_dict_from_files(runObj%user_inputs, radical, posinp_id, bigdft_mpi)
       else
          runObj%user_inputs => run_dict
       end if

       !this will fill atoms and inputs by parsing the input dictionary.
       call set_run_objects(runObj)

       !the user input is not needed anymore
       if (dict_from_files) call dict_free(runObj%user_inputs)

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
          if (associated(runObj%mm_rst)) call f_unref(runObj%mm_rst%refcnt)
          runObj%mm_rst => source%mm_rst
          call f_ref(runObj%mm_rst%refcnt)

       else
          ! Allocate persistent structures.
          allocate(runObj%rst)
          call nullify_QM_restart_objects(runObj%rst)
          !init and update the restart objects
          call init_QM_restart_objects(bigdft_mpi%iproc,&
               runObj%inputs,runObj%atoms,runObj%rst)

          allocate(runObj%mm_rst)
          call nullify_MM_restart_objects(runObj%mm_rst)
          call init_MM_restart_objects(runObj,runObj%mm_rst,bigdft_nat(runObj),runObj%run_mode)
       end if
       ! Start the signaling loop in a thread if necessary.
       if (runObj%inputs%signaling .and. bigdft_mpi%iproc == 0) then
          call bigdft_signals_init(runObj%inputs%gmainloop, 2, &
               & runObj%inputs%domain, len_trim(runObj%inputs%domain))
          call bigdft_signals_start(runObj%inputs%gmainloop, runObj%inputs%signalTimeout)
       end if

    else if (present(source)) then
       call run_objects_associate(runObj,&
            source%inputs,source%atoms,source%rst,source%mm_rst)
    end if

    if (bigdft_mpi%iproc==0 .and. runObj%run_mode /= 'QM_RUN_MODE') &
         call yaml_sequence_open('Initializing '//trim(f_str(runObj%run_mode)))
    call f_release_routine()

  END SUBROUTINE run_objects_init

  subroutine bigdft_init(options, with_taskgroups)
    use yaml_parse
    use dictionaries
    !use yaml_output, only: yaml_map
    use yaml_strings, only: f_strcpy,yaml_toa
    use module_defs, only: bigdft_mpi
    use module_input_dicts, only: merge_input_file_to_dict,set_dict_run_file
    use f_utils, only: f_file_exists
    use dynamic_memory
    implicit none
    !> dictionary of the options of the run
    !! on entry, it contains the options for initializing
    !! on exit, it contains in the key "BigDFT", a list of the
    !! dictionaries of each of the run that the local instance of BigDFT
    !! code has to execute.
    !! if this argument is not present, the code is only initialized
    !! in its normal mode: no taskgroups and default values of radical and posinp
    type(dictionary), pointer, optional :: options
    logical, intent(in), optional :: with_taskgroups
    !local variables
    logical :: exist_list,posinp_name
    integer :: ierr,mpi_groupsize,iconfig
    character(len=max_field_length) :: posinp_id,run_id,err_msg
    integer, dimension(4) :: mpi_info
    type(dictionary), pointer :: dict_run,opts
    logical :: uset

    call f_routine(id='bigdft_init')

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
    uset = .true.
    if (present(with_taskgroups)) uset = with_taskgroups

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
       if (uset) then
          do iconfig=1,bigdft_mpi%ngroup
             if (bigdft_mpi%ngroup > 1 .and. trim(run_id) /= "input") then
                call add(dict_run,trim(run_id) // trim(adjustl(yaml_toa(iconfig,fmt='(i3)'))))
             else
                call add(dict_run,trim(run_id))
             end if
          end do
       else
          call add(dict_run,trim(run_id))
       end if
    end if

    !call yaml_map('Dict of runs',dict_run)

    if (present(options)) then
       if (.not. associated(options)) call dict_init(options)
       !here the dict_run is given, and in each of the taskgroups a list of
       !runs for BigDFT instances has to be given
       do iconfig=0,dict_len(dict_run)-1
          if (modulo(iconfig,bigdft_mpi%ngroup)==bigdft_mpi%igroup .or. .not. uset) then
             run_id=dict_run//iconfig
             call set_dict_run_file(run_id,options)
          end if
       end do
    end if

    call dict_free(dict_run)

    call f_release_routine()

  end subroutine bigdft_init

  !>identify the options from command line
  !! and write the result in options dict
  subroutine bigdft_command_line_options(options)
    use yaml_parse
    use dictionaries
    use module_input_dicts, only: bigdft_options
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
  function bigdft_nat(runObj,filename) result(nat)
    use module_atoms, only: atomic_structure,nullify_atomic_structure,&
         set_astruct_from_file,deallocate_atomic_structure
    use module_base, only: bigdft_mpi
    implicit none
    type(run_objects), intent(in),optional :: runObj !> BigDFT run structure
    character(len=*), intent(in), optional :: filename
    integer :: nat !> Number of atoms
    !local
    type(atomic_structure) :: astruct

    if (present(runObj) .eqv. present(filename)) then
       call f_err_throw('Error in bigdft_nat: runObj *xor* filename '//&
            'should be present',err_name='BIGDFT_RUNTIME_ERROR')
    end if


    if(present(runObj))then
       if (associated(runObj%atoms)) then
          nat=runObj%atoms%astruct%nat
       else
          nat=-1
       end if
       if (f_err_raise(nat < 0 ,'Number of atoms uninitialized',&
            err_name='BIGDFT_RUNTIME_ERROR')) return
       return
    else if (present(filename)) then
       call nullify_atomic_structure(astruct)
       call set_astruct_from_file(trim(filename),bigdft_mpi%iproc,astruct)
       nat=astruct%nat
       call deallocate_atomic_structure(astruct)
       return
    end if

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
    !! max(bigdft_norb(runObj),1)
    real(gp), dimension(*), intent(out) :: eval
    !local variables
    integer :: norb

    norb=bigdft_norb(runObj)

    if(associated(runObj%rst%KSwfn%orbs%eval))then
        call f_memcpy(n=norb,src=runObj%rst%KSwfn%orbs%eval(1),dest=eval(1))
    else
       call f_err_throw('KSwfn%orbs%eval uninitialized',&
            err_name='BIGDFT_RUNTIME_ERROR')
    endif
  end subroutine bigdft_get_eval

  !BS
  function bigdft_get_units(runObj) result(units)
    implicit none
    type(run_objects), intent(in) :: runObj
    character(len=20) :: units

    units=repeat(' ',20)
    if (associated(runObj%atoms)) then
       units=runObj%atoms%astruct%units
    else
       call f_err_throw('Units uninitialized',&
            err_name='BIGDFT_RUNTIME_ERROR')
    endif
  end function bigdft_get_units
  subroutine bigdft_set_units(runObj,units)
    implicit none
    type(run_objects), intent(inout) :: runObj
    character(len=*), intent(in) :: units

    if (associated(runObj%atoms)) then
       runObj%atoms%astruct%units=units
    else
       call f_err_throw('Units uninitialized',&
            err_name='BIGDFT_RUNTIME_ERROR')
    endif

  end subroutine bigdft_set_units

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

  function bigdft_get_cell_ptr(runObj) result(cell)
    implicit none
    type(run_objects), intent(in) :: runObj
    real(gp), dimension(:), pointer :: cell

    nullify(cell)
    if (associated(runObj%atoms))then
        cell => runObj%atoms%astruct%cell_dim
    else
       call f_err_throw('Cell uninitialized',&
            err_name='BIGDFT_RUNTIME_ERROR')
    endif

  end function bigdft_get_cell_ptr

  !=====================================================================
  subroutine bigdft_state(runObj,outs,infocode)
    !IMPORTANT:
    !returns energies in hartree and
    !forces in hartree/bohr
    !(except for LJ)
    !receives distances in Bohr
    use module_lj
    use module_lenosky_si
    use public_enums
    use module_defs
    use dynamic_memory, only: f_memcpy,f_routine,f_release_routine
    use yaml_strings, only: yaml_toa, operator(+)
    use yaml_output
    use module_forces, only: clean_forces
    use module_morse_bulk
    use module_tersoff
    use module_BornMayerHugginsTosiFumi
    use module_cp2k
    use module_dftbp
    use f_enums, enum_int => int
    implicit none
    !parameters
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    integer, intent(inout) :: infocode
    !local variables
    logical :: write_mapping
    integer :: nat
    integer :: icc !for amber
    real(gp) :: maxdiff
    real(gp), dimension(3) :: alatint
    real(gp), dimension(:,:), pointer :: rxyz_ptr
    integer :: policy_tmp
!!integer :: iat , l
!!real(gp) :: anoise,tt
    call f_routine(id='bigdft_state')

    rxyz_ptr => bigdft_get_rxyz_ptr(runObj)
    nat=bigdft_nat(runObj)

    !Check the consistency between MPI processes of the atomic coordinates and broadcast them
    if (bigdft_mpi%nproc >1) then
       call mpibcast(rxyz_ptr,comm=bigdft_mpi%mpi_comm,&
            maxdiff=maxdiff)
       if (maxdiff > epsilon(1.0_gp)) then
          if (bigdft_mpi%iproc==0) then
             call yaml_warning('Input positions not identical! '//&
                  '(difference:'//trim(yaml_toa(maxdiff,fmt='(1pe12.5)'))//&
                  ' ), broadcasting from master node')
             call yaml_flush_document()
          end if
       end if
    end if

    !@NEW ####################################################
    ! Apply the constraints expressed in internal coordinates
    if (runObj%atoms%astruct%inputfile_format=='int') then
        call constraints_internal(runObj%atoms%astruct)
    end if
    !#########################################################

    call clean_state_properties(outs) !zero the state first

    !BS: is new document necessary (high overhead for FF)?
    !LG: unfortunately it it important to make testing possible. We should probably investigate 
    !    the reasons for such high overhead
    !    The new document has been substituted by sequence, not to have multiple documents for FF runs
    !    However this hybrid scheme has to be tested in the case of QM/MM runs
    !    In any case the verbosity value is used to (un)mute the output
    write_mapping= runObj%run_mode /= 'QM_RUN_MODE' .and. bigdft_mpi%iproc==0 .and. verbose > 0
    !open the document if the run_mode has not it inside
    if (write_mapping) then
       call yaml_sequence(advance='no')
       call yaml_mapping_open(trim(f_str(runObj%run_mode)),flow=.true.)
       !call yaml_new_document()
      end if
    infocode = 0
    !choose what to do by following the mode prescription
    select case(trim(f_str(runObj%run_mode)))
    case('LENNARD_JONES_RUN_MODE')
       call lenjon(nat,rxyz_ptr,outs%fxyz,outs%energy)
       !         if (bigdft_mpi%iproc == 0) then
       !            call yaml_map('LJ state, energy',outs%energy,fmt='(1pe24.17)')
       !         end if
    case('MORSE_SLAB_RUN_MODE')
        call morse_slab_wrapper(nat,bigdft_get_cell(runObj),rxyz_ptr, outs%fxyz, outs%energy)
    case('MORSE_BULK_RUN_MODE')
        call morse_bulk_wrapper(nat,bigdft_get_cell(runObj),rxyz_ptr, outs%fxyz, outs%energy)
    case('TERSOFF_RUN_MODE')
        call tersoff(nat,bigdft_get_cell(runObj),rxyz_ptr,outs%fxyz,outs%strten,outs%energy)
    case('BMHTF_RUN_MODE')
        call energyandforces_bmhtf(nat,rxyz_ptr,outs%fxyz,outs%energy)
    case('LENOSKY_SI_CLUSTERS_RUN_MODE')
       !else if(trim(adjustl(efmethod))=='LENSIc')then!for clusters
       call f_memcpy(src=rxyz_ptr,dest=runObj%mm_rst%rf_extra)
       !convert from bohr to angstroem
       call vscal(3*nat,Bohr_Ang,runObj%mm_rst%rf_extra(1,1),1)
       alatint=Bohr_Ang*bigdft_get_cell(runObj)
       call lenosky_si_shift(nat,runObj%mm_rst%rf_extra,outs%fxyz,&
            outs%energy)
       !convert energy from eV to Hartree
       outs%energy=ev_Ha*outs%energy
       !convert forces from eV/Angstroem to hartree/bohr
       call vscal(3*nat,eVAng_HaBohr,outs%fxyz(1,1),1)
    case('LENOSKY_SI_BULK_RUN_MODE')
       call f_memcpy(src=rxyz_ptr,dest=runObj%mm_rst%rf_extra)
       !convert from bohr to angstroem
       call vscal(3*nat,Bohr_Ang,runObj%mm_rst%rf_extra(1,1),1)
       alatint=Bohr_Ang*bigdft_get_cell(runObj)
       call lenosky_si(nat,alatint,runObj%mm_rst%rf_extra,outs%fxyz,outs%energy)
       !convert energy from eV to Hartree
       outs%energy=ev_Ha*outs%energy
       !convert forces from eV/Angstroem to hartree/bohr
       call vscal(3*nat,eVAng_HaBohr,outs%fxyz(1,1),1)
    case('AMBER_RUN_MODE')
       !else if(trim(adjustl(efmethod))=='AMBER')then
       icc=1
       call f_memcpy(src=rxyz_ptr,dest=runObj%mm_rst%rf_extra)
       !convert from bohr to angstroem
       call vscal(3*nat,Bohr_Ang,runObj%mm_rst%rf_extra(1,1),1)
!BS: is new document necessary (high overhead for FF)?
       !ATTENTION: call_nab_gradient returns gradient, not forces
       call call_nab_gradient(runObj%mm_rst%rf_extra,outs%fxyz,outs%energy,icc)
       outs%energy=kcalMol_Ha*outs%energy
       !convert from gradient in kcal_th/mol/angstrom to
       !force in hartree/bohr (minus before kcalMolAng_HaBohr)
       call vscal(3*nat,-kcalMolAng_HaBohr,outs%fxyz(1,1),1)
    case('QM_RUN_MODE') ! traditional BigDFT run
       !else if(trim(adjustl(efmethod))=='BIGDFT')then
       !the yaml document is created in the cluster routine
       call quantum_mechanical_state(runObj,outs,infocode)
       if (bigdft_mpi%iproc==0) call yaml_map('BigDFT infocode',infocode)
    case('CP2K_RUN_MODE') ! CP2K run mode
       call cp2k_energy_forces(nat,bigdft_get_cell(runObj),rxyz_ptr,&
            outs%fxyz,outs%energy,infocode)
       if (bigdft_mpi%iproc==0) call yaml_map('CP2K infocode',infocode)
    case('DFTBP_RUN_MODE') ! DFTB+ run mode
        call bigdft_get_input_policy(runObj, policy_tmp)
        call dftbp_energy_forces(policy_tmp,nat,bigdft_get_cell(runObj),&
             bigdft_get_astruct_ptr(runObj),bigdft_get_geocode(runObj),rxyz_ptr,&
             outs%fxyz,outs%strten,outs%energy,infocode)
       if (bigdft_mpi%iproc==0) call yaml_map('DFTB+ infocode',infocode)
    case default
       call f_err_throw('Following method for evaluation of '//&
            'energies and forces is unknown: '+ enum_int(runObj%run_mode)//&
            '('+f_str(runObj%run_mode)+')',err_name='BIGDFT_RUNTIME_ERROR')
    end select
!!         anoise=2.d-5
!!         if (anoise.ne.0.d0) then
!!         do iat=1,nat
!!         do l=1,3
!!          call random_number(tt)
!!          outs%fxyz(l,iat)=outs%fxyz(l,iat)+anoise*(tt-.5d0)
!!         enddo
!!         enddo
!!         endif
    call clean_forces(bigdft_mpi%iproc,bigdft_get_astruct_ptr(runObj),&
         rxyz_ptr,outs%fxyz,outs%fnoise,runObj%run_mode)

    !broadcast the state properties
    call broadcast_state_properties(outs)

    if (write_mapping) then
       !call yaml_release_document()
       call yaml_map('Energy',outs%energy)
       call yaml_mapping_close()
    end if

    call f_release_routine()

  end subroutine bigdft_state

  !> Routine to use BigDFT as a blackbox
  subroutine quantum_mechanical_state(runObj,outs,infocode)
    use module_base
    use yaml_output
    use module_atoms, only: astruct_dump_to_file,rxyz_inside_box
    use locregs, only: deallocate_locreg_descriptors
    use module_types, only: old_wavefunction_null
    use module_input_keys, only: set_inputpsiid
    !use communications_base
    implicit none
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    integer, intent(inout) :: infocode
    !local variables
    character(len=*), parameter :: subname='quantum_mechanical_state'
    logical :: exists
    integer :: istep,policy_tmp
    type(f_enumerator) :: inputPsiId_orig
    !integer :: iat
    external :: cluster
    !put a barrier for all the processes
    call mpibarrier(bigdft_mpi%mpi_comm)

    call f_routine(id=subname)
    !fill the rxyz array with the positions
    !wrap the atoms in the periodic directions when needed
    call rxyz_inside_box(runObj%atoms%astruct,rxyz=runObj%rst%rxyz_new)
    !assign the verbosity of the output
    !the verbose variables is defined in module_defs
    verbose=runObj%inputs%verbosity

    ! Use the restart for the linear scaling version... probably to be modified.
!!$    if(runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_WVL) then
!!$       if (runObj%rst%version == LINEAR_VERSION) then
!!$          runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_LINEAR
    if ((runObj%inputs%inputPsiId .hasattr. 'LINEAR') .and. &
         (runObj%inputs%inputPsiId .hasattr. 'MEMORY')) then
       if (any(runObj%inputs%lin%locrad_lowaccuracy /= &
            runObj%inputs%lin%locrad_highaccuracy)) then
          call f_err_throw('The radii for low and high accuracy must be the same '//&
               'when using the linear restart!',&
               err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       end if
    end if
!    end if
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

       if (runObj%inputs%inputPsiId == 'INPUT_PSI_LCAO' .and. associated(runObj%rst%KSwfn%psi)) then
          call f_free_ptr(runObj%rst%KSwfn%psi)
          call f_free_ptr(runObj%rst%KSwfn%orbs%eval)
          call deallocate_locreg_descriptors(runObj%rst%KSwfn%Lzd%Glr)
       end if

       !backdoor for hacking, finite difference method for calculating forces on particular quantities
       call f_file_exists('input.finite_difference_forces',exists)
       if (exists) then
          runObj%inputs%last_run=1 !do the last_run things nonetheless
          call set_inputpsiid(INPUT_PSI_LCAO,runObj%inputs%inputPsiId)
          !runObj%inputs%inputPsiId = INPUT_PSI_LCAO !the first run always restart from IG
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

       !the infocode should have some values defined as parameters
       select case(infocode)
       case(1) !wfn not yet converged
          call bigdft_get_input_policy(runObj,policy_tmp)
          if (policy_tmp == INPUT_POLICY_SCRATCH) &
               call bigdft_set_input_policy(INPUT_POLICY_MEMORY,runObj)
          !do not know why we do this, meaningful only if we did not exit
          exit loop_cluster
       case(2) !too big residue after restart, come back to scratch
          call bigdft_set_input_policy(INPUT_POLICY_SCRATCH,runObj)
       case(3) !even with scracth policy, residue is too big. Exit with failure
          !for the moment there is no crash for the linear scaling
          if (runObj%inputs%inputPsiId .hasattr. 'CUBIC') then
             if (bigdft_mpi%iproc == 0) then
                call astruct_dump_to_file(runObj%atoms%astruct,"posfail",&
                     'UNCONVERGED WF ',outs%energy,runObj%rst%rxyz_new)
             end if

             call f_free_ptr(runObj%rst%KSwfn%psi)
             call f_free_ptr(runObj%rst%KSwfn%orbs%eval)
             call deallocate_locreg_descriptors(runObj%rst%KSwfn%Lzd%Glr)

             !test if stderr works
             write(0,*) 'unnormal end'
             call mpibarrier(bigdft_mpi%mpi_comm)
             call f_err_throw('Convergence error (probably gnrm>4.0), cannot proceed. '//&
                  'Writing positions in file posfail.xyz',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          call set_inputpsiid(INPUT_PSI_RANDOM,runObj%inputs%inputPsiId)
       case default
          exit loop_cluster
       end select

!!$       !Check infocode in function of the inputPsiId parameters
!!$       !and change the strategy of input guess psi
!!$       if (runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_WVL .and. infocode==2) then
!!$          if (runObj%inputs%gaussian_help) then
!!$             runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_GAUSS
!!$          else
!!$             runObj%inputs%inputPsiId = INPUT_PSI_LCAO
!!$          end if
!!$       else if (runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_LINEAR .and. infocode==2) then
!!$          ! problems after restart for linear version
!!$          runObj%inputs%inputPsiId = INPUT_PSI_LINEAR_AO
!!$       else if ((runObj%inputs%inputPsiId == INPUT_PSI_MEMORY_WVL .or. &
!!$            runObj%inputs%inputPsiId == INPUT_PSI_LCAO) .and. infocode==1) then
!!$          runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_WVL
!!$          exit loop_cluster
!!$       else if (runObj%inputs%inputPsiId == INPUT_PSI_LCAO .and. infocode==3) then
!!$          if (bigdft_mpi%iproc == 0) then
!!$             call astruct_dump_to_file(runObj%atoms%astruct,"posfail",&
!!$                  'UNCONVERGED WF ',outs%energy,runObj%rst%rxyz_new)
!!$          end if
!!$
!!$          call f_free_ptr(runObj%rst%KSwfn%psi)
!!$          call f_free_ptr(runObj%rst%KSwfn%orbs%eval)
!!$          call deallocate_locreg_descriptors(runObj%rst%KSwfn%Lzd%Glr)
!!$
!!$          !test if stderr works
!!$          write(0,*) 'unnormal end'
!!$          call mpibarrier(bigdft_mpi%mpi_comm)
!!$          call f_err_throw('Convergence error (probably gnrm>4.0), cannot proceed. '//&
!!$               'Writing positions in file posfail.xyz',err_name='BIGDFT_RUNTIME_ERROR')
!!$       else
!!$          exit loop_cluster
!!$       end if

    end do loop_cluster

    !preserve the previous value
    runObj%inputs%inputPsiId=inputPsiId_orig

    !put a barrier for all the processes
    call f_release_routine()
    call mpibarrier(bigdft_mpi%mpi_comm)

  END SUBROUTINE quantum_mechanical_state


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
    use yaml_strings, only: yaml_toa
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

  !> Calculate atomic forces via finite differences (test purpose)
  subroutine forces_via_finite_differences(iproc,nproc,atoms,inputs,energy,fxyz,fnoise,rst,infocode)
    use module_base
    use module_types
    use module_atoms, only: move_this_coordinate
    use module_forces
    use module_input_keys, only: inputpsiid_set_policy
    implicit none
    integer, intent(in) :: iproc,nproc
    integer, intent(inout) :: infocode
    real(gp), intent(inout) :: energy,fnoise
    type(input_variables), intent(inout) :: inputs
    type(atoms_data), intent(inout) :: atoms
    type(QM_restart_objects), intent(inout) :: rst
    real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: fxyz
    !local variables
    character(len=*), parameter :: subname='forces_via_finite_differences'
    character(len=4) :: cc
    integer :: ik,km,n_order,iat,ii,i,k,order,iorb_ref
    real(gp) :: dd,alat,functional_ref,fd_alpha,energy_ref,pressure
    real(gp), dimension(3) :: fd_step
    real(gp), dimension(6) :: strten
    integer, dimension(:), allocatable :: kmoves
    real(gp), dimension(:), allocatable :: functional,dfunctional
    real(gp), dimension(:,:), allocatable :: rxyz_ref, fxyz_fake
    type(energy_terms) :: energs

    if (iproc == 0) then
       write(*,*)
       write(*,'(1x,a,59("="))') '=Forces via finite Difference '
    end if

    !read the file (experimental version)
    open(unit=79,file='input.finite_difference_forces',status='unknown')
    read(79,*) order,fd_alpha
    read(79,*) iorb_ref
    close(unit=79)

    !read the step size
    ! Initialize freq_step (step to move atoms)
    fd_step(1) = fd_alpha*inputs%hx
    fd_step(2) = fd_alpha*inputs%hy
    fd_step(3) = fd_alpha*inputs%hz

    !first, mark the reference energy
    energy_ref=energy

    !assign the reference
    functional_ref=functional_definition(iorb_ref,energy)

    if (order == -1) then
       n_order = 1
       kmoves = f_malloc(src=(/ -1 /),id='kmoves')
    else if (order == 1) then
       n_order = 1
       kmoves = f_malloc(src=(/ 1 /),id='kmoves')
    else if (order == 2) then
       n_order = 2
       kmoves = f_malloc(src=(/ -1, 1 /),id='kmoves')
    else if (order == 3) then
       n_order = 4
       kmoves = f_malloc(src=(/ -2, -1, 1, 2 /),id='kmoves')
    else
       print *, "Finite Differences: This order",order," is not implemented!"
       stop
    end if

    functional = f_malloc(n_order,id='functional')
    dfunctional = f_malloc0(3*atoms%astruct%nat,id='dfunctional')
    rxyz_ref = f_malloc(src=rst%rxyz_new,id='rxyz_ref')
    fxyz_fake = f_malloc((/ 3, atoms%astruct%nat /),id='fxyz_fake')

    do iat=1,atoms%astruct%nat

       do i=1,3 !a step in each of the three directions

          if (.not.move_this_coordinate(atoms%astruct%ifrztyp(iat),i)) then
             if (iproc == 0) write(*,"(1x,a,i0,a,i0,a)") '=F:The direction ',i,' of the atom ',iat,' is frozen.'
             cycle
          end if

          ii = i+3*(iat-1)
          if (i==1) then
             alat=atoms%astruct%cell_dim(1)
             cc(3:4)='*x'
          else if (i==2) then
             alat=atoms%astruct%cell_dim(2)
             cc(3:4)='*y'
          else
             alat=atoms%astruct%cell_dim(3)
             cc(3:4)='*z'
          end if
          km = 0
          functional=0.0_gp
          do ik=1,n_order
             k = kmoves(ik)
             !-1-> 1, 1 -> 2, y = ( x + 3 ) / 2
             km = km + 1
             write(cc(1:2),"(i2)") k
             !Displacement
             dd=real(k,gp)*fd_step(i)
             !We copy atomic positions (necessary?)
             !call vcopy(3*atoms%astruct%nat,rxyz_ref(1,1),1,rst%rxyz_new(1,1),1)
             call f_memcpy(src=rxyz_ref,dest=rst%rxyz_new)
             if (iproc == 0) then
                write(*,"(1x,a,i0,a,a,a,1pe20.10,a)") &
                     '=FD Move the atom ',iat,' in the direction ',cc,' by ',dd,' bohr'
             end if
             if (atoms%astruct%geocode == 'P') then
                rst%rxyz_new(i,iat)=modulo(rxyz_ref(i,iat)+dd,alat)
             else if (atoms%astruct%geocode == 'S') then
                rst%rxyz_new(i,iat)=modulo(rxyz_ref(i,iat)+dd,alat)
             else
                rst%rxyz_new(i,iat)=rxyz_ref(i,iat)+dd
             end if
             !inputs%inputPsiId=1
             call inputpsiid_set_policy(ENUM_MEMORY,inputs%inputPsiId)
             !here we should call cluster
             call cluster(nproc,iproc,atoms,rst%rxyz_new,energy,energs,fxyz_fake,strten,fnoise,pressure,&
                  rst%KSwfn,rst%tmb,&!psi,rst%Lzd,rst%gaucoeffs,rst%gbd,rst%orbs,&
                  rst%rxyz_old,inputs,rst%GPU,infocode)

             !assign the quantity which should be differentiated
             functional(km)=functional_definition(iorb_ref,energy)

          end do
          ! Build the finite-difference quantity if the calculation has converged properly
          if (infocode ==0) then
             !Force is -dE/dR
             if (order == -1) then
                dd = - (functional_ref - functional(1))/fd_step(i)
             else if (order == 1) then
                dd = - (functional(1) - functional_ref)/fd_step(i)
             else if (order == 2) then
                dd = - (functional(2) - functional(1))/(2.0_gp*fd_step(i))
             else if (order == 3) then
                dd = - (functional(4) + functional(3) - functional(2) - functional(1))/(6.d0*fd_step(i))
             else
                stop "BUG (FD_forces): this order is not defined"
             end if
             !if (abs(dd).gt.1.d-10) then
             dfunctional(ii) = dd
             !end if
          else
             if (iproc==0)&
                  write(*,*)&
                  'ERROR: the wavefunctions have not converged properly, meaningless result. Exiting. Infocode:',&
                  infocode
             stop
          end if

       end do
    end do

    !copy the final value of the energy and of the dfunctional
    if (.not. experimental_modulebase_var_onlyfion) then !normal case
       call vcopy(3*atoms%astruct%nat,dfunctional(1),1,fxyz(1,1),1)
    else
       call axpy(3*atoms%astruct%nat,2.0_gp*rst%KSwfn%orbs%norb,dfunctional(1),1,fxyz(1,1),1)
    end if
    !clean the center mass shift and the torque in isolated directions
    call clean_forces(iproc,atoms%astruct,rxyz_ref,fxyz,fnoise)
    if (iproc == 0) call write_forces(atoms%astruct,fxyz)

    energy=functional_ref

    if (iproc == 0) then
       write(*,"(1x,2(a,1pe20.10))") &
            '=FD Step done, Internal Energy:',energy_ref,' functional value:', functional_ref
    end if

    call f_free(kmoves)
    call f_free(functional)
    call f_free(dfunctional)
    call f_free(rxyz_ref)
    call f_free(fxyz_fake)

  contains

    function functional_definition(iorb_ref,energy)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iorb_ref
      real(gp), intent(in) :: energy
      real(gp) :: functional_definition
      !local variables
      real(gp) :: mu

      !chemical potential =1/2(e_HOMO+e_LUMO)= e_HOMO + 1/2 GAP (the sign is to be decided - electronegativity?)
      !definition which brings to Chemical Potential
      if (rst%KSwfn%orbs%HLgap/=UNINITIALIZED(rst%KSwfn%orbs%HLgap) .and. iorb_ref< -1) then
         mu=-abs(rst%KSwfn%orbs%eval(-iorb_ref)+ 0.5_gp*rst%KSwfn%orbs%HLgap)
      else
         mu=UNINITIALIZED(1.0_gp)
      end if

      !assign the reference
      if (iorb_ref==0) then
         functional_definition=energy
      else if (iorb_ref == -1) then
         if (rst%KSwfn%orbs%HLgap/=UNINITIALIZED(rst%KSwfn%orbs%HLgap)) then
            functional_definition=rst%KSwfn%orbs%HLgap !here we should add the definition which brings to Fukui function
         else
            stop ' ERROR (FDforces): gap not defined'
         end if
      else if(iorb_ref < -1) then      !definition which brings to the neutral fukui function (chemical potential)
         if (rst%KSwfn%orbs%HLgap/=UNINITIALIZED(rst%KSwfn%orbs%HLgap)) then
            functional_definition=mu!-mu*real(2*orbs%norb,gp)+energy
         else
            stop ' ERROR (FDforces): gap not defined, chemical potential cannot be calculated'
         end if
      else
         functional_definition=rst%KSwfn%orbs%eval(iorb_ref)
      end if

    end function functional_definition

  end subroutine forces_via_finite_differences


  !> Get the current run policy
  subroutine bigdft_get_input_policy(runObj, policy)
    use module_base
    use module_input_keys, only: inputpsiid_get_policy
    use public_enums
    implicit none
    ! Calling arguments
    type(run_objects),intent(in) :: runObj
    integer,intent(out) :: policy
    !local variables
    type(f_enumerator) :: policy_enum

    call inputpsiid_get_policy(runObj%inputs%inputPsiId,policy_enum)

    select case(trim(f_char(policy_enum)))
    case('SCRATCH')
       policy = INPUT_POLICY_SCRATCH
    case('MEMORY')
       policy = INPUT_POLICY_MEMORY
    case('FILE')
       policy = INPUT_POLICY_DISK
    case default
       call f_err_throw('The specified value of inputPsiId ('//&
            trim(f_char(runObj%inputs%inputPsiId))//') is not valid for policy '//&
            trim(f_char((policy_enum))), &
            err_name='BIGDFT_RUNTIME_ERROR')
    end select

!!$    select case(runObj%inputs%inputPsiId)
!!$    case(INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_LCAO, INPUT_PSI_LINEAR_AO)
!!$        ! Start the calculation from scratch
!!$        policy = INPUT_POLICY_SCRATCH
!!$    case(INPUT_PSI_MEMORY_WVL, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_MEMORY_LINEAR)
!!$        ! Start the calculation from data in memory
!!$        policy = INPUT_POLICY_MEMORY
!!$    case(INPUT_PSI_CP2K, INPUT_PSI_DISK_WVL, INPUT_PSI_DISK_GAUSS, INPUT_PSI_DISK_LINEAR)
!!$        ! Start the calculation from data on disk
!!$        policy = INPUT_POLICY_DISK
!!$    case default
!!$        call f_err_throw('The specified value of inputPsiId ('//&
!!$            &trim(yaml_toa(runObj%inputs%inputPsiId))//') is not valid', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$    end select
  end subroutine bigdft_get_input_policy


  !> Set the current run policy
  !! modify the values of inputpsiid according to
  !! the previously chosen policy
  subroutine bigdft_set_input_policy(policy, runObj)
    use module_base
    use module_input_keys, only: inputpsiid_set_policy
    use public_enums
    implicit none
    ! Calling arguments
    integer,intent(in) :: policy
    type(run_objects),intent(inout) :: runObj
    ! Local variables
    integer,parameter :: MODE_PLAIN    = 0
    integer,parameter :: MODE_CUBIC    = 201
    integer,parameter :: MODE_LINEAR   = 202
    integer,parameter :: MODE_GAUSSIAN = 203
    integer,parameter :: MODE_CP2K     = 204
    integer,parameter :: MODE_EMPTY    = 205
    integer,parameter :: MODE_RANDOM   = 206
    integer :: mode
    type(f_enumerator) :: policy_enum


    ! Set the new ID for the input guess
    select case(policy)
    case(INPUT_POLICY_SCRATCH)
       policy_enum=ENUM_SCRATCH
    case(INPUT_POLICY_MEMORY)
       policy_enum=ENUM_MEMORY
    case(INPUT_POLICY_DISK)
       policy_enum=ENUM_FILE
    end select

    call inputpsiid_set_policy(policy_enum,runObj%inputs%inputPsiId)
    
!!$
!!$    ! Check which type of run this is, based on the current value of inputPsiId
!!$    select case (runObj%inputs%inputPsiId)
!!$    case(INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, INPUT_PSI_DISK_WVL)
!!$        ! This is like a cubic run
!!$        mode = MODE_CUBIC
!!$    case(INPUT_PSI_LINEAR_AO, INPUT_PSI_MEMORY_LINEAR, INPUT_PSI_DISK_LINEAR)
!!$        ! This is a linear run
!!$        mode = MODE_LINEAR
!!$    case(INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS)
!!$        ! This is a run based on a Gaussian input guess
!!$        mode = MODE_GAUSSIAN
!!$    case(INPUT_PSI_CP2K)
!!$        ! This is a run based on an input guess coming from CP2K
!!$        mode = MODE_CP2K
!!$    case(INPUT_PSI_EMPTY)
!!$        ! This is a run based on an empty input guess
!!$        mode = MODE_EMPTY
!!$    case(INPUT_PSI_RANDOM)
!!$        ! This is a run based on a random input guess
!!$        mode = MODE_RANDOM
!!$    case default
!!$        call f_err_throw('The specified value of inputPsiId ('//&
!!$            &trim(yaml_toa(runObj%inputs%inputPsiId))//') is not valid', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$    end select
!!$
!!$    ! Set the new ID for the input guess
!!$    select case(policy)
!!$    case(INPUT_POLICY_SCRATCH)
!!$        ! Start the calculation from scratch
!!$        select case(mode)
!!$        case(MODE_CUBIC)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_LCAO
!!$        case(MODE_LINEAR)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_LINEAR_AO
!!$        case(MODE_GAUSSIAN)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_LCAO_GAUSS
!!$        case(MODE_CP2K)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_CP2K
!!$        case(MODE_EMPTY)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_EMPTY
!!$        case(MODE_RANDOM)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_RANDOM
!!$        case default
!!$        call f_err_throw('The specified value of mode ('//&
!!$            &trim(yaml_toa(mode))//') is not valid', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$        end select
!!$
!!$    case(INPUT_POLICY_MEMORY)
!!$        ! Start the calculation from data in memory
!!$        select case(mode)
!!$        case(MODE_CUBIC)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_WVL
!!$        case(MODE_LINEAR)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_LINEAR
!!$        case(MODE_GAUSSIAN)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_MEMORY_GAUSS
!!$        case(MODE_CP2K, MODE_EMPTY, MODE_RANDOM)
!!$        call f_err_throw('The specified value of mode ('//&
!!$            &trim(yaml_toa(mode))//') is not compatible with the input policy ('//&
!!$            &trim(yaml_toa(policy))//')', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$        case default
!!$        call f_err_throw('The specified value of mode ('//&
!!$            &trim(yaml_toa(mode))//') is not valid', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$        end select
!!$
!!$    case(INPUT_POLICY_DISK)
!!$        select case(mode)
!!$        case(MODE_CUBIC)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_DISK_WVL
!!$        case(MODE_LINEAR)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_DISK_LINEAR
!!$        case(MODE_GAUSSIAN)
!!$            runObj%inputs%inputPsiId = INPUT_PSI_DISK_GAUSS
!!$        case(MODE_CP2K, MODE_EMPTY, MODE_RANDOM)
!!$        call f_err_throw('The specified value of mode ('//&
!!$            &trim(yaml_toa(mode))//') is not compatible with the input policy ('//&
!!$            &trim(yaml_toa(policy))//')', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$        case default
!!$        call f_err_throw('The specified value of mode ('//&
!!$            &trim(yaml_toa(mode))//') is not valid', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$        end select
!!$
!!$    case default
!!$        call f_err_throw('The specified value of inputPsiId ('//&
!!$            &yaml_toa(runObj%inputs%inputPsiId)//') is not valid', &
!!$            err_name='BIGDFT_RUNTIME_ERROR')
!!$    end select
  end subroutine bigdft_set_input_policy


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
  !run_dict => dict_new('name' .is. radical, 'posinp' .is. posinp)

  !call bigdft_run_new(run_dict)
  nullify(run_dict)
  call bigdft_set_run_properties(run_dict,run_id=radical,posinp_id=posinp)

  call run_objects_init(runObj,run_dict)
  !call dict_free(run_dict) In this case the run_dict has not to be freed
END SUBROUTINE run_objects_init_from_run_name

subroutine run_objects_update(runObj, dict)
  use module_base, only: bigdft_mpi,f_int
  use bigdft_run, only: run_objects,init_QM_restart_objects,init_MM_restart_objects,set_run_objects,bigdft_nat
  use dictionaries!, only: dictionary, dict_update,dict_copy,dict_free,dict_iter,dict_next
  use yaml_output
  use module_input_dicts, only: create_log_file
  implicit none
  type(run_objects), intent(inout) :: runObj
  type(dictionary), pointer :: dict
  !local variables
  type(dictionary), pointer :: item
  logical :: dict_from_files

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

  call create_log_file(runObj%user_inputs,dict_from_files)

  ! Parse new dictionary.
  call set_run_objects(runObj)

  !init and update the restart objects
  call init_QM_restart_objects(bigdft_mpi%iproc,runObj%inputs,runObj%atoms,&
       runObj%rst)
  call init_MM_restart_objects(runObj,runObj%mm_rst,bigdft_nat(runObj),runObj%run_mode)
END SUBROUTINE run_objects_update

!> this routine should be used in memguess executable also
subroutine run_objects_system_setup(runObj, iproc, nproc, rxyz, shift, mem)
  use module_base, only: gp,f_memcpy,f_enumerator
  use bigdft_run
  use module_types
  use module_fragments
  use module_interfaces, only: system_initialization
  use psp_projectors_base, only: free_DFT_PSP_projectors
  use communications_base, only: deallocate_comms
  implicit none
  type(run_objects), intent(inout) :: runObj
  integer, intent(in) :: iproc, nproc
  real(gp), dimension(3,runObj%atoms%astruct%nat), intent(out) :: rxyz
  real(gp), dimension(3), intent(out) :: shift
  type(memory_estimation), intent(out) :: mem

  integer :: input_wf_format
  type(DFT_PSP_projectors) :: nlpsp
  type(f_enumerator) :: inputpsi
  type(system_fragment), dimension(:), pointer :: ref_frags
  character(len = *), parameter :: subname = "run_objects_estimate_memory"

  ! Copy rxyz since system_size() will shift them.
!!$  allocate(rxyz(3,runObj%atoms%astruct%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,rxyz,'rxyz',subname)
  call f_memcpy(src=runObj%atoms%astruct%rxyz,dest=rxyz)
  !call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, rxyz(1,1), 1)
  inputpsi=runObj%inputs%inputpsiid
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
