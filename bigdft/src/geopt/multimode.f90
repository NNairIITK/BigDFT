module bigdft_run_sections
  use bigdft_run, only: state_properties
  use module_atoms, only: atomic_structure
  implicit none

  type sub_state_properties
     type(state_properties) :: sub
     type(state_properties), pointer :: ref
     integer, dimension(:), pointer :: map
  end type sub_state_properties

contains
  subroutine sub_state_properties_from_asub(subouts, asub, outs)
    use bigdft_run, only: init_state_properties
    use public_keys, only: ASTRUCT_ATT_ORIG_ID
    use dictionaries
    use dynamic_memory
    use f_utils
    implicit none
    type(sub_state_properties), intent(out) :: subouts
    type(atomic_structure), intent(in) :: asub
    type(state_properties), target, intent(in) :: outs

    integer :: iat

    call init_state_properties(subouts%sub, asub%nat)
    subouts%map = f_malloc_ptr((/ asub%nat /), id = "map")
    call f_zero(subouts%map)
    do iat = 1, asub%nat
       if (ASTRUCT_ATT_ORIG_ID .in. asub%attributes(iat)%d) then
          subouts%map(iat) = asub%attributes(iat)%d // ASTRUCT_ATT_ORIG_ID
       end if
    end do
    subouts%ref => outs
  end subroutine sub_state_properties_from_asub

  subroutine deallocate_sub_state_properties(subouts)
    use bigdft_run, only: deallocate_state_properties
    use dynamic_memory
    implicit none
    type(sub_state_properties), intent(inout) :: subouts

    call deallocate_state_properties(subouts%sub)
    call f_free_ptr(subouts%map)
  end subroutine deallocate_sub_state_properties

end module bigdft_run_sections

subroutine multi_mode_extract(asub, runObj, section, passivation, buf, last)
  use module_atoms, only: atomic_structure, astruct_at_from_dict
  use bigdft_run
  use dictionaries
  use dynamic_memory
  implicit none
  ! Parameters
  type(atomic_structure), intent(out) :: asub
  type(run_objects), intent(in) :: runObj
  character(len = *), intent(in) :: section
  logical, intent(in) :: passivation, last
  integer, intent(in) :: buf
  ! Local variables
  integer :: iat
  logical, dimension(:), allocatable :: mask
  character(len = max_field_length) :: mode

  ! Generate the mask from the MODE atomic attribute.
  mask = f_malloc(runObj%atoms%astruct%nat, id = "mask")
  do iat = 1, runObj%atoms%astruct%nat
     call astruct_at_from_dict(runObj%atoms%astruct%attributes(iat)%d, mode = mode)
     mask(iat) = (trim(mode) == section) .or. (last .and. len_trim(mode) == 0)
  end do
  call astruct_from_subset(asub, runObj%atoms%astruct, runObj%atoms%astruct%rxyz, &
       & mask, passivation, buf, "yes")
  call f_free(mask)
END SUBROUTINE multi_mode_extract

subroutine multi_mode_state(runObj, outs, infocode)
  use module_atoms, only: atomic_structure, deallocate_atomic_structure, astruct_dump_to_file
  use bigdft_run
  use bigdft_run_sections
  use dynamic_memory
  use public_keys, only: ASTRUCT_ATT_ORIG_ID
  use module_defs, only: gp
  use dictionaries
  use module_f_objects, only: f_object_signal_prepare, f_object_has_signal, f_object_signal_emit
  implicit none
  ! Parameters
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  integer, intent(inout) :: infocode
  ! Local variables
  type(sub_state_properties), dimension(:), allocatable :: subouts
  integer :: ln, i
  type(atomic_structure) :: asub

  ln = size(runObj%sections)
  allocate(subouts(ln))
  infocode = 0

  ! Run subparts and store forces.
  do i = 1, ln
     !> @todo Add a signal to enable relabeling here.

     ! Need to fully re-extract here.
     call multi_mode_extract(asub, runObj, trim(runObj%sections(i)%label), &
          & runObj%inputs%multi_pass(i), runObj%inputs%multi_buf(i), (i == ln))
     !> @todo Handle the case where the number of atoms in the section
     !! vary because of movements.
     call bigdft_set_rxyz(runObj%sections(i), rxyz = asub%rxyz)

     call sub_state_properties_from_asub(subouts(i), asub, outs)
     call deallocate_atomic_structure(asub)
     call process_run(trim(runObj%sections(i)%label), runObj%sections(i), subouts(i)%sub)
  end do

  ! Mix state_properties, either default or by signal.
!!$  if (f_object_signal_prepare("run_objects", "mix")) then
!!$     call f_object_signal_add_arg("run_objects", "mix", runObj)
!!$     call f_object_signal_add_arg("run_objects", "mix", outs)
!!$     call f_object_signal_add_arg("run_objects", "mix", subouts)
!!$     call f_object_signal_emit("run_objects", "mix")
!!$  else
  call union_mix_subouts(runObj, outs, subouts)
!!$  end if

  do i = 1, ln
     call deallocate_sub_state_properties(subouts(i))
  end do
  deallocate(subouts)
END SUBROUTINE multi_mode_state

subroutine union_mix_subouts(runObj, outs, subouts)
  use bigdft_run
  use bigdft_run_sections
  use f_utils
  implicit none
  type(run_objects), intent(in) :: runObj
  type(state_properties), intent(inout) :: outs
  type(sub_state_properties), dimension(size(runObj%sections)), intent(in) :: subouts

  integer :: i, iat, jat

  outs%energy = 0.
  call f_zero(outs%fxyz)

  ! Update position, energy and forces.
  do i = 1, size(runObj%sections)
     ! Update positions and forces (positions may have been altered by the run).
     do iat = 1, size(subouts(i)%map)
        jat = subouts(i)%map(iat)
        if (jat > 0) then
           runObj%atoms%astruct%rxyz(:, jat) = runObj%sections(i)%atoms%astruct%rxyz(:, iat)
           outs%fxyz(:, jat) = outs%fxyz(:, jat) + subouts(i)%sub%fxyz(:, iat)
        end if
     end do
     outs%energy = outs%energy + subouts(i)%sub%energy
  end do
end subroutine union_mix_subouts
