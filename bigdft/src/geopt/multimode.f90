subroutine multi_mode_state(runObj, outs, infocode)
  use module_atoms, only: atomic_structure, deallocate_atomic_structure, astruct_dump_to_file
  use bigdft_run
  use dynamic_memory
  use module_defs, only: gp
  use dictionaries
  use module_f_objects
  implicit none
  ! Parameters
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  integer, intent(inout) :: infocode
  ! Local variables
  type(state_properties), dimension(:), allocatable :: subouts
  integer :: ln, i
  type(atomic_structure) :: asub
  type(signal_ctx) :: sig

  ln = size(runObj%sections)
  allocate(subouts(ln))
  infocode = 0

  ! Run subparts and store forces.
  do i = 1, ln
     !> @todo Add a signal to enable relabeling here.

     ! Need to fully re-extract here.
     call section_extract(asub, runObj%sections(i)%astruct_map, runObj, trim(runObj%sections(i)%label), &
          & runObj%inputs%multi_pass(i), runObj%inputs%multi_buf(i), (i == ln))
     !> @todo Handle the case where the number of atoms in the section
     !! vary because of movements.
     call bigdft_set_rxyz(runObj%sections(i), rxyz = asub%rxyz)

     call init_state_properties(subouts(i), asub%nat)
     call deallocate_atomic_structure(asub)
     call process_run(trim(runObj%sections(i)%label), runObj%sections(i), subouts(i))
  end do

  ! Mix state_properties, either default or by signal.
  if (f_object_signal_prepare("run_objects", "join", sig)) then
     call f_object_signal_add_arg(sig, runObj)
     call f_object_signal_add_arg(sig, outs)
     call f_object_signal_add_arg(sig, subouts)
     call f_object_signal_emit(sig)
  else
     call union_mix_subouts(runObj, outs, subouts)
  end if

  do i = 1, ln
     call deallocate_state_properties(subouts(i))
  end do
  deallocate(subouts)
END SUBROUTINE multi_mode_state

subroutine union_mix_subouts(runObj, outs, subouts)
  use bigdft_run
  use f_utils
  implicit none
  type(run_objects), intent(in) :: runObj
  type(state_properties), intent(inout) :: outs
  type(state_properties), dimension(size(runObj%sections)), intent(in) :: subouts

  integer :: i, iat, jat

  outs%energy = 0.
  call f_zero(outs%fxyz)

  ! Update position, energy and forces.
  do i = 1, size(runObj%sections)
     ! Update positions and forces (positions may have been altered by the run).
     do iat = 1, size(runObj%sections(i)%astruct_map)
        jat = runObj%sections(i)%astruct_map(iat)
        if (jat > 0) then
           runObj%atoms%astruct%rxyz(:, jat) = runObj%sections(i)%atoms%astruct%rxyz(:, iat)
           outs%fxyz(:, jat) = outs%fxyz(:, jat) + subouts(i)%fxyz(:, iat)
        end if
     end do
     outs%energy = outs%energy + subouts(i)%energy
  end do
end subroutine union_mix_subouts
