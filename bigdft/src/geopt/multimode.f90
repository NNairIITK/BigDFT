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
  use dynamic_memory
  use public_keys, only: ASTRUCT_ATT_ORIG_ID
  use module_defs, only: gp
  use dictionaries
  implicit none
  ! Parameters
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  integer, intent(inout) :: infocode
  ! Local variables
  type(state_properties) :: subouts
  integer :: ln, i, iat, nat, icode
  integer, dimension(:), allocatable :: map
  real(gp), dimension(:), allocatable :: coeffs
  type(atomic_structure) :: asub

  ln = size(runObj%sections)
  infocode = 0

  ! Run subparts and accumulate forces.
  do i = 1, ln
     ! Need to fully re-extract here.
     call multi_mode_extract(asub, runObj, trim(runObj%sections(i)%label), &
          & runObj%inputs%multi_pass(i), runObj%inputs%multi_buf(i), (i == ln))
     !@todo Handle the case where the number of atoms in the section
     !      vary because of movements.
     call bigdft_set_rxyz(runObj%sections(i), rxyz = asub%rxyz)
     call deallocate_atomic_structure(asub)

     nat = bigdft_nat(runObj%sections(i))
     call init_state_properties(subouts, nat)
     call bigdft_state(runObj%sections(i), subouts, icode)
     infocode = max(infocode, icode)

     map = f_malloc0((/ nat /), id = "maps")
     coeffs = f_malloc0((/ nat /), id = "coeffs")
     do iat = 1, nat
        if (ASTRUCT_ATT_ORIG_ID .in. runObj%sections(i)%atoms%astruct%attributes(iat)%d) then
           map(iat) = runObj%sections(i)%atoms%astruct%attributes(iat)%d // ASTRUCT_ATT_ORIG_ID
           !@todo Simple force mixing model, coefficients are unity.
           coeffs(iat) = 1._gp
        end if
     end do

     ! Mix the outs.
     call multi_fxyz_axpy(coeffs, subouts, outs, map)

     ! Update the positions that may have been altered by the run.
     do iat = 1, nat
        if (map(iat) > 0) then
           runObj%atoms%astruct%rxyz(:, map(iat)) = runObj%sections(i)%atoms%astruct%rxyz(:, iat)
        end if
     end do

     !@todo Hardcode the total energy to the first section only. To be improved.
     outs%energy = subouts%energy

     call f_free(coeffs)
     call f_free(map)

     call deallocate_state_properties(subouts)
  end do
END SUBROUTINE multi_mode_state

subroutine multi_fxyz_axpy(alpha, outx, outy, map)
  use bigdft_run
  use module_defs, only: gp
  implicit none
  type(state_properties), intent(inout) :: outy
  type(state_properties), intent(in) :: outx
  real(gp), dimension(outx%fdim), intent(in) :: alpha
  integer, dimension(outx%fdim), intent(in) :: map

  integer :: idim

  do idim = 1, outx%fdim
     if (map(idim) > 0) then
        outy%fxyz(:, map(idim)) = alpha(idim) * outx%fxyz(:, idim) + outy%fxyz(:, map(idim))
     end if
  end do
END SUBROUTINE multi_fxyz_axpy
