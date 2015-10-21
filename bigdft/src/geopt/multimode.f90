subroutine multi_mode_state(runObj, outs, infocode)
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

  ln = size(runObj%sections)
  infocode = 0

  ! Run subparts and accumulate forces.
  do i = 1, ln
     nat = bigdft_nat(runObj%sections(i))
     call init_state_properties(subouts, nat)
     call bigdft_state(runObj%sections(i), subouts, icode)
     infocode = max(infocode, icode)

     map = f_malloc((/ nat /), id = "maps")
     map = 1
     coeffs = f_malloc0((/ nat /), id = "coeffs")
     do iat = 1, nat
        if (ASTRUCT_ATT_ORIG_ID .in. runObj%sections(i)%atoms%astruct%attributes(iat)%d) then
           map(iat) = runObj%sections(i)%atoms%astruct%attributes(iat)%d // ASTRUCT_ATT_ORIG_ID
           ! Simple model, coefficients are unity.
           coeffs(iat) = 1._gp
        end if
     end do

     ! Mix the outs.
     call multi_fxyz_axpy(coeffs, subouts, outs, map)

     ! Update the positions that may have been altered by the run.
     do iat = 1, nat
        runObj%atoms%astruct%rxyz(:, map(iat)) = runObj%sections(i)%atoms%astruct%rxyz(:, iat)
     end do

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
     outy%fxyz(:, map(idim)) = alpha(idim) * outx%fxyz(:, idim) + outy%fxyz(:, map(idim))
  end do
END SUBROUTINE multi_fxyz_axpy
