subroutine init()
  use yaml_output
  use module_base, only: bigdft_mpi
  use module_f_objects
  integer :: sid
  type(kernel_ctx) :: kernel

  interface
     subroutine smooth_mix(runObj, outs, subs)
       use bigdft_run
       type(run_objects), intent(in) :: runObj
       type(state_properties), intent(inout) :: outs
       type(state_properties), dimension(size(runObj%sections)), intent(in) :: subs
     end subroutine smooth_mix
  end interface

  if (bigdft_mpi%iproc == 0) call yaml_map("mix plugin", .true.)

  call f_object_kernel_new(kernel, smooth_mix, 3)
  call f_object_signal_connect("run_objects", "join", kernel, sid)
end subroutine init

subroutine smooth_mix(runObj, outs, subs)
  use module_defs, only: gp
  use bigdft_run
  use f_utils
  use f_enums
  use public_enums, only: QM_RUN_MODE
  implicit none
  type(run_objects), intent(in) :: runObj
  type(state_properties), intent(inout) :: outs
  type(state_properties), dimension(size(runObj%sections)), intent(in) :: subs
  
  integer :: i, iat, jat
  real(gp) :: factor

  outs%energy = 0.
  call f_zero(outs%fxyz)

  ! Update position, energy and forces.
  do i = 1, size(runObj%sections)
     factor = 1._gp
     if (runObj%sections(i)%run_mode /= QM_RUN_MODE) factor = 0.5_gp
     ! Update positions and forces (positions may have been altered by the run).
     do iat = 1, size(runObj%sections(i)%astruct_map)
        jat = runObj%sections(i)%astruct_map(iat)
        if (jat > 0) then
           runObj%atoms%astruct%rxyz(:, jat) = runObj%sections(i)%atoms%astruct%rxyz(:, iat)
           outs%fxyz(:, jat) = outs%fxyz(:, jat) + factor * subs(i)%fxyz(:, iat)
        end if
     end do
     outs%energy = outs%energy + subs(i)%energy
  end do
end subroutine smooth_mix
