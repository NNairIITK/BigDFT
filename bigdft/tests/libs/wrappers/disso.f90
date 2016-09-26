module disso
  use module_base, only: gp
  use dictionaries, only: dictionary

  real(gp), dimension(:, :), allocatable, save :: en
  real(gp), dimension(3), save :: step, rxyz0
  integer, save :: iter, niter
  type(dictionary), pointer :: rerun

  character(len = 128) :: filename
end module disso

subroutine init(runObj)
  use module_base, only: bigdft_mpi, gp
  use module_f_objects
  use bigdft_run
  use dynamic_memory
  use yaml_output
  use public_keys, only: PLUGINS
  use dictionaries

  use disso

  implicit none

  type(run_objects), intent(inout) :: runObj

  integer :: sid
  type(kernel_ctx) :: kernel
  type(dictionary), pointer :: params
  real(gp), dimension(:,:), pointer :: rxyz

  external :: update, output, swap

  iter = 0
  niter = runObj%inputs%ncount_cluster_x
  en = f_malloc((/ 3, niter /), "en")
  params => runObj%user_inputs // PLUGINS // "disso" // "parameters"
  step = params // "delta"
  filename = params // "output"

  nullify(rerun)
  call dict_copy(rerun, params // "rerun")

  rxyz => bigdft_get_rxyz_ptr(runObj)
  rxyz0(:) = rxyz(:, 2)

  if (bigdft_mpi%iproc == 0) then
     call yaml_mapping_open("disso plugin")
     call yaml_map("step", step)
     call yaml_mapping_close()
  end if

  call f_object_kernel_new(kernel, update, 5)
  call f_object_kernel_add_arg(kernel, en)
  call f_object_kernel_add_arg(kernel, niter)
  call f_object_kernel_add_arg(kernel, iter)
  call f_object_signal_connect(RUN_OBJECTS_TYPE, POST_SCF_SIG, kernel)

  call f_object_kernel_new(kernel, swap, 3)
  call f_object_kernel_add_arg(kernel, niter)
  call f_object_kernel_add_arg(kernel, iter)
  call f_object_signal_connect(RUN_OBJECTS_TYPE, PRE_SCF_SIG, kernel)

  call f_object_kernel_new(kernel, output, 1)
  call f_object_signal_connect(RUN_OBJECTS_TYPE, DESTROY_SIG, kernel)
end subroutine init

subroutine swap(runObj, n, i)
  use bigdft_run
  use public_enums, only: ENUM_SCRATCH
  use module_input_keys, only: inputpsiid_set_policy
  use disso, only: rerun
  implicit none
  type(run_objects), intent(inout) :: runObj
  integer, intent(in) :: i, n
  
  if (i == n) then
     call run_objects_update(runObj, rerun)
     runObj%inputs%ncount_cluster_x = runObj%inputs%ncount_cluster_x * 2
     call inputpsiid_set_policy(ENUM_SCRATCH, runObj%inputs%inputPsiId)
  end if
end subroutine swap

subroutine update(runObj, outs, en, n, i)
  use module_base, only: gp
  use bigdft_run
  use wrapper_linalg
  use disso, only: step, rxyz0
  implicit none
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  integer, intent(in) :: n
  integer, intent(inout) :: i
  real(gp), dimension(3, n), intent(inout) :: en
  
  real(gp), dimension(:,:), pointer :: rxyz
  real(gp), dimension(3) :: d

  rxyz => bigdft_get_rxyz_ptr(runObj)
  d = rxyz(:, 1) - rxyz(:, 2)

  en(1, modulo(i, n) + 1) = nrm2(3, d(1), 1)
  en(2 + i / n, modulo(i, n) + 1) = outs%energy
  i = i + 1
  if (i == n) then
     rxyz(:, 2) = rxyz0(:)
     runObj%inputs%ncount_cluster_x = runObj%inputs%ncount_cluster_x * 2
  else
     rxyz(:, 2) = rxyz(:, 2) + step
  end if
end subroutine update

subroutine output(runObj)
  use module_base, only: gp
  use bigdft_run
  use dynamic_memory
  use disso
  use dictionaries
  implicit none
  type(run_objects), intent(in) :: runObj
  
  integer, parameter :: unit = 42
  
  open(unit = unit, file = trim(filename))
  write(unit, "(3F16.12)") en
  close(unit)
  call f_free(en)
  call dict_free(rerun)
end subroutine output
