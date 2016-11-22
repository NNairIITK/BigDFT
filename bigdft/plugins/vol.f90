module vol
  use module_base, only: gp

  real(gp), dimension(3), save :: alat0
  real(gp), dimension(:), pointer :: steps
  real(gp), dimension(:, :), pointer :: en

  character(len = 128) :: filename
end module vol

subroutine init(runObj)
  use module_base, only: bigdft_mpi, gp
  use module_f_objects
  use bigdft_run
  use dynamic_memory
  use yaml_output
  use public_keys, only: PLUGINS, GEOPT_VARIABLES, NCOUNT_CLUSTER_X, GEOPT_METHOD
  use dictionaries
  use vol

  implicit none

  type(run_objects), intent(inout) :: runObj

  type(kernel_ctx) :: kernel
  type(dictionary), pointer :: params

  external :: compute, output

  params => runObj%user_inputs // PLUGINS // "vol" // "parameters"
  steps = f_malloc_ptr(dict_len(params // "steps"), "steps")
  steps = params // "steps"
  en = f_malloc_ptr((/ 2, size(steps) /), "en")
  filename = params // "output"
  alat0 = runObj%atoms%astruct%cell_dim

  write(runObj%inputs%geopt_approach, "(A)") "LOOP"
  runObj%inputs%ncount_cluster_x = size(steps)
  call set(runObj%user_inputs // GEOPT_VARIABLES // GEOPT_METHOD, "LOOP")
  call set(runObj%user_inputs // GEOPT_VARIABLES // NCOUNT_CLUSTER_X, size(steps))

  if (bigdft_mpi%iproc == 0) then
     call yaml_mapping_open("volume plugin")
     call yaml_map("steps", steps)
     call yaml_map("export", filename)
     call yaml_mapping_close()
  end if

  call f_object_kernel_new(kernel, compute, 4)
  call f_object_signal_connect(PROCESS_RUN_TYPE, GEOPT_LOOP_SIG, kernel)

  call f_object_kernel_new(kernel, output, 1)
  call f_object_signal_connect(RUN_OBJECTS_TYPE, DESTROY_SIG, kernel)
end subroutine init

subroutine compute(outs, runObj, it, check)
  use bigdft_run
  use dictionaries
  use public_enums, only: ENUM_MEMORY
  use module_input_keys, only: inputpsiid_set_policy
  use vol
  implicit none
  type(state_properties), intent(inout) :: outs
  type(run_objects), intent(inout) :: runObj
  integer, intent(in) :: it
  integer, intent(out) :: check

  type(dictionary), pointer :: cell
  integer :: infocode

  call dict_init(cell)
  call set(cell // "posinp" // "cell", alat0 * steps(it + 1))
  call run_objects_update(runObj, cell)
  call dict_free(cell)

  call inputpsiid_set_policy(ENUM_MEMORY, runObj%inputs%inputPsiId)
  call bigdft_state(runObj, outs, infocode)

  en(1, it + 1) = product(alat0 * steps(it + 1))
  en(2, it + 1) = outs%energy

  check = 0
  if (it + 1 == size(steps)) check = 1
end subroutine compute

subroutine output(runObj)
  use module_base, only: bigdft_mpi
  use bigdft_run
  use dynamic_memory
  use yaml_output
  use vol
  implicit none
  type(run_objects), intent(in) :: runObj
  
  integer, parameter :: unit = 42
  integer :: i

  if (.not. associated(runObj%sections)) return
  
  if (bigdft_mpi%iproc == 0) then
     open(unit = unit, file = trim(filename))
     write(unit, "(2F16.12)") en
     close(unit)
     call yaml_sequence_open("Energy (Hartree) per volume (Bohr^3)")
     do i = 1, size(steps)
        call yaml_sequence(advance = "no")
        call yaml_mapping_open(flow = .true.)
        call yaml_map("volume", en(1, i))
        call yaml_map("energy", en(2, i))
        call yaml_mapping_close()
     end do
     call yaml_sequence_close()
  end if
  call f_free_ptr(en)
  call f_free_ptr(steps)
end subroutine output
