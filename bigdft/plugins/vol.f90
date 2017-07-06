module vol
  use module_base, only: gp

  integer, save :: iter
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

  external :: setCell, testLoop

  params => runObj%user_inputs // PLUGINS // "vol" // "parameters"
  steps = f_malloc_ptr(dict_len(params // "steps"), "steps")
  steps = params // "steps"
  en = f_malloc_ptr((/ 2, size(steps) /), "en")
  filename = params // "output"
  alat0 = runObj%atoms%astruct%cell_dim
  iter = 1

  write(runObj%inputs%geopt_approach, "(A)") "LOOP"
  runObj%inputs%ncount_cluster_x = size(steps) + 1
  call set(runObj%user_inputs // GEOPT_VARIABLES // GEOPT_METHOD, "LOOP")
  call set(runObj%user_inputs // GEOPT_VARIABLES // NCOUNT_CLUSTER_X, size(steps) + 1)

  if (bigdft_mpi%iproc == 0) then
     call yaml_mapping_open("volume plugin")
     call yaml_map("steps", steps)
     call yaml_map("export", filename)
     call yaml_mapping_close()
  end if

  call f_object_kernel_new(kernel, setCell, 1)
  call f_object_signal_connect(RUN_OBJECTS_TYPE, PRE_SCF_SIG, kernel)

  call f_object_kernel_new(kernel, testLoop, 3)
  call f_object_signal_connect(PROCESS_RUN_TYPE, GEOPT_CONDITIONAL_SIG, kernel)
end subroutine init

subroutine setCell(runObj)
  use bigdft_run
  use dictionaries
  use public_keys, only: DFT_VARIABLES, INPUTPSIID
  use vol
  implicit none
  type(run_objects), intent(inout) :: runObj

  type(dictionary), pointer :: cell

  call dict_init(cell)
  call set(cell // "posinp" // "cell", alat0 * steps(iter))
  call set(cell // DFT_VARIABLES // INPUTPSIID, 1)
  call run_objects_update(runObj, cell)
  call dict_free(cell)
end subroutine setCell

subroutine testLoop(outs, it, check)
  use module_base, only: bigdft_mpi
  use bigdft_run
  use dictionaries
  use dynamic_memory
  use yaml_output
  use vol
  implicit none
  type(state_properties), intent(in) :: outs
  integer, intent(in) :: it
  integer, intent(out) :: check

  integer, parameter :: unit = 42
  integer :: i

  en(1, it) = product(alat0 * steps(it))
  en(2, it) = outs%energy
  iter = it + 1

  check = 0
  if (it == size(steps)) then
     if (bigdft_mpi%iproc == 0) then
        open(unit = unit, file = trim(filename))
        write(unit, "(2F18.12)") en
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
     check = 1
  end if
end subroutine testLoop
