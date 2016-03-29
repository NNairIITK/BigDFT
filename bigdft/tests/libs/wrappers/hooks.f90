program test_hooks
  use module_base
  use bigdft_run
  use yaml_parse
  use yaml_output
  use dictionaries
  implicit none 
  !input variables
  type(run_objects) :: runObj
  type(run_objects_class_type), pointer :: klass
  !output variables
  type(state_properties) :: outs
  character(len=60) :: posinp_id
  integer :: ierr
  type(dictionary), pointer :: run, at

  interface
     subroutine bigdft_python_exec_dict(dict, status)
       use dictionaries
       implicit none
       type(dictionary), pointer :: dict
       integer, intent(out) :: status
     end subroutine bigdft_python_exec_dict
  end interface

  call f_lib_initialize()

  call bigdft_init()

  klass => run_objects_class()
  call f_bind_define(klass%hook_init, runInit, 1)
  call f_bind_define(klass%hook_pre, runPre, 1)
  call f_bind_define(klass%hook_post, runPost, 2)
  call f_bind_define(klass%hook_destroy, runDestroy, 1)

  call dict_init(at)
  call set(at // "H" // 0, 0.)
  call set(at // "H" // 1, 0.)
  call set(at // "H" // 2, 0.)
  call dict_init(run)
  call set(run // "posinp" // "positions" // 0, at)
  call set(run // "posinp" // "properties" // "format", "ascii")
  call set(run // "py_hooks" // "init0" // 0, 'print "hello"')
  call set(run // "py_hooks" // "pre0" // 0, 'print run.atoms.nat')

  call run_objects_init(runObj, run)
  call init_state_properties(outs, bigdft_nat(runObj))
  call bigdft_get_run_properties(run, posinp_id = posinp_id)
  call dict_free(run)

  call process_run(posinp_id, runObj, outs)

  call deallocate_state_properties(outs)
  call free_run_objects(runObj)

  call bigdft_finalize(ierr)

  call f_lib_finalize()

contains

  subroutine runInit(obj)
    use yaml_output
    type(run_objects), intent(in) :: obj

    integer :: ierr

    write(*,*) "init"
    if (associated(obj%py_hooks)) then
       call bigdft_python_init(bigdft_mpi%iproc, bigdft_mpi%nproc, &
            & bigdft_mpi%igroup, bigdft_mpi%ngroup)
       if ("init0" .in. obj%py_hooks) then
          call bigdft_python_exec_dict(obj%py_hooks // "init0", ierr)
          if (ierr /= 0) stop ! Should raise a proper error later.
       end if
    end if
  end subroutine runInit

  subroutine runPre(obj)
    type(run_objects), intent(in) :: obj

    write(*,*) "pre"
    if ("pre0" .in. obj%py_hooks) then
       call run_objects_to_python(obj, "run", 3)
       !call bigdft_python_exec_dict(obj%py_hooks // "pre", ierr)
       !if (ierr /= 0) stop
    end if
  end subroutine runPre

  subroutine runPost(obj, prop)
    type(run_objects), intent(in) :: obj
    type(state_properties), intent(in) :: prop

    write(*,*) "post"
  end subroutine runPost

  subroutine runDestroy(obj)
    type(run_objects), intent(in) :: obj

    write(*,*) "destroy"
  end subroutine runDestroy

end program test_hooks
