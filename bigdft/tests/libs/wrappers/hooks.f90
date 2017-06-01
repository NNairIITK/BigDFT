program test_hooks
  use module_base
  use bigdft_run
  use yaml_parse
  use yaml_output
  use dictionaries
  use module_f_objects
  implicit none 
  !input variables
  type(run_objects) :: runObj
  !output variables
  type(state_properties) :: outs
  character(len=60) :: posinp_id
  integer :: ierr, sid
  type(dictionary), pointer :: run, at

  call f_lib_initialize()

  call bigdft_init()

  call run_objects_type_init()
  call f_object_signal_connect("run_objects", "init", runInit, 1, sid)
  call f_object_signal_connect("run_objects", "pre", runPre, 1, sid)
  call f_object_signal_connect("run_objects", "post", runPost, 2, sid)
  call f_object_signal_connect("run_objects", "mix", runMix, 3, sid)
  call f_object_signal_connect("run_objects", "destroy", runDestroy, 1, sid)

  call dict_init(run)
  call dict_init(at)
  call set(at // "H" // 0, 0.)
  call set(at // "H" // 1, 0.)
  call set(at // "H" // 2, 0.)
  call add(run // "posinp" // "positions", at)
  call dict_init(at)
  call set(at // "H" // 0, 1.)
  call set(at // "H" // 1, 0.)
  call set(at // "H" // 2, 0.)
  call add(run // "posinp" // "positions", at)
  call set(run // "posinp" // "properties" // "format", "ascii")
  call add(run // "py_hooks" // "pre", 'print " number of atoms: %d" % run.nat(0)')
  call add(run // "py_hooks" // "post", 'print " forces: |\n%s" % prop.fxyz(futile.ARRAY)')
  call add(run // "py_hooks" // "post", 'print " energy: %s" % prop.energy(futile.SCALAR_R8)')

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
    use f_python
    type(run_objects), intent(in) :: obj

    integer :: ierr

    if (associated(obj%py_hooks)) then
       call f_python_initialize(bigdft_mpi%iproc, bigdft_mpi%nproc, &
            & bigdft_mpi%igroup, bigdft_mpi%ngroup)
       if ("init" .in. obj%py_hooks) then
          call f_python_execute_dict(obj%py_hooks // "init", ierr)
          if (ierr /= 0) stop ! Should raise a proper error later.
       end if
    end if
  end subroutine runInit

  subroutine runPre(obj)
    use f_python
    type(run_objects), intent(in) :: obj

    if ("pre" .in. obj%py_hooks) then
       call f_python_add_object("run_objects", "run", obj)
       call f_python_execute_dict(obj%py_hooks // "pre", ierr)
       if (ierr /= 0) stop
    end if
  end subroutine runPre

  subroutine runPost(obj, prop)
    use f_python
    type(run_objects), intent(in) :: obj
    type(state_properties), intent(in) :: prop

    if ("post" .in. obj%py_hooks) then
       call f_python_add_object("run_objects", "run", obj)
       call f_python_add_object("state_properties", "prop", prop)
       call f_python_execute_dict(obj%py_hooks // "post", ierr)
       if (ierr /= 0) stop
    end if
  end subroutine runPost

  subroutine runMix(obj, outs, subouts)
    type(run_objects), intent(in) :: obj
    type(state_properties), intent(out) :: outs
    type(state_properties), dimension(size(obj%sections)), intent(in) :: subouts

    !...
  end subroutine runMix

  subroutine runDestroy(obj)
    type(run_objects), intent(in) :: obj

    if ("destroy" .in. obj%py_hooks) then
       call f_python_add_object("run_objects", "run", obj)
       call f_python_execute_dict(obj%py_hooks // "destroy", ierr)
       if (ierr /= 0) stop
    end if
  end subroutine runDestroy

end program test_hooks
