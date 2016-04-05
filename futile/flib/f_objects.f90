module module_f_objects
  use dictionaries

  implicit none

  private

  type(dictionary), pointer :: class_library => null()

  public :: f_object_new, f_object_add_method, f_object_get_method
  public :: f_object_finalize
  public :: f_object_add_signal, f_object_signal_prepare, f_object_signal_add_arg
  public :: f_object_signal_emit, f_object_signal_connect

contains

  subroutine ensure_init()
    implicit none
    if (.not. associated(class_library)) call dict_init(class_library)
  end subroutine ensure_init

  subroutine f_object_finalize()
    implicit none
    if (associated(class_library)) call dict_free(class_library)
  end subroutine f_object_finalize

  subroutine f_object_new(obj_id, constructor_add, destructor_add)
    use f_precisions
    implicit none
    character(len = *), intent(in) :: obj_id
    integer(f_address), optional, intent(in) :: constructor_add, destructor_add
    
    call ensure_init()
    if (obj_id .in. class_library) stop

    call set(class_library // obj_id, "")
    if (present(constructor_add)) &
         & call f_object_add_method(obj_id, "constructor", constructor_add, 1)
    if (present(destructor_add)) &
         & call f_object_add_method(obj_id, "destructor", destructor_add, 0)
  end subroutine f_object_new
  
  subroutine f_object_add_method(obj_id, id, method_add, n_args)
    use f_precisions
    implicit none
    character(len = *), intent(in) :: obj_id, id
    integer(f_address), intent(in) :: method_add
    integer, intent(in) :: n_args

    call ensure_init()
    if (.not. (obj_id .in. class_library) .and. obj_id /= "class") stop
    
    call set(class_library // obj_id // "methods" // id // "address", method_add)
    call set(class_library // obj_id // "methods" // id // "n_args", n_args)
  end subroutine f_object_add_method

  subroutine f_object_get_method(obj_id, method_id, n_args, callback)
    use f_precisions
    implicit none
    character(len = *), intent(in) :: obj_id, method_id
    integer, intent(out) :: n_args
    integer(f_address), intent(out) :: callback

    n_args = 0
    callback = 0

    if (.not. (obj_id .in. class_library) .and. obj_id /= "class") return
    if (.not. (method_id .in. class_library // obj_id // "methods")) return

    n_args = class_library // obj_id // "methods" // method_id // "n_args"
    callback = class_library // obj_id // "methods" // method_id // "address"
  end subroutine f_object_get_method

  subroutine f_object_add_signal(obj_id, id, n_args)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    integer, intent(in) :: n_args

    call ensure_init()
    if (.not. (obj_id .in. class_library) .and. obj_id /= "class") stop
    
    call set(class_library // obj_id // "signals" // id // "n_args", n_args)
  end subroutine f_object_add_signal

  function ensure_signal(obj_id, id)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    logical :: ensure_signal

    call ensure_init()
    if (.not. (obj_id .in. class_library) .and. obj_id /= "class") stop
    if (.not. (id .in. class_library // obj_id // "signals")) stop
    
    ensure_signal = .true.
  end function ensure_signal

  function f_object_signal_prepare(obj_id, id) result(emit)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    logical :: emit

    if (.not. ensure_signal(obj_id, id)) return

    if (dict_len(class_library // obj_id // "signals" // id // "arguments") > 0) &
         & call dict_remove(class_library // obj_id // "signals" // id, "arguments")
    emit = (dict_len(class_library // obj_id // "signals" // id // "hooks") > 0)
  end function f_object_signal_prepare

  subroutine f_object_signal_add_arg(obj_id, id, arg_add)
    use f_precisions
    implicit none
    character(len = *), intent(in) :: obj_id, id
    integer(f_address), intent(in) :: arg_add

    if (.not. ensure_signal(obj_id, id)) return

    call add(class_library // obj_id // "signals" // id // "arguments", arg_add)
  end subroutine f_object_signal_add_arg

  subroutine f_object_signal_emit(obj_id, id)
    use f_precisions
    implicit none
    character(len = *), intent(in) :: obj_id, id
    
    type(dictionary), pointer :: iter, signal
    integer(f_address) :: callback
    integer(f_address), dimension(2) :: args
    integer :: n_args

    if (.not. ensure_signal(obj_id, id)) return

    signal => class_library // obj_id // "signals" // id
    n_args = signal // "n_args"
    if (dict_len(signal // "arguments") /= n_args) stop
    if (n_args > 0) args(1:n_args) = signal // "arguments"

    iter => dict_iter(signal // "hooks")
    do while (associated(iter))
       callback = iter // "address"
       select case(n_args)
       case (0)
          call call_external_c_fromadd(callback)
       case (1)
          call call_external_c_fromadd_data(callback, args(1))
       case (2)
          call call_external_c_fromadd_data_data(callback, args(1), args(2))
       end select
       iter => dict_next(iter)
    end do
  end subroutine f_object_signal_emit

  subroutine f_object_signal_connect(obj_id, id, hook_add, n_args, sid)
    use f_precisions
    implicit none
    character(len = *), intent(in) :: obj_id, id
    integer(f_address), intent(in) :: hook_add
    integer, intent(in) :: n_args
    integer, intent(out) :: sid

    type(dictionary), pointer :: hook
    integer :: n_args_signal

    sid = -1

    if (.not. ensure_signal(obj_id, id)) return

    n_args_signal = class_library // obj_id // "signals" // id // "n_args"
    if (n_args_signal /= n_args) stop

    ! Get the last hook, to retrieve its id.
    hook => class_library // obj_id // "signals" // id // "hooks"
    if (dict_len(hook) > 0) then
       hook => hook // (dict_len(hook) - 1)
       sid = hook // "id"
       sid = sid + 1
    else
       sid = 0
    end if

    call dict_init(hook)
    call set(hook // "id", sid)
    call set(hook // "address", hook_add)

    call add(class_library // obj_id // "signals" // id // "hooks", hook)
  end subroutine f_object_signal_connect

end module module_f_objects




subroutine f_object_new(id, constructor, destructor)
  use f_precisions
  use module_f_objects, only: wrapper_new => f_object_new
  character(len = *), intent(in) :: id
  external :: constructor, destructor
  
  call wrapper_new(id, f_loc(constructor), f_loc(destructor))
end subroutine f_object_new

subroutine f_object_add_method(obj_id, id, method, n_args)
  use f_precisions
  use module_f_objects, only: wrapper_add => f_object_add_method
  character(len = *), intent(in) :: obj_id, id
  integer, intent(in) :: n_args
  external :: method
  
  call wrapper_add(obj_id, id, f_loc(method), n_args)
end subroutine f_object_add_method

subroutine f_object_get_method(obj_id, method_id, n_args, callback)
  use f_precisions
  use module_f_objects, only: wrapper_get => f_object_get_method
  character(len = *), intent(in) :: obj_id, method_id
  integer, intent(out) :: n_args
  integer(f_address), intent(out) :: callback
  
  call wrapper_get(obj_id, method_id, n_args, callback)
end subroutine f_object_get_method

subroutine f_object_signal_add_arg(obj_id, id, arg)
  use f_precisions
  use module_f_objects, only: wrapper_add => f_object_signal_add_arg
  character(len = *), intent(in) :: obj_id, id
  external :: arg

  call wrapper_add(obj_id, id, f_loc(arg))
end subroutine f_object_signal_add_arg

subroutine f_object_signal_connect(obj_id, id, hook, n_args, sid)
  use f_precisions
  use module_f_objects, only: wrapper_connect => f_object_signal_connect
  character(len = *), intent(in) :: obj_id, id
  external :: hook
  integer, intent(out) :: sid

  call wrapper_connect(obj_id, id, f_loc(hook), n_args, sid)
end subroutine f_object_signal_connect
