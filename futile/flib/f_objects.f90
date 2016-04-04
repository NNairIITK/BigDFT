module module_f_objects

  use dictionaries

  private

  type(dictionary), pointer :: class_library => null()
  type(dictionary), pointer :: object_library => null()

  public :: f_object_new, f_object_add_method, f_object_get_method
  public :: f_object_finalize
!!$  public :: f_object_add_instance, f_object_remove_instance

contains

  subroutine ensure_init()
    if (.not. associated(class_library)) call dict_init(class_library)
    if (.not. associated(object_library)) call dict_init(object_library)
  end subroutine ensure_init

  subroutine f_object_finalize()
    call ensure_init()
    call dict_free(class_library)
    call dict_free(object_library)
  end subroutine f_object_finalize

  subroutine f_object_new(obj_id, constructor_add, destructor_add)
    use f_precisions
    character(len = *), intent(in) :: obj_id
    integer(f_address), intent(in) :: constructor_add, destructor_add
    
    call ensure_init()
    if (obj_id .in. class_library) stop

    call set(class_library // obj_id, "")
    call f_object_add_method(obj_id, "constructor", constructor_add, 1)
    call f_object_add_method(obj_id, "destructor", destructor_add, 0)
  end subroutine f_object_new
  
  subroutine f_object_add_method(obj_id, id, method_add, n_args)
    use f_precisions
    character(len = *), intent(in) :: obj_id, id
    integer(f_address), intent(in) :: method_add
    integer, intent(in) :: n_args

    call ensure_init()
    if (.not. (obj_id .in. class_library) .and. obj_id /= "class") stop
    
    call set(class_library // obj_id // id // "address", method_add)
    call set(class_library // obj_id // id // "n_args", n_args)
  end subroutine f_object_add_method

  subroutine f_object_get_method(obj_id, method_id, n_args, callback)
    use f_precisions
    character(len = *), intent(in) :: obj_id, method_id
    integer, intent(out) :: n_args
    integer(f_address), intent(out) :: callback

    n_args = 0
    callback = 0

    if (.not. (obj_id .in. class_library) .and. obj_id /= "class") return
    if (.not. (method_id .in. class_library // obj_id)) return

    n_args = class_library // obj_id // method_id // "n_args"
    callback = class_library // obj_id // method_id // "address"
  end subroutine f_object_get_method

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

