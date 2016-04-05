module my_object

  type my_object_type
     integer :: a_value
  end type my_object_type

contains

  subroutine my_object_type_init()
    use module_f_objects, only: f_object_new, f_object_add_signal

    call f_object_new("my_object")
    call f_object_add_signal("my_object", "serialize", 1)
    call f_object_add_signal("my_object", "add", 2)
    call f_object_add_signal("my_object", "value-changed", 1)
  end subroutine my_object_type_init
  
  subroutine add(obj, obj2)
    use module_f_objects, only: f_object_signal_prepare, f_object_signal_emit
    type(my_object_type), intent(inout) :: obj
    type(my_object_type), intent(in) :: obj2

    integer :: old_value

    old_value = obj%a_value

    if (f_object_signal_prepare("my_object", "add")) then
       call f_object_signal_add_arg("my_object", "add", obj)
       call f_object_signal_add_arg("my_object", "add", obj2)
       call f_object_signal_emit("my_object", "add")
    end if

    if (old_value /= obj%a_value .and. &
         & f_object_signal_prepare("my_object", "value-changed")) then
       call f_object_signal_add_arg("my_object", "value-changed", obj)
       call f_object_signal_emit("my_object", "value-changed")
    end if
  end subroutine add

  subroutine serialize(obj)
    use module_f_objects, only: f_object_signal_prepare, f_object_signal_emit
    type(my_object_type), intent(inout) :: obj
    
    if (f_object_signal_prepare("my_object", "serialize")) then
       call f_object_signal_add_arg("my_object", "serialize", obj)
       call f_object_signal_emit("my_object", "serialize")
    end if
  end subroutine serialize

end module my_object

subroutine objtoyaml(obj)
  use yaml_output
  use my_object
  type(my_object_type), intent(in) :: obj
  
  call yaml_mapping_open("obj")
  call yaml_map("value", obj%a_value)
  call yaml_mapping_close()
end subroutine objtoyaml

subroutine changed1(obj)
  use yaml_output
  use my_object
  type(my_object_type), intent(in) :: obj
  
  call yaml_map("obj changed", .true.)
end subroutine changed1

subroutine changed2(obj)
  use yaml_output
  use my_object
  type(my_object_type), intent(in) :: obj
  
  call yaml_map("new value", obj%a_value)
end subroutine changed2

subroutine objadd(obj, obj2)
  use my_object
  type(my_object_type), intent(inout) :: obj
  type(my_object_type), intent(in) :: obj2

  obj%a_value = obj%a_value + obj2%a_value
end subroutine objadd

program test_hooks
  use my_object
  type(my_object_type) :: o, o2
  external :: objtoyaml, objadd, changed1, changed2

  integer :: sid

  call f_lib_initialize()

  call my_object_type_init()

  call f_object_signal_connect("my_object", "serialize", objtoyaml, 1, sid)
  call f_object_signal_connect("my_object", "add", objadd, 2, sid)
  call f_object_signal_connect("my_object", "value-changed", changed1, 1, sid)
  call f_object_signal_connect("my_object", "value-changed", changed2, 1, sid)

  o%a_value = 42
  call serialize(o)
  
  o2%a_value = 41
  call add(o, o2)
  call serialize(o)

  call f_lib_finalize()
end program test_hooks
