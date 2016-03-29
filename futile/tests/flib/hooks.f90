module my_object
  use module_f_bind

  type my_object_type
     integer :: a_value
     type(f_bind) :: serialize, add
  end type my_object_type

contains
  
  subroutine add(obj, obj2)
    type(my_object_type), intent(inout) :: obj
    type(my_object_type), intent(in) :: obj2
    
    call f_bind_prepare(obj%add)
    call f_bind_add_arg(obj%add, obj)
    call f_bind_add_arg(obj%add, obj2)
    call f_bind_execute(obj%add)
  end subroutine add

  subroutine serialize(obj)
    type(my_object_type), intent(inout) :: obj
    
    call f_bind_prepare(obj%serialize)
    call f_bind_add_arg(obj%serialize, obj)
    call f_bind_execute(obj%serialize)
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

subroutine objadd(obj, obj2)
  use my_object
  type(my_object_type), intent(inout) :: obj
  type(my_object_type), intent(in) :: obj2

  obj%a_value = obj%a_value + obj2%a_value
end subroutine objadd

program test_hooks
  use my_object
  type(my_object_type) :: o, o2
  external :: objtoyaml, objadd

  call f_lib_initialize()

  o%a_value = 42
  call f_bind_define(o%serialize, objtoyaml, 1)
  call serialize(o)
  
  o2%a_value = 41
  call f_bind_define(o%add, objadd, 2)
  call add(o, o2)
  call serialize(o)

  call f_lib_finalize()
end program test_hooks
