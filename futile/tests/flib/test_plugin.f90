program test
  use yaml_output
  use module_f_objects, except => f_object_signal_add_str, except2 => f_object_signal_connect
  use dictionaries

  implicit none
  
  integer :: ierr
  character(len = max_field_length) :: pong
  character(len = 256) :: mess

  call f_lib_initialize()

  ! Define here some signals.
  call f_object_add_signal("class", "ping", 1)

  call f_object_signal_connect("class", "ping", me, 1, ierr)

  call plugin_load("pong", ierr)
  call yaml_map("status", ierr)
  if (ierr /= 0) then
     call plugin_error(mess)
     call yaml_map("error", mess)
  end if

  if (f_object_signal_prepare("class", "ping")) then
     call f_object_signal_add_str("class", "ping", pong)
     call f_object_signal_emit("class", "ping")
     call yaml_map("ping", pong)
  else
     call yaml_map("ping", "no listeners")
  end if

  call f_lib_finalize()

contains

  subroutine me(p)
    use dictionaries
    character(len = max_field_length), intent(out) :: p

    write(p, "(A)") "no answer"
  end subroutine me

end program test
