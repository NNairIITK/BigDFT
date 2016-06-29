!!! [init]
subroutine init()
!!! [init]
  use yaml_output

  integer :: sid

  interface
     subroutine pong(p)
       use dictionaries
       character(len = max_field_length), intent(out) :: p
     end subroutine pong
  end interface

  call yaml_map("plugin", "initialised")

  call f_object_signal_connect("class", "ping", pong, 1, sid)
end subroutine init

subroutine pong(p)
  use dictionaries
  character(len = max_field_length), intent(out) :: p

  write(p, "(A)") "pong"
end subroutine pong
