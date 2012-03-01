subroutine get_command_argument(count, value, status)
  implicit none
  integer, intent(in) :: count
  integer, intent(out) :: status
  character(len = *), intent(out) :: value

  write(value, "(A)") ""
  status = 1
end subroutine get_command_argument
