subroutine wrtout(unit, message, mpi_trt)
  implicit none
  
  integer, intent(in)              :: unit
  character(len = 500), intent(in) :: message
  character(len = *), intent(in)   :: mpi_trt
  
  write(unit, "(A)") trim(message)
end subroutine wrtout
