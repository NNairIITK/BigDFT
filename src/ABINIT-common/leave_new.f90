subroutine leave_new(mpi_trt)
  implicit none
  
  character(len=4),intent(in) :: mpi_trt

  write(6, "(A)") "leave_new : decision taken to exit ..."
  stop 1
end subroutine leave_new
