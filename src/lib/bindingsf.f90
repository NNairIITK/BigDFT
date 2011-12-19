subroutine memocc_report()
  use m_profiling, only: mreport => memocc_report

  call mreport()
end subroutine memocc_report

subroutine deallocate_double(array)
  use module_base
  implicit none

  double precision, dimension(:,:), pointer :: array
  integer :: i_all, i_stat

  if (associated(array)) then
     i_all=-product(shape(array))*kind(array)
     deallocate(array,stat=i_stat)
     call memocc(i_stat,i_all,'array',"deallocate_double")
  end if
end subroutine deallocate_double

subroutine glr_new(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  allocate(glr)
end subroutine glr_new
subroutine glr_free(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  deallocate(glr)
end subroutine glr_free
subroutine glr_get_n(glr, n)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: glr
  integer, dimension(3), intent(out) :: n

  n(1) = glr%d%n1
  n(2) = glr%d%n2
  n(3) = glr%d%n3
end subroutine glr_get_n
