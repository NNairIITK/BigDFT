  !then perform all the checks and profile the allocation procedure
  if (size(shape(array))==m%rank) then
     call pad_with_nan(array,m%rank,m%shape)
     !profile the array allocation
     call profile()
     ierror=SUCCESS
  else
     ierror=INVALID_RANK
     lasterror='rank not valid'
     deallocate(array,stat=ierror)
     call check_for_errors(ierror,m%try)
  end if

contains 

  subroutine profile()
    implicit none
    integer :: ierr
    character(len=info_length) :: address

    call getaddress(array,address,len(address),ierr)
    call profile_allocation(ierror,address,kind(array),m)

  end subroutine profile
