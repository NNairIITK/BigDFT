  call timing(0,'Init to Zero  ','IR') 
  !then perform all the checks and profile the allocation procedure
  if (size(shape(array))==m%rank) then
      call metadata_address(array,iadd)

!!$     if (m%rank ==1) then
!!$        call pad_array(array,m%put_to_zero,m%shape(1),m%shape(1)+ndebug)
!!$     else
!!$        call pad_array(array,m%put_to_zero,&
!!$             product(m%shape(1:m%rank)),&
!!$             product(m%shape(1:m%rank-1))*(m%shape(m%rank)+ndebug))
!!$     end if
     !if (ndebug /=0) call pad_with_nan(array,m%rank,m%shape)
     !profile the array allocation
     call profile(iadd)
     ierror=SUCCESS
  else
     ierror=INVALID_RANK
     lasterror='rank not valid'
     deallocate(array,stat=ierror)
     call check_for_errors(ierror,m%try)
  end if
  call timing(0,'Init to Zero  ','RS') 
contains 

  subroutine profile(iadd)
    implicit none
    integer(kind=8), intent(in) :: iadd
    integer :: ierr
    character(len=info_length) :: address
    
    !write the address of the first element in the address string
    call getaddress(array,address,len(address),ierr)
    call profile_allocation(ierror,iadd,address,kind(array),m)

  end subroutine profile

