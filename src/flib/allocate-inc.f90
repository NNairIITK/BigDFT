!  call timing(0,'Init to Zero  ','IR') 
  !then perform all the checks and profile the allocation procedure
  if (size(shape(array))==m%rank) then
      call metadata_address(array,iadd)
      call pad_array(array,m%put_to_zero,m%shape,ndebug)

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
!  call timing(0,'Init to Zero  ','RS') 
contains 

  subroutine profile(iadd)
    implicit none
    integer(kind=8), intent(in) :: iadd
    !local variables
    integer :: ierr,i,sizeof
    integer(kind=8) :: ilsize
    type(dictionary), pointer :: dict_tmp
    character(len=info_length) :: address
    
    !write the address of the first element in the address string
    call getaddress(array,address,len(address),ierr)

    !call profile_allocation(ierror,iadd,address,kind(array),m)

    !size
    sizeof=kind(array)
    do i=1,m%rank
       ilsize=int(sizeof*m%shape(i),kind=8)
    end do

    !create the dictionary array
    if (.not. associated(dict_routine)) then
       call dict_init(dict_routine)
    end if
    !add the array to the routine
    !call dict_array(m%routine_id,m%array_id,ilsize,dict_tmp)
    call dict_init(dict_tmp)
    call set(dict_tmp//arrayid,trim(m%array_id))
    call set(dict_tmp//sizeid,ilsize)
    call set(dict_tmp//routineid,trim(m%routine_id))
    call set(dict_tmp//metadatadd,iadd)

    call set(dict_routine//trim(address),dict_tmp)

    call check_for_errors(ierror,m%try)
    call memocc(ierror,product(m%shape(1:m%rank))*sizeof,m%array_id,m%routine_id)

  end subroutine profile

