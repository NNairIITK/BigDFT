!> @file
!! Include fortran file for allocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!  call timing(0,'Init to Zero  ','IR') 
  !then perform all the checks and profile the allocation procedure
  if (f_err_raise(ierror/=0,&
       'Allocation problem, error code '//trim(yaml_toa(ierror)),ERR_ALLOCATE)) return
  if (size(shape(array))==m%rank) then
     call pad_array(array,m%put_to_zero,m%shape,ndebug)
     !profile the array allocation
     if (m%profile) then
        call metadata_address(array,iadd)
        call profile(iadd)
     end if
     !ierror=SUCCESS
  else
     call f_err_throw('Rank specified by f_malloc ('//trim(yaml_toa(m%rank))//&
          ') is not coherent with the one of the array ('//trim(yaml_toa(size(shape(array))))//')',&
          ERR_INVALID_MALLOC)
     return
     !ierror=INVALID_RANK
     !lasterror='rank not valid'
     !call check_for_errors(ierror,m%try)
     !if (m%try) deallocate(array,stat=ierror)
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

    !size
    sizeof=kind(array)
    ilsize=int(sizeof*product(m%shape(1:m%rank)),kind=8)
    !do i=1,m%rank
    !   ilsize=int(sizeof*m%shape(i),kind=8)
    !end do

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
    !call set(dict_tmp//metadatadd,iadd)
    call set(dict_tmp//firstadd,trim(address))

    !call set(dict_routine//trim(address),dict_tmp)
    call set(dict_routine//trim(long_toa(iadd)),dict_tmp)

    !call check_for_errors(ierror,m%try)
    call memocc(ierror,product(m%shape(1:m%rank))*sizeof,m%array_id,m%routine_id)

  end subroutine profile

