!> @file
!! Include fortran file for deallocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  !local variables
  integer :: ierr
  character(len=info_length) :: address
  !local variables
  logical :: use_global
  integer(kind=8) :: ilsize,jlsize,iadd
  character(len=namelen) :: array_id,routine_id
  type(dictionary), pointer :: dict_add
  !      call timing(0,'AllocationProf','IR') 

  !address of first element (not needed for deallocation)
  !call getaddress(array,address,len(address),ierr)
  !address of the metadata 
  call metadata_address(array,iadd)
  address=repeat(' ',len(address))
  address=trim(long_toa(iadd))

  !here the size should be corrected with ndebug
  ilsize=int(product(shape(array))*kind(array),kind=8)
  !fortran deallocation
  deallocate(array,stat=ierror)
  !hopefully only address is necessary for the deallocation
  !in the include file this part can be expanded if needed
  !call profile_deallocation(ierror,ilsize,address) 

  call check_for_errors(ierror,.false.)
  !search in the dictionaries the address
  !a error event should be raised in this case
  dict_add=>find_key(dict_routine,trim(address))
  if (.not. associated(dict_add)) then
     dict_add=>find_key(dict_global,trim(address))
     if (.not. associated(dict_add)) then
        print *,'address:',trim(address)
        call f_malloc_dump_status()
        stop 'profile deallocations: address not present'
     else
        use_global=.true.
     end if
  else
     use_global=.false.
  end if
  !the global dictionary should be used instead
  array_id=dict_add//arrayid
  routine_id=dict_add//routineid
  jlsize=dict_add//sizeid
  if (ilsize /= jlsize) then
     print *,'array id',array_id,jlsize
     print *,'array id',ilsize,routine_id
     call f_malloc_dump_status()
     stop 'size error'
  end if

  call memocc(ierror,-int(ilsize),trim(array_id),trim(routine_id))

  if (use_global) then
     !call yaml_dict_dump(dict_global)
     call pop(dict_global,trim(address))
  else
     call pop(dict_routine,trim(address))
  end if


  !      call timing(0,'AllocationProf','RS') 
