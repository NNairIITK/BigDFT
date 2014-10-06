!> @file
!! Include fortran file for deallocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!$  !$ if (not_omp) then
!!$  call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
!!$  !$ end if
!!$
!!$  !here the size should be corrected with ndebug (or maybe not)
!!$  ilsize=int(kind(array),kind=8)*product(int(shape(array),kind=8))
!!$!  ilsize=int(product(shape(array))*kind(array),kind=8)
!!$  !retrieve the address of the first element if the size is not zero
!!$  iadd=int(0,kind=8)
!!$  if (ilsize /= int(0,kind=8)) call getlongaddress(array,iadd)
!!$  !fortran deallocation
!!$  deallocate(array,stat=ierror)
!!$
!!$  if (ierror/=0) then
!!$     !$ if (not_omp) then
!!$     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!!$     !$ end if
!!$     call f_err_throw('Deallocation problem, error code '//trim(yaml_toa(ierror)),&
!!$          ERR_DEALLOCATE)
!!$     return
!!$  end if
!!$
!!$  !profile address, in case of profiling activated
!!$  !  if (m%profile) then 
!!$  !address of first element (not needed for deallocation)
!!$  if (track_origins .and. iadd/=int(0,kind=8)) then
!!$     !hopefully only address is necessary for the deallocation
!!$
!!$     !search in the dictionaries the address
!!$     dict_add=>find_key(mems(ictrl)%dict_routine,long_toa(iadd))
!!$     if (.not. associated(dict_add)) then
!!$        dict_add=>find_key(mems(ictrl)%dict_global,long_toa(iadd))
!!$        if (.not. associated(dict_add)) then
!!$           !$ if (not_omp) then
!!$           call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!!$           !$ end if
!!$           call f_err_throw('Address '//trim(long_toa(iadd))//&
!!$                ' not present in dictionary',ERR_INVALID_MALLOC)
!!$           return
!!$        else
!!$           use_global=.true.
!!$        end if
!!$     else
!!$        use_global=.false.
!!$     end if
!!$
!!$     !transform the dict_add in a list
!!$     !retrieve the string associated to the database
!!$     array_info=dict_add
!!$     dict_add => yaml_a_todict(array_info)
!!$     !then retrieve the array information
!!$     array_id=dict_add//0
!!$     routine_id=dict_add//1
!!$     jlsize=dict_add//2
!!$
!!$     call dict_free(dict_add)
!!$     
!!$     if (ilsize /= jlsize) then
!!$        !$ if (not_omp) then
!!$        call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!!$        !$ end if
!!$        call f_err_throw('Size of array '//trim(array_id)//&
!!$             ' ('//trim(yaml_toa(ilsize))//') not coherent with dictionary, found='//&
!!$             trim(yaml_toa(jlsize)),ERR_MALLOC_INTERNAL)
!!$        return
!!$     end if
!!$     if (use_global) then
!!$        !call yaml_dict_dump(dict_global)
!!$        call dict_remove(mems(ictrl)%dict_global,long_toa(iadd))
!!$     else
!!$        call dict_remove(mems(ictrl)%dict_routine,long_toa(iadd))
!!$     end if
!!$  else
!!$     array_id(1:len(array_id))=arrayid
!!$     routine_id(1:len(routine_id))=routineid
!!$  end if
!!$
!!$  call memstate_update(memstate,-ilsize,trim(array_id),trim(routine_id))
!!$
!!$  !$ if (not_omp) then
!!$  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!!$  !$ end if
  
  !$ if (not_omp) then
  call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
  !$ end if

  !here the size should be corrected with ndebug (or maybe not)
  ilsize=product(int(shape(array),kind=8))
  !retrieve the address of the first element if the size is not zero
  iadd=int(0,kind=8)
  if (ilsize /= int(0,kind=8)) call getlongaddress(array,iadd)
  !fortran deallocation
  deallocate(array,stat=ierror)

  if (ierror/=0) then
     !$ if (not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('Deallocation problem, error code '//trim(yaml_toa(ierror)),&
          ERR_DEALLOCATE)
     return
  end if

  call f_purge_database(ilsize,kind(array),iadd)

  !$ if (not_omp) then
  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
  !$ end if
  
  
