!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!$  if (ierror/=0) then
!!$     !$ if(not_omp) then
!!$     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!!$     !$ end if
!!$     call f_err_throw('Allocation problem, error code '//trim(yaml_toa(ierror)),ERR_ALLOCATE)
!!$     return
!!$  end if
!!$  if (size(shape(array))==m%rank) then
!!$     call pad_array(array,m%put_to_zero,m%shape,ndebug)
!!$     !also fill the array with the values of the source if the address is identified in the source
!!$     if (m%srcdata_add /= 0) call c_memcopy(array,m%srcdata_add,product(shape(array))*kind(array))
!!$     !profile the array allocation
!!$     if (m%profile) then
!!$        sizeof=kind(array)
!!$        ilsize=max(int(sizeof,kind=8)*int(product(m%shape(1:m%rank)),kind=8),int(0,kind=8))
!!$        if (track_origins) then
!!$           !write the address of the first element in the address string
!!$           call getlongaddress(array,iadd)
!!$           !store information only for array of size /=0
!!$           if (ilsize /= int(0,kind=8)) then
!!$              !create the dictionary array
!!$              if (.not. associated(mems(ictrl)%dict_routine)) then
!!$                 call dict_init(mems(ictrl)%dict_routine)
!!$              end if
!!$              call set(mems(ictrl)%dict_routine//long_toa(iadd),&
!!$                   !dict_new(arrayid .is. trim(m%array_id),&
!!$                   !routineid .is. trim(m%routine_id),&
!!$                   !sizeid .is. trim(yaml_toa(ilsize)),&
!!$                   !'Rank' .is. trim(yaml_toa(m%rank))))
!!$              '[ '//trim(m%array_id)//', '//trim(m%routine_id)//', '//&
!!$               trim(yaml_toa(ilsize))//', '//trim(yaml_toa(m%rank))//']')
!!$           end if
!!$        end if
!!$        call memstate_update(memstate,ilsize,m%array_id,m%routine_id)
!!$     end if
!!$  else
!!$     !$ if(not_omp) then
!!$     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!!$     !$ end if
!!$     call f_err_throw('Rank specified by f_malloc ('//trim(yaml_toa(m%rank))//&
!!$          ') is not coherent with the one of the array ('//trim(yaml_toa(size(shape(array))))//')',&
!!$          ERR_INVALID_MALLOC)
!!$     return
!!$  end if
!!$  !$ if(not_omp) then
!!$  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!!$  !$ end if
  

  if (ierror/=0) then
     !$ if(not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('Allocation problem, error code '//trim(yaml_toa(ierror)),ERR_ALLOCATE)
     return
  end if
  if (size(shape(array))==m%rank) then
     call pad_array(array,m%put_to_zero,m%shape,ndebug)
     !also fill the array with the values of the source if the address is identified in the source
     if (m%srcdata_add /= 0) call c_memcopy(array,m%srcdata_add,product(shape(array))*kind(array))
     !profile the array allocation
     iadd=int(0,kind=8)
     !write the address of the first element in the address string
     if (m%profile .and. track_origins) call getlongaddress(array,iadd)

     call f_update_database(product(int(m%shape(1:m%rank),kind=8)),kind(array),m%rank,&
          iadd,m%array_id,m%routine_id)

  else
     !$ if(not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('Rank specified by f_malloc ('//trim(yaml_toa(m%rank))//&
          ') is not coherent with the one of the array ('//trim(yaml_toa(size(shape(array))))//')',&
          ERR_INVALID_MALLOC)
     return
  end if
  !$ if(not_omp) then
  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
  !$ end if
