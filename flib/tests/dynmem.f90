!> @file
!! Test the dynamic memory allocation of the flib library
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test the dynamic memory allocation (flib)
subroutine test_dynamic_memory()
   use yaml_output
   use dynamic_memory
   use dictionaries
   implicit none

   type :: dummy_type
      integer :: i
      double precision, dimension(:), pointer :: ptr
   end type dummy_type

   !logical :: fl
   integer :: i,id
   real(kind=8), dimension(:), allocatable :: density,rhopot,potential,pot_ion,xc_pot
   real(kind=8), dimension(:), pointer :: extra_ref
   real(kind=8), dimension(:,:), allocatable :: ab
   integer, dimension(:), allocatable :: i1_all,i1_src
   integer, dimension(:), pointer :: i1_ptr,ptr1
   integer, dimension(:), pointer :: ptr2

   integer,dimension(:,:,:),allocatable :: weight
   integer,dimension(:,:,:,:),allocatable :: orbital_id
   type(dummy_type) :: dummy_test
   logical :: within_openmp
   real(kind=8) :: total
   character(len=250) :: message
   external :: abort2
   integer ithread
   !$ logical :: omp_in_parallel, omp_get_nested
   !$ integer :: omp_get_thread_num
  ithread=0
   call yaml_comment('Routine-Tree creation example',hfill='~')
   !call dynmem_sandbox()

   i1_all=f_malloc(0,id='i1_all')
   call yaml_map('Address of first element',f_loc(i1_all))
!   call yaml_map('Address of first element, explicit',f_loc(i1_all(1)))

   nullify(ptr1,ptr2)
   call yaml_map('Associated',(/associated(ptr1),associated(ptr2)/))

   i1_ptr=f_malloc_ptr(0,id='i1_ptr')

   ptr1=>i1_ptr

   ptr2=>ptr1
   call yaml_map('Associated',(/associated(ptr1),associated(ptr2)/))

   call yaml_map('Address of first element',f_loc(i1_ptr))
   call yaml_map('Address of first element, explicit (1)',f_loc(ptr1))
   call yaml_map('Address of first element, explicit (2)',f_loc(ptr2))

   !also allocating the dummy structure
   dummy_test=dummy_init(10)

!call f_malloc_finalize(dump=.true.)
!stop
   !call init_dummy(10,dummy_test)
   call free_dummy(dummy_test)

   call f_free_ptr(ptr2) !this one would crash
   !call f_free_ptr(i1_ptr)

   call yaml_map('Associated',(/associated(ptr1),associated(ptr2),associated(i1_ptr)/))
!!$
!!$
!!$
!!$   call f_malloc_dump_status()
   call f_free(i1_all)


call yaml_comment('debug 1')
   call f_routine(id='PS_Check')
call yaml_comment('debug 2')
   call f_routine(id='Routine 0')
   !Density
   density=f_malloc(3*2,id='density')
   !Density then potential
   potential=f_malloc(3,id='potential')!0

   call f_release_routine()

   call f_routine(id='Routine A')
weight=f_malloc((/1.to.8,1.to.8,2.to.4/),id='weight')
weight(1,1,2)=5
call f_free(weight)
   call f_release_routine()

!!$   call yaml_open_map('TemporaryA')
!!$    call f_malloc_dump_status()
!!$    call yaml_close_map()
!!$
!!$
!!$   call f_routine(id='Routine A')
!!$weight=f_malloc((/1.to.8,1.to.8,2.to.4/),id='weight')
!!$weight(1,1,2)=5
!!$call f_free(weight)
!!$   call f_release_routine()
!!$
!!$   call yaml_open_map('TemporaryB')
!!$    call f_malloc_dump_status()
!!$    call yaml_close_map()


!!$   call f_release_routine()
!!$
!!$   call yaml_open_map('Temporary')
!!$    call f_malloc_dump_status()
!!$    call yaml_close_map()
!!$stop
!!$
!!$!   call f_malloc_dump_status()
!!$   call f_routine(id=subname)
!!$
!!$   ncommarr=f_malloc(lbounds=(/0,1/),ubounds=(/nproc-1,4/),id='ncommarr')
!!$   ncommarr=f_malloc((/0.to.nproc,1.to.4/),id='ncommarr')
!!$   ncommarr=f_malloc_ptr((/0.to.nproc-1,4/),id='ncommarr')
!!$   call f_release_routine()

   call f_routine(id='Routine D')
    call f_routine(id='SubCase 1')
    !test of copying data
    i1_src=f_malloc(63,id='i1_src')

    i1_src=(/(i+24,i=1,size(i1_src))/)
    call yaml_map('Source vector',i1_src)
    call yaml_map('Its shape',shape(i1_src))
    call yaml_map('Its rank',size(shape(i1_src)))
    i1_all=f_malloc(src=i1_src,id='i1_all')

    call yaml_map('Vector difference',i1_all-i1_src)

!!$    call yaml_map('Vector allocation',[allocated(i1_all),allocated(i1_src)])
!!$    call f_free(i1_src)
!!$    call yaml_map('Vector allocation',[allocated(i1_all),allocated(i1_src)])
!!$    call f_free(i1_all)
!!$    call yaml_map('Vector allocation',[allocated(i1_all),allocated(i1_src)])

    call f_free(i1_src,i1_all)
    call f_release_routine()
    call f_routine(id='Subcase 2')
     call f_routine(id='SubSubcase1')
     !test of copying data
     i1_ptr=f_malloc_ptr(34,id='i1_ptr')

     i1_ptr=(/(i+789,i=1,size(i1_ptr))/)
     call yaml_map('Source pointer',i1_ptr)
     call yaml_map('Its shape',shape(i1_ptr))
     call yaml_map('Its rank',size(shape(i1_ptr)))
     ptr1=f_malloc_ptr(src=i1_ptr,id='ptr1')

     call yaml_map('Pointer difference',i1_ptr-ptr1)
!!$     call yaml_map('Pointer association',[associated(i1_ptr),associated(ptr1)])
!!$     call f_free_ptr(i1_ptr)
!!$     call yaml_map('Pointer association',[associated(i1_ptr),associated(ptr1)])
!!$     call f_free_ptr(ptr1)
!!$     call yaml_map('Pointer association',[associated(i1_ptr),associated(ptr1)])

     !again, to see if the copy works
     i1_ptr=i1_ptr-1
     call yaml_map('Pointer difference again (should be -1)',i1_ptr-ptr1)
     !then copy
     !call f_memcpy(src=i1_ptr(1),dest=ptr1(1),n=size(i1_ptr))
     call f_memcpy(src=i1_ptr,dest=ptr1)
     call yaml_map('Pointer difference again (should be 0)',i1_ptr-ptr1)

     call f_free_ptr(i1_ptr,ptr1)
     call f_release_routine()
    call f_release_routine()
    call f_routine(id='SubCase 3')
      weight    =f_malloc((/1.to.1,1.to.1,-1.to.-1/),id='weight')
      orbital_id=f_malloc((/1.to.1,1.to.1,1.to.7,0.to.0/),id='orbital_id')
    call f_malloc_dump_status()
      call f_free(weight)
      call f_free(orbital_id)
    call f_release_routine()
   call f_release_routine()
   !repeat it to see if it gives errors
   !the point is what to do with a subroutine which is called like the parent routine
   call yaml_comment('Test for debug')
   call f_routine(id='PS_Check')
   call f_release_routine()
   call yaml_comment('End test for debug')
   call f_routine(id='PS_Check')
   call f_release_routine()
   call yaml_comment('End test for debug2')
   !call f_malloc_dump_status()

   call f_routine(id='Routine E')
    call f_free(density)
   call f_release_routine()
   call f_routine(id='Routine F')
   call f_release_routine()
   ! call f_malloc_dump_status()

   !Allocations, considering also spin density
!!$   !ionic potential
   pot_ion=f_malloc(0,id='pot_ion')
   !XC potential
   xc_pot=f_malloc(3*2,id='xc_pot')

   !   call f_malloc_dump_status()
   extra_ref=f_malloc_ptr(0,id='extra_ref')

   rhopot=f_malloc(3*2,id='rhopot')

    call f_malloc_dump_status()

   call f_free(rhopot)
!!$
!!$   !   call f_free(density,potential,pot_ion,xc_pot,extra_ref)
!!$!   call f_malloc_dump_status()
!!$   !stop
   call f_free(pot_ion)
   call f_free(potential)
   call f_free(xc_pot)
   !   call f_malloc_dump_status()
   call f_free_ptr(extra_ref)
!!$   !   call yaml_open_map('Last')
!!$   !   call f_malloc_dump_status()
!!$   !   call yaml_close_map()


  !Allocation in an omp section is not permissible
  !We check that.
  call f_err_severe_override(f_err_ignore)
  call yaml_comment('Start OpenMP and f_malloc test',hfill='-')
  call f_err_open_try()
  total=0.d0
  !$omp parallel private(ab,ithread,within_openmp)
  !$omp critical ! (allocate_critical)
  !$ ithread=omp_get_thread_num()
!  ab = f_malloc((/ 10, 10 /),id='ab'//trim(yaml_toa(ithread)))
!  allocate(ab(10,10))
!  ab=1
!  !$omp end critical (allocate_critical)
!total=total+sum(ab)
!deallocate(ab)
!  !$omp critical (deallocate_critical)
!  call f_free(ab)
  !$omp end critical ! (deallocate_critical)
  within_openmp=.false.
!$ within_openmp=omp_in_parallel() .or. omp_get_nested()
!$omp end parallel
  call f_err_close_try()
  if (within_openmp) then
     if (f_err_check()) then
        id = f_get_last_error(message)
        call yaml_map('Error when f_malloc inside OpenMP', message)
     else
        call yaml_map('Error when f_malloc inside OpenMP', .false.)
     end if
  else
     call yaml_map('Error when f_malloc inside OpenMP', .false.)
  end if

!call yaml_map('Total of the allocations',total)
   call f_routine(id='Routine A')
   call f_release_routine()

   call f_routine(id='Routine A')
   call f_release_routine()

   call f_routine(id='Routine A')
   call f_release_routine()


   call f_release_routine()

   call f_malloc_dump_status()

   contains
     
     function dummy_init(n) result(dummy)
       implicit none
       integer, intent(in) :: n
       type(dummy_type) :: dummy
       call init_dummy(n,dummy)
     end function dummy_init
     subroutine init_dummy(n,dummy)
       implicit none
       integer, intent(in) :: n
       type(dummy_type), intent(out) :: dummy
       
       dummy%i=n
       dummy%ptr=f_malloc_ptr(n,id='ptr')
       dummy%ptr=dble(n)
     end subroutine init_dummy
     subroutine free_dummy(dummy)
       implicit none
       type(dummy_type), intent(inout) :: dummy
       call f_free_ptr(dummy%ptr)
       dummy%i=0
     end subroutine free_dummy
end subroutine test_dynamic_memory

subroutine dynmem_sandbox()
  use yaml_output
  use dictionaries, dict_char_len=> max_field_length
  type(dictionary), pointer :: dict2,dictA
  character(len=dict_char_len) :: routinename

  call yaml_comment('Sandbox')  
   !let used to imagine a routine-tree creation
   nullify(dict2)
   call dict_init(dictA)
   dict2=>dictA//'Routine Tree'
!   call yaml_map('Length',dict_len(dict2))
   call add_routine(dict2,'Routine 0')
   call close_routine(dict2,'Routine 0')
   call add_routine(dict2,'Routine A')
   call close_routine(dict2,'Routine A')
   call add_routine(dict2,'Routine B')
   call close_routine(dict2,'Routine B')
   call add_routine(dict2,'Routine C')
   call close_routine(dict2,'Routine C')
   call add_routine(dict2,'Routine D')

   call open_routine(dict2)
   call add_routine(dict2,'SubCase 1')
   call close_routine(dict2,'SubCase 1')

   call add_routine(dict2,'Subcase 2')
   call open_routine(dict2)
   call add_routine(dict2,'SubSubCase1')
   call close_routine(dict2,'SubSubCase1')

   call close_routine(dict2,'SubSubCase1')
   
!   call close_routine(dict2)
   call add_routine(dict2,'SubCase 3')
   call close_routine(dict2,'SubCase 3')
   call close_routine(dict2,'SubCase 3')

   call add_routine(dict2,'Routine E')
   call close_routine(dict2,'Routine E')

   call add_routine(dict2,'Routine F')
!   call yaml_comment('Look Below',hfill='v')

   call yaml_open_map('Test Case before implementation')
   call yaml_dict_dump(dictA)
   call yaml_close_map()
!   call yaml_comment('Look above',hfill='^')

   call dict_free(dictA)

 contains

   subroutine open_routine(dict)
     implicit none
     type(dictionary), pointer :: dict
     !local variables
     integer :: ival
     type(dictionary), pointer :: dict_tmp

     !now imagine that a new routine is created
     ival=dict_len(dict)-1
     routinename=dict//ival

     !call yaml_map('The routine which has to be converted is',trim(routinename))

     call pop(dict,ival)

     dict_tmp=>dict//ival//trim(routinename)
     dict => dict_tmp
     nullify(dict_tmp)

   end subroutine open_routine

   subroutine close_routine(dict,name)
     implicit none
     type(dictionary), pointer :: dict
     character(len=*), intent(in), optional :: name
     !local variables
     logical :: jump_up
     type(dictionary), pointer :: dict_tmp

     if (.not. associated(dict)) stop 'ERROR, routine not associated' 

     !       call yaml_map('Key of the dictionary',trim(dict%data%key))

     if (present(name)) then
        !jump_up=(trim(dict%data%key) /= trim(name))
        jump_up=(trim(routinename) /= trim(name))
     else
        jump_up=.true.
     end if

     !       call yaml_map('Would like to jump up',jump_up)
     if (jump_up) then
        !now the routine has to be closed
        !we should jump at the upper level
        dict_tmp=>dict%parent 
        if (associated(dict_tmp%parent)) then
           nullify(dict)
           !this might be null if we are at the topmost level
           dict=>dict_tmp%parent
        end if
        nullify(dict_tmp)
     end if

     routinename=repeat(' ',len(routinename))
   end subroutine close_routine

   subroutine add_routine(dict,name)
     implicit none
     type(dictionary), pointer :: dict
     character(len=*), intent(in) :: name

     routinename=trim(name)
     call add(dict,trim(name))

   end subroutine add_routine

end subroutine dynmem_sandbox
