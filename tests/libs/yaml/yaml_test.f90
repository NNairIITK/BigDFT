!> @file
!! Test yaml_output module
!! @author
!!    Copyright (C) 2012-2012 BigDFT group
!!    This file is distributed oneder the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test yaml output module
program yaml_test
   use yaml_output
   use dictionaries
   use dynamic_memory
   implicit none
   logical :: fl
   integer :: i,l,j,d,ip,ic
!   integer :: i,l,j,d
   character(len=10), dimension(:), allocatable :: cv
   integer, dimension(:), allocatable :: iv
   real(kind=8), dimension(:), allocatable :: dv
   type(dictionary), pointer :: dict,dictA
   type(dictionary), pointer :: dict2,dictA2,dict_routine,dict_global,dict_array
   real(kind=8), dimension(:), allocatable :: density,rhopot,potential,pot_ion,xc_pot,extra_ref
   integer :: ival

   !First document
   call yaml_new_document()
   call yaml_open_map("Test")
      call yaml_map("Short sentence",.true.)
!      call yaml_stream_attributes(iflowlevel=i,ilevel=l,ilast=j,indent=d,flowrite=fl,indent_previous=ip,icursor=ic)
!      print *,'iflowlevel',i,'ilevel',l,'ilast',j,'indent',d,'flowrite',fl,'indent_previous',ip,'icursor',ic
      call yaml_open_map("Foo",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("two",2)
      call yaml_close_map()
 !     call yaml_stream_attributes(iflowlevel=i,ilevel=l,ilast=j,indent=d,flowrite=fl,indent_previous=ip,icursor=ic)
 !     print *,'iflowlevel',i,'ilevel',l,'ilast',j,'indent',d,'flowrite',fl,'indent_previous',ip,'icursor',ic
!      call yaml_stream_attributes()
!      call yaml_scalar("1.0")
      call yaml_open_map("toto",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("two",2)
      call yaml_close_map()
!      call yaml_stream_attributes(iflowlevel=i,ilevel=l,ilast=j,indent=d)
!      print *,'iflowlevel',i,'ilevel',l,'ilast',j,'indent',d
   call yaml_close_map()
   call yaml_release_document()

   !Second document
   call yaml_new_document()
   call yaml_open_map("Test")
      call yaml_map("I have a very long sentence in order to test if yaml_output fails to print that",.true.)
      call yaml_open_map("Foo",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("two",2)
      call yaml_close_map()
      !call yaml_comment('Bug at this level!: the indentation is not correct')
      !Works if the comment is uncommented!!
      call yaml_open_map("toto",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("two",2)
      call yaml_close_map()
   call yaml_close_map()
   call yaml_release_document()

   !Third document
   allocate(cv(0))
   allocate(iv(0))
!this raises a bug for a vector which is too long
   allocate(dv(11))
dv=3.d0
   call yaml_new_document()
   !Check calling twice yaml_new_document
   call yaml_new_document()
   call yaml_map('Vector of characters',cv)
   call yaml_map('Vector of integers',iv)
call yaml_open_sequence('Vector of double',flow=.true.)
do i=1,size(dv)
   call yaml_sequence(trim(yaml_toa(dv(i),fmt='(1pe12.5)')))
end do
call yaml_close_sequence()
   call yaml_map('Vector of real(kind=8)',dv,fmt='(f3.0)')
   call yaml_release_document()
   deallocate(cv)
   deallocate(iv)
   deallocate(dv)

   print *,'Dictionary test'

   call dict_init(dict)

   !Normal filling of the dictionary
   !this fills a last level
   call set(dict//'Number of Groups',1)

   !this fills a nested level
   call set(dict//'First'//'One',1)
   call set(dict//'First'//'Two',2)

   !alternative way of filling
   dict2 => dict//'First'
   call set(dict//'First'//'Three',3)
   call set(dict2//'Threeb','3b')

   !print dictionary status
   call yaml_dict_dump(dict,flow=.true.)

   !popping a term from the dictionary
   !only a normal pointer can be used
   !try with these examples
   call pop(dict2,'One')
   call pop(dict2,'Two')
!   call pop(dict2,'Three')
   !a further element can be added
   call set(dict//'First'//'Four',4)
   call yaml_open_map('After pop')
   call yaml_dict_dump(dict)
   call yaml_close_map()

   call pop(dict,'First')
   call yaml_open_map('Complete pop')
   call yaml_dict_dump(dict)
   call yaml_close_map()

   !note that we do not have a garbage collector!
   !a call to this will produce a crash due to association above
   !call set(dict2//'Five',5)

   !search for a key and point to it without modifying
   dict2=>find_key(dict,'Number of Gruops')
   call yaml_map('Key found',associated(dict2))
   !the key was wrong, try to find again
   dict2=>find_key(dict,'Number of Groups')
   call yaml_map('Second try, Key found',associated(dict2))
   ival=dict2
   call yaml_map('Value found',ival)
   !increase the value
   call set(dict//'Number of Groups',ival+1)  
   !retrieve it
   ival=dict//'Number of Groups'
   call yaml_map('Alternative way',ival)

  !test if a complete pop will disassociate the dictionry
  call yaml_map('Dictionary associated before last pop',associated(dict))
  call pop(dict,'Number of Groups')
  call yaml_map('Last pop done, still associated',associated(dict))


   call dict_init(dictA)

   call dict_init(dictA2)

   call set(dictA2//'Test1'//'Toto',5)
   call set(dictA2//'Test1'//'Titi',6)

   call set(dictA//'Stack'//0,5)
   call set(dictA//'Stack'//1,4)
   call set(dictA//'Stack'//2,2)
   call set(dictA//'Stack'//3,dictA2)

   call set(dictA//'Stack2',(/'1','2','3'/))

   call yaml_dict_dump(dictA)

   dict2=>find_key(dictA,'Stack')
   call pop(dict2)


   call pop(dict2)

   !  call push(dict2,'Element')
   !  call append(dictA,dictA2)
   call yaml_dict_dump(dictA)
   call dict_free(dictA)
   call yaml_map('Ended',.true.)
   call f_malloc_set_status(memory_limit=0.e0)
   call f_malloc_routine_id('PS_Check')
!!$
!!$   !Allocations, considering also spin density
   !Density
   density=f_malloc(3*2,id='density')
   !Density then potential
   potential=f_malloc(3,id='potential')
   !ionic potential
   pot_ion=f_malloc(3,id='pot_ion')
   !XC potential
   xc_pot=f_malloc(3*2,id='xc_pot')

   call f_malloc_dump_status()
   extra_ref=f_malloc(3,id='extra_ref')

   rhopot=f_malloc(3*2,id='rhopot')
   call f_free(rhopot)
   call yaml_map('Ended1',.true.)
   !call f_free(density,potential,pot_ion,xc_pot,extra_ref)
!!$   call f_malloc_dump_status()
   call f_free(potential)
   call yaml_map('Ended0',.true.)
   call f_free(pot_ion)
   call f_malloc_dump_status()
   call f_free(xc_pot)
   call f_free(extra_ref)
   call yaml_map('Ended2',.true.)
   call f_malloc_finalize()
end program yaml_test
