subroutine test_dictionaries0()
  use yaml_output
  use dictionaries
  implicit none
  type(dictionary), pointer :: dict1,dict2,dict3
  !local variables
  integer :: ival,nval
  character(len=2) :: val
  type(dictionary), pointer :: list
  !finding operations
!!$  print *,' Filling a linked list' 
!!$  call dict_init(list)
!!$
!!$  call add(list,'x1')
!!$  call add(list,'x2')
!!$  call add(list,'y1')
!!$  call add(list,'z1')
!!$  call add(list,'z1')
!!$ 
!!$  call yaml_dict_dump(list)
!!$
!!$  nval=dict_len(list)
!!$  print *,' Number of elements',dict_len(list)
!!$
!!$  do ival=0,nval-1
!!$     val=list//ival
!!$     print *,'value',ival,valget
!!$  end do
!!$  
!!$  call dict_free(list)

  call yaml_comment('Now we test dictionaries inside yaml.')
  !Test a dictionary
  !alternative way of initializing a dictionary
  !call dict_init(dict1)

  dict1=>dict_new()
  call f_err_open_try()
  ival=dict1//'Toto' 

  call yaml_map('ival not existing, fake value',ival)

  call yaml_map('An error has been raised',f_err_check())
  call yaml_map('Its error id is',f_get_last_error())
  !routine to retrieve the error
  call f_dump_last_error()
  call f_err_close_try()
  call yaml_map('Error pipe is still full',f_err_check())
 
  ! a single scalar
!!  call set(dict1//'',1)
!!$  !can be also set like that, should be avoided
!  call set(dict1,1
!  call yaml_dict_dump(dict1)
!stop
  call set(dict1//'toto',1)
!stop
  call set(dict1//'titi',1.d0)
  call set(dict1//'tutu',(/ '1', '2' /))
  call dict_init(dict2)
  call set(dict2//'a',0)

  !this had  a bug, now solved
  call set(dict1//'List',list_new((/.item. dict2,.item. '4',.item. '1.0'/)))

  !this works
!!$  call add(dict1//'List',dict2)
!!$  call add(dict1//'List',4)
!!$  call add(dict1//'List',1.0)

!!$  !this works too
!!$  list=>dict_new()
!!$  call add(list,dict2)
!!$  call add(list,4)
!!$  call add(list,1.0)
!!$  call set(dict1//'List',list)

!!$  !this also
!!$  list=>dict_new()
!!$  call set(list//'First',dict2)
!!$  call set(list//'Second',4)
!!$  call set(list//'Third',1.0)
!!$  call set(dict1//'List',list)

!!$  dict3=>dict1//'List'
!!$  call yaml_map('Elements of the new dictionary (elems, list)',&
!!$       (/dict3%data%nelems,dict3%data%nitems/))

  dict3=> dict1//'New key'
  call set(dict3//'Example',4)
  call yaml_dict_dump(dict1,flow=.true.)

  !test length functions of dictionaries
  call yaml_map('List length',dict_len(dict1//'List'))
  call yaml_map('Dictionary size',dict_size(dict1))
  call dict_free(dict1)

  !new test, build dictionary on-the-fly
  dict1=>dict_new((/'Key1' .is. 'One',&
       'Key2' .is. 'Two','Key3' .is. 'Three'/))

  call yaml_dict_dump(dict1)
  call dict_free(dict1)
!!$
!!$  !new test, build list on-the-fly
!!$  dict1=list_new((/ .item. 'Val1', .item. 'Val2', .item. 'Val3' ,&
!!$       .item. 'Val4'/))
!!$  call yaml_dict_dump(dict1)
!!$  call dict_free(dict1)

  
end subroutine test_dictionaries0

subroutine test_dictionaries1()
  use yaml_output
  use dictionaries
  implicit none
  !local variables
   integer :: ival,i
   type(dictionary), pointer :: dict2
   type(dictionary), pointer :: dict,dictA
   type(dictionary), pointer :: dictA2,dict_tmp

   !testing add
   call dict_init(dict)
!   call set(dict//0,1)
!   call set(dict//1,2)
!   call set(dict//2,3)
   call add(dict,'1')
   call add(dict,'2')
   call add(dict,'3')
   call yaml_open_map('List')
   call yaml_dict_dump(dict,flow=.true.)
   call yaml_close_map()
!after this call the document has to finish
   call yaml_release_document()

   call yaml_new_document()

   
   call yaml_map('Dictionary length',dict_len(dict))
   call yaml_map('Dictionary size',dict_size(dict))

   call dict_free(dict)
   
   call yaml_comment('Fortran Dictionary Test',hfill='~')

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
   call yaml_map('Size before popping',dict_size(dict2))
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
   call yaml_map('Size after popping',dict_size(dict))
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

   !try to see if extra information can be added after the value
   call set(dictA//'Test Field',6,fmt='(i6.6)')
   ival = dictA//'Test Field'
   call yaml_map('Retrieving Test Field',ival)
   call yaml_map('Retrieving actual value',dict_value(dictA//'Test Field'))


   call set(dictA//'Test Field','6   #extra comment # extra')
   ival = dictA//'Test Field'
   call yaml_map('Retrieving Test Field Again',ival)
   call yaml_map('Retrieving actual value',dict_value(dictA//'Test Field'))
   call yaml_map('Index of comment',index(dict_value(dictA//'Test Field'),'#'))

   call yaml_comment('Prepend dictionary example',hfill='~')

   call yaml_map('Size of dict A',dict_size(dictA))
   call yaml_open_map('Dict A')
   call yaml_dict_dump(dictA)
   call yaml_close_map()

stop
   call dict_init(dict2)
   call set(dict2//'Test1'//'Toto',5)
   call set(dict2//'Test1'//'Titi',6)
   call set(dict2//'Test2'//'Toto',4)
   call set(dict2//'Test2'//'Titi',2)

   call yaml_map('Size of dict 2',dict_size(dict2))
   call yaml_open_map('Dict 2')
   call yaml_dict_dump(dict2)
   call yaml_close_map()

   call yaml_map('Keys of first dict',dict_keys(dictA))
   call yaml_map('Keys of second dict',dict_keys(dict2))


   call prepend(dictA,dict2)
   call yaml_map('Size of prepended',dict_size(dictA))
   call yaml_open_map('Prepended')
   !call yaml_dict_dump2(dictA,verbatim=.true.)
   call yaml_dict_dump(dictA)
   call yaml_close_map()
   
   call yaml_map('Keys of prepended dict',dict_keys(dictA))

   !perform an iterator on dictA
   dict_tmp=>dict_next(dictA)
   do while(associated(dict_tmp))
      call yaml_map('Iterating in dictA',.true.)
      call yaml_map('Key of dictA',dict_key(dict_tmp))
      call yaml_map('Value of dictA',dict_value(dict_tmp))
      dict_tmp=>dict_next(dict_tmp)
   end do

   call dict_free(dictA)

   !fill a list and iterate over it
   dictA=>dict_new()
   do i=1,10
      call add(dictA,'Value'//adjustl(trim(yaml_toa(i))))
   end do

   !perform an iterator on dict
   dict_tmp=>dict_next(dictA)
   do while(associated(dict_tmp))
      call yaml_map('Item of dictA',dict_item(dict_tmp))
      call yaml_map('Value of dictA',dict_value(dict_tmp))
      dict_tmp=>dict_next(dict_tmp)
   end do
   call dict_free(dictA)

   !example which has a bug
   dictA=>list_new((/.item. '5',.item. '6',.item. list_new((/.item.'55',.item. '66'/))/))
!!$!call yaml_open_sequence("",flow=.false.)
!!$call yaml_sequence(advance="no")
!!$call yaml_open_map("SUCCESS",flow=.false.)
!!$call yaml_map("Id","0")
!!$call yaml_map("Message","Operation has succeeded")
!!$call yaml_map("Action","No action")
!!$call yaml_close_map()
!!$call yaml_sequence(advance="no")
!!$call yaml_open_map("GENERIC_ERROR",flow=.false.)
!!$call yaml_map("Id","1")
!!$call yaml_map("Message","UNSPECIFIED")
!!$call yaml_map("Action","UNKNOWN")
!!$call yaml_close_map()
!!$!call yaml_close_sequence()

   !what should be, also this writing has problem in the indentation
!!$    call yaml_sequence('5')
!!$    call yaml_sequence('6')
!!$    call yaml_sequence(advance='no')
!!$    call yaml_open_sequence()
!!$      call yaml_sequence('55')
!!$      call yaml_sequence('66')
!!$    call yaml_close_sequence()
!!$
   call yaml_open_sequence('List in a list')
   call yaml_dict_dump(dictA,verbatim=.true.)
   call yaml_dict_dump(dictA)
   call yaml_dict_dump(dictA,flow=.true.,verbatim=.true.)
   call yaml_dict_dump(dictA,flow=.true.)
   call yaml_close_sequence()

   !perform an iterator on dict
   dict_tmp=>dict_next(dictA)
   do while(associated(dict_tmp))
      call yaml_map('Item of dictA',dict_item(dict_tmp))
      call yaml_map('Key of dictA',dict_key(dict_tmp))
      call yaml_map('Value of dictA',dict_value(dict_tmp))
      dict_tmp=>dict_next(dict_tmp)
   end do
   call dict_free(dictA)

 end subroutine test_dictionaries1
