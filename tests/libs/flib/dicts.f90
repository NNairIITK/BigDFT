subroutine test_dictionaries0()
  use yaml_output
  use dictionaries
  implicit none
  type(dictionary), pointer :: dict1,dict2,dict3
  !local variables
  

  call yaml_comment('Now we test dictionaries inside yaml.')
  !Test a dictionary
  call dict_init(dict1)
  call set(dict1//'toto',1)
  call set(dict1//'titi',1.d0)
  call set(dict1//'tutu',(/ '1', '2' /))
  call dict_init(dict2)
  call set(dict2//'a',0)
  !call set(dict1//'dict2',dict2)
  call set(dict1//'List'//0,dict2)
  call set(dict1//'List'//1,4)
  call set(dict1//'List'//2,1.0)
  dict3=> dict1//'New key'
  call set(dict3//'Example',4)
  call yaml_dict_dump(dict1,flow=.true.)

  !test length functions of dictianries
  call yaml_map('List length',dict_len(dict1//'List'))
  call yaml_map('Dictionary size',dict_size(dict1))
  call dict_free(dict1)

end subroutine test_dictionaries0

subroutine test_dictionaries1()
  use yaml_output
  use dictionaries
  implicit none
  !local variables
   integer :: ival
   type(dictionary), pointer :: dict2
   type(dictionary), pointer :: dict,dictA
   type(dictionary), pointer :: dictA2

   !testing add
   call dict_init(dict)
!   call set(dict//0,1)
!   call set(dict//1,2)
!   call set(dict//2,3)
   call add(dict,'1')
   call add(dict,'2')
   call add(dict,'3')
   call yaml_dict_dump(dict,flow=.true.)
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

   call yaml_comment('Prepend dictionary example',hfill='~')

   call yaml_map('Size of dict A',dict_size(dictA))
   call yaml_open_map('Dict A')
   call yaml_dict_dump(dictA)
   call yaml_close_map()
   call dict_init(dict2)
   call set(dict2//'Test1'//'Toto',5)
   call set(dict2//'Test1'//'Titi',6)
   call set(dict2//'Test2'//'Toto',4)
   call set(dict2//'Test2'//'Titi',2)

   call yaml_map('Size of dict 2',dict_size(dict2))
   call yaml_open_map('Dict 2')
   call yaml_dict_dump(dict2)
   call yaml_close_map()



   call prepend(dictA,dict2)
   call yaml_map('Size of prepended',dict_size(dictA))
   call yaml_open_map('Prepended')
   call yaml_dict_dump(dictA)
   call yaml_close_map()
   
   call dict_free(dictA)

 end subroutine test_dictionaries1
