!> @file
!! Test the dictionaries of flib
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed oneder the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Routine testing the dictionary object of flib
subroutine test_dictionaries0()
  use yaml_output
  use dictionaries
  implicit none
  type(dictionary), pointer :: dict1,dict2,dict3
  !local variables
  integer :: ival!,nval
!!$  character(len=2) :: val
  character(len=30) :: val2
  type(dictionary), pointer :: dict_tmp
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

  dict1=>dict_new()

  call set(dict1//'hgrid',0.5,fmt='(1pe17.5)')
  call yaml_map('Length and size before',(/dict_len(dict1//'hgrid'),dict_size(dict1//'hgrid')/))
  !call add(dict1//'hgrid','new')
  call set(dict1//'hgrid'//0,'new')

  call yaml_open_map('There was a hidden problem here')
  call yaml_dict_dump(dict1)
  call yaml_close_map()

  call yaml_map('Value of dict1//hgrid',trim(dict_value(dict1//'hgrid')))

  !retrieve value
  val2=dict1//'hgrid' !dict_value(dict1//'hgrid')
  call yaml_map('Value retrieved with equal sign',trim(val2))
  
  !test value of the dictionary, explicitly
  dict_tmp=>dict1//'hgrid'
  call yaml_map('Value explicitly written in the dictionary',trim(dict_tmp%data%value))

  !test length and sizes of the dictionary
  call yaml_map('Length and size after',(/dict_len(dict_tmp),dict_size(dict_tmp)/))

  call dict_free(dict1)

!stop
  dict1=>dict_new()
  call set(dict1//'hgrid',dict_new((/'test1' .is. '1','test2' .is. '2'/)))
  call yaml_map('Length and size before',(/dict_len(dict1//'hgrid'),dict_size(dict1//'hgrid')/))
  call set(dict1//'hgrid'//0,'new')

  call yaml_open_map('Hidden problem here')
  call yaml_dict_dump(dict1)
  call yaml_close_map()

  call yaml_map('Value of dict1//hgrid',trim(dict_value(dict1//'hgrid')))

  !retrieve value
  val2=dict1//'hgrid' !dict_value(dict1//'hgrid')
  call yaml_map('Value retrieved with equal sign',trim(val2))
  
  !test value of the dictionary, explicitly
  dict_tmp=>dict1//'hgrid'
  call yaml_map('Verify that the child is still associated',associated(dict_tmp%child))

  !test length and sizes of the dictionary
  call yaml_map('Length and size after',(/dict_len(dict_tmp),dict_size(dict_tmp)/))

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
   integer :: ival,i,j
   type(dictionary), pointer :: dict2
   type(dictionary), pointer :: dict,dictA
   type(dictionary), pointer :: dictA2,dict_tmp,zero1,zero2
   double precision, dimension(3) :: tmp_arr

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

   !retrieve the value from the Stack2 key
   tmp_arr=dictA//'Stack2'

   call yaml_map('Values retrieved from the dict',tmp_arr,fmt='(1pg12.5)')

   dict2=>find_key(dictA,'Stack')
   call pop(dict2)


   call pop(dict2)

   !  call push(dict2,'Element')
   !  call append(dictA,dictA2)
   call yaml_dict_dump(dictA)

   !retrieve the value from the Stack key
   tmp_arr(1:2)=dictA//'Stack'
   call yaml_map('Two values from Stack key',tmp_arr,fmt='(1pg12.5)')

   !retrieve the value from the a scalar
   tmp_arr=dictA//'Stack'//0
   call yaml_map('Array filled with a scalar',tmp_arr,fmt='(1pg12.5)')

!!$   !try to see if extra information can be added after the value
!!$   call set(dictA//'Test Field',6,fmt='(i6.6)')
!!$   ival = dictA//'Test Field'
!!$   call yaml_map('Retrieving Test Field',ival)
!!$   call yaml_map('Retrieving actual value',dict_value(dictA//'Test Field'))
!!$
!!$
!!$   call set(dictA//'Test Field','6   #extra comment # extra')
!!$   ival = dictA//'Test Field'
!!$   call yaml_map('Retrieving Test Field Again',ival)
!!$   call yaml_map('Retrieving actual value',dict_value(dictA//'Test Field'))
!!$   call yaml_map('Index of comment',index(dict_value(dictA//'Test Field'),'#'))

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

   !verify the euqlity between dictionaries
   call yaml_map('Test1 and and Test2 are equal',dict2//'Test1' == dict2//'Test2')

   !now correct
   call set(dict2//'Test1'//'Toto',4)
   call set(dict2//'Test1'//'Titi',2)

   call yaml_map('Corrected version',dict2)
   
   !verify the equality between dictionaries
   call yaml_map('Test1 and and Test2 are equal',dict2//'Test1' == dict2//'Test2')

   !now add another element, written differently
   call set(dict2//'Test1'//'Tutu',4.d0,fmt='(1pe12.5)')
   call set(dict2//'Test2'//'Tutu','4.d0')

   call yaml_map('Added version',dict2)

   !verify the equality between dictionaries
   call yaml_map('Test1 and and Test2 are equal',dict2//'Test1' == dict2//'Test2')

   !also test the possibility two arrays filled with zeroes
   zero1=>list_new(.item. list_new(.item. '0.0000000000000000',&
        .item. '0.0000000000000000',.item. '0.0000000000000000'))
   zero2=>list_new(.item. list_new(.item. '0.',.item. '0.',.item. '0.'))

   call yaml_map('List of list of zeroes, first version',zero1)
   call yaml_map('List of list of zeroes, second version',zero2)

      !verify the equality between dictionaries
   call yaml_map('Zero1 and and Zero2 are equal',zero1==zero2)

   call dict_free(zero1)
   call dict_free(zero2)

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
   dict_tmp=>dict_iter(dictA)
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
   dict_tmp => list_new((/.item.'55',.item. '66'/))
   dictA=>list_new((/.item. '5',.item. '6',.item. dict_tmp/))
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
   call yaml_dict_dump(dictA,flow=.false.)
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

!!$   !try to steel a argument (does not work, should arrange routine set to be full-proof)
!!$   !fill a list and iterate over it
!!$   dictA=>dict_new()
!!$   do i=1,10
!!$      call add(dictA,trim(yaml_toa((/ (j,j=i,i+3) /))))
!!$   end do
!!$
!!$   call yaml_map('List before',dictA)
!!$
!!$   dict_tmp=>dict_new('ciao' .is. '1','hello' .is. '2')
!!$   dictA2=>dictA//3
!!$   call set(dict_tmp//'bonjour',dictA2)
!!$
!!$   call yaml_map('Thief dict',dict_tmp)
!!$
!!$   call yaml_map('List after',dictA)
!!$   call dict_free(dictA)
!!$   call dict_free(dict_tmp)



 end subroutine test_dictionaries1

 subroutine test_copy_merge()
   use dictionaries
   use yaml_output
   implicit none

   type(dictionary), pointer :: dict, cpy, subd

   dict => dict_new(&
         & "__comment__" .is. 'Grid shifts', &
         & "__cond__"    .is. dict_new("__master_key__" .is. "kpt_method", "__when__" .is. list_new( .item. "MPGrid")), &
         & "__default__" .is. list_new( .item."0.", .item."0.", .item."0.") )

   call yaml_open_map("test dict_copy")
   call yaml_open_map("original")
   call yaml_dict_dump(dict)
   call yaml_close_map()
   nullify(cpy)
   call dict_copy(cpy, dict)
   call yaml_open_map("copy")
   call yaml_dict_dump(cpy)
   call yaml_close_map()
   call dict_free(cpy)
   call yaml_close_map()

   subd => dict_new(  &
         & "__exclusive__" .is. dict_new( "123" .is. "operation 123", &
         &                                  "456" .is. "operation 456" ), &
         & "__default__"   .is. list_new(.item."1.", .item."2.", .item."3." ) )
   call yaml_open_map("test dict_update")
   call dict_update(dict, subd)
   call yaml_open_map("additional")
   call yaml_dict_dump(subd)
   call yaml_close_map()
   call yaml_open_map("after merge")
   call yaml_dict_dump(dict)
   call yaml_close_map()
   call yaml_close_map()
   call dict_free(subd)

   call dict_free(dict)
 end subroutine test_copy_merge

subroutine test_dictionary_for_atoms()
  use yaml_output
  implicit none

  character(len = 50) :: gu,fmts
  double precision, dimension(3) :: cell, xred, hgrids
  double precision :: tt


  call yaml_open_map("Atomic structure")

  cell = 20.345752999999998
  call yaml_map('Cell', cell)

  hgrids = cell / (/ 54, 40, 40 /)

  call yaml_open_sequence('Positions')

  call yaml_sequence(advance='no')
  xred = (/ 0.2516085125D-05,  0.5826606155D-05,  20.34574212d0 /)
  call print_one_atom('Si',xred,hgrids,1)

  call yaml_sequence(advance='no')
  xred = (/ 5.094032326d0,  5.153107111d0,  0.3047989908d-01 /)
  call print_one_atom('Si',xred,hgrids,2)
!!$  call yaml_map("Si", xred, fmt="(g18.10)", advance = "no")
!!$  xred = xred / hgrids
!!$  write(gu, "('[ 'F6.2', 'F6.2', 'F6.2'] 'I4.4)") xred, 2
!!$  call yaml_comment(gu)

  call yaml_sequence(advance='no')
  xred = (/ 0.3049344014d-01,  5.153107972d0,  5.094018600d0 /)
  call print_one_atom('Si',xred,hgrids,3)
!!$  call yaml_map("Si", xred, fmt="(g18.10)", advance = "no")
!!$  xred = xred / hgrids
!!$  write(gu, "('[ 'F6.2', 'F6.2', 'F6.2'] 'I4.4)") xred, 3
!!$  call yaml_comment(gu)

  call yaml_close_sequence()

  call yaml_close_map()

  !now print some double precision values to understand which is the best format
  tt=real(0.5e0,kind=8) !use a conversion from float

  call yaml_map('Real without format',tt)
  fmts(1:len(fmts))='(1pe25.17)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(1pe24.16)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(es23.16)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(es24.17)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(es25.18)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(es26.19)'
  call yaml_map('Real with format '//trim(fmts),tt+epsilon(1.d0),fmt=fmts)
  fmts(1:len(fmts))='(es27.20)'
  call yaml_map('Real with format '//trim(fmts),tt-epsilon(1.d0),fmt=fmts)
  fmts(1:len(fmts))='(es26.19)'
  call yaml_map('Real with format '//trim(fmts),epsilon(1.d0),fmt=fmts)
  fmts(1:len(fmts))='(es27.20)'
  call yaml_map('Real with format '//trim(fmts),-epsilon(1.d0),fmt=fmts)


  contains

    subroutine print_one_atom(atomname,rxyz,hgrids,id)
      implicit none
      integer, intent(in) :: id
      character(len=*), intent(in) :: atomname
      double precision, dimension(3), intent(in) :: rxyz,hgrids
      !local variables
      character(len=*), parameter :: fmtat='(g18.10)',fmtg='(F6.2)'
      integer :: i

      call yaml_open_sequence(atomname,flow=.true.)
      do i=1,3
         call yaml_sequence(yaml_toa(rxyz(i),fmt=fmtat))
      end do
      call yaml_close_sequence(advance='no')
      call yaml_comment(trim(yaml_toa(rxyz/hgrids,fmt=fmtg))//trim(yaml_toa(id))) !we can also put tabbing=

    end subroutine print_one_atom

end subroutine test_dictionary_for_atoms
