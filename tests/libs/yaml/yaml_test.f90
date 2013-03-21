!> @file
!! Test yaml_output module
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed oneder the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test yaml output module
program yaml_test
   use yaml_output
   use dictionaries
   implicit none
   !logical :: fl
   integer :: i!,l,j,ip,ic,d
   character(len=10), dimension(:), allocatable :: cv
   integer, dimension(:), allocatable :: iv
   real(kind=8), dimension(:), allocatable :: dv
   type(dictionary), pointer :: dict1,dict2,dict3

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
   !This raises a bug for a vector which is too long
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

   !Fourth document
   allocate(dv(5))
   dv=1.d0
   call yaml_new_document()
   !Check a comment
   call yaml_comment('This document checks the call yaml_comment().')
   call yaml_comment(trim(yaml_toa(dv, fmt='(f14.10)')))
   !Check a very long comment
   call yaml_comment('See if this very long comment is correctly treated:' // &
   & trim(yaml_toa(dv, fmt='(f14.10)')))
   call yaml_open_map('Map')
   call yaml_map('One',1)
   call yaml_comment('No blank characters'//repeat('x',500))
   call yaml_comment(repeat('y',200),hfill='-')
   call yaml_comment(repeat('y',200),tabbing=5,hfill='-')
   call yaml_close_map()
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
   call yaml_release_document()
   deallocate(dv)
   call dict_free(dict1)

end program yaml_test
