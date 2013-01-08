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
   implicit none
   integer :: i,l,j,d
   character(len=10), dimension(:), allocatable :: cv
   integer, dimension(:), allocatable :: iv
   real(kind=8), dimension(:), allocatable :: dv
   !First document
   call yaml_new_document()
   call yaml_open_map("Test")
      call yaml_map("Short sentence",.true.)
      call yaml_stream_attributes(iflowlevel=i,ilevel=l,ilast=j,indent=d)
      call yaml_open_map("Foo",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("two",2)
      call yaml_close_map()
      call yaml_stream_attributes(iflowlevel=i,ilevel=l,ilast=j,indent=d)
      call yaml_stream_attributes()
      !call yaml_scalar("1.0")
      call yaml_open_map("toto",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("two",2)
      call yaml_close_map()
      call yaml_stream_attributes(iflowlevel=i,ilevel=l,ilast=j,indent=d)
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
   allocate(dv(0))
   call yaml_new_document()
   !Check calling twice yaml_new_document
   call yaml_new_document()
   call yaml_map('Vector of characters',cv)
   call yaml_map('Vector of integers',iv)
   call yaml_map('Vector of real(kind=8)',dv)
   call yaml_release_document()
   deallocate(cv)
   deallocate(iv)
   deallocate(dv)
end program yaml_test
