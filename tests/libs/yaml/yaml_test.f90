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
   logical :: fl
   integer :: i,l,j,d,ip,ic
   character(len=10), dimension(:), allocatable :: cv
   integer, dimension(:), allocatable :: iv
   real(kind=8), dimension(:), allocatable :: dv
   !First document
   call yaml_new_document()
   call yaml_open_map('Test 1')
      call yaml_map('Short sentence',.true.)
      call yaml_stream_attributes(iflowlevel=i,ilevel=l,ilast=j,indent=d,flowrite=fl,indent_previous=ip,icursor=ic)
      call yaml_open_map('Foo',flow=.true.)
      call yaml_map('one',1)
      call yaml_map('two',2)
      call yaml_close_map()
      call yaml_open_map('toto',flow=.true.)
      call yaml_map('one',1)
      call yaml_map('two',2)
      call yaml_close_map()
   call yaml_close_map()
   call yaml_release_document()

   !Second document
   call yaml_new_document()
   call yaml_open_map('Test 2')
      call yaml_map('I have a very long sentence in order to test if yaml_output fails to print that',.false.)
      call yaml_open_map('Foo',flow=.true.)
      call yaml_map('one',1)
      call yaml_map('two',2)
      call yaml_close_map()
      call yaml_open_map('toto',flow=.true.)
      call yaml_map('one',1)
      call yaml_map('two',2)
      call yaml_close_map()
   call yaml_close_map()
   call yaml_release_document()

   !Third document
   allocate(cv(0))
   allocate(iv(0))
   allocate(dv(0))
   call yaml_new_document()
   call yaml_comment('Test 3')
   call yaml_comment('Check calling twice yaml_new_document')
   call yaml_new_document()
   call yaml_map('Vector of characters',cv)
   call yaml_map('Vector of integers',iv)
   call yaml_map('Vector of real(kind=8)',dv)
   call yaml_close_sequence()
   deallocate(cv)
   allocate(cv(3))
   do i=1,3
      cv(i)=yaml_toa(i)
   end do
   call yaml_comment('Check flow=.true.',hfill='-')
   call yaml_open_map('Open a sequence')
      call yaml_open_sequence(flow=.true.)
         call yaml_sequence(trim(yaml_toa(cv)))
         call yaml_sequence(trim(yaml_toa(cv)))
      call yaml_close_sequence()
   call yaml_close_map()
   call yaml_open_map('Open a sequence',flow=.true.)
      call yaml_open_sequence('Key')
         call yaml_sequence(trim(yaml_toa(cv)))
         call yaml_sequence(trim(yaml_toa(cv)))
      call yaml_close_sequence()
   call yaml_close_map()
   call yaml_open_sequence('Open a sequence',flow=.true.)
      call yaml_sequence(trim(yaml_toa(cv)))
      call yaml_sequence(trim(yaml_toa(cv)))
   call yaml_close_sequence()
   call yaml_release_document()
   deallocate(cv)
   deallocate(iv)
   deallocate(dv)
   call yaml_release_document()

   !Fourth document
   call yaml_new_document()
   call yaml_comment('This a document with a sequence and not a map')
   call yaml_sequence('One')
   call yaml_sequence(trim(yaml_toa( (/ 1, 2, 3 /))))
   call yaml_sequence(advance='no')
      call yaml_open_map(flow=.true.)
         call yaml_map('Example','Test')
         call yaml_map('Two','Test')
      call yaml_close_map()
   call yaml_sequence(advance='no')
      call yaml_open_sequence(flow=.true.)
         call yaml_sequence(trim(yaml_toa(1)))
         call yaml_sequence('Two')
         call yaml_sequence(advance='no')
            call yaml_scalar('1.0',advance='no')
      call yaml_close_sequence()
   call yaml_release_document()
end program yaml_test
