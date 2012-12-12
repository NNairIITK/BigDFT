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
   call yaml_new_document()
   call yaml_open_map("Test")
   call yaml_map("I have a very long sentence in order to test if yaml_output fails to print that",.true.)
      call yaml_open_map("Foo",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("deux",2)
      call yaml_close_map()
      call yaml_open_map("toto",flow=.true.)
      call yaml_map("one",1)
      call yaml_map("deux",2)
      call yaml_close_map()
   call yaml_close_map()
   call yaml_release_document()
end program yaml_test
