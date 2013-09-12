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
   use dictionaries, dict_char_len=> max_field_length
   use dynamic_memory
   implicit none
   !logical :: fl
   !First document
  
   call yaml_new_document()
   call yaml_comment('Yaml Output Module Test',hfill='~')
   call test_yaml_output1()
   call yaml_release_document()

   !Second document
   call yaml_new_document()
    call test_yaml_output2()
   call yaml_release_document()
   !Third document
   call yaml_new_document()
   !Check calling twice yaml_new_document
   call yaml_new_document()
    call test_yaml_output_sequences1()
   call yaml_release_document()

   !Fourth document
   call yaml_new_document()
    call test_yaml_output_sequences2()
   call yaml_release_document()

   !Fourth-A document
   call yaml_new_document()
   call yaml_invoice_example()
   call yaml_release_document()

   !Fourth-A2 document
   call yaml_new_document()
   call yaml_invoice_example_with_dictionaries()
   call yaml_release_document()

   !Fourth-B document, to be moved
   call yaml_new_document()
    call test_error_handling()
   call yaml_release_document()

   !Fifth document
   call yaml_new_document()
    call test_dictionaries0()
   call yaml_release_document()

   call yaml_new_document()
   call test_dictionaries1()
   call yaml_release_document()

   call yaml_new_document()
    call test_dynamic_memory()
   call yaml_release_document()

   call yaml_new_document()
    call test_copy_merge()
   call yaml_release_document()

   !prepare the finalization of the library
   call f_lib_finalize()
end program yaml_test
