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
   use yaml_parse
   implicit none
   type(dictionary), pointer :: dict_tmp
   !type(yaml_cl_parse) :: parser
   !logical :: fl

   call f_lib_initialize()
!!$   call yaml_comment('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
!!$   call yaml_comment('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')
!!$   call yaml_comment('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
!!$   call yaml_comment('DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD')

!!$   parser=yaml_cl_parse_null()
!!$   !set valid options
!!$   call yaml_cl_parse_option(parser,'test','1',&
!!$        'this is a valid test option','t',&
!!$        dict_new('Test' .is. 'long help'))
!!$
!!$   !verify the parsing
!!$   call yaml_cl_parse_cmd_line(parser)
!!$
!!$   call yaml_map('Parsed options',parser%options)
!!$   call yaml_map('Parsed info',parser%args)
!!$
!!$   call yaml_cl_parse_free(parser)


   !call profile_dictionary_usage()
!!$   call f_lib_finalize()
!!$   stop

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

   !Sixth document: test dictionaries
   call yaml_new_document()
   call test_dictionaries1()
   call yaml_release_document()

   !Seventh document: Test dynamic memory allocation
   call yaml_new_document()
   call test_dynamic_memory()
   call yaml_release_document()

   call f_malloc_dump_status(dict_summary=dict_tmp)
   call yaml_map('Summary',dict_tmp)
   call dict_free(dict_tmp)
!   call f_lib_finalize()
!stop

   call yaml_new_document()
   call test_copy_merge()
   call yaml_release_document()

   !test the yaml parsing
   call yaml_parse_file_and_string()

   call yaml_new_document()
   call test_dictionary_for_atoms()
   call yaml_release_document()

   call profile_dictionary_usage()

   !prepare the finalization of the library
   call f_lib_finalize()

end program yaml_test

subroutine yaml_parse_file_and_string()
  use dictionaries
  use yaml_parse
  use yaml_output
  implicit none
  type(dictionary), pointer :: dict_parse
  character(len=*), parameter :: cr=char(13)//char(10) !carriage return
  character(len=*), parameter ::stream=&
       "---"//cr//&
       "Key1: field1"//cr//&
       "Key2: "//cr//&
       " Subkey2-1: ciao"//cr//&
       " Subkey2-2: [0, hello, 2]"//cr//&
       "Key3:"//cr//&
       " - One"//cr//&
       " - Two"//cr//&
       " - Three"//cr//&
       "---"//cr//&
       "#example of the input variables"//cr//&
       "inputvars: #name of the variable as declared in the code"//cr//&
       " COMMENT: 'This is the description of the variable as will appear in the logfile'"//cr//&
       " RANGE: ['from','to'] #always put two numbers (also .inf) can be put"//cr//&
       " EXCLUSIVE: #here follows a list of allowed values (either RANGE or EXCLUSIVE)"//cr//&
       "  - Value1:  'comments of value1'"//cr//&
       "  - Value2:  'comment of value2'"//cr//&
       " CONDITION: #here the conditions for which the variable makes sense are written "//cr//&
       "   MASTER_KEY: foo #this means that inputvars makes sense only if foo is specified "//cr//&
       "   WHEN: #provide a list of allowed values of foo for which inputvars is meaningful "//cr//&
       "     - fooval1 "//cr//&
       "     - fooval2   "//cr//&
       "#then the profiles follows, which gives to the variables the allowed name "//cr//&
       " default: 'value of the default, written as a string' "//cr//&
       " profile1: 'value1' # if the user specifies inputvars: profile1 then inputvars will be value1 "//cr//&
       " profile2: 'value2'"//cr

  
  call yaml_comment('Yaml parsing',hfill='-')

!  call yaml_parse_from_file(dict_parse,'testfile.yaml')
!  call yaml_dict_dump_all(dict_parse)
!  call dict_free(dict_parse)

  call yaml_parse_from_string(dict_parse,stream)
  call yaml_dict_dump_all(dict_parse)
  call dict_free(dict_parse)

  
end subroutine yaml_parse_file_and_string

  subroutine help_screen()
    write(*,*)' Usage of the command line instruction'
    write(*,*)' --taskgroup-size=<mpi_groupsize>'
    write(*,*)' --runs-file=<list_posinp filename>'
    write(*,*)' --run-id=<name of the run>: it can be also specified as unique argument'
    write(*,*)' --help : prints this help screen'
  end subroutine help_screen

