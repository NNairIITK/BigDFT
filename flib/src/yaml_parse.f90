!> @file
!! Module to parse the yaml (flib library)
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module containing routines to parse a yaml input
module yaml_parse
  implicit none

  private

  integer :: STREAM_START, STREAM_END
  integer :: DOCUMENT_START, DOCUMENT_END
  integer :: SEQUENCE_START, SEQUENCE_END
  integer :: MAPPING_START, MAPPING_END
  integer :: ALIAS, SCALAR, ERROR

  integer, public :: YAML_PARSE_ERROR       = 0
  integer, public :: YAML_PARSE_UNSUPPORTED = 0

  public :: yaml_parse_from_file
  public :: yaml_parse_from_char_array
  public :: yaml_parse_from_string

  !for internal f_lib usage
  public :: yaml_parse_errors

contains

  subroutine yaml_parse_errors()
    use dictionaries, only: f_err_define
    implicit none

    call f_err_define(err_name='YAML_PARSE_ERROR',&
         err_msg='YAML parse error.',&
         err_action='modify your inputs.',&
         err_id=YAML_PARSE_ERROR)
    call f_err_define(err_name='YAML_PARSE_UNSUPPORTED',&
         err_msg='YAML standard not supported.',&
         err_action='kindly ask developers to finish implementation.',&
         err_id=YAML_PARSE_UNSUPPORTED)

  end subroutine yaml_parse_errors

  subroutine yaml_parse_from_file(dict, fname)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: fname
    
    integer(kind = 8) :: parser

    call yaml_parser_c_init(parser, fname, len(fname))
    dict => yaml_parse_(parser)
  end subroutine yaml_parse_from_file

  subroutine yaml_parse_from_char_array(dict, carr)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character, dimension(:), intent(in) :: carr
    
    integer(kind = 8) :: parser

    call yaml_parser_c_init_from_buf(parser, carr(1), size(carr))
    dict => yaml_parse_(parser)
  end subroutine yaml_parse_from_char_array

  subroutine yaml_parse_from_string(dict, str)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: str
    
    integer(kind = 8) :: parser

    call yaml_parser_c_init_from_buf(parser, str, len_trim(str))
    dict => yaml_parse_(parser)
  end subroutine yaml_parse_from_string

  function yaml_parse_(parser) result(output)
    use dictionaries
    implicit none
    integer(kind = 8), intent(in) :: parser

    type(dictionary), pointer :: dict, doc, output
    integer :: event, errid
    character(max_field_length) :: val

    ! Get event values from C.
    call yaml_parser_c_get_stream_start(STREAM_START)
    call yaml_parser_c_get_stream_end(STREAM_END)
    call yaml_parser_c_get_document_start(DOCUMENT_START)
    call yaml_parser_c_get_document_end(DOCUMENT_END)
    call yaml_parser_c_get_sequence_start(SEQUENCE_START)
    call yaml_parser_c_get_sequence_end(SEQUENCE_END)
    call yaml_parser_c_get_mapping_start(MAPPING_START)
    call yaml_parser_c_get_mapping_end(MAPPING_END)
    call yaml_parser_c_get_alias(ALIAS)
    call yaml_parser_c_get_scalar(SCALAR)
    call yaml_parser_c_get_error(ERROR)

    ! Initialise error if required.
    if (YAML_PARSE_ERROR == 0) then
       call f_err_define(err_name='YAML_PARSE_ERROR',&
            err_msg='YAML parse error.',&
            err_action='modify your inputs.',&
            err_id=YAML_PARSE_ERROR)
    end if
    if (YAML_PARSE_UNSUPPORTED == 0) then
       call f_err_define(err_name='YAML_PARSE_UNSUPPORTED',&
            err_msg='YAML standard not supported.',&
            err_action='kindly ask developers to finish implementation.',&
            err_id=YAML_PARSE_UNSUPPORTED)
    end if

    call f_err_open_try()
    call dict_init(dict)
    event = 0
    nullify(doc)
    val=repeat(' ',len(val))
    do while (event /= STREAM_END)
       call yaml_parser_c_next(parser, event, val, max_field_length)
       !print *,'event',event_toa(event),event,trim(val),'end'
       if (event == ERROR) then
          call f_err_throw(err_id = YAML_PARSE_ERROR, err_msg = trim(val))
          exit
       end if

       if (event == DOCUMENT_END) then
          if (.not.associated(doc)) call dict_init(doc) ! empty document case
          call add(dict, doc)
          nullify(doc)
       else if (event == MAPPING_START) then
          doc => build_map(parser)                      ! dictionary document case
       else if (event == SEQUENCE_START) then
          doc => build_seq(parser)                      ! list document case
       else if (event == SCALAR) then
          call dict_init(doc)                           ! scalar document case
          call set(doc, val)
       end if

       if (f_err_check()) exit
    end do
    errid = 0
    if (f_err_check()) then
       if (associated(doc)) call add(dict, doc)
       errid = f_get_last_error(val)
    end if
    call f_err_close_try()

    if (event /= STREAM_END) call yaml_parser_c_finalize(parser)

    output => dict

    if (errid /= 0) call f_err_throw(err_id = errid, err_msg = val)
  contains
    !>determine which is the event that has been recognized, to be used mostly for debugging purposes
    function event_toa(event) result(toa)
      implicit none
      integer, intent(in) :: event
      character(len=32) :: toa
      if(event==STREAM_START) then
         toa(1:len(toa))='STREAM_START'
      else if(event==STREAM_END) then
         toa(1:len(toa))='STREAM_END'
      else if(event==DOCUMENT_START) then
         toa(1:len(toa))='DOCUMENT_START'
      else if(event==DOCUMENT_END) then
         toa(1:len(toa))='DOCUMENT_END'
      else if(event==SEQUENCE_START) then
         toa(1:len(toa))='SEQUENCE_START'
      else if(event==SEQUENCE_END) then
         toa(1:len(toa))='SEQUENCE_END'
      else if(event==MAPPING_START) then
         toa(1:len(toa))='MAPPING_START'
      else if(event==MAPPING_END) then
         toa(1:len(toa))='MAPPING_END'
      else if(event==ALIAS) then
         toa(1:len(toa))='ALIAS'
      else if(event==SCALAR) then
         toa(1:len(toa))='SCALAR'
      else if(event==ERROR) then
         toa(1:len(toa))='ERROR'
      else
         toa(1:len(toa))='UNKNOWN'
      end if
    end function event_toa
  end function yaml_parse_

  recursive function build_map(parser) result(map)
    use dictionaries
    !use yaml_output
    implicit none
    integer(kind = 8), intent(in) :: parser

    type(dictionary), pointer :: m, sub, map
    integer :: event
    character(max_field_length) :: val, key

    call dict_init(m)
    map => m

    event = 0
    key(1:max_field_length) = " "
    do while (event /= STREAM_END)
       call yaml_parser_c_next(parser, event, val, max_field_length)
       !write(*,*) "map", event

       if (event == ERROR) then
          call f_err_throw(err_id = YAML_PARSE_ERROR, err_msg = trim(val))
          return
       end if

       if (event == MAPPING_END) then
          exit
       else if (event == MAPPING_START) then
          sub => build_map(parser)
          if (len_trim(key) == 0) stop "no key"
          call set(m // key, sub)
          key(1:max_field_length) = " "
       else if (event == SEQUENCE_START) then
          sub => build_seq(parser)
          if (len_trim(key) == 0) stop "no key"
          call set(m // key, sub)
          key(1:max_field_length) = " "
       else if (event == SCALAR) then
          if (len_trim(key) > 0) then
             ! This is a simple key / val entry.
             call set(m // key, val)
             key(1:max_field_length) = " "
          else
             ! We store a key for later usage.
             key = val
             !write(*,*) "set ", key
          end if
       else if (event == ALIAS) then
          call f_err_throw(err_id = YAML_PARSE_UNSUPPORTED, &
               & err_msg = "unsupported alias to " // trim(val))
          return
       end if
       
       if (f_err_check()) return
       
    end do
  end function build_map

  recursive function build_seq(parser) result(seq)
    use dictionaries
    implicit none
    integer(kind = 8), intent(in) :: parser

    type(dictionary), pointer :: s, sub, seq
    integer :: event
    character(max_field_length) :: val

    call dict_init(s)
    seq => s

    event = 0
    do while (event /= STREAM_END)
       call yaml_parser_c_next(parser, event, val, max_field_length)
       !write(*,*) "seq", event

       if (event == ERROR) then
          call f_err_throw(err_id = YAML_PARSE_ERROR, err_msg = trim(val))
          return
       end if

       if (event == SEQUENCE_END) then
          exit
       else if (event == MAPPING_START) then
          sub => build_map(parser)
          call add(s, sub)
       else if (event == SEQUENCE_START) then
          sub => build_seq(parser)
          call add(s, sub)
       else if (event == SCALAR) then
          call add(s, val)
       else if (event == ALIAS) then
          call f_err_throw(err_id = YAML_PARSE_UNSUPPORTED, &
               & err_msg = "unsupported alias to " // trim(val))
          return
       end if
       
       if (f_err_check()) return
       
    end do

  end function build_seq

end module yaml_parse
