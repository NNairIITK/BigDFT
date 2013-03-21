!> @file
!! Define the modules (yaml_strings and yaml_output) and the methods to write yaml output
!! yaml: Yet Another Markeup Language (ML for Human)
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module yaml_output
  use yaml_strings
  use dictionaries
  implicit none
  private 

  !> Yaml events for dump routine
  integer, parameter :: NONE           = -1000
  integer, parameter :: DOCUMENT_START = -1001
  integer, parameter :: DOCUMENT_END   = -1002
  integer, parameter :: MAPPING_START  = -1003
  integer, parameter :: MAPPING_END    = -1004
  integer, parameter :: SEQUENCE_START = -1005
  integer, parameter :: SEQUENCE_END   = -1006
  integer, parameter :: SCALAR         = -1007
  integer, parameter :: COMMENT        = -1008
  integer, parameter :: MAPPING        = -1009
  integer, parameter :: SEQUENCE_ELEM  = -1010
  integer, parameter :: NEWLINE        = -1011
  integer, parameter :: COMMA_TO_BE_PUT= 10
  integer, parameter :: STREAM_ALREADY_PRESENT=-1

  integer, parameter :: tot_max_record_length=95   !< Max record length by default
  integer, parameter :: tot_max_flow_events=500    !< Max flow events
  integer, parameter :: tot_streams=10             !< Max total number of streams
  integer, parameter :: tab=5

  integer :: active_streams=0  !< Number of active streams (stdout always active after init)
  integer :: default_stream=1  !< Id of the default stream

  !parameter of the document
  type :: yaml_stream
     logical :: document_closed=.true.  !< Put the starting of the document if new_document is called
     logical :: pp_allowed=.true.       !< Pretty printing allowed
     integer :: unit=6                  !< Unit for the stdout
     integer :: max_record_length=tot_max_record_length
     logical :: flowrite=.false.        !< Write in flow (.false.=no .true.=yes)
     integer :: Wall=-1                 !< Warning messages of level Wall stop the program (-1 : none)
     integer :: indent=1                !< Blank spaces indentations for Yaml output level identification
     integer :: indent_previous=0       !< Indent level prior to flow writing
     integer :: indent_step=2           !< Indentation level
     integer :: tabref=40               !< position of tabular in scalar assignment (single column output)
     integer :: icursor=1               !< running position of the cursor on the line
     integer :: itab_active=0           !< number of active tabbings for the line in flowrite
     integer :: itab=0                  !< Tabbing to have a look on
     integer :: ilevel=0                !< Number of opened levels
     integer :: iflowlevel=0            !< Levels of flowrite simoultaneously enabled
     integer :: ilast=0                 !< Last level with flow==.false.
     integer :: icommentline=0          !< Active if the line being written is a comment
     integer, dimension(tot_max_record_length/tab) :: linetab=0   !< Value of the tabbing in the line
     integer :: ievt_flow=0                                       !< Events which track is kept of in the flowrite
     integer, dimension(tot_max_flow_events) :: flow_events=0     !< Set of events in the flow
     type(dictionary), pointer :: dict_warning=>null()            !< dictionary of warnings emitted in the stream
  end type yaml_stream

  type(yaml_stream), dimension(tot_streams), save :: streams    !< Private array containing the streams
  integer, dimension(tot_streams) :: stream_units=6 !default units unless otherwise specified               
                                                  
  interface yaml_map                              
     module procedure yaml_map,yaml_map_i,yaml_map_f,yaml_map_d,yaml_map_l,yaml_map_iv,yaml_map_dv,yaml_map_cv
  end interface                                   
                                                  
  public :: yaml_new_document,yaml_release_document
  public :: yaml_map,yaml_open_map,yaml_close_map 
  public :: yaml_sequence,yaml_open_sequence,yaml_close_sequence
  public :: yaml_comment,yaml_warning,yaml_toa,yaml_newline
  public :: yaml_set_stream,yaml_get_default_stream,yaml_set_default_stream,yaml_stream_attributes
  public :: yaml_close_all_streams
  public :: yaml_date_and_time_toa,yaml_scalar,yaml_date_toa,yaml_dict_dump

contains                                          
          
  !> Initialize the stream to default values
  function stream_null() result(strm)
    implicit none
    type(yaml_stream) :: strm

    strm%document_closed=.true.  
    strm%pp_allowed=.true.
    strm%unit=6
    strm%max_record_length=tot_max_record_length
    strm%flowrite=.false.        
    strm%Wall=-1                 
    strm%indent=1                
    strm%indent_previous=0       
    strm%indent_step=2           
    strm%tabref=40               
    strm%icursor=1               
    strm%itab_active=0           
    strm%itab=0                  
    strm%ilevel=0                
    strm%iflowlevel=0            
    strm%ilast=0                 
    strm%icommentline=0          
    strm%linetab=0
    strm%ievt_flow=0                                     
    strm%flow_events=0
    nullify(strm%dict_warning)
  end function stream_null

                                        
  !> Set the default stream of the module. Return  a STREAM_ALREADY_PRESENT errcode if 
  !! The stream has not be initialized.
  subroutine yaml_set_default_stream(unit,ierr)
    implicit none
    integer, intent(in) :: unit
    integer, intent(out) :: ierr
    !local variables
    integer :: istream

    !check if the stream is present
    call get_stream(unit,istream,istat=ierr)
    if (ierr==0) then
       default_stream=istream
    end if

  end subroutine yaml_set_default_stream


  !> Get the default stream unit
  subroutine yaml_get_default_stream(unit)
    implicit none
    integer, intent(out) :: unit

    unit=stream_units(default_stream)

  end subroutine yaml_get_default_stream


  !> Set all the output from now on to the file indicated by stdout
  subroutine yaml_set_stream(unit,filename,istat,tabbing,record_length)
    implicit none
    integer, optional, intent(in) :: unit              !< File unit (by default 6)
    integer, optional, intent(in) :: tabbing           !< Indicate a tabbing for the stream (0 no tabbing, default)
    integer, optional, intent(in) :: record_length     !< Maximum length of a record
    character(len=*), optional, intent(in) :: filename !< Filename of the stream
    integer, optional, intent(out) :: istat            !< Status
    !local variables
    integer, parameter :: NO_ERRORS           = 0
    integer :: istream,unt,ierr

    if (present(istat)) istat=NO_ERRORS !so far

    if (present(unit)) then
       unt=unit
    else
       unt=6
    end if

    !check if unit has been already assigned
    do istream=1,active_streams
       if (unt==stream_units(istream)) then
          if (present(istat)) then
             istat=STREAM_ALREADY_PRESENT
             return
          else
             stop 'yaml_set_stream:unit already present'
          end if
       end if
    end do
    !assign the unit to the new stream
    active_streams=active_streams+1
    !initalize the stream
    streams(active_streams)=stream_null()
    streams(active_streams)%unit=unt
    stream_units(active_streams)=unt

    ! set last opened stream as default stream
    default_stream=active_streams

    !open fortran unit if needed
    if (present(filename) .and. unt /= 6) then
       open(unit=unt,file=trim(filename),status='unknown',position='append',iostat=ierr)
       if (present(istat)) then
          istat=ierr
       else
          if (ierr /=0) then
             stop 'error in file opening'
          end if
       end if
    end if

    !set stream non-default attributes
    streams(active_streams)%unit=unt

    if (present(tabbing)) then
       streams(active_streams)%tabref=tabbing
       if (tabbing==0) streams(active_streams)%pp_allowed=.false.
    end if
    if (present(record_length)) then
       streams(active_streams)%max_record_length=record_length
    end if
  end subroutine yaml_set_stream


  !> Get the attributes of the stream at present
  !! Display is dump=.true.
  subroutine yaml_stream_attributes(unit,stream_unit,&
       icursor,flowrite,itab_active,iflowlevel,ilevel,ilast,indent,indent_previous,&
       record_length)
    implicit none
    integer, intent(in) , optional :: unit          !< File unit to display
    integer, intent(in) , optional :: stream_unit   !< Stream Id 
    logical, intent(out), optional :: flowrite
    integer, intent(out), optional :: icursor,itab_active
    integer, intent(out), optional :: iflowlevel,ilevel,ilast
    integer, intent(out), optional :: indent,indent_previous,record_length
    !local variables
    logical :: dump,flowritet
    integer :: sunt,unt,strm,icursort,itab_activet
    integer :: iflowlevelt,ilevelt,ilastt,indentt
    integer :: indent_previoust,record_lengtht
    integer, dimension(tot_max_record_length/tab) :: linetab

    !writing unit
    unt=0
    if (present(unit)) unt=unit
    !stream to be analyzed
    sunt=0
    if (present(stream_unit)) sunt=unit
    call get_stream(sunt,strm)

    !copy the values
    icursort=streams(strm)%icursor
    flowritet=streams(strm)%flowrite
    iflowlevelt=streams(strm)%iflowlevel
    ilevelt=streams(strm)%ilevel
    ilastt=streams(strm)%ilast
    itab_activet=streams(strm)%itab_active
    linetab=streams(strm)%linetab
    indentt=streams(strm)%indent
    indent_previoust=streams(strm)%indent_previous
    record_lengtht=streams(strm)%max_record_length

    dump=.true.
    !check if the variables have to be imported or not
    if (present(icursor)) then
       icursor=icursort
       dump=.false.
    end if
    if (present(flowrite)) then
       flowrite=flowritet
       dump=.false.
    end if
    if (present(indent)) then
       indent=indentt
       dump=.false.
    end if
    if (present(indent_previous)) then
       indent_previous=indent_previoust
       dump=.false.
    end if
    if (present(itab_active)) then
       itab_active=itab_activet
       dump=.false.
    end if
    if (present(iflowlevel)) then
       iflowlevel=iflowlevelt
       dump=.false.
    end if
    if (present(ilevel)) then
       ilevel=ilevelt
       dump=.false.
    end if
    if (present(ilast)) then
       ilast=ilastt
       dump=.false.
    end if
    if (present(record_length)) then
       record_length=record_lengtht
       dump=.false.
    end if


    if (dump) then
       call yaml_newline(unit=unt)
       call yaml_open_map('Attributes of the Stream',unit=unt)
       call yaml_map('Cursor position',icursort,unit=unt)
       call yaml_map('Max. Record Length',record_lengtht,unit=unt)
       call yaml_map('Indent value',indentt,unit=unt)
       call yaml_map('Indent value Saved',indent_previoust,unit=unt)
       call yaml_map('Write in Flow',flowritet,unit=unt)
       call yaml_map('Flow Level',iflowlevelt,unit=unt)
         call yaml_map('Level',ilevelt,unit=unt)
         call yaml_map('Last Level (flow==.false.)',ilastt,unit=unt)
       call yaml_map('Active Tabulars',itab_activet,unit=unt)
       if (itab_activet>0) call yaml_map('Tabular Values',linetab(1:itab_activet),unit=unt)
       call yaml_close_map(unit=unt)
       call yaml_newline(unit=unt)
    end if
  end subroutine yaml_stream_attributes


  !> Create a new document
  !! Put document_closed to .false.
  !! Check if already used before yaml_release_document by testing document_closed
  !! In this case, do nothing
  subroutine yaml_new_document(unit)
    implicit none
    integer, optional, intent(in) :: unit
    !local variables
    integer :: unt,strm

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    if (streams(strm)%document_closed) then
       !Check all indentation
       if (streams(strm)%indent /= 1) then
          call yaml_warning("Indentation error. Yaml Document has not been closed correctly",unit=stream_units(strm))
          streams(strm)%indent=1
       end if
       call dump(streams(strm),'---',event=DOCUMENT_START)
       streams(strm)%flow_events=NONE
       streams(strm)%document_closed=.false.
    end if

  end subroutine yaml_new_document


  !> After this routine is called, the new_document will become effective again
  subroutine yaml_release_document(unit)
    implicit none
    integer, optional, intent(in) :: unit  !< Stream Identity number
    !local variables
    integer :: unt,strm,unit_prev

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    !here we should print the warnings which have been obtained
    if (associated(streams(strm)%dict_warning)) then
       call yaml_newline()
       call yaml_comment('Warnings obtained during the run, check their relevance!',hfill='-')
       call yaml_dict_dump(streams(strm)%dict_warning)
       call dict_free(streams(strm)%dict_warning)
    end if

    !Initialize the stream, keeping the file unit
    unit_prev=streams(strm)%unit
    streams(strm)=stream_null()
    streams(strm)%unit=unit_prev

  end subroutine yaml_release_document

  subroutine yaml_close_all_streams()
    implicit none
    
    !local variables
    integer :: istream,unt,unts

    do istream=1,active_streams
       unt=stream_units(istream)
       unts=streams(istream)%unit
       if (unts /= unt) stop 'YAML close streams: unit inconsistency'
       !close files which are not stdout
       if (unt /= 6) close(unt)
       !reset the stream information
       stream_units(istream)=6
       streams(istream)=stream_null()
    end do
    active_streams=1 !stdout is always kept active
    default_stream=1
  end subroutine yaml_close_all_streams


  !> Display a warning (yaml comment starting with '#WARNING: ')
  subroutine yaml_warning(message,level,unit)
    implicit none
    integer, optional, intent(in) :: unit
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: level
    !local variables
    integer :: unt,strm,item
    type(dictionary), pointer :: dict_tmp

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    call dump(streams(strm),' #WARNING: '//trim(message))
    !here we should add a collection of all the warning which are printed out in the code.
    if (.not. streams(strm)%document_closed) then
       if (.not. associated(streams(strm)%dict_warning)) then
          call dict_init(streams(strm)%dict_warning)
          call set(streams(strm)%dict_warning//'WARNINGS'//0,trim(message))
       else
          !add the warning as a list
          dict_tmp=>streams(strm)%dict_warning//'WARNINGS'
          item=dict_tmp%data%nitems
          call set(dict_tmp//item,trim(message))
       end if
    end if
    if (present(level)) then
       if (level <= streams(strm)%Wall) then
          call dump(streams(strm),' Critical warning level reached, aborting...')
          call yaml_release_document(unit=unt)
          stop
       end if
    end if
  end subroutine yaml_warning


  !> Write a yaml comment (#......)
  !! Split the comment if too long
  subroutine yaml_comment(message,advance,unit,hfill,tabbing)
    implicit none
    character(len=*), intent(in) :: message           !< The given comment (without #)
    character(len=*), optional, intent(in) :: advance !< Advance or not
    integer, optional, intent(in) :: unit             !< Unit of the stream (by default unit=0)
    character(len=1), optional, intent(in) :: hfill   !< If present fill the line with the given character
    integer, optional, intent(in) :: tabbing          !< Number of space for tabbing
    !Local variables
    integer :: unt,strm,msg_lgt,tb,ipos
    integer :: lstart,lend,lmsg,lspace,hmax
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    !comment to be written
    if (present(advance)) then
       adv=advance
    else
       adv='yes'
    end if

    !Beginning of the message
    lstart=1
    !Length of the message to write (without blank characters)
    lmsg=len_trim(message)

    !Split the message if too long
    do
       !Position of the cursor
       ipos=max(streams(strm)%icursor,streams(strm)%indent)

       msg_lgt=0
       if (present(tabbing)) then
          tb=max(tabbing-ipos-1,1)
          call buffer_string(towrite,len(towrite),repeat(' ',tb),msg_lgt)
          ipos=ipos+tb
       end if

       !Detect the last character of the message
       lend=len_trim(message(lstart:))
       if (lend+msg_lgt > streams(strm)%max_record_length) then
          !We have an error from buffer_string so we split it!
          !-1 to be less and -2 for the character '#'
          lend=streams(strm)%max_record_length-msg_lgt-2
          !We are looking for the first ' ' from the end
          lspace=index(message(lstart:lstart+lend-1),' ',back=.true.)
          if (lspace /= 0) then
             lend = lspace
          end if
       end if
       call buffer_string(towrite,len(towrite),message(lstart:lstart+lend-1),msg_lgt)

       !Check if possible to hfill
       hmax = max(streams(strm)%max_record_length-ipos-len_trim(message)-3,0)
       if (present(hfill) .and. hmax > 0) then
          !Fill with the given character and dump
          call dump(streams(strm),repeat(hfill,hmax)//' '//towrite(1:msg_lgt),advance=adv,event=COMMENT)
       else
          !Dump the string towrite into the stream
          call dump(streams(strm),towrite(1:msg_lgt),advance=adv,event=COMMENT)
       end if

       !Check if all the message is written
       !So we start from iend+1
       lstart=lstart+lend
       if (lstart>lmsg) then
          exit
       end if
    end do

  end subroutine yaml_comment


  !> Write a scalar variable, takes care of indentation only
  subroutine yaml_scalar(message,advance,unit,hfill)
    implicit none
    character(len=1), optional, intent(in) :: hfill
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: unit
    character(len=*), intent(in), optional :: advance
    !local variables
    integer :: unt,strm
    character(len=3) :: adv

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    !comment to be written
    if (present(advance)) then
       adv=advance
    else
       adv='yes'
    end if
    if (present(hfill)) then
       call dump(streams(strm),&
            repeat(hfill,&
            max(streams(strm)%max_record_length-&
            max(streams(strm)%icursor,streams(strm)%indent)-&
            len_trim(message)-3,0))//' '//trim(message),&
            advance=adv,event=COMMENT)
    else
       call dump(streams(strm),trim(message),advance=adv,event=SCALAR)
    end if

  end subroutine yaml_scalar


  !> Open a yaml map (dictionary) 
  subroutine yaml_open_map(mapname,label,flow,unit)
    implicit none
    integer, optional, intent(in) :: unit
    character(len=*), optional, intent(in) :: mapname
    logical, optional, intent(in) :: flow
    character(len=*), optional, intent(in) :: label
    !local variables
    logical :: doflow
    integer :: msg_lgt
    integer :: unt,strm
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    doflow=streams(strm)%flowrite
    !override if already active
    if (present(flow)) doflow=flow .or. doflow

    msg_lgt=0
    !put the message
    if (present(mapname)) then
       call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
       !put the semicolon
       call buffer_string(towrite,len(towrite),':',msg_lgt)
    end if
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label),msg_lgt)
    end if

    call open_level(streams(strm),doflow)

    if (doflow .or. msg_lgt==0) then
       adv='no '
    else
       adv='yes'
    end if

    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING_START)

  end subroutine yaml_open_map


  !> Close the map
  subroutine yaml_close_map(advance,unit)
    implicit none
    integer, optional, intent(in) :: unit
    character(len=*), optional, intent(in) :: advance
    !local variables
    integer :: unt,strm
    character(len=3) :: adv
    logical :: doflow

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    if (streams(strm)%iflowlevel > 1) then
       adv='no'
    else
       adv='yes'
    end if
    if (present(advance)) adv=advance

    call dump(streams(strm),' ',advance=trim(adv),event=MAPPING_END)

    doflow = (streams(strm)%flowrite)
    call close_level(streams(strm),doflow)

  end subroutine yaml_close_map


  !> Open a yaml sequence
  subroutine yaml_open_sequence(mapname,label,flow,advance,unit)
    use yaml_strings
    implicit none
    character(len=*), optional, intent(in) :: mapname !< Key of the sequence
    character(len=*), optional, intent(in) :: label   !< Add a label to be referenced as &xxx
    logical, optional, intent(in) :: flow             !< .true.  Add [ and represent the sequence as a stream
                                                      !! .false. Add a flow level (go to the line and indent)
    character(len=*), optional, intent(in) :: advance !< Same option as write
    integer, optional, intent(in) :: unit             !< File unit
    !local variables
    logical :: doflow
    integer :: msg_lgt
    integer :: unt,strm
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    doflow=streams(strm)%flowrite
    !override if already active
    if (present(flow)) doflow=flow .or. doflow

    msg_lgt=0
    !put the message
    if (present(mapname)) then
       call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
       !put the semicolon
       call buffer_string(towrite,len(towrite),':',msg_lgt)
    end if
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label),msg_lgt)
    end if

    call open_level(streams(strm),doflow)

    if (doflow .or. msg_lgt==0) then
       adv='no '
    else
       adv='yes'
       if (present(advance)) adv = advance
    end if

    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=SEQUENCE_START)

  end subroutine yaml_open_sequence


  !> Close a yaml sequence
  subroutine yaml_close_sequence(advance,unit)
    implicit none
    character(len=*), optional, intent(in) :: advance
    integer, optional, intent(in) :: unit
    !local variables
    integer :: unt,strm
    character(len=3) :: adv
    logical :: doflow

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    if (streams(strm)%iflowlevel > 1) then
       adv='no'
    else
       adv='yes'
       if (present(advance)) adv=advance
    end if

    call dump(streams(strm),' ',advance=trim(adv),event=SEQUENCE_END)

    doflow = (streams(strm)%flowrite)
    call close_level(streams(strm),doflow)

  end subroutine yaml_close_sequence


  !> Add a new line in the flow 
  !! This routine has a effect only if a flow writing is active
  subroutine yaml_newline(unit)
    implicit none
    integer, optional, intent(in) :: unit
    !local variables
    integer :: unt,strm

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    if (streams(strm)%icursor > 1) then
    !if (streams(strm)%flowrite .and. streams(strm)%icursor > 1) then
       call dump(streams(strm),' ',advance='yes',event=NEWLINE)
    end if
  end subroutine yaml_newline


  !> Write directly a yaml sequence
  subroutine yaml_sequence(seqvalue,label,advance,unit)
    implicit none
    integer, optional, intent(in) :: unit
    character(len=*), optional, intent(in) :: label,seqvalue,advance
    !local variables 
    integer :: msg_lgt,unt,strm
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if
    !put the value
    if (present(seqvalue)) &
         call buffer_string(towrite,len(towrite),trim(seqvalue),msg_lgt)

    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=SEQUENCE_ELEM)
  end subroutine yaml_sequence


  !> Do a yaml map
  subroutine yaml_map(mapname,mapvalue,label,advance,unit)
    implicit none
    character(len=*), intent(in) :: mapname             !< key
    character(len=*), intent(in) :: mapvalue            !< value
    character(len=*), optional, intent(in) :: label     !< label for reference (&xxx)
    character(len=*), optional, intent(in) :: advance   !< advance or not
    integer, optional, intent(in) :: unit               !< unit for strem
    !local variables
    logical :: cut,redo_line
    integer :: msg_lgt,strm,unt,icut,istr,ierr,msg_lgt_ck
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0

    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),'&',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label),msg_lgt)
    end if

    !while putting the message verify that the string is not too long
    msg_lgt_ck=msg_lgt
    call buffer_string(towrite,len(towrite),trim(mapvalue),msg_lgt,istat=ierr)
    if (ierr ==0) then
       call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING,istat=ierr)
    end if
    redo_line=ierr/=0
    if (redo_line) then
       if (streams(strm)%flowrite) then
          call dump(streams(strm),towrite(1:msg_lgt_ck),advance=trim(adv),event=SCALAR)
       else
          if (present(label)) then
             call yaml_open_map(mapname,label=label,unit=unt)
          else
             call yaml_open_map(mapname,unit=unt)
          end if
       end if
!       if (streams(strm)%flowrite) call yaml_newline(unit=unt)
       icut=len_trim(mapvalue)
       istr=1
       cut=.true.
       msg_lgt=0
       cut_line: do while(cut)
          !print *,'hereOUTPU',cut,icut
       !verify where the message can be cut
          cut=.false.
          cut_message :do while(icut > streams(strm)%max_record_length - max(streams(strm)%icursor,streams(strm)%indent))
             icut=index(trim((mapvalue(istr:istr+icut-1))),' ',back=.true.)
             cut=.true.
          end do cut_message
          call buffer_string(towrite,len(towrite),mapvalue(istr:istr+icut-1),msg_lgt)
          if (streams(strm)%flowrite .and. .not. cut) call buffer_string(towrite,len(towrite),',',msg_lgt)
          call dump(streams(strm),towrite(1:msg_lgt),advance='yes',event=SCALAR)
          istr=icut
          icut=len_trim(mapvalue)-istr+1
          !print *,'icut',istr,icut,mapvalue(istr:istr+icut-1),cut,istr+icut-1,len_trim(mapvalue)
          msg_lgt=0
       end do cut_line
       if (.not.streams(strm)%flowrite) call yaml_close_map(unit=unt)
    end if

  end subroutine yaml_map


  subroutine yaml_map_i(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    character(len=*), intent(in) :: mapname
    integer, intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,unt
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if
    !put the value
    if (present(fmt)) then
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
    else
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    end if
    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
  end subroutine yaml_map_i


  subroutine yaml_map_f(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    character(len=*), intent(in) :: mapname
    real, intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,unt
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if
    !put the value
    if (present(fmt)) then
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
    else
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    end if
    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
  end subroutine yaml_map_f


  subroutine yaml_map_d(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    character(len=*), intent(in) :: mapname
    real(kind=8), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,unt
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)


    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if

    !put the value
    if (present(fmt)) then
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
    else
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    end if
    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
  end subroutine yaml_map_d


  subroutine yaml_map_l(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    character(len=*), intent(in) :: mapname
    logical,  intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,unt
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)


    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if
    !put the value
    if (present(fmt)) then
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
    else
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    end if
    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
  end subroutine yaml_map_l


  subroutine yaml_map_dv(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    character(len=*), intent(in) :: mapname
    real(kind=8), dimension(:), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,unt,nl,nu,tmp_lgt,i
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    adv='def' !default value
    if (present(advance)) adv=advance

    nl=lbound(mapvalue,1)
    nu=ubound(mapvalue,1)

    msg_lgt=0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if

    !check whether the final message will be too long or not
    if (nu-nl > 0) then
       !change strategy if the remaining space is too low
       !template of an element
       if (present(fmt)) then
          tmp_lgt=len_trim(yaml_toa(mapvalue(nl),fmt=fmt))
       else
          tmp_lgt=len_trim(yaml_toa(mapvalue(nl)))
       end if
       tmp_lgt=tmp_lgt+3 !comma and spaces
       tmp_lgt=tmp_lgt*(nu-nl)
       if (max(streams(strm)%icursor+msg_lgt+1,streams(strm)%tabref)+tmp_lgt > &
            streams(strm)%max_record_length) then
          !implement the writing explicitly per element
          call yaml_open_sequence(mapname,flow=.true.,unit=unt)
          do i=nl,nu
             if (present(fmt)) then
                call yaml_sequence(trim(yaml_toa(mapvalue(i),fmt=fmt)),unit=unt)
             else
                call yaml_sequence(trim(yaml_toa(mapvalue(i))),unit=unt)
             end if
          end do
          call yaml_close_sequence(unit=unt)
          return
       end if
    end if

    !put the value
    if (present(fmt)) then
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
    else
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    end if
    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
  end subroutine yaml_map_dv


  !> Character vector
  subroutine yaml_map_cv(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    character(len=*), intent(in) :: mapname
    character(len=*), dimension(:), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,unt
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if
    !put the value
    if (present(fmt)) then
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    else
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    end if
    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
  end subroutine yaml_map_cv


  subroutine yaml_map_iv(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    character(len=*), intent(in) :: mapname
    integer, dimension(:), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,unt
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if
    !put the value
    if (present(fmt)) then
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
    else
       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
    end if
    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
  end subroutine yaml_map_iv


  !> Get the stream, initialize if not already present (except if istat present)
  subroutine get_stream(unt,strm,istat)
    implicit none
    integer, intent(in) :: unt
    integer, intent(out) :: strm
    integer, optional, intent(out) :: istat
    !local variables
    logical :: stream_found
    integer :: istream,prev_def

    if (present(istat)) istat=0

    if (unt==0) then
       strm=default_stream
    else
       !it is assumed that the unit exists
       stream_found=.false.
       do istream=1,active_streams
          if (stream_units(istream)==unt) then
             strm=istream
             stream_found=.true.
             exit
          end if
       end do
       if (.not. stream_found) then
          if (present(istat)) then
             istat=STREAM_ALREADY_PRESENT
          else
             !otherwise initialize it, no pretty printing
             prev_def=default_stream
             call yaml_set_stream(unit=unt,tabbing=0)
             strm=default_stream
             !but do not change default stream
             default_stream=prev_def
          end if
       end if
    end if

  end subroutine get_stream


  subroutine dump(stream,message,advance,event,istat)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: event
    integer, intent(out), optional :: istat
    !local variables
    logical :: ladv,change_line,reset_line,pretty_print,reset_tabbing,comma_postponed
    integer :: evt,indent_lgt,msg_lgt,shift_lgt,prefix_lgt
    integer :: towrite_lgt
    character(len=1) :: anchor
    character(len=3) :: adv
    character(len=5) :: prefix
    character(len=stream%max_record_length) :: towrite

    if(present(istat)) istat=0 !no errors

    if (present(event)) then
       evt=event
    else !default event: scalar value
       evt=SCALAR
    end if

    !decide whether to write advanced or not
    !decide if advanced output or not
    ladv=.not.stream%flowrite
    if (present(advance)) then
       if (trim(advance)=='no' .or. trim(advance)=='NO') then
          ladv=.false.
       else if (trim(advance)=='yes' .or. trim(advance)=='YES') then
          ladv=.true.
       end if
    end if
    if (ladv) then
       adv='yes'
    else
       adv='no '
    end if

    !decide whether the line has to be reset (no by default)
    reset_line=.false.

    !decide whether the line has to be continuated (no by default)
    change_line=.false.

    !possible indentation (depending of the event) and of the cursor
    indent_lgt=indent_value(stream,evt)
    !write(*,fmt='(a,i0,a)',advance="no") '(lgt ',indent_lgt,')'

    !calculate the number of objects to be written before
    !these objects should go to the active line in case of a new line
    !string length, and message body
    !initialize it
    towrite=repeat(' ',len(towrite))
    msg_lgt=0
    !a empty message is not written
    if (len_trim(message) > 0) &
         call buffer_string(towrite,len(towrite),message,msg_lgt)

    prefix_lgt=0
    !initialize it
    prefix=repeat(' ',len(prefix))
    !write(stdout,*)'NOT DONE',icomma,flowrite,iflowlevel
    !msg_lgt should be added to the function
    if (put_comma(stream,evt)) then!stream%icomma==1 .and. stream%flowrite) then
       call buffer_string(prefix,len(prefix),', ',prefix_lgt)
    end if
    !next time comma should be postponed
    comma_postponed=comma_not_needed(evt) .or. (flow_is_ending(evt) .and. stream%iflowlevel ==1)

    !no pretty printing by default
    pretty_print=.false.
    shift_lgt=0

    !reset_tabbing is disabled
    reset_tabbing=.false.

    !set module variables according to the event
    select case(evt)

    case(SEQUENCE_START)

       if (.not.stream%flowrite) then
          call open_indent_level(stream)
       else
          call buffer_string(towrite,len(towrite),' [',msg_lgt)
          !comma has to be written afterwards, if there is a message
          !stream%flowrite=-1
          stream%flowrite=.true.
          !added for pretty printing
          reset_tabbing=.true.
       end if

    case(SEQUENCE_END)
       !print *,'here',prefix_lgt,prefix,icomma,flowrite,iflowlevel

       if (.not.stream%flowrite) then
          call close_indent_level(stream)
       else
          if (stream%iflowlevel > 1 .and. ladv) then
             call buffer_string(prefix,len(prefix),']',prefix_lgt,back=.true.)
             !stream%flowrite=-1
             stream%flowrite=.true.
          else
             call buffer_string(prefix,len(prefix),']',prefix_lgt)
          end if
          reset_line=ladv
       end if

    case(MAPPING_START)

       if (.not.stream%flowrite) then
          call open_indent_level(stream)
       else
          !write(stdout,*)'here',prefix,'there',icomma,flowrite,iflowlevel
          call buffer_string(towrite,len(towrite),' {',msg_lgt)
          !stream%flowrite=-1
          stream%flowrite=.true.
          reset_tabbing=.true.
       end if
       !write(*,fmt='(a,i0,a,a,a)',advance='no') '|',stream%indent,'|',trim(towrite),'|'

       !pretty_print=.true. .and. stream%pp_allowed

    case(MAPPING_END)

       if (.not.stream%flowrite) then
          call close_indent_level(stream)
       else
          if (stream%iflowlevel > 1 .and. ladv) then
             call buffer_string(prefix,len(prefix),'}',prefix_lgt,back=.true.)
             !flowrite=-1
             reset_line=.true.
          else 
             call buffer_string(prefix,len(prefix),'}',prefix_lgt)
          end if
          reset_line=ladv
       end if

    case(COMMENT)
       if (stream%icommentline==0) then !no comment active
          call buffer_string(prefix,len(prefix),' #',prefix_lgt)
       end if
       if (.not. ladv) then
          stream%icommentline=1
       else
          reset_line=.true.
       end if

    case(MAPPING)

       pretty_print=.true. .and. stream%pp_allowed
       anchor=':'

    case(SEQUENCE_ELEM)

       if (.not.stream%flowrite) then
          !lower indent and update prefix
          indent_lgt=indent_lgt-2
          call buffer_string(prefix,len(prefix),'- ',prefix_lgt)
       else
          if (msg_lgt>0) comma_postponed=.false.
          !just added to change the line
          pretty_print=.true. .and. stream%pp_allowed
          anchor='.'
       end if

    case(SCALAR)

    case(NEWLINE)

       if (stream%flowrite) then
          !print *,'NEWLINE:',stream%flowrite
          change_line=.true.
          !stream%flowrite=-1
          stream%flowrite=.true.
          reset_line=ladv
          msg_lgt=0
       else
          change_line=.true.
          reset_line=.true.
          msg_lgt=0
       end if

    end select

    !adjust the towrite string to match with the closest tabular
    if (pretty_print) then
       call pretty_printing(stream%flowrite,anchor,towrite,&
            stream%icursor,indent_lgt,prefix_lgt,&
            msg_lgt,stream%max_record_length,shift_lgt,change_line)
    end if

    !standard writing,
    if (change_line) then
       !first write prefix, if needed
       if (prefix_lgt>0) then
          write(stream%unit,'(a)')prefix(1:prefix_lgt)
       else if (msg_lgt >0 .or. evt == NEWLINE) then
          !change line
          write(stream%unit,*)
       end if
       stream%icursor=1
       towrite_lgt=msg_lgt+shift_lgt
    else
       call shiftstr(towrite,prefix_lgt)
       if (prefix_lgt > 0)towrite(1:prefix_lgt)=prefix(1:prefix_lgt)
       towrite_lgt=prefix_lgt+msg_lgt+shift_lgt
    end if
    !print *,'adv',trim(adv),towrite_lgt,icursor,change_line,msg_lgt
    !here we should check whether the size of the string exceeds the maximum length
    if (towrite_lgt > 0) then
       if (towrite_lgt > stream%max_record_length) then
          if (present(istat)) then
             istat=-1
             return
          else
             !crop the writing 
             towrite_lgt=stream%max_record_length
             !stop 'ERROR (dump): writing exceeds record size'
          end if
       else
          !write(*,fmt='(a,i0,a)',advance="no") '(indent_lgt ',indent_lgt,')'
          write(stream%unit,'(a)',advance=trim(adv))repeat(' ',max(indent_lgt,0))//towrite(1:towrite_lgt)
       end if
    end if

    !if advancing i/o cursor is again one
    if (ladv) then
       stream%icursor=1
    else
       !cursor after writing
       stream%icursor=stream%icursor+indent_lgt+towrite_lgt
    end if

    if (reset_tabbing) then
       stream%itab_active=0
       stream%itab=0
    end if

    if (reset_line) call carriage_return(stream)

    !keep history of the event for a flowrite
    !needed for the comma
    if (stream%flowrite) then
       stream%ievt_flow=modulo(stream%ievt_flow,tot_max_flow_events)+1 !to avoid boundary problems
       if (comma_postponed) then
          stream%flow_events(stream%ievt_flow)=evt
       else
          stream%flow_events(stream%ievt_flow)=COMMA_TO_BE_PUT
       end if
    else
       stream%ievt_flow=0
    end if

  contains

    subroutine pretty_printing(rigid,anchor,message,icursor,&
         indent_lgt,prefix_lgt,msg_lgt,max_lgt,shift_lgt,change_line)
      implicit none
      logical, intent(in) :: rigid
      integer, intent(in) :: icursor,prefix_lgt,msg_lgt,max_lgt
      integer, intent(inout) :: indent_lgt
      character(len=*), intent(in) :: anchor
      character(len=*), intent(inout) :: message
      logical, intent(out) :: change_line
      integer, intent(out) :: shift_lgt
      !local variables
      integer :: iscpos,ianchor_pos,tabeff

      change_line=.false.
      iscpos=index(message,anchor)
      shift_lgt=0
      if (iscpos==0) return !no anchor, no pretty printing
      ianchor_pos=icursor+prefix_lgt+indent_lgt+iscpos-1
      call closest_tab(ianchor_pos,tabeff)
      !first attempt to see if the line enters
      shift_lgt=tabeff-ianchor_pos
      !print *, 'there',tabeff,itab,ianchor_pos,shift_lgt,msg_lgt,prefix_lgt,indent_lgt,icursor
      !see if the line enters
      if (icursor+msg_lgt+prefix_lgt+indent_lgt+shift_lgt >= max_lgt) then
         !restart again 
         change_line=.true.
         !reset newly created tab
         if (stream%itab==stream%itab_active .and. stream%itab > 1)&
              stream%itab_active=max(stream%itab_active-1,0)
         stream%itab=1
         if (indent_lgt==0) indent_lgt=1
         ianchor_pos=indent_lgt+iscpos
         call closest_tab(ianchor_pos,tabeff)

         shift_lgt=tabeff-ianchor_pos
      end if
      !print *, 'here',tabeff,itab,ianchor_pos,shift_lgt,change_line
      !at this point we know the size of the message.
      !we know also whether to write it or to pass to the following line
      !once the tabbing has been decided, adjust the message to the anchor
      call align_message(rigid,len(message),shift_lgt+iscpos,anchor,message)

    end subroutine pretty_printing


    !> Calculate the reference tabular value
    subroutine closest_tab(ianchor_pos,tabeff)
      implicit none
      integer, intent(in) :: ianchor_pos
      integer, intent(out) :: tabeff

      if (stream%flowrite) then
         !first check that the tabbing is already done, otherwise add another tab
         if (stream%itab < stream%itab_active) then
            !realign the value to the tabbing
            do 
               if (ianchor_pos <= stream%linetab(stream%itab) .or. &
                    stream%itab==stream%itab_active) exit
               stream%itab=modulo(stream%itab,tot_max_record_length/tab)+1
            end do
         end if

         if (stream%itab < stream%itab_active .and. stream%itab>0) then
            tabeff=stream%linetab(stream%itab)
         else
            tabeff=ianchor_pos
            stream%itab=modulo(stream%itab,tot_max_record_length/tab)+1
            stream%itab_active=modulo(stream%itab_active,tot_max_record_length/tab)+1
            stream%linetab(stream%itab_active)=tabeff
         end if
      else
         !for the moment do not check compatibility of the line
         tabeff=max(stream%tabref,ianchor_pos)
      end if
    end subroutine closest_tab

    function indent_value(stream,evt)
      implicit none
      integer, intent(in) :: evt
      type(yaml_stream), intent(in) :: stream
      integer :: indent_value

      !write(*,fmt='(a,i0,i3,a)',advance="no") '(stream_indent ',stream%indent,stream%flowrite,')'
      if (.not.stream%flowrite .and. stream%icursor==1) then
         indent_value=stream%indent!max(stream%indent,0) !to prevent bugs
      !if first time in the flow recuperate the saved indent
      else if (stream%icursor==1 .and. stream%iflowlevel==1 & 
           .and. stream%ievt_flow==0) then
         indent_value=stream%indent_previous
      else
         indent_value=0!1
         if (stream%icursor==1) indent_value=1
      end if

      if (evt==DOCUMENT_START) indent_value=0

!      if (stream%icursor > 1) then
!         indent_value=0
!      else if(.not.stream%flowrite) then !other conditions have to be added here
!         indent_value=max(stream%indent,0) !to prevent bugs
!      else
!         indent_value=stream%indent_previous
!      end if
      
    end function indent_value

    !> Decide whether comma has to be put
    function put_comma(stream,evt)
      implicit none
      integer, intent(in) :: evt
      type(yaml_stream), intent(inout) :: stream
      logical :: put_comma

      put_comma=stream%flowrite .and. stream%ievt_flow>0

      if (stream%ievt_flow > 0) then
         put_comma=stream%flow_events(stream%ievt_flow) == COMMA_TO_BE_PUT
         if (.not. put_comma .and. comma_potentially_needed(evt)) then
            !print *,'event'
            !control whether there is a ending flow
            !if (stream%iflowlevel > 1 .and. 
         end if
      end if
      !in any case the comma should not be put before a endflow
      if (flow_is_ending(evt)) put_comma=.false.

    end function put_comma

  end subroutine dump


  function flow_is_starting(evt)
    implicit none
    integer, intent(in) :: evt
    logical flow_is_starting
    
    flow_is_starting=(evt==MAPPING_START .or. evt == SEQUENCE_START)

  end function flow_is_starting


  function flow_is_ending(evt)
    implicit none
    integer, intent(in) :: evt
    logical flow_is_ending

    flow_is_ending=(evt==MAPPING_END .or. evt == SEQUENCE_END)

  end function flow_is_ending


  function comma_not_needed(evt)
    implicit none
    integer, intent(in) :: evt
    logical :: comma_not_needed

    comma_not_needed=evt==NONE           .or. &
                     evt==MAPPING_START  .or. &
                     evt==SEQUENCE_START .or. &
                     evt==SCALAR         .or. &
                     evt==COMMENT        .or. &
                     evt==SEQUENCE_ELEM  .or. &
                     evt==NEWLINE 
  end function comma_not_needed


  function comma_potentially_needed(evt)
    implicit none
    integer, intent(in) :: evt
    logical :: comma_potentially_needed

    comma_potentially_needed=evt==MAPPING_START  .or. &
                             evt==SEQUENCE_START .or. &
                             evt==SCALAR

  end function comma_potentially_needed


!!$  subroutine write_prefix(prefix_lgt,prefix,event)
!!$    implicit none
!!$    integer, intent(in) :: prefix_lgt,event
!!$    character(len=prefix_lgt), intent(out) :: prefix
!!$    !local variables
!!$    logical :: write_comma
!!$    integer :: iflw,nopen,nclose,ievt
!!$    if (flowrite==0) then
!!$       write_comma=.false.
!!$    else
!!$       nopen=0
!!$       nclose=0
!!$       do iflw=1,ievt_flow
!!$          ievt=flow_events(iflw)
!!$          if (ievt == MAPPING_START .or. ievt == SEQUENCE_START) nopen=nopen+1
!!$          if (ievt == MAPPING_END .or. ievt == SEQUENCE_END) nclose=nclose+1
!!$       end do
!!$       write_comma=(nopen > nclose)
!!$    end if
!!$
!!$  end subroutine write_prefix


  !> Reset the line control quantities, and reset the indentation
  subroutine carriage_return(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    !if a yaml_comment is called put the has in front
    stream%icommentline=0
    !beginining of the line
    stream%icursor=1
    !no tabbing decided yet
    stream%itab_active=0
    stream%itab=0
    !all needed commas are placed in the previous line
  end subroutine carriage_return  


  !> Open a level
  subroutine open_level(stream,doflow)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    logical, intent(in) :: doflow
    stream%ilevel = stream%ilevel + 1
    if(doflow) then
       call open_flow_level(stream)
    else
       stream%ilast = stream%ilevel
    end if
  end subroutine open_level


  !> Open a flow level (Indent more)
  subroutine open_flow_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    if (.not.stream%flowrite) then
       if (stream%iflowlevel==0) stream%indent_previous=stream%indent
       stream%indent=1
    end if
    stream%iflowlevel=stream%iflowlevel+1
    if (.not.stream%flowrite) stream%flowrite=.true. !start to write
  end subroutine open_flow_level


  !> Close a level
  subroutine close_level(stream,doflow)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    logical, intent(in) :: doflow
    stream%ilevel = stream%ilevel - 1
    if(doflow) then
       call close_flow_level(stream)
    else
       stream%ilast = min(stream%ilevel,stream%ilast)
    end if
  end subroutine close_level


  !> Close a flow level (Indent less)
  subroutine close_flow_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    !lower the flowlevel
    stream%iflowlevel=stream%iflowlevel-1
    if (stream%iflowlevel==0) then
       stream%indent=stream%indent_previous
       stream%flowrite=.false.
       !reset the events in the flow
       stream%flow_events=NONE
       stream%ievt_flow=0
    else
       stream%indent=1
       stream%flowrite=.true.
    end if

  end subroutine close_flow_level


  !> Increase the indentation of the strean without changing the flow level
  subroutine open_indent_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    stream%indent=stream%indent+stream%indent_step
  end subroutine open_indent_level


  !> Decrease the indentation of the strean without changing the flow level
  subroutine close_indent_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    stream%indent=max(stream%indent-stream%indent_step,0) !to prevent bugs
  end subroutine close_indent_level

  !> Dump a dictionary
  subroutine yaml_dict_dump(dict,flow)
   implicit none
   type(dictionary), intent(in) :: dict   !< Dictionary to dump
   logical, intent(in), optional :: flow  !< if .true. inline
   !local variables
   logical :: flowrite
     character(len=3) :: adv

   flowrite=.false.
   if (present(flow)) flowrite=flow

     !TEST (the first dictionary has no key)
     !if (.not. associated(dict%parent)) then
     if (associated(dict%child)) then
        call yaml_dict_dump_(dict%child,flowrite)
     else
        if (flowrite) then
           adv='no '
        else
           adv='yes'
        end if
        call yaml_scalar(dict%data%value,advance=adv)
     end if

   contains
     recursive subroutine yaml_dict_dump_(dict,flowrite)
         implicit none
         type(dictionary), intent(in) :: dict
         logical, intent(in) :: flowrite

         if (associated(dict%child)) then
            !see whether the child is a list or not
            !print *trim(dict%data%key),dict%data%nitems
            if (dict%data%nitems > 0) then
               call yaml_open_sequence(trim(dict%data%key),flow=flowrite)
               call yaml_dict_dump_(dict%child,flowrite)
               call yaml_close_sequence()
            else
               if (dict%data%item >= 0) then
                  call yaml_sequence(advance='no')
                  call yaml_dict_dump_(dict%child,flowrite)
               else
                  call yaml_open_map(trim(dict%data%key),flow=flowrite)
                  !call yaml_map('No. of Elems',dict%data%nelems)
                  call yaml_dict_dump_(dict%child,flowrite)
                  call yaml_close_map()
               end if
            end if
         else 
            !print *,'ciao',dict%key,len_trim(dict%key),'key',dict%value,flowrite
            if (dict%data%item >= 0) then
               call yaml_sequence(trim(dict%data%value))
            else
               call yaml_map(trim(dict%data%key),trim(dict%data%value))
            end if
         end if
         if (associated(dict%next)) then
            call yaml_dict_dump_(dict%next,flowrite)
         end if

       end subroutine yaml_dict_dump_
  end subroutine yaml_dict_dump

end module yaml_output
