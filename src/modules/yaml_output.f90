module yaml_strings

  integer :: max_value_length=95

  interface yaml_toa
     module procedure yaml_itoa,yaml_litoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_dvtoa,yaml_ivtoa
  end interface
  private :: yaml_itoa,yaml_litoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_dvtoa,yaml_ivtoa,max_value_lenght

contains

  !> Add a buffer to a string and increase its length
  subroutine buffer_string(string,string_lgt,buffer,string_pos,back,istat)
    implicit none
    integer, intent(in) :: string_lgt
    integer, intent(inout) :: string_pos
    character(len=*), intent(in) :: buffer
    character(len=string_lgt), intent(inout) :: string
    logical, optional, intent(in) :: back
    integer, optional, intent(out) :: istat
    !local variables
    integer :: lgt_add

    if (present(istat)) istat=0 !no errors

    lgt_add=len(buffer)
    !do not copy strings which are too long
    if (lgt_add+string_pos > string_lgt) then
       if (present(istat)) then
          istat=-1
          return
       else
          stop 'ERROR (buffer string): string too long'
       end if
    end if
       
    if (lgt_add==0) return
    if (present(back)) then
       if (back) then
          call shiftstr(string,lgt_add)
          string(1:lgt_add)=buffer(1:lgt_add)
       else
          string(string_pos+1:string_pos+lgt_add)=buffer(1:lgt_add)
       end if
    else
       string(string_pos+1:string_pos+lgt_add)=buffer(1:lgt_add)
    end if

    string_pos=string_pos+lgt_add

  end subroutine buffer_string

  !> add the spaces necessary to align the first occurrence of a given anchor
  !! into a tabular value. Can be done either by moving rigidly the message or 
  !! by adding spaces between the anchor and the rest of the message
  subroutine align_message(rigid,maxlen,tabval,anchor,message)
    implicit none
    logical, intent(in) :: rigid
    integer, intent(in) :: maxlen
    integer, intent(in) :: tabval
    character(len=*), intent(in) :: anchor
    character(len=maxlen), intent(inout) :: message
    !local variables
    integer :: iscpos,ishift

    !cannot align, tabular too far
    if (tabval>maxlen) return

    iscpos=index(message,anchor)      
    ishift=tabval-iscpos
    if (rigid) then
       call shiftstr(message,ishift)
    else
       message=message(1:iscpos-1)//repeat(' ',ishift)//anchor//&
            message(iscpos+1:maxlen-ishift)  ! shift right 
    end if

  end subroutine align_message


  !> Convert integer to character
  function yaml_itoa(i,fmt)
    implicit none
    integer, intent(in) :: i
    character(len=max_value_length) :: yaml_itoa
    character(len=*), optional, intent(in) :: fmt

    yaml_itoa=repeat(' ',max_value_length)
    if (present(fmt)) then
       write(yaml_itoa,fmt)i
    else
       write(yaml_itoa,'(i0)')i
    end if

    yaml_itoa=yaml_adjust(yaml_itoa)

  end function yaml_itoa

  !> Convert longinteger to character
  function yaml_litoa(i,fmt)
    implicit none
    integer(kind=8), intent(in) :: i
    character(len=max_value_length) :: yaml_litoa
    character(len=*), optional, intent(in) :: fmt

    yaml_litoa=repeat(' ',max_value_length)
    if (present(fmt)) then
       write(yaml_litoa,fmt)i
    else
       write(yaml_litoa,'(i0)')i
    end if

    yaml_litoa=yaml_adjust(yaml_litoa)

  end function yaml_litoa


!!$
  !> Convert float to character
  function yaml_ftoa(f,fmt)
    implicit none
    real, intent(in) :: f
    character(len=max_value_length) :: yaml_ftoa
    character(len=*), optional, intent(in) :: fmt

    yaml_ftoa=repeat(' ',max_value_length)
    if (present(fmt)) then
       write(yaml_ftoa,fmt)f
    else
       write(yaml_ftoa,'(1pe17.9)')f
    end if

    yaml_ftoa=yaml_adjust(yaml_ftoa)


  end function yaml_ftoa
!!$
  !> Convert double to character
  function yaml_dtoa(d,fmt)
    implicit none
    real(kind=8), intent(in) :: d
    character(len=max_value_length) :: yaml_dtoa
    character(len=*), optional, intent(in) :: fmt

    yaml_dtoa=repeat(' ',max_value_length)
    if (present(fmt)) then
       write(yaml_dtoa,fmt)d
    else
       write(yaml_dtoa,'(1pe25.17)')d
    end if
    yaml_dtoa=yaml_adjust(yaml_dtoa)

  end function yaml_dtoa

  !> Convert logical to character
  function yaml_ltoa(l,fmt)
    implicit none
    logical, intent(in) :: l
    character(len=max_value_length) :: yaml_ltoa
    character(len=*), optional, intent(in) :: fmt

    yaml_ltoa=repeat(' ',max_value_length)

    if (present(fmt)) then
       write(yaml_ltoa,fmt)l
    else
       if (l) then
          write(yaml_ltoa,'(a3)')'Yes'
       else
          write(yaml_ltoa,'(a3)')'No'
       end if
    end if

    yaml_ltoa=yaml_adjust(yaml_ltoa)
  end function yaml_ltoa

  !> Convert vector of double to character
  function yaml_dvtoa(dv,fmt)
    implicit none
    real(kind=8), dimension(:), intent(in) :: dv
    character(len=max_value_length) :: yaml_dvtoa
    character(len=*), optional, intent(in) :: fmt
    !local variables
    character(len=max_value_length) :: tmp
    integer :: nl,nu,i,length,pos

    tmp=repeat(' ',max_value_length)
    yaml_dvtoa=tmp

    nl=lbound(dv,1)
    nu=ubound(dv,1)

    yaml_dvtoa(1:2)='[ '
    pos=3
    do i=nl,nu
       if (present(fmt)) then
          tmp=yaml_dtoa(dv(i),fmt=fmt)
       else
          tmp=yaml_dtoa(dv(i))
       end if
       length=len(trim(tmp))-1
       if (pos+length > max_value_length) exit
       yaml_dvtoa(pos:pos+length)=tmp(1:length+1)
       if (i < nu) then
          yaml_dvtoa(pos+length+1:pos+length+2)=', '
       else
          yaml_dvtoa(pos+length+1:pos+length+2)=' ]'
       end if
       pos=pos+length+3
    end do

    yaml_dvtoa=yaml_adjust(yaml_dvtoa)

  end function yaml_dvtoa

  !> Convert vector of integer to character
  function yaml_ivtoa(iv,fmt)
    implicit none
    integer, dimension(:), intent(in) :: iv
    character(len=max_value_length) :: yaml_ivtoa
    character(len=*), optional, intent(in) :: fmt
    !local variables
    character(len=max_value_length) :: tmp
    integer :: nl,nu,i,length,pos

    tmp=repeat(' ',max_value_length)
    yaml_ivtoa=tmp

    nl=lbound(iv,1)
    nu=ubound(iv,1)

    yaml_ivtoa(1:2)='[ '
    pos=3
    do i=nl,nu
       if (present(fmt)) then
          tmp=yaml_itoa(iv(i),fmt=fmt)
       else
          tmp=yaml_itoa(iv(i))
       end if
       length=len(trim(tmp))-1
       if (pos+length > max_value_length) exit
       yaml_ivtoa(pos:pos+length)=tmp(1:length+1)
       if (i < nu) then
          yaml_ivtoa(pos+length+1:pos+length+2)=', '
       else
          yaml_ivtoa(pos+length+1:pos+length+2)=' ]'
       end if
       pos=pos+length+3
    end do

    yaml_ivtoa=yaml_adjust(yaml_ivtoa)

  end function yaml_ivtoa

  !> Yaml Spaced format for Date and Time
  function yaml_date_and_time_toa(values,zone)
    implicit none
    logical, optional, intent(in) :: zone
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_date_and_time_toa
    !local variables
    character(len=*), parameter :: &
         deffmt='i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3'
    logical :: zon
    integer :: zonhrs,zonmin
    integer, dimension(8) :: vals
    character(len=4) :: sgn

    zon=.false.
    if (present(zone)) zon=zone

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    if (zon) then
       zonmin=abs(mod(vals(4),60))
       zonhrs=abs(vals(4)/60)
       if (vals(4) < 0) then
          sgn='" -"'
       else
          sgn='" +"'
       end if
       write(yaml_date_and_time_toa,'('//deffmt//','//sgn//',i2.2,":",i2.2)')&
            vals(1:3),vals(5:8),zonhrs,zonmin

    else
       write(yaml_date_and_time_toa,'('//deffmt//')')vals(1:3),vals(5:8)
    end if

    yaml_date_and_time_toa=yaml_adjust(yaml_date_and_time_toa)

  end function yaml_date_and_time_toa

  !> Yaml Spaced format for Date
  function yaml_date_toa(values)
    implicit none
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_date_toa
    !local variables
    integer, dimension(8) :: vals

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    write(yaml_date_toa,'(i4.4,"-",i2.2,"-",i2.2)')vals(1:3)

    yaml_date_toa=yaml_adjust(yaml_date_toa)

  end function yaml_date_toa

  function yaml_time_toa(values)
    implicit none
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_time_toa
    !local variables
    integer, dimension(8) :: vals

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    write(yaml_time_toa,'(i2.2,":",i2.2,":",i2.2,".",i3.3)')vals(5:8)

    yaml_time_toa=yaml_adjust(yaml_time_toa)

  end function yaml_time_toa

  function yaml_adjust(str)
    implicit none
    character(len=*), intent(in) :: str
    character(len=max_value_length) :: yaml_adjust

    yaml_adjust=adjustl(str)

    !put a space if there is no sign
    if (yaml_adjust(1:1)/='-') then
       call shiftstr(yaml_adjust,1)
    else
       call shiftstr(yaml_adjust,0)
    end if
    

  end function yaml_adjust

  !> Shifts characters in in the string 'str' n positions (positive values
  !! denote a right shift and negative values denote a left shift). Characters
  !! that are shifted off the end are lost. Positions opened up by the shift 
  !! are replaced by spaces.
  !! This routine has been downloaded from the website http://gbenthien.net/strings/index.html
  subroutine shiftstr(str,n)
    implicit none
    integer, intent(in) :: n
    character(len=*), intent(inout) :: str
    !local variables
    integer :: lenstr,nabs

    lenstr=len(str)
    nabs=iabs(n)
    if(nabs>=lenstr) then
       str=repeat(' ',lenstr)
       return
    end if
    if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
    if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
    return

  end subroutine shiftstr

end module yaml_strings
!> Needed to control yaml indentation and to control output on stdout
module yaml_output
  use yaml_strings
  implicit none
  private 

  !yaml events for dump routine
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
  integer, parameter :: tot_max_record_length=95,tot_max_flow_events=500,tab=5,tot_streams=10
  integer :: active_streams=0,default_stream=1

  !parameter of the document
  type :: yaml_stream
     logical :: pp_allowed=.true. !< Pretty printing allowed
     integer :: unit=6 !<unit for the stdout
     integer :: max_record_length=tot_max_record_length
     integer :: flowrite=0 !< Write in flow (0=no -1=start  1=yes)
     integer :: Wall=-1 !< Warning messages of level Wall stop the program (-1 : none)
     integer :: indent=1 !<Blank spaces indentations for Yaml output level identification
     integer :: indent_previous=0 !< indent level prior to flow writing
     integer :: indent_step=2 !< indentation level
     integer :: tabref=40 !> position of tabular in scalar assignment (single column output)
     integer :: icursor=1 !> running position of the cursor on the line
     integer :: itab_active=0 !> number of active tabbings for the line in flowrite
     integer :: itab=0 !> tabbing to have a look on
     integer :: iflowlevel=0 !>levels of flowrite simoultaneously enabled
     integer :: icommentline=0 !> Active if the line being written is a comment
     integer, dimension(tot_max_record_length/tab) :: linetab !>value of the tabbing in the line
     integer :: ievt_flow=0 !>events which track is kept of in the flowrite
     integer, dimension(tot_max_flow_events) :: flow_events !> Set of events in the flow
  end type yaml_stream

  type(yaml_stream), dimension(tot_streams), save :: streams
  integer, dimension(tot_streams) :: stream_units

  interface yaml_map
     module procedure yaml_map,yaml_map_i,yaml_map_f,yaml_map_d,yaml_map_l,yaml_map_iv,yaml_map_dv
  end interface

  public :: yaml_map,yaml_sequence,yaml_new_document,yaml_set_stream,yaml_warning
  public :: yaml_newline,yaml_open_map,yaml_close_map,yaml_stream_attributes
  public :: yaml_open_sequence,yaml_close_sequence,yaml_comment,yaml_toa,yaml_set_default_stream
  public :: yaml_get_default_stream,yaml_date_and_time_toa,yaml_scalar
  
contains

  !> Set the default stream of the module. Return a STREAM_ALREADY_PRESENT errcode if 
  !! the stream has not be initialized
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

  subroutine yaml_get_default_stream(unit)
    implicit none
    integer, intent(out) :: unit

    unit=stream_units(default_stream)

  end subroutine yaml_get_default_stream


  !set all the output from now on to the file indicated by stdout
  subroutine yaml_set_stream(unit,filename,istat,tabbing,record_length)
    implicit none
    integer, optional, intent(in) :: unit,tabbing,record_length
    character(len=*), optional, intent(in) :: filename
    integer, optional, intent(out) :: istat
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
 
  !> print the attributes of the stream at present
  subroutine yaml_stream_attributes(stream_unit,unit,&
       icursor,flowrite,itab_active,iflowlevel,indent,indent_previous,&
       record_length)
    implicit none
    integer, intent(in) ,optional :: unit,stream_unit
    integer, intent(out), optional :: icursor,flowrite,itab_active,iflowlevel
    integer, intent(out), optional :: indent,indent_previous,record_length
    !local variables
    logical :: dump
    integer :: sunt,unt,strm,icursort,flowritet,itab_activet,iflowlevelt,indentt
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
         call yaml_map('Active Tabulars',itab_activet,unit=unt)
         if (itab_activet>0) call yaml_map('Tabular Values',linetab(1:itab_activet),unit=unt)
       call yaml_close_map(unit=unt)
       call yaml_newline(unit=unt)
    end if
  end subroutine yaml_stream_attributes

  subroutine yaml_new_document(unit)
    implicit none
    integer, optional, intent(in) :: unit
    !local variables
    integer :: unt,strm

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    !check all indentation
    if (streams(strm)%indent /= 1) then
       call yaml_warning("Indentation error. Yaml Document has not been closed correctly",unit=stream_units(strm))
       streams(strm)%indent=1
    end if
    call dump(streams(strm),'---',event=DOCUMENT_START)
    !write(stdout,'(3a)')'---'
    streams(strm)%flow_events=NONE
  end subroutine yaml_new_document

  subroutine yaml_warning(message,level,unit)
    implicit none
    integer, optional, intent(in) :: unit
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: level
    !local variables
    integer :: ierr
    integer :: unt,strm

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    call dump(streams(strm),' #WARNING:'//trim(message))
    if (present(level)) then
       if (level <= streams(strm)%Wall) then
          call dump(streams(strm),' Critical warning level reached, aborting...')
          stop
       end if
    end if
  end subroutine yaml_warning

  subroutine yaml_comment(message,advance,unit,hfill,tabbing)
    implicit none
    character(len=1), optional, intent(in) :: hfill
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: unit,tabbing
    character(len=*), intent(in), optional :: advance
    !local variables
    integer :: unt,strm,msg_lgt,tb,ipos
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

    ipos=max(streams(strm)%icursor,streams(strm)%indent)

    msg_lgt=0
    if (present(tabbing)) then
       tb=max(tabbing-ipos-1,1)
       call buffer_string(towrite,len(towrite),repeat(' ',tb),msg_lgt)
       ipos=ipos+tb
    end if

    call buffer_string(towrite,len(towrite),trim(message),msg_lgt)


    if (present(hfill)) then
       call dump(streams(strm),&
            repeat(hfill,&
            max(streams(strm)%max_record_length-ipos-&
            len_trim(message)-3,0))//' '//towrite(1:msg_lgt),&
            advance=adv,event=COMMENT)
    else
       call dump(streams(strm),towrite(1:msg_lgt),advance=adv,event=COMMENT)
    end if

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

  subroutine yaml_open_map(mapname,label,flow,unit)
    use yaml_strings
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

    doflow=streams(strm)%flowrite /=0
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

    if (doflow) call open_flow_level(streams(strm))

    if (doflow .or. msg_lgt==0) then
       adv='no '
    else
       adv='yes'
    end if

    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING_START)

  end subroutine yaml_open_map

  subroutine yaml_close_map(advance,unit)
    implicit none
    integer, optional, intent(in) :: unit
    character(len=*), optional, intent(in) :: advance
    !local variables
    integer :: unt,strm
    character(len=3) :: adv

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

    if (streams(strm)%flowrite /=0) call close_flow_level(streams(strm))
  end subroutine yaml_close_map

!!$  !> Routine useful only for a flow writing
!!$  subroutine yaml_open_sequence(flow,unit)
!!$    implicit none
!!$    logical, optional, intent(in) :: flow
!!$    integer, optional, intent(in) :: unit
!!$    !local variables
!!$    logical :: doflow
!!$    integer :: unt,strm
!!$    character(len=3) :: adv
!!$
!!$    unt=0
!!$    if (present(unit)) unt=unit
!!$    call get_stream(unt,strm)
!!$
!!$    !override it if flow has been already opened
!!$    doflow=streams(strm)%flowrite /=0
!!$    if (present(flow)) doflow=flow .or. doflow
!!$
!!$    if (doflow) then
!!$       adv='no '
!!$    end if
!!$
!!$    if (doflow) then
!!$       call open_flow_level(streams(strm))
!!$    end if
!!$
!!$    call dump(streams(strm),' ',advance=trim(adv),event=SEQUENCE_START)
!!$
!!$  end subroutine yaml_open_sequence

  subroutine yaml_open_sequence(mapname,label,flow,unit)
    use yaml_strings
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

    doflow=streams(strm)%flowrite /=0
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

    if (doflow) call open_flow_level(streams(strm))

    if (doflow .or. msg_lgt==0) then
       adv='no '
    else
       adv='yes'
    end if

    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=SEQUENCE_START)

  end subroutine yaml_open_sequence


  subroutine yaml_close_sequence(advance,unit)
    implicit none
    character(len=*), optional, intent(in) :: advance
    integer, optional, intent(in) :: unit
    !local variables
    integer :: unt,strm
    character(len=3) :: adv

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

    if (streams(strm)%flowrite /=0) call close_flow_level(streams(strm))

  end subroutine yaml_close_sequence

  !> Add a new line in the flow 
  !! this routine has a effect only if a flow writing is active
  subroutine yaml_newline(unit)
    implicit none
    integer, optional, intent(in) :: unit
    !local variables
    integer :: unt,strm

    unt=0
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    if (streams(strm)%icursor > 1) then
    !if (streams(strm)%flowrite /=0 .and. streams(strm)%icursor > 1) then
       call dump(streams(strm),' ',advance='yes',event=NEWLINE)
    end if
  end subroutine yaml_newline

  subroutine yaml_sequence(seqvalue,label,advance,unit)
    use yaml_strings
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

  subroutine yaml_map(mapname,mapvalue,label,advance,unit,fmt)
    use yaml_strings
    implicit none
    character(len=*), intent(in) :: mapname
    character(len=*), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    logical :: cut,redo_line
    integer :: msg_lgt,strm,istream,unt,icut,istr,ierr,msg_lgt_ck
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
    !print *,'ierr',ierr
    if (redo_line) then
       if (streams(strm)%flowrite/=0) then
          call dump(streams(strm),towrite(1:msg_lgt_ck),advance=trim(adv),event=SCALAR)
       else
          if (present(label)) then
             call yaml_open_map(mapname,label=label,unit=unt)
          else
             call yaml_open_map(mapname,unit=unt)
          end if
       end if
!       if (streams(strm)%flowrite/=0) call yaml_newline(unit=unt)
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
          if (streams(strm)%flowrite/=0 .and. .not. cut) call buffer_string(towrite,len(towrite),',',msg_lgt)
          call dump(streams(strm),towrite(1:msg_lgt),advance='yes',event=SCALAR)
          istr=icut
          icut=len_trim(mapvalue)-istr+1
          !print *,'icut',istr,icut,mapvalue(istr:istr+icut-1),cut,istr+icut-1,len_trim(mapvalue)
          msg_lgt=0
       end do cut_line
       if (streams(strm)%flowrite==0) call yaml_close_map(unit=unt)
    end if

  end subroutine yaml_map

  subroutine yaml_map_i(mapname,mapvalue,label,advance,unit,fmt)
    use yaml_strings
    implicit none
    character(len=*), intent(in) :: mapname
    integer, intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,istream,unt
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
    use yaml_strings
    implicit none
    character(len=*), intent(in) :: mapname
    real, intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,istream,unt
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
    use yaml_strings
    implicit none
    character(len=*), intent(in) :: mapname
    real(kind=8), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,istream,unt
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
    use yaml_strings
    implicit none
    character(len=*), intent(in) :: mapname
    logical,  intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,istream,unt
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
    use yaml_strings
    implicit none
    character(len=*), intent(in) :: mapname
    real(kind=8), dimension(:), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,istream,unt
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
  end subroutine yaml_map_dv

  subroutine yaml_map_iv(mapname,mapvalue,label,advance,unit,fmt)
    use yaml_strings
    implicit none
    character(len=*), intent(in) :: mapname
    integer, dimension(:), intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label,advance,fmt
    integer, optional, intent(in) :: unit
    !local variables
    integer :: msg_lgt,strm,istream,unt
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
    use yaml_strings
    implicit none
    type(yaml_stream), intent(inout) :: stream
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: event
    integer, intent(out), optional :: istat
    !local variables
    logical :: ladv,change_line,reset_line,pretty_print,reset_tabbing,comma_postponed
    integer :: lgt,evt,indent_lgt,msg_lgt,shift_lgt,prefix_lgt,iscpos
    integer :: ianchor_pos,tabeff,towrite_lgt
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
    ladv=(stream%flowrite==0)
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
    if (put_comma(stream,evt)) then!stream%icomma==1 .and. stream%flowrite==1) then
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

       if (stream%flowrite==0) then
          call open_indent_level(stream)
       else
          call buffer_string(towrite,len(towrite),' [',msg_lgt)
          !comma has to be written afterwards, if there is a message
          stream%flowrite=-1
       end if

    case(SEQUENCE_END)
       !print *,'here',prefix_lgt,prefix,icomma,flowrite,iflowlevel

       if (stream%flowrite==0) then
          call close_indent_level(stream)
       else
          if (stream%iflowlevel > 1 .and. ladv) then
             call buffer_string(prefix,len(prefix),']',prefix_lgt,back=.true.)
             stream%flowrite=-1
          else
             call buffer_string(prefix,len(prefix),']',prefix_lgt)
          end if
          reset_line=ladv
       end if

    case(MAPPING_START)

       if (stream%flowrite ==0) then
          call open_indent_level(stream)
       else
          !write(stdout,*)'here',prefix,'there',icomma,flowrite,iflowlevel
          call buffer_string(towrite,len(towrite),' {',msg_lgt)
          stream%flowrite=-1
          reset_tabbing=.true.
       end if

       !pretty_print=.true. .and. stream%pp_allowed

    case(MAPPING_END)

       if (stream%flowrite==0) then
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

    case(SEQUENCE_ELEM)

       if (stream%flowrite==0) then
          !lower indent and update prefix
          indent_lgt=indent_lgt-2
          call buffer_string(prefix,len(prefix),'- ',prefix_lgt)
       else
          if (msg_lgt>0) comma_postponed=.false.
       end if

    case(SCALAR)

    case(NEWLINE)
       if (stream%flowrite/=0) then
          !print *,'NEWLINE:',stream%flowrite
          change_line=.true.
          stream%flowrite=-1
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
       call pretty_printing((stream%flowrite/=0),':',towrite,&
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
          write(stream%unit,'(a)',advance=trim(adv))repeat(' ',indent_lgt)//towrite(1:towrite_lgt)
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
    if (stream%flowrite /=0) then
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
      use yaml_strings
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

    !>calculate the reference tabular value
    subroutine closest_tab(ianchor_pos,tabeff)
      implicit none
      integer, intent(in) :: ianchor_pos
      integer, intent(out) :: tabeff

      if (stream%flowrite /= 0) then
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

      if (stream%flowrite==0 .and. stream%icursor==1) then
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
!      else if(stream%flowrite/=-1) then !other conditions have to be added here
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
      !local variables
      integer :: ievt
      put_comma=stream%flowrite/=0 .and. stream%ievt_flow>0

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
      
    subroutine open_flow_level(stream)
      implicit none
      type(yaml_stream), intent(inout) :: stream
      if (stream%flowrite ==0) then
         if (stream%iflowlevel==0) stream%indent_previous=stream%indent
         stream%indent=1
      end if
      stream%iflowlevel=stream%iflowlevel+1
      if (stream%flowrite==0) stream%flowrite=-1 !start to write
    end subroutine open_flow_level

    subroutine close_flow_level(stream)
      implicit none
      type(yaml_stream), intent(inout) :: stream
      !lower the flowlevel
      stream%iflowlevel=stream%iflowlevel-1
      if (stream%iflowlevel==0) then
         stream%indent=stream%indent_previous
         stream%flowrite=0
         !reset the events in the flow
         stream%flow_events=NONE
      else
         stream%indent=1
         stream%flowrite=-1
      end if

    end subroutine close_flow_level

    subroutine open_indent_level(stream)
      implicit none
      type(yaml_stream), intent(inout) :: stream
      stream%indent=stream%indent+stream%indent_step
    end subroutine open_indent_level

    subroutine close_indent_level(stream)
      implicit none
      type(yaml_stream), intent(inout) :: stream
      stream%indent=max(stream%indent-stream%indent_step,0) !to prevent bugs
    end subroutine close_indent_level

    
end module yaml_output
