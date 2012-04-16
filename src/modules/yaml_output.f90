!> Needed to control yaml indentation and to control output on stdout
module yaml_output

  implicit none
  private 

  !yaml events for dump routine
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

  integer, parameter :: stdout=70 !<unit for the stdout
  integer, parameter :: max_value_length=50,max_record_length=95,tab=5 !> tabbing size
  integer :: flowrite=0 !< Write in flow (0=no -1=start  1=yes)
  integer :: Wall=-1 !< Warning messages of level Wall stop the program (-1 : none)
  integer :: yaml_indent=1 !<Blank spaces indentations for Yaml output level identification
  integer :: yaml_indent_previous=0 !< indent level prior to flow writing
  integer :: yaml_level=2 !< indentation level
  integer :: tabmain=40 !> position of tabular in scalar assignment (single column output)
  integer :: icursor=1 !> running position of the cursor on the line
  integer :: itab_active=0 !> number of active tabbings for the line in flowrite
  integer :: itab=0 !> tabbing to have a look on
  integer :: icomma=0 !> if zero, no comma has to be placed if the flow continues
  integer :: iflowlevel=0 !>levels of flowrite simoultaneously enabled
  integer :: icommentline=0 !> Active if the line being written is a comment
  integer, dimension(max_record_length/tab), save :: linetab !>value of the tabbing in the line
  interface yaml_toa
     module procedure yaml_itoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_dvtoa
  end interface

  public :: yaml_toa,yaml_map,yaml_indent_map,yaml_close_indent_map
  public :: yaml_flow_map,yaml_close_flow_map,yaml_flow_newline,yaml_flow_sequence,yaml_close_flow_sequence
  public :: yaml_sequence_element,yaml_close_sequence_element,yaml_comment

contains

  subroutine dump(message,fmt,advance,event)
    implicit none
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: event
    !local variables
    logical :: ladv,change_line,reset_line,pretty_print,reset_tabbing
    integer :: lgt,evt,indent_lgt,msg_lgt,shift_lgt,prefix_lgt,iscpos
    integer :: ianchor_pos,tabeff,towrite_lgt
    character(len=3) :: adv
    character(len=5) :: prefix
    character(len=max_record_length) :: towrite


    !some consistency check
    if (icomma /= 0 .and. flowrite ==0) then
       !comma should not be there or flowrite is not correct
       !switch comma out
       icomma=0
    end if


    if (present(event)) then
       evt=event
    else !default event: scalar value
       evt=SCALAR
    end if

    !decide whether to write advanced or not
    !decide if advanced output or not
    ladv=(flowrite==0)
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
    if (icursor > 1) then
       indent_lgt=0
    else !other conditions have to be added here
       indent_lgt=max(yaml_indent,0) !to prevent bugs
    end if


    !calculate the number of objects to be written before
    !these objects should go to the active line in case of a new line
    prefix_lgt=0
    !initialize it
    prefix=repeat(' ',len(prefix))

    !check if comma have to be written
    if (icomma==1 .and. flowrite==1) then
       call buffer_string(prefix,len(prefix),', ',prefix_lgt)
    end if
    !write comma next time
    if (icomma==1 .and. flowrite==-1) then
       flowrite=1
    end if
    !string length, and message body
    !initialize it
    towrite=repeat(' ',len(towrite))
    msg_lgt=0
    !a empty message is not written
    if (len_trim(message) > 0) &
         call buffer_string(towrite,len(towrite),message,msg_lgt)

    !no pretty printing by default
    pretty_print=.false.
    shift_lgt=0

    !reset_tabbing is disabled
    reset_tabbing=.false.

    !set module variables according to the event
    select case(evt)
    case(SEQUENCE_START)

       if (flowrite ==0) then
          call open_indent_level()
       else
          call buffer_string(towrite,len(towrite),' [',msg_lgt)
          !comma has to be written afterwards
          icomma=1
          flowrite=-1
       end if

    case(SEQUENCE_END)
       !print *,'here',prefix_lgt,prefix,icomma,flowrite,iflowlevel

       if (flowrite==0) then
          call close_indent_level()
       else
          if (iflowlevel > 1 .and. ladv) then
             call buffer_string(prefix,len(prefix),']',prefix_lgt,back=.true.)
             flowrite=-1
             !the comma will be erased by the end of the sequence
          else if (icomma==1 .and. prefix_lgt>0) then
             prefix_lgt=prefix_lgt-1
             prefix(prefix_lgt:prefix_lgt)=']'
          else
             call buffer_string(prefix,len(prefix),']',prefix_lgt)
          end if

          if (iflowlevel==1) then
             icomma=0
          end if
          reset_line=ladv
       end if

    case(MAPPING_START)

       if (flowrite ==0) then
          call open_indent_level()
       else
          call buffer_string(towrite,len(towrite),' {',msg_lgt)
          !comma has to be written afterwards
          icomma=1
          flowrite=-1
          reset_tabbing=.true.
       end if

       pretty_print=.true.

    case(MAPPING_END)

       if (flowrite==0) then
          call close_indent_level()
       else
          if (iflowlevel > 1 .and. ladv) then
             call buffer_string(prefix,len(prefix),'}',prefix_lgt,back=.true.)
             flowrite=-1
             reset_line=.true.
          !the comma will be erased by the end of the mapping
          else if (icomma==1 .and. prefix_lgt>0) then
             prefix_lgt=prefix_lgt-1
             prefix(prefix_lgt:prefix_lgt)='}'
          !terminate nonetheless the mapping
          else 
             call buffer_string(prefix,len(prefix),'}',prefix_lgt)
          end if

          if (iflowlevel==1) then
             icomma=0
          end if
          reset_line=ladv
       end if

    case(COMMENT)
       if (icommentline==0) then !no comment active
          call buffer_string(prefix,len(prefix),' #',prefix_lgt)
       end if
       if (.not. ladv) then
          icommentline=1
       else
          reset_line=.true.
       end if
    case(MAPPING)

       pretty_print=.true.

    case(SEQUENCE_ELEM)

       if (flowrite==0) then
          !lower indent and update prefix
          indent_lgt=indent_lgt-2
          call buffer_string(prefix,len(prefix),'- ',prefix_lgt)
       end if

    end select

    !adjust the towrite string to match with the closest tabular
    if (pretty_print) then
       call pretty_printing((flowrite/=0),':',towrite,icursor,indent_lgt,prefix_lgt,&
            msg_lgt,max_record_length,shift_lgt,change_line)
    end if

    !standard writing,
    if (change_line) then
       !first write prefix, if needed
       if (prefix_lgt>0) then
          write(stdout,'(a)')prefix(1:prefix_lgt)
       else if (msg_lgt >0) then
          !change line
          write(stdout,*)
       end if
       icursor=1
       towrite_lgt=msg_lgt+shift_lgt
    else
       call shiftstr(towrite,prefix_lgt)
       if (prefix_lgt > 0)towrite(1:prefix_lgt)=prefix(1:prefix_lgt)
       towrite_lgt=prefix_lgt+msg_lgt+shift_lgt
    end if
!print *,'adv',trim(adv),towrite_lgt,icursor,change_line,msg_lgt
    if (present(fmt)) then
       write(stdout,fmt,advance=trim(adv))repeat(' ',indent_lgt)//towrite(1:towrite_lgt)
    else if (towrite_lgt > 0) then
       write(stdout,'(a)',advance=trim(adv))repeat(' ',indent_lgt)//towrite(1:towrite_lgt)
    end if

    !if advancing i/o cursor is again one
    if (ladv) then
       icursor=1
    else
       !cursor after writing
       icursor=icursor+indent_lgt+towrite_lgt
    end if

    if (reset_tabbing) then
       itab_active=0
       itab=0
    end if

    if (reset_line) call carriage_return()

  end subroutine dump

  subroutine pretty_printing(rigid,anchor,message,icursor,&
       indent_lgt,prefix_lgt,msg_lgt,max_lgt,shift_lgt,change_line)
    implicit none
    logical, intent(in) :: rigid
    integer, intent(in) :: icursor,indent_lgt,prefix_lgt,msg_lgt,max_lgt
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
       if (itab==itab_active .and. itab > 1) itab_active=itab_active-1
       itab=1
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

      if (flowrite /= 0) then
         !first check that the tabbing is already done, otherwise add another tab
         if (itab<itab_active) then
            !realign the value to the tabbing
            do 
               if (ianchor_pos <= linetab(itab) .or. itab==itab_active) exit
               itab=itab+1
            end do
         end if

         if (itab<itab_active) then
            tabeff=linetab(itab)
         else
            tabeff=ianchor_pos
            itab=itab+1
            itab_active=itab_active+1
            linetab(itab_active)=tabeff
         end if
      else
         !for the moment do not check compatibility of the line
         tabeff=max(tabmain,ianchor_pos)
      end if
    end subroutine closest_tab


    !> Add a buffer to a string and increase its length
    subroutine buffer_string(string,string_lgt,buffer,string_pos,back)
      implicit none
      integer, intent(in) :: string_lgt
      integer, intent(inout) :: string_pos
      character(len=*), intent(in) :: buffer
      character(len=string_lgt), intent(inout) :: string
      logical, optional, intent(in) :: back
      !local variables
      integer :: lgt_add

      lgt_add=len(buffer)
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
    

    !> Reset the line control quantities, and reset the indentation
    subroutine carriage_return()
      implicit none
      !if a yaml_comment is called put the has in front
      icommentline=0
      !beginining of the line
      icursor=1
      !no tabbing decided yet
      itab_active=0
      itab=0
      !all needed commas are placed in the previous line
      !icomma=0
    end subroutine carriage_return
      
    subroutine yaml_new_document()
      implicit none
      !check all indentation
      if (yaml_indent /= 1) then
         call yaml_warning("Indentation error. Yaml Document has not been closed correctly")
         yaml_indent=1
      end if
      call dump('---',fmt='(3a)')
      !write(stdout,'(3a)')'---'
    end subroutine yaml_new_document
  
    subroutine yaml_warning(message,level)
      use module_base, only: MPI_COMM_WORLD
      implicit none
      character(len=*), intent(in) :: message
      integer, optional, intent(in) :: level
      !local variables
      integer :: ierr
      call dump('#WARNING:'//trim(message),fmt='(1x,a)')
      !write(stdout,'(1x,a,a)')'#WARNING:',trim(message)
      if (present(level)) then
         if (level <= Wall) then
            call dump('Critical warning level reached, aborting...',fmt='(1x,a)')
            !write(stdout,'(1x,a)')'Critical warning level reached, aborting...'
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
         end if
      end if
    end subroutine yaml_warning

    subroutine yaml_comment(message,advance)
      implicit none
      character(len=*), intent(in) :: message
      character(len=*), intent(in), optional :: advance

      !comment to be written
      if (present(advance)) then
         call dump(trim(message),advance=advance,event=COMMENT)
      else
         call dump(trim(message),advance='yes',event=COMMENT)
      end if

!!$      if (icommentline == 0) then !no comment active
!!$         if (present(advance)) then
!!$            write(stdout,'(a)',advance=advance)' #'//trim(message)
!!$            if (advance=='no') then
!!$               icursor=icursor+len(trim(message))+2
!!$               icommentline=1
!!$            else
!!$               call carriage_return()
!!$            end if
!!$         else
!!$            write(stdout,'(a)')' #'//trim(message)
!!$            call carriage_return()
!!$         end if
!!$      else
!!$         if (present(advance)) then
!!$            write(stdout,'(a)',advance=advance)trim(message)
!!$            if (advance=='no') then
!!$               icursor=icursor+len(trim(message))
!!$               icommentline=1
!!$            else
!!$               call carriage_return()
!!$            end if
!!$         else
!!$            write(stdout,'(a)')trim(message)
!!$            call carriage_return()
!!$         end if
!!$      end if
    end subroutine yaml_comment

    !> Open a indented map
    subroutine yaml_indent_map(mapname,label)
      implicit none
      character(len=*), intent(in) :: mapname
      character(len=*), optional, intent(in) :: label
      !character(len=3), optional, intent(in) :: verbatim
      !local variables       

      !tentative solution
      if (present(label)) then
         call dump(mapname//': &'//label,advance='yes',event=MAPPING_START)
      else
         call dump(mapname//':',advance='yes',event=MAPPING_START)
      end if
      
      
!!$      if (present(label)) then
!!$         call yaml_map(mapname,label=label)
!!$      else
!!$         call yaml_map(mapname)
!!$      end if
!!$      yaml_indent=yaml_indent+yaml_level

    end subroutine yaml_indent_map

    !Adjust the indentation for a indented map
    subroutine yaml_close_indent_map()

      call dump(' ',advance='yes',event=MAPPING_END)

!!$      call close_indent_level()
      !yaml_indent=yaml_indent-yaml_level
    end subroutine yaml_close_indent_map
      
    subroutine open_flow_level()
      implicit none
      if (flowrite ==0) then
         if (iflowlevel==0) yaml_indent_previous=yaml_indent
         yaml_indent=1
      end if
      iflowlevel=iflowlevel+1
      flowrite=-1 !start to write
    end subroutine open_flow_level

    subroutine close_flow_level()
      implicit none
      
      !lower the flowlevel
      iflowlevel=iflowlevel-1
      if (iflowlevel==0) then
         yaml_indent=yaml_indent_previous
         flowrite=0
      else
         yaml_indent=1
         flowrite=-1
      end if

    end subroutine close_flow_level

    !> Open a hash table written in flow format
    subroutine yaml_flow_map(mapname,label)
      implicit none
      character(len=*), optional, intent(in) :: mapname
      character(len=*), optional, intent(in) :: label
      !local variables

      call open_flow_level()
      if (present(mapname)) then
         if (present(label)) then
            call dump(mapname//': &'//label,advance='no',event=MAPPING_START)
         else
            call dump(mapname//':',advance='no',event=MAPPING_START)
         end if
      else
         call dump(' ',advance='no',event=MAPPING_START)
      end if

!!$
!!$      if (present(mapname))then
!!$         if (present(label)) then
!!$            call yaml_map(mapname,label=label,advance='no')
!!$         else
!!$            call yaml_map(mapname,advance='no')
!!$         end if
!!$      end if
!!$      write(stdout,'(a)',advance='no')' {'
!!$      icursor=icursor+2


    end subroutine yaml_flow_map

    subroutine yaml_close_flow_map(advance)
      implicit none
      character(len=*), optional, intent(in) :: advance
      !terminate mapping
      if (present(advance)) then
         call dump(' ',advance=advance,event=MAPPING_END)
      else
         call dump(' ',advance='yes',event=MAPPING_END)
      end if
      call close_flow_level()

!!$      write(stdout,'(a)',advance='no')'}'
!!$      icursor=icursor+1
!!$      call close_flow_level()
!!$      if (iflowlevel==0) write(stdout,'(a)')' '
!!$      if (present(advance)) then
!!$         if (advance=='yes') call yaml_flow_newline()
!!$      else
!!$         call yaml_flow_newline()
!!$      end if
      !write(stdout,*)'debug',iflowlevel,flowrite
    end subroutine yaml_close_flow_map

    !> Open a sequence written in flow format
    subroutine yaml_flow_sequence()
      implicit none
      !local variables
      call open_flow_level()
      call dump(' ',advance='no',event=SEQUENCE_START)

!!$      write(stdout,'(a)',advance='no')' ['
!!$      icursor=icursor+2
!!$      if (flowrite ==0) then
!!$         if (iflowlevel==0) yaml_indent_previous=yaml_indent
!!$         yaml_indent=1
!!$      end if
!!$      iflowlevel=iflowlevel+1
!!$      flowrite=-1 !start to write
!!$      call open_flow_level()
    end subroutine yaml_flow_sequence

    subroutine yaml_close_flow_sequence(advance)
      implicit none
      character(len=*), optional, intent(in) :: advance
      !local variables
      character(len=3) :: adv

      adv='yes' !default value
      if (present(advance)) adv=advance
      call dump(' ',advance=trim(adv),event=SEQUENCE_END)

      call close_flow_level()
!!$      !terminate mapping
!!$      if (present(advance)) then
!!$         write(stdout,'(a)',advance=advance)']'
!!$         icursor=icursor+1
!!$         call close_flow_level()
!!$         if (advance=='yes') then
!!$            call carriage_return()
!!$         end if
!!$      else
!!$         write(stdout,'(a)')']'
!!$         call close_flow_level()
!!$         call carriage_return()
!!$      end if
!!$      !no commas at the end
!!$      if (iflowlevel==0) icomma=0
    end subroutine yaml_close_flow_sequence

    !> Add a new line in the flow 
    !! this routine has a effect only if a flow writing is active
    subroutine yaml_flow_newline()
      implicit none
      if (flowrite ==-1) then !flowrite had just been opened
         write(stdout,'(a)')' '
      else if (flowrite==1) then !just close the last field
         if (icomma==1) then
            write(stdout,'(a)')', '
         else
            write(stdout,'(a)')' '
         end if
         flowrite=-1 !restart from line
         icomma=0
      end if
      !reset tabbing
      call carriage_return()
    end subroutine yaml_flow_newline

    subroutine open_indent_level()
      implicit none
      yaml_indent=yaml_indent+yaml_level
    end subroutine open_indent_level

    subroutine close_indent_level()
      yaml_indent=max(yaml_indent-yaml_level,0) !to prevent bugs
    end subroutine close_indent_level

    subroutine yaml_sequence_element(label,advance)
      implicit none
      character(len=*), optional, intent(in) :: label,advance
      !needed only for indented writings
      if (flowrite==0) then
         if(present(advance)) then
            if (present(label)) then
               write(stdout,'(a)',advance=advance)repeat(' ',yaml_indent)//'- &'//&
                    trim(adjustl(label))
               if (advance=='no')icursor=icursor+yaml_indent+1+len(trim(adjustl(label)))
            else
               write(stdout,'(a)',advance=advance)repeat(' ',yaml_indent)//'-'
               if (advance=='no')icursor=icursor+yaml_indent+1
            end if
         else
            if (present(label)) then
               write(stdout,'(a)')repeat(' ',yaml_indent)//'- &'//&
                    trim(adjustl(label))
            else
               write(stdout,'(a)')repeat(' ',yaml_indent)//'-'
            end if
         end if
         call open_indent_level()
         !yaml_indent=yaml_indent+yaml_level
      end if
      !comma should be placed after
      icomma=1
    end subroutine yaml_sequence_element
    
    subroutine yaml_close_sequence_element(advance)
      implicit none
      character(len=*), optional, intent(in) :: advance
      if (flowrite==0 .and. iflowlevel==0) then
         call close_indent_level()
         !yaml_indent=max(yaml_indent-yaml_level,0)
      else
         if(present(advance)) then
            if (advance=='no') then
               if (icomma==1) write(stdout,'(a)',advance='no')','
               icursor=icursor+1
               icomma=0
            else
               write(stdout,'(a)')','
               icursor=1
               itab=0
               icomma=0
            end if
         else
            if (icomma==1) write(stdout,'(a)',advance='no')','
            icursor=icursor+1
            icomma=0
         end if
      end if
    end subroutine yaml_close_sequence_element
    
    subroutine yaml_map(mapname,mapvalue,label,advance)
      implicit none
      character(len=*), intent(in) :: mapname
      character(len=*), optional, intent(in) :: label,mapvalue,advance
      !local variables
      integer :: msg_lgt
      character(len=3) :: adv
      character(len=max_record_length) :: towrite

      adv='def' !default value
      if (present(advance)) adv=advance
      
      msg_lgt=0
      !put the message
      call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
      !put the semicolon
      call buffer_string(towrite,len(towrite),':',msg_lgt)
      !put the optional name
      if (present(label)) then
         call buffer_string(towrite,len(towrite),' &',msg_lgt)
         call buffer_string(towrite,len(towrite),trim(label),msg_lgt)
      end if
      !put the value
      if (present(mapvalue)) &
           call buffer_string(towrite,len(towrite),trim(mapvalue),msg_lgt)

      call dump(towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
    end subroutine yaml_map

    
    !> fill the hash table value, if present.
    !! This should be the only routine which writes hash table elements 
    subroutine yaml_map_old(mapname,mapvalue,label,advance)
      implicit none
      character(len=*), intent(in) :: mapname
      character(len=*), optional, intent(in) :: mapvalue,label,advance
      !local variables
      integer :: ipos,lgt,tabeff,ish
      character(len=max_record_length) :: towrite

      !reinitialize itab if at beginning of the line
      if (flowrite /=0 .and. icursor==1) itab=0
      
      !write spaces or indentation
      if (icursor==1 .and. flowrite==0) then
         towrite(1:yaml_indent)=repeat(' ',yaml_indent)
         ipos=yaml_indent+1
      else
         !if (flowrite==0) then
         !   towrite(1:1)=' '
         !   ipos=2
         !else
            ipos=1
         !end if
      end if
      !then write name 
      lgt=len(trim(mapname))
      towrite(ipos:ipos+lgt-1)=mapname(1:lgt)
      ipos=ipos+lgt
      !put the semicolon here if name is absent
      if (.not. present(mapvalue)) then
         towrite(ipos:ipos)=':'
         ipos=ipos+1
         !then label (if present)
         if (present(label)) then
            lgt=len(trim(label))
            towrite(ipos:ipos+1)=' &'
            towrite(ipos+2:ipos+lgt+1)=label(1:lgt)
            ipos=ipos+lgt+2
         end if
      end if
      
      !then the value, possibly aligned
      if (present(mapvalue)) then
         !align the value only if flow is desactivated
         !otherwise "magnetize" the value to the closest multiple of tab
         if (flowrite /= 0) then
            !first check that the tabbing is already done, otherwise add another tab
            if (itab<itab_active) then
               !realign the value to the tabbing
               do 
                  if (icursor+ipos-1 <= linetab(itab) .or. itab==itab_active) exit
                  itab=itab+1
               end do
            end if
            if (itab<itab_active) then
               tabeff=linetab(itab)
               !itab=itab+1
            else
               tabeff=ipos+icursor-1!((ipos+icursor-2)/tab+1)*tab
               itab=itab+1
               itab_active=itab_active+1
               linetab(itab_active)=tabeff
!write(stdout,*)mapname,tabeff,ipos
!print *,'there',(ipos+icursor-1),tabeff
            end if
            lgt=tabeff-icursor-ipos+1!-1
            !add tabbing handling
            !shift the value at right 
            if (lgt /= 0) call shiftstr(towrite,lgt)
            ipos=ipos+lgt
!print *,'here',mapname,'a',towrite,'b',lgt,tabeff,(ipos+icursor-1),tab
         else
            lgt=max(tabmain-(ipos+icursor-1),1)
            !put the spaces
            towrite(ipos:ipos+lgt-1)=repeat(' ',lgt)
            ipos=ipos+lgt
         end if
         
         !add the semicolon
         towrite(ipos:ipos)=':'
         ipos=ipos+1

         !then label (if present)
         if (present(label)) then
            lgt=len(trim(label))
            towrite(ipos:ipos+1)=' &'
            towrite(ipos+2:ipos+lgt+1)=label(1:lgt)
            ipos=ipos+lgt+2
         end if

         !add the value
         lgt=len(trim(mapvalue))
         towrite(ipos:ipos+lgt-1)=mapvalue(1:lgt)
         ipos=ipos+lgt
      end if

      !decide whether to write advanced or not
      if (flowrite /= 0) then
         if (flowrite == 1) then
            write(stdout,'(a)',advance='no')', ' !separate previous field
            icursor=icursor+2
         end if
         !control if there is enough space in the record
         if (ipos+icursor-1 < max_record_length) then
            !write(stdout,*)'A',ipos,icursor,max_record_length,ipos+icursor-1
         else
            !close the line and carriage return
            !write(stdout,*)'B',ipos,icursor,max_record_length,ipos+icursor-1,mapname
            write(stdout,'(a)')' '
            icursor=1
            itab=1
            !realign the value to the tabbing
            do 
               if (index(towrite,':') < linetab(itab)) exit
               itab=itab+1
            end do
            !itab=itab-1
            ish=linetab(itab)-index(towrite,':')
            !print *,'a',mapname,ish,ipos,icursor,linetab(1:itab_active),itab_active
            !shift the string if it fits in the line
            if (ipos+ish-1 < max_record_length) then
               call shiftstr(towrite,ish)
               ipos=ipos+ish
            end if
         end if
         write(stdout,'(a)',advance='no')towrite(1:ipos-1)
         icursor=icursor+ipos-1
         !first time writing done
         if (flowrite==-1) flowrite=1
      else
         if (present(advance)) then
            write(stdout,'(a)',advance=advance)towrite(1:ipos-1)
            if (trim(advance)=='no') icursor=icursor+ipos
         else
            write(stdout,'(a)')towrite(1:ipos-1)
         end if
         icursor=1
      end if

    end subroutine yaml_map_old

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
         write(yaml_ltoa,'(l1)')l
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

    function yaml_adjust(str)
      implicit none
      character(len=*), intent(in) :: str
      character(len=max_value_length) :: yaml_adjust

      yaml_adjust=adjustl(str)

      !put a space if there is no sign
      if (yaml_adjust(1:1)/='-') then
         call shiftstr(yaml_adjust,2)
      else
         call shiftstr(yaml_adjust,1)
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
    
end module yaml_output
