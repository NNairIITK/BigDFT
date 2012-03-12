!> Needed to control yaml indentation and to control output on stdout
module yaml_output

  implicit none
  private 


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
  integer :: icomma=0 !> comma has been already written
  integer :: iflowlevel=0 !>levels of flowrite simoultaneously enabled
  integer, dimension(max_record_length/tab), save :: linetab !>value of the tabbing in the line
  interface yaml_toa
     module procedure yaml_itoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_dvtoa
  end interface

  public :: yaml_toa,yaml_map,yaml_indent_map,yaml_close_indent_map
  public :: yaml_flow_map,yaml_close_flow_map,yaml_flow_newline,yaml_flow_sequence,yaml_close_flow_sequence
  public :: yaml_sequence_element,yaml_close_sequence_element

  contains

    subroutine yaml_new_document()
      implicit none
      !check all indentation
      if (yaml_indent /= 1) then
         call yaml_warning("Indentation error. Yaml Document has not been closed correctly")
         yaml_indent=1
      end if
      write(stdout,'(3a)')'---'
    end subroutine yaml_new_document
  
    subroutine yaml_warning(message,level)
      use module_base, only: MPI_COMM_WORLD
      implicit none
      character(len=*), intent(in) :: message
      integer, optional, intent(in) :: level
      !local variables
      integer :: ierr
      write(stdout,'(1x,a,a)')'#WARNING:',trim(message)
      if (present(level)) then
         if (level <= Wall) then
            write(stdout,'(1x,a)')'Critical warning level reached, aborting...'
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
         end if
      end if
    end subroutine yaml_warning

    subroutine yaml_indent_map(mapname,label)
      implicit none
      character(len=*), intent(in) :: mapname
      character(len=*), optional, intent(in) :: label
      !character(len=3), optional, intent(in) :: verbatim
      !local variables

      if (present(label)) then
         call yaml_map(mapname,label=label)
      else
         call yaml_map(mapname)
      end if
      yaml_indent=yaml_indent+yaml_level

    end subroutine yaml_indent_map

    subroutine yaml_close_indent_map()
      yaml_indent=max(yaml_indent-yaml_level,0)
    end subroutine yaml_close_indent_map
      
    !> Open a hash table written in flow format
    subroutine yaml_flow_map(mapname,label)
      implicit none
      character(len=*), optional, intent(in) :: mapname
      character(len=*), optional, intent(in) :: label
      !local variables
      character(len=500) :: line


      if (present(mapname))then
         if (present(label)) then
            call yaml_map(mapname,label=label,advance='no')
         else
            call yaml_map(mapname,advance='no')
         end if
      end if
      write(stdout,'(a)',advance='no')' {'
      icursor=icursor+2

      flowrite=-1 !start to write
      if (iflowlevel==0) yaml_indent_previous=yaml_indent
      iflowlevel=iflowlevel+1

      yaml_indent=1

    end subroutine yaml_flow_map

    subroutine yaml_close_flow_map(advance)
      implicit none
      character(len=*), optional, intent(in) :: advance
      !terminate mapping
      if (present(advance)) then
         write(stdout,'(a)',advance=advance)'}'
         icursor=icursor+1
         if (advance=='yes') then
            icursor=1
            flowrite=0
            iflowlevel=iflowlevel-1
            if (iflowlevel==0) then
               yaml_indent=yaml_indent_previous
            else
               yaml_indent=1
            end if
            itab_active=0
            itab=0
         end if
      else
         write(stdout,'(a)')'}'
         icursor=1
         flowrite=0
         iflowlevel=iflowlevel-1
         if (iflowlevel==0) then
            yaml_indent=yaml_indent_previous
         else
            yaml_indent=1
         end if
         itab_active=0
         itab=0
      end if
    end subroutine yaml_close_flow_map

    !> Open a sequence written in flow format
    subroutine yaml_flow_sequence()
      implicit none
      !local variables
      character(len=500) :: line

      write(stdout,'(a)',advance='no')' ['
      icursor=icursor+2

      if (flowrite ==0) then
         if (iflowlevel==0) yaml_indent_previous=yaml_indent
         iflowlevel=iflowlevel+1
         yaml_indent=1
      end if

      flowrite=-1 !start to write
    end subroutine yaml_flow_sequence

    subroutine yaml_close_flow_sequence(advance)
      implicit none
      character(len=*), optional, intent(in) :: advance
      !terminate mapping
      if (present(advance)) then
         write(stdout,'(a)',advance=advance)']'
         icursor=icursor+1
         if (advance=='yes') then
            icursor=1
            flowrite=0
            iflowlevel=iflowlevel-1
            if (iflowlevel==0) then
               yaml_indent=yaml_indent_previous
            else
               yaml_indent=1
            end if
            itab_active=0
            itab=0
         end if
      else
         write(stdout,'(a)')']'
         icursor=1
         flowrite=0
         iflowlevel=iflowlevel-1
         if (iflowlevel==0) then
            yaml_indent=yaml_indent_previous
         else
            yaml_indent=1
         end if
         itab_active=0
         itab=0
      end if
    end subroutine yaml_close_flow_sequence

    subroutine yaml_sequence_element(advance)
      implicit none
      character(len=*), optional, intent(in) :: advance
      !needed only for indented writings
      if (flowrite==0) then
         if(present(advance)) then
            write(stdout,'(a)',advance=advance)repeat(' ',yaml_indent)//'-'
            if (advance=='no')icursor=icursor+yaml_indent
         else
            write(stdout,'(a)')repeat(' ',yaml_indent)//'-'
         end if
         yaml_indent=yaml_indent+yaml_level
      end if
      !comma should be placed after
      icomma=0
    end subroutine yaml_sequence_element
    
    subroutine yaml_close_sequence_element(advance)
      implicit none
      character(len=*), optional, intent(in) :: advance
      if (flowrite==0 .and. iflowlevel==0) then
         yaml_indent=max(yaml_indent-yaml_level,0)
      else
         if(present(advance)) then
            if (advance=='no') then
               if (icomma==0) write(stdout,'(a)',advance='no')','
               icursor=icursor+1
               icomma=1
            else
               write(stdout,'(a)')','
               icursor=1
               itab=0
            end if
         else
            if (icomma==0) write(stdout,'(a)',advance='no')','
            icursor=icursor+1
            icomma=1
         end if
      end if
    end subroutine yaml_close_sequence_element
    
    
    !> Add a new line in the flow 
    !! this routine has a effect only if flowrite is active
    subroutine yaml_flow_newline()
      implicit none
      if (flowrite ==-1) then !flowrite had just been opened
         write(stdout,'(a)')' '
      else if (flowrite==1) then !just close the last field
         if (icomma==0) then
            write(stdout,'(a)')', '
         else
            write(stdout,'(a)')' '
         end if
         flowrite=-1 !restart from line
         icomma=0
      end if
      !reset tabbing
      itab_active=0
      itab=0
      icursor=1
    end subroutine yaml_flow_newline

    !> fill the hash table value, if present.
    !! This should be the only routine which writes hash table elements 
    subroutine yaml_map(mapname,mapvalue,label,advance)
      implicit none
      character(len=*), intent(in) :: mapname
      character(len=*), optional, intent(in) :: mapvalue,label,advance
      !local variables
      integer :: ipos,lgt,tabeff,ish
      character(len=2) :: fmt
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

    end subroutine yaml_map

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
