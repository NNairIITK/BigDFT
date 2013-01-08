!> @file
!! Define the modules (yaml_strings and yaml_output) and the methods to write yaml output
!! yaml: Yet Another Markeup Language (ML for Human)
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Define all yaml strings for output
!! This module must be only used by the module yaml_output
module yaml_strings

  implicit none

  integer :: max_value_length=95 !< Not a parameter in order to be used by C bindings but constant

  character(len=*), parameter :: yaml_int_fmt  = '(i0)'       !< Default format for integer
  character(len=*), parameter :: yaml_real_fmt = '(1pe17.9)' !< Default format for integer
  character(len=*), parameter :: yaml_dble_fmt = '(1pe25.17)' !< Default format for integer

  interface yaml_toa             !< Convert into a character string yaml_toa(xxx,fmt)
     module procedure yaml_itoa,yaml_litoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_dvtoa,yaml_ivtoa,yaml_cvtoa
  end interface

  public ::  yaml_toa, buffer_string, align_message, shiftstr
  private :: yaml_itoa,yaml_litoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_dvtoa,yaml_ivtoa,max_value_length
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


  !> Add the spaces necessary to align the first occurrence of a given anchor
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
       write(yaml_itoa,fmt) i
    else
       write(yaml_itoa,yaml_int_fmt) i
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
       write(yaml_litoa,fmt) i
    else
       write(yaml_litoa,yaml_int_fmt) i
    end if

    yaml_litoa=yaml_adjust(yaml_litoa)

  end function yaml_litoa


  !> Convert float to character
  function yaml_ftoa(f,fmt)
    implicit none
    real, intent(in) :: f
    character(len=max_value_length) :: yaml_ftoa
    character(len=*), optional, intent(in) :: fmt

    yaml_ftoa=repeat(' ',max_value_length)
    if (present(fmt)) then
       write(yaml_ftoa,fmt) f
    else
       write(yaml_ftoa,yaml_real_fmt) f
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
       write(yaml_dtoa,fmt) d
    else
       write(yaml_dtoa,yaml_dble_fmt) d
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

    if (nl > nu) then
       !Special case for size 0 (nl is > nu!)
       yaml_dvtoa(1:2) = '[]'
    else
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
    end if

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

    if (nl > nu) then
       !Special case for size 0 (nl is > nu!)
       yaml_ivtoa(1:2) = '[]'
    else
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
    end if

    yaml_ivtoa=yaml_adjust(yaml_ivtoa)

  end function yaml_ivtoa


  !> Convert vector of characters to a chain of characters
  function yaml_cvtoa(cv,fmt)
    implicit none
    character(len=*), dimension(:), intent(in) :: cv
    character(len=max_value_length) :: yaml_cvtoa
    character(len=*), optional, intent(in) :: fmt
    !local variables
    character(len=max_value_length) :: tmp
    integer :: nl,nu,i,length,pos

    tmp=repeat(' ',max_value_length)
    yaml_cvtoa=tmp

    nl=lbound(cv,1)
    nu=ubound(cv,1)

    if (nl > nu) then
       !Special case for size 0 (nl is > nu!)
       yaml_cvtoa(1:2) = '[]'
    else
       yaml_cvtoa(1:2)='[ '
       pos=3
       do i=nl,nu
          tmp=trim(cv(i))
          length=len(trim(tmp))-1
          if (pos+length > max_value_length) exit
          yaml_cvtoa(pos:pos+length)=tmp(1:length+1)
          if (i < nu) then
             yaml_cvtoa(pos+length+1:pos+length+2)=', '
          else
             yaml_cvtoa(pos+length+1:pos+length+2)=' ]'
          end if
          pos=pos+length+3
       end do
    end if
    yaml_cvtoa=yaml_adjust(yaml_cvtoa)
  end function yaml_cvtoa


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

    !There is no - sign so we skip this step (TD)
    !yaml_date_and_time_toa=yaml_adjust(yaml_date_and_time_toa)

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

