!> @file
!! Define the modules (yaml_strings and yaml_output) and the methods to write yaml output
!! yaml: Yet Another Markeup Language (ML for Human)
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining all yaml strings for output.
!! This module must be only used by the module yaml_output
module yaml_strings

  implicit none

  private

  integer :: max_value_length=95 !< Not a parameter in order to be used by C bindings but constant

  character(len=*), parameter :: yaml_int_fmt  = '(i0)'       !< Default format for integer
  character(len=*), parameter :: yaml_real_fmt = '(1pe18.9)' !< Default format for single
  character(len=*), parameter :: yaml_dble_fmt = '(1pg26.17e3)'!'(1pe25.17)' !< Default format for double
  character(len=*), parameter :: yaml_char_fmt = '(a)' !< Default format for strings

  interface yaml_toa             !< Convert into a character string yaml_toa(xxx,fmt)
     module procedure yaml_itoa,yaml_litoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_ctoa
     module procedure yaml_dvtoa,yaml_ivtoa,yaml_cvtoa,yaml_ztoa,yaml_zvtoa,yaml_lvtoa,yaml_rvtoa
  end interface

  interface cnv_fmt
     module procedure fmt_i,fmt_r,fmt_d,fmt_a,fmt_li
  end interface

  !Public routines
  public ::  yaml_toa, buffer_string, align_message, shiftstr,yaml_date_toa
  public :: yaml_date_and_time_toa,yaml_time_toa,is_atoi,is_atof,is_atol

contains

  pure function fmt_li(data)
    implicit none
    integer(kind=8), intent(in) :: data
    character(len=len(yaml_int_fmt)) :: fmt_li
    fmt_li=yaml_int_fmt
  end function fmt_li

  pure function fmt_i(data)
    implicit none
    integer, intent(in) :: data
    character(len=len(yaml_int_fmt)) :: fmt_i
    fmt_i=yaml_int_fmt
  end function fmt_i

  pure function fmt_r(data)
    implicit none
    real, intent(in) :: data
    character(len=len(yaml_real_fmt)) :: fmt_r
    fmt_r=yaml_real_fmt
  end function fmt_r

  pure function fmt_d(data)
    implicit none
    double precision, intent(in) :: data
    character(len=len(yaml_dble_fmt)) :: fmt_d
    fmt_d=yaml_dble_fmt
  end function fmt_d

  pure function fmt_a(data)
    implicit none
    character(len=*), intent(in) :: data
    character(len=len(yaml_char_fmt)) :: fmt_a
    fmt_a=yaml_char_fmt
  end function fmt_a


  !write the strings as they were written by write
  pure subroutine string_assignment(stra,strb)
    implicit none
    character(len=*), intent(out) :: stra
    character(len=*), intent(in) :: strb
    !local variables
    integer :: i

    stra=repeat(' ',len(stra))
    
    do i=1,min(len(stra),len(strb))
       stra(i:i)=strb(i:i)
    end do
    
  end subroutine string_assignment

  !> Add a buffer to a string and increase its length
  pure subroutine buffer_string(string,string_lgt,buffer,string_pos,back,istat)
    implicit none
    integer, intent(in) :: string_lgt                   !< Length of the string towrite
    character(len=string_lgt), intent(inout) :: string  !< String towrite
    integer, intent(inout) :: string_pos                !< Position to add buffer into string and for the next.
    character(len=*), intent(in) :: buffer              !< Buffer to add
    logical, optional, intent(in) :: back               !< Add string from the end
    integer, optional, intent(out) :: istat             !< Error status (if present otherwise stops if error)
    !local variables
    integer :: lgt_add

    if (present(istat)) istat=0 !no errors

    lgt_add=len(buffer)
    !do not copy strings which are too long if istat is present
    if (lgt_add+string_pos > string_lgt) then
       !try to eliminate trailing spaces
       lgt_add=len_trim(buffer)
       if (present(istat)) then
          istat=-1
          return
       else if (lgt_add+string_pos > string_lgt) then
!          write(*,*)'#ERROR (buffer string): string too long'
!          write(*,*)'#Initial String: ',string(1:string_pos)
!          write(*,*)'#Buffer: ',trim(buffer)
!          write(*,*)'#String position: ',string_pos
!          write(*,*)'#Length of Buffer: ',lgt_add
!          write(*,*)'#String limit: ',string_lgt
          lgt_add=string_lgt-string_pos-1
!          write(*,*)'#Buffer shortened into: ',buffer(1:lgt_add)
          !stop 
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
  pure function yaml_itoa(data,fmt) result(str)
    implicit none
    integer, intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_itoa

  !> Convert longinteger to character
  pure function yaml_litoa(data,fmt) result(str)
    implicit none
    integer(kind=8), intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_litoa

  !> Convert float to character
  pure function yaml_ftoa(data,fmt) result(str)
    implicit none
    real, intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_ftoa

  !> Convert double to character
  pure function yaml_dtoa(data,fmt) result(str)
    implicit none
    real(kind=8), intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_dtoa

  !> character to character, only for genericity
  pure function yaml_ctoa(d,fmt)
    implicit none
    character(len=*), intent(in) :: d
    character(len=max_value_length) :: yaml_ctoa
    character(len=*), optional, intent(in) :: fmt

    yaml_ctoa(1:max_value_length)=trim(d)

  end function yaml_ctoa

  !> Convert double complex to character
  !! use python notation for yaml complex
  pure function yaml_ztoa(z,fmt)
    implicit none
    double complex, intent(in) :: z
    character(len=max_value_length) :: yaml_ztoa
    character(len=*), optional, intent(in) :: fmt
    !local variables
    integer :: ipos,rpos
    character(len=max_value_length) :: zr,zi
    double complex :: ztmp
    double precision, dimension(2) :: zeta

    yaml_ztoa=repeat(' ',max_value_length)
    ztmp=z
    zeta=transfer(ztmp,zeta)
    zr=yaml_ztoa
    zi=yaml_ztoa

    if (present(fmt)) then
       write(zr,fmt) zeta(1)
       write(zi,fmt) zeta(2)
    else
       write(zr,yaml_dble_fmt) zeta(1)
       write(zi,yaml_dble_fmt) zeta(2)
    end if
    
    zr=yaml_adjust(zr)
    zi=yaml_adjust(zi)
    rpos=len(trim(zr))
    ipos=min(len(trim(zi)),max_value_length-rpos-2)
    
    yaml_ztoa(1:rpos)=zr(1:rpos)
    if (zeta(2) >= 0.d0) then
       yaml_ztoa(rpos+1:rpos+2)='+'
       rpos=rpos+1
    end if
    yaml_ztoa(rpos+1:rpos+ipos)=zi(1:ipos)
    yaml_ztoa(rpos+ipos+1:rpos+ipos+2)='j'
  end function yaml_ztoa

  !> Convert logical to character
  pure function yaml_ltoa(l,fmt)
    implicit none
    logical, intent(in) :: l
    character(len=max_value_length) :: yaml_ltoa
    character(len=*), optional, intent(in) :: fmt

    yaml_ltoa=repeat(' ',max_value_length)

    if (present(fmt)) then
       write(yaml_ltoa,fmt) l
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
  pure function yaml_dvtoa(vec,fmt) result(vec_toa)
    implicit none
    real(kind=8), dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_dvtoa

  !> Convert vector of double complex to character
  pure function yaml_zvtoa(vec,fmt) result(vec_toa)
    implicit none
    double complex, dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_zvtoa

  !> Convert vector of integer to character
  pure function yaml_ivtoa(vec,fmt) result(vec_toa)
    implicit none
    integer, dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_ivtoa

  !> Convert vector of characters to a chain of characters
  pure function yaml_cvtoa(vec,fmt) result(vec_toa)
    implicit none
    character(len=*), dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_cvtoa

  pure function yaml_lvtoa(vec,fmt) result(vec_toa)
    implicit none
    logical, dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_lvtoa

  pure function yaml_rvtoa(vec,fmt) result(vec_toa)
    implicit none
    real, dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_rvtoa


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

  pure function yaml_adjust(str)
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

  !>find if a string is an integer
  !! use the portable mode described in 
  !! http://flibs.sourceforge.net/fortran_aspects.html#check_integers
  function is_atoi(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    integer :: ierr,ival
    character(len=20) :: form
    
    !fill the string describing the format to be used for reading
    !use the trimmed string and the yaml_toa function as i0 can add extra zeros in the specifications
    write(form,'(a20)')'(i'//adjustl(trim(yaml_itoa(len_trim(str),fmt='(i17)')))//')' 
    read(str,trim(form),iostat=ierr)ival
    yes=ierr==0
  end function is_atoi

  !>check if str contains a floating point number. 
  !!note that in principle this function gives positive answer also 
  !!if the number in str is an integer. Therefore is_atoi should be used to check before
  function is_atof(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    integer :: ierr,is,ie
    double precision :: rval

    ie=len(trim(str))
    is=max(scan(trim(str),' '),1)
    yes=scan(str(is:ie),' ') ==0 !there is no other space in the string
    if (yes) then
       read(str(is:ie),*,iostat=ierr)rval
       yes=ierr==0 .and. str(ie:ie)/='/' !the slash terminator is not allowed
    end if
  end function is_atof

  !> check if str contains a logical in yaml specification (Yes=.true. and No=.false.)
  function is_atol(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    integer :: is,ie
    !fill the string describing the format to be used for reading
    !use the trimmed string and the yaml_toa function as i0 can add extra zeros in the specifications
    ie=len(trim(str))
    is=max(scan(trim(str),' '),1)
    yes=scan(str(is:ie),' ') ==0 !there is no other space in the string
    if (yes) yes= (ie-is+1==3 .and. str(is:ie)=='Yes') .or. (ie-is+1==2 .and. str(is:ie)=='No')
  end function is_atol



  !> Shifts characters in in the string 'str' n positions (positive values
  !! denote a right shift and negative values denote a left shift). Characters
  !! that are shifted off the end are lost. Positions opened up by the shift 
  !! are replaced by spaces.
  !! This routine has been downloaded from the website http://gbenthien.net/strings/index.html
  pure subroutine shiftstr(str,n)
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

  end subroutine shiftstr

end module yaml_strings

