!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fstrings
!! NAME
!!  m_fstrings
!!
!! FUNCTION
!!  This module contains basic tools to operate on Fortran strings.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2013 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

!#include "abi_common.h"

MODULE m_fstrings

 use defs_basis, only : dp, std_out

 implicit none

 private 
 
 public :: is_letter       ! Returns .TRUE. if ch is a letter and .FALSE. otherwise
 public :: is_digit        ! Returns .TRUE. if ch is a digit (0,1,...,9) and .FALSE. otherwise
 public :: upper           ! Convert lower case letters to UPPER CASE
 public :: toupper         ! Convert lower case letters to UPPER CASE (function version)
 public :: lower           ! Convert UPPER CASE letters to lower case 
 public :: tolower         ! Convert UPPER CASE letters to lower case  (function version)
 public :: compact         ! Converts multiple spaces and tabs to single spaces; deletes control characters and initial spaces
 public :: removesp        ! Removes spaces, tabs, and control characters in string str
 public :: tovalue         ! Converts number string to a number
 public :: lstrip          ! Remove leading spaces from string
 public :: ljust           ! Return a left-justified string of length width.
 public :: write_num       ! Writes a number to a string using format fmt
 public :: trimzero        ! Deletes nonsignificant trailing zeroes from a number string. 
 public :: writeq          ! Writes a string of the form <name> = value to unit
 public :: match           ! Find the position of a delimiter in a string
 public :: OPERATOR(.jn.)  
 public :: strcat          ! Concatenate two strings (function version for buggy PGI)
 public :: itoa            ! Convert an integer into a string
 public :: basename        ! Returns the final component of a pathname.
 public :: starts_with     ! Returns .TRUE. is the first character in a string belongs to a gives set. 

 !TODO method to center a string

 interface tovalue
  module procedure value_int
  module procedure value_dp
 end interface 

 interface write_num
  module procedure write_rdp_0D
  module procedure write_int_0D
 end interface

 interface writeq
  module procedure writeq_rdp_0D
  module procedure writeq_int_0D
 end interface

 interface is_digit
  module procedure is_digit_0D
 end interface

 interface operator (.jn.)    
  module procedure strcat
 end interface 

  character(len=1),parameter :: BLANK=' ' 
  character(len=1),parameter :: DIR_SEPARATOR = '/'

  integer,parameter :: ASCII_A=ICHAR('A')
  integer,parameter :: ASCII_Z=ICHAR('Z')
  integer,parameter :: ASCII_aa=ICHAR('a')
  integer,parameter :: ASCII_zz=ICHAR('z')
  integer,parameter :: SHIFT=ASCII_aa-ASCII_A ! Capital letters have smaller Dec value in the ASCII table.
  integer,parameter :: ASCII_0=ICHAR('0')
  integer,parameter :: ASCII_9=ICHAR('9')

CONTAINS  !===========================================================
!!*** 

!!****f* m_fstrings/is_letter
!! NAME
!!  is_letter 
!!
!! FUNCTION
!!  Returns .TRUE. if ch is a letter and .FALSE. otherwise.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function is_letter(ch) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'is_letter'
!End of the abilint section

 character(len=1),intent(in) :: ch
 logical :: ans
! *********************************************************************

 select case (ICHAR(ch))
 case (ASCII_A:ASCII_Z,ASCII_aa:ASCII_zz)
  ans=.TRUE.
 case DEFAULT
  ans=.FALSE.
 end select

end function is_letter
!!***

!!****f* m_fstrings/is_digit_0D
!! NAME
!!  is_digit_0D
!!
!! FUNCTION
!!  Returns .TRUE. if ch is a digit (0,1,...,9) and .FALSE. otherwise.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
function is_digit_0D(ch) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'is_digit_0D'
!End of the abilint section

 character(len=1),intent(in) :: ch
 logical :: ans
! *********************************************************************

 select case (ICHAR(ch))
 case(ASCII_0:ASCII_9)
  ans=.TRUE.
 case default
  ans=.FALSE.
 end select

end function is_digit_0D
!!***

!!****f* m_fstrings/upper
!! NAME
!!  upper 
!!
!! FUNCTION
!!  Convert lower case letters to UPPER CASE.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine upper(str)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'upper'
!End of the abilint section

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: ic,iasc
! ********************************************************************* 

 do ic=1,LEN_TRIM(str)
  iasc=IACHAR(str(ic:ic))
  if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) str(ic:ic)=ACHAR(iasc-SHIFT)
 end do

end subroutine upper
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/toupper
!! NAME
!!  toupper 
!!
!! FUNCTION
!!  Convert lower case letters to UPPER CASE (function version).
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function toupper(str_in) result(str_out)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'toupper'
!End of the abilint section

 character(len=*),intent(in) :: str_in
 character(len=LEN_TRIM(str_in)) :: str_out

!Local variables-------------------------------
 integer :: ic,iasc
! ********************************************************************* 

 do ic=1,LEN_TRIM(str_in)
  iasc=IACHAR(str_in(ic:ic))
  if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) then 
   str_out(ic:ic)=ACHAR(iasc-SHIFT)
  else 
   str_out(ic:ic)=str_in(ic:ic)
  end if
 end do

end function toupper
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/lower
!! NAME
!!  lower 
!!
!! FUNCTION
!!  Convert UPPER CASE letters to lower case.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      fftprof,ioprof,lapackprof
!!
!! CHILDREN
!!
!! SOURCE

subroutine lower(str)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lower'
!End of the abilint section

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str)
  iasc=IACHAR(str(ic:ic))
  if (iasc>=ASCII_A.and.iasc<=ASCII_Z) str(ic:ic)=ACHAR(iasc+SHIFT)
 end do

end subroutine lower
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/tolower
!! NAME
!!  tolower 
!!
!! FUNCTION
!!  Convert UPPER CASE letters to lower case (function version).
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
function tolower(str_in) result(str_out)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tolower'
!End of the abilint section

 character(len=*),intent(in) :: str_in
 character(len=LEN_TRIM(str_in)) :: str_out

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str_in)
  iasc=IACHAR(str_in(ic:ic))
  if (iasc>=ASCII_A.and.iasc<=ASCII_Z) then 
   str_out(ic:ic)=ACHAR(iasc+SHIFT)
  else 
   str_out(ic:ic)=str_in(ic:ic)
  end if
 end do

end function tolower
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/compact
!! NAME
!!  compact  
!!
!! FUNCTION
!! Converts multiple spaces and tabs to single spaces; 
!! deletes control characters; removes initial spaces.
!!  
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine compact(str)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compact'
!End of the abilint section

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: isp,i,k,lenstr,ich
 character(len=1):: ch
 character(len=LEN_TRIM(str)):: outstr
! *********************************************************************

 str=ADJUSTL(str) ; lenstr=LEN_TRIM(str)

 outstr=BLANK ; isp=0 ; k=0
 do i=1,lenstr
  ch=str(i:i) ; ich=IACHAR(ch)
   
  select case(ich)
  case(9,32)     ! space or tab character
   if (isp==0) then
    k=k+1
    outstr(k:k)=' '
   end if
   isp=1
  case(33:)      ! not a space, quote, or control character
   k=k+1
   outstr(k:k)=ch
   isp=0
  end select
 end do

 str=ADJUSTL(outstr)

end subroutine compact
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/removesp
!! NAME
!!  removesp 
!!
!! FUNCTION
!!  Removes spaces, tabs, and control characters in string str.
!!  
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine removesp(str)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'removesp'
!End of the abilint section

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: i,k,lenstr,ich
 character(len=1):: ch
 character(len=LEN_TRIM(str)):: outstr
! *********************************************************************

 str=ADJUSTL(str) ; lenstr=LEN_TRIM(str)

 outstr=BLANK ; k=0
 do i=1,lenstr
  ch=str(i:i)
  ich=IACHAR(ch)
  select case(ich)    
  case(0:32)  ! space, tab, or control character
   CYCLE       
  case(33:)  
   k=k+1
   outstr(k:k)=ch
  end select
 end do

 str=ADJUSTL(outstr)

end subroutine removesp
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/value_dp
!! NAME
!!  value_dp
!!
!! FUNCTION
!!  Convert number string to a number
!!  
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      m_fstrings
!!
!! CHILDREN
!!
!! SOURCE
subroutine value_dp(str,rnum,ios)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'value_dp'
!End of the abilint section

 character(len=*),intent(in) ::str
 integer,intent(out) :: ios
 real(dp),intent(out) :: rnum

!Local variables-------------------------------
 integer :: ilen,ipos
! *********************************************************************

 ilen=LEN_TRIM(str) ; ipos=SCAN(str,'Ee')
 if (.not.is_digit(str(ilen:ilen)).and.ipos/=0) then
  ios=3
  RETURN
 end if
 read(str,*,iostat=ios) rnum

end subroutine value_dp
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/value_int
!! NAME
!!  value_int
!!
!! FUNCTION
!!  Convert number string to a number.
!!  
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine value_int(str,inum,ios)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'value_int'
!End of the abilint section

 character(len=*),intent(in) ::str
 integer,intent(out) :: inum,ios

!Local variables-------------------------------
 real(dp) :: rnum
! *********************************************************************

 call value_dp(str,rnum,ios)
 if (ABS(rnum)>HUGE(inum)) then
  ios=1
  RETURN
 end if
 inum=NINT(rnum)

end subroutine value_int 
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/lstrip
!! NAME
!!  lstrip 
!!
!! FUNCTION
!!  Removes leading spaces from the input string.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function lstrip(istr) result(ostr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lstrip'
!End of the abilint section

 character(len=*),intent(in) :: istr
 character(len=len(istr)) :: ostr

!Local variables-------------------------------
 integer :: ii,jj,lg
! *********************************************************************

 lg=LEN(istr)
 do ii=1,lg
   if (istr(ii:ii)/=BLANK) EXIT
 end do

 do jj=1,lg-ii+1
   ostr(jj:jj) = istr(ii:ii)
   ii=ii+1
 end do
 ostr(ii:lg)=BLANK

end function lstrip
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/ljust
!! NAME
!!  ljust 
!!
!! FUNCTION
!!  Return S left-justified in a string of length width. Padding is
!!  done using the specified fill character (default is a space).
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ljust(istr, width, fillchar) result(ostr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ljust'
!End of the abilint section

 character(len=*),intent(in) :: istr
 integer,intent(in) :: width 
 character(len=width) :: ostr
 character(len=1),optional,intent(in) :: fillchar

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 ostr = ADJUSTL(istr)

 if (PRESENT(fillchar)) then
   do ii=LEN_TRIM(ostr)+1,width
     ostr(ii:ii) = fillchar
   end do
 end if

end function ljust
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/write_rdp_0d
!! NAME
!!  write_rdp_0d
!!
!! FUNCTION
!!  Writes a number to a string using format fmt.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_rdp_0d(rnum,str,fmt)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_rdp_0d'
!End of the abilint section

 real(dp),intent(in) :: rnum
 character(len=*),intent(in) :: fmt
 character(len=*),intent(out) :: str

!Local variables-------------------------------
 character(len=LEN(fmt)+2) :: formt
! *********************************************************************

 formt='('//TRIM(fmt)//')'
 write(str,formt)rnum
 str=ADJUSTL(str)

end subroutine write_rdp_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/write_int_0d
!! NAME
!!  write_int_0d
!!
!! FUNCTION
!!  Writes a number to a string using format fmt.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine write_int_0D(inum,str,fmt)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_int_0D'
!End of the abilint section

 integer,intent(in) :: inum
 character(len=*),intent(in) :: fmt
 character(len=*),intent(out) :: str

!Local variables-------------------------------
 character(len=LEN(fmt)+2) :: formt
! *********************************************************************

 formt='('//TRIM(fmt)//')'
 write(str,formt) inum
 str=ADJUSTL(str)

end subroutine write_int_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/trimzero
!! NAME
!!  trimzero 
!!
!! FUNCTION
!! Deletes nonsignificant trailing zeroes from number string str. If number
!! string ends in a decimal point, one trailing zero is added.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      m_fstrings
!!
!! CHILDREN
!!
!! SOURCE
! NOT sure it will work

subroutine trimzero(str)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'trimzero'
!End of the abilint section

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: i,ipos,lstr
 character :: ch
 character(len=10) :: sexp
! *********************************************************************

 ipos=SCAN(str,'eE')
 if (ipos>0) then
  sexp=str(ipos:)
  str=str(1:ipos-1)
 end if
 lstr=LEN_TRIM(str)
 do i=lstr,1,-1
  ch=str(i:i)
  if (ch=='0') CYCLE          
  if (ch=='.') then
   str=str(1:i)//'0'
   if (ipos>0) str=TRIM(str)//TRIM(sexp)
   EXIT
  end if
  str=str(1:i)
  EXIT
 end do

 if (ipos>0) str=TRIM(str)//TRIM(sexp)

end subroutine trimzero
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/writeq_rdp_0D 
!! NAME
!!  writeq_rdp_0D
!!
!! FUNCTION
!!  Writes a string of the form <name> = value to unit.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine writeq_rdp_0D(unit,namestr,value,fmt)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'writeq_rdp_0D'
!End of the abilint section

 real(dp),intent(in) :: value
 integer,intent(in) :: unit
 character(len=*),intent(in) :: fmt
 character(len=*),intent(in) :: namestr

!Local variables-------------------------------
 character(len=32) :: tempstr
! *********************************************************************

 call write_num(value,tempstr,fmt)
 call trimzero(tempstr)
 write(unit,*)TRIM(namestr)//' = '//TRIM(tempstr)

end subroutine writeq_rdp_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/writeq_int_0D 
!! NAME
!!  writeq_int_0D
!!
!! FUNCTION
!!  Writes a string of the form <name> = value to unit.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine writeq_int_0D(unit,namestr,ivalue,fmt)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'writeq_int_0D'
!End of the abilint section

 integer,intent(in) :: ivalue
 integer,intent(in) :: unit
 character(len=*),intent(in) :: namestr
 character(len=*),intent(in) :: fmt

!Local variables-------------------------------
 character(len=32) :: tempstr
! *********************************************************************

 call write_num(ivalue,tempstr,fmt)
 call trimzero(tempstr)
 write(unit,*)TRIM(namestr)//' = '//TRIM(tempstr)

end subroutine writeq_int_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/match
!! NAME
!!  match 
!!
!! FUNCTION
!!  Sets imatch to the position in string of the delimiter matching the delimiter
!!  in position ipos. Allowable delimiters are (), [], {}, <>.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine match(str,ipos,imatch,ios)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'match'
!End of the abilint section

 character(len=*),intent(in) :: str
 integer,intent(in) :: ipos
 integer,intent(out) :: ios,imatch

!Local variables-------------------------------
 integer :: i,isum,lenstr,idelim2,istart,iend,inc
 character :: delim1,delim2,ch
! *********************************************************************

 lenstr=LEN_TRIM(str)
 delim1=str(ipos:ipos)
 select case(delim1)
  case('(')
   idelim2=IACHAR(delim1)+1
   istart=ipos+1
   iend=lenstr
   inc=1
  case(')')
   idelim2=IACHAR(delim1)-1
   istart=ipos-1
   iend=1
   inc=-1
  case('[','{','<')
   idelim2=IACHAR(delim1)+2
   istart=ipos+1
   iend=lenstr
   inc=1
  case(']','}','>')
   idelim2=IACHAR(delim1)-2
   istart=ipos-1
   iend=1
   inc=-1
  case default
   write(std_out,*)delim1,' is not a valid delimiter'
   RETURN
 end select

 if (istart<1 .or. istart>lenstr) then
  write(std_out,*) delim1,' has no matching delimiter'
  RETURN
 end if

 delim2=ACHAR(idelim2) ! matching delimiter

 isum=1
 do i=istart,iend,inc
  ch=str(i:i)
  if (ch/=delim1 .and. ch/=delim2) CYCLE
  if (ch==delim1) isum=isum+1
  if (ch==delim2) isum=isum-1
  if (isum==0) EXIT
 end do

 if (isum/=0) then
  write(std_out,*)delim1,' has no matching delimiter'
  ios=1
  RETURN
 end if   
 ios=0 ; imatch=i

end subroutine match
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat
!! NAME
!! strcat
!!
!! FUNCTION 
!!  Returns two concatenated strings.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
function strcat(str1,str2) result(cnct)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strcat'
!End of the abilint section

 character(len=*),intent(in) :: str1,str2
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)) :: cnct

!Local variables-------------------------------
! *********************************************************************

 cnct=TRIM(str1)//TRIM(str2)

end function strcat 
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/itoa
!! NAME
!! itoa
!!
!! FUNCTION 
!!  Convert an integer into a string
!!
!! INPUTS
!!  value=The integer
!!
!! PARENTS
!!
!! CHILDREN
!!

function itoa(value) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'itoa'
!End of the abilint section

 integer,intent(in) :: value
 character(len=22) :: itoa

! *********************************************************************

 ! len=22 is large enough to contain integer*8
 write(itoa,"(i0)")value
 itoa = ADJUSTL(itoa)

end function itoa 
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/basename
!! NAME
!! basename
!!
!! FUNCTION
!!  Returns the final component of a pathname.
!!
!! INPUTS
!!  string=The input string
!!
!! NOTES
!!  * If the input string in not a valid path to a file (i.e not in the form foo/name) 
!!    a blank strink is returned
!!  * We do a backward search becase we want to optimize the algorithm for Fortran strings.
!!
!! PARENTS
!!
!! CHILDREN
!!  
!!
!! SOURCE

function basename(string)
    
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'basename'
!End of the abilint section

 character(len=*),intent(in) :: string
 character(len=LEN_TRIM(string)) :: basename

!Local variables-------------------------------
 integer :: ic,nch_trim,nch
!************************************************************************

 nch     =LEN     (string)
 nch_trim=LEN_TRIM(string)

 ic = INDEX (TRIM(string), DIR_SEPARATOR, back=.TRUE.)
 !write(std_out,*)'DEBUG ',TRIM(string),ic

 if (ic >= 1 .and. ic <= nch_trim-1) then ! there is stuff after the separator.
  basename = string(ic+1:nch_trim)
  return
 else if (ic==0 .or. ic == nch_trim+1) then ! no separator in string or zero length string, 
  basename = TRIM(string)                   ! return trimmed string.
  return
 else              ! (ic == nch_trim) separator is the last char.
  basename= BLANK  ! This is not a valid path to a file, return blank.
  return
 end if

end function basename
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/starts_with
!! NAME
!! starts_with
!!
!! FUNCTION
!!  Returns .TRUE. is the first character of the string belongs to a given list.
!!
!! INPUTS
!!  string=The string whose first character has to be cheched
!!  char_list=The list of characters.
!!  [csens]=.TRUE. if comparison is done regardless of case. Defaults to .FALSE.
!!
!! PARENTS
!!
!! CHILDREN
!!  
!!
!! SOURCE

function starts_with(string,char_list,csens) result(ans)
    
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'starts_with'
!End of the abilint section

 logical :: ans
 logical,optional,intent(in) :: csens
 character(len=*),intent(in) :: string
 character(len=1),intent(in) :: char_list(:)

!Local variables-------------------------------
 integer :: ii
 logical :: my_csens
 character(len=1) :: first_ch
!************************************************************************

 my_csens=.FALSE.; if (PRESENT(csens)) my_csens = csens

 first_ch = string(1:1)

 ans=.FALSE.

 if (.not.my_csens) then
   do ii=1,SIZE(char_list)
     ans = ( first_ch == char_list(ii) ); if (ans) EXIT
   end do
 else 
   do ii=1,SIZE(char_list)
     ans = ( toupper(first_ch) == toupper(char_list(ii)) ); if (ans) EXIT
   end do
 end if

end function starts_with
!!***

END MODULE m_fstrings
!!***
