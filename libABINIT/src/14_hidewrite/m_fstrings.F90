!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fstrings
!! NAME
!!  m_fstrings
!!
!! FUNCTION
!!  This module contains basic tools to operate on Fortran strings.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MG, XG, MT)
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
#include "config.h"
#endif

#include "abi_common.h"

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
 public :: lstrip          ! Remove leading spaces from string
 public :: ljust           ! Return a left-justified string of length width.
 public :: lpad            ! Pad a string adding repeat characters fillchar on the left side.
 public :: quote           ! Return a new string enclosed by quotation marks.
 public :: write_num       ! Writes a number to a string using format fmt
 public :: trimzero        ! Deletes nonsignificant trailing zeroes from a number string.
 public :: writeq          ! Writes a string of the form <name> = value to unit
 public :: strcat          ! Concatenate strings (function version)
 public :: sjoin           ! Joins strings with a space separator.
 public :: itoa            ! Convert an integer into a string
 public :: atoi            ! Convert a string into a integer
 public :: basename        ! Returns the final component of a pathname.
 public :: starts_with     ! Returns .TRUE. is the first character in a string belongs to a gives set.
 public :: indent          ! Indent text
 public :: int2char4       ! Convert a positive integer number (zero included) to a character(len=*)
                           ! with trailing zeros if the number is <=9999
 public :: int2char10      ! Convert a positive integer number (zero included) to a character(len=10)
                           ! with trailing blanks

 !TODO method to center a string

 interface write_num
   module procedure write_rdp_0D
   module procedure write_int_0D
 end interface write_num

 interface writeq
   module procedure writeq_rdp_0D
   module procedure writeq_int_0D
 end interface writeq

 interface is_digit
   module procedure is_digit_0D
 end interface is_digit

 interface starts_with
   module procedure starts_with_0d
   module procedure starts_with_1d
 end interface starts_with

 interface sjoin
   module procedure sjoin_2
 end interface sjoin

 interface strcat
   module procedure strcat_2
   module procedure strcat_3
   module procedure strcat_4
 end interface strcat

  character(len=1),parameter :: BLANK=' '
  character(len=1),parameter :: NCHAR = char(10)
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

 ostr = " "
 do jj=1,lg-ii+1
   ostr(jj:jj) = istr(ii:ii)
   ii=ii+1
 end do

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

!!****f* m_fstrings/lpad
!! NAME
!!  lpad
!!
!! FUNCTION
!!  Pad a string adding repeat characters fillchar on the left side.
!!  Padding is done using the specified fill character (default is a blanck character).
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

function lpad(istr, repeat, fillchar) result(ostr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lpad'
!End of the abilint section

 character(len=*),intent(in) :: istr
 integer,intent(in) :: repeat
 character(len=LEN_TRIM(istr) + repeat) :: ostr
 character(len=1),optional,intent(in) :: fillchar

!Local variables-------------------------------
 integer :: ii
 character(len=1) :: ch
! *********************************************************************

 ostr(repeat+1:) = TRIM(istr)

 ch = " "; if (PRESENT(fillchar)) ch = fillchar
 do ii=1,repeat
   ostr(ii:ii) = ch
 end do

end function lpad
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/quote
!! NAME
!!  quote
!!
!! FUNCTION
!!  Return a new string enclosed by quotation marks.
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

function quote(istr) result(ostr)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'quote'
!End of the abilint section

 character(len=*),intent(in) :: istr
 character(len=LEN_TRIM(istr)+2) :: ostr

!Local variables-------------------------------
 integer :: ii
 character(len=1) :: qq
 character(len=LEN(istr)+2) :: tmp

! *********************************************************************

 do ii=1,LEN(istr)
   if (istr(ii:ii)/=BLANK) EXIT
 end do

 qq = istr(ii:ii)

 if (qq == "'" .or. qq == '"') then
   ! Don't add quotation marks if they already present.
   tmp = istr
   ii = LEN_TRIM(tmp)
   ! If the string is not closed, fix it.
   if (tmp(ii:ii) /= qq) tmp(ii+1:ii+1) = qq
   ostr = TRIM(tmp)

 else
   qq = '"'
   ostr(1:1) = qq
   ostr(2:) = TRIM(istr)
   ii = LEN_TRIM(ostr)+1
   ostr(ii:ii) = qq
 end if

end function quote
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

!!****f* m_fstrings/sjoin_2
!! NAME
!! sjoin_2
!!
!! FUNCTION
!!  Joins two strings with a space separator unless first string is empty.
!!

function sjoin_2(str1,str2) result(ostr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sjoin_2'
!End of the abilint section

 character(len=*),intent(in) :: str1,str2
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+1) :: ostr

! *********************************************************************

 if (len_trim(str1) > 0) then
   ostr=TRIM(str1)//" "//TRIM(str2)
 else
   ostr=TRIM(str2)
 end if

end function sjoin_2
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat_2
!! NAME
!! strcat_2
!!
!! FUNCTION
!!  Returns two concatenated strings.
!!

function strcat_2(str1,str2) result(ostr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strcat_2'
!End of the abilint section

 character(len=*),intent(in) :: str1,str2
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)) :: ostr

! *********************************************************************

 ostr=TRIM(str1)//TRIM(str2)

end function strcat_2
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat_3
!! NAME
!! strcat_3
!!
!! FUNCTION
!!  Concatenate 3 strings
!!

function strcat_3(str1, str2, str3) result(ostr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strcat_3'
!End of the abilint section

 character(len=*),intent(in) :: str1,str2,str3
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)) :: ostr

! *********************************************************************

 ostr = TRIM(str1)//TRIM(str2)//TRIM(str3)

end function strcat_3
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat_4
!! NAME
!! strcat_3
!!
!! FUNCTION
!!  Concatenate 4 strings
!!

function strcat_4(str1, str2, str3, str4) result(ostr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strcat_4'
!End of the abilint section

 character(len=*),intent(in) :: str1,str2,str3,str4
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+LEN_TRIM(str4)) :: ostr

! *********************************************************************

 ostr = TRIM(str1)//TRIM(str2)//TRIM(str3)//TRIM(str4)

end function strcat_4
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/atoi
!! NAME
!! atoi
!!
!! FUNCTION
!!  Convert a string into a integer
!!
!! PARENTS
!!
!! CHILDREN
!!

function atoi(string)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'atoi'
!End of the abilint section

 integer :: atoi
 character(len=*),intent(in) :: string

! *********************************************************************

 read(string,*,err=10)atoi
 return
 10 write(std_out,*)"Error while trying to convert string to integer. string: ",trim(string)

end function atoi
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
 !write(*,*)'DEBUG ',TRIM(string),ic

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

!!****f* m_fstring/starts_with_0d
!! NAME
!! starts_with_0d
!!
!! FUNCTION
!!   Return True if string starts with the specified character
!!
!! INPUTS
!!  string=The string whose first character has to be cheched
!!  ch=Character
!!  [csens]=.TRUE. if comparison is done regardless of case. Defaults to .FALSE.
!!
!! PARENTS
!!
!! CHILDREN
!!
!!
!! SOURCE

function starts_with_0d(string,ch,csens) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'starts_with_0d'
!End of the abilint section

 logical :: ans
 logical,optional,intent(in) :: csens
 character(len=*),intent(in) :: string
 character(len=1),intent(in) :: ch

!Local variables-------------------------------
 logical :: my_csens
!************************************************************************

 my_csens=.FALSE.; if (PRESENT(csens)) my_csens = csens

 if (.not.my_csens) then
   ans = ( string(1:1) == ch)
 else
   ans = ( toupper(string(1:1)) == toupper(ch))
 end if

end function starts_with_0d
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/starts_with_1d
!! NAME
!! starts_with_1d
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

function starts_with_1d(string,char_list,csens) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'starts_with_1d'
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

end function starts_with_1d
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/indent
!! NAME
!!  indent
!!
!! FUNCTION
!!  Indent text
!!
!! INPUTS
!!   istr=Input string
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function indent(istr) result(ostr)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indent'
!End of the abilint section

 character(len=*),intent(in) :: istr
 character(len=len(istr)*4+4) :: ostr

!Local variables-------------------------------
 integer,parameter :: n=4 ! ostr is large enough to allocate all the possible indentations.
 integer :: ii,jj,kk
 character(len=1) :: ch

! *********************************************************************

 ostr = " "
 jj = n
 do ii=1,LEN_TRIM(istr)
   ch = istr(ii:ii)
   jj = jj + 1
   if (ch == NCHAR) then
      ostr(jj:jj) = NCHAR
      do kk=jj+1,jj+n
        ostr(kk:kk) = " "
      end do
      jj = jj+n
   else
     ostr(jj:jj) = ch
   end if
 end do
 !ostr(jj+1:) = "H"

end function indent
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/int2char4
!! NAME
!! int2char4
!!
!! FUNCTION
!! Convert an integer number to ("2") a character(len=*)
!! with trailing zeros if the number is <=9999.
!! Exemple : 123 will be mapped to "0123" ; 12345 will be mapped to "12345"
!! Makes sure that the integer fits the string length
!! (ex.: between 0 and 99999 if the string is a character(len=5)).
!!
!! INPUTS
!!  iint=integer to be converted
!!
!! OUTPUT
!!  string=character string ('####...' if error)
!!
!! PARENTS
!!      aim,anaddb,driver,dtfil_init1,gaus_dos,get_all_gkq,iofn1,m_atprj
!!      m_green,m_io_redirect,m_phonon_supercell,m_self,mrgscr,optic,pawmkaewf
!!      prtfatbands,read_wfrspa,scfcv,tddft,tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

subroutine int2char4(iint,string)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int2char4'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iint
 character(len=*),intent(out) :: string

!Local variables-------------------------------
 integer :: lenstr

! *************************************************************************

 lenstr=min(len(string),25)
 if(iint<0 .or. iint>10._dp**(lenstr-1))then
   string=repeat('#',lenstr)
   return
 end if
 if(iint<10)then
   write(string,'("000",i1)')iint
 else if(iint<100)then
   write(string,'("00",i2)')iint
 else if(iint<1000)then
   write(string,'("0",i3)')iint
 else if(iint<10000)then
   write(string,'(i4)')iint
 else if(iint<1.0d5)then
   write(string,'(i5)')iint
 else if(iint<1.0d6)then
   write(string,'(i6)')iint
 else if(iint<1.0d7)then
   write(string,'(i7)')iint
 else if(iint<1.0d8)then
   write(string,'(i8)')iint
 else if(iint<1.0d9)then
   write(string,'(i9)')iint
 else if(iint<1.0d9)then
   write(string,'(i10)')iint
 else
   string=repeat('#',lenstr)
 end if

end subroutine int2char4
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/int2char10
!! NAME
!! int2char10
!!
!! FUNCTION
!! Convert a positive integer number (zero included) to a character(len=10),
!! with blanks to COMPLETE the string.
!! Exemple : 1234 will be mapped to "1234      "
!! Makes sure that the integer is between 0 and 9 999 999 999
!! Should be enough for integer*4
!!
!! INPUTS
!!  iint=integer to be converted
!!
!! OUTPUT
!!  string=character string ('##########' if error)
!!
!! PARENTS
!!      handle_ncerr,m_bands_sym,m_dfti,m_dyson_solver,m_qparticles,m_wfs
!!      prt_cif,wffile
!!
!! CHILDREN
!!
!! SOURCE

subroutine int2char10(iint,string)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int2char10'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iint
 character(len=10),intent(out) :: string

! *************************************************************************

!Note the use of floating numbers instead of large integers, for portability
 if(iint<0 .or. iint>=1.d10)then
   string='####'
   return
 end if
 if(iint<10)then
   write(string,'(i1,9x)')iint
 else if(iint<100)then
   write(string,'(i2,8x)')iint
 else if(iint<1.0d3)then
   write(string,'(i3,7x)')iint
 else if(iint<1.0d4)then
   write(string,'(i4,6x)')iint
 else if(iint<1.0d5)then
   write(string,'(i5,5x)')iint
 else if(iint<1.0d6)then
   write(string,'(i6,4x)')iint
 else if(iint<1.0d7)then
   write(string,'(i7,3x)')iint
 else if(iint<1.0d8)then
   write(string,'(i8,2x)')iint
 else if(iint<1.0d9)then
   write(string,'(i9,1x)')iint
 else
   write(string,'(i10)')iint
 end if

end subroutine int2char10
!!***

!----------------------------------------------------------------------

!!!
!!! !!****f* m_fstrings/n2ch10
!!! !! NAME
!!! !!  n2ch10
!!! !!
!!! !! FUNCTION
!!! !!
!!! !!
!!! !! INPUTS
!!! !!   istr=Input string
!!! !!
!!! !! PARENTS
!!! !!
!!! !! CHILDREN
!!! !!
!!! !! SOURCE
!!!
!!! function n2ch10(istr) result(ostr)
!!!
!!! !Arguments ------------------------------------
!!! !scalars
!!!
!!! !This section has been created automatically by the script Abilint (TD).
!!! !Do not modify the following lines by hand.
!!! #undef ABI_FUNC
!!! #define ABI_FUNC 'n2ch10'
!!! !End of the abilint section
!!!
!!!  character(len=*),intent(in) :: istr
!!!  character(len=len(istr)*2) :: ostr
!!!
!!! !Local variables-------------------------------
!!!  integer :: ii,jj,lentrim
!!!  integer :: in_escape
!!!  character(len=1) :: ch,next_ch
!!!
!!! ! *********************************************************************
!!!
!!!  ostr = ""
!!!
!!!  lentrim = LEN_TRIM(istr)
!!!
!!!  in_escape = 0
!!!  if (istr(1:1) == "\") in_escape = 1
!!!  if (in_escape==0) ostr(1:1) = istr(1:1)
!!!
!!!  jj = 0
!!!  do ii=2,lentrim
!!!    ch = istr(ii:ii)
!!!    jj = jj + 1
!!!    if (in_escape==0) then
!!!      ostr(jj:jj) = ch
!!!    else if (in_escape==1) then
!!!      next_ch =  istr(ii+1:ii+1)
!!!      if next_ch == "n"
!!!        ostr(jj:jj) = ch10
!!!        in_escape = in_escape + 1
!!!      else
!!!    end if
!!!  end do
!!!
!!! end function n2ch10
!!! !!***
!!!
!!! !----------------------------------------------------------------------

END MODULE m_fstrings
!!***
