!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_libpaw_tools
!! NAME
!!  m_libpaw_tools
!!
!! FUNCTION
!!  Several libPAW tools: message printing, error handling, string handling...
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2014 ABINIT group (MT, MG, ...)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  This module comes directly from hide_write & hide_leave src files delivered with ABINIT.
!!
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/42_??libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_libpaw_tools
    
 USE_DEFS
 USE_MPI_WRAPPERS

#if defined HAVE_LIBPAW_BIGDFT
  use yaml_output
#endif

 implicit none

 private

!PUBLIC FUNCTIONS - MESSAGE HANDLING
 public  :: libpaw_wrtout        ! Parallel output of messages
 public  :: libpaw_msg_hndl      ! Basic error handler
 public  :: libpaw_flush         ! Wrapper for the standard flush routine

!PUBLIC FUNCTIONS - STRING HANDLING
 public  :: libpaw_basename      ! String, base name extraction from path
 public  :: libpaw_to_upper      ! String conversion to uppercase
 public  :: libpaw_lstrip        ! String right blanks removal
 public  :: libpaw_indent        ! String indentation

!PRIVATE FUNCTIONS
 private :: libpaw_wrtout_myproc ! Sequential output of messages
 private :: libpaw_write_lines   ! OS-compatible string output
 private :: libpaw_leave         ! Clean exit of F90 routines
 private :: libpaw_die           ! Clean exit
!!***

CONTAINS !===========================================================

!!****f* m_libpaw_tools/libpaw_wrtout
!! NAME
!! libpaw_wrtout
!!
!! FUNCTION
!!  Organizes the sequential or parallel version of the write intrinsic
!!
!! INPUTS
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!   "INIT" to change the rank of the master node that prints the message if "COLL" is used.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the WRTOUT routine delivered with ABINIT.
!!
!! SOURCE

subroutine libpaw_wrtout(unit,msg,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_wrtout'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: mode_paral

!Local variables ------------------------------
 integer :: comm,me,nproc
 integer,save :: master=0
 character(len=len(msg)+50) :: string
 character(len=500) :: my_mode_paral

!******************************************************************

 if (unit == dev_null) RETURN

 my_mode_paral = "COLL"; if (PRESENT(mode_paral)) my_mode_paral = mode_paral

!Communicator is xpaw_mpi_world by default
 comm=xpaw_mpi_world

!Determine who I am in COMM
 nproc = xpaw_mpi_comm_size(comm)
 me    = xpaw_mpi_comm_rank(comm)

 if ((my_mode_paral=='COLL').or.(nproc==1)) then
   if (me==master) then
     call libpaw_wrtout_myproc(unit,msg)
   end if
 else if (my_mode_paral=='PERS') then
   call libpaw_write_lines(unit,msg)
 else if (my_mode_paral=='INIT') then
   master=unit
 else
   write(string,'(7a)')ch10,&
&   'libpaw_wrtout: ERROR -',ch10,&
&   '  Unknown write mode: ',my_mode_paral,ch10,&
&   '  Continuing anyway ...'
   write(unit,'(A)') trim(string)
 end if

end subroutine libpaw_wrtout
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_wrtout_myproc
!! NAME
!!  libpaw_wrtout_myproc
!!
!! FUNCTION
!!  Do the output for one proc.
!!
!! INPUTS
!!  unit=unit number for writing
!!  msg=(character(len=*)) message to be written
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the WRTOUT_MYPROC routine delivered with ABINIT.
!!
!! SOURCE

subroutine libpaw_wrtout_myproc(unit,msg,mpicomm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_wrtout_myproc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 integer,intent(in),optional :: mpicomm
 character(len=*),intent(in) :: msg

!Local variables ------------------------------
!scalars
 integer,save :: iexit=0,ncomment=0,nwarning=0
 integer :: ierr
 logical :: print_std_err
!arrays
 integer :: buf(3)

!******************************************************************

!When I/O are redirected, it is sometimes necessary to reduce counters (saved) values;
!This can be done by passing mpicomm optional argument to the routine.
!In that case, no printing is done.
 if (present(mpicomm)) then
   buf(1)=iexit;buf(2)=ncomment;buf(3)=nwarning
   call xpaw_mpi_sum(buf,mpicomm,ierr)
   iexit=buf(1);ncomment=buf(2);nwarning=buf(3)
   if (iexit/=0) iexit=1
   return
 end if

 print_std_err=(unit==std_out.and.(index(trim(msg),'BUG')/=0.or.index(trim(msg),'ERROR')/=0))

 call libpaw_write_lines(unit,msg)
 if (print_std_err) then
   call libpaw_write_lines(std_err,msg)
 end if

 if (index(trim(msg),'BUG')/=0) then
   write(unit,'(a)') '  Action : contact libPAW developers.'
   if (print_std_err) write(std_err, '(a)' ) '  Action : contact libPAW developers.'
   write(unit,*); if (print_std_err) write(std_err,*)
   if(iexit/=0) write(unit, '(a)' ) ' Note: exit requested by the user.'
 end if

 if (index(trim(msg),'Exit')/=0) iexit=1

!Count the number of warnings and comments. Only take into
!account unit std_out, in order not to duplicate these numbers.
 if (index(trim(msg),'WARNING') /= 0 .and. unit==std_out) nwarning=nwarning+1
 if (index(trim(msg),'COMMENT') /= 0 .and. unit==std_out) ncomment=ncomment+1

end subroutine libpaw_wrtout_myproc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_write_lines
!! NAME
!!  libpaw_write_lines
!!
!! FUNCTION
!!  This routine receives a string, split the message in lines according to the 
!!  ch10 character and output the text to the specified unit. 
!!  Allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!! INPUTS
!!  unit=unit number for writing
!!  msg=(character(len=*)) message to be written
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the WRITE_LINES routine delivered with ABINIT.
!!
!! SOURCE

subroutine libpaw_write_lines(unit,msg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_write_lines'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg

!Local variables ------------------------------
!scalars
 integer :: msg_size,ii,jj,rtnpos

!******************************************************************

 msg_size=len_trim(msg)

#if defined HAVE_LIBPAW_BIGDFT
 if (msg_size>0 .and. unit==std_out) then
    call yaml_comment(msg)
 end if
 return
#endif

 if (msg_size==0) then
   write(unit,*) ; return 
 end if

!Here, split the message, according to the char(10) characters (carriage return). 
!This technique is portable accross different OS.
 rtnpos=index(msg,ch10)
 if (rtnpos==0) then
   write(unit,"(a)")msg(1:msg_size) ; return
 end if 

 ii=1; jj=rtnpos
 do 
   if (ii==jj) then
     write(unit,*)
   else
     write(unit,'(a)') msg(ii:jj-1)
   end if
   ii=jj+1 ; if (ii>msg_size) exit
   jj=index(msg(ii:msg_size),ch10) 
   if (jj==0) then 
     jj=msg_size+1
   else
     jj=jj+ii-1
   end if
 end do

 if (msg(msg_size:msg_size)==ch10) write(unit,*)

end subroutine libpaw_write_lines
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_msg_hndl
!! NAME
!!  libpaw_msg_hndl
!!
!! FUNCTION
!!  Basic error handler.
!!
!! INPUTS 
!!  msg=string containing additional information on the nature of the problem
!!  level=string defining the type of problem. Possible values are:
!!   COMMENT, WARNING, ERROR,BUG
!!  mode_paral=Either "COLL" or "PERS".
!!  [line]=line number of the file where problem occurred (optional)
!!  [file]=name of the f90 file containing the caller (optional)
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the MSG_HNDL routine delivered with ABINIT.
!!
!! SOURCE

subroutine libpaw_msg_hndl(msg,level,mode_paral,file,line)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_msg_hndl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: level,msg,mode_paral
 character(len=*),optional,intent(in) :: file

!Local variables ------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=LEN(msg)) :: my_msg 
 character(len=MAX(4*LEN(msg),2000)) :: sbuf

! *********************************************************************

 if (PRESENT(line)) f90line=line
 if (PRESENT(file)) f90name=libpaw_basename(file)

 my_msg=libpaw_lstrip(msg)
 write(sbuf,'(12a,i0,3a)')ch10,&
  "--- !",TRIM(libpaw_to_upper(level)),ch10,&
  "message: |",ch10,TRIM(libpaw_indent(my_msg)),ch10,&
  "src_file: ",TRIM(f90name),ch10,&
  "src_line: ",f90line,ch10,&
  "...",ch10

 select case (libpaw_to_upper(level))
 case ('COMMENT','WARNING')
   call libpaw_wrtout(std_out,sbuf,mode_paral) 
 case ('ERROR','BUG')
   call libpaw_wrtout(std_out,sbuf,mode_paral) 
   call libpaw_leave(mode_paral)
 case default 
   write(sbuf,'(4a)') ch10,' lbpaw_msg_hndl: BUG**2 - ',ch10,' Wrong value for level'
   call libpaw_die(sbuf,&
&  __FILE__,__LINE__)
 end select

end subroutine libpaw_msg_hndl
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_leave
!! NAME
!!  libpaw_leave
!!
!! FUNCTION
!!  Routine for clean exit of f90 code, taking into account possible parallelization.
!!
!! INPUTS
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same msg to be written once only
!!   'PERS' if the procs are calling the routine with different msgs each to be written,
!!          or if one proc is calling the routine
!!  [exit_status]=(optional, default=1 or -1, see below) the return code of the routine
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the LEAVE_NEW routine delivered with ABINIT.
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! NOTES
!!  This routine comes directly from the LEAVE_NEW routine delivered with ABINIT.
!!
!! SOURCE

subroutine libpaw_leave(mode_paral,exit_status)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_leave'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in),optional :: exit_status
 character(len=4),intent(in) :: mode_paral

!Local variables ------------------------------

! **********************************************************************

 call libpaw_wrtout(std_out,ch10//' leave_new : decision taken to exit ...','PERS')

!Caveat: Do not use MPI collective calls!
 if (mode_paral=="COLL") then
   call libpaw_wrtout(std_out,"Why COLL? Are you sure that ALL the processors are calling leave_new?")
 end if

 if (present(exit_status)) then
   call xpaw_mpi_abort(exit_status=exit_status)
 else
   call xpaw_mpi_abort()
 end if

end subroutine libpaw_leave
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_die
!! NAME
!!  libpaw_die
!!
!! FUNCTION
!!  Stop smoothly the execution in case of unexpected events reporting the
!!  line number and the file name where the error occurred as well as the 
!!  MPI rank of the processor.
!!
!! INPUTS 
!!  msg=String containing additional information on the nature of the problem
!!  [file]=Name of the f90 file containing the caller
!!  [line]=Line number of the file where problem occurred
!!
!! NOTES
!!  This routine comes directly from the DIE routine delivered with ABINIT.
!!
!! SOURCE

subroutine libpaw_die(message,file,line)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_die'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables ------------------------------
 integer :: rank 
 integer :: f90line=0
 character(len=10) :: lnum,strank
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=500) :: msg

! *********************************************************************

 if (PRESENT(line)) f90line=line
 if (PRESENT(file)) f90name= libpaw_basename(file)

 rank=xpaw_mpi_comm_rank(xpaw_mpi_world) !Determine my rank inside MPI_COMM_WORLD

 write(lnum,"(i0)") f90line
 write(strank,"(i0)") rank
 msg=TRIM(f90name)//':'//TRIM(lnum)//' P'//TRIM(strank)
 write(msg,'(a,2x,2a,2x,a)') ch10,TRIM(msg),ch10,TRIM(message)

 call libpaw_wrtout(std_out,msg,'PERS') 
 call libpaw_leave('PERS')

end subroutine libpaw_die
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_flush
!! NAME
!!  libpaw_flush
!!
!! FUNCTION
!!  Wrapper for the standard flush routine
!!  Available only if the compiler implements this intrinsic procedure.
!!
!! INPUTS
!!  unit=Fortran logical Unit number
!!
!! NOTES
!!  This routine comes directly from the FLUSH_UNIT routine delivered with ABINIT.
!!
!! SOURCE

subroutine libpaw_flush(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_flush'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit

!Local variables ------------------------------
 integer, parameter :: dev_null=-1
 logical :: isopen

!************************************************************************

 if (unit==dev_null) return

!FLUSH on unconnected unit is illegal: F95 std., 9.3.5.
 inquire(unit=unit,opened=isopen)

#if defined HAVE_FC_FLUSH
 if (isopen) then
   call flush(unit)
 endif
#elif defined HAVE_FC_FLUSH_
 if (isopen) then
   call flush_(unit)
  end if
#endif

end subroutine libpaw_flush
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_basename
!! NAME
!! libpaw_basename
!!
!! FUNCTION
!!  Returns the final component of a pathname (function version).
!!
!! INPUTS
!!  string=The input string
!!
!! NOTES
!!  This routine comes directly from the BASENAME routine delivered with ABINIT.
!!  If the input string in not a valid path to a file, a blank strink is returned
!!
!! SOURCE

pure function libpaw_basename(istr) result(ostr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_basename'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=LEN_TRIM(istr)) :: ostr

!Local variables ------------------------------
 integer :: ic,nch_trim,nch
 character(len=1),parameter :: BLANK=' '
 character(len=1),parameter :: DIR_SEPARATOR = '/'

!************************************************************************

 nch     =LEN     (istr)
 nch_trim=LEN_TRIM(istr)

 ic = INDEX (TRIM(istr), DIR_SEPARATOR, back=.TRUE.)
 if (ic >= 1 .and. ic <= nch_trim-1) then ! there is stuff after the separator.
   ostr = istr(ic+1:nch_trim)
 else if (ic==0 .or. ic == nch_trim+1) then ! no separator in string or zero length string,
   ostr = TRIM(istr)     ! return trimmed string.
 else                    ! (ic == nch_trim) separator is the last char.
   ostr = BLANK ! This is not a valid path to a file, return blank.
 end if
 return

end function libpaw_basename
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/to_upper
!! NAME
!!  libpaw_to_upper
!!
!! FUNCTION
!!  Convert a string to UPPER CASE (function version).
!!
!! INPUTS
!!   istr=Input string
!!
!! NOTES
!!  This routine comes directly from the TOUPPER routine delivered with ABINIT.
!!
!! SOURCE

pure function libpaw_to_upper(istr) result(ostr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_to_upper'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=LEN_TRIM(istr)) :: ostr

!Local variables ------------------------------
 integer,parameter :: ASCII_aa=ICHAR('a')
 integer,parameter :: ASCII_zz=ICHAR('z')
 integer,parameter :: SHIFT=ICHAR('a')-ICHAR('A')
 integer :: ic,iasc

! *********************************************************************

 do ic=1,LEN_TRIM(istr)
   iasc=IACHAR(istr(ic:ic))
   if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) then
     ostr(ic:ic)=ACHAR(iasc-SHIFT)
   else
     ostr(ic:ic)=istr(ic:ic)
   end if
 end do

end function libpaw_to_upper
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_lstrip
!! NAME
!!  libpaw_lstrip
!!
!! FUNCTION
!!  Removes leading spaces from the input string.
!!
!! NOTES
!!  This routine comes directly from the LSTRIP routine delivered with ABINIT.
!!
!! SOURCE

pure function libpaw_lstrip(istr) result(ostr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_lstrip'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=len(istr)) :: ostr

!Local variables ------------------------------
 integer :: ii,jj,lg
 character(len=1),parameter :: BLANK=' '

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

end function libpaw_lstrip
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_indent
!! NAME
!!  libpaw_indent
!!
!! FUNCTION
!!  Indent text (function version).
!!
!! INPUTS
!!   istr=Input string
!!
!! NOTES
!!  This routine comes directly from the INDENT routine delivered with ABINIT.
!!
!! SOURCE

pure function libpaw_indent(istr) result(ostr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_indent'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=len(istr)*4+4) :: ostr

!Local variables-------------------------------
 character(len=1),parameter :: NCHAR = char(10)
 integer,parameter :: n=4
 integer :: ii,jj,kk
 character(len=1) :: ch

! *********************************************************************

 ostr=" "
 jj=n
 do ii=1,LEN_TRIM(istr)
   ch=istr(ii:ii)
   jj=jj+1
   if (ch==NCHAR) then
      ostr(jj:jj)=NCHAR
      do kk=jj+1,jj+n
        ostr(kk:kk)=" "
      end do
      jj=jj+n
   else
     ostr(jj:jj)=ch
   end if
 end do

end function libpaw_indent
!!***

!----------------------------------------------------------------------

end module m_libpaw_tools
!!***
