!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_wrtout
!! NAME
!!  abi_wrtout
!!
!! FUNCTION
!!  Organizes the sequential or parallel version of the write intrinsic
!!  Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
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
!!   The routine uses optional arguments, therefore the interface must be explicit.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_wrtout(unit,msg,mode_paral)

 use abi_defs_basis
 use abi_interfaces_lowlevel, except_this_one => abi_wrtout
 use abi_m_xmpi

#undef ABI_FUNC
#define ABI_FUNC 'abi_wrtout'

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: mode_paral

!Local variables-------------------------------
 integer :: comm,me,nproc
 integer,save :: master=0
 character(len=len(msg)+50) :: string
 character(len=500) :: my_mode_paral

!******************************************************************

 if ((unit == std_out).and.(.not.do_write_log)) RETURN
 if (unit == dev_null) RETURN

 my_mode_paral = "COLL"; if (PRESENT(mode_paral)) my_mode_paral = mode_paral

!Communicator is MPI_COMM_WORLD by default
 if (abinit_comm_output/=-1) then
   comm=abinit_comm_output
 else
   comm=abi_xmpi_comm_world
 end if

!Determine who I am in COMM_WORLD
 nproc = abi_xmpi_comm_size(comm)
 me    = abi_xmpi_comm_rank(comm)

 if( (my_mode_paral=='COLL') .or. (nproc==1) ) then
   if (me==master) then
     call abi_wrtout_myproc(unit, msg)
   end if

 else if (my_mode_paral=='PERS') then
   call abi_write_lines(unit,msg)

 else if (my_mode_paral=='INIT') then
   master=unit

 else
   write(string,'(7a)')ch10,&
&   'abi_wrtout: ERROR -',ch10,&
&   '  Unknown write mode: ',my_mode_paral,ch10,&
&   '  Continuing anyway ...'
   write(unit, '(A)' ) trim(string)
 end if

end subroutine abi_wrtout
!!***

!!****f* ABINIT/abi_wrtout_myproc
!! NAME
!!  abi_wrtout_myproc
!!
!! FUNCTION
!!  Do the output for one proc. For parallel or sequential output use abi_wrtout()
!!  instead. Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!!  Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR)
!! INPUTS
!!  unit=unit number for writing
!!  message=(character(len=*)) message to be written
!!  [mpicomm]= Optional argument. If present, no printing is done
!!             Variables iexit, nwarning and ncomment are
!!             summed over the mpicomm communicator
!!
!! OUTPUT
!!  (only writing)
!!
!! SOURCE

subroutine abi_wrtout_myproc(unit,message,mpicomm) ! optional argument

 use abi_defs_basis
 use abi_interfaces_lowlevel, except_this_one => abi_wrtout_myproc
 use abi_m_xmpi

#undef ABI_FUNC
#define ABI_FUNC 'abi_wrtout_myproc'

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: message
 integer,intent(in),optional :: mpicomm

!Local variables-------------------------------
!scalars
 integer,save :: iexit=0,ncomment=0,nwarning=0
 integer :: ierr
 logical :: print_std_err
!arrays
 integer :: buf(3)

!******************************************************************

!When I/O are redirected, it is sometimes necessary to reduce counters (saved) values;
!this can be done by passing mpicomm optional argument to the routine In that case, no printing is done.
 if (present(mpicomm)) then
   buf(1)=iexit;buf(2)=ncomment;buf(3)=nwarning
   call abi_xmpi_sum(buf,mpicomm,ierr)
   iexit=buf(1);ncomment=buf(2);nwarning=buf(3)
   if (iexit/=0) iexit=1
   return
 end if

 print_std_err=(unit==std_out.and.(index(trim(message),'BUG')/=0.or.index(trim(message),'ERROR')/=0))

 call abi_write_lines(unit,message)
 if (print_std_err) then
   call abi_write_lines(std_err,message)
 end if

 if( index(trim(message),'BUG') /= 0 )then
   write(unit, '(a)' ) '  Action : contact ABINIT group.'
   if (print_std_err) write(std_err, '(a)' ) '  Action : contact ABINIT group.'
   write(unit,*)
   if (print_std_err) write(std_err,*)
 end if

 if( index(trim(message),'BUG') /= 0   .or. index(trim(message),'Calculation completed') /= 0 )then
   if(nwarning<10000 .and. ncomment<1000)then
     write(unit, '(a,i5,a,i4,a)' ) '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   else
     write(unit, '(a,i6,a,i6,a)' ) '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   end if
   if(iexit/=0)then
     write(unit, '(a)' ) ' Note : exit requested by the user.'
   end if
 end if

 if( index(trim(message),'Exit') /= 0 )then
   iexit=1
 end if

!Count the number of warnings and comments. Only take into
!account unit std_out, in order not to duplicate these numbers.
 if( index(trim(message),'WARNING') /= 0 .and. unit==std_out )then
   nwarning=nwarning+1
 end if
 if( index(trim(message),'COMMENT') /= 0 .and. unit==std_out )then
   ncomment=ncomment+1
 end if

end subroutine abi_wrtout_myproc
!!***

!!****f* ABINIT/abi_write_lines
!! NAME
!!  abi_write_lines
!!
!! FUNCTION
!!  This routine receives a string, split the message in lines according to the 
!!  ch10 character and output the text to the specified unit 
!!
!! INPUTS
!!  unit=unit number for writing
!!  message=(character(len=*)) message to be written
!!
!! OUTPUT
!!  Only writing.
!!
!! SOURCE

subroutine abi_write_lines(unit,message)

 use abi_defs_basis

#undef ABI_FUNC
#define ABI_FUNC 'abi_write_lines'

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: message

!Local variables-------------------------------
!scalars
 integer :: msg_size,ii,jj,rtnpos

!******************************************************************

 msg_size = len_trim(message)

 if (msg_size == 0) then
   write(unit,*)
   return 
 end if

 ! Here, split the message, according to the char(10) characters (carriage return). 
 ! This technique is portable accross different OS.
 rtnpos = index(message,ch10)

 if (rtnpos == 0) then
   write(unit,"(a)")message(1:msg_size)
   return
 end if 

 ii = 1; jj = rtnpos
 do 
   if (ii == jj) then
     write(unit,*)
   else
     write(unit, '(a)' ) message(ii:jj-1)
   end if
   ii = jj + 1
   if (ii > msg_size) exit
   jj = index(message(ii:msg_size),ch10) 
   if (jj == 0) then 
     ! Will write the last line at the next iteration and exit .
     jj = msg_size + 1
   else
     jj = jj + ii - 1
   end if
   !write(*,*)"ii, jj, msg_size",ii, jj, msg_size
 end do

 ! This is needed to preserve the od behaviour: a ch10 at the 
 ! end of the string was causing an extra newline!
 if (message(msg_size:msg_size) == ch10) write(unit,*)

end subroutine abi_write_lines
!!***
