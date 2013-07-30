!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_errors
!! NAME
!!  m_errors
!!
!! FUNCTION
!!  This module contains low-level procedures to check assertions and handle errors.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2013 ABINIT group (MG,YP,NCJ,MT)
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


MODULE m_errors

 use defs_basis
 use m_xmpi
!#ifdef HAVE_TRIO_NETCDF
! use netcdf
!#endif
!#ifdef HAVE_TRIO_ETSF_IO
! use etsf_io_low_level
! use etsf_io
!#endif
#ifdef HAVE_MPI2
 use mpi
#endif
!#ifdef FC_NAG
! use f90_unix
!#endif

 use m_build_info,      only : dump_config
 use m_io_tools,        only : flush_unit
 use m_fstrings,        only : toupper, basename

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 private

 public :: assert_eq        ! Report and die gracefully if integers not all equal (used for size checking).
 public :: assert           ! Report and die if any logical is false (used for argument range checking).
 public :: sentinel         ! Announce the entering or the exiting from a procedure.
 public :: die              ! Stop execution in case of unexpected events.
 public :: memerr           ! Stop execution when memory allocation failed.
 public :: msg_hndl         ! Basic Error handlers.
 public :: netcdf_check     ! Stop execution after a NetCDF I/O error
 public :: check_mpi_ierr   ! Erro handler for MPI routines.
 public :: io_hndl          ! Error handler for IO operations on external files.
 public :: unused_var       ! Helper function used to silence compiler warnings due to unused variables.
!#if defined HAVE_TRIO_ETSF_IO
! public :: abietsf_msg_hndl ! Error handler for ETSF-IO routines.
! public :: abietsf_warn     ! Write warnings reported by ETSF-IO routines.
!#endif

 interface assert_eq  
   module procedure assert_eq2
   module procedure assert_eq3
   module procedure assert_eq4
   module procedure assert_eqn
 end interface assert_eq

 interface assert 
   module procedure assert1
   module procedure assert2
   module procedure assert3
   module procedure assert4
   module procedure assert_v
 end interface assert

 interface unused_var
   module procedure unused_int
   module procedure unused_int_d1
   module procedure unused_real_dp
   module procedure unused_real_sp_d1
   module procedure unused_real_dp_d1
   module procedure unused_cplx_dpc
   module procedure unused_cplx_dpc_d1
   module procedure unused_cplx_spc
   module procedure unused_cplx_spc_d1
   module procedure unused_logical_d0
   module procedure unused_ch_d0
 end interface unused_var

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eq2
!! NAME
!!  assert_eq2
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! INPUTS 
!!  l1,l2,.. Integers to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! OUTPUT
!!
!! PARENTS 
!! 
!! CHILDREN
!!
!! SOURCE

function assert_eq2(l1,l2,message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert_eq2'
!End of the abilint section

 integer,intent(in) :: l1,l2 
 integer,optional,intent(in) :: line
 integer :: assert_eq2
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'

! *************************************************************************

 if (l1==l2) then
  assert_eq2=l1
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,line)
 end if

end function assert_eq2
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eq3
!! NAME
!!  assert_eq3
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! INPUTS 
!!  l1,l2,.. Integers to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! OUTPUT
!!
!! PARENTS 
!! 
!! CHILDREN
!!
!! SOURCE

function assert_eq3(l1,l2,l3,message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert_eq3'
!End of the abilint section

 integer,intent(in) :: l1,l2,l3 
 integer,optional,intent(in) :: line
 integer :: assert_eq3
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (l1==l2.and.l2==l3) then
  assert_eq3=l1
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,line)
 end if

end function assert_eq3
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eq4
!! NAME
!!  assert_eq4
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! INPUTS 
!!  l1,l2,.. Integers to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! OUTPUT
!!
!! PARENTS 
!! 
!! CHILDREN
!!
!! SOURCE

function assert_eq4(l1,l2,l3,l4,message,file,line)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert_eq4'
!End of the abilint section

 integer,intent(in) :: l1,l2,l3,l4 
 integer,optional,intent(in) :: line
 integer :: assert_eq4
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (l1==l2.and.l2==l3.and.l3==l4) then
  assert_eq4=l1
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,line)
 end if

end function assert_eq4
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eqn
!! NAME
!!  assert_eqn
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
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

function assert_eqn(nn,message,file,line)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert_eqn'
!End of the abilint section

 integer,optional,intent(in) :: line
 integer :: assert_eqn
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file 
!arrays
 integer,intent(in) :: nn(:)

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (ALL(nn(2:)==nn(1))) then
  assert_eqn=nn(1)
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,line)
 end if

end function assert_eqn
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert1
!! NAME
!!  assert1
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if 
!!  any logical is false (used for arg range checking).
!!
!! INPUTS 
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additiona information
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine assert1(l1,message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert1'
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.l1) then
   if (PRESENT(line)) f90line=line
   if (PRESENT(file)) f90name= basename(file)
   call msg_hndl(message,'ERROR','COLL',f90name,f90line)
 end if

end subroutine assert1
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert2
!! NAME
!!  assert2
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if 
!   any logical is false (used for arg range checking).
!!
!! INPUTS 
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE
subroutine assert2(l1,l2,message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert2'
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.(l1.and.l2)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,f90line)
 end if

end subroutine assert2
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert3
!! NAME
!!  assert3
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if 
!!  any logical is false (used for arg range checking).
!!
!! INPUTS 
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE
subroutine assert3(l1,l2,l3,message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert3'
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2,l3

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.(l1.and.l2.and.l3)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,f90line)
 end if

end subroutine assert3
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert4
!! NAME
!!  assert4
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if 
!!  any logical is false (used for arg range checking).
!!
!! INPUTS 
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE
subroutine assert4(l1,l2,l3,l4,message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert4'
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2,l3,l4

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.(l1.and.l2.and.l3.and.l4)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,f90line)
 end if

end subroutine assert4
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_v
!! NAME
!!  assert_v
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if 
!!  any logical is false (used for arg range checking).
!!
!! INPUTS 
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE
subroutine assert_v(n,message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assert_v'
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: n(:)

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.ALL(n)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','COLL',f90name,f90line)
 end if

end subroutine assert_v
!!***

!----------------------------------------------------------------------

!!****f* m_errors/memerr
!! NAME
!!  memerr
!!
!! FUNCTION
!!  Reports errors occurring when allocating memory.
!!
!! INPUTS 
!!  array_name= name of the array not properly allocated.
!!  array_size= size of the array.
!!  file_name= name of the file where the allocation was performed.
!!  file_line= line number in the file where the allocation was performed.
!!
!! NOTES
!!  This routine is usually interfaced with the macros defined in abi_common.h
!!  and uses this information to define a line offset.
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine memerr(array_name,array_size,file_name,file_line)

 use defs_basis

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memerr'
!End of the abilint section

 character(len=*),intent(in) :: array_name,file_name
 integer(i8b),intent(in) :: array_size
 integer,intent(in) :: file_line

!Local variables-------------------------------
 character(len=500) :: msg

 write(msg,'(a,a,a,a,a,i10)') &
&  '  Memory allocation failed',ch10, &
&  '  for array ',trim(array_name),' of size ',array_size
 call die(msg,file_name,file_line)

end subroutine memerr
!!***

!----------------------------------------------------------------------

!!****f* m_errors/netcdf_check
!! NAME
!!  netcdf_check
!!
!! FUNCTION
!!  Error handler for Netcdf calls.
!!
!! INPUTS 
!!  ncerr=Status error returned by the Netcdf library.
!!  msg=User-defined string with info on the action that was performed
!!  file= name of the file.
!!  line= line number.
!!
!! NOTES
!!  This routine is usually interfaced with the macros defined in abi_common.h
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine netcdf_check(ncerr,msg,file,line)

 use defs_basis
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'netcdf_check'
!End of the abilint section

 integer,intent(in) :: ncerr
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: file
 integer,optional,intent(in) :: line

!Local variables-------------------------------
 integer :: f90line
 character(len=500) :: f90name
 character(len=1024) :: nc_msg 
 character(len=2048) :: my_msg

! *************************************************************************

!#ifdef HAVE_TRIO_NETCDF
! if (ncerr /= NF90_NOERR) then
!   if (PRESENT(line)) then
!     f90line=line
!   else 
!     f90line=0
!   end if
!   if (PRESENT(file)) then 
!     f90name = basename(file)
!   else
!     f90name='Subroutine Unknown'
!   end if
!   !
!   ! Append Netcdf string to user-defined message.
!   write(nc_msg,'(a,3x,a)')'NetCDF library returned:',TRIM(nf90_strerror(ncerr))
!   my_msg = TRIM(msg) // TRIM(nc_msg)
!
!   call msg_hndl(my_msg,"ERROR","COLL",f90name,f90line)
! end if
!#endif
!
end subroutine netcdf_check
!!***

!----------------------------------------------------------------------

!!****f* m_errors/sentinel
!! NAME
!!  sentinel
!!
!! FUNCTION
!!  Announce the entering and the exiting from a function. Useful for poor-man debugging.
!!
!! INPUTS 
!!  level=1 when entering, 2 for exit.
!!  mode_paral= ['COLL'|'PERS'|'COLL_SILENT|PERS_SILENT'] 
!!   'COLL' and 'PERS' refer to the output mode used in wrtout to report the message.
!!   'COLL_SILENT' and 'PERS_SILENT' can be used if the procedure is called several times inside a loop.
!!   In this case sentinel will report only the first entry and the first exit using either 'COLL' or 'PERS' mode.
!!  funcname=Name of the procedure to be tested (TODO should be passed through ABI_func)
!!  [lineno]=Line number. Defaults to 0.
!!
!! NOTES
!!  This routine is usually interfaced with the macros defined in abi_common.h
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine sentinel(level,mode_paral,funcname,lineno)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sentinel'
 use interfaces_14_hidewrite
!End of the abilint section

 integer,intent(in) :: level
 integer,optional,intent(in) :: lineno
 character(len=*),intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: funcname

!Local variables-------------------------------
 integer,save :: level_save=0 
 integer :: ii
 integer :: f90line
 character(len=500),save :: funcname_save
 character(len=4) :: my_mode
 character(len=10) :: lnum
 character(len=500) :: my_funcname='Function Unknown'
 character(len=500) :: msg

! *********************************************************************

 if (toupper(mode_paral)=='COLL_SILENT'.or.toupper(mode_paral)=='PERS_SILENT') then
   ! * Silent mode, check if we are inside a loop.
   if (level==level_save .and. funcname==funcname_save) RETURN
   ii = index( toupper(mode_paral), '_SILENT')
   my_mode=toupper(mode_paral(1:ii-1))
 else ! * Normal mode.
   my_mode=mode_paral
 end if

 level_save   =level
 funcname_save=funcname

 if (my_mode/='COLL'.or.my_mode/='PERS') my_mode='COLL'
 if (PRESENT(funcname)) my_funcname = basename(funcname)

 f90line=0; if (PRESENT(lineno)) f90line=lineno
 write(lnum,"(i0)")f90line

 my_funcname= TRIM(my_funcname)//":"//TRIM(lnum)

 if (level==1) then 
  msg = ' '//TRIM(my_funcname)//' : enter'//ch10
 else if (level==2) then
  msg = ' '//TRIM(my_funcname)//' : exit '//ch10
 else 
  call die('Wrong level',__FILE__,__LINE__)
 end if

 call wrtout(std_out,msg,my_mode) 
 call flush_unit(std_out)

end subroutine sentinel
!!***

!----------------------------------------------------------------------

!!****f* m_errors/die
!! NAME
!!  die
!!
!! FUNCTION
!!  Stop smoothly the execution in case of unexpected events reporting the
!!  line number and the file name where the error occurred as well as the 
!!  MPI rank of the processor. This routine is usually interfaced through 
!!  some macro defined in abi_common.h
!!
!! INPUTS 
!!  message=String containing additional information on the nature of the problem
!!  line=Line number of the file where problem occurred
!!  f90name=Name of the f90 file containing the caller
!!
!! PARENTS
!!      m_errors,m_xc_vdw
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine die(message,file,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'die'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: rank 
 integer :: f90line=0
 character(len=10) :: lnum,strank
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=500) :: msg

! *********************************************************************

 if (PRESENT(line)) f90line=line
 write(lnum,"(i0)")f90line

 ! === Determine my rank inside MPI_COMM_WORLD ===
 rank = xcomm_rank(xmpi_world)
 write(strank,"(i0)")rank

 if (PRESENT(file)) f90name= basename(file)
 msg=TRIM(f90name)//':'//TRIM(lnum)//' P'//TRIM(strank)

 write(msg,'(a,2x,2a,2x,a)')ch10,&
& TRIM(msg),ch10,&
& TRIM(message)

 call wrtout(std_out,msg,'COLL') 
 call leave_new('COLL')

end subroutine die
!!***

!----------------------------------------------------------------------

!!****f* m_errors/msg_hndl
!! NAME
!!  msg_hndl
!!
!! FUNCTION
!!  Basic error handler for abinit. This routine is usually interfaced through some macro defined in abi_common.h
!!
!! INPUTS 
!!  message=string containing additional information on the nature of the problem
!!  level=string defining the type of problem. Possible values are
!!   COMMENT
!!   WARNING
!!   ERROR
!!   BUG
!!  line=line number of the file where problem occurred
!!  file=name of the f90 file containing the caller
!!  mode_paral=Either "COLL" or "PERS".
!!  NODUMP= (optional) if present dump config before stopping
!!  NOSTOP= (optional) if present don't stop even in the case of an error or a bug
!!
!! OUTPUT
!!
!! PARENTS
!!      m_errors
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine msg_hndl(message,level,mode_paral,file,line,NODUMP,NOSTOP)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'msg_hndl'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 integer,optional,intent(in) :: line
 logical,optional,intent(in) :: NODUMP,NOSTOP
 character(len=*),intent(in) :: level,message
 character(len=*),optional,intent(in) :: file
 character(len=*),intent(in) :: mode_paral

!Local variables-------------------------------
 integer :: f90line
 character(len=10) :: lnum
 character(len=500) :: f90name
 character(len=500) :: my_msg,sbuf

! *********************************************************************

 if (PRESENT(line)) then
   f90line=line
 else 
   f90line=0
 end if
 write(lnum,"(i0)")f90line

 if (PRESENT(file)) then 
   f90name = basename(file)
 else
   f90name='Subroutine Unknown'
 end if

 my_msg=TRIM(f90name)//":"//TRIM(lnum)//":"

 select case (toupper(level))

 case ('COMMENT','WARNING')
   write(sbuf,'(a,1x,3a,1x,a)')ch10,&
&    TRIM(my_msg),toupper(level),ch10,&
&    TRIM(message)
   call wrtout(std_out,sbuf,mode_paral) 

 case ('ERROR','BUG')

   if ((.not.present(NOSTOP)).and.(.not.present(NODUMP))) then
     call print_kinds()
     call xmpi_show_info()
     call dump_config()
   end if

   write(sbuf,'(a,1x,3a,1x,a)')ch10,&
&    TRIM(my_msg),toupper(level),ch10,&
&    TRIM(message)
   call wrtout(std_out,sbuf,mode_paral) 

   if (.not.present(NOSTOP)) then
     call leave_new(mode_paral)
     !call leave_new(mode_paral,print_config=.FALSE.)
     !call leave_new(mode_paral,print_config=.TRUE.)
   end if

 case default 
   write(sbuf,'(4a)')ch10,&
&    ' msg_hndl: BUG**2 - ',ch10,' Wrong value for level '
   call die(sbuf,&
&  __FILE__,__LINE__)
 end select

end subroutine  msg_hndl
!!***

!----------------------------------------------------------------------

!!****f* m_errors/io_hndl
!! NAME
!!  io_hndl
!!
!! FUNCTION
!!  Basic error handler for I/O operations on external files. This routine is usually interfaced 
!!  through some macro defined in abi_common.h
!!
!! INPUTS 
!!  unit=Fortran unit number
!!  ios=IO status
!!  f90name=name of the subroutine where error occurs
!!  line=line number in the f90name file
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine io_hndl(ios,unit,f90name,line)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_hndl'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 integer,intent(in) :: unit,ios,line
 character(len=*),intent(in) :: f90name

!Local variables-------------------------------
 character(len=fnlen) :: fname
 character(len=10) :: lnum,s_ios,s_unt
 character(len=500) :: msg
 logical :: lexist,lopened,lnamed
! *********************************************************************

 write(lnum,"(i0)")line
 msg=TRIM(f90name)//':'//TRIM(lnum)

 write(s_ios,"(i0)")ios
 write(s_unt,"(i0)")unit

 write(msg,'(8a)')ch10,&
& TRIM(msg),' I/O ERROR- ',ch10,&
& ' while operating on unit ',TRIM(s_unt),', iostat = ',TRIM(s_ios)
 call wrtout(std_out,msg,'COLL') 

 inquire(unit=unit,exist=lexist,named=lnamed,opened=lopened)
 fname='None' ; if (lnamed) inquire(unit=unit,name=fname)

 write(msg,'(2a,2(a,l7,a),2a)')&
& ' Inquire reports : ',ch10,&
& '  exist  = ',lexist,ch10,&
& '  opened = ',lopened,ch10,&
& '  name   = ',TRIM(fname)
 call wrtout(std_out,msg,'COLL') 
 
 call leave_new('COLL')

end subroutine io_hndl
!!***

!----------------------------------------------------------------------

!!****f* m_errors/check_mpi_ierr
!! NAME
!!  check_mpi_ierr
!!
!! FUNCTION
!!  Basic error handler for MPI calls. This routine is usually interfaced through some macro defined in abi_common.h
!!
!! INPUTS 
!!  ierr=Exit status reported by an MPI call.
!!  line=line number of the file where problem occurred
!!  file=name of the f90 file containing the caller
!!  mode_paral=Either "COLL" or "PERS".
!!
!! OUTPUT
!!  Write error message thep stop execution.
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine check_mpi_ierr(ierr,msg,mode_paral,file,line)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_mpi_ierr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ierr
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: msg,mode_paral
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer,parameter :: mpi_msg_len=1000
 integer :: f90line=0 
 integer :: ilen,ierr2
 character(len=10) :: lnum
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=mpi_msg_len) :: mpi_msg_error
 character(len=500) :: my_msg
! *********************************************************************

#ifdef HAVE_MPI
 if (ierr==MPI_SUCCESS) RETURN
 call MPI_ERROR_STRING(ierr, mpi_msg_error, ilen, ierr2)
#else
 ilen=0; ierr2=0
 mpi_msg_error = " Check_mpi_ierr should not be called in non-MPI mode!"
 if (ierr==0) RETURN
#endif

 if (ilen>mpi_msg_len) write(std_out,*)" Warning_ MPI message has been truncated!"
 if (ierr2/=0) write(std_out,*)" Warning: MPI_ERROR_STRING returned ierr2= ",ierr2

 if (PRESENT(line)) f90line=line
 write(lnum,"(i0)")f90line

 if (PRESENT(file)) f90name = basename(file)
 my_msg=" "//TRIM(f90name)//":"//TRIM(lnum)//":"//TRIM(msg)

 call print_kinds()
 call xmpi_show_info()
 call dump_config()

 call wrtout(std_out,my_msg,mode_paral) 
 call wrtout(std_out,mpi_msg_error,mode_paral) 

 !call leave_new('PERS',print_config=.FALSE.)
 call leave_new('PERS')

end subroutine check_mpi_ierr
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_int
!! NAME
!!  unused_int
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS 
!!  var=Scalar integer value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_int(var) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_int'
!End of the abilint section

 integer,intent(in) :: var

!Local variables-------------------------------
 integer :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_int
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_int_d1
!! NAME
!!  unused_int_d1
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro. Target: one-dimensional integer vector.
!!
!! INPUTS 
!!  var_arr=Vector of integer values.
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_int_d1(var_arr) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_int_d1'
!End of the abilint section

 integer,intent(in) :: var_arr(:)

!Local variables-------------------------------
 integer :: dummy(SIZE(var_arr))
! *********************************************************************

 dummy = var_arr

end subroutine unused_int_d1
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_real_dp
!! NAME
!!  unused_real_dp
!!
!! FUNCTION
!!  Helper function used to silence warning messages due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS 
!!  var=Scalar real value.
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_real_dp(var) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_real_dp'
!End of the abilint section

 real(dp),intent(in) :: var

!Local variables-------------------------------
 real(dp) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_real_dp
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_real_sp_d1
!! NAME
!!  unused_real_sp_d1
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro. Target: one-dimensional real(dp) vector.
!!
!! INPUTS 
!!  var_arr=Vector of real(dp) values.
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_real_sp_d1(var_arr) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_real_sp_d1'
!End of the abilint section

 real(sp),intent(in) :: var_arr(:)

!Local variables-------------------------------
 real(sp) :: dummy(SIZE(var_arr))
! *********************************************************************

 dummy = var_arr

end subroutine unused_real_sp_d1
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_real_dp_d1
!! NAME
!!  unused_real_dp_d1
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro. Target: one-dimensional real(dp) vector.
!!
!! INPUTS 
!!  var_arr=Vector of real(dp) values.
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_real_dp_d1(var_arr) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_real_dp_d1'
!End of the abilint section

 real(dp),intent(in) :: var_arr(:)

!Local variables-------------------------------
 real(dp) :: dummy(SIZE(var_arr))
! *********************************************************************

 dummy = var_arr

end subroutine unused_real_dp_d1
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_cplx_spc
!! NAME
!!  unused_cplx_spc
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS 
!!  var=Scalar complex value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_cplx_spc(var) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_cplx_spc'
!End of the abilint section

 complex(spc),intent(in) :: var

!Local variables-------------------------------
 complex(spc) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_cplx_spc
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_cplx_spc_d1
!! NAME
!!  unused_cplx_spc_d1
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro. Target: one-dimensional complex(spc) vector.
!!
!! INPUTS 
!!  var_arr=Vector of complex(spc) values.
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_cplx_spc_d1(var_arr) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_cplx_spc_d1'
!End of the abilint section

 complex(spc),intent(in) :: var_arr(:)

!Local variables-------------------------------
 complex(spc) :: dummy(SIZE(var_arr))
! *********************************************************************

 dummy = var_arr

end subroutine unused_cplx_spc_d1
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_cplx_dpc
!! NAME
!!  unused_cplx_dpc
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS 
!!  var=Scalar complex value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_cplx_dpc(var) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_cplx_dpc'
!End of the abilint section

 complex(dpc),intent(in) :: var

!Local variables-------------------------------
 complex(dpc) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_cplx_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_cplx_dpc_d1
!! NAME
!!  unused_cplx_dpc_d1
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro. Target: one-dimensional complex(dpc) vector.
!!
!! INPUTS 
!!  var_arr=Vector of complex(dpc) values.
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_cplx_dpc_d1(var_arr) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_cplx_dpc_d1'
!End of the abilint section

 complex(dpc),intent(in) :: var_arr(:)

!Local variables-------------------------------
 complex(dpc) :: dummy(SIZE(var_arr))
! *********************************************************************

 dummy = var_arr

end subroutine unused_cplx_dpc_d1
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_logical_d0
!! NAME
!!  unused_logical_d0
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS 
!!  var=Scalar logical value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_logical_d0(var) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_logical_d0'
!End of the abilint section

 logical,intent(in) :: var

!Local variables-------------------------------
 logical :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_logical_d0
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_ch_d0
!! NAME
!!  unused_ch_d0
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.  
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS 
!!  var=Scalar character value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

subroutine unused_ch_d0(var) 

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unused_ch_d0'
!End of the abilint section

 character(len=*),intent(in) :: var

!Local variables-------------------------------
 character(len=LEN(var)) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_ch_d0
!!***

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abietsf_msg_hndl
!! NAME
!!  abietsf_msg_hndl
!!
!! FUNCTION
!!  Wrapper to interface the abinint error handlers with the error handling routines used in etsf-io.
!!  It is usually interfaced via the macro ETSF_* defined in abi_common.h
!!
!! INPUTS 
!!  lstat=Logical flag returned by etsf-io routines.
!!  Error_data<ETSF_io_low_error>=Structure storing the error returned by etsf-io calls.
!!  [line]=line number of the file where the problem occurred
!!  [file]=name of the f90 file containing the caller
!!  mode_paral=Either "COLL" or "PERS".
!!
!! OUTPUT
!!  Only writing, then the code is stopped.
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,msg_hndl
!!
!! SOURCE

!#if defined HAVE_TRIO_ETSF_IO
!
!subroutine abietsf_msg_hndl(lstat,Error_data,mode_paral,file,line)
!
!
!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'abietsf_msg_hndl'
!!End of the abilint section
!
! implicit none
!
!!Arguments ------------------------------------
! integer,optional,intent(in) :: line
! character(len=*),optional,intent(in) :: file
! character(len=*),intent(in) :: mode_paral
! logical,intent(in) :: lstat
! type(ETSF_io_low_error),intent(in) :: Error_data
!
!!Local variables-------------------------------
! integer :: f90line=0
! character(len=500) :: f90name='Subroutine Unknown'
! character(len=etsf_io_low_error_len) :: errmess
!! *********************************************************************
!
! if (lstat) RETURN
!
! if (PRESENT(line)) f90line=line
! if (PRESENT(file)) f90name = file
! call etsf_io_low_error_to_str(errmess,Error_data)
!
! call msg_hndl(errmess,"ERROR",mode_paral,f90name,f90line)
!
!end subroutine abietsf_msg_hndl
!!!***
!
!!----------------------------------------------------------------------
!
!!!****f* m_abi_etsf/abietsf_warn
!!! NAME
!!!  abietsf_warn
!!!
!!! FUNCTION
!!!  Wrapper to write warning messages, only used for ETSF_IO routines 
!!!  It is usually interfaced via the macro ETSF_WARN defined in abi_common.h
!!!
!!! INPUTS 
!!!  lstat=status error.
!!!  Error_data<ETSF_io_low_error>=Structure storing the error returned by etsf-io calls.
!!!  [line]=line number of the file where the problem occurred
!!!  [file]=name of the f90 file containing the caller
!!!  mode_paral=Either "COLL" or "PERS".
!!!
!!! OUTPUT
!!!  Only writing.
!!!
!!! PARENTS
!!!
!!! CHILDREN
!!!      etsf_io_low_error_to_str,msg_hndl
!!!
!!! SOURCE
!
!
!subroutine abietsf_warn(lstat,Error_data,mode_paral,file,line)
!
!
!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'abietsf_warn'
!!End of the abilint section
!
! implicit none
!
!!Arguments ------------------------------------
! integer,optional,intent(in) :: line
! logical,intent(in) :: lstat
! character(len=*),optional,intent(in) :: file
! character(len=*),intent(in) :: mode_paral
! type(ETSF_io_low_error),intent(in) :: Error_data
!
!!Local variables-------------------------------
! integer :: f90line=0
! character(len=500) :: f90name='Subroutine Unknown'
! character(len=etsf_io_low_error_len) :: errmess
!! *********************************************************************
!
! if (lstat) RETURN 
!
! if (PRESENT(line)) f90line=line
! if (PRESENT(file)) f90name = file
! call etsf_io_low_error_to_str(errmess,Error_data)
!
! call msg_hndl(errmess,"WARNING",mode_paral,f90name,f90line)
!
!end subroutine abietsf_warn
!!!***
!
!#endif

!----------------------------------------------------------------------

END MODULE m_errors
!!***
