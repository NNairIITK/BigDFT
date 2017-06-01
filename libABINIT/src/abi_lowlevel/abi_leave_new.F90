!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_leave_new
!! NAME
!!  abi_leave_new
!!
!! FUNCTION
!!  Routine for clean exit of f90 code, taking into account possible parallelization.
!!
!! INPUTS
!!  exit_status=(optional, default=1 or -1, see below) the return code of the routine
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same message to be
!!     written once only or
!!   'PERS' if the procs are calling the routine with different mesgs
!!     each to be written, or if one proc is calling the routine
!!
!! OUTPUT
!!  (only writing, then stop)
!!
!! NOTES
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_leave_new(mode_paral,exit_status)

 use abi_defs_basis
 use abi_interfaces_lowlevel, except_this_one => abi_leave_new
 use m_abi_xmpi

#undef ABI_FUNC
#define ABI_FUNC 'abi_leave_new'

 implicit none

!Arguments ------------------------------------
 character(len=4),intent(in) :: mode_paral
 integer,intent(in),optional :: exit_status

!Local variables-------------------------------
 !character(len=500) :: msg

! **********************************************************************

 call abi_wrtout(std_out,ch10//' abi_leave_new : decision taken to exit ...','PERS')

! Caveat: Do not use MPI collective calls!
 if (mode_paral == "COLL") then
   call abi_wrtout(std_out,"Why are you using COLL? Are you sure that ALL the processors are calling abi_leave_new?")
 end if

 if (present(exit_status)) then
   call abi_xmpi_abort(exit_status=exit_status)
 else
   call abi_xmpi_abort()
 end if

end subroutine abi_leave_new
!!***
