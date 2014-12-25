!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_leave_new
!! NAME
!!  abi_leave_new
!!
!! FUNCTION
!!  Routine for clean exit of f90 code, taking into account possible parallelization.
!!
!!  Note the this routine is private and should never be called explicitly. 
!!  Please, use the macros:
!!    MSG_ERROR, MSG_BUG 
!!  defined in abi_common.h to abort the execution.
!!  XG : this is not true, in very rare cases, ABINIT has to exit without giving an error (e.g. for non-zero prtkpt )
!!
!! INPUTS
!!  exit_status=(optional, default=1 or -1, see below) the return code of the routine
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same message to be
!!     written once only or
!!   'PERS' if the procs are calling the routine with different mesgs
!!     each to be written, or if one proc is calling the routine
!!  print_config=(optional, default=true)
!!       if true print out several informations before leaving
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

#include "abi_common.h"


subroutine abi_leave_new(mode_paral,exit_status,print_config)

 use defs_basis
 use m_xmpi

 !use m_build_info,      only : dump_config

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_leave_new'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=4),intent(in) :: mode_paral
 integer,intent(in),optional :: exit_status
 logical,intent(in),optional :: print_config

!Local variables-------------------------------
 logical :: print_config_
 !character(len=500) :: msg

! **********************************************************************

 call abi_wrtout(std_out,ch10//' abi_leave_new : decision taken to exit ...','PERS')

! Caveat: Do not use MPI collective calls!
 if (mode_paral == "COLL") then
   call abi_wrtout(std_out,"Why are you using COLL? Are you sure that ALL the processors are calling abi_leave_new?")
 end if

!Dump configuration before exiting
 print_config_=.true.; if (present(print_config)) print_config_=print_config
 if (print_config_) then
   call print_kinds()
   call xmpi_show_info()
   !call dump_config(std_out)
 end if

 if (present(exit_status)) then
   call xmpi_abort(exit_status=exit_status)
 else
   call xmpi_abort()
 end if

end subroutine abi_leave_new
!!***
