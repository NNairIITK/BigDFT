!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi_xmpi
!! NAME
!!  m_abi_xmpi
!!
!! FUNCTION
!!  This module provides a few MPI named constants and generic interfaces
!!  for inquiring the MPI environment.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (MG, MT, ...)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE m_abi_xmpi

 use abi_defs_basis
#ifdef HAVE_MPI2
 use mpi
#endif

 implicit none

 private
!!***

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

!MPI constants
#ifdef HAVE_MPI
 integer,public,parameter :: abi_xmpi_comm_world = MPI_COMM_WORLD
 integer,public,parameter :: abi_xmpi_comm_self  = MPI_COMM_SELF
 integer,public,parameter :: abi_xmpi_comm_null  = MPI_COMM_NULL
#else
 ! Fake replacements for the sequential version.
 integer,public,parameter :: abi_xmpi_comm_world = 0
 integer,public,parameter :: abi_xmpi_comm_self  = 0
 integer,public,parameter :: abi_xmpi_comm_null  = 0
#endif

!----------------------------------------------------------------------
!!***

! Public procedures.
 public :: abi_xmpi_abort                 ! Wrapper for MPI_ABORT
 public :: abi_xmpi_comm_rank             ! Wrapper for MPI_COMM_RANK
 public :: abi_xmpi_comm_size             ! Wrapper for MPI_COMM_SIZE

 public :: abi_xmpi_sum                   ! Wrapper for MPI_ALLREDUCE(SUM)
interface abi_xmpi_sum
  module procedure abi_xmpi_sum_int1d
  module procedure abi_xmpi_sum_dp1d
  module procedure abi_xmpi_sum_dp3d
end interface abi_xmpi_sum
!!***

CONTAINS  !===========================================================
!!***

!!****f* m_abi_xmpi/abi_xmpi_abort
!! NAME
!!  abi_xmpi_abort
!!
!! FUNCTION
!!  Wrapper for MPI_ABORT
!!
!! INPUTS
!!  [comm]=communicator of tasks to abort.
!!  [mpierr]=Error code to return to invoking environment.
!!  [msg]=User message
!!  [exit_status]=optional, shell return code, default 1
!!
!! SOURCE

subroutine abi_xmpi_abort(comm,mpierr,msg,exit_status)

#undef ABI_FUNC
#define ABI_FUNC 'abi_xmpi_abort'

 implicit none

!Arguments-------------------------
 integer,optional,intent(in) :: comm,mpierr,exit_status
 character(len=*),optional,intent(in) :: msg

!Local variables-------------------
 integer :: ierr,my_comm,my_errorcode,ilen,ierr2
 logical :: testopen
#ifdef HAVE_MPI
 character(len=MPI_MAX_ERROR_STRING) :: mpi_msg_error
#endif

! *************************************************************************

 ierr=0
 my_comm = abi_xmpi_comm_world; if (PRESENT(comm)) my_comm = comm

 if (PRESENT(msg)) then
   write(std_out,'(2a)')"User message: ",TRIM(msg)
 end if

 ! Close std_out and ab_out
 inquire(std_out,opened=testopen)
 if (testopen) close(std_out)

 inquire(ab_out,opened=testopen)
 if (testopen) close(ab_out)

#ifdef HAVE_MPI
 my_errorcode=MPI_ERR_UNKNOWN; if (PRESENT(mpierr)) my_errorcode=mpierr
 call MPI_ERROR_STRING(my_errorcode, mpi_msg_error, ilen, ierr2)
 call MPI_ABORT(my_comm,my_errorcode,ierr)
#endif

#if defined FC_NAG
 call exit(exit_status)
#elif defined HAVE_FC_EXIT
 call exit(exit_status)
#else
 if (exit_status== 0) stop  "0"
 if (exit_status== 1) stop  "1"
 if (exit_status==-1) stop "-1"
#endif
 stop "1"

end subroutine abi_xmpi_abort
!!***

!----------------------------------------------------------------------

!!****f* m_abi_xmpi/abi_xmpi_comm_rank
!! NAME
!!  abi_xmpi_comm_rank
!!
!! FUNCTION
!!  Wrapper for MPI_COMM_RANK
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  abi_xmpi_comm_rank=The rank of the node inside comm
!!
!! SOURCE

function abi_xmpi_comm_rank(comm)

#undef ABI_FUNC
#define ABI_FUNC 'abi_xmpi_comm_rank'

 implicit none

!Arguments-------------------------
 integer,intent(in) :: comm
 integer :: abi_xmpi_comm_rank

!Local variables-------------------
 integer :: mpierr

! *************************************************************************

 mpierr=0
#ifdef HAVE_MPI
 abi_xmpi_comm_rank=-1  ! Return non-sense value if the proc does not belong to the comm
 if (comm/=abi_xmpi_comm_null) then
   call MPI_COMM_RANK(comm,abi_xmpi_comm_rank,mpierr)
 end if
#else
 abi_xmpi_comm_rank=0
#endif

end function abi_xmpi_comm_rank
!!***

!----------------------------------------------------------------------

!!****f* m_abi_xmpi/abi_xmpi_comm_size
!! NAME
!!  abi_xmpi_comm_size
!!
!! FUNCTION
!!  Wrapper for MPI_COMM_SIZE
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  abi_xmpi_comm_size=The number of processors inside comm.
!!
!! SOURCE

function abi_xmpi_comm_size(comm)

#undef ABI_FUNC
#define ABI_FUNC 'abi_xmpi_comm_size'

 implicit none

!Arguments-------------------------
 integer,intent(in) :: comm
 integer :: abi_xmpi_comm_size

!Local variables-------------------------------
!scalars
 integer :: mpierr

! *************************************************************************

 mpierr=0; abi_xmpi_comm_size=1
#ifdef HAVE_MPI
 if (comm/=abi_xmpi_comm_null) then
   call MPI_COMM_SIZE(comm,abi_xmpi_comm_size,mpierr)
 end if
#endif

end function abi_xmpi_comm_size
!!***

!----------------------------------------------------------------------

!!****f* m_abi_xmpi/abi_xmpi_sum_int1d
!! NAME
!!  abi_xmpi_sum_int1d
!!
!! FUNCTION
!!  Wrapper for MPI_ALLREDUCE(SUM)
!!  Target: 1D integer array
!!
!! SOURCE

subroutine abi_xmpi_sum_int1d(xval,comm,mpierr)

#undef ABI_FUNC
#define ABI_FUNC 'abi_xmpi_sum_int1d'

 implicit none

!Arguments ------------------------------------
 integer,intent(inout) :: xval(:)
 integer,intent(in) :: comm
 integer,intent(out) :: mpierr

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: xsum(size(xval))
#endif

! *************************************************************************
 mpierr=0
#if defined HAVE_MPI
 if (comm /= abi_xmpi_comm_self .and. comm /= abi_xmpi_comm_null) then
   call MPI_ALLREDUCE(xval,xsum,size(xval),MPI_INTEGER,MPI_SUM,comm,mpierr)
   xval = xsum
 end if
#endif
end subroutine abi_xmpi_sum_int1d
!!***

!!****f* m_abi_xmpi/abi_xmpi_sum_dp1d
!! NAME
!!  abi_xmpi_sum_dp1d
!!
!! FUNCTION
!!  Wrapper for MPI_ALLREDUCE(SUM)
!!  Target: 1D double precision array
!!
!! SOURCE

subroutine abi_xmpi_sum_dp1d(xval,comm,mpierr)

#undef ABI_FUNC
#define ABI_FUNC 'abi_xmpi_sum_dp1d'

 implicit none

!Arguments ------------------------------------
 real(dp),intent(inout) :: xval(:)
 integer,intent(in) :: comm
 integer,intent(out) :: mpierr

!Local variables-------------------------------
#if defined HAVE_MPI
 real(dp) :: xsum(size(xval,dim=1))
#endif

! *************************************************************************
 mpierr=0
#if defined HAVE_MPI
 if (comm /= abi_xmpi_comm_self .and. comm /= abi_xmpi_comm_null) then
   call MPI_ALLREDUCE(xval,xsum,size(xval),MPI_DOUBLE_PRECISION,MPI_SUM,comm,mpierr)
   xval = xsum
 end if
#endif
end subroutine abi_xmpi_sum_dp1d
!!***

!!****f* m_abi_xmpi/abi_xmpi_sum_dp3d
!! NAME
!!  abi_xmpi_sum_dp3d
!!
!! FUNCTION
!!  Wrapper for MPI_ALLREDUCE(SUM)
!!  Target: 3D double precision array
!!
!! SOURCE

subroutine abi_xmpi_sum_dp3d(xval,comm,mpierr)

#undef ABI_FUNC
#define ABI_FUNC 'abi_xmpi_sum_dp3d'

 implicit none

!Arguments ------------------------------------
 real(dp),intent(inout) :: xval(:,:,:)
 integer,intent(in) :: comm
 integer,intent(out) :: mpierr

!Local variables-------------------------------
#if defined HAVE_MPI
 real(dp) :: xsum(size(xval,dim=1),size(xval,dim=2),size(xval,dim=3))
#endif

! *************************************************************************
 mpierr=0
#if defined HAVE_MPI
 if (comm /= abi_xmpi_comm_self .and. comm /= abi_xmpi_comm_null) then
   call MPI_ALLREDUCE(xval,xsum,size(xval),MPI_DOUBLE_PRECISION,MPI_SUM,comm,mpierr)
   xval = xsum
 end if
#endif
end subroutine abi_xmpi_sum_dp3d
!!***

!----------------------------------------------------------------------

END MODULE m_abi_xmpi
!!***
