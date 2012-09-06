!{\src2tex{textfont=tt}}
!!****f* ABINIT/xalltoallv_mpi
!! NAME
!!  xalltoallv_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the   CPP flags.
!!  xalltoallv_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group (AR,XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xalltoallv_mpi_dp2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:)
 real(dp),intent(inout) :: recvbuf(:,:)
 integer ,intent(in) :: sendcnts(:),sdispls(:),rdispls(:),recvcnts(:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
! *********************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_DOUBLE_PRECISION,recvbuf,&
&   recvcnts,rdispls,MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xalltoallv_mpi_dp2d
!!***

!!****f* ABINIT/xalltoallv_mpi_int2d
!! NAME
!!  xalltoallv_mpi_int2d
!!
!! FUNCTION
!!  Sends data from all to all processes.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  sendcnts= number of sent elements
!!  sdispls= postions of values sent by the processor
!!  rdispls= positions of values received by the processor
!!  recvcnts= number of received elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xalltoallv_mpi_int2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(in) :: xval(:,:)
 integer,intent(inout) :: recvbuf(:,:)
 integer ,intent(in) :: sendcnts(:),sdispls(:),rdispls(:),recvcnts(:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
! *********************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_INTEGER,recvbuf,&
&   recvcnts,rdispls,MPI_INTEGER,spaceComm,ier)
 end if
#endif
end subroutine xalltoallv_mpi_int2d
!!***

!!****f* ABINIT/xalltoallv_mpi_dp1d
!! NAME
!!  xalltoallv_mpi_dp1d
!!
!! FUNCTION
!!  Sends data from all to all processes.
!!  Target: double precision one-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  sendcnts= number of sent elements
!!  sdispls= postions of values sent by the processor
!!  recvcnts= number of received elements
!!  spaceComm= MPI communicator
!!  rdispls= positions of values received by the processor
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xalltoallv_mpi_dp1d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:)
 real(dp),intent(inout) :: recvbuf(:)
 integer ,intent(in) :: sendcnts(:),sdispls(:),recvcnts(:)
 integer ,intent(in) :: spaceComm, rdispls
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer, allocatable :: rdispls_on(:)
#endif

! *********************************************************************

 ier=0
 if(.false.)write(6,*)rdispls
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   allocate(rdispls_on(size(sendcnts)))
   rdispls_on = 0
!  allgather xval on all proc. in spaceComm
   call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_DOUBLE_PRECISION,recvbuf,&
&   recvcnts,rdispls_on,MPI_DOUBLE_PRECISION,spaceComm,ier)
   deallocate(rdispls_on)
 end if
#endif
end subroutine xalltoallv_mpi_dp1d
!!***

!!****f* ABINIT/xalltoallv_mpi_dp1d2
!! NAME
!!  xalltoallv_mpi_dp1d2
!!
!! FUNCTION
!!  Sends data from all to all processes.
!!  Target: double precision one-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  sendcnts= number of sent elements
!!  sdispls= postions of values sent by the processor
!!  recvcnts= number of received elements
!!  spaceComm= MPI communicator
!!  rdispls= positions of values received by the processor
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xalltoallv_mpi_dp1d2(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:)
 real(dp),intent(inout) :: recvbuf(:)
 integer ,intent(in) :: sendcnts(:),sdispls(:),recvcnts(:),rdispls(:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *********************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_DOUBLE_PRECISION,recvbuf,&
&   recvcnts,rdispls,MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xalltoallv_mpi_dp1d2
!!***


