!{\src2tex{textfont=tt}}
!!****f* ABINIT/xallgatherv_mpi_int2d
!! NAME
!!  xallgatherv_mpi_int2d
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the MPI CPP flags.
!!  xallgatherv_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group (AR,XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allgatherv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xallgatherv_mpi_int2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

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
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHERV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
&   MPI_INTEGER,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_int2d
!!***

!!****f* ABINIT/xallgatherv_mpi_int
!! NAME
!!  xallgatherv_mpi_int
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
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
!!      mpi_allgatherv
!!
!! SOURCE
subroutine xallgatherv_mpi_int(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(in) :: xval(:)
 integer,intent(inout)   :: recvbuf(:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHERV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
&   MPI_INTEGER,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_int
!!***

!!****f* ABINIT/xallgatherv_mpi_dp
!! NAME
!!  xallgatherv_mpi_dp
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
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
!!      mpi_allgatherv
!!
!! SOURCE
subroutine xallgatherv_mpi_dp(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

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
 real(dp),intent(inout)   :: recvbuf(:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_dp
!!***

!!****f* ABINIT/xallgatherv_mpi_dp2d
!! NAME
!!  xallgatherv_mpi_dp2d
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
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
!!      mpi_allgatherv
!!
!! SOURCE
subroutine xallgatherv_mpi_dp2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

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
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_dp2d
!!***

!!****f* ABINIT/xallgatherv_mpi_dp3d
!! NAME
!!  xallgatherv_mpi_dp3d
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: double precision three-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
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
!!      mpi_allgatherv
!!
!! SOURCE
subroutine xallgatherv_mpi_dp3d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:,:)
 real(dp),intent(inout)   :: recvbuf(:,:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_dp3d
!!***

!!****f* ABINIT/xallgatherv_mpi_dp4d
!! NAME
!!  xallgatherv_mpi_dp4d
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: double precision four-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
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
!!      mpi_allgatherv
!!
!! SOURCE
subroutine xallgatherv_mpi_dp4d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:,:,:)
 real(dp),intent(inout)   :: recvbuf(:,:,:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_dp4d
!!***
