!{\src2tex{textfont=tt}}
!!****f* ABINIT/xrecv_mpi
!! NAME
!!  xrecv_mpi
!!
!! FUNCTION
!!  This module contains functions that call MPI routine MPI_RECV,
!!  to receive data on one processor sent by another,
!!  if we compile the code using the MPI CPP flags.
!!  xrecv_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE
!!***

!!****f* ABINIT/xrecv_mpi_intv
!! NAME
!!  xrecv_mpi_intv
!!
!! FUNCTION
!!  Receives data from one processor sent by another.
!!  Target: single integer.
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_RECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xrecv_mpi_intv(xval,source,tag,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
  integer,intent(inout) :: xval
  integer,intent(in) :: source,tag,spaceComm
  integer,intent(out)   :: ier
!Local variables-------------------
#if defined HAVE_MPI
  integer :: my_tag
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   my_tag=MOD(tag,xmpi_tag_ub+1)
   call MPI_RECV(xval,1,MPI_INTEGER,source,my_tag,spaceComm,MPI_STATUS_IGNORE,ier)
 end if
#endif

 end subroutine xrecv_mpi_intv
!!***

!!****f* ABINIT/xrecv_mpi_dp2d
!! NAME
!!  xrecv_mpi_dp2d
!!
!! FUNCTION
!!  Receives data from one proc sent by another.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_RECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xrecv_mpi_dp2d(xval,source,tag,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:)
 integer ,intent(in) :: source,tag,spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,my_tag
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   my_tag=MOD(tag,xmpi_tag_ub+1)
   call MPI_RECV(xval,n1*n2,MPI_DOUBLE_PRECISION,source,my_tag,spaceComm,MPI_STATUS_IGNORE,ier)
 end if
#endif

end subroutine xrecv_mpi_dp2d
!!***

!!****f* ABINIT/xrecv_mpi_dp3d
!! NAME
!!  xrecv_mpi_dp3d
!!
!! FUNCTION
!!  Receives data from one proc sent by another.
!!  Target: double precision three-dimensional arrays.
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_RECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xrecv_mpi_dp3d(xval,source,tag,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: source,tag,spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,my_tag
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   my_tag=MOD(tag,xmpi_tag_ub+1)
   call MPI_RECV(xval,n1*n2*n3,MPI_DOUBLE_PRECISION,source,my_tag,spaceComm,MPI_STATUS_IGNORE,ier)
 end if
#endif

end subroutine xrecv_mpi_dp3d
!!***
