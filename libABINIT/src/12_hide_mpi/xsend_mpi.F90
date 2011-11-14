!{\src2tex{textfont=tt}}
!!****f* ABINIT/xsend_mpi
!! NAME
!!  xsend_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine MPI_SEND,
!!  to send data from one processor to another,
!!  if we compile the code using the MPI CPP flags.
!!  xsend_mpi is the generic function.
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

!!****f* ABINIT/xsend_mpi_intv
!! NAME
!!  xsend_mpi_intv
!!
!! FUNCTION
!!  Sends data from one processor to another.
!!  Target: single integer.
!!
!! INPUTS
!!  dest :: rank of destination process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xsend_mpi_intv(xval,dest,tag,spaceComm,ier)

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
  integer,intent(in) :: dest,tag,spaceComm
  integer,intent(out)   :: ier
!Local variables-------------------
#if defined HAVE_MPI
 integer :: my_tag
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   my_tag = MOD(tag,xmpi_tag_ub+1)
   call MPI_SEND(xval,1,MPI_INTEGER,dest,my_tag,spaceComm,ier)
 end if
#endif

 end subroutine xsend_mpi_intv
!!***

!!****f* ABINIT/xsend_mpi_dp2d
!! NAME
!!  xsend_mpi_dp2d
!!
!! FUNCTION
!!  Sends data from one proc to another.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  dest :: rank of destination process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsend_mpi_dp2d(xval,dest,tag,spaceComm,ier)

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
 integer ,intent(in) :: dest,tag,spaceComm
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
   my_tag = MOD(tag,xmpi_tag_ub+1)
   call MPI_SEND(xval,n1*n2,MPI_DOUBLE_PRECISION,dest,my_tag,spaceComm,ier)
 end if
#endif

end subroutine xsend_mpi_dp2d
!!***

!!****f* ABINIT/xsend_mpi_dp3d
!! NAME
!!  xsend_mpi_dp3d
!!
!! FUNCTION
!!  Sends data from one proc to another.
!!  Target: double precision three-dimensional arrays.
!!
!! INPUTS
!!  dest :: rank of destination process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsend_mpi_dp3d(xval,dest,tag,spaceComm,ier)

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
 integer ,intent(in) :: dest,tag,spaceComm
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
   my_tag = MOD(tag,xmpi_tag_ub+1)
   call MPI_SEND(xval,n1*n2*n3,MPI_DOUBLE_PRECISION,dest,my_tag,spaceComm,ier)
 end if
#endif

end subroutine xsend_mpi_dp3d
!!***
!
