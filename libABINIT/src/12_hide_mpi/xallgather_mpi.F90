!{\src2tex{textfont=tt}}
!!****f* ABINIT/xallgather_mpi
!! NAME
!!  xallgather_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the  MPI CPP flags.
!!  xallgather_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group (AR,XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!!***

!!****f* ABINIT/xallgather_mpi_int
!! NAME
!!  xallgather_mpi_int
!!
!! FUNCTION
!!  Gathers data from all tasks and distributes it to all.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xallgather_mpi_int(xval,recvcounts,spaceComm,ier)

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
 integer,intent(inout) :: recvcounts(:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHER(xval,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,spaceComm,ier)
 else
#endif
   recvcounts(1)=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xallgather_mpi_int
!!***



!!****f* ABINIT/xallgather_mpi_char
!! NAME
!!  xallgather_mpi_char
!!
!! FUNCTION
!!  Gathers data from all tasks and distributes it to all.
!!  Target: one-dimensional character(20) arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  charval= buffer array
!!  recvcounts= number of received elements
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xallgather_mpi_char(charval,recvcounts,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(in)  :: spaceComm
 integer,intent(out) :: ier
 character(20),intent(inout) :: charval
 character(20),intent(inout) :: recvcounts(:)

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHER(charval,20,MPI_CHARACTER,recvcounts,20,MPI_CHARACTER,spaceComm,ier)
 else
#endif
   recvcounts = charval
#if defined HAVE_MPI
 end if
#endif
end subroutine xallgather_mpi_char
!!***


!!****f* ABINIT/xallgather_mpi_int1d
!! NAME
!!  xallgather_mpi_int1d
!!
!! FUNCTION
!!  Gathers data from all tasks and distributes it to all.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvcounts= number of received elements
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xallgather_mpi_int1d(xval,nelem,recvcounts,spaceComm,ier)

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
 integer,intent(inout) :: recvcounts(:)
 integer ,intent(in) :: nelem,spaceComm
 integer ,intent(out)   :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHER(xval,nelem,MPI_INTEGER,recvcounts,nelem,MPI_INTEGER,spaceComm,ier)
 else
#endif
   recvcounts(1:nelem)=xval(1:nelem)
#if defined HAVE_MPI
 end if
#endif
end subroutine xallgather_mpi_int1d
!!***

!!****f* ABINIT/xallgather_mpi_dp1d
!! NAME
!!  xallgather_mpi_dp1d
!!
!! FUNCTION
!!  Gathers data from all tasks and distributes it to all.
!!  Target: double precision one-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvcounts= number of received elements
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xallgather_mpi_dp1d(xval,nelem,recvcounts,spaceComm,ier)

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
 real(dp),intent(inout) :: recvcounts(:)
 integer ,intent(in) :: nelem,spaceComm
 integer ,intent(out)   :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvcounts,nelem,MPI_DOUBLE_PRECISION,spaceComm,ier)
 else
#endif
   recvcounts(1:nelem)=xval(1:nelem)
#if defined HAVE_MPI
 end if
#endif
end subroutine xallgather_mpi_dp1d
!!***

!!****f* ABINIT/xallgather_mpi_dp2d
!! NAME
!!  xallgather_mpi_dp2d
!!
!! FUNCTION
!!  Gathers data from all tasks and distributes it to all.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvcounts= number of received elements
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xallgather_mpi_dp2d(xval,nelem,recvcounts,spaceComm,ier)

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
 real(dp),intent(inout) :: recvcounts(:,:)
 integer ,intent(in) :: nelem,spaceComm
 integer ,intent(out)   :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvcounts,nelem,MPI_DOUBLE_PRECISION,spaceComm,ier)
 else
#endif
   recvcounts(:,:)=xval(:,:)
#if defined HAVE_MPI
 end if
#endif
end subroutine xallgather_mpi_dp2d
!!***

!!****f* ABINIT/xallgather_mpi_dp3d
!! NAME
!!  xallgather_mpi_dp3d
!!
!! FUNCTION
!!  Gathers data from all tasks and distributes it to all.
!!  Target: double precision three-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvcounts= number of received elements
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xallgather_mpi_dp3d(xval,nelem,recvcounts,spaceComm,ier)

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
   real(dp),intent(inout) :: recvcounts(:,:,:)
   integer ,intent(in) :: nelem,spaceComm
   integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvcounts,nelem,MPI_DOUBLE_PRECISION,spaceComm,ier)
 else
#endif
   recvcounts(:,:,:)=xval(:,:,:)
#if defined HAVE_MPI
 end if
#endif
end subroutine xallgather_mpi_dp3d
!!***

!!****f* ABINIT/xallgather_mpi_dp4d
!! NAME
!!  xallgather_mpi_dp4d
!!
!! FUNCTION
!!  Gathers data from all tasks and distributes it to all.
!!  Target: double precision four-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvcounts= number of received elements
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xallgather_mpi_dp4d(xval,nelem,recvcounts,spaceComm,ier)

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
   real(dp),intent(inout) :: recvcounts(:,:,:,:)
   integer ,intent(in) :: nelem,spaceComm
   integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  allgather xval on all proc. in spaceComm
   call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvcounts,nelem,MPI_DOUBLE_PRECISION,spaceComm,ier)
 else
#endif
   recvcounts(:,:,:,:)=xval(:,:,:,:)
#if defined HAVE_MPI
 end if
#endif
end subroutine xallgather_mpi_dp4d
!!***
