!{\src2tex{textfont=tt}}
!!****f* ABINIT/xmin_mpi
!! NAME
!!  xmin_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the MPI  CPP flags.
!!  xmin_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group (AR,XG,MB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xmin_mpi_intv(xval,xmin,spaceComm,ier)

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer ,intent(inout) :: xval,xmin
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_ALLREDUCE(xval,xmin,1,MPI_INTEGER,MPI_MIN,spaceComm,ier)
 else
#endif
   xmin=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xmin_mpi_intv
!!***

!!****f* ABINIT/xmin_mpi_dpv
!! NAME
!!  xmin_mpi_dpv
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: scalar double precisions.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!  xmin= number of elements in send buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xmin_mpi_dpv(xval,xmin,spaceComm,ier)
 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval,xmin
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_ALLREDUCE(xval,xmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,spaceComm,ier)
 else
#endif
   xmin=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xmin_mpi_dpv
!!***
