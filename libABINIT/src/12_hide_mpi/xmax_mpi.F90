!{\src2tex{textfont=tt}}
!!****f* ABINIT/xmax_mpi
!! NAME
!!  xmax_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the MPI  CPP flags.
!!  xmax_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group (AR,XG,MB)
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

subroutine xmax_mpi_intv(xval,xmax,spaceComm,ier)

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer ,intent(inout) :: xval,xmax
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_ALLREDUCE(xval,xmax,1,MPI_INTEGER,MPI_MAX,spaceComm,ier)
 else
#endif
   xmax=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xmax_mpi_intv
!!***

!!****f* ABINIT/xmax_mpi_dpv
!! NAME
!!  xmax_mpi_dpv
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
!!  xmax= number of elements in send buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xmax_mpi_dpv(xval,xmax,spaceComm,ier)
 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval,xmax
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_ALLREDUCE(xval,xmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,spaceComm,ier)
 else
#endif
   xmax=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xmax_mpi_dpv
!!***
